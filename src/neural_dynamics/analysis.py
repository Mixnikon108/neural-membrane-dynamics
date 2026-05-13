"""Dynamical systems analysis tools for neural models.

Provides: numerical Jacobian, fixed point finding, eigenvalue classification,
nullcline computation, bifurcation sweeps, and f-I curve computation.
"""

import numpy as np
from scipy.optimize import fsolve
from .hodgkin_huxley import (
    hh_rhs_parameterized, HH_DEFAULT_PARAMS,
    alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h,
    n_inf, m_inf, h_inf,
    make_hh_rhs,
)
from .solvers import solve_ode
import numba as nb


def compute_jacobian_hh(y, I_ext, params=HH_DEFAULT_PARAMS, eps=1e-6):
    """Compute the Jacobian matrix of the HH system numerically.

    Uses central finite differences for accuracy.
    """
    n = len(y)
    J = np.zeros((n, n))
    C_m, g_Na, g_K, g_L, E_Na, E_K, E_L = params

    for j in range(n):
        y_plus = y.copy()
        y_minus = y.copy()
        y_plus[j] += eps
        y_minus[j] -= eps
        f_plus = hh_rhs_parameterized(0.0, y_plus, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
        f_minus = hh_rhs_parameterized(0.0, y_minus, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
        J[:, j] = (f_plus - f_minus) / (2.0 * eps)

    return J


def find_fixed_point(I_ext, params=HH_DEFAULT_PARAMS, V0_guess=-65.0):
    """Find a fixed point of the HH system for a given I_ext.

    Uses scipy.optimize.fsolve starting from steady-state guess.
    """
    C_m, g_Na, g_K, g_L, E_Na, E_K, E_L = params

    def equations(y):
        return hh_rhs_parameterized(0.0, y, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)

    y_guess = np.array([V0_guess, n_inf(V0_guess), m_inf(V0_guess), h_inf(V0_guess)])
    y_fp, info, ier, msg = fsolve(equations, y_guess, full_output=True)

    if ier != 1:
        raise RuntimeError(f"Fixed point solver did not converge: {msg}")

    return y_fp


def classify_fixed_point(y_fp, I_ext, params=HH_DEFAULT_PARAMS):
    """Classify a fixed point by its eigenvalues.

    Returns
    -------
    classification : str
        One of: "stable node", "unstable node", "saddle",
        "stable spiral", "unstable spiral", "center"
    eigenvalues : np.ndarray
    """
    J = compute_jacobian_hh(y_fp, I_ext, params)
    eigenvalues = np.linalg.eigvals(J)
    real_parts = np.real(eigenvalues)
    imag_parts = np.imag(eigenvalues)

    has_complex = np.any(np.abs(imag_parts) > 1e-8)
    all_negative = np.all(real_parts < 0)
    all_positive = np.all(real_parts > 0)
    mixed_sign = not all_negative and not all_positive

    if mixed_sign:
        classification = "saddle"
    elif has_complex and all_negative:
        classification = "stable spiral"
    elif has_complex and all_positive:
        classification = "unstable spiral"
    elif has_complex:
        classification = "center"
    elif all_negative:
        classification = "stable node"
    else:
        classification = "unstable node"

    return classification, eigenvalues


def compute_nullclines(V_range, I_ext, params=HH_DEFAULT_PARAMS, n_points=500):
    """Compute V-nullcline and n-nullcline for a 2D reduction (V, n).

    The V-nullcline uses the quasi-steady-state approximation m = m_inf(V), h = h_inf(V).

    Returns
    -------
    V_arr : np.ndarray
    n_V_nullcline : np.ndarray
        Values of n where dV/dt = 0
    n_n_nullcline : np.ndarray
        n_inf(V) — the n-nullcline
    """
    C_m, g_Na, g_K, g_L, E_Na, E_K, E_L = params
    V_arr = np.linspace(V_range[0], V_range[1], n_points)
    n_V_null = np.full(n_points, np.nan)
    n_n_null = np.zeros(n_points)

    for i, V in enumerate(V_arr):
        m = m_inf(V)
        h = h_inf(V)
        n_n_null[i] = n_inf(V)

        # dV/dt = 0 => I_ext - g_Na*m^3*h*(V-E_Na) - g_K*n^4*(V-E_K) - g_L*(V-E_L) = 0
        # Solve for n: g_K * n^4 * (V - E_K) = I_ext - g_Na*m^3*h*(V-E_Na) - g_L*(V-E_L)
        rhs_val = I_ext - g_Na * m**3 * h * (V - E_Na) - g_L * (V - E_L)
        denom = g_K * (V - E_K)

        if abs(denom) > 1e-10:
            n4 = rhs_val / denom
            if n4 >= 0:
                n_V_null[i] = n4 ** 0.25

    return V_arr, n_V_null, n_n_null


def _count_spikes(V, threshold=0.0):
    """Count spikes by detecting upward threshold crossings."""
    count = 0
    above = False
    for i in range(len(V)):
        if V[i] > threshold and not above:
            count += 1
            above = True
        elif V[i] < threshold:
            above = False
    return count


def _make_constant_current(I_val):
    """Create a njit-compiled constant-current function for use with make_hh_rhs."""
    # We need a njit function — use a closure approach compatible with Numba.
    # Since Numba does not support closures over non-constant Python scalars in
    # all contexts, we encode the value as a global via a module-level factory
    # that returns a fresh function each time.
    I_const = float(I_val)

    @nb.njit
    def _const_current(t):
        return I_const

    return _const_current


def compute_fi_curve_hh(I_range, n_points, t_sim, dt, params=HH_DEFAULT_PARAMS):
    """Compute the frequency-current (f-I) curve for the HH model.

    Parameters
    ----------
    I_range : tuple (I_min, I_max)
    n_points : int
    t_sim : float, simulation time in ms
    dt : float, time step

    Returns
    -------
    currents : np.ndarray
    rates : np.ndarray, firing rates in Hz
    """
    currents = np.linspace(I_range[0], I_range[1], n_points)
    rates = np.zeros(n_points)

    for i, I_ext in enumerate(currents):
        I_fn = _make_constant_current(I_ext)
        rhs = make_hh_rhs(I_fn, params=params)

        y0 = np.array([-65.0, n_inf(-65.0), m_inf(-65.0), h_inf(-65.0)])
        t, y = solve_ode(rhs, y0, t_span=(0.0, t_sim), dt=dt, method="rk4")

        # Skip transient (first 50 ms)
        skip = int(50.0 / dt)
        V = y[skip:, 0]
        n_spikes = _count_spikes(V)
        t_window = (t_sim - 50.0) / 1000.0  # convert ms to s
        rates[i] = n_spikes / t_window if t_window > 0 else 0.0

    return currents, rates


def bifurcation_sweep(I_range, n_points, params=HH_DEFAULT_PARAMS):
    """Sweep I_ext and classify the fixed point at each value.

    Returns
    -------
    currents : np.ndarray
    V_fps : np.ndarray, fixed point V at each current
    real_eigenvalues : np.ndarray, shape (n_points, 4), real parts of eigenvalues
    classifications : list of str
    """
    currents = np.linspace(I_range[0], I_range[1], n_points)
    V_fps = np.zeros(n_points)
    real_eigs = np.zeros((n_points, 4))
    classifications = []

    V_guess = -65.0
    for i, I_ext in enumerate(currents):
        try:
            y_fp = find_fixed_point(I_ext, params, V0_guess=V_guess)
            classification, eigs = classify_fixed_point(y_fp, I_ext, params)
            V_fps[i] = y_fp[0]
            real_eigs[i] = np.real(eigs)
            classifications.append(classification)
            V_guess = y_fp[0]  # warm start
        except RuntimeError:
            V_fps[i] = np.nan
            real_eigs[i] = np.nan
            classifications.append("no convergence")

    return currents, V_fps, real_eigs, classifications
