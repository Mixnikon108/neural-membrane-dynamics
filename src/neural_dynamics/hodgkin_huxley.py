"""Hodgkin-Huxley model of the squid giant axon (1952).

Original parameters at 6.3°C. Voltage convention: absolute membrane potential
(rest ≈ -65 mV), not the shifted variable used in the original paper.

State vector: y = [V, n, m, h]
    V : membrane potential (mV)
    n : K+ activation gating variable
    m : Na+ activation gating variable
    h : Na+ inactivation gating variable
"""

import numpy as np
import numba as nb

# Original HH parameters (modern voltage convention)
HH_DEFAULT_PARAMS = (
    1.0,     # C_m   : membrane capacitance (µF/cm²)
    120.0,   # g_Na  : max Na+ conductance (mS/cm²)
    36.0,    # g_K   : max K+ conductance (mS/cm²)
    0.3,     # g_L   : leak conductance (mS/cm²)
    50.0,    # E_Na  : Na+ reversal potential (mV)
    -77.0,   # E_K   : K+ reversal potential (mV)
    -54.4,   # E_L   : leak reversal potential (mV)
)


@nb.njit
def alpha_n(V):
    """K+ activation opening rate."""
    dV = V + 55.0
    if abs(dV) < 1e-7:
        return 0.1
    return 0.01 * dV / (1.0 - np.exp(-dV / 10.0))


@nb.njit
def beta_n(V):
    """K+ activation closing rate."""
    return 0.125 * np.exp(-(V + 65.0) / 80.0)


@nb.njit
def alpha_m(V):
    """Na+ activation opening rate."""
    dV = V + 40.0
    if abs(dV) < 1e-7:
        return 1.0
    return 0.1 * dV / (1.0 - np.exp(-dV / 10.0))


@nb.njit
def beta_m(V):
    """Na+ activation closing rate."""
    return 4.0 * np.exp(-(V + 65.0) / 18.0)


@nb.njit
def alpha_h(V):
    """Na+ inactivation opening rate (recovery from inactivation)."""
    return 0.07 * np.exp(-(V + 65.0) / 20.0)


@nb.njit
def beta_h(V):
    """Na+ inactivation closing rate."""
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))


@nb.njit
def n_inf(V):
    return alpha_n(V) / (alpha_n(V) + beta_n(V))

@nb.njit
def m_inf(V):
    return alpha_m(V) / (alpha_m(V) + beta_m(V))

@nb.njit
def h_inf(V):
    return alpha_h(V) / (alpha_h(V) + beta_h(V))

@nb.njit
def tau_n(V):
    return 1.0 / (alpha_n(V) + beta_n(V))

@nb.njit
def tau_m(V):
    return 1.0 / (alpha_m(V) + beta_m(V))

@nb.njit
def tau_h(V):
    return 1.0 / (alpha_h(V) + beta_h(V))


@nb.njit
def hh_rhs_parameterized(t, y, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L):
    """Full HH ODE system with explicit parameters.

    Parameters
    ----------
    y : array [V, n, m, h]
    I_ext : float, injected current (µA/cm²)

    Returns
    -------
    dydt : array [dV/dt, dn/dt, dm/dt, dh/dt]
    """
    V, n, m, h = y[0], y[1], y[2], y[3]

    I_Na = g_Na * m**3 * h * (V - E_Na)
    I_K = g_K * n**4 * (V - E_K)
    I_L = g_L * (V - E_L)

    dVdt = (I_ext - I_Na - I_K - I_L) / C_m
    dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n
    dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h

    return np.array([dVdt, dndt, dmdt, dhdt])


@nb.njit
def hh_rhs(t, y):
    """HH ODE with default parameters and I_ext = 0."""
    return hh_rhs_parameterized(t, y, 0.0, 1.0, 120.0, 36.0, 0.3, 50.0, -77.0, -54.4)


def make_hh_rhs(I_ext_fn, params=HH_DEFAULT_PARAMS):
    """Factory: create an ODE rhs with a time-dependent current function.

    Parameters
    ----------
    I_ext_fn : callable(t) -> float (must be njit-compiled)
    params : tuple of HH parameters

    Returns
    -------
    rhs : njit-compiled function(t, y) -> dydt
    """
    C_m, g_Na, g_K, g_L, E_Na, E_K, E_L = params

    @nb.njit
    def rhs(t, y):
        I_ext = I_ext_fn(t)
        return hh_rhs_parameterized(t, y, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)

    return rhs


def hh_initial_state(V0=-65.0):
    """Return initial state [V, n_inf, m_inf, h_inf] at given voltage."""
    return np.array([V0, n_inf(V0), m_inf(V0), h_inf(V0)])


@nb.njit
def _rk4_step_hh(t, y, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L, dt):
    """Single RK4 step for HH with explicit current value."""
    k1 = hh_rhs_parameterized(t, y, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
    k2 = hh_rhs_parameterized(t + 0.5*dt, y + 0.5*dt*k1, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
    k3 = hh_rhs_parameterized(t + 0.5*dt, y + 0.5*dt*k2, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
    k4 = hh_rhs_parameterized(t + dt, y + dt*k3, I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L)
    return y + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)


@nb.njit
def hh_solve_current_array(I_ext_array, dt, C_m=1.0, g_Na=120.0, g_K=36.0, g_L=0.3, E_Na=50.0, E_K=-77.0, E_L=-54.4):
    """Simulate HH model with an arbitrary current waveform.

    Parameters
    ----------
    I_ext_array : np.ndarray
        External current at each time step (µA/cm²). Length determines simulation duration.
    dt : float
        Time step (ms).

    Returns
    -------
    t : np.ndarray, shape (n_steps,)
    y : np.ndarray, shape (n_steps, 4) — columns are [V, n, m, h]
    """
    n_steps = len(I_ext_array)
    t = np.empty(n_steps)
    y = np.empty((n_steps, 4))

    V0 = -65.0
    y[0] = np.array([V0, n_inf(V0), m_inf(V0), h_inf(V0)])
    t[0] = 0.0

    for i in range(n_steps - 1):
        t[i + 1] = t[i] + dt
        y[i + 1] = _rk4_step_hh(t[i], y[i], I_ext_array[i], C_m, g_Na, g_K, g_L, E_Na, E_K, E_L, dt)

    return t, y


@nb.njit
def hh_solve_constant(I_ext, t_end, dt, C_m=1.0, g_Na=120.0, g_K=36.0, g_L=0.3, E_Na=50.0, E_K=-77.0, E_L=-54.4):
    """Simulate HH with constant injected current.

    Parameters
    ----------
    I_ext : float
        Constant current (µA/cm²).
    t_end : float
        Duration (ms).
    dt : float
        Time step (ms).

    Returns
    -------
    t : np.ndarray
    y : np.ndarray, shape (n_steps+1, 4)
    """
    n_steps = int(t_end / dt)
    t = np.empty(n_steps + 1)
    y = np.empty((n_steps + 1, 4))

    V0 = -65.0
    y[0] = np.array([V0, n_inf(V0), m_inf(V0), h_inf(V0)])
    t[0] = 0.0

    for i in range(n_steps):
        t[i + 1] = t[i] + dt
        y[i + 1] = _rk4_step_hh(t[i], y[i], I_ext, C_m, g_Na, g_K, g_L, E_Na, E_K, E_L, dt)

    return t, y
