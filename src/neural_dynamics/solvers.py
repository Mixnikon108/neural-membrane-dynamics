"""Numerical ODE integrators with Numba JIT compilation.

Provides Euler and 4th-order Runge-Kutta methods for solving systems of ODEs.
All core stepping functions are compiled with Numba for performance.
"""

import numpy as np
import numba as nb


@nb.njit
def euler_step(rhs, t, y, dt):
    """Single forward Euler step: y_{n+1} = y_n + dt * f(t_n, y_n)."""
    return y + dt * rhs(t, y)


@nb.njit
def rk4_step(rhs, t, y, dt):
    """Single 4th-order Runge-Kutta step.

    Computes the weighted average of four slope evaluations:
        k1 = f(t, y)
        k2 = f(t + dt/2, y + dt/2 * k1)
        k3 = f(t + dt/2, y + dt/2 * k2)
        k4 = f(t + dt, y + dt * k3)
        y_{n+1} = y_n + (dt/6)(k1 + 2*k2 + 2*k3 + k4)
    """
    k1 = rhs(t, y)
    k2 = rhs(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = rhs(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = rhs(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


@nb.njit
def _solve_euler(rhs, y0, t_start, t_end, dt):
    n_steps = int((t_end - t_start) / dt)
    n_vars = y0.shape[0]
    t = np.empty(n_steps + 1)
    y = np.empty((n_steps + 1, n_vars))
    t[0] = t_start
    y[0] = y0.copy()
    for i in range(n_steps):
        y[i + 1] = euler_step(rhs, t[i], y[i], dt)
        t[i + 1] = t[i] + dt
    return t, y


@nb.njit
def _solve_rk4(rhs, y0, t_start, t_end, dt):
    n_steps = int((t_end - t_start) / dt)
    n_vars = y0.shape[0]
    t = np.empty(n_steps + 1)
    y = np.empty((n_steps + 1, n_vars))
    t[0] = t_start
    y[0] = y0.copy()
    for i in range(n_steps):
        y[i + 1] = rk4_step(rhs, t[i], y[i], dt)
        t[i + 1] = t[i] + dt
    return t, y


def solve_ode(rhs, y0, t_span, dt, method="rk4"):
    """Integrate a system of ODEs.

    Parameters
    ----------
    rhs : callable
        Right-hand side function f(t, y) -> dy/dt. Must be Numba-compatible.
    y0 : np.ndarray
        Initial state vector.
    t_span : tuple
        (t_start, t_end).
    dt : float
        Time step.
    method : str
        "euler" or "rk4".

    Returns
    -------
    t : np.ndarray, shape (n_steps+1,)
    y : np.ndarray, shape (n_steps+1, n_vars)
    """
    t_start, t_end = t_span
    y0 = np.asarray(y0, dtype=np.float64)
    if method == "euler":
        return _solve_euler(rhs, y0, t_start, t_end, dt)
    elif method == "rk4":
        return _solve_rk4(rhs, y0, t_start, t_end, dt)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'euler' or 'rk4'.")
