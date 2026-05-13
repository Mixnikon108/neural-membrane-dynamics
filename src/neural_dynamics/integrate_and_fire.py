"""Leaky Integrate-and-Fire (LIF) neuron model.

The simplest biologically motivated spiking neuron model:
    tau_m * dV/dt = -(V - V_rest) + R_m * I_ext

When V >= V_th: record spike, reset V = V_reset, enforce refractory period.

State: single variable V (membrane potential).
"""

import numpy as np
import numba as nb

# (V_rest, tau_m, R_m, V_th, V_reset, t_ref)
LIF_DEFAULT_PARAMS = (
    -65.0,   # V_rest  : resting potential (mV)
    10.0,    # tau_m   : membrane time constant (ms)
    10.0,    # R_m     : membrane resistance (MΩ)
    -50.0,   # V_th    : spike threshold (mV)
    -65.0,   # V_reset : reset potential (mV)
    2.0,     # t_ref   : absolute refractory period (ms)
)


@nb.njit
def _lif_simulate_core(I_ext_const, t_end, dt, V_rest, tau_m, R_m, V_th, V_reset, t_ref):
    n_steps = int(t_end / dt)
    t = np.empty(n_steps + 1)
    V = np.empty(n_steps + 1)
    spike_indices = np.empty(n_steps, dtype=nb.int64)
    n_spikes = 0

    t[0] = 0.0
    V[0] = V_rest
    t_last_spike = -1e9

    for i in range(n_steps):
        t[i + 1] = t[i] + dt

        if (t[i] - t_last_spike) < t_ref:
            V[i + 1] = V_reset
            continue

        dVdt = (-(V[i] - V_rest) + R_m * I_ext_const) / tau_m
        V[i + 1] = V[i] + dt * dVdt

        if V[i + 1] >= V_th:
            spike_indices[n_spikes] = i + 1
            n_spikes += 1
            V[i + 1] = V_reset
            t_last_spike = t[i + 1]

    return t, V, spike_indices[:n_spikes]


def lif_simulate(I_ext, t_end, dt, params=None):
    """Simulate a LIF neuron with constant injected current.

    Parameters
    ----------
    I_ext : float
        Constant injected current (nA).
    t_end : float
        Simulation duration (ms).
    dt : float
        Time step (ms).
    params : tuple, optional
        (V_rest, tau_m, R_m, V_th, V_reset, t_ref). Defaults to LIF_DEFAULT_PARAMS.

    Returns
    -------
    t : np.ndarray
    V : np.ndarray
    spike_indices : np.ndarray of int
    """
    if params is None:
        params = LIF_DEFAULT_PARAMS
    V_rest, tau_m, R_m, V_th, V_reset, t_ref = params
    return _lif_simulate_core(I_ext, t_end, dt, V_rest, tau_m, R_m, V_th, V_reset, t_ref)


@nb.njit
def lif_simulate_array(I_ext_array, dt, V_rest, tau_m, R_m, V_th, V_reset, t_ref):
    """Simulate LIF with a time-varying current array.

    Parameters
    ----------
    I_ext_array : np.ndarray
        Current at each time step.

    Returns
    -------
    V : np.ndarray
    spike_indices : np.ndarray
    """
    n_steps = len(I_ext_array)
    V = np.empty(n_steps)
    spike_indices = np.empty(n_steps, dtype=nb.int64)
    n_spikes = 0

    V[0] = V_rest
    t_last_spike = -1e9
    t = 0.0

    for i in range(1, n_steps):
        t += dt
        if (t - t_last_spike) < t_ref:
            V[i] = V_reset
            continue

        dVdt = (-(V[i - 1] - V_rest) + R_m * I_ext_array[i - 1]) / tau_m
        V[i] = V[i - 1] + dt * dVdt

        if V[i] >= V_th:
            spike_indices[n_spikes] = i
            n_spikes += 1
            V[i] = V_reset
            t_last_spike = t

    return V, spike_indices[:n_spikes]
