import numpy as np
import numba as nb
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from neural_dynamics.hodgkin_huxley import (
    alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h,
    n_inf, m_inf, h_inf,
    hh_rhs, hh_rhs_parameterized,
    HH_DEFAULT_PARAMS,
    hh_initial_state,
)
from neural_dynamics.solvers import solve_ode


class TestRateFunctions:
    def test_alpha_beta_positive(self):
        """Rate functions must be non-negative for all voltages."""
        voltages = np.linspace(-100.0, 50.0, 100)
        for V in voltages:
            assert alpha_n(V) >= 0
            assert beta_n(V) >= 0
            assert alpha_m(V) >= 0
            assert beta_m(V) >= 0
            assert alpha_h(V) >= 0
            assert beta_h(V) >= 0

    def test_steady_state_bounds(self):
        """Gating steady states must be in [0, 1]."""
        voltages = np.linspace(-100.0, 50.0, 100)
        for V in voltages:
            assert 0.0 <= n_inf(V) <= 1.0
            assert 0.0 <= m_inf(V) <= 1.0
            assert 0.0 <= h_inf(V) <= 1.0

    def test_m_inf_sigmoid_shape(self):
        """m_inf should be near 0 for very negative V and near 1 for positive V."""
        assert m_inf(-80.0) < 0.05
        assert m_inf(20.0) > 0.95

    def test_h_inf_sigmoid_shape(self):
        """h_inf should be near 1 for very negative V and near 0 for positive V."""
        assert h_inf(-80.0) > 0.90
        assert h_inf(20.0) < 0.05


class TestHHResting:
    def test_resting_potential_stability(self):
        """At resting potential with no input, dV/dt should be approximately zero."""
        V_rest = -65.0
        y0 = hh_initial_state(V_rest)
        dydt = hh_rhs(0.0, y0)
        assert abs(dydt[0]) < 0.5

    def test_resting_gating_stability(self):
        """At rest with steady-state gating, dn/dt, dm/dt, dh/dt should be ~0."""
        V_rest = -65.0
        y0 = hh_initial_state(V_rest)
        dydt = hh_rhs(0.0, y0)
        assert abs(dydt[1]) < 0.01
        assert abs(dydt[2]) < 0.01
        assert abs(dydt[3]) < 0.01


class TestHHSpike:
    def test_suprathreshold_generates_spike(self):
        """A strong current pulse should generate an action potential."""
        y0 = hh_initial_state(-65.0)

        @nb.njit
        def rhs_with_current(t, y):
            I_ext = 10.0 if t < 1.0 else 0.0
            return hh_rhs_parameterized(t, y, I_ext, 1.0, 120.0, 36.0, 0.3, 50.0, -77.0, -54.4)

        t, y = solve_ode(rhs_with_current, y0, t_span=(0.0, 5.0), dt=0.01, method="rk4")
        V = y[:, 0]
        assert np.max(V) > 0.0
        assert np.min(V) < -70.0

    def test_subthreshold_no_spike(self):
        """A weak current should not produce a spike."""
        y0 = hh_initial_state(-65.0)

        @nb.njit
        def rhs_weak(t, y):
            I_ext = 2.0 if t < 1.0 else 0.0
            return hh_rhs_parameterized(t, y, I_ext, 1.0, 120.0, 36.0, 0.3, 50.0, -77.0, -54.4)

        t, y = solve_ode(rhs_weak, y0, t_span=(0.0, 5.0), dt=0.01, method="rk4")
        V = y[:, 0]
        assert np.max(V) < -40.0


class TestHHParams:
    def test_default_params_values(self):
        """Verify original HH 1952 parameter values."""
        C_m, g_Na, g_K, g_L, E_Na, E_K, E_L = HH_DEFAULT_PARAMS
        assert C_m == 1.0
        assert g_Na == 120.0
        assert g_K == 36.0
        assert g_L == 0.3
        assert E_Na == 50.0
        assert E_K == -77.0
        assert abs(E_L - (-54.4)) < 0.1
