import numpy as np
import numba as nb
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from neural_dynamics.solvers import euler_step, rk4_step, solve_ode


@nb.njit
def exponential_rhs(t, y):
    """dy/dt = -y, solution: y(t) = y0 * exp(-t)"""
    return -y


@nb.njit
def harmonic_rhs(t, y):
    """Simple harmonic oscillator: dy0/dt = y1, dy1/dt = -y0"""
    return np.array([y[1], -y[0]])


class TestEulerStep:
    def test_single_step_exponential(self):
        y0 = np.array([1.0])
        dt = 0.01
        y1 = euler_step(exponential_rhs, 0.0, y0, dt)
        expected = np.array([1.0 + dt * (-1.0)])
        np.testing.assert_allclose(y1, expected, rtol=1e-12)

    def test_preserves_shape(self):
        y0 = np.array([1.0, 0.0])
        y1 = euler_step(harmonic_rhs, 0.0, y0, 0.01)
        assert y1.shape == y0.shape


class TestRK4Step:
    def test_single_step_accuracy(self):
        """RK4 should be much more accurate than Euler for smooth ODEs."""
        y0 = np.array([1.0])
        dt = 0.1
        y_rk4 = rk4_step(exponential_rhs, 0.0, y0, dt)
        y_exact = np.array([np.exp(-0.1)])
        np.testing.assert_allclose(y_rk4, y_exact, rtol=1e-7)

    def test_fourth_order_convergence(self):
        """Error should decrease as O(h^4)."""
        y0 = np.array([1.0])
        errors = []
        dts = [0.1, 0.05, 0.025]
        for dt in dts:
            y_rk4 = rk4_step(exponential_rhs, 0.0, y0, dt)
            y_exact = np.array([np.exp(-dt)])
            errors.append(abs(y_rk4[0] - y_exact[0]))
        ratio1 = errors[0] / errors[1]
        ratio2 = errors[1] / errors[2]
        assert ratio1 > 14
        assert ratio2 > 14


class TestSolveODE:
    def test_exponential_decay(self):
        y0 = np.array([1.0])
        t, y = solve_ode(exponential_rhs, y0, t_span=(0.0, 1.0), dt=0.001, method="rk4")
        np.testing.assert_allclose(y[-1], np.exp(-1.0), rtol=1e-6)
        assert t.shape[0] == y.shape[0]

    def test_harmonic_oscillator_energy(self):
        """Energy (y0^2 + y1^2) should be conserved for harmonic oscillator."""
        y0 = np.array([1.0, 0.0])
        t, y = solve_ode(harmonic_rhs, y0, t_span=(0.0, 10.0), dt=0.001, method="rk4")
        energy_initial = y[0, 0]**2 + y[0, 1]**2
        energy_final = y[-1, 0]**2 + y[-1, 1]**2
        np.testing.assert_allclose(energy_final, energy_initial, rtol=1e-4)

    def test_euler_method_flag(self):
        y0 = np.array([1.0])
        t, y = solve_ode(exponential_rhs, y0, t_span=(0.0, 1.0), dt=0.001, method="euler")
        np.testing.assert_allclose(y[-1], np.exp(-1.0), rtol=1e-2)
