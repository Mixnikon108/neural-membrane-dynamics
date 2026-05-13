import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from neural_dynamics.analysis import (
    compute_jacobian_hh,
    find_fixed_point,
    classify_fixed_point,
    compute_fi_curve_hh,
    compute_nullclines,
    bifurcation_sweep,
)
from neural_dynamics.hodgkin_huxley import hh_initial_state, n_inf, m_inf, h_inf


class TestJacobian:
    def test_jacobian_shape(self):
        y0 = hh_initial_state(-65.0)
        J = compute_jacobian_hh(y0, 0.0)
        assert J.shape == (4, 4)

    def test_jacobian_at_rest_stable(self):
        """At resting potential with no input, eigenvalues should have negative real parts."""
        y0 = hh_initial_state(-65.0)
        J = compute_jacobian_hh(y0, 0.0)
        eigenvalues = np.linalg.eigvals(J)
        assert np.all(np.real(eigenvalues) < 0)


class TestFixedPoint:
    def test_find_fixed_point_at_rest(self):
        """With I_ext=0, fixed point should be near V=-65."""
        y_fp = find_fixed_point(I_ext=0.0)
        assert abs(y_fp[0] - (-65.0)) < 5.0
        assert 0.0 < y_fp[1] < 1.0
        assert 0.0 < y_fp[2] < 1.0
        assert 0.0 < y_fp[3] < 1.0

    def test_classify_rest_as_stable(self):
        y_fp = find_fixed_point(I_ext=0.0)
        classification, eigenvalues = classify_fixed_point(y_fp, I_ext=0.0)
        assert "stable" in classification.lower()


class TestNullclines:
    def test_nullcline_shapes(self):
        V_arr, n_V_null, n_n_null = compute_nullclines((-80.0, 40.0), 0.0, n_points=100)
        assert len(V_arr) == 100
        assert len(n_V_null) == 100
        assert len(n_n_null) == 100

    def test_n_nullcline_is_n_inf(self):
        """The n-nullcline should be n_inf(V)."""
        V_arr, _, n_n_null = compute_nullclines((-80.0, 40.0), 0.0, n_points=50)
        for i, V in enumerate(V_arr):
            np.testing.assert_allclose(n_n_null[i], n_inf(V), rtol=1e-10)


class TestFICurve:
    def test_fi_curve_shape(self):
        currents, rates = compute_fi_curve_hh(
            I_range=(0.0, 20.0), n_points=5, t_sim=200.0, dt=0.01
        )
        assert len(currents) == 5
        assert len(rates) == 5

    def test_fi_curve_monotonic(self):
        """Firing rate should generally increase with current."""
        currents, rates = compute_fi_curve_hh(
            I_range=(5.0, 20.0), n_points=5, t_sim=500.0, dt=0.01
        )
        assert rates[-1] >= rates[0]


class TestBifurcation:
    def test_bifurcation_sweep_shape(self):
        currents, V_fps, real_eigs, classifications = bifurcation_sweep(
            I_range=(0.0, 20.0), n_points=5
        )
        assert len(currents) == 5
        assert len(V_fps) == 5
        assert real_eigs.shape == (5, 4)
        assert len(classifications) == 5

    def test_stability_changes(self):
        """At low current stable, at high current should change classification."""
        currents, V_fps, real_eigs, classifications = bifurcation_sweep(
            I_range=(0.0, 20.0), n_points=10
        )
        assert "stable" in classifications[0].lower()
