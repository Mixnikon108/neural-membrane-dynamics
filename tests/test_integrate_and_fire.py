import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from neural_dynamics.integrate_and_fire import lif_simulate, LIF_DEFAULT_PARAMS


class TestLIFBasic:
    def test_subthreshold_no_spike(self):
        """Weak current should not reach threshold."""
        t, V, spikes = lif_simulate(
            I_ext=1.0, t_end=100.0, dt=0.01
        )
        assert len(spikes) == 0
        assert np.all(V < LIF_DEFAULT_PARAMS[3])  # V_th

    def test_suprathreshold_fires(self):
        """Strong constant current should produce spikes."""
        t, V, spikes = lif_simulate(
            I_ext=15.0, t_end=100.0, dt=0.01
        )
        assert len(spikes) > 0

    def test_firing_rate_increases_with_current(self):
        """Higher current -> higher firing rate."""
        _, _, spikes_low = lif_simulate(I_ext=12.0, t_end=500.0, dt=0.01)
        _, _, spikes_high = lif_simulate(I_ext=20.0, t_end=500.0, dt=0.01)
        assert len(spikes_high) > len(spikes_low)

    def test_reset_after_spike(self):
        """After each spike, V should reset to V_reset."""
        t, V, spikes = lif_simulate(I_ext=15.0, t_end=100.0, dt=0.01)
        V_reset = LIF_DEFAULT_PARAMS[4]
        for s_idx in spikes:
            if s_idx + 1 < len(V):
                np.testing.assert_allclose(V[s_idx + 1], V_reset, atol=1.0)

    def test_membrane_decay_to_rest(self):
        """With no input, V should stay at V_rest."""
        t, V, _ = lif_simulate(
            I_ext=0.0, t_end=100.0, dt=0.01
        )
        V_rest = LIF_DEFAULT_PARAMS[0]
        np.testing.assert_allclose(V[-1], V_rest, atol=0.1)
