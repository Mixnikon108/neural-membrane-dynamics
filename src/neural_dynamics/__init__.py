"""Neural Membrane Dynamics: Hodgkin-Huxley and Integrate-and-Fire models."""

from .hodgkin_huxley import (
    hh_rhs,
    hh_rhs_parameterized,
    make_hh_rhs,
    hh_initial_state,
    hh_solve_current_array,
    hh_solve_constant,
    alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h,
    n_inf, m_inf, h_inf, tau_n, tau_m, tau_h,
    HH_DEFAULT_PARAMS,
)
from .integrate_and_fire import (
    lif_simulate,
    lif_simulate_array,
    LIF_DEFAULT_PARAMS,
)
from .solvers import euler_step, rk4_step, solve_ode
from .analysis import (
    compute_jacobian_hh,
    find_fixed_point,
    classify_fixed_point,
    compute_nullclines,
    compute_fi_curve_hh,
    bifurcation_sweep,
)
