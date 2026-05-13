# Neural Membrane Dynamics
### From Biophysics to Dynamical Systems

A first-principles exploration of the Hodgkin-Huxley model, deriving the mathematics from membrane biophysics, implementing efficient Numba-accelerated simulations, analyzing dynamical stability through phase plane analysis and bifurcation theory, and comparing against the Leaky Integrate-and-Fire model.

---

## Overview

This project develops the full theory of neural excitability from the ground up. Starting with the physical chemistry of ion transport across lipid bilayers, we build the Hodgkin-Huxley equations step by step, implement numerically efficient solvers using Numba JIT compilation, and subject the resulting dynamical system to a rigorous mathematical analysis. The work is structured as a self-contained computational notebook that emphasizes both mathematical derivation and scientific computing best practices.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Neuron as an Electrical Circuit](#2-the-neuron-as-an-electrical-circuit)
3. [Deriving the Hodgkin-Huxley Equations](#3-deriving-the-hodgkin-huxley-equations)
4. [Numerical Implementation](#4-numerical-implementation)
5. [Simulating Neural Dynamics](#5-simulating-neural-dynamics)
6. [Phase Plane Analysis](#6-phase-plane-analysis)
7. [Dynamical Stability and Bifurcations](#7-dynamical-stability-and-bifurcations)
8. [The Leaky Integrate-and-Fire Model](#8-the-leaky-integrate-and-fire-model)
9. [Quantitative Comparison: HH vs LIF](#9-quantitative-comparison-hh-vs-lif)
10. [Conclusions and References](#10-conclusions-and-references)

---

## Key Features

- **Full mathematical derivations from first principles** — Nernst equation, Goldman-Hodgkin-Katz current equation, and the complete Hodgkin-Huxley formalism, derived from membrane biophysics rather than assumed
- **Numba JIT-compiled solvers** — Euler and Runge-Kutta 4 integrators compiled with `@njit` for high-performance simulation without leaving Python
- **Phase plane analysis** — Nullcline geometry, vector fields, and phase trajectories revealing the qualitative structure of neural excitability
- **Bifurcation analysis** — Hopf bifurcation detection, Type I vs Type II excitability classification, and F-I curve characterization
- **Quantitative HH vs LIF comparison** — Side-by-side benchmarks of accuracy, computational cost, and physiological fidelity

---

## Project Structure

```
neural-membrane-dynamics/
├── notebooks/
│   └── hodgkin_huxley_dynamics.ipynb
├── src/
│   └── neural_dynamics/
│       ├── __init__.py
│       ├── hodgkin_huxley.py
│       ├── integrate_and_fire.py
│       ├── solvers.py
│       └── analysis.py
├── tests/
├── figures/
├── requirements.txt
└── README.md
```

---

## Installation

```bash
git clone https://github.com/Mixnikon108/neural-membrane-dynamics.git
cd neural-membrane-dynamics
pip install -r requirements.txt
```

**Dependencies:** NumPy, Numba, Matplotlib, SciPy, Jupyter, pytest.

---

## Usage

Open the main notebook to follow the full derivation and analysis interactively:

```bash
jupyter notebook notebooks/hodgkin_huxley_dynamics.ipynb
```

---

## Running Tests

```bash
PYTHONPATH=src pytest tests/ -v
```

The test suite covers numerical solver accuracy, gating variable consistency, spike detection, and LIF threshold behaviour.

---

## References

Hodgkin, A. L. & Huxley, A. F. (1952). A quantitative description of membrane current and its application to conduction and excitation in nerve. *Journal of Physiology*, 117(4), 500-544.

Goldman, D. E. (1943). Potential, impedance, and rectification in membranes. *Journal of General Physiology*, 27(1), 37-60.

Rinzel, J. & Ermentrout, B. (1989). Analysis of neural excitability and oscillations. In *Methods in Neuronal Modeling*, MIT Press.

Izhikevich, E. M. (2007). *Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting*. MIT Press.

Dayan, P. & Abbott, L. F. (2001). *Theoretical Neuroscience: Computational and Mathematical Modeling of Neural Systems*. MIT Press.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
