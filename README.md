# Neural Membrane Dynamics: From Biophysics to Dynamical Systems

*A first-principles exploration of the Hodgkin-Huxley model*

## Why Model Neurons Mathematically?

The nervous system is arguably the most complex structure in the known universe. A single human brain contains roughly $10^{11}$ neurons, each forming thousands of synaptic connections. Yet the fundamental electrical event underlying all neural computation (the **action potential**) can be described by a system of just four ordinary differential equations. This remarkable compression of biological complexity into mathematical form is one of the great triumphs of 20th-century biophysics.

Mathematical models of neurons serve multiple purposes. First, they provide a **quantitative framework** for understanding how ionic currents flowing through protein channels in the cell membrane give rise to the rapid, all-or-nothing voltage spikes that carry information along nerve fibers. Without equations, we are limited to qualitative descriptions; with them, we can predict precisely how a neuron will respond to any given stimulus. Second, neural models reveal **universal dynamical principles** (bifurcations, limit cycles, excitability thresholds) that transcend the specific biological details and connect neuroscience to the broader theory of nonlinear dynamical systems.

Third, and perhaps most importantly, mathematical neuron models are the **building blocks of computational neuroscience**. From small circuits of a few neurons to large-scale brain simulations involving millions of units, every computational model of neural activity rests on some mathematical description of individual cells. Understanding the Hodgkin-Huxley model is therefore not merely an exercise in biophysics: it is the foundation upon which modern theoretical neuroscience is built.

### Historical Context

In the late 1940s and early 1950s, Alan Hodgkin and Andrew Huxley carried out a series of extraordinary experiments on the **giant axon of the Atlantic squid** (*Loligo*). The squid giant axon, with a diameter of up to 1 mm, was large enough to allow insertion of electrodes directly into the interior of the nerve fiber; a feat impossible with the much thinner axons of vertebrate neurons at the time.

Using the **voltage clamp technique**, Hodgkin and Huxley were able to hold the membrane potential at a fixed value and measure the ionic currents flowing through the membrane at that voltage. By systematically varying the clamp voltage and using ion substitution experiments (replacing external sodium with choline), they separated the total membrane current into its sodium, potassium, and leak components. They then fitted empirical equations to describe how the conductance of each ion channel depended on voltage and time.

The result, published in a landmark series of five papers in 1952, was a mathematical model that could **quantitatively reproduce** the shape, amplitude, duration, and propagation velocity of the action potential: all from first principles, without any free parameters adjusted after the fact. For this work, Hodgkin and Huxley were awarded the **Nobel Prize in Physiology or Medicine in 1963**.

### Roadmap

In this notebook, we will:

1. **Derive the electrical circuit model** of the neuronal membrane from physical principles
2. **Develop the Hodgkin-Huxley equations** step by step, motivating every term
3. **Simulate action potentials** and explore how they depend on stimulus parameters
4. **Analyze the ionic currents** underlying each phase of the action potential
5. **Investigate the dynamical systems structure**: fixed points, stability, and bifurcations
6. **Construct phase portraits** and nullcline diagrams
7. **Compute the frequency-current (f-I) relationship** and explore neural coding
8. **Compare with simplified models** such as the leaky integrate-and-fire neuron

---

# 2. The Neuron as an Electrical Circuit

To understand how a neuron generates electrical signals, we must first understand the **physical properties of the cell membrane** and how they give rise to electrical behavior. In this section, we build up the equivalent circuit model of a small patch of neuronal membrane, deriving each component from first principles.

## 2.1 The Cell Membrane as a Capacitor

The neuronal cell membrane is a **lipid bilayer**: a thin sheet (approximately $d \approx 7\text{-}8 \; \text{nm}$ thick) composed of phospholipid molecules whose hydrophobic tails face inward, forming a nearly impermeable barrier to ions. On either side of this insulating layer lies an electrically conducting aqueous solution: the **intracellular fluid** (cytoplasm) and the **extracellular fluid**.

This arrangement (two conductors separated by an insulating dielectric) is precisely the definition of a **parallel-plate capacitor**. From elementary electrostatics, the capacitance of a parallel-plate capacitor with plate area $A$, plate separation $d$, and dielectric constant $\varepsilon$ is:

$$C = \varepsilon \frac{A}{d} = \varepsilon_0 \varepsilon_r \frac{A}{d}$$

where $\varepsilon_0 \approx 8.85 \times 10^{-12} \; \text{F/m}$ is the permittivity of free space and $\varepsilon_r \approx 2$ is the relative permittivity of the lipid core. Substituting typical values:

$$c_m = \frac{C}{A} = \frac{\varepsilon_0 \varepsilon_r}{d} \approx \frac{(8.85 \times 10^{-12})(2)}{7.5 \times 10^{-9}} \approx 2.4 \times 10^{-3} \; \text{F/m}^2 \approx 0.24 \; \mu\text{F/cm}^2$$

Experimentally, the specific membrane capacitance is measured to be:

$$\boxed{C_m \approx 1 \; \mu\text{F/cm}^2}$$

The discrepancy (roughly a factor of 4) arises because the effective dielectric constant of the full membrane (including polar head groups, bound water molecules, and membrane proteins) is higher than that of the pure lipid core. The value $C_m \approx 1 \; \mu\text{F/cm}^2$ is remarkably universal across cell types and organisms.

The charge stored on the membrane capacitor is $Q = C_m V_m$, where $V_m = V_\text{in} - V_\text{out}$ is the **membrane potential**. The capacitive current (current flowing to charge or discharge the capacitor) is therefore:

$$\boxed{I_C = C_m \frac{dV_m}{dt}}$$

This tells us that the membrane potential can change only when current flows onto or off of the membrane capacitor. This is the starting point for all electrical models of neurons.

## 2.2 Ion Channels as Conductances

While the lipid bilayer is an excellent insulator, the cell membrane is studded with **ion channel proteins**: large transmembrane proteins that form aqueous pores through which specific ions can flow. The key channel types for the action potential are:

- **Voltage-gated $	ext{Na}^+$ channels**: selectively permeable to sodium ions, open rapidly upon depolarization
- **Voltage-gated $	ext{K}^+$ channels**: selectively permeable to potassium ions, open more slowly upon depolarization
- **Leak channels**: always open, primarily permeable to $	ext{K}^+$ and $	ext{Cl}^-$, responsible for the resting potential

Each population of ion channels can be characterized by a **conductance** $g_\text{ion}$ (the inverse of resistance), measured in units of $\text{mS/cm}^2$ (millisiemens per square centimeter of membrane). The conductance represents how easily ions can flow through the channels: high conductance means many channels are open, low conductance means few are open.

The current through a given ion channel population obeys an **Ohmic** (linear) relationship:

$$\boxed{I_\text{ion} = g_\text{ion} \cdot (V_m - E_\text{ion})}$$

where $E_\text{ion}$ is the **reversal potential** (or equilibrium potential) of that ion species. The driving force $(V_m - E_\text{ion})$ is the difference between the actual membrane potential and the potential at which the net current through those channels would be zero. When $V_m = E_\text{ion}$, the electrical and chemical driving forces on the ion are exactly balanced, and no net current flows: hence the name "reversal potential."

The critical insight that distinguishes the Hodgkin-Huxley model from simpler models is that the conductances $g_\text{Na}$ and $g_\text{K}$ are **not constant**: they are functions of both voltage and time. Understanding how to describe this voltage- and time-dependence is the central challenge of the model, and we will address it in Section 3.

## 2.3 The Nernst Equation (Full Derivation)

Where does the reversal potential $E_\text{ion}$ come from? It arises from the fact that ion concentrations differ dramatically between the inside and outside of the cell. For example, $	ext{K}^+$ is roughly 20 times more concentrated inside the cell than outside, while $	ext{Na}^+$ is roughly 10 times more concentrated outside. These concentration gradients are maintained by metabolic pumps (primarily the $	ext{Na}^+$/$	ext{K}^+$-ATPase).

When a channel selective for a particular ion opens, two forces act on that ion:

1. **The chemical force** (diffusion): ions tend to flow down their concentration gradient, from high to low concentration.
2. **The electrical force** (drift): ions are charged particles and are attracted or repelled by the electric field across the membrane.

At equilibrium, these two forces balance exactly. Let us derive the voltage at which this occurs.

### Derivation from Thermodynamic First Principles

Consider a single ion species with valence $z$ (e.g., $z = +1$ for $	ext{Na}^+$ and $	ext{K}^+$, $z = -1$ for $	ext{Cl}^-$). The **electrochemical potential** of this ion on each side of the membrane has two contributions:

$$\mu = \mu_\text{chemical} + \mu_\text{electrical}$$

The chemical potential of an ideal solute at concentration $[\text{ion}]$ is:

$$\mu_\text{chemical} = \mu^0 + RT \ln[\text{ion}]$$

where $R = 8.314 \; \text{J/(mol} \cdot \text{K)}$ is the gas constant, $T$ is absolute temperature, and $\mu^0$ is a reference chemical potential.

The electrical potential energy of one mole of ions with charge $z$ at voltage $V$ is:

$$\mu_\text{electrical} = zFV$$

where $F = 96{,}485 \; \text{C/mol}$ is Faraday's constant. Therefore the total electrochemical potential is:

$$\mu = \mu^0 + RT \ln[\text{ion}] + zFV$$

At **thermodynamic equilibrium**, the electrochemical potential must be equal on both sides of the membrane:

$$\mu_\text{out} = \mu_\text{in}$$

$$\mu^0 + RT \ln[\text{ion}]_\text{out} + zF V_\text{out} = \mu^0 + RT \ln[\text{ion}]_\text{in} + zF V_\text{in}$$

The reference potentials $\mu^0$ cancel. Rearranging to solve for the transmembrane potential $E_\text{ion} = V_\text{in} - V_\text{out}$:

$$zF(V_\text{in} - V_\text{out}) = RT \ln[\text{ion}]_\text{out} - RT \ln[\text{ion}]_\text{in}$$

$$E_\text{ion} = V_\text{in} - V_\text{out} = \frac{RT}{zF} \ln \frac{[\text{ion}]_\text{out}}{[\text{ion}]_\text{in}}$$

This is the **Nernst equation**:

$$\boxed{E_\text{ion} = \frac{RT}{zF} \ln \frac{[\text{ion}]_\text{out}}{[\text{ion}]_\text{in}}}$$

At body temperature ($T = 310 \; \text{K}$) or squid axon temperature ($T \approx 279 \; \text{K}$ for $6.3°\text{C}$), the prefactor $RT/F \approx 24\text{-}27 \; \text{mV}$ for monovalent ions. Using the ionic concentrations of the squid giant axon:

| Ion | $[\text{ion}]_\text{in}$ (mM) | $[\text{ion}]_\text{out}$ (mM) | $z$ | $E_\text{ion}$ (mV) |
|-----|------|------|---|------|
| $	ext{K}^+$ | 400 | 20 | +1 | $\approx -77$ |
| $	ext{Na}^+$ | 50 | 440 | +1 | $\approx +50$ |
| $	ext{Cl}^-$ | 40 | 560 | $-1$ | $\approx -66$ |

The large difference between $E_\text{Na} \approx +50 \; \text{mV}$ and $E_\text{K} \approx -77 \; \text{mV}$ is what makes the action potential possible: the membrane can swing between these two "battery" voltages by selectively opening $	ext{Na}^+$ or $	ext{K}^+$ channels.

<p align="center">
  <img src="figures/fig01_nernst_equation.png" width="700" alt="Nernst equilibrium potential as a function of concentration ratio for different ion valences.">
  <br><em>Figure 1. Nernst equilibrium potential as a function of concentration ratio for different ion valences.</em>
</p>

## 2.4 The Goldman-Hodgkin-Katz Equation

The Nernst equation gives the equilibrium potential for a single ion species. But at rest, the membrane is permeable to **multiple ion species simultaneously**: primarily $	ext{K}^+$, $	ext{Na}^+$, and $	ext{Cl}^-$. The resting membrane potential is therefore not equal to any single Nernst potential, but rather a weighted average that depends on the **relative permeabilities** of the membrane to each ion.

The **Goldman-Hodgkin-Katz (GHK) voltage equation** generalizes the Nernst equation to multiple ions. It is derived from the Nernst-Planck equation for electrodiffusion under the assumption of a constant electric field across the membrane (the "constant field" approximation). For $	ext{Na}^+$, $	ext{K}^+$, and $	ext{Cl}^-$:

$$\boxed{V_m = \frac{RT}{F} \ln \frac{P_\text{K}[\text{K}^+]_\text{out} + P_\text{Na}[\text{Na}^+]_\text{out} + P_\text{Cl}[\text{Cl}^-]_\text{in}}{P_\text{K}[\text{K}^+]_\text{in} + P_\text{Na}[\text{Na}^+]_\text{in} + P_\text{Cl}[\text{Cl}^-]_\text{out}}}$$

where $P_\text{K}$, $P_\text{Na}$, $P_\text{Cl}$ are the membrane permeabilities to each ion. Note that chloride concentrations are "flipped" (inside in numerator, outside in denominator) because $	ext{Cl}^-$ has negative valence.

At rest, the squid axon membrane is much more permeable to $	ext{K}^+$ than to $	ext{Na}^+$, with typical permeability ratios $P_\text{K} : P_\text{Na} : P_\text{Cl} \approx 1 : 0.04 : 0.45$. Substituting the squid axon concentrations:

$$V_\text{rest} \approx -65 \; \text{mV}$$

This value is close to $E_\text{K}$ (since $	ext{K}^+$ permeability dominates at rest) but is pulled slightly positive by the small $	ext{Na}^+$ permeability. During an action potential, when $	ext{Na}^+$ channels open and $P_\text{Na}$ increases dramatically, the GHK equation predicts that $V_m$ will swing toward $E_\text{Na} \approx +50 \; \text{mV}$.

## 2.5 The Equivalent Circuit

We can now assemble the complete electrical equivalent circuit for a patch of neuronal membrane. The membrane consists of:

1. A **capacitor** $C_m$ (the lipid bilayer)
2. Three parallel **conductance branches**, one for each ion species, each consisting of a variable resistor (the channels) in series with a battery (the Nernst potential)
3. An optional **external current source** $I_\text{ext}$ (representing synaptic input or an experimentalist's electrode)

By **Kirchhoff's current law**, the total current flowing into the membrane patch must be zero. All current that enters as $I_\text{ext}$ must either charge the capacitor or flow out through the ion channels:

$$I_\text{ext} = I_C + I_\text{Na} + I_\text{K} + I_L$$

Substituting the expressions for each current:

$$I_\text{ext} = C_m \frac{dV_m}{dt} + g_\text{Na}(V_m - E_\text{Na}) + g_\text{K}(V_m - E_\text{K}) + g_L(V_m - E_L)$$

Rearranging to express the dynamics of the membrane potential:

$$\boxed{C_m \frac{dV_m}{dt} = -g_\text{Na}(V_m - E_\text{Na}) - g_\text{K}(V_m - E_\text{K}) - g_L(V_m - E_L) + I_\text{ext}}$$

**This is THE fundamental equation of the Hodgkin-Huxley model.** Everything that follows (the gating variables, the rate functions, the full system of ODEs) is devoted to specifying how the conductances $g_\text{Na}$ and $g_\text{K}$ depend on voltage and time. The leak conductance $g_L$ is constant (representing channels that are always open).

The physical content of this equation is transparent:
- The left side is the rate of change of charge on the membrane capacitor.
- The right side is the net current: external current minus the sum of all ionic currents.
- Each ionic current is proportional to its conductance times its driving force.
- The membrane potential $V_m$ evolves toward whichever reversal potential has the largest conductance at that moment.

With constant conductances, this equation describes a simple leaky integrator; the membrane potential exponentially relaxes toward a weighted average of the reversal potentials. The action potential emerges only when we allow $g_\text{Na}$ and $g_\text{K}$ to depend on voltage, creating a powerful positive feedback loop. This is the subject of the next section.

---

# 3. Deriving the Hodgkin-Huxley Equations

In Section 2 we established that the membrane potential obeys $C_m \, dV/dt = -\sum I_\text{ion} + I_\text{ext}$, with ionic currents of the form $I_\text{ion} = g_\text{ion}(V_m - E_\text{ion})$. The remaining challenge is to describe how the conductances $g_\text{Na}$ and $g_\text{K}$ depend on voltage and time. This is where Hodgkin and Huxley's experimental genius and mathematical creativity come together.

## 3.1 Voltage-Dependent Conductances

Hodgkin and Huxley's voltage clamp experiments revealed a crucial fact: **the membrane conductances are not constant.** When the membrane is depolarized (voltage stepped to a more positive value):

- The **$	ext{Na}^+$ conductance** rises rapidly to a peak and then decays back toward zero, even though the voltage is held constant. This transient behavior reflects two distinct processes: fast **activation** (opening) followed by slower **inactivation** (a conformational change that blocks the open channel).
- The **$	ext{K}^+$ conductance** rises more slowly to a sustained plateau and remains elevated as long as the depolarization is maintained. $	ext{K}^+$ channels activate but do not inactivate on the timescale of the action potential.

These are the experimental signatures that the model must capture. The conductances are functions of both voltage (because the rate of opening/closing depends on $V_m$) and time (because the channels take time to respond to voltage changes).

## 3.2 Gating Variables

Hodgkin and Huxley's key insight was to decompose each conductance into the product of a **maximum conductance** $\bar{g}$ (when all channels are open) and one or more **gating variables** that represent the fraction of channels in the open state. Each gating variable takes values between 0 (all gates closed) and 1 (all gates open).

### The $	ext{K}^+$ Channel: Four Identical Gates ($n^4$)

The $	ext{K}^+$ conductance rises sigmoidally upon depolarization, not exponentially. Hodgkin and Huxley found that the time course of $	ext{K}^+$ activation could be well described by the **fourth power** of a single gating variable $n$:

$$\boxed{g_\text{K}(t) = \bar{g}_\text{K} \cdot n(t)^4}$$

Why $n^4$? The physical picture is that the $	ext{K}^+$ channel has **four identical, independent subunits** (gates), each of which must be in the "open" configuration for the channel to conduct. If each gate has probability $n$ of being open (independently), then the probability that all four are open simultaneously is $n^4$. This explains the sigmoidal rise: at early times, $n$ is increasing from a small value, and $n^4$ increases even more steeply. We now know from molecular biology that voltage-gated $	ext{K}^+$ channels are indeed tetramers (four identical $\alpha$-subunits arranged around a central pore) a remarkable confirmation of HH's purely mathematical deduction.

### The $	ext{Na}^+$ Channel: Three Activation Gates and One Inactivation Gate ($m^3 h$)

The $	ext{Na}^+$ conductance has a more complex time course: it rises rapidly (activation) and then decays (inactivation), even at constant voltage. Hodgkin and Huxley modeled this by introducing **two types of gates**:

- **$m$**: an activation gate. Three identical $m$-gates must be open for the channel to conduct (similar logic to $n$ for $	ext{K}^+$).
- **$h$**: an inactivation gate. This single gate must also be open. Unlike $m$, the $h$-gate is open at rest and *closes* upon depolarization.

The $	ext{Na}^+$ conductance is therefore:

$$\boxed{g_\text{Na}(t) = \bar{g}_\text{Na} \cdot m(t)^3 \cdot h(t)}$$

Why $m^3 h$? The $m^3$ factor produces a rapid, sigmoidal rise in conductance (since $m$ activates quickly upon depolarization). The $h$ factor produces a slower decay (since $h$ decreases upon depolarization, eventually shutting off the channel). The product $m^3 h$ thus naturally produces a transient conductance that rises fast and falls slowly: exactly as observed experimentally.

Molecularly, we now know that $	ext{Na}^+$ channels are single large proteins with four homologous domains (I-IV). Three of these contribute to activation (analogous to the three $m$-gates), while a cytoplasmic loop between domains III and IV acts as the "inactivation ball" that plugs the open channel (analogous to the $h$-gate).

## 3.3 First-Order Kinetics

Each gating variable $x \in \{n, m, h\}$ represents the probability that a single gate is in the open state. A gate can transition between two states:

$$\text{Closed} \underset{\beta_x(V)}{\overset{\alpha_x(V)}{\rightleftharpoons}} \text{Open}$$

where $\alpha_x(V)$ is the voltage-dependent **opening rate** (transition rate from closed to open, in $\text{ms}^{-1}$) and $\beta_x(V)$ is the **closing rate** (open to closed). If $x$ is the fraction of gates in the open state (and $1-x$ is the fraction closed), the rate of change is:

$$\boxed{\frac{dx}{dt} = \alpha_x(V)(1 - x) - \beta_x(V) \, x}$$

This is a **first-order linear ODE** in $x$ (for fixed $V$). The first term describes gates opening (proportional to the fraction that are closed, $1-x$), and the second describes gates closing (proportional to the fraction that are open, $x$).

### Equivalent Formulation: Steady-State and Time Constant

The rate equation can be rewritten in a more physically transparent form. Defining:

$$x_\infty(V) = \frac{\alpha_x(V)}{\alpha_x(V) + \beta_x(V)} \qquad \text{and} \qquad \tau_x(V) = \frac{1}{\alpha_x(V) + \beta_x(V)}$$

we can verify by substitution that the ODE becomes:

$$\boxed{\frac{dx}{dt} = \frac{x_\infty(V) - x}{\tau_x(V)}}$$

This form makes the physical meaning completely clear:

- $x_\infty(V)$ is the **steady-state value** of $x$ at voltage $V$. If the voltage is held constant, $x$ will exponentially relax toward $x_\infty(V)$.
- $\tau_x(V)$ is the **time constant** of this relaxation. A small $\tau_x$ means the gate responds quickly to voltage changes; a large $\tau_x$ means it responds slowly.

The solution at constant voltage is:

$$x(t) = x_\infty - (x_\infty - x_0) \, e^{-t/\tau_x}$$

where $x_0 = x(0)$ is the initial value. This exponential relaxation is the signature behavior of HH gating variables.

## 3.4 The Rate Functions

Hodgkin and Huxley determined the rate functions $\alpha_x(V)$ and $\beta_x(V)$ by fitting empirical expressions to their voltage clamp data. The original paper used a shifted voltage variable $v = V - V_\text{rest}$, but we present the equations in the **modern convention** where $V$ is the absolute membrane potential (rest $\approx -65 \; \text{mV}$):

### $	ext{K}^+$ activation ($n$)

$$\alpha_n(V) = \frac{0.01 \, (V + 55)}{1 - \exp\!\big(-(V+55)/10\big)} \qquad \beta_n(V) = 0.125 \, \exp\!\big(-(V+65)/80\big)$$

### $	ext{Na}^+$ activation ($m$)

$$\alpha_m(V) = \frac{0.1 \, (V + 40)}{1 - \exp\!\big(-(V+40)/10\big)} \qquad \beta_m(V) = 4 \, \exp\!\big(-(V+65)/18\big)$$

### $	ext{Na}^+$ inactivation ($h$)

$$\alpha_h(V) = 0.07 \, \exp\!\big(-(V+65)/20\big) \qquad \beta_h(V) = \frac{1}{1 + \exp\!\big(-(V+35)/10\big)}$$

### Singularity Handling

Note that $\alpha_n$ and $\alpha_m$ have the form $f(V) = a(V - V_0) / (1 - e^{-(V-V_0)/k})$, which is indeterminate ($0/0$) at $V = V_0$ (i.e., at $V = -55 \; \text{mV}$ for $\alpha_n$ and $V = -40 \; \text{mV}$ for $\alpha_m$). Applying **L'Hopital's rule**:

$$\lim_{V \to V_0} \frac{a(V - V_0)}{1 - e^{-(V-V_0)/k}} = \lim_{V \to V_0} \frac{a}{e^{-(V-V_0)/k} / k} = a \cdot k$$

So $\alpha_n(-55) = 0.01 \times 10 = 0.1$ and $\alpha_m(-40) = 0.1 \times 10 = 1.0$. Our numerical implementation handles this by checking whether $|V - V_0| < \epsilon$ and returning the limiting value directly.

<p align="center">
  <img src="figures/fig02_rate_functions.png" width="700" alt="All six Hodgkin-Huxley rate functions plotted against membrane voltage.">
  <br><em>Figure 2. All six Hodgkin-Huxley rate functions plotted against membrane voltage.</em>
</p>

<p align="center">
  <img src="figures/fig03_steady_state_gating.png" width="700" alt="Steady-state gating variables m_inf, h_inf, n_inf as functions of voltage.">
  <br><em>Figure 3. Steady-state gating variables m_inf, h_inf, n_inf as functions of voltage.</em>
</p>

<p align="center">
  <img src="figures/fig04_time_constants.png" width="700" alt="Voltage-dependent time constants for the three gating variables.">
  <br><em>Figure 4. Voltage-dependent time constants for the three gating variables.</em>
</p>

## 3.5 The Complete Model

We can now assemble the full Hodgkin-Huxley model. The state of the membrane patch at any instant is described by four variables: the membrane potential $V$ and the three gating variables $n$, $m$, $h$. Their dynamics are governed by the following system of **four coupled ordinary differential equations**:

$$\boxed{\begin{aligned}
C_m \frac{dV}{dt} &= -\bar{g}_\text{Na} \, m^3 h \, (V - E_\text{Na}) - \bar{g}_\text{K} \, n^4 \, (V - E_\text{K}) - g_L \, (V - E_L) + I_\text{ext} \\[6pt]
\frac{dn}{dt} &= \alpha_n(V)(1 - n) - \beta_n(V) \, n \\[4pt]
\frac{dm}{dt} &= \alpha_m(V)(1 - m) - \beta_m(V) \, m \\[4pt]
\frac{dh}{dt} &= \alpha_h(V)(1 - h) - \beta_h(V) \, h
\end{aligned}}$$

The coupling is bidirectional and nonlinear: the gating variables affect $V$ through the conductances $m^3 h$ and $n^4$, while $V$ affects the gating variables through the rate functions $\alpha_x(V)$ and $\beta_x(V)$. This feedback is what produces the rich dynamics of the model, including the action potential.

### Parameter Table

The following are the **original parameters** from Hodgkin and Huxley (1952), measured at $6.3\,°\text{C}$:

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Membrane capacitance | $C_m$ | 1.0 | $\mu\text{F/cm}^2$ |
| Max $	ext{Na}^+$ conductance | $\bar{g}_\text{Na}$ | 120.0 | $\text{mS/cm}^2$ |
| Max $	ext{K}^+$ conductance | $\bar{g}_\text{K}$ | 36.0 | $\text{mS/cm}^2$ |
| Leak conductance | $g_L$ | 0.3 | $\text{mS/cm}^2$ |
| $	ext{Na}^+$ reversal potential | $E_\text{Na}$ | +50.0 | mV |
| $	ext{K}^+$ reversal potential | $E_\text{K}$ | $-77.0$ | mV |
| Leak reversal potential | $E_L$ | $-54.4$ | mV |

Note the striking asymmetry: $\bar{g}_\text{Na}$ is more than three times larger than $\bar{g}_\text{K}$, which in turn is 120 times larger than $g_L$. When the $	ext{Na}^+$ channels open fully, they can drive an enormous inward current that rapidly depolarizes the membrane.

These seven parameters, together with the six rate functions defined in Section 3.4, completely specify the model. In the next sections, we will simulate this system numerically and explore its remarkably rich dynamical behavior.

---

# 4. Numerical Implementation

We have derived a beautiful system of four coupled nonlinear ODEs, the Hodgkin-Huxley equations. But beauty alone does not produce action potentials. To actually *solve* these equations, we must turn to **numerical methods**.

Why can’t we simply write down an analytical solution? The answer lies in the **nonlinear coupling** between the equations. The membrane potential $V$ appears inside exponential functions (through the rate functions $\alpha_x(V)$ and $\beta_x(V)$), and these nonlinear functions of $V$ are multiplied by the gating variables, which themselves depend on $V$. There is no closed-form solution to this system. Even the existence and uniqueness of solutions must be established through the Picard-Lindelöf theorem, which guarantees local existence but tells us nothing about the solution’s form.

We therefore need **numerical integration**: we approximate the continuous ODEs by stepping forward in discrete time increments $\Delta t$, computing the state at each step from the state at the previous step. The art of numerical methods lies in choosing schemes that are accurate (small error per step), stable (errors don’t grow catastrophically), and efficient (few function evaluations per step).

A particular challenge for the HH system is **stiffness**. The fast $	ext{Na}^+$ activation variable $m$ has a time constant as small as $\tau_m \approx 0.05 \; \text{ms}$ near the spike peak, while the slow $	ext{K}^+$ activation $n$ has $\tau_n \approx 5 \; \text{ms}$; a ratio of 100:1. This disparity in timescales means that explicit methods must use time steps small enough to resolve the fastest dynamics, even during periods when only the slow dynamics are active. We will see how this affects the choice of method and step size.

## 4.1 The Euler Method

The simplest numerical method for ODEs is the **forward Euler method**, derived directly from the Taylor expansion of the solution. Suppose we know the state $\mathbf{y}(t)$ and wish to approximate $\mathbf{y}(t + \Delta t)$. Expanding in a Taylor series:

$$\mathbf{y}(t + \Delta t) = \mathbf{y}(t) + \Delta t \, \mathbf{y}'(t) + \frac{\Delta t^2}{2} \, \mathbf{y}''(t) + \mathcal{O}(\Delta t^3)$$

Since the ODE tells us that $\mathbf{y}'(t) = \mathbf{f}(t, \mathbf{y}(t))$, we can drop all terms of order $\Delta t^2$ and higher to obtain the **Euler approximation**:

$$\boxed{\mathbf{y}_{n+1} = \mathbf{y}_n + \Delta t \, \mathbf{f}(t_n, \mathbf{y}_n)}$$

### Error Analysis

The **local truncation error** (error per step, assuming exact input) is the first neglected term in the Taylor series:

$$e_\text{local} = \frac{\Delta t^2}{2} \, \mathbf{y}''(\xi) = \mathcal{O}(\Delta t^2)$$

for some $\xi \in [t_n, t_{n+1}]$ (by the mean value theorem). Over a total integration time $T$, we take $N = T / \Delta t$ steps, so the **global error** accumulates as:

$$e_\text{global} \leq N \cdot e_\text{local} = \frac{T}{\Delta t} \cdot \mathcal{O}(\Delta t^2) = \mathcal{O}(\Delta t)$$

The Euler method is therefore a **first-order method**: halving the step size halves the global error (but doubles the computational cost).

### Stability Considerations

For the HH system, the Euler method requires $\Delta t \lesssim 0.01 \; \text{ms}$ to maintain stability, primarily because of the fast $	ext{Na}^+$ activation dynamics. With a time constant $\tau_m \approx 0.05 \; \text{ms}$, the stability condition $\Delta t < 2\tau_m$ gives $\Delta t < 0.1 \; \text{ms}$, and in practice a safety margin is needed. This is acceptable for short simulations but becomes expensive for long runs.

## 4.2 The 4th-Order Runge-Kutta Method (RK4)

Can we achieve higher accuracy without computing higher derivatives of $\mathbf{f}$? The **Runge-Kutta family** of methods answers this question affirmatively: instead of evaluating $\mathbf{f}$ at a single point (as in Euler), we evaluate it at **multiple carefully chosen points** within the interval $[t_n, t_n + \Delta t]$ and combine the results.

The most widely used member of this family is the **classical 4th-order Runge-Kutta method** (RK4). It computes four "slope estimates":

$$\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n, \, \mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{\Delta t}{2}, \, \mathbf{y}_n + \tfrac{\Delta t}{2} \, \mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + \tfrac{\Delta t}{2}, \, \mathbf{y}_n + \tfrac{\Delta t}{2} \, \mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{f}\!\left(t_n + \Delta t, \, \mathbf{y}_n + \Delta t \, \mathbf{k}_3\right)
\end{aligned}$$

and then advances the solution using their weighted average:

$$\boxed{\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)}$$

### Understanding Each Slope

The intuition behind the four slopes is as follows:

- $\mathbf{k}_1$: the slope at the **beginning** of the interval (same as Euler)
- $\mathbf{k}_2$: the slope at the **midpoint**, using an Euler half-step with slope $\mathbf{k}_1$ to predict the midpoint state. This is a first correction: we use the initial slope to estimate where we’ll be at the midpoint, then evaluate the slope there.
- $\mathbf{k}_3$: the slope at the **midpoint** again, but now using $\mathbf{k}_2$ (the corrected midpoint slope) for the half-step prediction. This is a second correction: we refine our midpoint estimate.
- $\mathbf{k}_4$: the slope at the **endpoint**, using a full step with slope $\mathbf{k}_3$.

### The 1:2:2:1 Weighting and Connection to Simpson’s Rule

The weights $\frac{1}{6}(1, 2, 2, 1)$ are not arbitrary. They correspond to **Simpson’s rule** for numerical integration. Recall that integrating $\mathbf{y}' = \mathbf{f}$ from $t_n$ to $t_n + \Delta t$ gives:

$$\mathbf{y}(t_n + \Delta t) - \mathbf{y}(t_n) = \int_{t_n}^{t_n + \Delta t} \mathbf{f}(\tau, \mathbf{y}(\tau)) \, d\tau$$

Simpson’s rule approximates an integral using values at the endpoints and midpoint with weights $\frac{\Delta t}{6}(1, 4, 1)$. The RK4 weights $\frac{\Delta t}{6}(1, 2, 2, 1)$ arise because the midpoint is sampled *twice* (as $\mathbf{k}_2$ and $\mathbf{k}_3$), splitting the weight of 4 into $2 + 2$.

### Error Analysis

The local truncation error of RK4 is $\mathcal{O}(\Delta t^5)$, and the global error is:

$$e_\text{global} = \mathcal{O}(\Delta t^4)$$

This is a dramatic improvement over Euler. Halving the step size reduces the error by a factor of $2^4 = 16$, at a cost of only doubling the number of steps (each step requires 4 function evaluations instead of 1, so the total cost increases by a factor of 8, but the error drops by 16). For the HH system, RK4 with $\Delta t = 0.01 \; \text{ms}$ gives results accurate to $\sim 10^{-8}$, far more than sufficient for any practical purpose.

## 4.3 Numba JIT Compilation

Python’s flexibility comes at a cost: as an interpreted language, every arithmetic operation incurs significant overhead from type checking, dynamic dispatch, and memory management. For a tight numerical loop (stepping an ODE forward in time for tens of thousands of iterations) this overhead can slow execution by a factor of 100 or more compared to compiled languages like C or Fortran.

**Numba** solves this problem by providing a **just-in-time (JIT) compiler** that translates Python functions directly into optimized machine code using the LLVM compiler framework. By adding the `@nb.njit` decorator (“nopython JIT”), we tell Numba to compile the function without any Python interpreter involvement: all operations are translated to native machine instructions, types are inferred at compile time, and the result runs at speeds comparable to hand-written C code.

The `@nb.njit` mode imposes restrictions: no Python objects (lists, dicts, classes), no dynamic typing, and limited library support. But for numerical code using NumPy arrays and scalar arithmetic (exactly the kind of code in an ODE integrator) these restrictions are easily met, and the payoff is enormous.

In our `neural_dynamics` package, every performance-critical function is decorated with `@nb.njit`:
- The HH rate functions (`alpha_n`, `beta_n`, etc.)
- The ODE right-hand side (`hh_rhs_parameterized`)
- The Euler and RK4 stepping functions (`euler_step`, `rk4_step`)
- The integration loop (`_solve_euler`, `_solve_rk4`)

Let us quantify the speedup with a concrete benchmark.

## 4.4 Validation: Reproducing Hodgkin & Huxley’s Results

Before exploring the model’s dynamics, we must first verify that our numerical implementation is correct. The gold standard is to reproduce the results from the original 1952 paper: specifically, the shape, amplitude, and duration of the action potential in response to a brief current stimulus.

We inject a **step current pulse** of $10 \; \mu\text{A/cm}^2$ for 1 ms into a resting membrane. If our implementation is correct, we should observe:

1. A rapid **depolarization** from the resting potential ($\approx -65 \; \text{mV}$) to a peak near $+40 \; \text{mV}$
2. A swift **repolarization** driven by $	ext{K}^+$ channel opening and $	ext{Na}^+$ channel inactivation
3. An **undershoot** (hyperpolarization) below the resting potential, caused by the lingering $	ext{K}^+$ conductance
4. A slow return to rest as the $	ext{K}^+$ channels close

<p align="center">
  <img src="figures/fig05_single_action_potential.png" width="700" alt="Single action potential evoked by a 1 ms current pulse, reproducing the classic HH waveform.">
  <br><em>Figure 5. Single action potential evoked by a 1 ms current pulse, reproducing the classic HH waveform.</em>
</p>

<p align="center">
  <img src="figures/fig06_four_state_variables.png" width="700" alt="All four state variables (V, m, h, n) during the action potential.">
  <br><em>Figure 6. All four state variables (V, m, h, n) during the action potential.</em>
</p>

**Validation against the original paper.** The waveform above closely matches the classic action potential shape from Hodgkin & Huxley (1952), Figure 12. Key features reproduced:

- **Peak amplitude** $\approx +40 \; \text{mV}$, approaching $E_\text{Na} = +50 \; \text{mV}$
- **Spike duration** $\approx 1.5 \; \text{ms}$ (measured at half-maximum amplitude)
- **Undershoot** to $\approx -80 \; \text{mV}$, transiently exceeding $E_\text{K} = -77 \; \text{mV}$ due to the slow deactivation of $	ext{K}^+$ channels
- **Temporal ordering** of gating variables: $m$ rises first (fast $	ext{Na}^+$ activation), then $h$ falls ($	ext{Na}^+$ inactivation) and $n$ rises ($	ext{K}^+$ activation)

The gating variable dynamics reveal the **biophysical mechanism** of the action potential:
1. The stimulus depolarizes the membrane past threshold
2. $m$ gates open rapidly $\rightarrow$ $	ext{Na}^+$ influx $\rightarrow$ positive feedback (more depolarization $\rightarrow$ more $m$ opening)
3. $h$ gates close (slower) $\rightarrow$ $	ext{Na}^+$ channels inactivate, terminating the inward current
4. $n$ gates open (slow) $\rightarrow$ $	ext{K}^+$ efflux $\rightarrow$ repolarization and undershoot
5. $n$ gates slowly close, $h$ gates slowly recover $\rightarrow$ return to rest

---

# 5. Simulating Neural Dynamics

Having validated our implementation, we can now use it as a *computational microscope* to explore the rich dynamical behavior of the Hodgkin-Huxley model. We will dissect the action potential into its constituent ionic currents, probe the threshold behavior that makes neural signaling digital, investigate the refractory periods that limit firing rates, and map out how the neuron encodes stimulus intensity into spike frequency.

## 5.1 Anatomy of an Action Potential

The action potential is not a single event but a sequence of **distinct biophysical phases**, each dominated by a different ionic current. By computing the individual currents $I_\text{Na}$, $I_\text{K}$, and $I_L$ from the simulation data, we can precisely identify which current drives each phase of the voltage waveform.

<p align="center">
  <img src="figures/fig07_action_potential_anatomy.png" width="700" alt="Detailed anatomy of the action potential with annotated phases.">
  <br><em>Figure 7. Detailed anatomy of the action potential with annotated phases.</em>
</p>

## 5.2 Threshold Behavior

One of the most striking features of the action potential is its **all-or-none** character. A stimulus that is too weak produces only a small, passive depolarization that decays back to rest. But once the stimulus exceeds a critical **threshold**, the full action potential is triggered with the same amplitude regardless of how far above threshold the stimulus is.

This threshold arises from the positive feedback loop in $	ext{Na}^+$ activation: depolarization opens $	ext{Na}^+$ channels $\rightarrow$ $	ext{Na}^+$ influx $\rightarrow$ further depolarization. Below threshold, the outward leak and $	ext{K}^+$ currents are strong enough to counteract the $	ext{Na}^+$ current and pull the membrane back to rest. Above threshold, the $	ext{Na}^+$ current wins, and the regenerative cycle drives the membrane to the peak of the spike.

Let us visualize this by injecting current pulses of increasing amplitude.

<p align="center">
  <img src="figures/fig08_threshold_behavior.png" width="700" alt="Threshold behavior: membrane responses to current pulses of varying amplitude.">
  <br><em>Figure 8. Threshold behavior: membrane responses to current pulses of varying amplitude.</em>
</p>

## 5.3 Refractory Periods

Immediately after an action potential, the neuron enters a **refractory period** during which it is either impossible or difficult to elicit a second spike. There are two phases:

- **Absolute refractory period** ($\sim 1\text{-}2 \; \text{ms}$): No stimulus, no matter how strong, can trigger another spike. This is because the $	ext{Na}^+$ inactivation gate $h$ is near zero; the channels are physically blocked, and cannot be opened until $h$ recovers.

- **Relative refractory period** ($\sim 3\text{-}5 \; \text{ms}$): A second spike *can* be triggered, but only by a **stronger-than-normal** stimulus. During this period, $h$ is partially recovered but $n$ is still elevated ($	ext{K}^+$ channels are still partly open), so the threshold is effectively raised.

We can demonstrate this with a **paired-pulse protocol**: deliver a first stimulus to trigger a spike, then deliver a second identical stimulus at varying delays.

<p align="center">
  <img src="figures/fig09_refractory_periods.png" width="700" alt="Paired-pulse protocol demonstrating absolute and relative refractory periods.">
  <br><em>Figure 9. Paired-pulse protocol demonstrating absolute and relative refractory periods.</em>
</p>

## 5.4 Repetitive Firing and Frequency-Current Relationship

When a neuron receives sustained suprathreshold input, it fires **repetitive action potentials**: a spike train. The rate of firing encodes information about the stimulus intensity. Understanding the relationship between injected current and firing rate is central to neural coding theory.

The **frequency-current (f-I) curve** maps the firing rate (in Hz) as a function of the constant injected current $I_\text{ext}$. For the Hodgkin-Huxley model, this curve has a distinctive property: firing onset is **discontinuous**. Below threshold, the firing rate is zero; at threshold, it jumps to a finite frequency (typically $\sim 50\text{-}70 \; \text{Hz}$). This is the hallmark of **Type II excitability**, in contrast to Type I neurons (such as those modeled by the Connor-Stevens equations) where firing can begin at arbitrarily low frequencies.

The biophysical origin of Type II behavior in the HH model is the coexistence of the resting state with a **limit cycle** (periodic orbit) at the bifurcation point: specifically, a subcritical Hopf bifurcation. We will analyze this in detail in the dynamical systems sections; here, we demonstrate it computationally.

<p align="center">
  <img src="figures/fig10_repetitive_firing.png" width="700" alt="Repetitive firing under constant suprathreshold current injection.">
  <br><em>Figure 10. Repetitive firing under constant suprathreshold current injection.</em>
</p>

<p align="center">
  <img src="figures/fig11_fi_curve.png" width="700" alt="Frequency-current (f-I) relationship of the Hodgkin-Huxley model.">
  <br><em>Figure 11. Frequency-current (f-I) relationship of the Hodgkin-Huxley model.</em>
</p>

## 5.5 Anode Break Excitation

One of the most counterintuitive predictions of the Hodgkin-Huxley model is **anode break excitation**: a neuron can fire an action potential *at the end* of a hyperpolarizing (inhibitory) stimulus, even though no depolarizing current is applied.

The mechanism is elegantly simple when understood through the gating variables:

1. **During hyperpolarization**: The membrane is held below rest. At these negative voltages, the $	ext{Na}^+$ inactivation gate $h$ has a *higher* steady-state value than at rest (recall from Section 3 that $h_\infty$ increases at negative voltages). So $h$ gradually increases. $	ext{Na}^+$ channels are being **de-inactivated**, priming them to open.

2. **Upon release**: When the hyperpolarizing current is removed, the membrane begins to return toward rest. But now $h$ is *larger than its resting value*: more $	ext{Na}^+$ channels are available to open. As $V$ rises past threshold (even passively, driven by the leak current), the enhanced $	ext{Na}^+$ availability triggers the positive feedback loop: $m$ gates open $\rightarrow$ $	ext{Na}^+$ influx $\rightarrow$ depolarization $\rightarrow$ more $m$ gates open.

3. **Result**: A full action potential fires, triggered not by a depolarizing stimulus but by the *removal* of a hyperpolarizing one.

This phenomenon is also called **post-inhibitory rebound** and plays important roles in neural circuits, particularly in central pattern generators for rhythmic motor behaviors.

<p align="center">
  <img src="figures/fig12_anode_break_excitation.png" width="700" alt="Anode break excitation: post-inhibitory rebound spike.">
  <br><em>Figure 12. Anode break excitation: post-inhibitory rebound spike.</em>
</p>

---

# 6. Phase Plane Analysis

Sections 4 and 5 explored the Hodgkin-Huxley model through *time-domain simulations*: we injected current, recorded voltage traces, and watched gating variables evolve. This approach reveals *what* the neuron does, but it does not fully explain *why* it does it. Why is there a sharp threshold? Why does repetitive firing start abruptly at a specific current? Why does the spike always have the same shape?

To answer these questions, we need a different perspective: not following a single trajectory through time, but mapping the **entire landscape of possible states** and understanding how the system flows through that landscape. This is the perspective of **dynamical systems theory**, and its primary visual tool is the **phase portrait**.

## 6.1 Motivation: Dimensionality Reduction

The HH model lives in a **four-dimensional state space** $(V, n, m, h)$. Visualizing trajectories in 4D is impossible for creatures confined to a 3D visual cortex. We need a principled way to reduce the dimensionality.

Fortunately, the HH system contains a natural **separation of timescales** that makes reduction possible:

1. **$m(t)$ is much faster than $n(t)$ and $h(t)$.** Recall from Section 3 that $\tau_m \approx 0.05\text{-}0.5 \; \text{ms}$, while $\tau_n \approx 1\text{-}8 \; \text{ms}$ and $\tau_h \approx 1\text{-}9 \; \text{ms}$. The fast variable $m$ reaches its steady-state $m_\infty(V)$ quasi-instantaneously relative to the slow dynamics of $V$, $n$, and $h$. We can therefore replace $m(t)$ with $m_\infty(V)$ everywhere; the **quasi-steady-state approximation**.

2. **$h$ and $n$ are approximately linearly related.** Hodgkin and Huxley themselves noted that, empirically, $h(t) + n(t) \approx 0.83$ throughout the action potential. This is not a coincidence: it reflects the fact that $h_\infty(V) + n_\infty(V)$ varies only weakly with $V$. We can therefore eliminate $h$ using the approximation $h \approx 0.83 - n$.

These two reductions collapse the 4D system to a **2D system in the $(V, n)$ plane**:

$$C_m \frac{dV}{dt} = -\bar{g}_\text{Na} \, m_\infty(V)^3 \, (0.83 - n) \, (V - E_\text{Na}) - \bar{g}_\text{K} \, n^4 \, (V - E_\text{K}) - g_L \, (V - E_L) + I_\text{ext}$$

$$\frac{dn}{dt} = \alpha_n(V)(1 - n) - \beta_n(V) \, n$$

In two dimensions, we can draw **phase portraits** (complete pictures of the system’s dynamics) and apply the full machinery of planar dynamical systems theory: nullclines, fixed points, stability analysis, and bifurcation diagrams.

## 6.2 FitzHugh-Nagumo as Conceptual Bridge

Before diving into the full HH reduction, it is instructive to consider a **simplified caricature** that captures the essential topology of excitability.

In 1961, Richard FitzHugh proposed a two-variable model inspired by the van der Pol oscillator, and in 1962, Jin-ichi Nagumo built an electronic circuit implementing the same equations:

$$\frac{dv}{dt} = v - \frac{v^3}{3} - w + I_\text{ext} \qquad \text{(fast excitatory variable)}$$

$$\frac{dw}{dt} = \varepsilon\,(v + a - bw) \qquad \text{(slow recovery variable)}$$

Here $v$ is analogous to the membrane potential, $w$ is analogous to a slow recovery variable (like $n$), and $\varepsilon \ll 1$ enforces the timescale separation. The cubic $v - v^3/3$ term produces the **N-shaped** nullcline that is the geometric hallmark of excitability.

The FitzHugh-Nagumo model demonstrates that **excitability is a topological property**, not a biological one. Any 2D system with a cubic-like fast nullcline and a monotone slow nullcline will exhibit threshold behavior, all-or-none spikes, and refractory periods. The HH model is simply a biophysically grounded realization of this universal geometry.

We will now construct the full $(V, n)$ phase portrait for the HH reduction, which has the same qualitative structure as FitzHugh-Nagumo but uses the actual biophysical rate functions.

## 6.3 Nullclines

A **nullcline** is the set of points in the phase plane where one of the state variables has zero rate of change. For our 2D reduction:

- **$V$-nullcline** ($dV/dt = 0$): The curve $n = f(V)$ where the membrane potential is instantaneously stationary. Using $m = m_\infty(V)$ and $h = h_\infty(V)$:

$$0 = -\bar{g}_\text{Na} \, m_\infty(V)^3 \, h_\infty(V) \, (V - E_\text{Na}) - \bar{g}_\text{K} \, n^4 \, (V - E_\text{K}) - g_L(V - E_L) + I_\text{ext}$$

  Solving for $n$:

$$n^4 = \frac{I_\text{ext} - \bar{g}_\text{Na} \, m_\infty^3 \, h_\infty \, (V - E_\text{Na}) - g_L(V - E_L)}{\bar{g}_\text{K}(V - E_\text{K})} \;\;\Rightarrow\;\; n = \left(\frac{\cdots}{\bar{g}_\text{K}(V-E_\text{K})}\right)^{1/4}$$

  This curve has a characteristic **cubic (N-shaped)** profile, the geometric origin of excitability.

- **$n$-nullcline** ($dn/dt = 0$): The curve where $n$ is stationary, which is simply $n = n_\infty(V)$. This is a **monotone sigmoid** that increases from 0 to 1.

**Fixed points** of the system are the **intersections** of the two nullclines, where both $dV/dt = 0$ and $dn/dt = 0$ simultaneously. The number and position of these intersections, and how they change with $I_\text{ext}$, determine the qualitative behavior of the system.

<p align="center">
  <img src="figures/fig13_nullclines_I0.png" width="700" alt="V- and n-nullclines at I_ext = 0, showing the resting state geometry.">
  <br><em>Figure 13. V- and n-nullclines at I_ext = 0, showing the resting state geometry.</em>
</p>

<p align="center">
  <img src="figures/fig14_nullclines_shift_Iext.png" width="700" alt="Nullcline shift with increasing I_ext, illustrating the transition to oscillation.">
  <br><em>Figure 14. Nullcline shift with increasing I_ext, illustrating the transition to oscillation.</em>
</p>

## 6.4 Linear Stability Analysis

Finding a fixed point tells us *where* the system can be stationary, but not *how* it behaves near that point. Does the system return to the fixed point after a small perturbation (stable), or does it diverge (unstable)? Does it spiral or approach monotonically?

To answer these questions, we **linearize** the system around the fixed point. Let $\mathbf{y}^*$ be a fixed point where $\mathbf{f}(\mathbf{y}^*) = 0$, and consider a small perturbation $\delta\mathbf{y} = \mathbf{y} - \mathbf{y}^*$. Taylor-expanding to first order:

$$\frac{d(\delta\mathbf{y})}{dt} = \mathbf{f}(\mathbf{y}^* + \delta\mathbf{y}) \approx \underbrace{\mathbf{f}(\mathbf{y}^*)}_{= \, 0} + \mathbf{J} \, \delta\mathbf{y}$$

where $\mathbf{J}$ is the **Jacobian matrix** evaluated at the fixed point:

$$J_{ij} = \frac{\partial f_i}{\partial y_j}\bigg|_{\mathbf{y} = \mathbf{y}^*}$$

The linearized system $\delta\mathbf{y}\,' = \mathbf{J}\,\delta\mathbf{y}$ has the general solution $\delta\mathbf{y}(t) = \sum_k c_k \, \mathbf{v}_k \, e^{\lambda_k t}$, where $\lambda_k$ and $\mathbf{v}_k$ are the eigenvalues and eigenvectors of $\mathbf{J}$. The **eigenvalues** therefore determine everything about the local dynamics:

| Eigenvalue type | Condition | Behavior | Classification |
|:---|:---|:---|:---|
| All real, all $\text{Re}(\lambda) < 0$ | $\lambda_k \in \mathbb{R}$, $\lambda_k < 0$ | Monotone decay | **Stable node** |
| All real, any $\text{Re}(\lambda) > 0$ | $\lambda_k \in \mathbb{R}$, $\exists \, \lambda_k > 0$ | Monotone growth | **Unstable node** |
| Real, mixed signs | $\lambda_1 < 0 < \lambda_2$ | Saddle dynamics | **Saddle point** |
| Complex, $\text{Re}(\lambda) < 0$ | $\lambda = \sigma \pm i\omega$, $\sigma < 0$ | Damped oscillations | **Stable spiral** |
| Complex, $\text{Re}(\lambda) > 0$ | $\lambda = \sigma \pm i\omega$, $\sigma > 0$ | Growing oscillations | **Unstable spiral** |
| Complex, $\text{Re}(\lambda) = 0$ | $\lambda = \pm i\omega$ | Sustained oscillations | **Center** (marginal) |

<p align="center">
  <img src="figures/fig15_eigenvalue_tracking.png" width="700" alt="Eigenvalue real parts as a function of I_ext, identifying the Hopf bifurcation.">
  <br><em>Figure 15. Eigenvalue real parts as a function of I_ext, identifying the Hopf bifurcation.</em>
</p>

## 6.5 Vector Field and Trajectories

Nullclines partition the phase plane into regions where $dV/dt$ and $dn/dt$ have definite signs. The **vector field** $\big(dV/dt, \, dn/dt\big)$ at each point $(V, n)$ shows the direction and magnitude of the system’s flow. Overlaying trajectories on this vector field reveals the geometry of spikes, threshold, and limit cycles in their full glory.

<p align="center">
  <img src="figures/fig16_phase_portrait_I0.png" width="700" alt="Phase portrait in the (V, n) plane at I_ext = 0 (resting state with perturbation).">
  <br><em>Figure 16. Phase portrait in the (V, n) plane at I_ext = 0 (resting state with perturbation).</em>
</p>

<p align="center">
  <img src="figures/fig17_phase_portrait_I10.png" width="700" alt="Phase portrait in the (V, n) plane at I_ext = 10 (stable limit cycle).">
  <br><em>Figure 17. Phase portrait in the (V, n) plane at I_ext = 10 (stable limit cycle).</em>
</p>

---

# 7. Dynamical Stability and Bifurcations

In Section 6, we saw that the fixed point of the HH system changes its stability as $I_\text{ext}$ increases: at rest ($I_\text{ext} = 0$), the fixed point is a **stable spiral**; at higher currents, eigenvalues cross the imaginary axis and the fixed point becomes an **unstable spiral**, with a limit cycle (periodic firing) appearing in its place.

This qualitative change in the system’s behavior at a critical parameter value is called a **bifurcation**. Bifurcation theory is the mathematical framework for classifying and predicting such transitions. In this section, we analyze the specific type of bifurcation in the HH model, its consequences for neural coding, and the remarkable phenomenon of hysteresis that it produces.

## 7.1 What is a Bifurcation?

A **bifurcation** occurs when a smooth, small change in a parameter causes a **topological change** in the system’s phase portrait: the number of fixed points changes, a fixed point changes stability, or a periodic orbit appears or disappears.

The key idea is that we are not just studying one dynamical system, but a **family** of dynamical systems parameterized by some control parameter (here, $I_\text{ext}$). As we tune this parameter, the vector field deforms continuously, but the *qualitative structure* of the flow can change discontinuously at isolated critical values. These critical values are the **bifurcation points**.

In the Hodgkin-Huxley model, $I_\text{ext}$ is the bifurcation parameter, and we have already seen evidence of a bifurcation: the transition from a stable resting state (no firing) to periodic spiking (repetitive firing) as current increases. The precise nature of this bifurcation has profound consequences for how the neuron encodes information.

## 7.2 Hopf Bifurcation in the HH Model

As $I_\text{ext}$ increases from zero, the fixed point of the HH system undergoes the following sequence of changes:

1. **Low current** ($I_\text{ext} \lesssim 6$ $\mu\text{A/cm}^2$): The fixed point is a **stable spiral**. The eigenvalues are complex with negative real parts: $\lambda = \sigma \pm i\omega$ with $\sigma < 0$. Perturbations produce **damped oscillations** that decay back to rest. The neuron does not fire.

2. **Critical current** ($I_\text{ext} \approx 6\text{-}10$ $\mu\text{A/cm}^2$): The real part of the complex eigenvalues crosses zero: $\sigma = 0$. At this instant, the fixed point is a **center**: perturbations neither grow nor decay.

3. **Above critical current**: $\sigma > 0$, and the fixed point becomes an **unstable spiral**. Perturbations grow, and the trajectory is captured by a **stable limit cycle**: the periodic spike train.

This scenario (a fixed point losing stability through complex eigenvalues crossing the imaginary axis, accompanied by the birth of a limit cycle) is called a **Hopf bifurcation**.

### Supercritical vs. Subcritical Hopf

There are two flavors of Hopf bifurcation, with dramatically different physical consequences:

- **Supercritical Hopf**: The limit cycle is born *at* the bifurcation with **zero amplitude** and grows smoothly as the parameter moves past the critical value. Oscillations appear gradually.

- **Subcritical Hopf**: An *unstable* limit cycle exists below the bifurcation point, shrinks as the parameter approaches the critical value, and collides with the fixed point at the bifurcation. Beyond the critical point, the system jumps to a **large-amplitude** limit cycle that already existed further from the fixed point. The onset of oscillations is **sudden and discontinuous**.

The HH model exhibits a **subcritical Hopf bifurcation**. This has two critical consequences:

1. **Discontinuous onset of firing**: The neuron jumps from silence to finite-frequency spiking, with no intermediate regime of infinitesimally small oscillations.

2. **Hysteresis**: Because the large-amplitude limit cycle exists even below the bifurcation point (in a range where both the stable fixed point and the limit cycle coexist), the system exhibits **bistability**. Once firing is established, it can persist even when the current is reduced below the onset threshold.

<p align="center">
  <img src="figures/fig18_bifurcation_diagram.png" width="700" alt="Bifurcation diagram showing fixed point branch and limit cycle extrema.">
  <br><em>Figure 18. Bifurcation diagram showing fixed point branch and limit cycle extrema.</em>
</p>

## 7.3 Type I vs. Type II Excitability

The type of bifurcation at firing onset has a profound consequence for how the neuron encodes stimulus intensity into firing rate. Neurons are classified into two fundamental types based on their **frequency-current (f-I) relationship**:

### Type I Excitability (SNIC Bifurcation)

In a **Saddle-Node on an Invariant Circle (SNIC)** bifurcation, a stable fixed point and a saddle point collide and annihilate on a periodic orbit. After the bifurcation, the trajectory must traverse the “ghost” of the vanished fixed points, which slows it down. The result:

- Firing rate starts from **arbitrarily close to zero** at threshold
- $f(I) \propto \sqrt{I - I_\text{th}}$ near threshold (a square-root law)
- The neuron can fire at **any frequency** above zero
- Continuous onset of firing

### Type II Excitability (Hopf Bifurcation)

In a **Hopf bifurcation** (especially subcritical), a stable fixed point loses stability to an oscillation with a characteristic frequency $\omega$. The result:

- Firing starts at a **nonzero minimum frequency** $f_\text{min} > 0$
- **Discontinuous jump** from silence to $f_\text{min}$ at threshold
- The neuron **cannot** fire below $f_\text{min}$
- The HH model is **Type II**, as we demonstrated in Section 5.4

This classification has deep implications for neural coding. Type I neurons are better suited for **rate coding** (smoothly encoding stimulus intensity). Type II neurons are better at **coincidence detection** (responding to synchronous inputs, since they prefer inputs at their resonant frequency).

<p align="center">
  <img src="figures/fig19_fi_typeI_vs_typeII.png" width="700" alt="f-I curve comparison: Type II (HH) vs Type I (theoretical SNIC) excitability.">
  <br><em>Figure 19. f-I curve comparison: Type II (HH) vs Type I (theoretical SNIC) excitability.</em>
</p>

## 7.4 Hysteresis and Bistability

The subcritical nature of the HH Hopf bifurcation creates a parameter regime where **two stable attractors coexist**: the resting fixed point and the spiking limit cycle. This bistability produces **hysteresis**: the system’s behavior depends not only on the current value of $I_\text{ext}$, but also on its *history*.

Concretely:
- When $I_\text{ext}$ is **ramped up** from zero, firing begins at a current $I_\text{on}$.
- When $I_\text{ext}$ is **ramped down** from a high value, firing persists below $I_\text{on}$ and ceases at a lower current $I_\text{off} < I_\text{on}$.

The difference $I_\text{on} - I_\text{off}$ is the **hysteresis width**, and its existence is a direct consequence of the subcritical bifurcation structure. This phenomenon has been observed experimentally in real neurons and has implications for bistable neural circuits and memory.

<p align="center">
  <img src="figures/fig20_hysteresis.png" width="700" alt="Hysteresis demonstration via slow current ramp up and down.">
  <br><em>Figure 20. Hysteresis demonstration via slow current ramp up and down.</em>
</p>

### Summary of Dynamical Systems Analysis

We have now established the complete dynamical systems picture of the Hodgkin-Huxley model:

1. **Phase plane reduction**: By exploiting the timescale separation ($m$ fast, $n$ and $h$ slow) and the approximate linear relationship $h \approx 0.83 - n$, we reduced the 4D HH system to a 2D system in the $(V, n)$ plane.

2. **Nullclines and fixed points**: The $V$-nullcline has an N-shaped (cubic-like) profile, and the $n$-nullcline is a monotone sigmoid. Their intersection determines the fixed point. As $I_\text{ext}$ increases, the $V$-nullcline shifts upward, moving the fixed point along the nullclines.

3. **Linear stability**: The Jacobian eigenvalues reveal that the resting state is a **stable spiral** at low currents. As $I_\text{ext}$ increases, the real parts of complex eigenvalue pairs cross zero; a **Hopf bifurcation**.

4. **Subcritical Hopf**: The bifurcation is subcritical, meaning firing onset is **discontinuous** (Type II excitability) and the system exhibits **hysteresis** (bistability between rest and spiking).

5. **Type II excitability**: The f-I curve starts at a nonzero frequency, in contrast to Type I neurons (SNIC bifurcation) that can fire arbitrarily slowly. This reflects the biophysical properties of the HH model’s ionic currents.

These results illustrate a profound principle: **the qualitative behavior of a neural model is determined by its bifurcation structure**, not by its specific biophysical details. Different ionic mechanisms can produce the same type of bifurcation (and hence the same firing phenomenology), while changes in a single parameter can switch the system between fundamentally different dynamical regimes.

---

# 8. The Leaky Integrate-and-Fire Model

## 8.1 From Biophysics to Simplicity

The Hodgkin-Huxley model is a triumph of mathematical biophysics, but it comes at a cost: **four coupled nonlinear ODEs** that require careful numerical integration with small time steps. For many applications in computational neuroscience (large-scale network simulations, machine learning with spiking neurons, theoretical analysis of neural coding) this level of detail is neither necessary nor computationally feasible.

The **Leaky Integrate-and-Fire (LIF)** model occupies the opposite end of the realism-tractability spectrum. It retains the essential feature of neural computation (spike generation through threshold crossing) while discarding all voltage-dependent conductance dynamics. The result is a **single linear ODE** with a discontinuous reset rule: orders of magnitude cheaper to simulate, and often analytically solvable.

The fundamental trade-off:

| Property | Hodgkin-Huxley | LIF |
|:---|:---|:---|
| State variables | 4 ($V, n, m, h$) | 1 ($V$) |
| ODE type | Nonlinear, coupled | Linear, scalar |
| Spike mechanism | Emergent from ion channels | Imposed threshold rule |
| Computational cost per spike | High | Very low |
| Analytical tractability | Limited | Extensive |

## 8.2 Derivation from RC Circuit

Strip the membrane down to its simplest form: a **passive RC circuit**. No voltage-dependent conductances, no gating variables: just a capacitor (the lipid bilayer) in parallel with a resistor (the combined leak conductance):

$$C_m \frac{dV}{dt} = -\frac{V - V_\text{rest}}{R_m} + I_\text{ext}$$

where $R_m$ is the membrane resistance (inverse of the total leak conductance) and $V_\text{rest}$ is the resting potential. Defining the **membrane time constant** $\tau_m = R_m C_m$, we obtain the standard LIF equation:

$$\tau_m \frac{dV}{dt} = -(V - V_\text{rest}) + R_m I_\text{ext}$$

This is a **linear first-order ODE** (unlike the HH system, which is nonlinear!). For constant $I_\text{ext}$, the analytical solution exists in closed form:

$$V(t) = V_\text{rest} + R_m I_\text{ext}\left(1 - e^{-t/\tau_m}\right)$$

The voltage exponentially approaches the asymptotic value $V_\infty = V_\text{rest} + R_m I_\text{ext}$ with time constant $\tau_m$. If $V_\infty < V_\text{th}$, the neuron never fires (subthreshold regime). If $V_\infty \geq V_\text{th}$, the voltage will reach threshold in finite time, and we need the fire-and-reset rule.

## 8.3 The Fire-and-Reset Rule

The RC circuit alone does not spike: it merely charges toward an asymptote. To produce discrete spike events, we impose a **discontinuous threshold rule**:

> When $V \geq V_\text{th}$: **(1)** record a spike at time $t_\text{spike}$, **(2)** reset $V \to V_\text{reset}$, **(3)** enforce an absolute refractory period $t_\text{ref}$ during which the voltage is clamped at $V_\text{reset}$.

The spike itself has **no shape** in the LIF model: it is a mathematical point event (a Dirac delta function in the spike train). This is the most fundamental difference from the HH model, where the action potential waveform emerges naturally from the ionic current dynamics.

**Default parameters:**

| Parameter | Symbol | Value | Units |
|:---|:---|:---|:---|
| Resting potential | $V_\text{rest}$ | $-65$ | mV |
| Membrane time constant | $\tau_m$ | $10$ | ms |
| Membrane resistance | $R_m$ | $10$ | M$\Omega$ |
| Spike threshold | $V_\text{th}$ | $-50$ | mV |
| Reset potential | $V_\text{reset}$ | $-65$ | mV |
| Refractory period | $t_\text{ref}$ | $2$ | ms |

The threshold current (rheobase) is $I_\text{rheo} = (V_\text{th} - V_\text{rest}) / R_m = 15/10 = 1.5$ nA.

## 8.4 LIF Implementation and Simulation

Let us simulate the LIF model using the `lif_simulate` function from our package and visualize the characteristic **sawtooth waveform** of the LIF membrane potential.

<p align="center">
  <img src="figures/fig21_lif_simulation.png" width="700" alt="Leaky Integrate-and-Fire model simulation.">
  <br><em>Figure 21. Leaky Integrate-and-Fire model simulation.</em>
</p>

<p align="center">
  <img src="figures/fig22_hh_vs_lif_sidebyside.png" width="700" alt="Side-by-side comparison of HH and LIF responses to the same stimulation protocol.">
  <br><em>Figure 22. Side-by-side comparison of HH and LIF responses to the same stimulation protocol.</em>
</p>

## 8.5 What LIF Discards

The simplification from HH to LIF comes at a significant cost in biophysical detail:

1. **No spike shape.** In HH, the action potential waveform (rapid depolarization, overshoot to +40 mV, repolarization, undershoot to -80 mV) emerges from the interplay of $	ext{Na}^+$ and $	ext{K}^+$ currents. In LIF, the spike is a dimensionless point event; the voltage jumps discontinuously from $V_\text{th}$ to $V_\text{reset}$.

2. **No $	ext{Na}^+$ dynamics or inactivation.** The rapid activation and subsequent inactivation of sodium channels (the mechanism that terminates the action potential) is entirely absent. There is no gating variable $m$ or $h$.

3. **No natural refractory period.** In HH, the refractory period arises naturally from $	ext{Na}^+$ channel inactivation ($h \to 0$) and $	ext{K}^+$ channel activation ($n \to 1$). In LIF, refractoriness must be imposed artificially as a hard-coded time window.

4. **No subthreshold nonlinearities.** The HH model exhibits rich subthreshold dynamics: resonance, subthreshold oscillations, and a voltage-dependent effective time constant. The LIF subthreshold dynamics are purely linear.

5. **No anode break excitation.** The post-inhibitory rebound spike we observed in the HH model (Section 5.5) relies on the differential recovery rates of $	ext{Na}^+$ and $	ext{K}^+$ channels. LIF cannot exhibit this phenomenon.

Despite these limitations, the LIF model captures the essential feature of neural spiking (**threshold-triggered discrete events**) and does so with remarkable efficiency. In the next section, we quantify exactly how the two models compare.

---

# 9. Quantitative Comparison: HH vs LIF

Having introduced both models, we now compare them head-to-head on four quantitative dimensions: firing rate, spike shape, refractory behavior, and computational cost.

## 9.1 Firing Rate Comparison: f-I Curves

The frequency-current (f-I) curve is the single most important input-output characterization of a neuron model. It reveals the neuron's **excitability type**: Type I neurons can fire at arbitrarily low rates (continuous onset), while Type II neurons begin firing at a nonzero minimum frequency (discontinuous onset).

<p align="center">
  <img src="figures/fig23_fi_curves_overlaid.png" width="700" alt="f-I curves overlaid: HH (Type II) vs LIF (Type I).">
  <br><em>Figure 23. f-I curves overlaid: HH (Type II) vs LIF (Type I).</em>
</p>

## 9.2 Spike Shape Comparison

A single spike reveals the most striking qualitative difference between the two models.

<p align="center">
  <img src="figures/fig24_spike_shape_comparison.png" width="700" alt="Single spike shape comparison between HH and LIF models.">
  <br><em>Figure 24. Single spike shape comparison between HH and LIF models.</em>
</p>

## 9.3 Refractory Period Comparison

In the HH model, the refractory period emerges naturally from the slow recovery of $	ext{Na}^+$ inactivation ($h$) and the slow deactivation of $	ext{K}^+$ channels ($n$). The recovery is **continuous**: at short interpulse intervals, the second spike is smaller and slower; at longer intervals, it gradually recovers to full amplitude.

In the LIF model, the refractory period is a **binary switch**: during $t_\text{ref}$, the neuron is completely unable to fire; after $t_\text{ref}$, it fires with exactly the same characteristics as before.

<p align="center">
  <img src="figures/fig25_refractory_comparison.png" width="700" alt="Paired-pulse refractory period comparison: HH vs LIF.">
  <br><em>Figure 25. Paired-pulse refractory period comparison: HH vs LIF.</em>
</p>

## 9.4 Computational Cost Benchmark

For network simulations involving thousands of neurons, computational cost per neuron is critical. Here we benchmark 10 seconds of simulation time for each model and method combination.

<p align="center">
  <img src="figures/fig26_computational_cost_benchmark.png" width="700" alt="Computational cost benchmark: HH vs LIF.">
  <br><em>Figure 26. Computational cost benchmark: HH vs LIF.</em>
</p>

## 9.5 When to Use Which?

The choice between HH and LIF depends on the scientific question:

| Criterion | Hodgkin-Huxley | Leaky Integrate-and-Fire |
|:---|:---:|:---:|
| Spike shape matters | ✓ | ✗ |
| Network of >1000 neurons | ✗ | ✓ |
| Ion channel pharmacology | ✓ | ✗ |
| Refractory mechanisms | Natural (continuous) | Artificial (binary) |
| Computational cost | High | Very low |
| Phase plane analysis | Rich structure | Trivial (1D) |
| Analytical tractability | Limited | Extensive |
| Subthreshold dynamics | Nonlinear, resonant | Linear, passive |
| Anode break excitation | ✓ | ✗ |
| Excitability type | Type II (Hopf) | Type I (threshold) |

**Rule of thumb:** Use HH (or similar conductance-based models) when the question involves **ionic mechanisms, spike waveforms, or dynamical systems structure**. Use LIF (or similar threshold models) when the question involves **network-level computation, coding, or learning** and the internal spike dynamics are not the focus.

Many intermediate models exist (the Izhikevich model, the adaptive exponential integrate-and-fire (AdEx), the FitzHugh-Nagumo model), each offering different trade-offs along the realism-tractability axis.

---

# 10. Conclusions and References

## 10.1 Summary

This notebook has traced a path from the **biophysics of ion channels** to the **mathematics of dynamical systems**, using the Hodgkin-Huxley model as our guide. The key insights are:

1. **The neuron is an electrical circuit with voltage-dependent conductances.** The lipid bilayer acts as a capacitor, ion channels as variable resistors, and concentration gradients as batteries. This equivalent circuit naturally gives rise to the membrane equation $C_m \, dV/dt = -\sum I_\text{ion} + I_\text{ext}$.

2. **The HH model captures spike generation through four coupled nonlinear ODEs.** The four state variables $(V, n, m, h)$ describe the membrane potential and the three gating variables governing $	ext{Na}^+$ and $	ext{K}^+$ channel kinetics. Despite its apparent simplicity (four equations, seven parameters), this system reproduces the full richness of neural excitability: threshold behavior, all-or-nothing spikes, refractory periods, repetitive firing, and anode break excitation.

3. **Phase plane analysis reveals the topological structure of excitability.** By exploiting timescale separation, we reduced the 4D system to a 2D phase portrait in the $(V, n)$ plane. The N-shaped $V$-nullcline and sigmoid $n$-nullcline organize the dynamics: their intersection determines the fixed point, and the shape of the nullclines explains the existence of threshold behavior and all-or-nothing responses.

4. **Hopf bifurcation separates resting from oscillating regimes.** As the applied current $I_\text{ext}$ increases, the fixed point loses stability through a (subcritical) Hopf bifurcation, and a stable limit cycle (repetitive firing) emerges. The subcritical nature explains the discontinuous onset of firing (Type II excitability) and the hysteresis observed in current ramp experiments.

5. **Simpler models (LIF) trade realism for computational efficiency.** The leaky integrate-and-fire model reduces the neuron to a single linear ODE with a threshold rule, achieving orders-of-magnitude speedups at the cost of biophysical detail. The choice between detailed and simplified models depends on the scientific question.

## 10.2 References

1. Hodgkin, A.L. & Huxley, A.F. (1952). "A quantitative description of membrane current and its application to conduction and excitation in nerve." *J. Physiol.* **117**(4): 500-544.

2. FitzHugh, R. (1961). "Impulses and Physiological States in Theoretical Models of Nerve Membrane." *Biophysical Journal* **1**(6): 445-466.

3. Nagumo, J., Arimoto, S., & Yoshizawa, S. (1962). "An Active Pulse Transmission Line Simulating Nerve Axon." *Proc. IRE* **50**(10): 2061-2070.

4. Izhikevich, E.M. (2007). *Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting.* MIT Press.

5. Dayan, P. & Abbott, L.F. (2001). *Theoretical Neuroscience: Computational and Mathematical Modeling of Neural Systems.* MIT Press.

6. Gerstner, W., Kistler, W.M., Naud, R., & Paninski, L. (2014). *Neuronal Dynamics: From Single Neurons to Networks and Models of Cognition.* Cambridge University Press.

7. Ermentrout, G.B. & Terman, D.H. (2010). *Mathematical Foundations of Neuroscience.* Springer.

8. Koch, C. (1999). *Biophysics of Computation: Information Processing in Single Neurons.* Oxford University Press.

---

*Notebook created as part of the Neural Membrane Dynamics project. All simulations use custom Numba-compiled solvers for performance.*

---

## Project Structure

```
neural-membrane-dynamics/
├── notebooks/
│   └── hodgkin_huxley_dynamics.ipynb
├── src/neural_dynamics/
│   ├── hodgkin_huxley.py
│   ├── integrate_and_fire.py
│   ├── solvers.py
│   └── analysis.py
├── tests/
├── figures/
├── requirements.txt
└── README.md
```

## Installation

```bash
git clone https://github.com/Mixnikon108/neural-membrane-dynamics.git
cd neural-membrane-dynamics
pip install -r requirements.txt
```

**Dependencies:** NumPy, Numba, Matplotlib, SciPy, Jupyter, pytest.

## Running Tests

```bash
PYTHONPATH=src pytest tests/ -v
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.