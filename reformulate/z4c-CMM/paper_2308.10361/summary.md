# summary.md — arXiv:2308.10361

## Motivation
Cauchy-characteristic evolution (CCE) extracts unambiguous waveforms at future null
infinity but is one-way: the Cauchy evolution never receives information from the
characteristic system, so nonlinear backscatter off spacetime curvature beyond the
outer boundary is lost. Production SpEC simulations freeze psi0 = 0 at the boundary,
eliminating all backscatter and admitting spurious reflections.

## Goal
Implement and test fully relativistic, three-dimensional Cauchy-characteristic
matching (CCM) for the *physical* degrees of freedom — using the characteristic
system's psi0 to supply the boundary value w-|_BC of the physical Bjorhus boundary
condition of the GH Cauchy system in SpECTRE — without symmetry assumptions,
perturbative matching, or linearization.

## Result scope
- Matching of the two physical incoming modes only; the four incoming gauge modes
  keep Sommerfeld conditions (future work). Constraint-preserving BCs are untouched
  (matching them is neither necessary nor well-motivated).
- Two equivalent construction choices: from partially flat Bondi-like coordinates
  (Choice 1: Type I transformation of m-hat plus Type II of psi0-hat) or from
  Bondi-like coordinates (Choice 2: trivial Type I for m plus Type II of psi0).
  Choice 2 is used in production (single Lorentz transformation).
- Tests: Teukolsky wave on flat background (X = 1e-5 perturbative, X = 2 nonlinear),
  Teukolsky-wave perturbation of a Kerr BH (chi = 0.5), and injection of a GW pulse
  from the characteristic grid into an initially Minkowski Cauchy domain.

## Conclusions
- No numerical instabilities observed in any test (smooth data, >1000 code units).
- CCM leaves GH gauge/three-index constraint violations unchanged relative to CCE,
  but systematically reduces Bondi-gauge constraint violations and brings nonlinear
  (X = 2) waveforms about one order of magnitude closer to a causally disconnected
  reference run.
- The characteristic-pulse test confirms the boundary is transparent to incoming
  radiation under CCM (vs perfectly reflecting under CCE).
- CCM permits smaller Cauchy domains (cheaper simulations) at fixed waveform quality.

## Key challenge
The Cauchy (GH) and characteristic (CCE) systems use different coordinates and
different null tetrads. Feeding psi0 back requires a chain of Jacobians across four
coordinate systems (Cauchy -> null-radius -> Bondi-like -> partially flat
Bondi-like), the classification of allowed tetrad transformations (Type I null
rotations leave psi0 invariant; Type II boosts rescale it by A^2), and
time-dependent angular interpolation between the characteristic and Cauchy grids.

## Method innovation
- Observation that the characteristic outgoing null vector l-hat is *by construction*
  proportional to the Cauchy l at the worldtube, so only a Type II boost
  (computable from alpha, beta^i, s^j, beta-hat) relates psi0 to psi0', and all
  Type-I ambiguities can be dropped (the "approx" calculus).
- Phase freedom Theta of m provably cancels in w-, so no angular-tetrad alignment
  is needed.
- Cheap inverse angular map: evolve x^A(u-hat, x-hat) by d_uhat x^A = U^(0)A
  (Cartesian, spin-weighted form) alongside the forward map instead of inverting
  numerically.

## Possible bottlenecks
- The full Bondi-like characteristic system is only weakly hyperbolic, so CCM is
  not well-posed; instabilities may appear for non-smooth/high-frequency data or
  BBH mergers (untested). [OPEN]
- Gauge boundary conditions remain Sommerfeld; spurious reflection (tertiary pulse
  at t = 2 r_out + (r_out - r_c)) persists and limits how far inward the boundary
  can be moved. [FUTURE]
- The perturbative tests cannot resolve CCE-vs-CCM differences because backscatter
  (psi0 ~ r^-5) is tiny at the chosen boundary radii.
