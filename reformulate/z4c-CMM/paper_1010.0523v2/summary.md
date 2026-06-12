# summary.md — arXiv:1010.0523v2 (Ruiz, Hilditch, Bernuzzi 2010)

## Motivation
Numerical relativity introduces an artificial timelike outer boundary; the boundary
conditions must yield a well-posed IBVP and must not inject constraint violations.
Z4c propagates every constraint at the speed of light, which reduces bulk constraint
violation but makes it acutely sensitive to boundary treatment: Sommerfeld conditions
(standard with BSSN) are not constraint preserving and produce non-convergent boundary
artifacts that can masquerade as physics (e.g. permanently shifting a star's central density).

## Goal
Specify high order absorbing constraint preserving boundary conditions (CPBCs) for Z4c
coupled to the moving-puncture gauge family, prove well-posedness where tractable, and
demonstrate their efficacy numerically.

## Result scope
- CPBCs: (r^2 l-ring^a d_a)^L Theta = 0, (r^2 l-ring^a d_a)^L Z_i = 0 on the boundary, plus six gauge/metric conditions (eq:general_CPBCs - eq:BCs_lastII), generalizing Bona et al. (L = 0).
- Well-posedness (boundary stability) of the constraint subsystem proven in the frozen coefficient approximation on a four dimensional compact manifold, for any order L, via Laplace-Fourier / Kreiss theory.
- Well-posedness of the remaining (gauge + metric) system shown only for a spherical reduction with 1+log slicing, mu_S = 1, and (for the cascade structure) the new "asymptotically harmonic" shift (eps_alpha, eps_chi) = (1, 1/2).
- Numerical evidence in spherical symmetry (flat perturbed, stable star, puncture and Kerr-Schild black holes): 2nd order CPBCs (L = 1) absorb outgoing constraint violation, retain 4th order convergence in reflections, improve the constraint monitor by 2-4 orders of magnitude over Sommerfeld; 1st order CPBCs reflect everything (R ~ 1) and are discarded.

## Conclusions
CPBCs for puncture-gauge Z4c are practical, nearly Sommerfeld-like to implement, and
essential for reliable Z4c evolutions. The constraint subsystem IBVP is boundary stable
for all orders L. First demonstrated use of the asymptotically harmonic shift on puncture data.

## Key challenge
The system is not symmetric hyperbolic and the CPBCs are not maximally dissipative, so
energy methods fail; the Kreiss-Winicour cascade method also fails for the standard
puncture gauge (gauge and metric sectors stay coupled), forcing the frozen coefficient
approximation and, beyond the constraint subsystem, a spherical reduction.

## Method innovation
- New eps_alpha, eps_chi terms in the Gamma-driver shift (asymptotically harmonic shift) producing a cascade structure that decouples the gauge sector in the bulk.
- First order pseudo-differential reduction of fully second order wave problems and iteration of the boundary operator to algebraic form (A^{L+1}), enabling boundary-stability estimates for arbitrary-order BCs.
- Implementation of CPBCs as Sommerfeld-like replacement evolution equations at the boundary (Appendix C), with an experimental reflection coefficient diagnostic.

## Possible bottleneck
Results beyond the constraint subsystem hold only in spherical symmetry and the
high-frequency/frozen-coefficient limit; no symmetrizer is exhibited for the standard
puncture gauge; 3D tests (BAM) indicate discarded tangential terms are required for
stability; numerical stability of the discretization is not proven (no SBP energy).
