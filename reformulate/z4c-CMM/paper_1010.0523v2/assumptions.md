# assumptions.md — arXiv:1010.0523v2

Enumerated symmetries / limits / regimes. Markers per `_common/contracts/markers.md`.

1. [ASSUMPTION] Vacuum or perfect-fluid GR with geometric units c = G = 1; signature (-,+,+,+); matter enters only through (rho, S^i, S_ij).
2. [ASSUMPTION] Constraint damping discarded: the Gundlach et al. damping scheme of Z4 is switched off throughout the analysis and the numerics (kappa_1 = kappa_2 = 0).
3. [ASSUMPTION] Algebraic constraints D = ln det(gamma-tilde) = 0 and T = gamma-tilde^ij A-tilde_ij = 0 are imposed continuously (projected) during evolution.
4. [ASSUMPTION] Puncture gauge family eq:punc_alpha / eq:punc_beta without the auxiliary B^i field; analysis fixes 1+log slicing mu_L = 2/alpha and mu_S = 1; standard numerical gauge uses mu_L = 2/alpha, mu_S = 3/4, eps_alpha = eps_chi = 0.
5. [ASSUMPTION] Asymptotically harmonic shift (mu_S, eps_alpha, eps_chi) = (1, 1, 1/2) for the cascade (decoupled gauge sector) analytic results; lambda != 0 required, else only weak hyperbolicity.
6. [ASSUMPTION] Frozen coefficient (high-frequency) approximation: coefficients of equations of motion and boundary operators frozen at a point p; linear constant-coefficient half-space problem, flat boundary x = 0 with frozen normal shift beta-ring; well-posedness here is necessary, and only expected to transfer to the nonlinear problem.
7. [ASSUMPTION] Background metric for boundary constructions: conformally flat 3-metric (eq:auto-16); normal vectors n-ring, s-ring defined against the background, not the physical metric.
8. [ASSUMPTION] Physical and background metrics sufficiently close to flat at the boundary so that exactly ten characteristic variables are incoming (fixes the number of BCs).
9. [ASSUMPTION] Linear regime around the background for constraint preservation/absorption statements; principal-part (\simeq) equalities only for the wave-equation structure.
10. [ASSUMPTION] Trivial initial data in all IBVP proofs (general data reducible by the u-bar = u - g(t) f(x) transformation, Gustafsson et al.).
11. [ASSUMPTION] L2-class solutions in Laplace-Fourier space with Re(s) = eta > 0; eigenvectors normalized to remain finite as omega -> 0, +-infinity and |s| -> infinity.
12. [ASSUMPTION] Strong hyperbolicity of puncture-gauge Z4c assumed generic ("except a handful of special cases" not discussed); Kreiss theory strictly applies to strictly hyperbolic systems, Agranovich extension for constant multiplicity.
13. [ASSUMPTION] Spherical symmetry for the gauge/metric-sector well-posedness analysis and for all numerical tests; spherical reduction of Z4c per Bernuzzi-Hilditch Appendix A.
14. [ASSUMPTION] Boundary conditions eq:BCs-alpha - eq:BCs_lastII are not tailored to absorb all outgoing fields (absorption of gravitational radiation not addressed); tangential boundary terms discarded (3D tests show they matter).
15. [ASSUMPTION] Numerical setup: spherical code of Bernuzzi-Hilditch, 4th order finite differences, Kreiss-Oliger dissipation, sixth order ghostzone extrapolation (eq:auto-36); BCs implemented by replacing evolution equations at the boundary, linearized around flat space (eq:sph_bc_1st - eq:sph_bc_last); reference solution with causally disconnected boundary (r'_out ~ 1000 M) as truth proxy.
16. [ASSUMPTION] Norm comparisons restricted to r in (0, r_out); staggered grids used in applications although only non-staggered grids show clean 4th order boundary convergence.
