# assumptions.md — arXiv:2308.10361

Enumerated symmetries / limits / regimes used by the paper. Each carries a marker
per `_common/contracts/markers.md`.

1. [ASSUMPTION] Vacuum general relativity: R_μ'ν' = 0. No matter sources anywhere
   in the matched system.
2. [ASSUMPTION] GH evolution system of Lindblom et al. 2005 (first-order, 50
   variables) with constraint damping γ₂; characteristic-field expressions and
   constraint-preserving BCs imported unmodified.
3. [ASSUMPTION] "Typical shift values" at the outer boundary: exactly forty
   incoming characteristic fields. The counting can change if β^i s_i is large.
4. [ASSUMPTION] Outer boundary only: inner (excision) boundaries are out of scope
   (paper footnote in the introduction).
5. [ASSUMPTION] Gauge boundary conditions are NOT matched: Sommerfeld conditions
   (Rinne et al. 2007) are used for the four incoming gauge modes. CCM here covers
   only the two physical modes. Final accuracy claims are conditional on this.
6. [ASSUMPTION] Characteristic region free of caustics of the outgoing null
   congruence: the worldtube must be far enough out that Bondi-like coordinates
   exist (null rays from the worldtube remain non-crossing).
7. [ASSUMPTION] Falloff regime: hatted metric functions satisfy
   eq:falloff_partial_inertial (Ŵ, Û^Â = O(r̂^-2), ĥ = q + O(r̂^-1)); β̂ need not
   fall off (partially flat, not true Bondi-Sachs). Waveforms at scri are computed
   after transformation to true Bondi-Sachs coordinates (up to BMS).
8. [ASSUMPTION] The characteristic (and hence CCM) system is only weakly
   hyperbolic; well-posedness of the continuum matched problem is NOT established.
   All stability statements are empirical, for smooth data only.
9. [ASSUMPTION] Worldtube = Cauchy outer boundary (r'_out), i.e. matching surface
   coincides with the boundary where the Bjørhus condition is imposed.
10. [ASSUMPTION] Smooth test data: Teukolsky waves (Gaussian profiles), Kerr +
    small pulse (X = 0.01), characteristic pulse small enough not to collapse
    (Z = 1e-3). No high-frequency or non-smooth robustness tests.
11. [ASSUMPTION] Perturbative regime for analytic comparisons: X = 1e-5 such that
    the linearized Teukolsky solution (eq:Teukolsky_metric) is the reference;
    nonlinear terms O(X²) neglected in the reference only.
12. [ASSUMPTION] For the X = 2 test, the "exact" reference is a SpEC CCE run with
    causally disconnected outer boundary (r'_ref = 900) — itself a numerical
    solution of a weakly hyperbolic system, justified by prior accuracy studies.
13. [ASSUMPTION] Type I tetrad equivalence ("≈" calculus): terms proportional to
    the outgoing null vector l are dropped because ψ₀ and w− are invariant under
    Type I transformations; validity requires l^a' ∝ l^â, proven at the worldtube
    only (eq:l_GH_and_l_CCE_transform).
14. [ASSUMPTION] Phase gauge freedom Θ of m is physically irrelevant: w− invariant
    under m → m e^{iΘ} because ψ₀ co-transforms (eq:lorentz_psi0_ii).
15. [ASSUMPTION] Angular maps x̂^Â(u, x^A) and x^A(û, x̂^Â) remain smooth,
    invertible diffeomorphisms of S² during evolution (ω̂ > 0, b̂b̄̂ − ââ̄ > 0).
16. [ASSUMPTION] l = 2, m = 0 axisymmetric even-parity initial data in all tests;
    the algorithm itself is 3D and symmetry-free, but evidence is from
    axisymmetric data.
17. [ASSUMPTION] Specific numerical configuration: 3rd-order Adams-Bashforth,
    dt = 0.001, CCE angular resolution l = 24, 12 radial points, stated domain
    partitions; convergence claims are relative to these settings.
18. [ASSUMPTION] Kerr test interprets BH mass = 1 (bare mass), spin χ = 0.5,
    spherical Kerr-Schild conformal data, excision at 1.8 inside r₊ = 1.87.
19. [ASSUMPTION] Quasinormal-mode comparison relies on BH perturbation theory
    (qnm package) as external reference.
