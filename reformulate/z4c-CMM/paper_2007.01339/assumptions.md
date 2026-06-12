# Assumptions — arXiv:2007.01339

Enumerated symmetries, limits, and regimes used by the paper's calculations. Every item is `[ASSUMPTION]` unless noted.

1. `[ASSUMPTION]` **Asymptotic flatness.** The spacetime admits a Bondi-Sachs frame with the Minkowski falloffs `eq:BondiSachsFalloffs`; scri+ exists and is reachable by the compactified domain.
2. `[ASSUMPTION]` **Vacuum exterior.** Outside the worldtube the Einstein field equations hold in vacuum; Section VI explicitly sets R_{μν} = 0 for the Weyl-scalar identities (eq:auto-44).
3. `[ASSUMPTION]` **Bondi-like gauge attainability.** The Cauchy worldtube data admits the local gauge conditions g_rr = g_rA = 0, det(h_AB) = det(q_AB), and is "sufficiently well behaved" that metric components are asymptotically finite (`eq:BondiLikeFalloff` is not actually enforced by the standard algorithm — paper's own caveat).
4. `[ASSUMPTION]` **Null foliation regularity.** Outgoing null hypersurfaces emanating from the worldtube reach scri+ without caustics or crossings in the exterior region (implicit in the characteristic method; worldtube at ~100–1000 Schwarzschild radii).
5. `[ASSUMPTION]` **Regular (polynomial in 1/r) asymptotic expansions.** All spin-weighted scalars are assumed expandable as F = F^(0) + l F^(1) + l² F^(2) + …; the regularity analysis (Section V) shows closure of this ansatz in partially flat gauge but assumes no half-integer or log terms in the input data.
6. `[ASSUMPTION]` **Partially flat gauge conditions.** Section V formally imposes Ĵ^(0) = Û^(0) = 0; all regularity conclusions are conditional on these plus Ĵ^(2) = 0 on the initial hypersurface.
7. `[ASSUMPTION]` **Initial-data falloff for Ĵ.** Initialization prescribes Ĵ with pure r^{-1} falloff at scri+ — terms with r^0 or r^{-2} dependence are forbidden (computational roadmap, initialization step 2).
8. `[ASSUMPTION]` **Perturbative solvability near scri+.** eq:uRequation and eq:bsxA are solved as power series in l̂; only the first orders are needed for asymptotic quantities. No claim of convergence or numerical stability of a non-perturbative radial integration (paper explicitly disclaims this).
9. `[ASSUMPTION]` **Dyad and normalization conventions.** Complex dyad q^A q̄_A = 2 with q^A = {−1, −i/sinθ}; spin-weighted transform libraries capped at |s| ≤ 2 suffice (the equations were arranged to respect this).
10. `[ASSUMPTION]` **BMS residual freedom.** All "gauge-invariant" statements (news, Weyl scalars) are invariant only up to BMS transformations fixed by the metric on the initial CCE hypersurface; tetrad ambiguity additionally affects Weyl scalars.
11. `[ASSUMPTION]` **First-derivative worldtube data.** The Cauchy code supplies metric components and their first derivatives on the worldtube; identities eq:auto-14/eq:auto-15 (using one Einstein equation component) remove any need for second derivatives.
12. `[ASSUMPTION]` **Smooth worldtube radius function.** R(u, θ, φ) = r|_Γ is smooth, and the worldtube remains a y̆ = −1 surface (not constant Bondi-like r) throughout.
13. `[ASSUMPTION]` **Hierarchy well-posedness.** On each hypersurface the radial ODE hierarchy eq:BondiHierarchy has unique solutions given worldtube boundary values; u-stepping is treated as a standard ODE problem for J.
14. `[ASSUMPTION]` **Peeling.** Weyl scalars obey Ψ_n ~ r^{n−5} (eq:auto-49), as required for the asymptotic Weyl-scalar formulas.
15. `[ASSUMPTION]` **Initial angular-coordinate coincidence.** x̂^Â(u=0) = x^A, hence ω̂(u=0) = 1 — used to simplify initialization and the regularity argument.

`[FUTURE]` The paper does not treat matter at the worldtube, caustic formation, or Cauchy-characteristic matching (CCM); these are explicitly out of scope.
