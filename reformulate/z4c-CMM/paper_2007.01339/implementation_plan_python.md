# Implementation plan (Python) — arXiv:2007.01339

Partition of work per logic.md node clusters C1–C6. Target stack: `sympy` for symbolic verification (primary stage-2 goal), `numpy`/`scipy` + a spin-weighted spherical-harmonic library (e.g. `spinsfast`/`spherical`) for numerical prototypes. Each module closes specific obligations from obligations.md.

## M1 — Spin-weighted calculus core (cluster C1)
Nodes: eq:auto-4, eq:ComplexDyadChoice, eq:DyadRotation, eq:auto-5, eq:auto-6, eq:auto-7, eq:QDef, WDef, eq:KDef.
- `swsh.py`: dyad algebra, spin-weight tracking type, ð/ð̄ operators in coordinate form (eq:auto-7) and spectral form; K(J) helper with the stable 1−K = JJ̄/(1+K) evaluation.
- Symbolic twin in sympy: dyad contraction utilities used by all later verification modules.
- Tests: ð on spherical harmonics (known eigenvalues); spin-weight closure |s| ≤ 2 (O-8).

## M2 — Worldtube boundary extraction (cluster C2)
Nodes: eq:auto-8 … eq:auto-16, eq:BondiLikeRadius, eq:auto-13 group.
- `boundary.py`: ADM → null-radius metric (eq:auto-12), Bondi-like radius (eq:BondiLikeRadius), Bondi-like scalars (eq:auto-13), first-derivative identities (eq:auto-14, eq:auto-15), worldtube H (eq:auto-16).
- Verification: sympy check that eq:auto-13 follows from the up-index metric; analytic Schwarzschild-in-Kerr-Schild worldtube data as regression fixture.

## M3 — Gauge transformations to partially flat / Bondi-Sachs (cluster C3)
Nodes: eq:UpIndexBondiLikeMetric … eq:auto-23, eq:ifcQ, eq:auto-52, eq:ifcH, eq:auto-53 … eq:auto-55.
- `gauge.py`: angular-coordinate evolution (eq:XIfcEom, eq:auto-38), Jacobians â, b̂ (eq:auto-19), ω̂ (eq:ifc_omega), scalar transformations (eq:ifcBeta group, eq:ifcQ, eq:auto-52, eq:ifcH), ∂_ûω̂ (eq:auto-20).
- `bondi_sachs.py` (symbolic only): verify eq:uRequation, eq:bsxA from the transformation ansatz; reproduce Appendix C orders (eq:auto-53 … eq:auto-55) and eq:auto-23 — closes O-2.

## M4 — Compactified hypersurface solvers (cluster C4)
Nodes: eq:auto-24 … eq:auto-31, eq:form1, eq:form2, eq:form3, eq:auto-26, all eq:Hypersurface* and *numeric equations.
- `equations.py`: source bundles Λ_Q̆ (eq:auto-28), Λ_W̆, 𝒜/ℬ/𝒞/𝒟/ℒ_H̆ (eq:auto-30, eq:auto-31) with conjugate-caching as the paper recommends.
- `radial_solvers.py`: form-1 direct integration, form-2/form-3 regular singular solvers on y̆ ∈ [−1,1] (Legendre/Chebyshev collocation), matrix real/imag decomposition (eq:auto-26).
- Symbolic verification: derive each printed equation from the Einstein tensor of eq:BondiLikeMetric (sympy; heavy — may stage through the companion Mathematica package, O-1, O-7). Reconcile eq:KDef sign (O-6).

## M5 — Regularity and evolution loop (cluster C5)
Nodes: eq:QWabstract … eq:Hregularity, eq:auto-32 … eq:auto-38, eq:UHat, eq:SubleadingQRegularity, eq:auto-33.
- `regularity.py` (symbolic): order-by-order expansion of each hypersurface equation under Ĵ^(0)=Û^(0)=0, Ĵ^(2)=0; confirm eq:auto-32 … eq:auto-36 — closes O-3.
- `evolution.py`: the roadmap loop — boundary transforms, hierarchy pass with 𝒰 trick (eq:UHat, eq:auto-37), J step via H̆, angular-coordinate step; RK4/Dormand-Prince stepping.
- Tests: linearized solutions (Teukolsky waves) and log-detection diagnostic (fit residual r^{-n} ln r component must converge to zero spectrally).

## M6 — News and Weyl outputs (cluster C6)
Nodes: eq:BondiNews, eq:auto-17, eq:NewsDefinition*, eq:auto-39 … eq:auto-51, eq:NPCovariantDerivatives.
- `scri.py`: news via eq:NewsDefinitionIncompletelyFlat; optional arbitrary-gauge news eq:NewsDefinitionBondiLike; asymptotic Weyl scalars eq:auto-50, Bondi-frame conversion eq:auto-51.
- Symbolic verification of spin coefficients eq:auto-43 and Weyl identities eq:auto-44 → eq:auto-47/eq:auto-48 from the tetrad eq:auto-39 — closes O-4, O-5.

## Sequencing and gating
1. M1 → M2/M3(symbolic) → M4(symbolic) → M5(symbolic) — verification chain; each module's CAS outputs are the admission evidence promoting claims C-1 … C-5 from `[UNCHECKED]`.
2. Numerical prototypes (M4/M5/M6 numeric parts) only after the corresponding symbolic check passes — avoids implementing transcription errors (C-7).
3. All verification scripts emit machine-readable pass/fail; results feed result_seed.md → stage-3 validation and the knowledge DB (status promotion solid only after verifier output is pasted).
