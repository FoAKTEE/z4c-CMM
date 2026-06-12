# claims.md — arXiv:1010.0523v2

Schema per `_common/contracts/research_admission_contract.md`; markers per markers.md.
All claims below are currently transcriptions of the paper's claims: evidence type
"literature grounding" unless and until re-derived/re-run by implementation tasks
(implementation_plan_python.md). Status here is the status *within our pipeline*, not the
paper's own confidence.

## CL1 — Z4c conformal system equivalence
- Working context: vacuum/matter GR, Z4 formulation, conformal variables, assumptions 1–3.
- Claim: equations eq:Z4_decomp_first, eq:auto-6, eq:auto-7 with constraints eq:Conf_Constr_1/2 are equivalent to the Z4 system eq:Z4_ADM_1–eq:Z4_ADM_2; ADM/GR recovered when Theta = Z_i = 0.
- Evidence type: literature grounding (Bernuzzi:2009ex) → target: symbolic derivation (T1).
- Dependencies: C1 nodes. Status: [UNCHECKED] (transcribed; sympy re-derivation pending).

## CL2 — Constraint subsystem is wave-like
- Working context: principal part, linear regime (assumption 9).
- Claim: Theta, Z_i, H, M_i all satisfy wave equations at light speed (eq:sys-theta-Z, eq:sys-ham-mom) given eq:ham-Z4, eq:mom-Z4.
- Evidence type: literature grounding → symbolic derivation (T2).
- Dependencies: CL1. Status: [UNCHECKED].

## CL3 — Hyperbolicity of puncture-gauge Z4c
- Working context: gauge eq:punc_alpha/eq:punc_beta, generic parameters, assumption 12.
- Claim: the second order system is strongly hyperbolic with characteristic variables eq:auto-11–eq:auto-15 and speeds ±(sqrt(mu_L), lambda, 1, sqrt(mu_S)); weakly hyperbolic iff lambda = 0.
- Evidence type: literature grounding → symbolic derivation (T3).
- Dependencies: CL1, CL2. Status: [UNCHECKED]; the "handful of special cases" is [OPEN] (O5).

## CL4 — Constraint subsystem IBVP boundary stable for all L (headline analytic result)
- Working context: frozen coefficient approximation, half-space, trivial initial data, L2 class (assumptions 6, 10, 11); four dimensional compact manifold via local freezing.
- Claim: the IBVP for eq:sys-theta-Z with CPBCs eq:general_CPBCs of any order L ≥ 0 is boundary stable (eq:fcpbcTheta) and satisfies the estimate eq:esTheta; with zero data the constraints remain zero ⇒ constraint preservation.
- Evidence type: literature grounding (paper's symbolic derivation via Kreiss theory) → targets: symbolic + numerical certification (T6, T7).
- Dependencies: CL2, C6/C7 machinery, |a_+| ≥ delta (eq:bc_est-bc). Status: [PRELIMINARY] as paper result, [UNCHECKED] in our pipeline.

## CL5 — Gauge subsystem well-posed in spherical reduction
- Working context: spherical reduction, mu_L = 2/alpha, mu_S = 1, frozen coefficients (assumptions 4–6, 13).
- Claim: the lapse–shift subsystem with first and high order BCs is boundary stable (integration constants eq:auto-25, eq:auto-29 bounded since Re(s') > 0); similar arguments cover the metric wave problems.
- Evidence type: literature grounding → symbolic + numerical certification (T8).
- Dependencies: CL3, C5, C6. Status: [UNCHECKED]; relies on unspecified shorthands F,G,H,J ([OPEN] O4); "similar arguments apply to metric components" is asserted, not shown ([OPEN] O6).

## CL6 — Cascade property of asymptotically harmonic shift
- Working context: frozen coefficients, (eps_alpha, eps_chi) = (1, 1/2) (assumption 5).
- Claim: with the asymptotically harmonic shift the gauge sector couples to the metric only through the BCs, enabling sequential (cascade) analysis; with (0,0) everything couples in the bulk.
- Evidence type: literature grounding → symbolic derivation (T5).
- Dependencies: eq:full-sec-z4-alpha/-beta/-metric. Status: [UNCHECKED].

## CL7 — 2nd order CPBCs absorb constraint violation; 1st order do not
- Working context: spherical numerics, assumptions 13–16; perturbed flat spacetime test.
- Claim: reflection coefficient R ≈ 1 for L = 0 CPBCs (pure reflection); L = 1 CPBCs absorb the largest portion immediately, behave qualitatively like Sommerfeld in R(k) but constraint-preservingly, and keep 4th order convergence of reflections (Sommerfeld reflections are 1st order accurate).
- Evidence type: literature grounding (paper's numerical simulation) → numerical simulation (T9, T10).
- Dependencies: eq:refcoef, eq:theta _modes, eq:Cmonitor, C10. Status: [UNCHECKED] (empirical in paper).

## CL8 — CPBCs necessary for physical reliability (star test)
- Working context: stable TOV-type star, Z4c + GR hydro, close boundary r_out ≈ 20 M.
- Claim: Sommerfeld BCs permanently shift the star's central density (non-convergent artifact); 2nd order CPBCs remove the artifact and bring the constraint monitor within 2–4 orders of magnitude of the reference solution.
- Evidence type: literature grounding (numerical) → numerical simulation (T10, qualitative).
- Dependencies: CL7. Status: [UNCHECKED] (empirical in paper).

## CL9 — Robustness on black hole spacetimes; asymptotically harmonic shift on punctures
- Working context: single puncture (trumpet endpoint) and excised Kerr–Schild data (eq:auto-30).
- Claim: CPBCs perform comparably across gauges and data; Sommerfeld with close boundary crashes the Kerr–Schild run; first use of the asymptotically harmonic shift on puncture data.
- Evidence type: literature grounding (numerical). Status: [UNCHECKED]; full reproduction [FUTURE].

## CL10 — Constraint damping is not a substitute for CPBCs
- Working context: boundary-induced violations, damping scheme of Gundlach et al.
- Claim: damping coefficients required to suppress boundary perturbations are large because such perturbations are low frequency; indiscriminate damping acts like artificial dissipation.
- Evidence type: literature grounding (qualitative numerical experience). Status: [CONJECTURAL]/[UNCHECKED]; no quantitative data shown.
