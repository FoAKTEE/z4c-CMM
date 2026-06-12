# obligations.md — arXiv:1010.0523v2

Open obligations with expected evidence type to close each. Markers per markers.md.

- O1 [OPEN] Frozen-coefficient gap: well-posedness shown only in the high-frequency limit; transfer to the variable-coefficient/nonlinear IBVP is "expected", not proven. Close by: pseudo-differential symmetrizer construction (paper cites Taylor99b) or alternative proof. Owner: future analytic work (paper's own stated gap).
- O2 [OPEN] Standard puncture gauge (eps_alpha = eps_chi = 0): no cascade, no symmetrizer, no well-posedness proof; only numerical evidence. Close by: new mathematical approach (paper acknowledges); for our pipeline, numerical mode analysis (T8 extension).
- O3 [OPEN] Sign/consistency check of the frozen normal vectors eq:auto-18: tex prints s̊^a ∂_a = −β̊ ∂_x, which has no unit spatial part; expected s̊^a ∂_a ∝ −∂_x (s^i = −ê_x per Sec. III.B). Close by: symbolic re-derivation in T5; suspected typo in source tex.
- O4 [OPEN] Shorthands F(β̊), G(β̊), H(β̊), J(β̊) in eq:auto-26/eq:auto-27/eq:auto-29 never given explicitly; also γ_a vs γ_α and a^{L+1} vs a_+^{L+1} notational drift. Close by: reconstructing them via the A^{L+1} iteration in T8 (symbolic derivation).
- O5 [OPEN] "Handful of special cases" where puncture-gauge Z4c fails strong hyperbolicity is not enumerated. Close by: characteristic analysis over parameter space (T3).
- O6 [OPEN] "Similar arguments apply to the wave problem of the metric components" (end of Sec. IV.B) — the metric-sector boundary stability (coupled eq:Z4_gammasbc_2+1 / eq:Z4_gammaAbc_2+1 problem) is asserted, not derived. Close by: explicit 2x2 coupled boundary-determinant computation (T8 extension).
- O7 [OPEN] Number of incoming characteristics (ten) assumed from near-flatness (assumption 8); no quantitative closeness criterion. Close by: eigenvalue monitoring in implementation; dimensional/frame check of boundary operators.
- O8 [OPEN] Absorption properties of the gauge/metric conditions eq:BCs-alpha–eq:BCs_lastII and their relation to incoming gravitational radiation: explicitly deferred by the paper. [FUTURE] for our scope.
- O9 [OPEN] Discrete stability: no SBP/energy proof for the implementation; stability supported only by toy problems and long-run tests. Close by: numerical stability experiments (T9) or semi-discrete analysis. 
- O10 [OPEN] 3D validity: numerical results are spherically symmetric; BAM 3D tests indicate discarded tangential terms in the BCs are required for stability. [BLOCKING] for any 3D reuse of eq:sph_bc_1st–eq:sph_bc_last; unblocked by adding tangential terms (paper's follow-up work). 
- O11 [OPEN] eq:sph_bc_last third equation has Ã_rr (not ∂_t Ã_rr) on the left-hand side — dimensional/frame check needed (suspected missing ∂_t in source tex). Close by: re-derivation from T4 operators (symbolic).
- O12 [OPEN] Approximation remainders: linearization around flat space in Appendix C and around the background in eq:general_CPBCs carry uncontrolled remainders for strong-field boundaries. Close by: controlled approximation check in T9 (vary r_out).
- O13 [UNCHECKED] All 116 DAG node transcriptions are literature-grounded; none re-derived yet. Close by: T1–T8 symbolic verifications, promoting nodes per result_seed.md.
- O14 [OPEN] Constraint damping interaction (CL10) unquantified. Close by: numerical simulation with damping terms restored. [FUTURE].
