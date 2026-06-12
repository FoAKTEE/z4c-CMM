# Iteration 16 — N3 upgraded to EXACT after user review; docs/package/paper retargeted

1. **Paper anchor** — DAG node N3 (exact scope); P3 eq:wab pattern made Z4c-native; trigger: user review of paper eqs. 16-17 ("still perturbative — use the full precise calculation").
2. **What shipped** — `scripts/verify_n3_exact.py` (generic-Weyl exact identity with derived constants (σ,κ)=(1,−1); ψ̄₄ partner; dyad completeness; linearized on-shell limit recovered; Kerr-Schild Schwarzschild two-route autodiff check on 4 GPUs); `zccm.psi0_to_u_minus` + `zccm.u_minus_scalar` (exact constructors); `zccm.boundary` retargeted to U⁻; N6 verifier + package test updated and re-PASS; BC doc + DAG doc + RESEARCH_STATE updated; paper §II.D/§II.E rewritten exact (eqs. 16-19 new), recompiled clean (4 pp). Pushed through `58aab4c` earlier this iteration.
3. **Next-3 roadmap** — N9 (nonlinear constraint-damping addition — now the natural next node since the exact object is in place), oblique-incidence analysis (N7 obligation 2), full-GR worldtube test with E/B from a real Z4c state including spatial derivatives.
4. **Simplification flag** — psi0_to_physical_datum retained as documented linearized-limit helper (model1d uses the flat-frame operator form); no duplication.
5. **Verifier outputs** (certified): n3_exact_check.txt — 9/9 PASS, wall 21 s:
   generic Weyl 10 dof ✓ · E,B sym/traceless ✓ · EXACT identity (σ,κ)=(1,−1) ✓ · ψ̄₄ partner (−1) ✓ · dyad completeness ✓ · linearized on-shell limit (σ_B=+1) ✓ · Schwarzschild E 7.68e-17, B 3.93e-17, identity 8.46e-16 (1024 pts) ✓
   n6_composite_check.txt — 5/5 PASS (exact targets, wire-through 1.47e-16); n1_varmap_gpu_test.txt — 4/4 PASS (incl. exact κ=−1 constructor).
