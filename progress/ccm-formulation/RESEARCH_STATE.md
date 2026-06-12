# z4c-CCM — New-Formulation Mission — Verified Research State

Mission 2: derive a NEW FORMULATION of Cauchy-characteristic matching for Z4c,
combining arXiv:1010.0523v2 (Z4c CPBC), 2007.01339 (CCE), 2308.10361 (GH-CCM).
Organizing principle: physical dof + characteristic modes + boundary conditions.
End product constraint: ONE DAG (single connected acyclic graph) — enforced by
the verifier's connected-component check.
Phase: mission complete — N1–N6 solid, N7 sketch solid, N8 final product solid;
N9 remains a documented [FUTURE] obligation (not in the completion gate).
Branch: `main`. Loop state: `.claude/ralph-loop.local.md` (mission z4c-CCM;
terminal archive: `progress/ccm-formulation/loop_state_terminal.md`).

## Source Library

| ID | Source | Kind | Status | Notes |
|---|---|---|---|---|
| 1010.0523v2 | ref-paper/arxiv-1010.0523v2/ | paper | [SOLID] decomposed (mission 1) | Z4c CPBC: 10-mode split, CPBCs, gauge BCs, LF/Kreiss machinery |
| 2007.01339 | ref-paper/arxiv-2007.01339/ | paper | [SOLID] decomposed (mission 1) | CCE engine: worldtube data → hierarchy → psi0/news |
| 2308.10361 | ref-paper/arxiv-2308.10361/ | paper | [SOLID] decomposed (mission 1) | GH-CCM: psi0 → Type-II boost → w−|_BC Bjørhus injection |

## Working Context

The single formulation DAG: `reformulate/z4c-CMM/synthesis_z4c-CCM/formulation_dag.md`
(36 nodes, 53 edges; sectors A/B = P1, C = P2, D/E = P3, N = new work program).
Rendered figure: `results/figures/zccm_formulation_dag.html`
(generator `scripts/render_zccm_dag.py`).

## Active Claims

(none — the mission claim is discharged: the Z4c-CCM BC set is explicit,
verified, and implemented; see Accepted Results Log. Successor work program:
[FUTURE] items in Next Work Steps.)

## Binding user constraints (2026-06-12)

1. No test run > 10 minutes (`timeout 600` + internal `signal.alarm(540)`).
2. At most 4 GPUs (`CUDA_VISIBLE_DEVICES=0,1,2,3`; host has 8× A100-80GB).
3. FINAL PRODUCT = GPU-runnable code: JAX float64 Z4c-CCM implementation under
   `packages/zccm/`, multi-GPU capable (≤4 devices). N8 promoted into the
   completion gate.

## Accepted Results Log

| Claim | Evidence type | Evidence / verifier | Assumptions / deps | Status | Open obligations |
|---|---|---|---|---|---|
| ONE connected acyclic formulation DAG combining the three papers, every source node grounded in the per-paper equation ledgers | automated check | `python3 scripts/zccm_dag_check.py` → OVERALL: PASS; results/numerical/zccm_dag_check.txt | grounding = label resolution in mission-1 ledgers (transcription-level, not re-derived) | [SOLID] | N-node derivations remain (below) |
| N2: characteristic l̂ ∝ Cauchy (n+s)/√2 with factor √2/2·e^{2β̂}/(α−γ_ij β^i s^j) for fully general 3+1 data ⇒ paper-3 Type-II boost valid verbatim on the Z4c frame | exact rational-point identity testing + 4-GPU float64 sweep | `scripts/verify_n2_boost.py` → results/numerical/n2_boost_check.txt (OVERALL: PASS, 17.6 s) | worldtube at fixed coordinate radius; u=t on tube; outgoing root | [SOLID] | none |
| N3 EXACT (user-review upgrade, iter 16): (E + ε(s)B)(m,m) = −ψ₀ and (E − ε(s)B)(m,m) = −ψ̄₄ as pure Weyl algebra (generic 10-dof tensor); U⁻_AB = TT₂[E + ε(s)B]_AB = −(ψ₀m̄m̄+c.c.); Z4c realization E = R + KK − KK, B = +εDK verified NONPERTURBATIVELY on Schwarzschild (route agreement 7.7e-17 / 3.9e-17, identity 8.5e-16); linearized form demoted to on-shell flat corollary | exact generic-Weyl algebra + autodiff Schwarzschild two-route + 4-GPU | `scripts/verify_n3_exact.py` → results/numerical/n3_exact_check.txt (OVERALL: PASS, 21 s); corollary: n3_dictionary_check.txt | on-shell vacuum modulo Z4c constraint fields (CPBC-controlled, N4) | [SOLID] (exact) | nonlinear constraint-damping term (N9) |

| N1: Z4c ↔ ADM ↔ 4-metric maps exact; kinematic closure forces D_iβ^i = full covariant divergence (convention resolved, c=1 residual for naive reading); `packages/zccm` z4c_vars module (with N2 boost + N3 datum constructors) 4-GPU tested | exact full-symbolic (rational) + 4-GPU float64 package test | `scripts/verify_n1_varmap.py` + `packages/zccm/tests/test_z4c_vars.py` → results/numerical/n1_varmap_{check,gpu_test}.txt (both OVERALL: PASS) | det γ̃ = 1, tr_γ̃ Ã = 0 constraints; spatial derivatives supplied by evolution code | [SOLID] | none |
| N4: CCM/TT channel is constraint-orthogonal (H = M_i = 0 off-shell) ⇒ P1 order-L CPBCs retained verbatim; boundary-channel taxonomy TT/gauge(h_xx)/vector(h_xA)/trace derived with negative controls | exact symbolic + negative controls + 4-GPU sweep | `scripts/verify_n4_cpbc_compat.py` → results/numerical/n4_cpbc_compat_check.txt (PASS, 12.4 s) | linearized flat boundary frame; nonlinear coupling + 3D tangential terms + corners → obligations | [SOLID] | nonlinear coupling, corner compatibility (N7) |
| N5: gauge sector decouples from the CCM channel (trace and Γ̃^i sources vanish off-shell), and on P1's eq:eomgammasA model the boundary-adapted operator has R = 0 vs Sommerfeld R = −b(1−b)/((1+b)(2−b)) — P3's bottleneck removed at model level | exact symbolic + exact reflection algebra + 4-GPU sweep | `scripts/verify_n5_gauge.py` → results/numerical/n5_gauge_check.txt (PASS, 12.2 s) | normal incidence model; oblique/curvature → N7 | [SOLID] | oblique incidence + curvature (N7); numerical R tests (N8) |
| N6: composite scheme — explicit 10-mode BC table (z4c_ccm_boundary_conditions.md) + zccm.boundary GPU module; mode count exact, ψ₀→0 reduces to P1 CPBC scheme, N1+N2+N3 wire-through 2.9e-16 over 1M states, outgoing-transparent / incoming-constrained exactly | structural + 4-GPU float64 | `scripts/verify_n6_composite.py` → results/numerical/n6_composite_check.txt (PASS, 22.2 s) | linearized flat boundary frame for the dictionary rows | [SOLID] | O-N6-1..4 routed to N7/N8/N9 |
| N7: well-posedness sketch — CCM changes data not operator; D_L(s) = (2s)^{L+1} uniform Kreiss at all L, datum-independent; P1 stability transfers verbatim to inhomogeneous rows; full proof program documented as [FUTURE] obligations (coupled feedback, oblique/curvature, corners, nonlinear) | exact LF algebra + 4-GPU sweep | `scripts/verify_n7_lf_sketch.py` → results/numerical/n7_lf_sketch_check.txt (PASS, 13.5 s) | frozen coefficients, normal-incidence model, ψ₀′ as given data | [SOLID] (sketch scope) | full proof [FUTURE] — human sign-off requested (sketch doc items 1–4) |

| N8 (FINAL PRODUCT): packages/zccm complete — z4c_vars + boundary + model1d; dynamical verification: Sommerfeld reflects at the analytic coefficient (0.03% dev), P1/CCM operator absorbs at 3.2e-15, CCM injection two-way transparent to 1.1e-4 | 4-GPU dynamical evolution + package tests | `packages/zccm/tests/test_model1d.py` → results/numerical/n8_model1d_test.txt (PASS, 18.8 s) | model problem (P1 eq:eomgammasA); full-GR 3D tests are successor work | [SOLID] (model-level) | full-GR 3D test campaign [FUTURE] |

## Next Work Steps (successor missions — all [FUTURE])

- `[FUTURE] N9` — κ₁/κ₂ damping vs injected-data error (Z4c analog of γ₂c³).
- `[FUTURE]` N7 full-proof program: coupled Cauchy↔characteristic feedback (weak hyperbolicity), oblique/curvature, corners (sketch doc items 1–4; human sign-off requested).
- `[FUTURE]` full-GR 3D tests: Teukolsky wave, Kerr perturbation, pulse injection with a production Z4c code + SpECTRE-style CCE.

## Appendix

### Abandoned Methods

(none)

### Audit references

- z4c-CCM ledger: `knowledge-database/paper_z4c-CCM/nodes.jsonl` (N1–N9 + S2-synthesis-dag).
- Mission-1 equation DAGs: `knowledge-database/paper_arxiv-<id>/nodes.jsonl` (325 equation nodes).
- alignment.md anomaly (see mission-1 note): upstream commit 663804c replaced the protocol; binding rules taken from the settings.json wrapper.
