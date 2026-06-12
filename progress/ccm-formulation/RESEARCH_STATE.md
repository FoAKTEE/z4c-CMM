# z4c-CCM — New-Formulation Mission — Verified Research State

Mission 2: derive a NEW FORMULATION of Cauchy-characteristic matching for Z4c,
combining arXiv:1010.0523v2 (Z4c CPBC), 2007.01339 (CCE), 2308.10361 (GH-CCM).
Organizing principle: physical dof + characteristic modes + boundary conditions.
End product constraint: ONE DAG (single connected acyclic graph) — enforced by
the verifier's connected-component check.
Phase: roadmap admitted (S2-synthesis-dag [SOLID]); N-node derivation work next.
Branch: `main`. Loop state: `.claude/ralph-loop.local.md` (mission z4c-CCM, iteration 7).

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

| Claim | Needed evidence | Priority | Owner |
|---|---|---|---|
| Z4c-CCM BC set: 2 physical ← CCE ψ0 (Type-II boost), 4 constraint ← CPBC, 4 gauge ← P1 gauge BCs is consistent and explicitly writable | N2→N3→N6 derivations, symbolically verified | [PRELIMINARY] | loop agent |

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
| N3 (linearized core): (∂_t+∂_s)² γ_AB^TF = 4(ψ₀ m̄_A m̄_B + c.c.), ψ₀ = ¼(∂_t+∂_s)² h_mm = −C(l,m,l,m), OFF-SHELL on the P1 frame ⇒ CCE ψ₀′ (N2 boost) supplies the h_AB^TF datum of eq:BCs_lastII; L=0 Bjørhus form fixed | exact symbolic (arbitrary functions) + 4-GPU off-shell Fourier sweep | `scripts/verify_n3_dictionary.py` → results/numerical/n3_dictionary_check.txt (OVERALL: PASS, 20.1 s) | linear TT perturbation on flat boundary frame; nonlinear/conformal corrections → N6 obligations | [SOLID] | nonlinear corrections, constraint-term replacement (N6) |

## Next Work Steps

- `[PRELIMINARY] N1` — Z4c → worldtube 4-metric variable map; verifier: round-trip identity g(α,β,χ,γ̃); seeds `packages/zccm/` JAX module.
- `[PRELIMINARY] N4` — CPBC retention argument written out with corner-compatibility obligation.
- `[HYPOTHESIS] N5` — gauge-sector replacement statement + reflection analysis plan.
- `[HYPOTHESIS] N6` — composite scheme assembly (depends N3, N4, N5).
- `[FUTURE] N7/N8/N9` — well-posedness sketch, tests, damping interplay.

## Appendix

### Abandoned Methods

(none)

### Audit references

- z4c-CCM ledger: `knowledge-database/paper_z4c-CCM/nodes.jsonl` (N1–N9 + S2-synthesis-dag).
- Mission-1 equation DAGs: `knowledge-database/paper_arxiv-<id>/nodes.jsonl` (325 equation nodes).
- alignment.md anomaly (see mission-1 note): upstream commit 663804c replaced the protocol; binding rules taken from the settings.json wrapper.
