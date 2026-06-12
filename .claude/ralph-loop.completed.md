---
active: false
iteration: 60
iteration_base: 17
session_id: f14fd41d-d596-4dcc-9717-d73ef2b4c448
max_iterations: 1000
no_progress_limit: 8
stuck_counter_limit: 3
paper: "z4c-CCM"
completion_promise: "Z4C_CCM_LIVE_CCM_COMPLETE"
completion_status: "complete"
started_at: "2026-06-12T17:45:00Z"
---

Mission 3 (live CCM in AthenaK): (0) reference = analytic [N12 SOLID,
iter 18]; (1) constraint-preserving boundary sector in AthenaK [N13
PRELIMINARY, iter 19 — gated on O-N13-1]; (2) Cauchy and characteristic
evolutions running SIMULTANEOUSLY inside AthenaK, no postprocessing [N14].
Project note: `progress/live-ccm/RESEARCH_STATE.md`. Sources: CCM 2023 DAG
(knowledge-database/paper_arxiv-2308.10361), CCE DAG (paper_arxiv-2007.01339),
Z4c CPBC DAG (paper_arxiv-1010.0523v2), synthesis
(reformulate/z4c-CMM/synthesis_z4c-CCM/).

Each iteration: PLAN one obligation/node -> EDIT one artifact -> VERIFY with
a battery/verifier (transcript under results/numerical/) -> error-DB rows
for every failed attempt with root cause -> ledger row on admission ->
COMMIT (detailed, NO Claude attribution) + push -> update RESEARCH_STATE ->
schedule the next wake (ScheduleWakeup is the loop driver; the plugin
stop-hook promise detection is ledgered broken). Honest-admission invariant:
no gate is weakened to pass; partial results are admitted as [PRELIMINARY]
with the gap stated.

USER CONSTRAINTS (binding): NATIVE in-AthenaK implementation of the
characteristic solver (USER DIRECTIVE 2026-06-12 iter 43: "directly
implement all needed in athenak, instead of using bridge") — the C++ port
in athenak/src/z4c/ccm/bondi_solver.* is the production path, gated
against the ZccmJl-admitted numbers; ZccmJl remains the verification and
derivation harness (solver R&D stays Julia-first). FULL MIGRATION of the
Python toolchain to
Julia (packages/ZccmJl): GPU via CUDA.jl, symbolics via Symbolics.jl,
arbitrary precision via BigFloat/Rational{BigInt}; further development
ONLY in Julia once each ported component passes honest unit tests that
reproduce the admitted Python-verified numbers (no gate weakening; the
Python artifacts remain in-repo as the provenance baseline). AthenaK (C++)
is the production evolution code and stays. Plus: no run > 10 min
(timeout 600 + alarm 540);
<= 4 GPUs (mpirun -np 4 --mca pml ucx --mca osc ucx); 64-core builds with
PATH=/usr/local/cuda-12.8/bin:$PATH; commits without Claude attribution;
push to origin main.

Iteration queue (re-planned each wake; stage C parts 1+2 DONE iter 43-44:
live chain end-to-end, beta-corrected unringed solver, retarded-cone
timing, conditioning finding — see RESEARCH_STATE + n14_native_check.txt):
- iter 45: datum-FIDELITY test with real incoming content — the
  2308.10361 Sec V.C characteristic-pulse-injection run through the LIVE
  solver: ingoing l=2 J-pulse in the characteristic initial data (ZccmJl
  derivation + gates first: pulse ID, psi0 at the tube O(1) of data,
  conditioning benign), solver hands psi0 to the Cauchy side (mode 5 with
  quiescent Cauchy ID), gate vs the paper's Sec V.C response and vs a
  mode-2 prescribed-datum battery.
- iter 46+: SpECTRE CharacteristicExtract cross-oracle on identical
  worldtube data (ref-code/spectre-cce).
- iter 47+: N13 LF/GKS analysis revisit (O-N13-1); X=2 config
  (diss/kappa1 or solved ID); O(dt) coupling refinement (linear-in-time
  BC interpolation between worldtube samples).
- then: figures + method-paper refresh on corrected-solution batteries,
  N11-figure regeneration.

Exit conditions (ALL must hold):
- N13 [SOLID]: characteristic-exact constraint row stable + measurable
  absorption (ratio < 0.95 at X=0.05), or a documented LF/GKS analysis
  showing why the v1 damped row is the correct admissible realization.
- N14 [SOLID]: in-process characteristic module evolving in lockstep with
  the Cauchy evolution, worldtube data live, psi0 fed back through the
  Type-II boost (datum mode 4); verifiers: standalone linear solve matches
  the N12 exact strain/psi0 chain; live-CCM Teukolsky battery passes vs
  analytic truth; cross-oracle vs SpECTRE CharacteristicExtract on
  identical worldtube data documented.
- Figures + method paper refreshed from corrected-solution batteries.
- Ledger + RESEARCH_STATE consistent, no [BLOCKING] markers, all commits
  pushed.

Completion promise (output ONLY when unequivocally TRUE):
<promise>Z4C_CCM_LIVE_CCM_COMPLETE</promise>
