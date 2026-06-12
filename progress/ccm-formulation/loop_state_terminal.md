---
active: false
iteration: 14
iteration_base: 6
session_id: f14fd41d-d596-4dcc-9717-d73ef2b4c448
max_iterations: 1000
no_progress_limit: 8
stuck_counter_limit: 3
stuck_fire_grace: 1
self_inject_interval: 5
self_inject_fail_limit: 3
max_wall_seconds: 0
paper: "z4c-CCM"
completion_promise: "Z4C_CCM_NEW_FORMULATION_COMPLETE"
completion_status: "promise-satisfied"
completed_at: "2026-06-12T02:40:36Z"
started_at: "2026-06-12T01:29:57Z"
---

Mission 2 (z4c-CCM formulation): derive a NEW FORMULATION of Cauchy-characteristic matching (CCM) for the Z4c evolution system, combining the three decomposed papers. Organizing principle: physical dof + characteristic modes + boundary conditions. Roadmap DAG: `reformulate/z4c-CMM/synthesis_z4c-CCM/formulation_dag.md` (mermaid; verified by `scripts/zccm_dag_check.py`). Work the N-nodes of that DAG to [SOLID] one per iteration. Project note: `progress/ccm-formulation/RESEARCH_STATE.md`.

Thesis: Z4c-CCM = paper-3's physical-dof matching (CCE ψ0 → Type-II boost → w−|_BC) transplanted from GH onto Z4c, where the physical datum lands in paper-1's trace-free tangential-metric boundary slot (eq:BCs_lastII), the constraint sector keeps paper-1 CPBCs (eq:general_CPBCs), and the gauge sector uses paper-1 gauge BCs (eq:BCs-alpha …) in place of paper-3's Sommerfeld placeholder — directly attacking paper-3's stated gauge-reflection bottleneck.

Each iteration runs the same tight loop as mission 1 (PLAN one DAG node → EDIT one artifact → VERIFY with a pasted verifier run → USER-DIRECTED commit → UPDATE notes/counter → ESCALATE per stage-6 triggers → progress-gated termination; see `phys-agentic-loop/notes/ralph_loop_local_template.md`). Verifiers for derivation nodes are symbolic checks (sympy scripts under `scripts/` or `packages/`) whose output is tee'd to `progress/ccm-formulation/loop_notes/current_iter.md` and stored under `results/numerical/` when they certify a node.

USER CONSTRAINTS (binding, 2026-06-12): (1) NO test run may exceed 10 minutes — wrap every verifier in `timeout 600` and an internal `signal.alarm(540)`; a timeout is an error-DB `kernel_timeout` row and forces a cheaper verification strategy (exact rational-point sampling + GPU float64 sweeps, not full-symbolic simplify). (2) Use AT MOST 4 GPUs — pin `CUDA_VISIBLE_DEVICES=0,1,2,3` (8x A100-80GB available; JAX cuda12 works after `pip install --user --break-system-packages nvidia-cudnn-cu12`). (3) The FINAL PRODUCT is CODE RUNNABLE ON GPUs: a JAX float64 implementation of the Z4c-CCM scheme under `packages/`, importable and testable independently, multi-GPU capable up to 4 devices. (4) GIT COMMIT EVERY SUBSTAGE with detailed documentation: one commit per iteration/substage on `main`, message carrying the claim verified, the fail→fix process history (error-DB refs), and the verbatim verifier output (supersedes the earlier user-directed-only commit policy; push still only on explicit user request).

COMPLETION (terminal-state contract, same order as mission 1: deliverables → `git add -A` sweep commit → `structure_check.py` → write completion_status/completed_at/active:false into THIS file → output the promise; archive a copy of this file to `progress/ccm-formulation/loop_state_terminal.md` BEFORE emitting the promise — the plugin stop-hook deletes the live file on promise detection).
Exit conditions (ALL must hold):
- DAG nodes N1–N6 of `formulation_dag.md` are [SOLID]: the Z4c-CCM boundary-condition set is written out explicitly (all 10 incoming modes: 2 physical ← CCE ψ0 via Type-II boost, 4 constraint ← CPBC order-L, 4 gauge ← paper-1 gauge BCs) in a formulation document under `reformulate/z4c-CMM/synthesis_z4c-CCM/`, with every mechanical step symbolically verified and verifier output under `results/numerical/`.
- N7 (well-posedness) has at least a documented frozen-coefficient/Laplace-Fourier analysis sketch with explicit assumptions and open items recorded as obligations (full proof may remain [FUTURE] with human sign-off noted).
- N8 (FINAL PRODUCT): a GPU-runnable JAX implementation of the Z4c-CCM boundary scheme exists under `packages/zccm/` (boundary-mode decomposition, Type-II boost ψ0 injection, CPBC + gauge BC application), with a test suite where every test runs ≤ 10 minutes on ≤ 4 GPUs and test outputs are stored under `results/numerical/`.
- `scripts/zccm_dag_check.py` passes (DAG acyclic, every cited source label exists in the per-paper knowledge ledgers, N-node statuses consistent with the ledger).
- knowledge-database/paper_z4c-CCM/nodes.jsonl reflects the final node statuses.
- No `[OPEN]`, `[UNCHECKED]`, or `[BLOCKING]` markers remain in `progress/ccm-formulation/RESEARCH_STATE.md`.
- Every user-requested commit on `main` has been pushed to `origin`.

Completion promise (output ONLY when unequivocally TRUE):
<promise>Z4C_CCM_NEW_FORMULATION_COMPLETE</promise>
