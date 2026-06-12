# result_seed.md — initial result-log entries, arXiv:1010.0523v2

Markers per `_common/contracts/markers.md`; statuses per research_admission_contract.md.
Provenance for all: ref-paper/arxiv-1010.0523v2/src/Z4.tex (stage-0 import node
S0-import-1010.0523v2) + this decomposition.

## RS1 — Equation DAG transcription
- name: eqdag-1010.0523v2; working_context: stage-1 decomposition, whole paper.
- claim: all 94 equation environments (80 tex labels + 36 auto ids = 116 nodes) transcribed into derivation.md and knowledge DB with acyclic predecessor edges.
- evidence_type: citation/transcription + automated structural check.
- verifier_result: `python3 scripts/eqdag_check.py --paper 1010.0523v2` → OVERALL: PASS (acyclicity 117 nodes incl. stage-0; edge closure; 80/80 tex labels; registry complete).
- status: checked (structural coverage only — NOT a check of physical correctness). [SOLID] for coverage; content [UNCHECKED] (O13).

## RS2 — Z4c system & constraint propagation (CL1, CL2)
- claim: conformal Z4c equations and wave-type constraint subsystem as in C1/C2 nodes.
- evidence_type: literature grounding. status: unchecked. [UNCHECKED] → close via T1/T2 symbolic derivation.

## RS3 — Characteristic structure (CL3)
- claim: characteristic variables/speeds eq:auto-11–eq:auto-15; weak hyperbolicity iff lambda = 0.
- evidence_type: literature grounding. status: unchecked. [UNCHECKED] → T3; special-case enumeration [OPEN] (O5).

## RS4 — Constraint-subsystem boundary stability, all orders L (CL4)
- claim: IBVP for Box Theta ≃ 0, Box Z_i ≃ 0 with eq:general_CPBCs is boundary stable; estimate eq:esTheta; constraints preserved for trivial data.
- evidence_type: literature grounding (paper: symbolic derivation in frozen-coefficient approximation).
- assumptions: 6, 10, 11 (assumptions.md). status: conditional (on frozen-coefficient limit, O1). [PRELIMINARY] in paper, [UNCHECKED] here → T6/T7.

## RS5 — Gauge-subsystem boundary stability, spherical reduction (CL5, CL6)
- claim: cascade property of asymptotically harmonic shift; boundary stability of lapse–shift subsystem for first and high order BCs.
- evidence_type: literature grounding. status: conditional (spherical reduction, mu_S = 1, frozen coefficients; F,G,H,J unreconstructed O4; metric sector existence_only per O6). [UNCHECKED] → T5/T8.

## RS6 — Numerical efficacy of 2nd order CPBCs (CL7, CL8, CL9)
- claim: R ≈ 1 for 1st order CPBCs; 2nd order CPBCs absorb outgoing constraint violations, keep 4th order convergence, improve constraint monitor by 2–4 orders vs Sommerfeld on flat/star/BH tests; Sommerfeld artifacts permanently alter star central density.
- evidence_type: literature grounding (paper: numerical simulation/empirical measurement).
- status: empirical (paper), unchecked (pipeline). [UNCHECKED] → T9/T10 reproduction; 3D validity [BLOCKING] per O10 if reused beyond spherical symmetry.

## RS7 — Constraint damping inadequacy at boundaries (CL10)
- claim: damping scheme cannot efficiently cure low-frequency boundary violations.
- evidence_type: conjecture / qualitative experience. status: conjectural. [HYPOTHESIS] (O14, [FUTURE]).

## RS8 — Source-tex anomalies
- claim: suspected typos: eq:auto-18 s-ring expression (O3); missing ∂_t in eq:sph_bc_last third row (O11); label `eq:theta _modes` contains a space and `sol:lf`, `sol:lf2` are nonstandard — preserved verbatim everywhere.
- evidence_type: textual inspection. status: unchecked. [OPEN] O3/O11 close via symbolic re-derivation.
