# result_seed.md — initial result-log entries for arXiv:2308.10361

Markers per `_common/contracts/markers.md`; statuses per the research admission
contract. These seed entries are inputs to pipelines/4-result-log.

## R1 — Stage-1 decomposition complete
- id: r-2308.10361-decomp-1
- working_context: stage-1 decomposition of arXiv:2308.10361 into
  reformulate/z4c-CMM/paper_2308.10361/; equation DAG mirrored in
  knowledge-database/paper_arxiv-2308.10361/nodes.jsonl
  (task_id eq-dag-stage1-2308.10361).
- claim: all 106 equation environments (74 tex labels + 38 auto labels) are
  transcribed in derivation.md and registered as 106 acyclic DAG nodes.
- evidence_type: verifier output.
- evidence: `python3 scripts/eqdag_check.py --paper 2308.10361` →
  acyclicity PASS (107 nodes incl. stage-0 import), edge closure PASS,
  74/74 tex labels covered, all DB labels registered in derivation.md,
  OVERALL: PASS (run 2026-06-12).
- status: checked (for the registration claim only). [SOLID]

## R2 — Paper equation content
- id: r-2308.10361-content-1
- working_context: as R1.
- claim: the equations and dependency edges recorded in derivation.md/logic.md
  faithfully transcribe the paper's stated mathematical structure.
- evidence_type: literature grounding (citation/transcription).
- status: unchecked. [UNCHECKED] — transcription is NOT verified derivation;
  upgrades via implementation_plan_python.md T1-T8.

## R3 — CCM exactness (paper claim C1)
- id: r-2308.10361-ccm-exact
- claim: ψ₀ → w−|_BC mapping is approximation-free in the continuum.
- evidence_type required: symbolic derivation. status: unchecked. [OPEN]
  obligations 1-4.

## R4 — Empirical test results (paper claims C6-C12)
- id: r-2308.10361-tests
- claim set: stability, unchanged GH constraints, reduced Bondi-gauge violations,
  ~10x waveform-deviation reduction at X = 2, transparent boundary for injected
  characteristic pulses.
- evidence_type: numerical simulation (external; data at CCMData GitHub repo).
- status: empirical, unreproduced here. [UNCHECKED] — obligation 6.

## R5 — Well-posedness caveat
- id: r-2308.10361-wellposedness
- claim: CCM as implemented is not well-posed (weakly hyperbolic characteristic
  system); stability evidence is empirical and limited to smooth data.
- evidence_type: literature grounding (Giannakopoulos et al.).
- status: conditional. [ASSUMPTION] carried by all downstream results;
  obligation 7 [FUTURE]/[BLOCKING] for continuum-convergence claims.

## R6 — Typo correction usable downstream
- id: r-2308.10361-K-typo
- claim: K = √(1+JJ̄) (eq:CCE_all_variables:K) corrects Moxon:2020gha Eq. (10e).
- evidence_type required: symbolic derivation (short). status: unchecked.
  [OPEN] obligation 4 / plan T3.
