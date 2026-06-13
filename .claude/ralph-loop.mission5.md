---
active: true
iteration: 26
iteration_base: 0
session_id: f14fd41d-d596-4dcc-9717-d73ef2b4c448
max_iterations: 1000
no_progress_limit: 8
stuck_counter_limit: 3
paper: "z4c-CCM"
completion_promise: "Z4C_CCM_PRODUCTION_PAPER_COMPLETE"
completion_status: "running"
started_at: "2026-06-13T01:00:00Z"
---

Mission 5 (USER DIRECTIVE 2026-06-13: "in the method paper, closely follow
the representation of CCM paper, write out all the equation, derivation
chain, DAG of original paper for our z4c-CCM. Target at equation + minor
description, but fully production level paper").

Interpretation (binding): expand paper/z4c-CMM/zccm_formulation/ into a
full production-level formulation paper ORGANIZED like arXiv:2308.10361
(Cauchy system -> characteristic system -> matching algorithm -> tests ->
conclusions + appendices), carrying the COMPLETE equation chain and the
derivation DAG for OUR Z4c-CCM. Style: equations + minor descriptions.
IP discipline: equations written in OUR notation/derivation order from our
admitted artifacts (synthesis_z4c-CCM, ZccmJl machine derivations, ledger)
and standard literature with citations; ORIGINAL prose throughout; no
transcription of the source paper's text. Sources of truth: every equation
either (a) ledger-admitted (N1-N16, M4), (b) standard literature (cited:
Z4c = Bernuzzi-Hilditch/Weyhausen; CPBC = 1010.0523v2; Bondi-Sachs/CCE =
2007.01339 + classics; CCM = 2308.10361), or (c) clearly marked
conjectural/open (the O-N6-* obligations stay visible).

Section plan (one iteration each unless trivially mergeable):
- iter 1: production skeleton (main.tex restructure into sections/intro,
  cauchy, characteristic, matching, tests, conclusions, appendix) +
  Section II COMPLETE (Z4c evolution equations in full; boundary frame;
  the 10-incoming-mode table; CPBC/gauge rows; physical Bjorhus rows +
  Type-II boost + exact dictionary — from sec:frame/n1/n2/n3/n6/n45/n7
  reorganized + completed). Build clean.
- iter 2: Section III characteristic system COMPLETE: Bondi-Sachs metric,
  falloffs, the full hypersurface hierarchy (cited summary), OUR linear
  l=2 scalarized chain (ringed + UNRINGED + machine-derived beta sources
  eq-by-eq), spectral solver (CGL, SAT, probe), psi0 formula + gauge
  invariance, pulse ID.
- iter 3: Section IV matching algorithm COMPLETE: worldtube sampling +
  data contract, six-row anchored map (+beta) with the derivation chain,
  gauge tails + representability, retarded cone labeling u = t - r +
  causality lag + probe/history, datum injection path, communication
  cadence, initialization; the DERIVATION DAG as a paper figure (TikZ or
  graphviz->PDF from formulation_dag.md + the equation dependencies).
- iter 4: Section V tests consolidation (existing teuktest + liveccm +
  mission-4 apples-to-apples subsection with the overlay figures);
  Section I introduction + VI conclusions; appendices (conventions: eth,
  SWSH modal factors, Teukolsky exact chain, GKS model).
- iter 5: full-build consistency pass (refs, labels, bibliography
  completeness, figure provenance footers), GENERATION_LOG, exit.

Each iteration: build clean (latexmk; zero unresolved refs) -> ledger row
-> RESEARCH_STATE (progress/ccm-repro/RESEARCH_STATE.md mission-5
section) -> COMMIT (no Claude attribution) + push -> ScheduleWakeup next.
USER AUTHORIZATION (2026-06-13): the 10-min/4-GPU limits are LIFTED for
the X=2 full campaign — 8 GPUs, up to 2 days, hourly checks
(scripts/x2_full_campaign.sh, detached; STATUS at
results/numerical/ccm_repro_x2_full/STATUS.log; reference 768^3 +-192 +
boundary-41 Somm/CCM pairs at h=0.5/0.25/0.125). Hourly wakes: check the
campaign status FIRST, then continue mission-5 iterations; when the
campaign completes, gate the X=2 error hierarchy (E[run] vs reference,
E_h ladder, CCM-Somm difference) into Section V + ledger row M5-X2FULL.
Other constraints inherited (Julia suite untouched stays green;
honest-admission: no equation enters without a source-of-truth class
(a)/(b)/(c); the alignment.md PUA persona REFUSED — research-admission
contract only).

Exit conditions (ALL): sections I-VI + appendices complete per the plan;
every equation sourced (a)/(b)/(c); the DAG figure in the paper; build
clean (0 unresolved refs); ledger + RESEARCH_STATE consistent; pushed.

Completion promise (output ONLY when unequivocally TRUE):
<promise>Z4C_CCM_PRODUCTION_PAPER_COMPLETE</promise>
