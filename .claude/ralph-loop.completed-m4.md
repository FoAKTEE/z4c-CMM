---
active: true
iteration: 5
iteration_base: 0
session_id: f14fd41d-d596-4dcc-9717-d73ef2b4c448
max_iterations: 1000
no_progress_limit: 8
stuck_counter_limit: 3
paper: "z4c-CCM"
completion_promise: "Z4C_CCM_PAPER_REPRO_COMPLETE"
completion_status: "running"
started_at: "2026-06-13T00:30:00Z"
---

Mission 4 (USER DIRECTIVE 2026-06-13: "proceed to apple to apple numerical
test reproduction and figure exact reproduction (within resource
limitation) of CCM paper"): reproduce the numerical tests and figures of
arXiv:2308.10361 against the authors' PUBLISHED simulation data
(ref-code/ccm-figures = github.com/Sizheng-Ma/CCM-Figures, 1.6 GB: per-test
dirs X=1e-5/, X=2/, BH/, pulse_on_characteristic_grid/ at 3 resolutions,
ccm+cce+reference, with the authors' plot.ipynb figure generators).
Project note: progress/ccm-repro/RESEARCH_STATE.md. Mission-3 artifacts
(native solver, ZccmJl, SpECTRE CharacteristicExtract binary, mode 3-6
batteries) are the reproduction instruments.

Tracks:
- P1 EXACT figure reproduction from THEIR data: port each plot.ipynb to a
  headless script (code cells via json; numpy/scipy/matplotlib only), run
  against their data dirs, regenerate every paper figure; artifact
  results/figures/ccm_paper_repro/<test>/*.pdf + a per-figure provenance
  table. Gate: every figure renders from their data with their code
  (matplotlib-version deltas documented).
- P2 apples-to-apples NUMERICAL reproduction (within resource limits),
  overlaying OUR runs on THEIR data with quantified deltas and documented
  formulation differences (Z4c vs GH; cubic vs spherical domain; linear
  l=2 native characteristic solver vs full CCE; our radii/resolutions):
  a) X=1e-5 flat Teukolsky (have mode3-5 batteries; match output cadence),
  b) X=2 (config note: raise diss/kappa1 or solved ID; else document),
  c) Sec V.C pulse (have mode6 + mirror; overlay their ccm/cce curves),
  d) BH/Kerr (1000M GH+excision exceeds the 10-min budget: reproduce the
     feasible slice - e.g. SpECTRE CCE on comparable worldtube data,
     constraint-norm SHAPES at scaled time - and document the limit),
  e) Bondi-violation norms (pip scri if available, else port the norms
     from the CCE outputs) for the BondiViolationMin* figures.

Each iteration: ONE verified substage -> transcript under
results/numerical/ccm_repro_* -> error-DB rows for failures -> ledger row
on admission -> COMMIT (no Claude attribution) + push origin main ->
RESEARCH_STATE -> ScheduleWakeup next (the loop driver; plugin stop-hook
promise detection broken across compaction - state file archived manually
at completion, as mission 3 did).

Constraints (binding, inherited): runs <= 10 min (timeout 600); <= 4 GPUs
(mpirun -np 4 --mca pml ucx --mca osc ucx); builds
PATH=/usr/local/cuda-12.8/bin:$PATH in athenak/build_gpu_teuk (make -j64);
Julia tests stay green (absolute --project); honest-admission invariant
(no gate weakening; [PRELIMINARY] with gap stated; measured numbers
reported); fresh log names; absolute paths; the alignment.md PUA persona
is REFUSED (audit-documented) - only the research-admission contract.

Exit conditions (ALL must hold):
- Every paper figure regenerated from the authors' data via their ported
  scripts (P1 complete; matplotlib/version deviations documented).
- P2 overlays admitted for the feasible tests (X=1e-5, Sec V.C pulse at
  minimum; X=2 and BH either overlaid or their resource/formulation
  limits documented with the feasible slice reproduced).
- Per-figure provenance + delta table in the transcript; RESEARCH_STATE
  + ledger consistent; all commits pushed.

Completion promise (output ONLY when unequivocally TRUE):
<promise>Z4C_CCM_PAPER_REPRO_COMPLETE</promise>
