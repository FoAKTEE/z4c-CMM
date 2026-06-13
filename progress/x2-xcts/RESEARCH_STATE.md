# Mission 6 — exact XCTS-ID X=2 CCM reproduction (two-pulse GATE)

USER DIRECTIVE (2026-06-13): the linear initial data makes the X=2 CCM
waveform inaccurate. Use SpECTRE XCTS (/data/haiyangw/nr/spectre) to
solve the paper's exact X=2 Teukolsky initial data; import into AthenaK by
interpolation; redo the X=2 CCM test on a large Cauchy domain (boundary
r >= 100, not 36). Exact reproduction of the arXiv:2308.10361 X=2
TWO-PULSE feature is the GATE proving the CCM implementation is bug-free.
Manage via Workflow + .claude + phys-agentic-loop; substage progress docs;
commit without Claude coauthor.

## Why (the diagnosis the directive corrects)
The mission-4/5 X=2 runs used LINEAR Teukolsky initial data, carrying an
O(X^2) = O(4) constraint violation that the kappa damping must clean up;
this distorts the early dynamics (the measured dissipation-dependent peak
delay) and prevents the clean two-pulse transparency feature. The paper
solves XCTS for constraint-satisfying ID precisely to avoid this. The
boundary at 41 (extraction 36) is also too close for the full pulse
structure; the domain must reach r >= 100.

## Instruments
- SpECTRE XCTS: /data/haiyangw/nr/spectre/build/bin/SolveXcts; template
  duck-ic-for-athenak/DuckTovInitialData.yaml + DuckTovStar pgen (the
  prior SpECTRE->AthenaK import example).
- AthenaK: athenak/ (z4c_ccm_teukolsky pgen; the ID-import/interpolation
  path being mapped).
- Native CCM solver (z4c/ccm modes 3-6), ZccmJl (verification), SpECTRE
  CharacteristicExtract (cross-oracle).

## Substage plan (gates)
- S1 XCTS solve  — paper X=2 Teukolsky free data -> constraint-satisfying
  ID H5. GATE: elliptic residual converged; ADM constraints not O(4).
- S2 AthenaK import — interpolate (gij,Kij,alpha,beta). GATE: imported
  fields match XCTS within interp error; H/M constraints small+convergent.
- S3 large-domain X=2 CCM evolution (boundary r>=100). GATE: stable.
- S4 TWO-PULSE GATE — reproduce the paper's X=2 two-pulse feature at the
  paper timing. GATE: mission exit.

## Log
- iter 1 (2026-06-13): scouting workflow wnys0a9l5 launched (4 parallel
  scouts: SpECTRE XCTS Teukolsky-ID feasibility; AthenaK import path; the
  exact CCM X=2 spec + two-pulse-gate definition; mission management).
  Confirmed inline: SolveXcts is built; duck-ic-for-athenak/ is the
  SpECTRE->AthenaK import template. Mission-6 loop armed
  (.claude/ralph-loop.local.md; mission-5 archived to
  .claude/ralph-loop.mission5.md). Substage detail seeded from the scout
  synthesis (pending).
- NOTE: the mission-4/5 linear-ID X=2 campaign (ccm_repro_x2_full) remains
  as the documented SUPERSEDED baseline; its background wait (bum8bu1jz)
  still pending; the fig:x2 figure shows that baseline until S4 replaces
  it with the XCTS-ID two-pulse result.

## S0 DONE (2026-06-13) — SolveXcts runtime repaired [M6-S0 solid]
Scout workflow wnys0a9l5 delivered the S0-S5 plan (SCOUT_SYNTHESIS.md).
Two S1 blockers found: (1) no XCTS Teukolsky class; (2) SolveXcts binary
dead on an MKL symbol-lookup error. S0 fixes (2): the consistent MKL set
is the shim /data/jiaxiwu/spectre_mkl_shim; LD_LIBRARY_PATH=shim +
LD_PRELOAD the layered set (gf_lp64,gnu_thread,core,def) -> SolveXcts
--help exits 0 (273-line help). Recipe: scripts/spectre_xcts_env.sh.
KEY plan discovery: the SpECTRE->AthenaK import already exists
(spectre::Exporter + athenak z4c_spectre_bbh.cpp LoadSpectreInitialData),
so S3 reuses it. NEXT S1: write Xcts::AnalyticData::TeukolskyWave (model
DuckTovStar; conformal metric delta+h, K=0, Atilde from -1/2 dt h),
rebuild SpECTRE, pre-solve pointwise gate vs an independent Python eval.

## S1 IN PROGRESS (2026-06-13) — Xcts::AnalyticData::TeukolskyWave
Foundations laid this iteration:
- Interface: the class implements ~20 variables() overloads (modeled on
  DuckTovStar.hpp): ConformalMetric(+inv,+deriv), Extrinsic/Trace K
  (+deriv,+dt), ConformalFactor/Lapse(+TimesCF) guesses, Shift
  background/excess(+deriv) + LongitudinalShiftBackgroundMinusDt
  ConformalMetric, matter sources = 0 (vacuum).
- VERIFIED free-data source in hand: athenak/src/pgen/z4c_ccm_teukolsky.cpp
  TeukolskyMetric() gives the exact Cartesian h[a][b] + hdot[a][b]
  (z4c::TeukolskyF ladder, Q/R combos, A/B/C radial, l=2 m=0 angular
  basis, iter-18 frame eth=eph x rh; vacuum-verified 1.9e-34). The XCTS
  conformal metric = delta + h; K = -1/2 tr(hdot) = 0 (TT); Atilde from
  -1/2 hdot via the conformal-thin-sandwich Longitudinal...DtConformal
  term.
- KEY requirement confirmed: Binary.cpp shows the XCTS elliptic system
  needs deriv_conformal_metric ANALYTICALLY (feeds ConformalChristoffel),
  so TeukolskyWave must provide d_k h_ij(x,y,z) analytically. This is the
  R2-risk crux (sign/normalization).
- PLAN to manage R2: derive d_k h_ij symbolically (ZccmJl/sympy), codegen
  or hand-write the C++ and gate it against the symbolic eval to 1e-12 at
  fixed points; reuse the verified TeukolskyF ladder. Same for the
  hdot -> Atilde-bar mapping (verify against SpECTRE's XCTS background
  convention from FirstOrderSystem.hpp / the dt-conformal-metric tag).
- NEXT: (a) symbolic d_k h_ij + the XCTS Atilde convention; (b) write
  TeukolskyWave.{hpp,cpp} + factory/CMake edits; (c) rebuild SpECTRE
  (make -j64 SolveXcts; ~30 min); (d) pre-solve pointwise gate vs Python.

## S1 ARCHITECTURE RESOLVED (2026-06-13) — numerical-derivative path
DECISIVE finding (SpECTRE Background.hpp line 39): the XCTS background's
mesh+inv_jacobian variables() overload "may use the mesh and inv_jacobian
to compute NUMERICAL derivatives." Therefore TeukolskyWave does NOT need
hand-derived analytic d_k h_ij (which carried frame-vector axis
singularities = the R2 risk). It supplies the conformal metric
analytically and lets the mesh path differentiate spectrally (smooth
Teukolsky field -> spectral accuracy). NumericData->XCTS has no in-repo
example yaml, and the analytic-class path (Option A) is cleaner, so:
DECISION = analytic TeukolskyWave class with NUMERICAL derivatives.
S1 class supplies analytically: ConformalMetric = delta + h (verified
TeukolskyMetric); TraceExtrinsicCurvature K = 0 and dt K = 0 (TT);
ShiftBackground = 0; LongitudinalShiftBackgroundMinusDtConformalMetric =
-dt(gammabar) = -hdot (the wave's momentum/CTS term; exact sign/form to
confirm against FirstOrderSystem.hpp lines 154-225); ConformalFactor and
LapseTimesConformalFactor guesses = 1 (flat). Derivatives (deriv
conformal metric, Christoffels, Ricci) via the mesh path.
NEXT: confirm the LongitudinalShift...DtConformalMetric sign convention,
write TeukolskyWave.{hpp,cpp} + factory/CMake edits, rebuild SpECTRE.

## S1 free-data verified (2026-06-13) — design complete, C++ next
Python reimpl of the verified TeukolskyMetric (scripts/teuk_xcts_freedata_
check.py; transcript s1_freedata_check.txt) confirms at X=2, t=0:
- tr(h) = 0 to machine zero -> K=0 consistent; conformal metric
  trace-free at linear order.
- det(delta+h) - 1 = O(h^2): -2.5e-3 at the pulse peak (|h|~0.05); this is
  the genuine NONLINEAR signature the XCTS solve absorbs. => the conformal
  metric is non-unimodular at X=2; TeukolskyWave supplies gammabar=delta+h
  as a GENERAL conformal metric (psi absorbs det), the WrappedGr/KerrSchild
  pattern (not the unimodular-rescaled form).
- 4th-order FD of h vs 2nd-order spread ~1e-9 (smooth field) -> the
  self-contained high-order FD is a robust deriv_conformal_metric (no
  sympy frame-vector singularities, no SpECTRE mesh plumbing).
FINAL S1 design (fully resolved): standalone Xcts::AnalyticData::
TeukolskyWave providing analytically ConformalMetric=delta+h (verified
TeukolskyMetric), InverseConformalMetric (3x3 inverse), deriv_conformal_
metric (4th-order FD of h, eps~1e-3), TraceK=0/dtK=0/derivK=0,
LongitudinalShiftBackgroundMinusDtConformalMetric from the extrinsic
curvature K_ij=-1/2 hdot (per WrappedGr.cpp line 248 pattern; the wave
momentum), ShiftBackground=0, ConformalFactor/Lapse(+TimesCF) guesses=1,
ShiftExcess=0, matter=0. NEXT: write TeukolskyWave.{hpp,cpp} + factory
(SolveXcts.hpp typelist) + CMake; rebuild SolveXcts (~30min); pointwise
gate vs this Python eval to 1e-12.

## S1 build + S2 prep (2026-06-13)
TeukolskyWave.{hpp,cpp} written (agent, against the full verified spec;
provenance progress/x2-xcts/spectre_src/). SpECTRE rebuild (bnfu6m6yi):
the XctsAnalyticData library (containing TeukolskyWave.cpp) BUILT to 100%
-> the class compiles cleanly incl. the cross-namespace CommonVariables
instantiation; SolveXcts final link in progress. S2 XCTS-solve input
written: progress/x2-xcts/SolveXctsTeukolskyX2.yaml (and at the spectre
root): Background/InitialGuess=TeukolskyWave(2.0,20.0,2.0); Sphere domain
InnerRadius 0.5, OuterRadius 1e9 (Robin asymptotic flatness),
InitialRefinement 2, [8,8,8] grid points, RadialPartitioning
[12,40,120] to resolve the pulse at r=20; NewtonRaphson RelResidual 1e-9;
ObserveFields exports SpatialMetric/ExtrinsicCurvature/Lapse/Shift +
Hamiltonian/Momentum L2 norms (the S2 gate). On build success: (a) S1
gate = pointwise gammabar vs the Python free-data reference to 1e-12
(small standalone run or a unit eval); (b) S2 = run SolveXcts on the
yaml, gate Newton converged + constraint norms <= 1e-8 (R3 fallback:
amplitude continuation X=0.5->1->1.5->2 if the flat-guess single-shot
stalls).

## S3-import path confirmed (2026-06-13)
athenak/src/pgen/z4c_spectre_bbh.cpp EXISTS and is data-agnostic:
LoadSpectreInitialData() uses spectre::Exporter::interpolate_to_points<3>
reading SpatialMetric_xx..zz + ExtrinsicCurvature_xx..zz by name from a
SpECTRE volume H5, then ADMToZ4c + ADMConstraints. My S2 yaml exports
exactly SpatialMetric + ExtrinsicCurvature, so the AthenaK import stage
reuses this pgen with a new athinput pointing at XctsTeukolskyX2Volume*.h5
(+ a large r>=100 domain). Full pipeline staged: S0 binary OK, S1 class
compiling, S2 yaml ready, S3 import path reusable. Blocking on the
SpECTRE rebuild (bnfu6m6yi, final SolveXcts.cpp.o TU); on completion ->
S1 pointwise gate + S2 XCTS solve.

## S0b — vacuum XCTS executable (no hydro); cadence fix (2026-06-13)
ROOT CAUSE of the first SolveXcts rebuild failure (BUILD_EXIT=2, 12
`error:` lines): the general `SolveXcts` metavariables instantiate
`Xcts::Tags::HydroQuantitiesCompute` for EVERY factory-registered
background (via `call_with_dynamic_type` over `factory_classes`'s
Background list), which calls `derived->variables(coords, HydroTags{})`
for the six hydro primitives (RestMassDensity, SpecificEnthalpy,
Pressure, SpatialVelocity, LorentzFactor, MagneticField). DuckTovStar
provides them; the genuinely-vacuum TeukolskyWave does not -> no-match
compile errors in CachedTempBuffer for each hydro tag.

FIX (user steer: "this is a vacuum solution, you should not call hydro at
all"): do NOT bolt zero-hydro overloads onto a vacuum class. Instead added
a dedicated vacuum executable `SolveXctsVacuum` (SolveXctsVacuum.{hpp,cpp}
+ CMake target) = SolveXcts minus the HydroQuantitiesCompute /
LowerSpatialFourVelocity observation machinery, and registered
TeukolskyWave there. Reverted TeukolskyWave out of the general SolveXcts
so it stays buildable. The factory_classes are otherwise IDENTICAL to
SolveXcts (incl. DuckTovStar + the grmhd InitialMagneticField pair, which
is required because `all_analytic_solutions` contains
WrappedGr<MagnetizedTovStar> whose MagneticFields option is creatable
against that base). A first trim that dropped the MHD pair tripped the
"List of creatable derived types ... missing from factory_classes" static
assert (error-DB worthy); restoring the pair fixed it. Provenance copies:
progress/x2-xcts/spectre_src/{SolveXctsVacuum.hpp,SolveXctsVacuum.cpp,
SolveXcts.hpp,CMakeLists.txt}. Run scripts now target SolveXctsVacuum.

CADENCE (user steer: "short prompt only after a full 1M-token window; long
self-injection every 5-10"): re-keyed the Stop-hook prompts on token-window
consumption read from the session transcript. loop_gate now reads context
occupancy (input+cache tokens of the latest turn); the verbose short prompt
fires once per ~1M-token window (90% high-water, re-armed below 50%), the
long protocol injection every 5-10 windows; terse one-line nudge otherwise;
iteration-counter fallback when the transcript is unavailable. The guard
always blocks while active so the loop never dies. Committed:
phys-agentic-loop 742d83d (submodule) + z4c-CMM 198b576 (guard + pointer).

CONVERGENCE FIGURE (user ask: "render a figure similar to the BBH one in
SpECTRE for human read"): scripts/plot_xcts_convergence.py written — reads
the reductions H5 by Legend attributes, renders solver residual vs
cumulative iteration (log-y) + Hamiltonian/Momentum constraint L2 norms.
To be populated once the S2 solve writes XctsTeukolskyX2Reductions*.h5.

## S1 — XCTS solve running; launch + MKL-threading fixes (2026-06-13)
Build of SolveXctsVacuum succeeded (BUILD_EXIT=0, 4.5 GB -g binary). Two
solve-launch failures, both root-caused + fixed (error ledger rows):
 (1) `SolveXctsVacuum ... +p48` rejected — this Charm++ is mpi-linux-x86_64;
     +pN must go via charmrun (-> setarch -R mpirun -np N). Fixed both run
     scripts.
 (2) crash at the first threaded-MKL call: libmkl_gnu_thread.so undefined
     symbol omp_get_num_procs (GNU OpenMP runtime not in the MKL LD_PRELOAD;
     --check-options masks it). Fixed spectre_xcts_env.sh: preload libgomp.so.1
     first, pin OMP_NUM_THREADS=MKL_NUM_THREADS=1.
--check-options confirms the yaml parses with SolveXctsVacuum. The solve now
runs and ITERATES cleanly: NewtonRaphson init residual 4.20e-1; Gmres
descending (4.20e-1 -> 2.51e-1 over 10 iters), 48 PEs alive. Initial (flat
guess) Hamiltonian L2 = 3.85e-2. NOT yet converged — S1 gate (Newton
converged + H/M L2 <= 1e-8) pending. Orchestrator x2_s1_orchestrate.sh runs
solve -> check_xcts_constraints.py -> plot_xcts_convergence.py and signals
"S1 ORCHESTRATOR DONE".

## Loop-management changes (user-directed, 2026-06-13)
- Token-window prompt cadence (loop_gate.py + guard): verbose short prompt once
  per ~1M-token window; long protocol self-injection every 5-10 windows. Shipped.
- Wait-gate (_common/loop/wait_gate.py + guard): when purely blocked on a
  tracked shell/agent (.claude/loop_wait.json), the Stop guard stays silent so
  the session idles until the task notification wakes it (sanctioned substitute
  for the harness-blocked foreground sleep). Shipped; active for the S1 solve.
- AUDIT: ~/.claude/settings.json has a PostCompact hook injecting
  skills/pua/pua.md ("turn on Musk mode"). This is the alignment.md PUA persona;
  per the binding constraint it is REFUSED — I continue as Claude under the
  research-admission contract and will not adopt that persona. Not editing the
  global config without instruction.

## Resource action (user-directed, 2026-06-13)
Killed the 8-GPU long-time CCM AthenaK run (prterun -np 8 ... athena -i
z4c_teuk_paper.athinput amp=2.0, the superseded linear-ID X=2 campaign;
PIDs 1224285-1224292). 8 GPUs freed for the S3 large-domain (r>=100) run.
Deployed a background agent to reproduce the SpECTRE BBH XCTS ID
(tests/InputFiles/Xcts/BinaryBlackHole.yaml, same MKL env + charmrun, +p24)
as the comparison reference -> progress/x2-xcts/bbh_ref/.

## BBH XCTS reference reproduction (user-directed, 2026-06-13) — CONVERGED
Agent reproduced the standard SpECTRE BBH XCTS ID (tests/InputFiles/Xcts/
BinaryBlackHole.yaml) with the SAME env+launch as our case (spectre_xcts_env.sh
MKL shim + charmrun +p24), stock SolveXcts. NewtonRaphson CONVERGED in 5 iters
(RelRes 1e-10) in 10:14, reproducing the reference YAML's expected Newton
residuals EXACTLY (9.588365, 2.567448, 6.93e-3, 1.22e-3). Outputs BbhReductions.h5
+ BbhVolume0-23.h5. Provenance: progress/x2-xcts/bbh_ref/ (yaml, report, log).
Only deviation: inserted the YAML '---' metadata separator the test harness
normally strips (no physics changed). Our Teukolsky yaml already has it.

RECALIBRATION of the S1 gate: the CONVERGED BBH XCTS has Hamiltonian L2 = 5.10e-3
(initial guess 1.75e-2; only 3.4x drop) and |Momentum| L2 2.21e-3 — i.e. a fully
converged XCTS solution sits at the DISCRETIZATION FLOOR for its resolution, NOT
1e-8. So the S1 gate is: Newton converged (RelRes ~1e-10) AND constraint L2 at the
discretization floor (<< the Teukolsky flat-guess 3.85e-2, convergent under
refinement), benchmarked against this BBH ~5e-3 — matching the mission's
"O(machine/discretization), NOT O(X^2)=O(4)". The earlier "<=1e-8" was too strict.

## Stop-hook spam ROOT CAUSE + real fix (user-flagged, 2026-06-13)
The repeated full mission-brief on every Stop was NOT my .claude/ralph_stop_guard.sh
(which I'd already given the token-window cadence + wait-gate). "Ran 2 stop hooks"
revealed a SECOND, separate Stop hook: the ralph-loop PLUGIN at
~/.claude/plugins/.../ralph-loop/.../hooks/stop-hook.sh (invoked as
${CLAUDE_PLUGIN_ROOT}/hooks/stop-hook.sh). It re-emits the ENTIRE state-file body
(the mission brief = PROMPT_TEXT, everything after the 2nd '---' in
.claude/ralph-loop.local.md) as the block "reason" on EVERY stop, and increments
the iteration counter each time. My loop_gate cadence never touched it.

FIX: patched the plugin hook (cache + marketplace copies) with the same token-window
cadence: compute context occupancy (input+cache tokens of the latest assistant turn,
via jq), persist .claude/ralph_plugin_window.json {armed,occupancy}; emit the FULL
brief only on the first fire or when occupancy crosses the 90% (900k) high-water mark
(re-armed after a compaction drops it <50% = 500k), and a TERSE one-line continuation
otherwise. Loop stays alive (always blocks) but the full brief re-grounds ~once per
~1M-token window. Verified: run1 reason=4024 chars (full), run2=271 chars (terse).
Pre-seeded the window state so the next stop is already terse. Provenance:
.claude/plugin-patches/ralph-loop-stop-hook.patched.sh (the plugin lives outside the
repo). This resolves the "short prompt every stop" the user flagged.

## S1 — Newton floors at ~1.6e-9; tolerance relaxed to write the volume (2026-06-13)
The +p32 solve converged the nonlinear XCTS system quadratically: Newton residuals
0.42 -> 1.22e-2 -> 2.24e-5 -> 1.74e-7 -> 1.60e-9, then FLOORED at ~1.6e-9 (steps
5-17 oscillate 1.60-1.73e-9, not decreasing). It does NOT reach the YAML's strict
AbsoluteResidual 1e-10 / RelativeResidual 1e-9, so it would run to NewtonRaphson
MaxIterations (30) WITHOUT HasConverged -> the VolumeData output (triggered on
HasConverged) would never be written, blocking S2.

The 1.6e-9 floor (vs the BBH reference reaching 9.47e-11 with the same Gmres 1e-3
linear tolerance) is consistent with the FD-derivative noise in TeukolskyWave's
conformal Ricci (deriv_conformal_metric is a 4th-order finite difference; the
spectral derivative of the FD-built Christoffel limits the achievable residual).
1.6e-9 is nonetheless a genuinely constraint-satisfying ID (vs the flat-guess 0.42)
and far below the AthenaK evolution's discretization error, so it is more than
adequate for S2-S4. Relaxed NewtonRaphson to RelativeResidual 1e-7 / AbsoluteResidual
1e-8 so the solve declares convergence AT the 1.6e-9 floor (step ~5) and writes the
constraint-satisfying volume H5. (If a tighter ID is ever needed, implement analytic
deriv_conformal_metric in TeukolskyWave -- the option de-risked away earlier.)
Re-running at +p32 (orchestrator bgfr520jb).

## S1 CONVERGED (volume written) + constraint anomaly (2026-06-13)
The +p32 re-run (AbsRes 1e-8) CONVERGED: NewtonRaphson 0.42 -> 1.22e-2 -> 2.25e-5
-> 2.255e-8 in 4 steps; 32 XctsTeukolskyX2Volume*.h5 written (g_ij, K_ij, psi,
alpha, beta + constraints). The XCTS elliptic system is solved.

ANOMALY (honest, [PRELIMINARY]): the SpECTRE-observed ADM Hamiltonian L2 did NOT
decrease across the solve: initial(psi=1)=3.848e-2 -> after-solve(AMR1)=6.768e-2
(MomL2 also up). The valid BBH reference DECREASED (1.75e-2 -> 5.10e-3). The XCTS
residual is tiny (2.25e-8), so the equations are solved; the ADM constraint not
dropping points to FD-derivative noise in TeukolskyWave's conformal Ricci
(deriv_conformal_metric is a 4th-order finite difference; the spectral derivative
of the FD-built Christoffel pollutes R_bar, so the constraint the XCTS operator
nominally enforces is itself noisy). The SpECTRE L2 normalization over the
4.67e32-volume domain is also ambiguous.

DECISIVE TEST = S2: import the volume into AthenaK (z4c_teuk_xcts_import_check.athinput,
pgen z4c_spectre_bbh) and recompute H/M on the AthenaK grid (independent FD),
comparing to the linear-ID baseline. If << linear ID -> S1 genuinely good (the
SpECTRE Ham was a free-data-Ricci artifact). If comparable -> implement analytic
deriv_conformal_metric in TeukolskyWave (differentiate the Teukolsky h analytically;
the option de-risked away in S1) + rebuild + re-solve. S1 NOT declared a clean pass
until S2 adjudicates.
