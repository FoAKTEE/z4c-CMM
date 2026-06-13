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
