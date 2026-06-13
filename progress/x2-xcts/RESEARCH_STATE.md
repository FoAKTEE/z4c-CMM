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
