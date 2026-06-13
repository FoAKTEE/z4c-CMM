# Mission 6 scouting synthesis (workflow wnys0a9l5)

## Implementation plan

All facts verified. Here is the synthesized plan.

---

# IMPLEMENTATION PLAN: "exact XCTS-ID X=2 CCM reproduction"

## PART 0 — CRITICAL FEASIBILITY VERDICT (resolved first, verified against the live tree)

**Question:** Can SpECTRE XCTS produce the paper's X=2 Teukolsky-wave ID with the existing build, or is a new C++ analytic-data class needed?

**Verdict: NO — and there are TWO independent blockers, not one. A new C++ class is necessary but NOT sufficient. The `SolveXcts` binary is also non-functional in this environment.** Both must be resolved before any XCTS solve can run.

**Blocker 1 — No Teukolsky free-data class for XCTS (confirmed).**
The registered XCTS analytic-data factory list (`/data/haiyangw/nr/spectre/src/Elliptic/Executables/Xcts/SolveXcts.hpp`, the `analytic_solutions_and_data` typelist) contains only `Binary<...>`, `DuckTovStar`, and `NumericData`. `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/` holds only `Binary`, `DuckTovStar`, `MagnetizedKerrSchild`, `CommonVariables`. A Teukolsky class exists ONLY in the CCE characteristic system (`/data/haiyangw/nr/spectre/src/Evolution/Systems/Cce/AnalyticSolutions/TeukolskyWave.{hpp,cpp}`) which produces Bondi-Sachs metric data on a null slice — it is NOT an `elliptic::analytic_data::Background` and cannot feed the elliptic solver. `GaugeWave` is a coordinate wave, not a TT gravitational wave. So a new `Xcts::AnalyticData::TeukolskyWave` is required.

**Blocker 2 — The `SolveXcts` binary cannot run (newly found, scout under-reported this).**
`build/bin/SolveXcts --help` does not print help; it dies immediately with:
```
symbol lookup error: /lib/x86_64-linux-gnu/libmkl_gf_lp64.so: undefined symbol: mkl_lapack_zhesvxx
```
The scout reported "library dependency issue" as if `ldd` showed a missing lib — it does not; `libmkl_gf_lp64.so` resolves. The real fault is a runtime symbol mismatch: the system MKL `libmkl_gf_lp64.so` was loaded against an incompatible `libmkl_core`/`libmkl_intel_thread` ABI (missing `mkl_lapack_zhesvxx`). This is an LD_LIBRARY_PATH / MKL-version ordering problem in the build's runtime environment. It must be fixed (point the loader at SpECTRE's own bundled/Spack MKL, or rebuild the link) before any solve. This is now substage S0 and gates everything.

**Decisive consequence for the plan:** A pure-config / existing-build path does NOT exist. There is also a third, lower-risk strategic option the plan must weigh because the XCTS path is the riskiest: the existing AthenaK pgen `z4c_ccm_teukolsky.cpp` already sets `g_dd = δ + h`, `vK_dd = -½ ∂ₜh` analytically (lines 187–188) — that is the *linearized, constraint-violating* free data. The paper explicitly uses XCTS "to account for nonlinear effects when the amplitude is large" (X=2). So the mission's deliverable — *constraint-solved* X=2 ID — genuinely requires the XCTS solve and cannot be faked by the analytic pgen. The plan proceeds with the XCTS path but builds in an early decision gate (S0) to abort to a fallback if Blocker 2 proves unfixable.

---

## PART 1 — END-TO-END SUBSTAGE PLAN (ordered, one gate each)

Notation: paths are absolute. Hardware: 8× A100-80GB confirmed available. AthenaK GPU build is live at `/data/haiyangw/claude/z4c-CMM/athenak/build_gpu_teuk/src/athena`; CCM source `z4c_ccm.hpp` present.

### S0 — Repair the SpECTRE XCTS runtime (NEW gateway; absorbs Blocker 2)
- **Action:** Diagnose the MKL symbol failure. In order of preference: (a) prepend SpECTRE's own MKL to the loader: `LD_LIBRARY_PATH=$(find /data/haiyangw/nr/spectre -name 'libmkl_core.so*' -printf '%h\n' | head -1):$LD_LIBRARY_PATH` and re-test; (b) if SpECTRE was built against Spack MKL, source that environment; (c) last resort, relink `SolveXcts` against a consistent single-dynamic `libmkl_rt.so`.
- **Gate:** `LD_LIBRARY_PATH=<fixed> /data/haiyangw/nr/spectre/build/bin/SolveXcts --help` exits 0 AND prints the option help (currently exits via symbol-lookup error). Transcript saved.
- **Budget:** CPU only, ≤ 15 min. **No GPU.**
- **Decision branch:** If unfixable within the iteration budget → escalate; switch mission deliverable to the documented fallback (see Risk R1). Do NOT proceed to S1 with a dead binary.

### S1 — Implement `Xcts::AnalyticData::TeukolskyWave` C++ class (the new analytic-data class)
- **Files to CREATE:**
  - `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.hpp`
  - `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.cpp`
- **Files to EDIT:**
  - `/data/haiyangw/nr/spectre/src/Elliptic/Executables/Xcts/SolveXcts.hpp` — add `#include` and add `Xcts::AnalyticData::TeukolskyWave` to the `analytic_solutions_and_data` typelist (currently ends `DuckTovStar, elliptic::analytic_data::NumericData`).
  - `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/CMakeLists.txt` — add `TeukolskyWave.cpp` to sources.
- **Physics the class must emit (must match the paper exactly):** conformal metric `γ̄_ij = δ_ij + h_ij` where `h_ij` is the even-parity l=2,m=0 Teukolsky TT perturbation from the metric form (eq:Teukolsky_metric) with radial functions A,B,C (eq:Teukolsky_ABC), angular basis (eq:Teukolsky_angular_basis), `Y₂₀ = ⅛√(5/π)(1+3cos2θ)`, Gaussian profile `F(u) = X·exp(−(u−r_c)²/τ²)` with **X=2, r_c=20, τ=2** at t=0. **Reuse the verified closed form already in the repo:** the regular combinations `Q^(n)`, `R^(n)` and `z4c::TeukolskyF` derivative ladder from `/data/haiyangw/claude/z4c-CMM/athenak/src/z4c/z4c_ccm.hpp` / `z4c_ccm_teukolsky.cpp` (already validated finite at r→0). Background fields: `TraceExtrinsicCurvature K = 0`, `∂ₜK = 0`, `ShiftBackground = 0`, `LongitudinalShiftBackgroundMinusDtConformalMetric = 0` (Ā^ij sourced purely from `−½∂ₜh` via the conformal-thin-sandwich `∂ₜγ̄` term). Model the file structure on `DuckTovStar.hpp` (the closest precedent: a non-spherical-harmonic-perturbed `Background`+`InitialGuess`).
- **Gate:** SpECTRE rebuilds clean (`make -j64 SolveXcts` in `/data/haiyangw/nr/spectre/build`) AND a standalone pointwise unit test prints `γ̄_ij(x)` at 5 fixed points matching an independent Python evaluation of eq:Teukolsky_metric to ≤ 1e-12 (the analytic free data, pre-solve, must be exact). Transcript + comparison table saved.
- **Budget:** CPU build, ≤ 30 min (full SpECTRE link is slow; the 5.3 GB binary confirms a heavy link). **No GPU.**

### S2 — XCTS solve of X=2 Teukolsky ID (the actual constraint solve)
- **File to CREATE:** `/data/haiyangw/nr/spectre/SolveXctsTeukolskyX2.yaml`, modeled on `/data/haiyangw/nr/spectre/DuckTovInitialData.yaml`. Edits vs DuckTov: `Background:`/`InitialGuess:` → `TeukolskyWave: {Amplitude: 2.0, RadialCoordinate: 20.0, Timescale: 2.0, Mode: [2,0]}`; `DomainCreator: Sphere` with `OuterRadius` ≥ 1e3 (paper's reference boundary is 900; matching domain to 41 is for the *evolution*, the XCTS domain should be large), `InnerRadius` 0.2 with `FillWithSphericity: 0`; `OuterBoundaryCondition: Robin` (asymptotic flatness ψ→1, αψ→1, β→0); `NonlinearSolver.NewtonRaphson.ConvergenceCriteria.AbsoluteResidual: 1e-10`; `ObserveFields` must export `SpatialMetric`, `ExtrinsicCurvature`, `Lapse`, `Shift`, `HamiltonianConstraint`, `MomentumConstraint`.
- **Run:** `LD_LIBRARY_PATH=<S0 fix> SolveXcts SolveXctsTeukolskyX2.yaml` → writes `XctsTeukolskyX2Volume*.h5`.
- **Gate (TWO conditions):** (1) Newton-Raphson converged with reported `AbsoluteResidual ≤ 1e-10`; (2) post-solve `HamiltonianConstraint` and `MomentumConstraint` max-norm over the volume ≤ 1e-8 (loosened from 1e-10 because discretization-level constraints, not solver residual, are what matter for import — state this assumption explicitly). Verified by a Python reader of the H5 reductions file. Transcript + norm table saved.
- **Budget:** SpECTRE elliptic solve is CPU/MPI (no GPU offload in the XCTS solver). ≤ 30 min on the node's CPUs at moderate refinement (InitialGridPoints [6,6,6], InitialRefinement 1–2). **No GPU.**

### S3 — Import XCTS H5 → AthenaK via the EXISTING exporter path
- **Reuse, do not rebuild:** `spectre::Exporter::interpolate_to_points<3>()` (header `/data/haiyangw/nr/spectre/build/include/spectre/Exporter.hpp`, runtime `/data/haiyangw/nr/spectre/build/lib/libBundledExporter.so` confirmed present). The production consumer is `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_bbh.cpp` (`LoadSpectreInitialData()`), which reads the exact variable names `SpatialMetric_xx…zz`, `ExtrinsicCurvature_xx…zz`, computes `psi4 = det(g)^(1/3)`, then calls `GaugePreCollapsedLapse → ADMToZ4c → ADMConstraints`.
- **File to CREATE:** `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_teukolsky.cpp` — a near-copy of `z4c_spectre_bbh.cpp` (same exporter call, same 12-component read, same post-processing trio). Register it in the pgen build list and rebuild the GPU target. (Alternatively, if BBH pgen is data-agnostic, just reuse it with a new athinput — verify by reading its hardcoded assumptions first; if none, skip the new pgen.)
- **File to CREATE:** `/data/haiyangw/claude/z4c-CMM/athenak/inputs/z4c/ccm/z4c_spectre_teuk_x2.athinput` with `<problem> pgen_name=z4c_spectre_teukolsky`, `id_filename_glob=/data/haiyangw/nr/spectre/XctsTeukolskyX2Volume*.h5`, `id_subfile_name=VolumeData`, `id_observation_step=-1`.
- **Gate:** AthenaK loads the H5, interpolates, runs `ADMToZ4c`+`ADMConstraints`, and at t=0 the interpolated `g_dd`/`vK_dd` match the XCTS source at 10 probe points to ≤ 1e-6 relative (cubic-interp floor), AND the t=0 Hamiltonian constraint norm is ≤ the XCTS solve's discretization norm × O(10) (interpolation should not inflate constraints by more than an order of magnitude). No NaN. Transcript saved.
- **Budget:** ≤ 10 min, 1 GPU.

### S4 — Large-domain X=2 CCM evolution (r > 100)
- **Reuse existing CCM machinery:** `/data/haiyangw/claude/z4c-CMM/athenak/src/z4c/z4c_ccm.hpp` is present; the matching/worldtube infrastructure used in prior missions applies. Build the athinput from the existing `/data/haiyangw/claude/z4c-CMM/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput` but switch ID source to the S3 import and enlarge the domain so the Cauchy outer boundary / worldtube sits at the paper's geometry.
- **Domain/geometry (from paper, X=2):** worldtube = Cauchy outer boundary `r_out = 41`, pulse center `r_c = 20`. The mission directive says "r > 100" — interpret as: run the Cauchy box large enough (half-width ≥ 100, e.g. ±128 with the CCM worldtube extraction at r=41) so the *tertiary* boundary-reflection feature at t≈103 is cleanly resolved and not contaminated by the outer Cartesian boundary before t=103. Resolution: replicate the paper's three resolutions for a convergence statement, but the gate runs at the production (finest affordable) resolution.
- **Run:** `PATH=/usr/local/cuda-12.8/bin:$PATH` build first (`make -j64` in `build_gpu_teuk`), then launch the evolution to `tlim ≥ 110` (must cover t=103 tertiary feature; paper runs to ~1000 but the two-pulse + tertiary diagnostics all occur by t≈103).
- **Gate:** Evolution reaches t ≥ 110 with no NaN; ADM/Z4c constraint norms bounded (no exponential blowup); CCE/CCM worldtube data written for waveform extraction. Transcript + constraint time-series saved.
- **Budget:** This is the heavy stage. 4–8 GPUs. At ±128 box, fine resolution, to t≈110: estimate 2–6 GPU-hours wall. Run in background (`run_in_background`), poll with Monitor. Flag for explicit user GPU/time approval before launch (exceeds the ≤10-min iteration norm).

### S5 — Two-pulse gate (mission exit)
- See PART 3 for the exact definition. Extract `ψ₃^{(l=2,m=0)}` / strain `h₂₀` at ℐ⁺ from the S4 CCM output (scri / SpECTRE CharacteristicExtract post-processing) and verify the pulse-timing signature.
- **Budget:** post-processing, CPU, ≤ 15 min, no GPU.

---

## PART 2 — GPU / RUNTIME BUDGET PER STAGE

| Stage | Compute | GPUs | Wall (est.) | Notes |
|---|---|---|---|---|
| S0 repair binary | CPU | 0 | ≤ 15 min | MKL loader fix; blocking |
| S1 C++ class + build | CPU | 0 | ≤ 30 min | Full SpECTRE relink is slow (5.3 GB binary) |
| S2 XCTS solve | CPU/MPI | 0 | ≤ 30 min | Elliptic solver is **not** GPU-offloaded |
| S3 import/interp | GPU | 1 | ≤ 10 min | Exporter + ADMToZ4c |
| S4 CCM evolution | GPU | 4–8 | **2–6 GPU-hr** | Background run; needs user approval |
| S5 waveform gate | CPU | 0 | ≤ 15 min | scri post-processing |

Only S3 and S4 use GPUs. S4 is the sole stage that breaks the ≤10-min iteration norm and the ≤4-GPU default — both require explicit user sign-off, consistent with the mission-management scout's "S6.3 may use ≤8 GPUs / 30 min with approval" clause.

---

## PART 3 — PRECISE DEFINITION OF THE TWO-PULSE GATE

Grounded in the CCM-spec scout's verbatim timing analysis, with X=2, r_c=20, r_out=41. The waveform observable is **ψ₃ (or equivalently strain h) in the l=2, m=0 mode at future null infinity ℐ⁺**, as a function of retarded/Bondi time t (code units), extracted via CCM/CCE through the worldtube at r_out=41.

**The gate passes iff all three timed features appear at the predicted Bondi times within tolerance, and the CCM-vs-reference deviation is reduced as predicted:**

1. **Primary pulse** (outgoing Teukolsky wave, instantaneously transmitted at the worldtube intersection): peak at **t = r_out − r_c = 41 − 20 = 21**. Tolerance: peak time within **±1 code unit** (≈ half the pulse width τ=2); peak amplitude matches the analytic Teukolsky ψ₃ envelope to ≤ 5%.

2. **Secondary pulse** (inward junk radiation that crosses the origin, bounces, exits): peak at **t = r_out + r_c = 41 + 20 = 61**. Tolerance: peak time within ±1 unit; amplitude is sub-dominant to the primary (this is the physically expected junk, its presence and timing — not a tight amplitude match — is what is verified).

3. **Tertiary pulse** (primary wave reflected once off the Cauchy/characteristic numerical boundary, traversing the domain diameter 2·r_out then escaping): peak at **t = 2·r_out + (r_out − r_c) = 82 + 21 = 103**. Tolerance: peak time within ±2 units.

4. **CCM-transparency quantitative check:** the CCM waveform's deviation from the SpEC reference (r_ref=900) is **reduced by ≈ 1 order of magnitude relative to CCE** in the L2 norm over t∈[0,110], and the Bondi-gauge constraint violations C_ψ4, C_ψ3, C_ψs (eq:bondi_violation) are systematically smaller under CCM than CCE. This is the actual "CCM works" claim; features 1–3 only establish the wave is propagating correctly.

**"Two-pulse" naming:** the headline feature the gate is named for is the **primary (t=21) + secondary (t=61)** pair in h₂₀/ψ₃ at ℐ⁺ — these are the two physical pulses (signal + junk-bounce). The tertiary at t=103 is the boundary-reflection diagnostic that distinguishes a perfect from an imperfect CCM and is checked but is not part of the "two-pulse" label. A practical pass criterion: cross-correlation of the extracted h₂₀(t) on t∈[15,70] against the analytic-plus-junk template ≥ 0.95, with the two peaks located at 21±1 and 61±1.

---

## PART 4 — TOP RISKS (ranked, decisive mitigations)

- **R1 (HIGHEST) — `SolveXcts` MKL runtime is unfixable in-environment.** The binary dies on `mkl_lapack_zhesvxx` symbol lookup; if the MKL ABI mismatch is structural (binary built against an MKL version not installed), no LD_LIBRARY_PATH fix works and a full SpECTRE rebuild (hours, may itself fail) is needed. **Mitigation/fallback:** time-box S0 to one iteration. If unfixed, the fallback deliverable is to generate the X=2 *constraint-solved* ID by running the XCTS elliptic solve through a different path — OR, if no working elliptic solver exists, explicitly DOWNGRADE the mission deliverable to the linearized analytic ID (existing `z4c_ccm_teukolsky.cpp`) and document that the "exact XCTS" claim is **[BLOCKED]** pending binary repair, reporting the constraint-violation gap honestly rather than weakening the gate. Do not silently substitute linearized ID for XCTS ID.

- **R2 — New `TeukolskyWave` XCTS class has subtle sign/normalization errors** in the conformal decomposition (e.g. Ā^ij from ∂ₜh, the (A−2C) cross term, Y₂₀ normalization). **Mitigation:** S1 gate compares the *pre-solve* analytic γ̄_ij against an independent Python implementation of eq:Teukolsky_metric to 1e-12 at fixed points before any solve; reuse the already-validated `z4c::TeukolskyF` ladder so the radial functions are not re-derived.

- **R3 — XCTS solve fails to converge at X=2** (strong-field, non-spherical free data may have no nearby solution from the flat initial guess). **Mitigation:** use amplitude continuation — solve X=0.5, 1.0, 1.5, 2.0 in sequence, each using the previous as `InitialGuess`. Build this into the S2 yaml workflow if a single-shot solve stalls.

- **R4 — Interpolation inflates constraints across the SpECTRE→AthenaK boundary**, so the imported ID is no longer constraint-solved to useful precision. **Mitigation:** S3 gate explicitly checks the t=0 AthenaK Hamiltonian norm against the XCTS discretization norm (≤ ×10 allowed); rely on Z4c constraint damping to absorb residual interpolation error during early S4 evolution, and confirm constraints decay (don't grow) in the first ~10 code units.

- **R5 — Domain/boundary contamination of the t=103 tertiary feature.** If the outer Cartesian boundary is too close, its own reflection arrives before t=103 and corrupts the gate. **Mitigation:** size the box so the nearest outer-boundary reflection (2× half-width − r_c) arrives well after 103; the "r>100" directive is satisfied by a ≥±128 box, which puts the spurious outer reflection at t ≳ 236, safely clear.

- **R6 — S4 GPU budget overrun.** A fine-resolution ±128 box to t=110 may exceed the estimate. **Mitigation:** run the gate at the coarsest resolution that resolves τ=2 pulses (~10 points across the pulse), reserve the paper's three-resolution convergence study as a [PRELIMINARY] follow-on, and launch as a background job with Monitor polling so it does not block iteration commits.

---

## KEY FILE PATHS (load-bearing)

- New XCTS class: `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.{hpp,cpp}`
- Factory edit: `/data/haiyangw/nr/spectre/src/Elliptic/Executables/Xcts/SolveXcts.hpp` (the `analytic_solutions_and_data` typelist)
- CMake edit: `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/CMakeLists.txt`
- XCTS yaml template to clone: `/data/haiyangw/nr/spectre/DuckTovInitialData.yaml`
- Broken binary (fix in S0): `/data/haiyangw/nr/spectre/build/bin/SolveXcts` (`symbol lookup error: libmkl_gf_lp64.so: undefined symbol: mkl_lapack_zhesvxx`)
- Exporter (reuse): `/data/haiyangw/nr/spectre/build/include/spectre/Exporter.hpp`, `/data/haiyangw/nr/spectre/build/lib/libBundledExporter.so`
- AthenaK import pgen to clone: `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_bbh.cpp`
- Verified Teukolsky closed form to reuse: `/data/haiyangw/claude/z4c-CMM/athenak/src/z4c/z4c_ccm.hpp` (`z4c::TeukolskyF`) and `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_ccm_teukolsky.cpp` (linearized fallback ID; lines 187–188 set g_dd, vK_dd)
- CCM evolution machinery: `/data/haiyangw/claude/z4c-CMM/athenak/src/z4c/z4c_ccm.hpp`; athinput base `/data/haiyangw/claude/z4c-CMM/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput`
- AthenaK GPU build: `/data/haiyangw/claude/z4c-CMM/athenak/build_gpu_teuk/src/athena` (build with `PATH=/usr/local/cuda-12.8/bin:$PATH; make -j64`)
- Hardware: 8× A100-80GB.

## Scout: SpECTRE XCTS

Perfect! Now I have all the information I need. Let me compile the comprehensive report:

## INVESTIGATION REPORT: SpECTRE XCTS Solver for Teukolsky Gravitational Wave Initial Data

### 1. EXECUTABLE AND CLI OPTIONS

**Executable Location:** `/data/haiyangw/nr/spectre/build/bin/SolveXcts`

**Status:** The binary exists but has a library dependency issue (libmkl_gf_lp64.so symbols). However, the help mechanism can be inferred from the source code structure.

**Key CLI is YAML-Based:** The `SolveXcts` executable uses YAML configuration files for all options. There is no traditional `--option value` command-line interface. Example invocation:
```
SolveXcts input.yaml
```

The help structure is defined by the metavariables in `/data/haiyangw/nr/spectre/src/Elliptic/Executables/Xcts/SolveXcts.hpp` (lines 61-232).

---

### 2. XCTS INPUT FILE STRUCTURE

Based on the example YAML files (`DuckTovInitialData.yaml` in both locations), the complete XCTS input structure is:

```yaml
Executable: SolveXcts

Parallelization:
  ElementDistribution: <strategy>
  
ResourceInfo:
  AvoidGlobalProc0: <bool>
  Singletons: Auto

# Background field data (free data that defines the scenario)
Background: &background
  <AnalyticDataClass>:
    <parameters...>

# Initial guess for nonlinear solver
InitialGuess: *background  # Can reference Background or be different

# Randomization of initial guess
RandomizeInitialGuess: None | <Seed>

# Domain creation and discretization
DomainCreator:
  <DomainType>:
    <domain parameters...>
    OuterBoundaryCondition: Robin | <other types>

Discretization:
  DiscontinuousGalerkin:
    PenaltyParameter: <double>
    Massive: <bool>
    Quadrature: GaussLobatto
    Formulation: WeakInertial

# Output specification
Observers:
  VolumeFileName: "<filename>"
  ReductionFileName: "<filename>"

# Nonlinear solver (Newton-Raphson with globalization)
NonlinearSolver:
  NewtonRaphson:
    ConvergenceCriteria:
      MaxIterations: <int>
      RelativeResidual: <double>
      AbsoluteResidual: <double>
    SufficientDecrease: <double>
    MaxGlobalizationSteps: <int>
    DampingFactor: <double>
    Verbosity: <Verbose|Quiet|Silent>

# Linear solver (GMRES + Multigrid + Schwarz smoother)
LinearSolver:
  Gmres:
    ConvergenceCriteria:
      MaxIterations: <int>
      RelativeResidual: <double>
      AbsoluteResidual: <double>
    Verbosity: <Quiet|Silent>
  Multigrid:
    Iterations: <int>
    InitialCoarseLevels: Auto
    PreSmoothing: <bool>
    UseBottomSolver: <bool>
    PostSmoothingAtBottom: <bool>
    Verbosity: <Silent>
    OutputVolumeData: <bool>
  SchwarzSmoother:
    MaxOverlap: <int>
    Iterations: <int>
    Verbosity: <Silent>
    SubdomainSolver:
      Gmres: <...>
    # ...more config

# Radial coordinate compression
RadiallyCompressedCoordinates:
  InnerRadius: <double>
  OuterRadius: <double>
  Compression: Inverse

# Events and triggers at each iteration
EventsAndTriggersAtIterations:
  - Trigger: Always | HasConverged | <other>
    Events:
      - ObserveNorms:
          SubfileName: <name>
          TensorsToObserve:
            - Name: <VariableName>
              NormType: Max | L2Norm | Min
              Components: Individual
      - ObserveFields:
          SubfileName: <name>
          VariablesToObserve:
            - ConformalFactor
            - Lapse
            - Shift
            - SpatialMetric
            - ExtrinsicCurvature
            - HamiltonianConstraint
            - MomentumConstraint
          # ... more variables
          InterpolateToMesh: None
          CoordinatesFloatingPointType: Double
          FloatingPointTypes: [Double]
          BlocksToObserve: All

# Adaptive Mesh Refinement
Amr:
  Verbosity: Quiet
  Criteria: []
  Policies: { ... }
  MaxCoarseLevels: Auto
  Iterations: <int>

PhaseChangeAndTriggers: []

BuildMatrix:
  MatrixSubfileName: None
  Verbosity: Verbose
  EnableDirectSolve: <bool>
  SkipResets: <bool>
```

---

### 3. AVAILABLE ANALYTIC DATA AND SOLUTION CLASSES

The SolveXcts executable has access to the following analytic-data and solution classes (from `SolveXcts.hpp` lines 126-131):

**Analytic Solutions (XCTS::Solutions):**
1. **Flatness** — Flat spacetime (conformal factor ψ=1, all sources zero)
2. **WrappedGr<KerrSchild>** — Kerr-Schild metric wrapped for XCTS
3. **WrappedGr<SphericalKerrSchild>** — Spherical Kerr-Schild
4. **WrappedGr<HarmonicSchwarzschild>** — Harmonic Schwarzschild
5. **Schwarzschild** — Direct XCTS Schwarzschild implementation
6. **TovStar** — Tolman-Oppenheimer-Volkoff star solution
7. **WrappedGrMhd<RotatingStar>** — Rotating star with GRMHD
8. **WrappedGrMhd<CcsnCollapse>** — CCSN collapse
9. **WrappedGrMhd<MagnetizedTovStar>** — Magnetized TOV star

**Analytic Data (XCTS::AnalyticData):**
1. **Binary<AnalyticSolution, all_analytic_solutions>** — Binary superposition with window functions
2. **DuckTovStar** — Duck-shaped deformed TOV star
3. **NumericData** — Numeric/external data

**Key Implementation Files:**
- `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticSolutions/Xcts/Factory.hpp` (lines 21-34)
- `/data/haiyangw/nr/spectre/src/Elliptic/Executables/Xcts/SolveXcts.hpp` (lines 126-131)

---

### 4. TEUKOLSKY GRAVITATIONAL WAVE SUPPORT: CRITICAL FINDING

**Result: NO NATIVE TEUKOLSKY CLASS IN XCTS**

Search of the entire SpECTRE codebase reveals:

1. **Teukolsky exists only in CCE system:**
   - `/data/haiyangw/nr/spectre/src/Evolution/Systems/Cce/AnalyticSolutions/TeukolskyWave.hpp`
   - `/data/haiyangw/nr/spectre/src/Evolution/Systems/Cce/AnalyticSolutions/TeukolskyWave.cpp`
   - Used for characteristic extraction code (not XCTS)

2. **GR solutions available in SpECTRE do NOT include Teukolsky:**
   - `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticSolutions/GeneralRelativity/` contains only:
     - GaugeWave (coordinate wave in flat spacetime)
     - GaugePlaneWave
     - KerrSchild
     - Schwarzschild variants
     - Minkowski

3. **Wave-related analytic data:** Only for wave equation systems, not GR constraint equations.

4. **Brill wave:** Not found anywhere in the codebase.

---

### 5. HOW XCTS FREE DATA IS SPECIFIED

XCTS free data (the conformal metric $\bar{\gamma}_{ij}$ and trace-free extrinsic curvature $\bar{A}^{ij}$) is specified through **two mechanisms**:

#### Mechanism A: Background Field from Analytic Class
The `Background:` YAML key instantiates an `elliptic::analytic_data::Background` class that provides:

**Key Conformal Variables** (from `CommonVariables.hpp` lines 22-54):
- `Tags::ConformalMetric<DataType, 3, Frame::Inertial>` — The conformal spatial metric $\bar{\gamma}_{ij}$
- `gr::Tags::TraceExtrinsicCurvature<DataType>` — The trace $K$ (not the trace-free part)
- `::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataType>>` — Time derivative $\partial_t K$
- `Tags::ShiftBackground<DataType, 3, Frame::Inertial>` — Background shift $\beta^i_\mathrm{background}$
- Derivatives of the above

**Trace-free extrinsic curvature calculation:**
The trace-free part $\bar{A}^{ij}$ is computed from the background shift using:
```
Ā^ij = (1/2α̃)[L̄β - ū^ij]
```
where:
- $L̄\beta$ = conformal Killing operator (longitudinal operator applied to shift)
- $\bar{u}^{ij}$ = `LongitudinalShiftBackgroundMinusDtConformalMetric` (another background field)

**See:** 
- `/data/haiyangw/nr/spectre/src/Elliptic/Systems/Xcts/FirstOrderSystem.hpp` (lines 154-210) lists all background_fields
- `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/CommonVariables.hpp` (lines 22-54) defines common_tags

#### Mechanism B: Binary Superposition
The `Binary` class (lines 308-421 in Binary.hpp) superposes two isolated objects:
```yaml
Background:
  Binary:
    XCoords: [-5., 5.]
    CenterOfMassOffset: [0., 0.]
    ObjectLeft: <IsolatedObjectBase>    # e.g., WrappedGr<KerrSchild>
    ObjectRight: <IsolatedObjectBase>   # e.g., WrappedGr<KerrSchild>
    AngularVelocity: <Omega>
    Expansion: <adot0>
    LinearVelocity: [vx, vy

## Scout: AthenaK import

Perfect! Now I have comprehensive information. Let me create the report.

## Summary Report: AthenaK SpECTRE Initial Data Ingestion via Interpolation

### 1. Pgens Reading External Files

I found **4 pgens that read external initial data files**:

1. **`z4c_spectre_bbh.cpp`** - SpECTRE BBH volume HDF5 files
   - Path: `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_bbh.cpp`
   - Reads: HDF5 volume data (BbhVolume*.h5)
   - Key function: `LoadSpectreInitialData()`

2. **`kadath_bns.cpp`** - Kadath (FUKa) binary neutron star data
   - Path: `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/kadath_bns.cpp`
   - Reads: `.info` config file + `.dat` spectral coefficient file
   - Inline spectral interpolation from Kadath

3. **`lorene_bns.cpp`** - LORENE binary neutron star data
   - Path: `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/lorene_bns.cpp`
   - Reads: LORENE binary files
   - Query interface via C++ library

4. **`elliptica_bns.cpp`** - Elliptica BNS/BHNS data
   - Path: `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/elliptica_bns.cpp`
   - Reads: Elliptica format files (requires `initial_data_file` parameter)
   - Uses: `elliptica_id_reader_lib.h` C library

---

### 2. Existing SpECTRE Integration: `z4c_spectre_bbh.cpp`

**YES** - There is an existing production pgen for SpECTRE XCTS output.

**File path:** `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_bbh.cpp` (227 lines)

**Input file format:**
```
<problem>
pgen_name              = z4c_spectre_bbh
id_filename_glob       = /path/to/data/BbhVolume*.h5
id_subfile_name        = VolumeData
id_observation_step    = -1  # -1 loads last step
```

**Example athinput:** `/data/haiyangw/claude/z4c-CMM/athenak/inputs/z4c/spectre_bbh/z4c_spectre_bbh.athinput`

**SpECTRE Variable Names** (exact names read from H5 files):
- `SpatialMetric_xx`, `SpatialMetric_yx`, `SpatialMetric_yy`
- `SpatialMetric_zx`, `SpatialMetric_zy`, `SpatialMetric_zz`
- `ExtrinsicCurvature_xx`, `ExtrinsicCurvature_yx`, `ExtrinsicCurvature_yy`
- `ExtrinsicCurvature_zx`, `ExtrinsicCurvature_zy`, `ExtrinsicCurvature_zz`

(Lines 163-170 in z4c_spectre_bbh.cpp list these exact variable names)

**Interpolation Pipeline:**

1. **SpECTRE Exporter Library** (line 27):
   ```cpp
   #include <spectre/Exporter.hpp>
   ```
   Uses `spectre::Exporter::interpolate_to_points<3>()` for cubic interpolation from SpECTRE mesh to AthenaK cell centers.

2. **Per-meshblock workflow** (lines 140-220):
   - Collect target points (cell centers of AthenaK including ghost zones) in flat arrays
   - Call SpECTRE interpolator: `spectre::Exporter::interpolate_to_points<3>(...)`
   - Extract spatial metric (gxx, gyx, gyy, gzx, gzy, gzz) and extrinsic curvature components
   - Compute conformal factor: `psi4 = det(g_ij)^(1/3)` to enforce unit determinant on conformal metric (lines 206-218)
   - Copy to ADM variables on device

3. **Post-interpolation Z4c setup** (lines 64-91 in z4c_spectre_bbh.cpp):
   ```cpp
   pmbp->pz4c->GaugePreCollapsedLapse(pmbp, pin);  // set lapse from conformal factor
   pmbp->pz4c->ADMToZ4c<ng>(pmbp, pin);           // convert ADM to Z4c variables
   pmbp->pz4c->ADMConstraints<ng>(pmbp);          // compute constraint values
   ```

---

### 3. SpECTRE→AthenaK Import Procedure in hengrui_instructions.md

**File:** `/data/haiyangw/claude/4.1/multipole-cow/production-code/hengrui_instructions.md`

**Duck Star Implementation (TOV-based, not SpECTRE):**
- This document describes a duck-shaped TOV star pgen (analytic initial data + radial deformation)
- **Not a SpECTRE import procedure**, but illustrates the workflow for setting up initial data in AthenaK:
  1. Set spatial metric g_ij and extrinsic curvature K_ij
  2. Compute conformal factor (psi^4)
  3. Convert to Z4c variables via ADMToZ4c
  4. Run ADMConstraints to compute initial Hamiltonian/momentum violations
  5. Let Z4c constraint damping relax the violation

---

### 4. Teukolsky Pgen Structure (z4c_ccm_teukolsky.cpp)

**File:** `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_ccm_teukolsky.cpp` (262 lines)

**Analytic initial data setup** (lines 144-207):

```cpp
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  // Parameters from input file
  const Real X  = pin->GetOrAddReal("problem", "amp", 1e-5);
  const Real rc = pin->GetOrAddReal("problem", "r_c", 4.0);
  const Real tau = pin->GetOrAddReal("problem", "tau", 1.0);

  // Per-cell Kokkos parallel loop (lines 170-194)
  par_for("pgen_teukolsky", DevExeSpace(), 0, nmb-1, ksg, keg, jsg, jeg, isg, ieg,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Compute (x,y,z) coordinates
    Real x1v = CellCenterX(...);  // compute x, y, z
    
    // Evaluate Teukolsky metric analytically
    Real h[3][3], hdot[3][3];
    TeukolskyMetric(0.0, x1v, x2v, x3v, X, rc, tau, h, hdot);

    // Set ADM variables: g_ij = delta_ij + h_ij, K_ij = -1/2 dh_ij/dt
    adm.g_dd(m,a,b,k,j,i) = (a == b ? 1.0 : 0.0) + h[a][b];
    adm.vK_dd(m,a,b,k,j,i) = -0.5*hdot[a][b];
    
    // Gauge setup
    adm.alpha(m,k,j,i) = 1.0;
    z4c.alpha(m,k,j,i) = 1.0;
    z4c.vTheta(m,k,j,i) = 0.0;
  });

  // Convert ADM → Z4c variables
  pmbp->pz4c->ADMToZ4c<ng>(pmbp, pin);
  
  // Compute initial constraints
  pmbp->pz4c->ADMConstraints<ng>(pmbp);
}
```

**Hook point for imported ID:** Replace the `TeukolskyMetric()` computation with a call to SpECTRE interpolator within the Kokkos lambda, or pre-interpolate all points before the loop.

---

### 5. ADM and Z4c Variable Layout

**ADM Variables** (`/data/haiyangw/claude/z4c-CMM/athenak/src/coordinates/adm.hpp` lines 30-36):

```cpp
enum {
  I_ADM_GXX, I_ADM_GXY, I_ADM_GXZ, I_ADM_GYY, I_ADM_GYZ, I_ADM_GZZ,  // spatial metric
  I_ADM_KXX, I_ADM_KXY, I_ADM_KXZ, I_ADM_KYY, I_ADM_KYZ, I_ADM_KZZ,  // extrinsic curvature
  I_ADM_PSI4,                    // conformal factor
  I_ADM_ALPHA, I_ADM_BETAX, I_ADM_BETAY, I_ADM_BETAZ,  // lapse & shift
  nadm
};
```

**ADM Tensor structs** (lines 40-55):
```cpp
struct ADM_vars {
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> alpha;     // lapse
  AthenaTensor<Real, TensorSymm::NONE, 3, 1> beta_u;    // shift
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> psi4;      // conformal factor
  AthenaTensor<Real, TensorSymm::SYM2, 3, 2> g_dd;      // spatial metric
  AthenaTensor<Real, TensorSymm::SYM2, 3, 2> vK_dd;     // extrinsic curvature
};
```

**Z4c Variables** (`/data/haiyangw/claude/z4c-CMM/athenak/src/z4c/z4c.hpp` lines 44-54):
```cpp
enum {
  I_Z4C_CHI,                         // conformal factor (chi = psi^-2 or psi^4)
  I_Z4C_GXX, I_Z4C_GXY, I_Z4C_GXZ, I_Z4C_GYY, I_Z4C_GYZ, I_Z4C_GZZ,  // conformal metric
  I_Z4C_KHAT,                        // trace of extrinsic curvature
  I_Z4C_AXX, I_Z4C_AXY, I_Z4C_AXZ, I_Z4C_AYY, I_Z4C_AYZ, I_Z4C_AZZ,  // traceless extr. curv.
  I_Z4C_GAMX, I_Z4C_GAMY, I_Z4C_GAMZ,  // BSSN Gamma functions
  I_Z4C_THETA,                       // Z4c Theta variable
  I_Z4C_ALPHA, I_Z4C_BETAX, I_Z4C_BETAY, I_Z4C_BETAZ,  // lapse & shift
  nz4c
};
```

**Z4c Tensor structs** (lines 109-119):
```cpp
struct Z4c_vars {
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> chi;       // conformal factor
  AthenaTensor<Real, TensorSymm::SYM2, 3, 2> g_dd;      // conformal metric
  AthenaTensor<Real, TensorSymm::SYM2, 3, 2> vA_dd;     // traceless extrinsic curvature
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> vKhat;     // trace of extrinsic curvature
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> vTheta;    // Z4c constraint var
  AthenaTensor<Real, TensorSymm::NONE, 3, 1> vGam_u;    // Gamma functions
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> alpha;     // lapse
  AthenaTensor<Real, TensorSymm::NONE, 3, 1> beta_u;    // shift
};
```

---

## Concrete Import Path for SpECTRE XCTS → AthenaK Z4c

### Existing Mechanism (Production-Ready)

**Use `z4c_spectre_bbh.cpp` directly:**

1. Generate SpECTRE XCTS solution → HDF5 volume files (BbhVolume*.h5)
2. Create athinput file:
   ```
   <problem>
 

## Scout: CCM X=2 spec + two-pulse gate

Perfect. Now I have the complete picture. Let me compile the comprehensive specification.

## EXACT SPECIFICATION OF THE X=2 TEUKOLSKY-WAVE TEST FROM THE CCM PAPER

### 1. TEUKOLSKY-WAVE FREE DATA: EXACT RADIAL FUNCTIONS (A, B, C), GAUSSIAN PROFILE, AND PARAMETERS

**Metric form (eq:Teukolsky_metric, lines 1194-1198):**
```latex
ds^2 = -dt'^2 + (1+Af_{r'r'})dr'^2 + 2Bf_{r'\theta'}r'dr'd\theta'
     + 2Bf_{r'\phi'}r'\sin\theta'dr'd\phi'
     + (1+Cf^{(1)}_{\theta'\theta'}+Af^{(2)}_{\theta'\theta'})r'^2 d\theta'^2
     + 2(A-2C)f_{\theta'\phi'}r'^2 \sin\theta' d\theta' d\phi'
     + (1+Cf^{(1)}_{\phi'\phi'}+Af^{(2)}_{\phi'\phi'})r'^2 \sin^2\theta' d\phi'^2
```

**Radial amplitude functions A, B, C (eq:Teukolsky_ABC, lines 1202-1207):**
```latex
A = 3[F^(2)/r'^3 + 3F^(1)/r'^4 + 3F/r'^5],
B = -[F^(3)/r'^2 + 3F^(2)/r'^3 + 6F^(1)/r'^4 + 6F/r'^5],
C = (1/4)[F^(4)/r' + 2F^(3)/r'^2 + 9F^(2)/r'^3 + 21F^(1)/r'^4 + 21F/r'^5]
```

**Gaussian pulse profile (eq:Gaussian_pulse_F, lines 1224-1225):**
```latex
F(u') = X exp(-(u' - r'_c)^2/\tau^2)
```

**where u' = t' - r' is retarded time.**

**Retarded-time derivatives (eq:Teukolsky_wave_outgoing, lines 1229-1230):**
```latex
F^(n) ≡ [d^n F(u')/du'^n]_{u'=t'-r'}
```

**Angular basis functions (eq:Teukolsky_angular_basis, lines 1210-1217):**
- f_{r'r'} = 4√(π/5) Y_{20}(θ', φ')
- f_{r'θ'} = 2√(π/5) ∂_{θ'} Y_{20}(θ', φ')
- f_{r'φ'} = 0
- f^{(2)}_{θ'θ'} = -1
- f_{θ'φ'} = 0
- f^{(1)}_{θ'θ'} = 2√(π/5)[∂²_{θ'} - cot(θ')∂_{θ'} - ∂²_{φ'}/sin²(θ')] Y_{20}(θ', φ')
- f^{(1)}_{φ'φ'} = -f^{(1)}_{θ'θ'}
- f^{(2)}_{φ'φ'} = 1 - f_{r'r'}

**Spherical harmonic Y_{20} (lines 1219-1221):**
```latex
Y_{20} = (1/8)√(5/π)(1 + 3cos(2θ'))
```

**X=2 PARAMETERS (lines 1240, 1293):**
- **X = 2** (amplitude)
- **r'_c = 20** (center of pulse at t'=0)
- **τ = 2** (width)

---

### 2. INITIAL-DATA CONSTRUCTION: XCTS FORMULATION, CONFORMAL METRIC, TRACE-FREE EXTRINSIC CURVATURE, BOUNDARY CONDITIONS

**Initial data method (lines 1192-1193):**
> "The initial data of the Cauchy system is constructed utilizing the Extended Conformal Thin Sandwich (XCTS) formulation [York:1998hy, Pfeiffer:2002iy], which accounts for nonlinear effects that arise when the amplitude of the Teukolsky wave is large."

**XCTS solving statement (line 1193):** The paper states XCTS is "solved" but does **not detail** the conformal metric decomposition, trace-free extrinsic curvature Ã^{ij}, or XCTS boundary conditions in the Teukolsky-wave test section. These are referenced only generically as the formulation used.

**Boundary conditions and lapse/shift:**
- The paper follows the GH (generalized harmonic) formulation for Cauchy evolution (lines 204-220).
- The boundary conditions are imposed via the **Bjørhus method** (line 204).
- For the Teukolsky-wave free data, the **metric-perturbation form h_{ij}** is constructed directly from the (A, B, C) radial functions as part of the analytically specified perturbative metric (eq:Teukolsky_metric).
- The lapse and shift are **evolved dynamically** as part of the GH system; they are not explicitly specified in the test description.
- **Gauge:** GH harmonic gauge conditions are used; the lapse α and shift β^i are evolved according to the gauge source function H^μ (lines 1450 in commented text).

---

### 3. NUMERICAL SETUP FOR X=2 TEST

**Worldtube and domain parameters (line 1240, Fig. 1105-1106 caption):**
- **Cauchy outer boundary: r'_{out} = 41** (coincides with worldtube for CCM/CCE)
- **Reference system outer boundary: r'_{ref} = 900** (SpEC code, causally disconnected)
- **Worldtube location for CCE wave extraction: r'_{wt} = 41** (same as Cauchy boundary)

**Evolution parameters (lines 1250, 1293-1295):**
- **Evolution time: spanning over 1000 code units** (line 1250)
  - Specifically: "The duration of our simulation greatly exceeds the timescale of the physical process of interest, which involves the propagation of the Teukolsky wave from its initial location to null infinity within the first 50 code units."
- **Numerical resolutions: three different resolutions** (lines 1235, 1250, 1289)
  - The paper mentions "three numerical resolutions" but does not explicitly state the grid point counts.

**Gauge:** GH (generalized harmonic) evolution with dynamical lapse and shift.

**Extraction radii:** Waveforms extracted at:
- **Future null infinity ℐ^+** (characteristic grid, via CCE/CCM)
- **Worldtube r'_{out} = 41** (for Cauchy-characteristic matching)

---

### 4. THE GATE — TWO-PULSE FEATURE: TIMING ANALYSIS AND CCM TRANSPARENCY

**Definition of the waveform feature (lines 1312-1313):**
The "two-pulse feature" refers to the **primary outgoing Teukolsky pulse AND the secondary junk-radiation pulse** observed in **ψ₃^{(l=2,m=0)} and strain h₂₀ at future null infinity ℐ^+**.

**PRIMARY PULSE (main GW signal):**
> "The primary outgoing Teukolsky pulse, initially located at r'_c (Fig. 1), reaches the Cauchy outer boundary r'_{out} after a time interval of r'_{out} - r'_c = 21. Because of the null slicing of the characteristic system, any outgoing null GW is **instantaneously** transmitted to future null infinity as soon as it intersects the worldtube. Consequently, we observe that the main pulse also appears at ℐ^+ at a time of **t = 21**."

**SECONDARY PULSE (junk radiation):**
> "Since our constructed NR initial data is not a perfect error-free solution for a solely outgoing Teukolsky wave, there is a wave component that travels inward, commonly known as junk radiation. As the evolution progresses, this component falls towards the coordinate center, crosses it, and bounces back out. The junk wave eventually reaches the Cauchy outer boundary (which coincides with the worldtube) and ℐ^+ at a time of r'_{out} + r'_c = 61. This junk wave can be seen as a **secondary pulse at time 61** in the top left panel of Fig. large_amp_strain_only."

**TERTIARY PULSE (spurious reflection at boundary):**
> "Moreover, we find the appearance of a **tertiary wave at a time of 100**. This arises because a portion of the primary Teukolsky wave **reflects off the numerical boundary between the Cauchy and characteristic systems**. Subsequently, this reflected wave traverses the entire Cauchy region, with a propagation time of 2r'_{out}, namely the domain's diameter. At t = 2r'_{out} + (r'_{out} - r'_c) = **103**, the reflected wave escapes to future null infinity."

**Numerical values for X=2, r'_c = 20, r'_{out} = 41:**
- Main pulse arrival at ℐ^+: **t = r'_{out} - r'_c = 41 - 20 = 21**
- Junk/secondary pulse at ℐ^+: **t = r'_{out} + r'_c = 41 + 20 = 61**
- Tertiary (reflected) pulse at ℐ^+: **t = 2r'_{out} + (r'_{out} - r'_c) = 2(41) + 21 = 103**

**CCM vs. CCE/reflective-boundary distinction (lines 1315-1317):**
> "The presence of the tertiary reflected wave indicates that our current CCM algorithm cannot entirely eliminate the spurious reflection at the numerical boundary. We attribute this issue to the gauge part of the boundary conditions... Nevertheless, upon comparing the results obtained using CCE and CCM to the reference ones (in blue) in Fig. large_amp_strain_only, we find that the new physical boundary conditions already lead to noticeable improvements. Specifically, the differences between the waveforms shown in the bottom panels... clearly demonstrate that the **CCM method systematically reduces the deviation from the reference system by approximately one order of magnitude**."

**Bondi-gauge diagnostics for X=2 (lines 1320-1333, equations eq:bondi_violation, eq:im_psi2):**
The paper uses waveform-quality constraints to assess CCM performance:
```latex
C_{ψ_4} ≡ ψ_4 + ḧ = 0
C_{ψ_3} ≡ ψ_3 - (1/√2) ð ḣ = 0
C_{ψ_s} ≡ ψ̇_s + (1/√2) ð ψ_{s+1} - (3-s)/4 h̄ ψ_{s+2} = 0  (s=0,1,2)
C_{Im ψ_2} ≡ Im ψ_2 + Im[(1/2)ð²h + (1/4)h̄ḣ] = 0
```
> "We see the CCM method systematically reduces the constraint violations compared to CCE, with the exception of the constraints... for which CCM and CCE are similar."

---

### SUMMARY: CCM TRANSPARENCY WITH NO IMPLEMENTATION BUG

The X=2 Teukolsky-wave test demonstrates **CCM transparency** through:

1. **Amplitude:** X=2 (nonperturbative regime requiring XCTS)
2. **Free data:** Gaussian-profile Teukolsky wave centered at r'_c=20, width τ=2
3. **Domain:** Cauchy outer boundary r'_{out}=41 (close enough to reveal backscatter)
4. **Evolution:** Full nonlinear GH dynamics over 1000 code units
5. **Key waveform features at ℐ^+ (rpsi₄ ~ F^{(6)}, News ~ F^{(5)}):**
   - **t≈21:** Primary outgoing pulse (instantaneously transmitted via worldtube)
   - **t≈61:** Secondary junk-radiation pulse (inward → inward-to-outward bounce)
   - **t≈100-103:** Tertiary reflected pulse (boundary reflection; gauge-condition artifact)
6. **CCM advantage:** ~1 order of magnitude reduction in deviation from SpEC reference, Bondi-gauge-constraint suppression, demonstrating physical boundary conditions encode backscatter correctly.

**The two-pulse feature in the rpsi₄/News waveform at ℐ^+ (primary at t=21, secondary at t=61) 

## Scout: mission management

Excellent. Now I have all the context I need. Let me compile the comprehensive proposal:

## DESIGN PROPOSAL: Mission 6 Management Scaffolding

### PART 1: RALPH-LOOP STATE STRUCTURE AND PROGRESS TRACKING

#### 1.1 Ralph-Loop State File Format (`.claude/ralph-loop.local.md`)

The loop state file is a **nine-field frontmatter contract** read by `ralph_stop_guard.sh` and `loop_gate.py`:

```yaml
---
active: true                    # Boolean: loop disarm/rearm flag
iteration: 0                    # Integer: current iteration counter
iteration_base: 0               # Integer: mission's starting iteration (for multi-mission repos)
session_id: "{UNIQUE_UUID}"     # String: mission ID; archive on change (loop_gate auto-archives)
max_iterations: 1000            # Integer: hard iteration cap
no_progress_limit: 8            # Integer: halt after N iterations without SOLID progress
stuck_counter_limit: 3          # Integer: halt after N charged same-iteration stop fires
stuck_fire_grace: 1             # Integer: tolerated extra stop fires per iteration
paper: "z4c-CCM"                # String: paper/project code (ledger key)
completion_promise: "TAG"       # String: unambiguous completion token (all-caps)
completion_status: ""           # EMPTY while running; "promise-satisfied" on exit
started_at: "2026-06-13T..."    # UTC ISO 8601; mission start timestamp
---
```

**Key contracts:**
- `session_id` MUST change between missions. When `loop_gate.py decide` detects a change, it **archives** the previous mission's gate state to `.claude/runs/{old-session-id}.gate_state.json` and resets counters.
- `completion_promise` and the `<promise>` tag at the bottom of the file MUST match exactly; the gate detects mismatch as a halt error.
- `completion_status` MUST be empty while `active: true`. At disarm (final commit), the disarming agent sets `completion_status: promise-satisfied` and `active: false` together; the gate then reports `halt:completed` and freezes the record.
- All loop-state changes **commit to git** (no mutation of the file outside the ralph-loop protocol).

**Progress signal (loop_gate.py):**
The gate monitors three signals aggregated from **knowledge-database** and **result-database** JSONL ledgers:
- **SOLID knowledge nodes** (`status: "solid"` in `knowledge-database/paper_z4c-CCM/nodes.jsonl`)
- **Admitted/classified results** (`status: "admitted" | "classified"` in `result-database/paper_z4c-CCM/results.jsonl`)
- **Discharged errors** (`status: "discharged"` in `error-database/paper_z4c-CCM/trials.jsonl`)

If none of these advance for `no_progress_limit` consecutive iterations, the gate halts with `halt:no_progress` and writes `.claude/HUMAN_REVIEW_REQUIRED.md`.

#### 1.2 Self-Injection Cadence and Protocol Re-emission

Every `self_inject_interval` iterations (default 5), the Stop-hook block (`ralph_stop_guard.sh`) re-emits the binding protocol into the loop reason:
- `.claude/ralph-loop.local.md` (the state file itself)
- `phys-agentic-loop/alignment.md` (always-on agent discipline)
- `phys-agentic-loop/_common/contracts/research_admission_contract.md` (plain-language result gate)

If any of these files cannot be resolved or re-injected, the count increments. After `self_inject_fail_limit` (default 3) consecutive failed windows, the gate issues `escalation_required`, forcing a pipelines/6-escalation run immediately.

---

### PART 2: MISSION 6 SUBSTAGE DECOMPOSITION

**Mission 6 Directive**: "SpECTRE XCTS solve of the paper's X=2 Teukolsky ID → AthenaK import/interpolation → large-domain (r>100) X=2 CCM evolution → exact two-pulse-feature gate"

The mission decomposes into **ordered substages** (each is one iteration per plan-edit-verify-commit cycle):

#### S6.0: XCTS Exact ID Generation (gateway)
**Claim:** Produce a numerically exact X=2 Teukolsky constraint-solved initial data (XCTS) in SpECTRE native format (5 metric + extrinsic curvature components, HDF5 volumetric).

**Gate:** SpECTRE XCTS solver converges to constraint violations ≤ 1e-10 (uniform norm) on a 64³ pilot grid; intermediate-domain convergence verified.

**Evidence type:** Numerical convergence study + SpECTRE solver log + HDF5 metadata.

**Artifact path:** `results/numerical/xcts_x2_exact/xcts_64_final.h5` + convergence table `xcts_64_convergence.txt`.

**Verifier:** `python scripts/validate_xcts.py results/numerical/xcts_x2_exact/xcts_64_final.h5 --constraint-norm 1e-10 --report` → outputs `xcts_validation.txt` with constraint norms and grid metadata.

**Dependencies:** SpECTRE CharacteristicExtract binary (ref-code/spectre-cce), Python scri library.

**Open obligation after this:** O-M6-1 (XCTS→AthenaK metric interpretation).

---

#### S6.1: Metric Interpolator (SpECTRE HDF5 → AthenaK Cartesian)
**Claim:** Build a three-map interpolator chain: (1) SpECTRE (global spectral) → nodal Cartesian grid sampler; (2) HDF5 read + 3D Lagrange interp to AthenaK mesh points; (3) inverse-shift + lapse reconstruction (ℓ→α via ℓ²(r) exact formula at t=0).

**Gate:** Interpolated metric at 10 test points matches the source XCTS to relative error ≤ 1e-6 (infinity norm). Inverse-shift reconstruction yields α within 1e-5 of the ℓ-based formula.

**Evidence type:** Unit test + roundtrip validation script + interpolation error table.

**Artifact path:** `packages/zccm/xcts_interp.jl` (native Julia, reuses existing lapse-solving infrastructure from ZccmJl).

**Verifier:** `julia --project=packages/ZccmJl -e "include(\"packages/zccm/test_xcts_interp.jl\"); test_xcts_interpolation();"` → outputs `xcts_interp_validation.txt` with 3 error norms.

**Dependencies:** Julia HDF5.jl, Interpolations.jl, existing ZccmJl lapse-solving suite.

**Open obligation after this:** O-M6-2 (AthenaK restart-file format generation).

---

#### S6.2: AthenaK Restart-File Generator
**Claim:** Write the interpolated metric + extrinsic curvature into AthenaK's native HDF5 restart format, with proper volumetric field grouping, ghost-cell padding, and mesh-block top
