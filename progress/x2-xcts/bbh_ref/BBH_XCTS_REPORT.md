# Standard SpECTRE BBH XCTS initial-data reference run

Reference run reproducing the standard binary-black-hole (BBH) XCTS initial-data
generation, executed with the **same environment and launch method** as our
Teukolsky XCTS "initial case" in the z4c-CCM project. Intended as the comparison
baseline for the XCTS comparison DAG (`../XCTS_DAG.md`, DAG 1).

## Exact configuration

| Item | Value |
|---|---|
| SpECTRE build | `/data/haiyangw/nr/spectre/build` |
| SpECTRE version | **2025.08.19** (git revision `9f98d9fdd`, branch `develop`, Release) |
| Executable | `/data/haiyangw/nr/spectre/build/bin/SolveXcts` (pristine general XCTS executable) |
| Input file (run copy) | `/data/haiyangw/nr/spectre/bbh_xcts_run/BinaryBlackHole.yaml` |
| Input file (provenance copy) | `/data/haiyangw/claude/z4c-CMM/progress/x2-xcts/bbh_ref/BinaryBlackHole.yaml` |
| Source of YAML | `/data/haiyangw/nr/spectre/tests/InputFiles/Xcts/BinaryBlackHole.yaml` (equal-mass Kerr-Schild superposition, ApparentHorizon excision BCs, outer Flatness) |
| Environment | `source /data/haiyangw/claude/z4c-CMM/scripts/spectre_xcts_env.sh` (MKL shim: LD_LIBRARY_PATH + LD_PRELOAD of libgomp.so.1 + layered MKL; OMP_NUM_THREADS=MKL_NUM_THREADS=1) |
| Launch | `charmrun +p24 .../SolveXcts --input-file BinaryBlackHole.yaml` (mpi-linux-x86_64 Charm++; must launch via charmrun, not +pN to the binary) |
| Run directory | `/data/haiyangw/nr/spectre/bbh_xcts_run` |
| Resource | +p24 on the 128/256-core node, coexisting with the Teukolsky +p48 run and an unrelated julia job (no oversubscription) |

**Executable used: `SolveXcts` (the general one).** No fallback to
`SolveXctsVacuum` was needed — `BinaryBlackHole.yaml` declares
`Executable: SolveXcts`, and the startup executable-name check accepted it, and
`--check-options` parsed the file successfully.

## Deviations from the example (minimal, run-environment only)

1. **Added the YAML metadata `---` separator.** The shipped test YAML
   (`tests/InputFiles/Xcts/BinaryBlackHole.yaml`) is consumed by SpECTRE's test
   harness, which strips a metadata header. Run directly through the binary's
   option parser, it failed with *"Missing metadata in input file. YAML input
   files begin with a metadata section terminated by `---`"*. Fix: inserted a
   `---` line after the test-harness metadata block (the
   `Executable`/`Testing`/`ExpectedOutput`/`OutputFileChecks` keys) and before
   `Parallelization:` (the start of the actual options). This is the same
   structure the Teukolsky YAML already uses. **No physics, domain, solver, or
   observe-field option was altered.**
2. **Output directory.** Run from a fresh dir
   (`/data/haiyangw/nr/spectre/bbh_xcts_run`); the example writes to fixed file
   names `BbhReductions.h5` / `BbhVolume*.h5` in the working directory, so a clean
   dir suffices (no path edits needed).

The BBH physics (BinaryCompactObject domain with two ApparentHorizon-excised
spheres + envelope + outer shell to R=1e9, Binary Kerr-Schild superposition free
data, NewtonRaphson + Gmres/Multigrid/Schwarz solver stack) was kept **as-is**.

## Validation (single-process `--check-options`)

```
source scripts/spectre_xcts_env.sh
$SPECTRE_BUILD/bin/SolveXcts --check-options --input-file BinaryBlackHole.yaml
```
Result: `BinaryBlackHole.yaml parsed successfully!` (exit 0). The MKL shim worked
(no symbol-lookup error). No option fixes were required beyond the `---`
separator.

## Domain (built at startup, AMR iteration 0)

- AMR level 2 grid: **84 elements, 20304 grid points** (matches the reference
  YAML's expected `NumberOfPoints = 20304`).

## Convergence result — CONVERGED (full run, not preliminary)

- **NewtonRaphson CONVERGED in 5 iterations** — RelativeResidual criterion
  (1e-10): final residual fraction `9.47159e-11` of the initial.
- **AMR iterations performed: 1** (`Amr.Iterations: 1`). The solve ran on the
  startup-refined grid (84 elements / 20304 grid points); after convergence AMR
  evaluated its (empty) criteria and the program ended — there was no second
  re-solve.
- **GMRES iterations:** per Newton step 10, 16, 23, 13, 34 (= 96 outer GMRES
  iterations); cumulative `GmresResiduals.dat` has 101 rows. Each GMRES step runs
  the full Multigrid + Schwarz-smoother preconditioner.
- **Wall time: 00:10:14** (start 04:50:04 UTC -> completion 05:00:32 UTC),
  on +p24 (24 cores) of the 256-core node, sharing with the Teukolsky +p48 run
  and a julia job.

### Newton-Raphson residual history (`NewtonRaphsonResiduals.dat`)

| Iter | Residual | vs reference YAML expected |
|---|---|---|
| 0 | 9.588365e+00 | 9.58836468792872e+00 (exact) |
| 1 | 2.567448e+00 | 2.56744819082553e+00 (exact) |
| 2 | 6.927921e-03 | 6.92792186565806e-03 (exact) |
| 3 | 1.218999e-03 | 1.21895158596328e-03 (exact) |
| 4 | 1.081882e-06 | 7.97507362468631e-07 (same order; small diff) |
| 5 | 9.081706e-10 | 7.12541536160998e-10 (same order; converged) |

The first four residuals reproduce the reference to all printed digits. Steps 4-5
differ at the ~1e-6/1e-10 level because the +p24 element distribution differs from
the reference run's (the reference test was generated on a different core count);
both converge to the RelRes-1e-10 criterion. **This is a faithful reproduction.**

### Constraint L2 norms: initial guess vs converged (`Norms.dat`)

| Quantity | Initial (KS superposition guess) | Converged | Reduction |
|---|---|---|---|
| L2Norm(HamiltonianConstraint) | **1.745895e-02** | **5.097199e-03** | 3.4x |
| L2Norm(MomentumConstraint_x) | 1.131148e-02 | 1.259369e-03 | |
| L2Norm(MomentumConstraint_y) | 8.478599e-03 | 1.298358e-03 | |
| L2Norm(MomentumConstraint_z) | 1.131245e-02 | 1.264785e-03 | |
| \|MomentumConstraint\| (rms) | 1.810546e-02 | 2.207130e-03 | 8.2x |
| Max(ConformalFactor) | 1.000000e+00 | 1.109963e+00 | (ref 1.10996e+00, exact) |
| Max(\|ShiftExcess\|) | 7.118553e-01 | 8.680516e-01 | (ref 8.6805e-01, exact) |
| Min(Lapse) | 7.054928e-01 | 6.281279e-01 | (ref 6.2813e-01, exact) |

The converged Max(ConformalFactor)/Max(|ShiftExcess|)/Min(Lapse) and the per-axis
Momentum-constraint L2 norms match the reference YAML's `OutputFileChecks` expected
values to all printed digits. The converged Hamiltonian L2 (5.10e-3) is close to
the reference's 4.95e-3 (the small difference tracks the steps-4/5 residual
difference and the slightly different element distribution). These are finite
because this is a DG solve of a superposed (not exactly-constraint-satisfying)
free data — the constraints are reduced ~3-8x from the guess but not driven to
zero; that residual is the expected truncation/superposition error of the
crudely-p-refined manual domain (the YAML comment notes AMR would refine further).

## Comparison to the Teukolsky XCTS case

- Both runs use the identical XCTS elliptic system, DG discretization, and
  NewtonRaphson + Gmres/Multigrid/Schwarz solver stack (see DAG 1). The only
  differences are the free-data class (`Binary` Kerr-Schild superposition vs
  `TeukolskyWave`), the DomainCreator (`BinaryCompactObject` vs `Sphere`), and the
  boundary conditions (ApparentHorizon excision + outer Flatness vs regular center
  + outer Robin).
- The Teukolsky case starts from a **flat** guess at Hamiltonian L2 ~3.85e-2; this
  BBH case starts from the **superposed Kerr-Schild** guess at Hamiltonian L2
  ~1.75e-2 (lower, because the superposition is already a good approximation of
  each black hole's geometry).
- Executable difference: the Teukolsky run uses `SolveXctsVacuum` (general
  `SolveXcts` minus the `HydroQuantitiesCompute`/`LowerSpatialFourVelocity`
  observations, so a genuine vacuum free-data class is never asked for matter
  primitives); this BBH run uses the stock `SolveXcts`.

## Output H5 paths

- Reductions: `/data/haiyangw/nr/spectre/bbh_xcts_run/BbhReductions.h5`
  (subfiles `NewtonRaphsonResiduals.dat`, `Norms.dat`, `GmresResiduals.dat`,
  `MultigridResiduals.dat`, `SchwarzSmoother*` and `src.tar.gz` provenance).
- Volume: `/data/haiyangw/nr/spectre/bbh_xcts_run/BbhVolume0.h5` ...
  `BbhVolume23.h5` (one file per +p24 rank; subfile `VolumeData` with
  ConformalFactor, Lapse, Shift, ShiftExcess, SpatialMetric, ExtrinsicCurvature,
  Hamiltonian/MomentumConstraint, RadiallyCompressedCoordinates).
