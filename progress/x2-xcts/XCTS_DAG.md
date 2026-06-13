# SpECTRE XCTS initial-data dependency DAGs

Documentation deliverable for the z4c-CCM mission. Two dependency graphs:

1. **DAG 1 — BBH-vs-Teukolsky XCTS comparison.** What is *literally shared*
   between the standard binary-black-hole (BBH) XCTS solve and our new X=2
   Teukolsky-wave XCTS solve, and exactly where they diverge.
2. **DAG 2 — XCTS generation usage graph.** A practical "how to generate XCTS
   initial data with SpECTRE" pipeline, from authoring a free-data class to
   consuming the output in AthenaK.

Every node is annotated with the concrete file path and/or C++/YAML symbol so it
is directly actionable later. All paths are absolute.

Reference input files used:

- **BBH reference:** `/data/haiyangw/nr/spectre/tests/InputFiles/Xcts/BinaryBlackHole.yaml`
  (also the pipeline template `/data/haiyangw/nr/spectre/support/Pipelines/Bbh/InitialData.yaml`).
- **Teukolsky reference:** `/data/haiyangw/nr/spectre/SolveXctsTeukolskyX2.yaml`.

---

## The shared XCTS elliptic system (common to both)

Both problems solve the *same* extended-conformal-thin-sandwich (XCTS) elliptic
system. The C++ system type is identical in both executables:

```cpp
using system = Xcts::FirstOrderSystem<
    Xcts::Equations::HamiltonianLapseAndShift,   // solve psi, alpha*psi, beta^i
    Xcts::Geometry::Curved,                       // general conformal metric
    /*conformal_matter_scale=*/0>;
```

Defined in `/data/haiyangw/nr/spectre/src/Elliptic/Systems/Xcts/FirstOrderSystem.hpp`.

- **Primal (solved-for) fields** — `system::primal_fields`:
  - `Xcts::Tags::ConformalFactorMinusOne<DataVector>`               (psi - 1)
  - `Xcts::Tags::LapseTimesConformalFactorMinusOne<DataVector>`     (alpha*psi - 1)
  - `Xcts::Tags::ShiftExcess<DataVector, 3, Frame::Inertial>`       (beta^i_excess)
- **Free-data / background fields** the system demands (`system::background_fields`,
  curved + HamiltonianLapseAndShift specialization):
  - `Xcts::Tags::ConformalMetric<DataVector, 3, Frame::Inertial>`   (gammabar_ij)
  - `Xcts::Tags::InverseConformalMetric<...>`                       (gammabar^ij)
  - `::Tags::deriv<Xcts::Tags::ConformalMetric<...>>`               (d_k gammabar_ij)
  - `Xcts::Tags::ConformalChristoffel{First,Second}Kind`, `ConformalRicciTensor`,
    `ConformalRicciScalar`, `ConformalChristoffelContracted`
  - `gr::Tags::TraceExtrinsicCurvature<DataVector>`                 (K)
  - `::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataVector>>`     (dt K)
  - `::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataVector>>`  (d_i K)
  - `Xcts::Tags::ShiftBackground<DataVector, 3, Frame::Inertial>`   (beta^i_background)
  - `::Tags::deriv<Xcts::Tags::ShiftBackground<...>>`
  - `Xcts::Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<...>`
    ( (Lbar beta_bg)^ij - ubar^ij )  and its `::Tags::div<...>`
  - matter sources, scaled by conformal_matter_scale=0:
    `gr::Tags::Conformal<gr::Tags::EnergyDensity<DataVector>, 0>`,
    `gr::Tags::Conformal<gr::Tags::StressTrace<DataVector>, 0>`,
    `gr::Tags::Conformal<gr::Tags::MomentumDensity<DataVector, 3>, 0>`
- **Fluxes / sources computers** (identical):
  `Xcts::Fluxes<Eqs, Geom>`, `Xcts::Sources<Eqs, Geom, Scale>`,
  `Xcts::LinearizedSources<Eqs, Geom, Scale>`
  (`/data/haiyangw/nr/spectre/src/Elliptic/Systems/Xcts/Equations.hpp`).
- **ADM output** (computed identically via
  `Xcts::Tags::SpacetimeQuantitiesCompute`): gamma_ij = psi^4 gammabar_ij,
  K_ij, lapse alpha = (alpha*psi)/psi, shift beta^i = beta^i_bg + beta^i_excess,
  plus `HamiltonianConstraint` / `MomentumConstraint`.

The free data is the **only physics input**; everything downstream of it (the
discretization, the solver stack, the ADM reconstruction, the observers) is the
same code path for BBH and Teukolsky.

---

## DAG 1 — BBH vs Teukolsky XCTS (shared spine, branched leaves)

ASCII box diagram. The **center column is the shared spine** (identical C++ /
solver / observer code). The **left column** is BBH-specific, the **right
column** is Teukolsky-specific. `[I]` = identical component, `[S]` = swapped.

```
                  BBH branch                         SHARED SPINE                       Teukolsky branch
            (tests/.../BinaryBlackHole.yaml)      (same code both sides)            (SolveXctsTeukolskyX2.yaml)
 ============================================  =================================  ============================================

 FREE-DATA SOURCE [S]                                                             FREE-DATA SOURCE [S]
 +--------------------------------------+                                         +--------------------------------------+
 | Xcts::AnalyticData::Binary<...>      |                                         | Xcts::AnalyticData::TeukolskyWave    |
 | Binary.hpp                           |                                         | TeukolskyWave.hpp                    |
 | superpose 2x WrappedGr<KerrSchild>   |                                         | flat delta_ij + verified Teukolsky   |
 | (M=0.4229 each, q=1, a=0)            |                                         | h_ij(l=2,m=0,even); X=2,rc=20,tau=2  |
 +-------------------+------------------+                                         +------------------+-------------------+
                     |                                                                              |
 (a) gammabar_ij      = superposed boosted KS metric                              (a) gammabar_ij    = delta_ij + h_ij (trace-free, non-unimodular)
 (d) K                = nonzero (boosted KS slicing)         <--- DIFF (d) --->   (d) K              = 0  (and dt K = 0)
 (e) ubar=dt gammabar = 0 (quasi-equilibrium, conformally    <--- DIFF (e) --->   (e) Long.ShiftBg-DtConfMetric = -gammabar^ik gammabar^jl hdot_kl
                          stationary; Long.ShiftBg=Omega-driven)                       (propagating wave-momentum); ShiftBackground = 0
 sources             = vacuum 0                                                   sources           = vacuum 0  [I-in-value]
                     |                                                                              |
                     +-----------------------+                +-----------------------+------------+
                                             v                v
                                   +-------------------------------------------+
                                   | FREE DATA: gammabar, ubar, K, dtK, sources|   <-- the only physics input
                                   | (system::background_fields)               |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | XCTS ELLIPTIC SYSTEM                  [I]  |
                                   | Xcts::FirstOrderSystem<HamiltonianLapse-  |
                                   |   AndShift, Curved, 0>                     |
                                   | solve for psi, alpha*psi, beta^i_excess   |
                                   | Fluxes / Sources / LinearizedSources      |
                                   +---------------------+---------------------+
                                                         |
 BOUNDARY CONDITIONS [S]                                 |                       BOUNDARY CONDITIONS [S]
 +--------------------------------------+                |                       +--------------------------------------+
 | (c) excision = ApparentHorizon<Curv> |                |                       | (c) regular center (no inner BC);    |
 |     Theta=0 Neumann on psi;          |                |                       |     no excision                      |
 |     Lapse Dirichlet from Kerr;       |                |                       | outer = Robin<Eqs>                   |
 |     beta = Omega x r + a_dot r       |                |                       |   d_r(r psi)=0, d_r(r alpha psi)=0,  |
 | outer = Flatness (or Robin)          |                |                       |   d_r(r beta_excess)=0  (radiative)  |
 +-------------------+------------------+                |                       +------------------+-------------------+
                     +----------------------+            |            +----------+
                                            v            v            v
                                   +-------------------------------------------+
                                   | DOMAIN [S]                                |
                                   |  BBH: DomainCreator BinaryCompactObject   |
                                   |       (2 excised spheres + envelope +     |
                                   |        outer shell R=1e9)                 |
                                   |  Teuk: DomainCreator Sphere (InnerR=0.5,  |
                                   |        regular center fill, OuterR=1e9,   |
                                   |        RadialPartitioning [12,40,120] to  |
                                   |        resolve pulse at r=20)             |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | DG DISCRETIZATION                     [I] |
                                   | DiscontinuousGalerkin, WeakInertial,      |
                                   | Massive, GaussLobatto, Penalty=1          |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | NONLINEAR SOLVER                      [I] |
                                   | NewtonRaphson (RelRes 1e-9..1e-10,        |
                                   | SufficientDecrease 1e-4, globalization)   |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | LINEAR SOLVER STACK                   [I] |
                                   | Gmres + Multigrid + SchwarzSmoother       |
                                   | (subdomain Gmres + MinusLaplacian/        |
                                   |  ExplicitInverse preconditioner)          |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | ADM RECONSTRUCTION                    [I] |
                                   | SpacetimeQuantitiesCompute:               |
                                   | gamma_ij=psi^4 gammabar_ij, K_ij,         |
                                   | alpha, beta^i; H & M constraints          |
                                   +---------------------+---------------------+
                                                         |
                                                         v
                                   +-------------------------------------------+
                                   | OBSERVERS                             [I] |
                                   | dg::Events::field_observations +          |
                                   | ObserveNorms + ObserveAdmIntegrals        |
                                   +---------------------+---------------------+
                                                         |
                          +------------------------------+------------------------------+
                          v                                                             v
              +-----------------------------+                            +-----------------------------+
              | Reductions H5  [I-mechanism]|                            | Volume H5  [I-mechanism]    |
              | BbhReductions / XctsTeukol- |                            | BbhVolume / XctsTeukolskyX2 |
              | skyX2Reductions: NewtonRaph-|                            | Volume: ConformalFactor,    |
              | son residuals + Norms (H,M  |                            | Lapse, Shift, SpatialMetric,|
              | L2)                         |                            | ExtrinsicCurvature, H, M    |
              +-----------------------------+                            +-----------------------------+

 CONTROL LOOP [S, BBH-only]                                             CONTROL LOOP [S]
 +--------------------------------------+                               +--------------------------------------+
 | (f) BBH ID control system:           |                               | (f) NONE. Fixed analytic free data;  |
 |  PostprocessId / ControlId.py adjust |                               |     a single XCTS solve. Nothing to  |
 |  MassA/B, spins, CoM, ADM linear     |                               |     iterate (the wave is prescribed).|
 |  momentum, then re-solve XCTS        |                               |                                      |
 +--------------------------------------+                               +--------------------------------------+
```

**Executable note.** The Teukolsky run uses the new `SolveXctsVacuum`
executable (`/data/haiyangw/claude/z4c-CMM/progress/x2-xcts/spectre_src/SolveXctsVacuum.{hpp,cpp}`),
which is `SolveXcts` minus the `HydroQuantitiesCompute` observation tag (and the
`LowerSpatialFourVelocity` observation). The generic `SolveXcts` metavariables
instantiate `HydroQuantitiesCompute` for *every* registered background, forcing
each background to supply six hydro primitives; a genuine vacuum free-data class
like `TeukolskyWave` has no matter and must never be asked for them. The shipped
`SolveXctsTeukolskyX2.yaml` header still reads `Executable: SolveXcts` — update
it to `SolveXctsVacuum` to use the vacuum metavariables. The elliptic system,
solver stack, and observer machinery inside the two executables are byte-for-byte
the same (see CMake target list — only the `.cpp` and the registered backgrounds
differ).

### Prose: reading DAG 1

The diagram is deliberately a single shared spine with two thin branches. The
spine — free data -> elliptic system -> DG discretization -> Newton-Raphson ->
linear-solver stack -> ADM reconstruction -> observers -> H5 — is *the same C++*
in both runs. The two problems are distinguished entirely by three swapped
factory-creatable objects plus the option values fed to them:

- the **AnalyticData class** that manufactures the free data
  (`Binary` vs `TeukolskyWave`);
- the **DomainCreator** (`BinaryCompactObject` vs `Sphere`);
- the **BoundaryConditions** (`ApparentHorizon` + outer `Flatness`/`Robin`
  vs regular-center + outer `Robin`).

Everything else — the elliptic operator, Newton-Raphson, Gmres/Multigrid/Schwarz,
the ADM-quantity compute tags, the observer events — is shared verbatim. This is
exactly why a Teukolsky XCTS solve is a faithful, apples-to-apples test bed for
the same machinery that produces production BBH initial data: only the free-data
physics and the geometry/BC scaffolding change, not the solver.

The physics differences enter only at the free-data leaves:

- **(a) free-data source:** superposed boosted Kerr-Schild black holes vs a flat
  background plus the verified even-parity Teukolsky perturbation h_ij.
- **(b) domain:** `BinaryCompactObject` with two excised spheres, an envelope,
  and an outer shell to R=1e9, vs a single `Sphere` with a *regular center*
  (no excision) and a radial partitioning chosen to resolve the outgoing pulse
  near r=20.
- **(c) boundary conditions:** apparent-horizon excision BCs carrying the orbital
  shift (Omega x r + expansion a_dot r) vs no inner BC and an outer **Robin**
  asymptotic-flatness BC appropriate for a radiative configuration.
- **(d) trace K:** nonzero (the boosted Kerr-Schild slicing) vs **K = 0** (with
  dt K = 0), the maximal/transverse-traceless slice of the wave.
- **(e) wave-momentum ubar term:** zero in the quasi-equilibrium binary (the
  `LongitudinalShiftBackgroundMinusDtConformalMetric` is driven only by the
  prescribed orbital `ShiftBackground`), vs the genuinely propagating
  `-gammabar^ik gammabar^jl hdot_kl` for the Teukolsky wave (with
  `ShiftBackground = 0`).
- **(f) control loop:** the BBH pipeline wraps the XCTS solve in a control loop
  (`PostprocessId` / `ControlId.py`) that adjusts masses, spins, center-of-mass,
  and ADM linear momentum and re-solves; the Teukolsky run has **no control loop**
  — the free data is fixed analytic input and a single solve suffices.

### Summary table: identical vs different components (DAG 1)

| Component | BBH | Teukolsky | Status |
|---|---|---|---|
| Elliptic system type | `Xcts::FirstOrderSystem<HamiltonianLapseAndShift, Curved, 0>` | same | **identical** |
| Fluxes / Sources / Linearized | `Xcts::Fluxes/Sources/LinearizedSources` | same | **identical** |
| Primal fields | psi-1, (alpha psi)-1, ShiftExcess | same | **identical** |
| Executable metavariables | `SolveXcts` (with hydro obs) | `SolveXctsVacuum` (no hydro obs) | near-identical (hydro obs removed) |
| AnalyticData free-data class | `Xcts::AnalyticData::Binary<...>` (`Binary.hpp`) | `Xcts::AnalyticData::TeukolskyWave` (`TeukolskyWave.hpp`) | **different (swapped)** |
| (a) Conformal metric gammabar_ij | superposed boosted Kerr-Schild | delta_ij + Teukolsky h_ij | **different** |
| (d) Trace K, dt K | nonzero | K = 0, dt K = 0 | **different** |
| (e) ubar / Long.Shift-DtConfMetric | 0 (orbital ShiftBackground only) | -gammabar^ik gammabar^jl hdot_kl; ShiftBackground=0 | **different** |
| Matter sources | vacuum 0 | vacuum 0 | identical in value |
| DomainCreator | `BinaryCompactObject` (2 excisions, envelope, R=1e9) | `Sphere` (regular center, InnerR=0.5, R=1e9) | **different (swapped)** |
| (c) Inner BC | `ApparentHorizon<Curved>` (Theta=0, orbital shift) | none (regular center) | **different** |
| (c) Outer BC | `Flatness` (or `Robin`) | `Robin` | different |
| DG discretization | WeakInertial, Massive, GaussLobatto, Penalty 1 | same | **identical** |
| Nonlinear solver | `NewtonRaphson` (globalized) | same | **identical** |
| Linear solver stack | `Gmres`+`Multigrid`+`SchwarzSmoother`(+`MinusLaplacian`/`ExplicitInverse`) | same | **identical** |
| ADM reconstruction | `SpacetimeQuantitiesCompute` | same | **identical** |
| Observers | `field_observations`+`ObserveNorms`+`ObserveAdmIntegrals` | same (minus `LowerSpatialFourVelocity`) | near-identical |
| (f) Control loop | BBH ID control (`PostprocessId`/`ControlId.py`): mass/spin/CoM/P | none (fixed analytic free data) | **different** |
| Output | Reductions H5 + Volume H5 | same mechanism, different file names | identical mechanism |

---

## DAG 2 — XCTS generation usage graph (author -> consume)

A practical recipe graph. Each node names the concrete file/symbol to touch.

```
 (1) AUTHOR / CHOOSE A FREE-DATA CLASS
 +---------------------------------------------------------------------------+
 | A subclass of elliptic::analytic_data::AnalyticSolution (or Background +   |
 | InitialGuess) providing the XCTS free data via ~20 variables() overloads. |
 | EXISTING: Xcts::AnalyticData::Binary   (Binary.hpp)        -> BBH          |
 |           Xcts::AnalyticData::TeukolskyWave (TeukolskyWave.hpp) -> wave    |
 |           Xcts::AnalyticData::DuckTovStar (DuckTovStar.hpp) -> TOV import  |
 | NEW: copy TeukolskyWave.{hpp,cpp} as a template; implement the overloads.  |
 +-------------------------------------+-------------------------------------+
                                       | must implement these variables() tags:
                                       v
   +-----------------------------------------------------------------------+
   | REQUIRED variables() OVERLOADS (the free data + flat guesses):        |
   |  Conformal geometry:                                                  |
   |   - Tags::ConformalMetric<DataType,3,Inertial>          gammabar_ij   |
   |   - Tags::InverseConformalMetric<...>                   gammabar^ij   |
   |   - ::Tags::deriv<Tags::ConformalMetric<...>>           d_k gammabar  |
   |  Extrinsic-curvature trace sector:                                    |
   |   - gr::Tags::ExtrinsicCurvature<DataType,3>            K_ij          |
   |   - gr::Tags::TraceExtrinsicCurvature<DataType>         K             |
   |   - ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature>    d_i K         |
   |   - ::Tags::dt<gr::Tags::TraceExtrinsicCurvature>       dt K          |
   |  Conformal factor (initial guess + deriv):                            |
   |   - Tags::ConformalFactor / ConformalFactorMinusOne                   |
   |   - ::Tags::deriv<Tags::ConformalFactorMinusOne>                      |
   |  Lapse * conformal factor (initial guess + derivs):                   |
   |   - gr::Tags::Lapse / ::Tags::deriv<Lapse>                            |
   |   - Tags::LapseTimesConformalFactor(MinusOne)                         |
   |   - ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne>            |
   |  Shift sector:                                                        |
   |   - Tags::ShiftBackground<...> + ::Tags::deriv<...>                   |
   |   - Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<...>      |
   |       (the ubar / wave-momentum term)                                 |
   |   - Tags::ShiftExcess<...> + ::Tags::deriv<...>  (guess = 0)          |
   |  Matter sources (conformal_matter_scale=0); vacuum -> all 0:          |
   |   - gr::Tags::Conformal<gr::Tags::EnergyDensity, 0>                   |
   |   - gr::Tags::Conformal<gr::Tags::StressTrace, 0>                     |
   |   - gr::Tags::Conformal<gr::Tags::MomentumDensity<3>, 0>              |
   +-----------------------------------+-----------------------------------+
                                       |
                                       v
 (2) REGISTER IN AN EXECUTABLE FACTORY + CMAKE
 +---------------------------------------------------------------------------+
 | Pick the executable metavariables:                                        |
 |  - SolveXcts        (SolveXcts.hpp)        : has HydroQuantitiesCompute    |
 |                                              -> matter backgrounds         |
 |  - SolveXctsVacuum  (SolveXctsVacuum.hpp)  : NO hydro obs                  |
 |                                              -> vacuum backgrounds (USE    |
 |                                                 THIS for TeukolskyWave)    |
 | Add the class to factory_creation::analytic_solutions_and_data AND to the |
 | <Background>, <InitialGuess>, <AnalyticSolution> factory_classes pairs.   |
 | CMake: add_spectre_executable(SolveXctsVacuum ... SolveXctsVacuum.cpp)     |
 |        + target_link_libraries(... XctsAnalyticData ...) (CMakeLists.txt). |
 +-----------------------------------+---------------------------------------+
                                     |
                                     v
 (3) WRITE THE YAML INPUT FILE                (SolveXctsTeukolskyX2.yaml)
 +---------------------------------------------------------------------------+
 | Executable: SolveXctsVacuum                                               |
 | Background:  <YourClass>: {options...}    (e.g. TeukolskyWave: Amp/rc/tau)|
 | InitialGuess: *background  (reuse, or flat)                               |
 | DomainCreator: Sphere | BinaryCompactObject | ...                         |
 | Discretization: DiscontinuousGalerkin {WeakInertial, Massive, GaussLob.}  |
 | NonlinearSolver: NewtonRaphson {ConvergenceCriteria, globalization}       |
 | LinearSolver: Gmres + Multigrid + SchwarzSmoother                         |
 | (BoundaryConditions are attached to the domain: OuterBoundaryCondition:   |
 |  Robin; excisions carry ApparentHorizon)                                  |
 | Observers: VolumeFileName / ReductionFileName                             |
 | EventsAndTriggersAtIterations: ObserveNorms (H,M L2) + ObserveFields      |
 |   (ConformalFactor, Lapse, Shift, SpatialMetric, ExtrinsicCurvature, H,M) |
 +-----------------------------------+---------------------------------------+
                                     |
                                     v
 (4) BUILD                                                                   
 +---------------------------------------------------------------------------+
 | cmake --build <build> --target SolveXctsVacuum   (Charm++ build)          |
 | (run_xcts_teukolsky_x2.sh wraps the build+run)                            |
 +-----------------------------------+---------------------------------------+
                                     |
                                     v
 (5) RUN
 +---------------------------------------------------------------------------+
 | ./bin/SolveXctsVacuum --input-file SolveXctsTeukolskyX2.yaml +pN          |
 | (Charm++ "+pN" = N PEs; NewtonRaphson -> Gmres/Multigrid/Schwarz solve)   |
 +-----------------------------------+---------------------------------------+
                                     |
                +--------------------+--------------------+
                v                                         v
 (6a) REDUCTIONS H5                            (6b) VOLUME H5
 +-------------------------------+             +--------------------------------+
 | XctsTeukolskyX2Reductions.h5  |             | XctsTeukolskyX2Volume*.h5      |
 |  - NewtonRaphsonResiduals     |             |  subfile "VolumeData":         |
 |    (convergence history)      |             |   ConformalFactor, Lapse,      |
 |  - Norms: HamiltonianConstr., |             |   Shift, SpatialMetric_xx..zz, |
 |    MomentumConstraint (L2)     |             |   ExtrinsicCurvature_xx..zz,   |
 +-------------------------------+             |   Hamiltonian/MomentumConstr.  |
                                              +---------------+----------------+
                                                              |
                                                              v
 (7) CONSUME -> ATHENAK IMPORT
 +---------------------------------------------------------------------------+
 | spectre::Exporter::interpolate_to_points<3>(                              |
 |     volume_files_or_glob, subfile_name="VolumeData",                      |
 |     ObservationStep{step},                                                |
 |     {"SpatialMetric_xx",...,"ExtrinsicCurvature_zz"},                     |
 |     target_points, ...)                                                   |
 |   (Exporter.hpp interpolate_to_points)                                    |
 |   -> athenak/src/pgen/z4c_spectre_bbh.cpp :: LoadSpectreInitialData(      |
 |        pmbp, filename_glob, id_subfile_name, observation_step)            |
 |   reads SpatialMetric_xx..zz + ExtrinsicCurvature_xx..zz                  |
 |   -> pmbp->pz4c->ADMToZ4c<NGHOST>(pmbp, pin)  (ADM -> Z4c vars)           |
 |   -> Z4c evolution                                                        |
 +---------------------------------------------------------------------------+
```

### Prose: reading DAG 2

DAG 2 is the operational checklist for producing XCTS initial data and feeding it
into AthenaK. The graph is mostly linear with one fan-out at the output:

1. **Author/choose a free-data class.** The free data is supplied by a C++ class
   deriving from `elliptic::analytic_data::AnalyticSolution` (which also satisfies
   the `Background` and `InitialGuess` interfaces). It must implement the ~20
   `variables()` overloads listed above — exactly the
   `system::background_fields` plus flat initial guesses for the primal fields.
   `TeukolskyWave.hpp` is the cleanest *vacuum* template; `Binary.hpp` is the
   superposition template for compact objects.
2. **Register it.** Add the class to the chosen executable's
   `factory_creation::analytic_solutions_and_data` and to the `Background`,
   `InitialGuess`, and `AnalyticSolution` factory-class lists, and link
   `XctsAnalyticData` in CMake. Use `SolveXctsVacuum` for vacuum free data so the
   `HydroQuantitiesCompute` observation tag never asks the background for matter
   primitives.
3. **Write the YAML.** Choose `Background`/`InitialGuess` (your class), the
   `DomainCreator`, the solver stack, and the observed fields. Boundary
   conditions attach to the domain (outer `Robin` for radiative/asymptotically
   flat problems; `ApparentHorizon` on excisions for BHs).
4. **Build** the executable (Charm++).
5. **Run** with `+pN` PEs; Newton-Raphson drives the Gmres/Multigrid/Schwarz
   linear solves.
6. **Outputs:** a Reductions H5 (convergence history + constraint L2 norms) and a
   Volume H5 (the ADM fields psi, alpha, beta^i, gamma_ij, K_ij + constraints).
7. **Consume:** `spectre::Exporter::interpolate_to_points<3>` reads the Volume H5
   `VolumeData` subfile and interpolates `SpatialMetric_xx..zz` and
   `ExtrinsicCurvature_xx..zz` to AthenaK grid points;
   `z4c_spectre_bbh.cpp::LoadSpectreInitialData` fills the ADM fields, then
   `ADMToZ4c<NGHOST>` converts to Z4c variables for evolution.

---

## File / symbol index (for later reuse)

| Role | Path / symbol |
|---|---|
| Teukolsky free data | `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.{hpp,cpp}` (copy in `progress/x2-xcts/spectre_src/`) |
| BBH free data | `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/Binary.{hpp,cpp}` |
| TOV import free data | `/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/DuckTovStar.{hpp,cpp}` |
| XCTS elliptic system | `/data/haiyangw/nr/spectre/src/Elliptic/Systems/Xcts/FirstOrderSystem.hpp`, `Equations.hpp`, `Tags.hpp` |
| XCTS boundary conditions | `/data/haiyangw/nr/spectre/src/Elliptic/Systems/Xcts/BoundaryConditions/{Factory,ApparentHorizon,Robin,Flatness}.hpp` |
| Vacuum executable | `/data/haiyangw/claude/z4c-CMM/progress/x2-xcts/spectre_src/SolveXctsVacuum.{hpp,cpp}` (and `CMakeLists.txt`) |
| Standard executable | `/data/haiyangw/claude/z4c-CMM/progress/x2-xcts/spectre_src/SolveXcts.hpp` |
| Teukolsky YAML | `/data/haiyangw/nr/spectre/SolveXctsTeukolskyX2.yaml` |
| BBH YAML (reference) | `/data/haiyangw/nr/spectre/tests/InputFiles/Xcts/BinaryBlackHole.yaml` |
| BBH pipeline + control | `/data/haiyangw/nr/spectre/support/Pipelines/Bbh/{InitialData.yaml,PostprocessId.py,ControlId.py}` |
| Exporter | `/data/haiyangw/nr/spectre/src/IO/Exporter/Exporter.hpp` (`interpolate_to_points<3>`) |
| AthenaK importer | `/data/haiyangw/claude/z4c-CMM/athenak/src/pgen/z4c_spectre_bbh.cpp` (`LoadSpectreInitialData`, `ADMToZ4c`) |

---

## Companion figure

`xcts_dag.tex` (this directory) renders both DAGs as a single standalone,
paper-quality TikZ figure matching the style of the Z4c-CCM derivation-chain
figure in `paper/z4c-CMM/zccm_formulation/sections/matching.tex` (black `box`
nodes for shared/foundational components, dark-blue `newbox` nodes for the new /
swapped Teukolsky-specific content, `Latex`-tipped gray arrows). It compiles with
`pdflatex xcts_dag.tex`.
