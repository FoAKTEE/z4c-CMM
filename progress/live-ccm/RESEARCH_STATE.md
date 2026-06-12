# z4c-CCM — Live-CCM-in-AthenaK Mission (mission 3) — Verified Research State

User directive (2026-06-12, verbatim core): "first put constraint preserving
to Athena. In the mean time you should make Cauchy and characteristic can run
simultaneously in athena, instead of postprocessing. What is needed, you
should use CCM 2023 paper, CCM 2023 DAG and CCE paper (ref-paper,
reformulation). But even on step before this, reference should be exactly the
same as the analytical results. use phys-agentic-loop for management."

Diagnosis the directive corrects: the mission-2 AthenaK implementation drove
only the two physical rows (Sommerfeld everywhere + additive self-datum
injection); its reference run carried a 39% truncation deficit against the
analytic Teukolsky waveform (2nd-order stencils at ~4 points per tau); and
the CCM datum was analytic/postprocessed, never a live characteristic
evolution. All three are corrected by the stages below, in order.

## Stage decomposition (ledger node IDs continue the z4c-CCM DAG)

| Node | Stage | Claim | Status |
|---|---|---|---|
| N12 | reference = analytic | Exact linear psi4(t,r,theta) of the Teukolsky solution at FINITE radius derived and verified (sympy; limits: scri formula, psi0 datum, linearized vacuum); AthenaK runs at 6th-order stencils (nghost=4) converge to it with measured order, Richardson limit consistent; checker gate replaces the order-of-magnitude W1 | [HYPOTHESIS] |
| N13 | CPBC in AthenaK | The constraint-preserving (Theta, Z_s, Z_A) and paper-1 gauge rows of the Z4c-CCM boundary operator (formulation eq:bc-cpbc/eq:bc-gauge) implemented in AthenaK as Bjørhus rows replacing plain Sommerfeld; verifier: boundary-sourced constraint growth suppressed vs Sommerfeld at equal resolution, Minkowski preserved, convergence unchanged | [HYPOTHESIS] |
| N14 | simultaneous CCM | A characteristic Bondi-Sachs module (2007.01339 hierarchy; worldtube data map; compactified radial grid; spin-weighted angular basis, axisymmetric m=0 first) evolving INSIDE AthenaK in lockstep with the Cauchy evolution; worldtube boundary data from the live Cauchy state each step; psi0 at the worldtube fed back through the Type-II boost into the N13 physical rows (datum mode 4 = live). No postprocessing anywhere in the loop. Verifiers: standalone linear Teukolsky characteristic solve vs analytic strain at scri; cross-oracle vs SpECTRE CharacteristicExtract on identical worldtube data; full live-CCM Teukolsky battery | [HYPOTHESIS] |

Stage order is binding: N12 gates N13 (convergence claims need the exact
reference); N13 gates N14 (live datum lands on the full boundary operator).

## Source library (all decomposed, mission 1)

| ID | Use here |
|---|---|
| 1010.0523v2 (Z4c CPBC) | N13: explicit CPBC + gauge rows; knowledge-database/paper_arxiv-1010.0523v2 |
| 2007.01339 (CCE) | N14: full hypersurface hierarchy beta/Q/U/W + dJ/du, worldtube transformation; knowledge-database/paper_arxiv-2007.01339 |
| 2308.10361 (GH-CCM) | N13/N14: Bjørhus physical rows, Type-II boost, matching protocol; knowledge-database/paper_arxiv-2308.10361 |
| synthesis | reformulate/z4c-CMM/synthesis_z4c-CCM/ (ONE DAG + explicit 10-mode BC table) |

## Binding user constraints (inherited + this mission)

1. No test run > 10 minutes; <= 4 GPUs; 64-core builds.
2. Commits per substage, detailed docs, NO Claude attribution; push to origin.
3. Reference runs must agree with analytic results to quantified truncation
   (measured convergence order + Richardson consistency), not order of
   magnitude.
4. Cauchy and characteristic evolutions run simultaneously in AthenaK —
   in-process lockstep coupling; file/postprocessing round-trips are
   disqualified.

## Active claims

N12 in progress (this iteration). N13, N14 queued.

## N12 findings (iteration 18)

1. **Teukolsky pgen frame bug (root cause).** The pgen's
   `eth = rh x eph = -theta-hat` flipped the h_{r theta} cross term: the
   evolved "Teukolsky wave" violated the linearized vacuum equations at 5%
   of |Ricci| (FD residual 1.8e-4 vs 3.9e-3 at 30-digit precision; 1.9e-34
   after the flip). Invisible to constraint monitors (O(X^2)+truncation
   floor) and to run-vs-reference protocols (consistent on both sides) —
   caught only by the exact-analytic-reference program. Error-db row
   2026-06-12T16:16. Fix: `eth = eph x rh` (pgen) + `eph x s` (injection
   dyad, Im psi0 channel only).
2. **Exact finite-radius waveform (N12a, OVERALL PASS).** In the AthenaK
   extraction convention: rpsi4(2,0) = sqrt(6pi/5)[F6 + 2F5/r + 3F4/r^2
   + 3F3/r^3 + (3/2)F2/r^4](t-r) - sqrt(6pi/5)(3/2)F2(t+r)/r^4. Gates:
   traceless, linearized vacuum coefficient-exact, sin^2-factorization,
   peeling (k+j=6), scri limit, independent FD cross-check (1.4e-40).
   Convention map: rpsi4_AthenaK = (-1) x the paper formula
   -sqrt(6pi/5)F6.
3. **psi0 datum channel exact + orientation map.** In the concrete
   right-handed (s, e_th, e_ph) frame with eps = +[ikl]:
   Re[(E - eps(s)B)(m,m)] = -sqrt(27pi/10)F2(t-r)/r^5 * 2Y20 EXACTLY (no
   other retarded term; FD-adjudicated r^-5 scaling at r=36/72). The N3
   abstract identity (E + eps B)(m,m) = -psi0 corresponds under
   eps -> -eps (m <-> mbar); the AthenaK injection w_ab is quadratic in
   e_th and unaffected for real datums.
4. **Waveform time-label bug.** The Weyl/waveform tasks are Task_End tasks
   (state at t + dt) but the writer labeled rows with the pre-advance
   pmesh->time: data exactly one step early vs label (measured best-fit
   shifts = dt at all three resolutions; raw E 12-18% with peak ratio 0.98).
   Fixed: label = time + dt. After shift removal the ladder agreed with the
   exact waveform to 3.4%/0.38%/0.31% (h = 0.5/0.331/0.277, 6th order).

## Accepted results log

| Claim | Evidence type | Evidence / verifier | Status |
|---|---|---|---|
| N12a: exact linear finite-radius rpsi4(2,0) and psi0 of the Teukolsky solution in the AthenaK convention, 7 gates incl. linearized vacuum (coefficient-exact) and an independent mpmath-FD cross-check (1.4e-40); convention map AthenaK = -1 x paper; datum channel = (E - eps(s)B)(m,m) in the concrete frame | exact symbolic + independent FD | `scripts/verify_n12_exact_psi4.py` -> results/numerical/n12_exact_psi4_check.txt, coefficients n12_psi4_exact_coeffs.json | [SOLID] |
| N12b: AthenaK 6th-order ladder converges to the exact waveform — E = 4.99e-2 / 5.36e-3 / 2.02e-3 at h = 0.5/0.331/0.277, measured order 5.39/5.52, finest amplitude 1±2e-4, shift 0.000 | dynamical convergence vs exact truth | `scripts/test_athenak_teuk_ana.sh` + `check_athenak_teuk_ana.py` -> results/numerical/athenak_teuk_ana/gpu6_checker.txt | [SOLID] |
| Two production bugs root-caused: pgen frame (-theta-hat, h_{r theta} sign) violating linearized vacuum at 5% of |Ricci|; waveform time label one step behind the Task_End state | error ledger | error-database rows 2026-06-12T16:16; N11 amended | [SOLID] |

N12 stage table status: [SOLID].

## N13 status (iteration 19): [PRELIMINARY]

v1 shipped: `z4c/cpbc` runtime flag; Gamma-tilde row = stock dissipative
advection + causal-rate constraint sink `-2 Z_a/r`,
`Z_a = (Gamma_a - d_j gt_aj)/2`. Holds: stability, Minkowski preservation
(2.2e-14), waveform fidelity, 6th-order convergence (p = 5.39), exact
linear-amplitude consistency (cpbc == somm when Z ~ roundoff). NOT yet
held: measurable absorption (late-window int L2-H ratio 0.9999 at X = 0.05
vs gate 0.95). Failed attempts ledgered: L=1 advected row (face-stencil
ghost reads, NaN t=19.75); L=0 metric-transparent row (non-dissipative
face loop, NaN t=43.25).

- **O-N13-1 (gating the solid claim):** characteristic-exact Z transport
  with INWARD one-sided first-derivative face stencils (no outflow-ghost
  reads), GKS/LF-checked (mission-2 model1d closure lesson).
- **Config note:** X = 2 with linear ID is under-dissipated at nghost=4
  (NaN ~ t=50 even for stock Sommerfeld; raise diss/kappa1 or solve ID —
  N14-era X=2 batteries).

## N14 next-iteration plan

Physical-row datum path is independent of the Z row — N14 proceeds.
Design: host-side axisymmetric (m = 0) Bondi-Sachs module in
athenak/src/z4c/ccm/ — hierarchy beta/Q/U/W + dJ/du from the 2007.01339
equation ledger; compactified radius y = 1 - 2 r_wt/r on a Chebyshev or
high-order FD grid; Gauss-Legendre theta collocation; worldtube Bondi data
from the live Cauchy state via the cce/ sphere-interpolation machinery
(in-memory); psi0 at the worldtube -> Type-II boost -> Z4cCCMDatumPsi0
mode 4 (live), lockstep du = Cauchy dt. Verifiers: standalone linear
Teukolsky characteristic solve vs the N12 exact module; SpECTRE
CharacteristicExtract cross-oracle on identical worldtube data; live-CCM
battery vs the N12 analytic truth.

## Binding directive added (iter 29, 2026-06-12)

FULL MIGRATION of the Python toolchain to Julia (packages/ZccmJl): GPU via
CUDA.jl, symbolic/exact via Rational{BigInt} series + Symbolics.jl,
arbitrary precision via BigFloat. Further development Julia-only once each
ported component passes honest unit tests reproducing the admitted
Python-verified numbers. Python artifacts stay as provenance baseline;
AthenaK (C++) unaffected. N15 [PRELIMINARY]: v0.1.0 baseline admitted
(teukolsky/linseries/psi0chain, all tests pass, refs auto-generated from
the N12 module). Sympy hit its kernel_timeout wall on the Bondi-transform
derivation (error-db iter 28) — that derivation is the next Julia unit
(LinSeries x angular-shape exact algebra).
