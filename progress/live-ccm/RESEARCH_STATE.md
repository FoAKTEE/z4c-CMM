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

## Accepted results log

(empty — populated on verifier admission only)
