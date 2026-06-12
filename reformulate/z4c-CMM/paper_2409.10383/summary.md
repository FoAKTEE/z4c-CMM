# Summary — arXiv:2409.10383

**Title:** Performance-Portable Numerical Relativity with AthenaK
**Authors:** Zhu, Fields, Zappa, Radice, Stone, Rashti, Cook, Bernuzzi, Daszuta
**Source:** `ref-paper/arxiv-2409.10383/src/main.tex` (681 lines, single-file LaTeX)
**Evidence type for everything below:** literature grounding (transcription from source).

## Scope

Methods/validation paper for the **vacuum numerical-relativity (Z4c) module of AthenaK**, a Kokkos-based, performance-portable extension of Athena++ targeting exascale GPU systems. Part of a three-paper AthenaK series (with Stone et al. and Fields et al., cited as `Stone.6.24`, `Fields.6.24`). Code: https://github.com/IAS-Astrophysics/athenak.

Paper structure (verbatim section labels):
- `sec:z4c` — Z4c formulation choice (summary only; equations deferred, see below).
- `sec:num` — numerical implementation: Mesh/MeshBlock/MeshBlockPack, modified refinement-boundary communication (prolongation/restriction via Lagrange interpolation), cell-centered (CC) finite differencing vs GR-Athena++'s vertex-centered (VC), new chi-based AMR criteria.
- `sec:test` — linear wave convergence, single spinning puncture vs BAM, BBH calibration run (convergence + cross-code comparison with GR-Athena++).
- `sec:scale` — performance portability table and Frontier/Perlmutter scaling.
- `sec:conclude`.

## What the AthenaK Z4c module provides (as stated in this paper)

1. **Formulation:** conformally-decomposed Z4 ("Z4c") for vacuum Einstein equations; puncture-compatible (no excision); admits natural constraint damping. Continuum equations and gauge are stated to be **"exactly the same as in GR-Athena++"** (Daszuta et al. 2021, Astrophys. J. Supp. 257, 25 = arXiv:2101.08289, cite key `Daszuta:2021ecf`); the paper points to **"Eqn. 8-13"** (evolution system) and **"Eqn. 22"** (gauge) of that reference rather than reprinting them.
2. **Gauge:** Bona-Masso lapse (`Bona:1994a`) + gamma-driver shift (`Alcubierre:2003hr`), again by reference to `Daszuta:2021ecf` Eqn. 22.
3. **Discretization:** cell-centered high-order finite differencing (default 2nd order with 2 ghosts, 6th order with 4 ghosts), explicit Runge-Kutta family (RK4 in all production tests), oct-tree block AMR with 2:1 refinement and MeshBlockPacks; Kokkos for portability.
4. **AMR criteria:** new puncture-agnostic chi-min and dchi-max criteria, plus a *radius* criterion that keeps all MeshBlocks within a given radius of a point above a specified refinement level (used to protect wave-zone resolution).
5. **Diagnostics:** puncture trackers (advected by integrating the shift vector) and wave extraction (Psi4 via a coordinate tetrad, interpolated to geodesic spheres at various radii, decomposed in spin-weighted spherical harmonics).
6. **CCE status:** worldtube data output for Cauchy-Characteristic Extraction is explicitly **"ongoing work"** (together with a horizon finder) — not a delivered feature in this paper.

## What the paper does NOT contain (important for the CCM task)

- **No explicit Z4c evolution equations** — only the cross-reference to `Daszuta:2021ecf` Eqns. 8-13 and 22. [OPEN obligation: transcribe from arXiv:2101.08289.]
- **No enumeration of the evolved variable set** (chi, gtilde_ij, Khat, Atilde_ij, Theta, Gamtilde^i). Only the conformal factor chi is discussed by name (AMR criteria: "The conformal factor chi goes to zero towards the punctures"). The preamble defines a macro `\defG = \widehat{\Gamma}` that is never used in the body.
- **No Sommerfeld / outer-boundary-condition equations or discussion.** The physical outer boundary treatment is never specified; only periodic BCs (linear wave test) and ghost-zone communication mechanics appear. [OPEN: Sommerfeld BCs are in `Daszuta:2021ecf`; constraint-preserving alternatives in Ruiz-Hilditch-Bernuzzi arXiv:1010.0523, already at `ref-paper/arxiv-1010.0523v2`.]
- **No Psi4 equations** and no statement that the computation routes through electric/magnetic Weyl parts — only the one-sentence coordinate-tetrad description. [OPEN: the E/B-Weyl route must be verified in `Daszuta:2021ecf` and/or the AthenaK source (z4c Weyl-scalar kernel), not in this paper.]
- Only **one labeled equation** in the entire paper: `eq:qfactor` (convergence Q-factor). All other displayed equations are unlabeled.

## Key validated numbers (for reproduction/calibration)

- BBH calibration run (BAM calibration problem `Bruegmann:2006ulg`; TwoPunctures initial data `Ansorg:2004ds`): KO dissipation sigma = 0.5, kappa1 = 0.02, kappa2 = 0, CFL = 0.25, 6th-order FD + RK4; Mesh extends -1024 to +1024 M (M = ADM mass of initial data), root grid 4^3 MeshBlocks, 10 AMR levels with chi-min threshold 0.3; puncture resolutions 0.03125 / 0.015625 / 0.0078125 M; Psi4 extracted at coordinate radius 100 M; observed 4th-order convergence of amplitude and phase.
- Single spinning puncture (a = 0.5, TwoPunctures ID with secondary mass 1e-12 M at separation 1e-5 M): SMR grid to +/-1024 M, dx = 0.08333 M at puncture, 0.66667 M in wave zone; Psi4 (l=2, m=0) at r = 50 M consistent with BAM (which used 4th-order FD).
- Performance: >= 200x GPU vs single CPU core; 80% weak-scaling efficiency at 65,536 MI250X GPUs on Frontier (relative to 4); 84% / 77% strong scaling on Frontier / Perlmutter at 32x resources.
