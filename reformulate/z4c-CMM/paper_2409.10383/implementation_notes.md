# Implementation Notes — arXiv:2409.10383 (AthenaK Z4c module)

All items are literature-grounded transcriptions from `ref-paper/arxiv-2409.10383/src/main.tex`; section labels verbatim. [OPEN] marks facts this paper does not settle.

## 1. Boundary conditions

- Ghost-zone mechanics (Sec. `sec:num`): each MeshBlock = active region + ghost region; "the ghost zones are filled by either the active region of neighboring MeshBlocks or the boundary conditions" at **each stage** of the explicit RK time integration. => A new physical BC (e.g. CCM injection) plugs in at the per-stage ghost-fill point.
- Ghost-cell counts: N_g = 2 (default, 2nd-order FD) or N_g = 4 (6th-order FD). A CCM boundary module must fill N_g ghost layers per RK stage.
- Physical outer BC type: **unspecified in this paper** (no Sommerfeld discussion; see convention.md Sec. 2). Tests place the outer boundary at +/-1024 M. Periodic BCs used only for the linear-wave test.
- Refinement (2:1) interior boundaries: prolongation/restriction via Lagrange polynomial interpolation, stencil size n = N_g + 1 per dimension => error converges at 3rd order (2 ghosts) / 5th order (4 ghosts); footnote states P/R error never dominates spatial-differencing error in practice.
- CC-scheme-specific communication change: coarse buffers must also be exchanged between *same-level* neighbors at a refinement boundary (MeshBlocks B and C in Fig. `fig:p_and_r`) because the prolongation stencil overlaps those nodes. Restriction stencils are kept entirely within the active region ("to improve locality"); asymmetric stencil, "yet to find any sacrifice in accuracy".

## 2. Wave extraction implementation

- Verbatim (Sec. 3.4): Psi4 from a **coordinate tetrad**, interpolated onto **geodesic spheres** at various radii, decomposed in **spin-weighted spherical harmonics**. No equations; no E/B-Weyl statement in this paper. [OPEN: confirm E_ij/B_ij route in Daszuta:2021ecf / AthenaK source.]
- Extraction radii in tests: 50 M (single puncture), 100 M (BBH calibration).
- Puncture trackers: initialized at puncture locations of initial data, advected by integrating the shift vector each time step (relevant if CCM worldtube must track sources).

## 3. CCE / worldtube hooks

- Verbatim (Sec. 3.4): "Ongoing work includes a horizon finder and **outputting worldtube data for Cauchy Characteristic Extraction** [Bishop:1996gt, Reisswig:2009rx, Bishop:2016lgv, Barkett:2019uae, Moxon:2020gha, Moxon:2021gbv]."
- I.e., at paper time (Sep 2024) AthenaK had **no shipped worldtube output**; CCE (one-way) was the stated target, not CCM (two-way). Cited CCE stack: Bishop et al. 1996 (PRD 54, 6153); Reisswig et al. 2010 (CQG 27, 075014); Bishop & Rezzolla 2016 (Living Rev. Rel. 19, 2); Barkett et al. 2020 (PRD 102, 024004); Moxon, Scheel & Teukolsky 2020 (PRD 102, 044052); Moxon et al. 2023 (PRD 107, 064013) — the SpECTRE CCE line.
- Companion sources already in repo: arXiv:2007.01339 (improved CCE system), arXiv:2308.10361 (fully relativistic 3D CCM, SpECTRE) at `ref-paper/`.

## 4. Discretization / time integration

- Spatial: cell-centered finite differencing (CC; GR-Athena++ used vertex-centered VC). Default pairs: 2nd-order FD with 2 ghosts; 6th-order FD with 4 ghosts. Motivation for CC: Z4c variables co-located with finite-volume fluid variables, removing metric<->matter interpolation cost.
- Time: "a family of explicit Runge-Kutta methods"; "the time integrator is kept at 4th order" (RK4) in all tests; CFL = 0.25 (BBH calibration).
- Kreiss-Oliger dissipation: sigma = 0.5 (BBH calibration; GR-Athena++ comparison runs used 0.02). Macro `\KO = sigma`.
- Convergence evidence: linear wave (amplitude 1e-8, periodic, diagonal propagation) converges at expected orders; N_g=4 + RK4 shows 6th-order (FD error dominates); plateau at high resolution traced to quadratic (nonlinear) effects with coefficient O(10). BBH waveform: clean 4th-order convergence (Q-factor, label `eq:qfactor`); footnote: order limited by sharp early-time lapse features, cf. Etienne:2024ncu.

## 5. Mesh / AMR (relevant to placing a worldtube/boundary)

- Oct-tree block AMR, refinement by factors of 2, 2:1 neighbor condition; MeshBlocks grouped into MeshBlockPacks (fewer kernel launches; key GPU optimization).
- AMR criteria: (i) chi-min — refine MeshBlock if min(chi) < threshold (0.3 in calibration run); (ii) dchi-max — refine if |grad chi| exceeds threshold (handles mass ratios; q=20 demo gets 4 extra levels on the small hole); (iii) **radius criterion** — "all meshblocks within certain radius of a given point are kept above a specified refinement level" — designed to keep wave-zone (r >~ 100 r_g) resolution; this is the natural lever for guaranteeing resolution at a CCM worldtube radius.
- Cost model stated: +1 max refinement level => ~2x cost (CFL-driven, ~10% more gridpoints) vs global 2x refinement => 16x; wave-zone resolution then unchanged at fine levels — inspiral accuracy preserved, merger-ringdown (higher frequency content) degraded.

## 6. Performance portability (Kokkos)

- Kokkos abstracts machine-specific backends; demonstrated on NVIDIA (A100, Grace Hopper, RTX 3070 Ti), AMD (MI250X), Intel/AMD/Apple CPUs. ZCPS table (`tab:portability`): A100 1.36e7 ZCPS; >= 200x vs single CPU core.
- Weak scaling (Frontier, single BH, 10 SMR levels): 80% efficiency from 4 to 65,536 GPUs with 64^3 MeshBlocks. Strong scaling: 84% (Frontier) / 77% (Perlmutter) at 32x resources; 64^3 -> 128^3 MeshBlocks gives ~2x performance on MI250X.
- Implication for CCM: any boundary/worldtube code must be expressible as Kokkos kernels over MeshBlockPacks and must not serialize per-stage ghost-fill communication (cf. non-locality remark used to argue against excision/horizon finding in Sec. `sec:z4c`).

## 7. Units
- M = ADM mass of initial data; geometric units implied (not stated). Resolutions and radii quoted in M.

## CCM-relevant facts (digest)

1. AthenaK Z4c continuum system == GR-Athena++ system: Eqns. 8-13 + gauge Eqn. 22 of Daszuta:2021ecf (ApJS 257, 25). This paper reprints none of it. [OPEN]
2. No Sommerfeld/outer-BC equations appear in this paper; physical outer BC type is unspecified here; tests use far boundaries at +/-1024 M (and periodic for linear waves). [OPEN -> Daszuta:2021ecf; constraint-preserving BCs in arXiv:1010.0523.]
3. Ghost zones (2 or 4 layers) are refilled at *every RK stage* from neighbors or BCs — the hook point for CCM injection.
4. FD is cell-centered; CCM boundary data must be provided at cell centers, N_g layers deep, compatible with 6th-order stencils.
5. Time integration: explicit RK family, RK4 in practice, CFL 0.25; BC application is per-stage, not per-step.
6. Production dissipation/damping: KO sigma = 0.5, kappa1 = 0.02, kappa2 = 0 (Z4c damping active in kappa1 only).
7. Psi4: coordinate tetrad -> geodesic-sphere interpolation -> SWSH decomposition; extraction at 50-100 M; E/B-Weyl route not documented here. [OPEN]
8. CCE worldtube output was "ongoing work" (not delivered); cited target stack is SpECTRE CCE (Moxon et al.). No CCM (two-way) infrastructure exists in the paper.
9. The radius AMR criterion can pin a minimum refinement level inside any chosen radius — usable to guarantee uniform resolution at a CCM worldtube/boundary sphere.
10. All new code must be Kokkos-kernel over MeshBlockPacks to preserve the demonstrated 80% weak scaling at 65,536 GPUs.
