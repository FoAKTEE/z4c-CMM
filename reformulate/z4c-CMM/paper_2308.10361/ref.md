# ref.md — references used directly/semi-directly in the calculations of arXiv:2308.10361

Citation keys are from `ref-paper/arxiv-2308.10361/src/References.bib` / `paper.bbl`.

## Load-bearing (equations imported verbatim or relied upon)

- `Lindblom:2005qh` — Lindblom, Scheel, Kidder, Owen, Rinne (2006): first-order GH
  evolution system. Source of FOSH variables, characteristic fields u^0̂/u^1̂±/u^2̂
  (their Eqs. 32-34), constraint-preserving BCs (Eqs. 63-65), constraint fields
  c³ (Eq. 57), gauge constraint C_a (Eq. 40), three-index constraint (Eq. 26).
- `Kidder:2004rw` — Kidder et al. (2005): Bjørhus-style boundary-condition
  implementation; classification into constraint-preserving/physical/gauge subsets;
  ψ₀-freezing physical BC that CCM replaces.
- `doi:10.1137/0916035` — Bjørhus (1995): the boundary method of replacing time
  derivatives of incoming characteristic fields (Bjorhus_bc).
- `Moxon:2020gha` — Moxon, Scheel, Teukolsky (2020): SpECTRE CCE formulation.
  Source of the Bondi-Sachs form (eq:BS_PFB), coordinate-system chain (their Table I),
  the tetrad eq:CCE_tetrad, the ψ₀ expression eq:psi0_CCE, β̂ = −½ln(∂_λ̲ r̂)
  (their Eqs. 19a, 33a), the Ĵ/K̂ → J/K transformation (eq:J_BS_like), and the
  spin-weighted scalar definitions (with a typo in their Eq. 10e corrected here).
- `Moxon:2021gbv` — Moxon et al. (2021): SpECTRE CCE implementation; worldtube data,
  numerical settings, spin-weighted Clenshaw interpolation used in Sec. IV.D.
- `PhysRevD.26.745` — Teukolsky (1982): linearized quadrupolar Teukolsky-wave metric
  (eq:Teukolsky_metric, eq:Teukolsky_ABC, eq:Teukolsky_h_20, appendix analytic
  Weyl/News/strain expressions).
- `Gomez:1996ge` — Gómez et al.: eth formalism / spin-weighted derivatives on the sphere.

## Semi-direct (methods, diagnostics, well-posedness context)

- `O-M-Moreschi_1986`, `Iozzo:2020jcu` — Bondi-gauge Bianchi-identity constraints
  (eq:bondi_violation) used as waveform-quality diagnostics.
- `Flanagan:2015pxa` — reality condition on the Bondi mass aspect (eq:im_psi2).
- `Rinne:2007ui` — Sommerfeld gauge boundary conditions adopted for the gauge subset.
- `Giannakopoulos:2020dih`, `Giannakopoulos:2021pnh`, `Giannakopoulos:2023zzm` —
  weak hyperbolicity of the Bondi-like characteristic systems; CCM not well-posed.
- `York:1998hy`, `Pfeiffer:2002iy` — XCTS initial-data formulation used for test setups.
- `Chen:2021rtb` — spherical Kerr-Schild coordinates for the Kerr conformal metric.
- `1974ApJ...193..443T` — Teukolsky (1974): r^-5 falloff of the ingoing ψ₀ (Table I),
  explaining the weak backscatter at large boundary radius.
- `Barkett:2019uae` — SpEC CCE; Cauchy-to-characteristic worldtube data conventions.
- `Handmer:2014qha`, `Bishop:1997ik` — sources of the radial hypersurface-equation
  hierarchy source functions (eq:cce_com, commented-out appendix).
- `scri`, `mike_boyle_2020_4041972`, `Boyle:2013nka`, `Boyle:2014ioa`, `Boyle:2015nqa`
  — `scri` package used to compute Bondi-constraint norms.
- `Stein:2019mop` — `qnm` package used to verify quasinormal-mode frequencies.
- `SpECTREDomain`, `spectrecode`, `Kidder:2016hev` — SpECTRE code and domain
  infrastructure of the simulations.
- `SpECwebsite`, `SXSWebsite`, `SXSCatalog` — SpEC reference simulation (X=2 test).
- `doi:10.1063/1.525904`, `winicour2012characteristic` — CCE/CCM historical context
  (background only).
- `Frittelli:1999yr`, `Gomez:2003ew`, `Frittelli:2004pk` — existence/uniqueness and
  hyperbolicity context for characteristic systems.
- `CCMData` — GitHub repository with the simulation data of the paper's tests.

Status note: all of the above are [ASSUMPTION]-level imports (literature grounding);
none has been independently re-derived in this stage.
