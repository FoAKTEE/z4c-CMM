# Conventions and Equations — arXiv:2409.10383 (AthenaK NR module)

Binding note: this paper contains almost **no formulation-level equations**. Its Section 2 (`sec:z4c`) delegates the entire continuum system to the GR-Athena++ method paper. Everything below is transcribed from `ref-paper/arxiv-2409.10383/src/main.tex` with verbatim labels where they exist; absences are recorded explicitly as [OPEN] obligations rather than filled in from memory.

## 1. Formulation and variable conventions

Verbatim (Sec. 2, `sec:z4c`):

> "As in GR-Athena++, we evolve the Einstein equation in the conformally-decomposed Z4 formulation [Bona:2003fj], or Z4c [Bernuzzi:2009ex, Hilditch:2012fp]. The continuum form of our equations and gauge choices are exactly the same as in GR-Athena++ [Daszuta:2021ecf], to which we refer interested readers for details."

> "We refer readers to Eqn. 8-13 of [Daszuta:2021ecf] for the exact form of equations used in our implementation. In addition to these equations, one must also specify the gauge conditions to close the system. Same as in GR-Athena++, we use the Bona-Masso lapse [Bona:1994a] and the gamma-driver shift [Alcubierre:2003hr] (Eqn. 22 in [Daszuta:2021ecf])."

Resolved citations (from `main.bbl`):
- `Daszuta:2021ecf` = Daszuta, Zappa, Cook, Radice, Bernuzzi, Morozova 2021, Astrophys. J. Supp. 257, 25, doi:10.3847/1538-4365/ac157b (= GR-Athena++, arXiv:2101.08289).
- `Bernuzzi:2009ex` = Bernuzzi & Hilditch 2010, Phys. Rev. D81, 084003.
- `Hilditch:2012fp` = Hilditch, Bernuzzi, Thierfelder, Cao, Tichy, Bruegmann 2013, Phys. Rev. D88, 084057.
- `Bona:1994a` = Bona, Masso, Stela, Seidel 1996 (MG7 proceedings).
- `Alcubierre:2003hr` = Alcubierre et al. 2003, Class. Quant. Grav. 20, 3951.

### Evolved variables
The paper never lists the Z4c state vector. Variables explicitly named in the text:
- **chi** — "conformal factor"; property used: "goes to zero towards the punctures" (Sec. 3.3, AMR criteria). Note chi -> 0 at punctures fixes the GR-Athena++ chi-convention (chi = det-normalized conformal factor, not the BSSN phi).
- **kappa1, kappa2** — Z4c constraint-damping parameters; values used in BBH calibration: "constraint damping parameters kappa_1 = 0.02, kappa_2 = 0" (Sec. 4.3).
- Preamble macro `\defG` = widehat(Gamma) (the Z4c conformal connection function widehat-Gamma^i) is defined but **unused in the body**.

[OPEN] Full state vector {chi, gtilde_ij, Khat, Atilde_ij, Theta, Gamtilde^i (= widehat-Gamma^i)}, the evolution equations, and the kappa1/kappa2 damping terms must be transcribed from `Daszuta:2021ecf` Eqns. 8-13 (and the constraint definitions there). Not admissible from this paper alone.

### Gauge equations
[OPEN] Bona-Masso lapse + gamma-driver shift = Eqn. 22 of `Daszuta:2021ecf`. Not reprinted in this paper.

## 2. Sommerfeld / outer boundary conditions

**Absent.** The string "Sommerfeld" does not occur in `main.tex`. The only boundary-condition statements are:

- Sec. 3 (`sec:num`): "Each MeshBlock consists of an active region and a ghost region. Data in the ghost region must be communicated between MeshBlocks and is used to set boundary conditions. At each stage of the time integration process, where we employ a family of explicit Runge-Kutta methods, the ghost zones are filled by either the active region of neighboring MeshBlocks or the boundary conditions."
- Sec. 4.1 (linear wave test): "We set periodic boundary condition..."
- Outer boundary placement in tests: +/-1024 M (single puncture, Sec. 4.2; BBH calibration, Sec. 4.3.1) — far boundary, causally disconnected strategy.

[OPEN] The Sommerfeld radiative BC equations used by the GR-Athena++/AthenaK lineage are Eqns. (23)-(25)-type conditions in `Daszuta:2021ecf` (exact numbering to be verified there). Constraint-preserving boundary alternatives for Z4c: Ruiz, Hilditch & Bernuzzi, arXiv:1010.0523 (`ref-paper/arxiv-1010.0523v2`).

## 3. Wave extraction / Weyl scalar Psi4

Complete verbatim statement (Sec. 3.4, "Other Diagnostics"):

> "To extract gravitational waves, we compute the outgoing Weyl scalar Psi_4 using a coordinate tetrad, then interpolate Psi_4 onto geodesic spheres at various radii and decompose it into spin-weighted spherical harmonics."

- **No equations given.** No tetrad definition, no statement on whether Psi4 is assembled from the electric/magnetic parts of the Weyl tensor (E_ij, B_ij). [OPEN: route through E/B Weyl parts is the GR-Athena++ approach (`Daszuta:2021ecf` Sec. on wave extraction); confirm against AthenaK source.]
- Extraction radii used: r = 50 M (single puncture, l=2 m=0 mode), r = 100 M (BBH, l=m=2 mode).
- Wave-zone caveat (Sec. 3.3.3): "the gravitational waveform is usually extracted at large radii (>~ 100 r_g), and a low resolution at the wave zone could significantly impact its accuracy" — motivates the *radius* derefinement-protection criterion.

## 4. Equations actually present in the paper (all numerics, with labels as in source)

1. Lagrange interpolation weights, **unlabeled** (Sec. 3.1.1 Restriction):
   `ell_i(x) = prod_{0<=m<=k, m != j} (x - x_m)/(x_i - x_m)`
   (as printed; note the index mismatch m != j vs subscript i is in the source).
2. Interpolant, **unlabeled**: `f(x) = sum_{0<=i<=k} l_i(x) f(x_i)`.
3. Linear-wave error norm, **unlabeled** (Sec. 4.1):
   `L^1_RMS = sqrt( sum_{a,b in {x,y,z}} ( int_mesh dv |g_ab(T) - g_ab(0)| / Vol. )^2 )`.
4. Waveform error model, **unlabeled** (Sec. 4.3.1): `psi_NR(t, Dx) = psi_cont(t) + xi (Dx)^n`.
5. Residual between resolutions, **unlabeled**: `epsilon(Dx_1, Dx_2) = psi(t, Dx_1) - psi(t, Dx_2) = xi ((Dx_1)^n - (Dx_2)^n)`.
6. Convergence Q-factor, label **`eq:qfactor`** (the only labeled equation in the paper):
   `epsilon(Dx_1, Dx_2) / epsilon(Dx_2, Dx_3) ~= ((Dx_1)^n - (Dx_2)^n) / ((Dx_2)^n - (Dx_3)^n) =: Q_n`.

## 5. Units
- Masses/lengths in M, with "M measured in the ADM mass of the initial data" (Sec. 4.3.1). Geometric units G = c = 1 implied throughout (never stated explicitly in this paper).
