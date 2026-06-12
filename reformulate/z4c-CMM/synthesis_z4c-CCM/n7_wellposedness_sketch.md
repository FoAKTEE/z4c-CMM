# N7 — Well-posedness sketch for Z4c-CCM (frozen-coefficient Laplace–Fourier)

Status: **[SOLID] as a documented sketch with model-level verification**
(`scripts/verify_n7_lf_sketch.py` → `results/numerical/n7_lf_sketch_check.txt`);
the full proof program remains **[FUTURE]** — explicit open items below
require human sign-off as research obligations.

## Setup

Half-space frozen-coefficient IBVP at the worldtube (P1 Sec. IV machinery:
eq:lapl-four, eq:bc_general, sol:lf, sol:lf2; Kreiss theory). Laplace in t
(`s`, Re s > 0), Fourier tangentially (`ω_A`); interior x < 0, boundary x = 0.
Mode ansatz `u = ũ e^{st + λx}` (sol:lf); admissible interior modes have
Re λ > 0; the boundary condition is well-posed in the Kreiss sense when its
determinant on the admissible modes has no zeros with Re s > 0, uniformly
(uniform Kreiss condition ⇒ strong well-posedness with inhomogeneous data
estimated by the datum norm).

## The central reduction (verified at model level)

**CCM modifies only the inhomogeneous datum of rows 1–2 of the BC table
(`z4c_ccm_boundary_conditions.md`), never the boundary operator.** The
Kreiss/LF determinant condition is a property of the *operator + interior
symbol* alone — the datum appears on the right-hand side of the estimate.
Therefore every boundary-stability statement P1 proves for the homogeneous
problem transfers verbatim to Z4c-CCM with ψ₀-driven data:

| sector | P1 status (homogeneous) | Z4c-CCM consequence |
|---|---|---|
| constraint (Θ, Z), any order L | boundary stable for ALL L (P1 Sec. IV, four-dim compact manifold, frozen coeff.) | inherited UNCHANGED (CCM datum never enters: N4 orthogonality) |
| gauge + metric, spherical reduction, asympt.-harmonic shift | boundary stable (P1 Sec. V cascade) | inherited; CCM adds an inhomogeneous physical datum g(t) with estimate ‖u‖ ≤ C‖g‖ |
| physical rows at normal incidence (model) | determinant D(s) = 2s, uniform Kreiss constant (verified) | same operator ⇒ same D(s); CCM datum on the RHS |

Model-level verification (normal incidence, advected wave
[∂_t² − 2b∂_t∂_x − (1−b²)∂_x²]u = 0, the P1 eq:eomgammasA reduction):

- interior roots λ₊ = s/(1+b) (admissible, Re λ₊ > 0 ⟺ Re s > 0),
  λ₋ = −s/(1−b) (inadmissible);
- physical/CPBC operator (∂_t + (1+b)∂_x)^{L+1}:
  **D_L(s) = (2s)^{L+1} ≠ 0 for Re s > 0, all L; |D_L(s)|/|s|^{L+1} = 2^{L+1}
  uniformly** — uniform Kreiss condition at every absorbing order;
- naive Sommerfeld (∂_t + ∂_x): D(s) = s(2+b)/(1+b) ≠ 0 — also well-posed
  (its defect is reflection, N5, not ill-posedness — honest note).

## Assumptions (explicit)

1. Frozen coefficients; high-frequency limit (P1's own regime — no
   symmetrizer is known for the standard puncture gauge).
2. Linearization around the boundary state; smooth data.
3. Subluminal shift at the boundary: α − γ_ij β^i s^j > 0 (boost positivity;
   numerically enforced/verified in the zccm package tests).
4. Absorbing order L ≥ 1 (P1: L = 0 reflects, R ≈ 1).
5. ψ₀′(t) treated as GIVEN data on the Cauchy side (see open item 1).
6. Asymptotically harmonic shift (ε_α, ε_χ) = (1, ½) for the gauge-sector
   cascade (P1 Sec. V).

## Open items — [FUTURE], human sign-off requested

1. **Coupled Cauchy ↔ characteristic system.** ψ₀′ is a functional of the
   characteristic solution, itself driven by Z4c worldtube data: a feedback
   loop with finite causal delay (worldtube ↔ boundary separation). Sketch:
   iterate the maps over one light-crossing of the matching shell and seek a
   contraction; obstruction: the Bondi-like characteristic system is only
   WEAKLY hyperbolic (P3's own caveat), so the required characteristic-side
   estimate is not available. OPEN even for GH-CCM.
2. **Oblique incidence + curvature** for the physical and gauge rows
   (O-N6-1): P1's full-system results are spherical-reduction only; the 3D
   discarded-tangential-terms caveat (P1 "possible bottleneck") applies.
3. **Corner/edge compatibility** between BC sectors on realistic boundary
   topologies (O-N6-4).
4. **Nonlinear constraint feedback** of the injected datum (O-N6-2 / N9):
   the N4 orthogonality is linear; the κ₁/κ₂-damped absorption argument for
   quadratic remainders is conjectural.

## Verifier

`scripts/verify_n7_lf_sketch.py`: exact sympy LF roots and determinants for
the model rows (admissibility, D_L(s) = (2s)^{L+1} for L = 0..3, Sommerfeld
comparison, datum-independence of D), plus a 4-GPU sweep of |D(s)|/|s|^{L+1}
over random s in the right half-plane confirming the uniform bound 2^{L+1}.
