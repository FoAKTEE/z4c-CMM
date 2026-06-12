# Z4c-CCM: the explicit boundary-condition set (formulation-DAG node N6)

The NEW FORMULATION. Worldtube `x = const` with outward unit normal `s^i`
(unit w.r.t. γ_ij), Cauchy normal `n^μ = (1, −β^i)/α`, tangential dyad
`m^A = (e_2 + i e_3)/√2`; null legs `l̊ = (n + s)/√2` (outgoing),
`k̊ = (n − s)/√2` (incoming) — paper-1 frame (P1 eq:nullvector-k,
eq:nullvector-m). `l̊·∂` annihilates outgoing waves and extracts incoming
content (verified, N3). Z4c variables per P1 eq:Z4_decomp_first
(χ, γ̃_ij, Ã_ij, K̂, Θ; Z_i carried by Γ̃^i; gauge α, β^i with 1+log +
gamma-driver, asymptotically harmonic shift ε_α = 1, ε_χ = 1/2).

Ten incoming characteristic fields at the boundary (P1 Sec. II; linearized
taxonomy verified in N4/N5: TT = physical, h_ss = gauge, h_sA = constraint-
vector, tangential trace = constraint-scalar):

| # | sector | incoming field | speed | boundary condition | datum source |
|---|---|---|---|---|---|
| 1–2 | physical | γ_AB^TF (TT tangential metric) | 1 | (l̊·∂)² γ_AB^TF ≐ 2(ψ₀′ m̄_A m̄_B + ψ̄₀′ m_A m_B) | **CCM: ψ₀′ = A² ψ₀_CCE** |
| 3 | constraint | Θ | 1 | (r² l̊·∂)^{L+1} Θ ≐ 0 | CPBC (P1 eq:general_CPBCs) |
| 4 | constraint | Z_s | 1 | (r² l̊·∂)^{L+1} Z_s ≐ 0 | CPBC |
| 5–6 | constraint | Z_A | 1 | (r² l̊·∂)^{L+1} Z_A ≐ 0 | CPBC |
| 7 | gauge | α (lapse) | √μ_L | (r² m̊·∂)^{L+1} α ≐ h_α | P1 eq:BCs-alpha (default h_α = 0) |
| 8 | gauge | β_s (normal shift) | √(ν_s) | (r² k̊_{ν_s}·∂)^{L+1} β_s ≐ h_s | P1 eq:auto-9/eq:Z4_betas_2+1BC form (default 0) |
| 9–10 | gauge | β_A (tangential shift) | √μ_S-sector | (r² j̊·∂)^{L+1} β_A ≐ h_A | P1 eq:BCs_lastII row 1 (default 0) |

with `m̊, k̊_{ν_s}, j̊` the per-speed outgoing characteristic vectors of P1
(`X̊^a = (n̊^a + √v · s̊^a)/√2` at sector speed v), and L ≥ 1 in production
(P1: L = 0 reflects, R ≈ 1).

## The physical rows (1–2) — the new content

Chain (all three links verified):

1. **CCE side (P2):** worldtube data from the Z4c 4-metric (node N1:
   `zccm.z4c_vars.four_metric` — exact maps; the χ-equation divergence is the
   full covariant `D_iβ^i`, resolved in `results/numerical/n1_varmap_check.txt`)
   feeds the hierarchy β→Q→U→W→H; ψ₀ on the characteristic grid from
   P2 eq:auto-48 / P3 eq:psi0_CCE.
2. **Frame transport (N2, verified):** the Type-II boost is built purely from
   3+1/characteristic data,
   `A = (α − γ_ij β^i s^j) e^{−2β̂}`, `ψ₀′ = A² ψ₀` —
   the tetrad correspondence `l̊ ∝ l̂` holds for arbitrary (α, β^i, γ_ij)
   with exact factor √2/2·e^{2β̂}/(α − γ_ij β^i s^j)
   (`results/numerical/n2_boost_check.txt`).
3. **Injection dictionary (N3, verified off-shell):**
   `(∂_t + ∂_s)² γ_AB^TF = 4(ψ₀ m̄_A m̄_B + c.c.)` with
   `ψ₀ = ¼ (∂_t + ∂_s)² h_mm`; since `(l̊·∂)² = ½(∂_t + ∂_s)²` on the flat
   frame, the BC row reads `(l̊·∂)²γ_AB^TF ≐ 2(ψ₀′ m̄m̄ + c.c.)`
   (`results/numerical/n3_dictionary_check.txt`).

**Bjørhus replacement form (implementation):** with the first-order incoming
field `V_AB := (l̊·∂) γ_AB^TF`, the boundary update is
`(l̊·∂) V_AB |_boundary = 2(ψ₀′ m̄_A m̄_B + c.c.)` — the Z4c analog of
P3 eq:bc_bjorhus, with GH's `−γ₂ s^i c³` constraint-damping addition replaced
by nothing at linear order (the CCM channel is constraint-orthogonal, N4);
the nonlinear analog (κ-weighted Θ, Z terms) is obligation O-N6-2.

## Limits and consistency (verified in scripts/verify_n6_composite.py)

- **ψ₀′ → 0** reproduces paper-1's freezing-ψ₀ / absorbing CPBC scheme
  exactly (Z4c-CCM ⊃ Z4c-CPBC).
- **Transparency:** pure outgoing waves satisfy the physical rows with zero
  datum — no spurious reflection of outgoing radiation (N3/N5).
- **Sector orthogonality at linear order:** the CCM datum drives rows 1–2
  only — rows 3–6 (N4) and 7–10 (N5) are blind to it.

## Communication cadence (from P3, sectors D5/E)

Each Cauchy RHS evaluation at the boundary: interpolate ψ₀ from the
characteristic grid (angular maps P3 eq:duhat_x / eq:du_xhat_cartesian +
time interpolation), boost (N2), build the datum (N3), apply rows 1–2;
rows 3–10 are local. Each CCE step: worldtube data from Z4c (N1).
Initialization: first CCE step uses the previous Cauchy slice's worldtube
data; ψ₀ injection starts one step delayed (P3's ordering carries over; the
β̂-circularity check is obligation O-N6-3).

## Open obligations (N7/N8/N9 program)

- O-N6-1: oblique incidence + curvature corrections to rows 1–2, 7–10
  (frozen-coefficient LF analysis — N7).
- O-N6-2: nonlinear constraint-damping addition to the physical Bjørhus form
  (κ₁/κ₂-weighted; Z4c analog of γ₂c³ — N9).
- O-N6-3: initialization order / β̂ circularity at t = 0 (P3 resolution
  transplant — N8).
- O-N6-4: corner/edge compatibility between sectors on a cubed-sphere or
  spherical-patch boundary (N7/N8).
