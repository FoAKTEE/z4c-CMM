# convention.md — arXiv:1010.0523v2

Conventions resolved across the single source file `ref-paper/arxiv-1010.0523v2/src/Z4.tex`.
Geometric units c = G = 1; for matter spacetimes M_sun = 1, for black holes M_bh = 1.
Latin indices i,j,k,... are spatial (1–3); a,b are spacetime; upper case Latin A,B,C denote
quantities projected tangentially to the boundary 2-surface. Symmetrization with weight one:
`(ij)`; `|...|` excludes enclosed indices. `[ASSUMPTION]` signature (-,+,+,+).

## Macros
- `\p` — partial derivative ∂.
- `\hateq` — equality holding only on the boundary surface T = [0,T] × ∂Σ.
- `\Lie` — Lie derivative L.
- `\mbeta` = β̊ — normal component of the background shift at the frozen point p.
- `\re`, `\im` — real and imaginary parts.

## ADM / Z4 variables
- γ_ij — spatial 3-metric; γ = det γ_ij.
- K_ij — extrinsic curvature; K = γ^ij K_ij its trace.
- α — lapse function; β^i — shift vector.
- D_i — covariant derivative compatible with γ_ij; R_ij, R — its Ricci tensor/scalar.
- Θ — Z4 scalar constraint field (projection of Z_a along the normal); vanishes on GR solutions.
- Z_i — Z4 vector constraint field; vanishes on GR solutions.
- H — Hamiltonian constraint; M^i (M_i) — momentum constraint.
- ρ (ρ_ADM), S^i, S_ij, S — energy density, momentum density, stress tensor and its trace of matter.
- ∂_0 = (∂_t − β^i ∂_i)/α — derivative along the unit normal to the slices.
- □ — wave operator (d'Alembertian) in the principal part.

## Z4c conformal variables
- χ = γ^{-1/3} — conformal factor.
- γ̃_ij = γ^{-1/3} γ_ij — conformal metric (det γ̃ = 1 imposed).
- K̂ = γ^ij K_ij − 2Θ — modified trace of the extrinsic curvature.
- Ã_ij = γ^{-1/3}(K_ij − γ_ij K/3) — conformal trace-free extrinsic curvature.
- Γ̃^i = 2 γ̃^ij Z_j + γ̃^ij γ̃^kl γ̃_jk,l — conformal connection variable (absorbs Z_i).
- Γ̃_d^i = γ̃^ij γ̃^kl γ̃_jk,l — purely metric part of Γ̃^i.
- D̃_i — covariant derivative of γ̃_ij; Γ̃^k_ij its Christoffels; R̃_ij, R^χ_ij — conformal/χ parts of R_ij.
- [...]^tf — trace-free part with respect to γ_ij.
- D ≡ ln det γ̃ = 0, T ≡ γ̃^ij Ã_ij = 0 — algebraic constraints, imposed continuously.
- M̃^i — momentum constraint in conformal variables.

## Gauge
- μ_L — lapse-condition parameter (1+log slicing: μ_L = 2/α).
- μ_S — shift-condition parameter (standard puncture gauge: μ_S = 3/4; analysis uses μ_S = 1).
- ε_α, ε_χ — new shift parameters multiplying α∇α and γ̃^ij ∂_j χ terms; (μ_S, ε_α, ε_χ) = (1, 1, 1/2) = asymptotically harmonic shift; standard choice ε_α = ε_χ = 0.
- η — shift damping parameter (≈ 2/M for equal-mass evolutions).
- ν_s, ν_T — characteristic speeds of the shift equation in the scalar / vector sector.
- λ = sqrt(2(2γ^{1/3}μ_S − ε_χ)/(3α^2)) — scalar gauge characteristic speed; λ = 0 ⇒ weak hyperbolicity.

## 2+1 boundary decomposition
- s^i — outward spatial unit normal to the boundary (frozen analysis: s^i = −ê_x).
- q^i_j = δ^i_j − s^i s_j — projector tangential to the boundary surface.
- γ_ss, γ_qq, γ_sA, γ^TF_AB, β_s, β_A — scalar/vector/tensor projections of γ_ij, β_i.
- Λ = γ_ss + γ_qq — trace variable of the projected metric.

## Boundary geometry and conditions
- M = [0,T] × Σ; ∂Σ smooth boundary of the compact 3-manifold Σ; T = [0,T] × ∂Σ timelike boundary; Σ_t spacelike slices; S_t = {t} × ∂Σ.
- α̊, β̊_i, γ̊_ij, g̊_ab — background metric quantities; ψ̊ — background conformal factor; r — background isotropic radius; dΩ² — unit 2-sphere metric.
- n̊^a, s̊^a — background unit timelike normal to Σ_t and outgoing unit normal to S_t.
- l̊^a, k̊^a, j̊^a, m̊^a — background outgoing characteristic vectors at speeds 1, sqrt(ν_s), sqrt(ν_T), sqrt(μ_L).
- L — integer order of the constraint preserving boundary conditions (L = 0: first order; numerical tests use L = 1, "2nd order CPBCs").
- h_α, h_s, h_A, h^TF_AB, h, q — given boundary data; CPBC = constraint preserving boundary condition; IBVP = initial boundary value problem.

## Laplace–Fourier / Kreiss machinery
- s — Laplace frequency dual to t; η = Re(s) > 0; ω_A (ω_y, ω_z) — Fourier frequencies tangential to the boundary; ω² = ω_y² + ω_z².
- κ = sqrt(|s|² + ω²) (gauge subsystem section uses κ = |s|); s' = s/κ, ω' = ω/κ — normalized frequencies.
- ~ (tilde over a field) — Laplace–Fourier transform of the field.
- u, B, C^A — generic first order system state vector and coefficient matrices; Λ^I, Λ^II — incoming/outgoing blocks of B; m — number of incoming modes/BCs.
- M(s,ω) — companion (symbol) matrix of the reduced ODE in x; L(s,ω) — boundary matrix; g̃ — transformed boundary data; D(s,ω) — m×m boundary determinant matrix; Det D ≠ 0 — determinant condition.
- τ_±, ê_± (λ_i, e_i) — eigenvalues/eigenvectors of M; σ, σ_i — complex integration constants.
- DΘ̃, Dα̃, Dβ̃_s, DŨ — first order pseudo-differential reduction variables.
- γ = 1/sqrt(1 − β̊²), γ_α = 1/sqrt(2 − β̊²), γ_μ = 1/sqrt(μ² − β̊²) — boost-type normalization factors.
- λ² = s'² + γ^{-2} ω'² (context-dependent γ-factor) — boundary-symbol quantity (distinct from the gauge speed λ).
- a_± = s' ± λ (constraint sector); a_± = μ(s' ± λ) (toy model); a_+ = 2√2 s', b_+ = 2 s' (gauge subsystem) — eigenvalues of boundary matrices A.
- 𝓛, 𝓛_μ = (μ − β̊)s' − ∂_x/(κ γ_μ²) — pseudo-differential boundary operators; A, A^{L+1} — their algebraic matrix forms; F(β̊), G(β̊), H(β̊), J(β̊) — shorthand coefficient functions of β̊.
- H(s',ω') — Kreiss symmetrizer; C, C', C_T, δ, δ_2 — positive constants in estimates.
- U(t,x^i), F, μ, m — toy model field, source, wave speed, BC order.
- W, W̃ — reduced state vectors; f̃ — reduced source.

## Numerics
- C = sqrt(H² + M^i M_i + Θ² + Z^i Z_i) — constraint monitor; ||·||_2 — radial 2-norm.
- U^±_Θ — outgoing/incoming Θ characteristics; R = |Ũ^-_Θ(k)|/|Ũ^+_Θ(k)| — experimental reflection coefficient at wavenumber k.
- r_out, r'_out — outer boundary radii of test and reference runs; Δr — grid spacing; N — boundary grid index; f_i — grid function.
- γ̃_rr, γ̃_T, Ã_rr, Ã_T — spherical-symmetry metric/curvature variables.
- M — ADM (or black hole) mass; Kerr–Schild metric as in eq:auto-30.
