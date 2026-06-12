# convention.md — arXiv:2308.10361 symbol table

All meanings are literature-grounded transcriptions from `ref-paper/arxiv-2308.10361/src/paper.tex`.

## Index and coordinate-frame conventions

| Symbol | Meaning |
|---|---|
| `i, j, k, ...` | 3D spatial coordinate indices |
| `μ, ν, ...` | 4D spacetime coordinate indices |
| `a, b, ...` | abstract (frame-independent) indices, used only where stated |
| primed (`x'^μ'`, `t'`, `r'`, `x'^A'`, `θ'`, `φ'`) | Cauchy coordinates used by the interior GH evolution |
| underlined (`u̲`, `λ̲`, `x̲^A̲`) | null-radius coordinates: intermediate system sliced by outgoing null congruences from the worldtube |
| unhatted Bondi (`u`, `r`, `x^A`) | Bondi-like coordinates: metric in Bondi-Sachs form, no asymptotic falloff requirements; angles equal Cauchy angles |
| hatted (`û`, `r̂`, `x̂^Â`, `θ̂`, `φ̂`) | partially flat Bondi-like coordinates used by the characteristic evolution; all metric components except β̂ asymptote to Minkowski form |
| `A, B` (capital, angular) | angular coordinate indices on the 2-sphere ({θ, φ} pairs) in each frame |

## Cauchy (GH) evolution — Sec. II

| Symbol | Meaning |
|---|---|
| `g_μ'ν'` | spacetime metric in Cauchy coordinates |
| `α` | lapse function |
| `β^i'` | shift vector |
| `γ_i'j'` | spatial metric of the 3+1 split |
| `R_μ'ν'` | Ricci tensor (vacuum: R = 0) |
| `u^α'` | collection of GH dynamical variables {g, Π, Φ} (50 components) |
| `Π_μ'ν'` | `α^{-1}(β^i' ∂_i' g_μ'ν' − ∂_t' g_μ'ν')`, time-derivative variable of the metric |
| `Φ_i'μ'ν'` | `∂_i' g_μ'ν'`, spatial-derivative variable |
| `A^k'α'_β'` | FOSH characteristic matrix |
| `F^α'` | FOSH lower-order source terms |
| `e^α̂'_β'` | left eigenvectors of `s_k' A^k'` defining the characteristic fields `u^α̂' = e^α̂'_β' u^β'` |
| `v_(α̂')` | characteristic speed (eigenvalue) of field α̂'; incoming if v < 0 |
| `s^k'` | outward unit normal to the outer boundary (spatial, `s^t' = 0`) |
| `n^μ'`, `n_μ'` | future-directed unit normal to the spatial hypersurface |
| `d_t' u^α̂'`, `d_⊥ u^α̂'` | characteristic-projected time and normal derivatives |
| `D_t' u^α̂'` | characteristic projection of the full FOSH right-hand side (bulk value of d_t') |
| `u^0̂, u^1̂±, u^2̂` | GH characteristic field groups; `u^1̂−` is the incoming group carrying physical + gauge modes |
| `c³_i'ρ'τ'` | GH three-index constraint field `∂_i' g_ρ'τ' − Φ_i'ρ'τ'` |
| `γ₂` | GH constraint-damping parameter |
| `P_μ'ν'` | projector onto the 2D subspace orthogonal to n and s |
| `P^P ρ'τ'_μ'ν'` | physical (TT) projection operator |
| `w−_ρ'τ'` | inward-propagating components of the Weyl tensor at the boundary |
| `w−|_BC` | desired boundary value of w−; zero in standard SpEC, CCM-supplied here |
| `C_μ'η'ν'α'` | Weyl tensor |
| `ψ₀'` | Cauchy-frame Weyl scalar `C l m l m` (ingoing transverse radiation) |
| `l^μ'` | Cauchy outgoing null tetrad vector `(n+s)/√2` (unique) |
| `k^μ'` | Cauchy ingoing null tetrad vector `(n−s)/√2` (unique) |
| `m^μ'`, `m̄^μ'` | Cauchy complex angular null tetrad vector; fixed up to phase `e^{iΘ}` |
| `Θ` | residual phase (spin-rotation) gauge freedom of m |
| `C_a` | GH gauge constraint (Eq. 40 of Lindblom et al. 2005) |
| `C_iab` | GH three-index constraint (Eq. 26 of Lindblom et al. 2005) |

## Characteristic (CCE) system — Sec. III

| Symbol | Meaning |
|---|---|
| `β̂, V̂, Û^Â, ĥ_ÂB̂, Ŵ` | hatted Bondi-Sachs metric functions in partially flat Bondi-like coordinates; `Ŵ = (V̂ − r̂)/r̂²` |
| `β, V, U^A, h_AB, W` | same metric functions in (unhatted) Bondi-like coordinates |
| `q_ÂB̂`, `q_AB` | unit-sphere metric (det = sin²θ); `q_AB = (q_A q̄_B + q̄_A q_B)/2` |
| `q^A`, `q_A` (and hatted versions) | complex dyad on the unit sphere: `q^A ∂_A = −∂_θ − (i/sinθ) ∂_φ`, `q_A dx^A = −dθ − i sinθ dφ` |
| `U` | spin-weighted shift scalar `U^A q_A` |
| `J` | complex spin-weight-2 scalar `½ q^A q^B h_AB` (strain-like angular metric component) |
| `K` | `½ q^A q̄^B h_AB = √(1+JJ̄)` (fixed by the determinant condition) |
| `Ĵ, K̂, Û` | hatted counterparts evolved by the characteristic code |
| `Q̂` | auxiliary radial variable `r̂² e^{−2β̂} q^Â ĥ_ÂB̂ ∂_r̂ Û^B̂` |
| `Ĥ` | `∂_û Ĵ`, time derivative of Ĵ on the null slice |
| `S_β, S_Q, S_U, S_W, S_H, L_H, L_H̄` | source/linear operators of the radial hypersurface hierarchy (Bishop et al. / Handmer et al.) |
| `U^(0)A` | asymptotic (r→∞) value of U^A; removed to define the hatted frame |
| `m^μ, k^μ, l^μ` (unhatted) | CCE null tetrad in Bondi-like coordinates; `l = ∂_r/√2` outgoing, k ingoing, m angular |
| `m^μ̂, l^â` etc. | same tetrad with hatted variables in partially flat coordinates |
| `ψ₀` / `ψ̂₀` | Weyl scalar assembled from Bondi-like / partially flat metric functions (eq:psi0_CCE) |
| `ℐ⁺` (`\mathscr{I}^+`) | future null infinity |

## Coordinate maps and Jacobians — Sec. IV

| Symbol | Meaning |
|---|---|
| `λ̲` | affine parameter along the outgoing null rays from the worldtube |
| `∂_λ̲` | outgoing null generator at the worldtube, `(n+s)/(α − γ_i'j' β^i' s^j')` |
| `r` (Bondi-like radius) | areal radius `[det(g_A̲B̲)/det(q_A̲B̲)]^{1/4}` |
| `ω̂(u, x^A)` | conformal factor of the Bondi-like → partially flat map, `r̂ = r ω̂`; `ω̂ = ½√(b̂b̄̂ − ââ̄)` |
| `ω(û, x̂^Â)` | conformal factor of the inverse map; `ωω̂ = 1` |
| `â` | spin-weight-2 Jacobian factor `q̂^Â ∂_Â x^A q_A` |
| `b̂` | spin-weight-0 Jacobian factor `q̄̂^Â ∂_Â x^A q_A` |
| `a`, `b` | inverse-map factors `q^A ∂_A x̂^Â q̂_Â` (sw 2), `q̄^A ∂_A x̂^Â q̂_Â` (sw 0); `a = −â/ω̂²`, `b = b̄̂/ω̂²` |
| `∂_Â x^A`, `∂_A x̂^Â` | forward/inverse angular Jacobians, expanded on the dyads |
| `κ` | spin-weight-1 complex parameter of Type I (null-rotation) Lorentz transformations |
| `A`, `Â` | real boost parameters of Type II Lorentz transformations (unhatted: Choice 2; hatted: Choice 1) |
| `M̂_θ', M̂_φ'` | components of the Choice-1 Cauchy m vector on (∂_θ', ∂_φ') |
| `M_θ', M_φ'` | components of the Choice-2 Cauchy m vector |
| `R'_wt` | Cauchy coordinate radius of the worldtube |
| `x^i` / `x̂^î` | Cartesian unit-sphere embeddings of the angular coordinates (spin-weight 0) |
| `ð, ð̄` (`\eth`) | spin-weighted angular derivatives `q^B D_B`, `q̄^B D_B` |
| `D_A` | covariant derivative of the unit-sphere metric q_AB |
| `𝒰^(0)`, `𝒰^(0)Â` | auxiliary hatted-frame advection variable with `𝒰^(0)Â ∂_Â x^B = U^(0)B` |
| `C_a, C_b` (CCM) | proposed Jacobian-consistency constraints (commented-out in tex; distinct from the GH gauge constraint C_a) |
| `≈` | equality up to a Type I Lorentz transformation (l-proportional terms dropped) |

## Numerical tests — Sec. V and appendices

| Symbol | Meaning |
|---|---|
| `A, B, C` (Teukolsky) | radial amplitude functions of the perturbative Teukolsky wave |
| `f_r'r', f_r'θ', f^(1,2)_θ'θ', ...` | angular basis functions of the (l=2, m=0) Teukolsky metric |
| `F(u')` | free Gaussian pulse profile; `u' = t' − r'` retarded time |
| `F(v')` | ingoing profile; `v' = t' + r'` advanced time (Kerr test) |
| `F^(n)` | n-th derivative of F evaluated at retarded time |
| `X` | Teukolsky wave amplitude (1e-5 perturbative; 2 nonlinear) |
| `r'_c`, `τ` | pulse center and width |
| `r'_out` | Cauchy outer boundary = worldtube radius (41 or 150) |
| `r'_ref` | outer boundary of the causally disconnected SpEC reference run (900) |
| `Y_20`, `ₛY₂₀` | (spin-weighted) l=2, m=0 spherical harmonics |
| `h`, `h_20` | GW strain and its (2,0) harmonic |
| `N` | Bondi News function |
| `ψ₀...ψ₄` | Newman-Penrose Weyl scalars at scri |
| `C_ψs`, `C_Imψ2` | Bondi-gauge constraint violations (Bianchi-identity and mass-aspect reality checks) |
| `χ` | dimensionless Kerr spin (0.5 in tests) |
| `r₊` | Kerr outer-horizon coordinate radius `1+√(1−χ²)` (spherical Kerr-Schild) |
| `ŷ` | compactified characteristic radial coordinate `1 − 2R̂/r̂` |
| `R̂` | partially flat Bondi-like worldtube radius |
| `𝒥(ŷ)`, `Z`, `ŷ_c`, `ŷ_min/max` | injected-pulse radial profile, amplitude, center, support |
| `t₁, t₂` | pulse arrival times at the Cauchy boundary (t₂ = t₁ + 2R crossing) |

## Naming notes (cross-tex resolution)

- Single tex file; no sibling tex conventions to merge. Ref. Lindblom:2005qh uses
  `ψ_ab` for the spacetime metric and `g_ij` for the spatial metric; this paper
  renames them to `g_μ'ν'` and `γ_i'j'` (footnote in Sec. II).
- The paper notes a typo in Eq. (10e) of Moxon:2020gha; the corrected K is
  eq:CCE_all_variables:K.
- `eq:CCE_all_variables:K` appears twice in the tex (live Sec. III and a
  commented-out appendix block); both occurrences are registered.
