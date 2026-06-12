# Conventions — arXiv:2007.01339 (Moxon, Scheel, Teukolsky; improved CCE system)

Source: `ref-paper/arxiv-2007.01339/src/characteristic_formulation.tex` (single main tex; conventions converged against the in-paper Appendix A coordinate glossary, Table `tab:indices`). All content here is literature-grounded transcription of the paper's notation.

## Coordinate-system adornments (paper macros)

| Adornment | Macro | Coordinate system | Coordinates | Meaning |
|---|---|---|---|---|
| prime `'` | `\ca` | Cauchy coordinates | {t', r', x'^{A'}} | Input coordinates on the worldtube from the interior Cauchy simulation. |
| underline | `\rn` | Radial null coordinates | {u̲, λ̲, x̲^{A̲}} | Intermediate system with affine radial parameter λ along outgoing null rays; g_{λλ}=g_{λA}=0. |
| (none) | `\bl` | Bondi-like coordinates | {u, r, x^A} | Evolution coordinates of standard CCE: g_{rr}=g_{rA}=0, det(h_{AB})=det(q_{AB}); only finiteness (not Minkowski falloff) at scri. |
| hat `^` | `\ifc` | Partially flat ("incompletely flat") Bondi-like coordinates | {û, l̂ (or r̂), x̂^Â} | New intermediate gauge of this paper: Bondi-like conditions plus all metric components except β̂ asymptotically Minkowski (Ĵ^(0)=Û^(0)=0). |
| ring `∘` | `\bs` | Bondi-Sachs coordinates | {ů, r̊, x̊^Å} | True Bondi-Sachs gauge: Bondi-like conditions plus full Minkowski falloffs `eq:BondiSachsFalloffs`; unique up to BMS. |
| breve `˘` | `\na` | Numerically adapted coordinates | {ŭ, y̆, θ̆, φ̆} | Any Bondi-like system with compactified radial coordinate y̆ = 1 − 2R/r ∈ [−1, 1]. |

Index conventions: Greek = spacetime 4-indices; Latin i…n = spatial 3-indices; capital Latin = angular 2-indices; field symbols carry the same adornments as their coordinate system.

## Metric variables (Bondi-Sachs / Bondi-like)

| Symbol | Meaning |
|---|---|
| ds² | Spacetime line element in Bondi-Sachs form `eq:BondiSachsMetric` / Bondi-like form `eq:BondiLikeMetric`. |
| u | Retarded time labelling outgoing null hypersurfaces. |
| r | Areal radial coordinate (det of angular block fixed to unit sphere). |
| l | Inverse radial coordinate l ≡ 1/r; l = 0 at scri+ (used for asymptotic expansions). |
| x^A | Angular coordinate pair (θ, φ where explicit). |
| β | Bondi metric function: lapse-like factor in g_{ur} = −e^{2β} (NOT the Cauchy shift). |
| V | Bondi metric function in g_{uu}; V/r → 1 at scri in Bondi-Sachs gauge ("mass aspect" container). |
| W | Modified mass-aspect variable, r²W ≡ V − r (`WDef`); spin-weight 0. |
| U^A | Angular shift vector of the Bondi-like metric. |
| h_{AB} | Conformal 2-metric of the angular block, det(h_{AB}) = det(q_{AB}). |
| q_{AB} | Round unit-sphere metric. |
| g^{μν} | Up-index Bondi-like metric, `eq:UpIndexBondiLikeMetric`. |
| Γ | Worldtube: inner boundary 2+1 surface supplying Cauchy boundary data. |
| R(u,θ,φ) | Bondi-like radius of the worldtube, R = r|_Γ. |
| scri+, 𝓘⁺ | Future null infinity, outer boundary of the compactified domain. |

## Dyad and spin-weighted machinery

| Symbol | Meaning |
|---|---|
| q^A | Complex dyad on the unit sphere, normalization q^A q̄_A = 2; explicit choice q^A = {−1, −i/sinθ} (`eq:ComplexDyadChoice`). |
| spin-weight s | Quantity v with v → v e^{isψ} under dyad rotation q^A → q^A e^{iψ} (`eq:DyadRotation`). |
| ð, ð̄ | Spin-weight raising/lowering angular derivatives built from the unit-sphere covariant derivative D_A (`eq:DefEthEthbar`). |
| D_A | Covariant derivative of the unit-sphere metric q_{AB}. |
| U | U ≡ U^A q_A, spin-weight-1 shift scalar. |
| Q | Q ≡ r² e^{−2β} q^A h_{AB} ∂_r U^B (`eq:QDef`); spin-weight-1 auxiliary reducing the U equation to first order in r. |
| J | J ≡ ½ q^A q^B h_{AB}; spin-weight-2 scalar carrying the two gravitational degrees of freedom. |
| K | K ≡ ½ q^A q̄^B h_{AB} = √(1 + JJ̄) (`eq:KDef`); spin-weight 0 (determinant condition). |
| H | H ≡ ∂_u J, the sole evolved time-derivative quantity (worldtube/hypersurface variable). |
| v̄ (overbar) | Complex conjugation; flips spin-weight sign. |

## Boundary-transformation (Section II.B) quantities

| Symbol | Meaning |
|---|---|
| α | Lapse of the Cauchy ADM decomposition. |
| β^{i'} | Cauchy shift vector (distinct from Bondi β). |
| g_{i'j'} | Cauchy spatial 3-metric. |
| s^{α'} | Unit outward normal to the constant-r' extraction surface S_{r'}. |
| n^{α'} | Unit timelike hypersurface normal of the Cauchy slicing. |
| l^{α'} | Outgoing null vector generating the null cone at the worldtube, l = (n+s)/(α − g_{ij}β^i s^j). |
| λ̲ | Affine parameter along the null generators l^{α'}. |

## Gauge-transformation (Section III) quantities

| Symbol | Meaning |
|---|---|
| x̂^Â | Asymptotically inertial angular coordinates, evolved by `eq:XIfcEom`. |
| U^{(0)A} / U^{(R)A} | Asymptotic constant part / remainder of U^A in the expansion `eq:AsymptoticExpansionOfU`. |
| F^{(n)} | Coefficient of l^n in the asymptotic expansion F = F^(0) + l F^(1) + l² F^(2) + … near scri+. |
| ω̂ | Conformal factor of the partially flat radial rescaling, l̂ = l/ω̂ (`eq:conformal_factor_determinant`). |
| â, b̂ | Spin-weighted angular Jacobian factors â = q̂^Â ∂_Â x^A q_A, b̂ = q̄̂^Â ∂_Â x^A q_A (`eq:SpinWeightedJacobiansAB`); ω̂ = ½√(b̂b̄̂ − ââ̄). |
| 𝒰^{(0)} | Transformed asymptotic shift, 𝒰^(0) ≡ (b̄̂U^(0) − â Ū^(0))/(2ω̂²). |
| 𝒰 | Auxiliary full-shift quantity `eq:UHat` with ∂_r𝒰 = ∂_r Û; used to extract U^(0) during evolution. |
| ů^{(0)}, ů^{(R)}, ů^{(1)}, ů^{(2)} | Leading/remainder/expansion coefficients of the Bondi-Sachs time ů (step 3 of transformation). |
| x̊^{(R)Å}, x̊^{(1)Å}, x̊^{(2)Å} | Subleading angular coordinate corrections (step 4). |
| ω̊, ω̊^{(1)} | Conformal factor of the final areal-radius restoration l̊ = l̂/ω̊ (step 5) and its leading coefficient. |

## Compactified evolution (Section IV) quantities

| Symbol | Meaning |
|---|---|
| y̆ | Compactified radial coordinate y̆ = 1 − 2R/r ∈ [−1, 1]; worldtube at y̆ = −1, scri+ at y̆ = +1. |
| F̆ | Any spin-weighted scalar in numerical coordinates; F̆ = F for F ∈ {J, β, Q, U, W}, while H̆ = H + ∂_u R ∂_r J (`eq:numerical-h-conversion`). |
| S_β, S_Q, S_U, S_W, S_H, L_H | Schematic source/linear operators of the Bondi hierarchy `eq:BondiHierarchy`. |
| F_1, S_1 | Form-1 hypersurface equation (β, U): ∂_y̆F₁ = S₁ (`eq:form1`). |
| F_2, S_2^P, S_2^R | Form-2 equation (Q, W): (1−y̆)∂_y̆F₂ + 2F₂ = S₂^P + (1−y̆)S₂^R (`eq:form2`). |
| F_3, S_3^P, S_3^R, L_3^G, L_3^J | Form-3 equation (H) with conjugate coupling (`eq:form3`). |
| Λ_Q, Λ_W | Nonlinear source bundles in the Q and W hypersurface equations. |
| Λ_Q̆, Λ_W̆ | Their numerical-coordinate counterparts. |
| 𝒜_H, ℬ_H, 𝒞_H, 𝒟_H, ℒ_H | Source/linear bundles in the H hypersurface equation (`eq:script-b` block); 𝒞_H = ðβ − Q/2, 𝒟_H = ðŪ − ð̄U. |
| 𝒜_H̆, ℬ_H̆, 𝒞_H̆, 𝒟_H̆, ℒ_H̆ | Numerical-coordinate versions (replacement {J,β,Q,U,W} → breve quantities). |

## Asymptotics / news / Weyl quantities

| Symbol | Meaning |
|---|---|
| N_{ÅB̊} | Bondi news tensor, lim r̊ ∂_ů h_{ÅB̊} (`eq:BondiNews`). |
| N | Spin-weight-2 news scalar, N = ½ q̄^Å q̄^B̊ N_{ÅB̊} = lim r̊ ∂_ů J̊; primary waveform observable. |
| l^μ, n^μ, m^μ | Newman-Penrose null tetrad adapted to Bondi-like coordinates (`eq:NP_tetrads`). |
| D, Δ, δ | NP directional derivatives l^μ∇_μ, n^μ∇_μ, m^μ∇_μ (`eq:NPCovariantDerivatives`). |
| κ, ρ, σ, τ, ν, μ, λ, π, ε, β_NP, γ, α | NP spin coefficients (β_NP distinguishes the NP coefficient from Bondi β). |
| Θ | Unit-sphere dyad connection coefficient Θ = q^A q̄_B ∇_A q^B. |
| β_NP^SW, α^SW, γ^SW, Δ^SW, δ^SW | Spin-weighted versions with Θ set to zero / ð-substituted derivatives (coordinate-dependence isolated). |
| Ψ_0 … Ψ_4 | NP Weyl scalars; peeling Ψ_n ~ r^{n−5}. |
| Ψ_n^PF, Ψ_n^Bondi | Leading asymptotic Weyl scalars in partially flat / Bondi-Sachs gauges; Ψ_n^{PF(k)} = coefficient of r̂^{−k}. |
| R_{μν} | Ricci tensor; vacuum R_{μν} = 0 assumed for Weyl-scalar identities. |

## Other

| Symbol | Meaning |
|---|---|
| BMS | Bondi–Metzner–Sachs asymptotic symmetry group; residual gauge freedom of Bondi-Sachs coordinates. |
| CCE | Cauchy-characteristic evolution (extraction + characteristic evolution). |
| CCM | Cauchy-characteristic matching (future work). |
| Σ_u | Null hypersurface of constant retarded time u. |
| 𝒪(l^n) | Asymptotic order symbol in inverse radius near scri+. |
