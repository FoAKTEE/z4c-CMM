# derivation.md — arXiv:2308.10361 (Ma et al., CCM for physical degrees of freedom)

Source: `ref-paper/arxiv-2308.10361/src/paper.tex`. Every equation environment
(equation/align/eqnarray/gather/multline, starred or not) is listed in document
order. Existing `\label{...}` names are preserved VERBATIM. Equation
environments with no label carry sequential ids `eq:auto-1` … `eq:auto-38`.
All content here is a literature-grounded transcription (evidence type:
citation/transcription), NOT an independent derivation. Equations inside
`\begin{comment}` blocks and `\mycomment{}` bodies in the tex are included
(they are equation environments in the source and are covered for
completeness; marked "[commented-out in tex]").

Index convention: `(env N, tex ~LXXXX)` gives the environment's document-order
index and approximate source line.

---

## Section II — Summary of the Cauchy evolution and its boundary conditions (`sec:GH`)

### eq:3+1_metric (env 1, ~L209)
```latex
ds^2=g_{\mu^\prime \nu^\prime}dx^{\prime\, \mu^\prime}dx^{\prime\, \nu^\prime}=(-\alpha^2+\beta^{i^\prime}\beta^{i^\prime}\gamma_{i^\prime j^\prime})dt^{\prime 2}
+2\beta^{i^\prime}\gamma_{i^\prime j^\prime}dx^{\prime \,j^\prime}dt^\prime+\gamma_{i^\prime j^\prime}dx^{\prime \,i^\prime}dx^{\prime \,j^\prime}
```
3+1 (ADM) decomposition of the spacetime metric in Cauchy (primed) coordinates: lapse α, shift β^i', spatial metric γ_i'j'.
Dependencies: none (root; standard ADM split).

### FOSH (env 2, ~L215)
```latex
\partial_{t^\prime} u^{\alpha^\prime}+\tensor[]{A}{^{k^\prime \alpha^\prime}_{\beta^\prime}}\partial_{k^\prime} u^{\beta^\prime}=F^{\alpha^\prime}
```
Vacuum Einstein equations R_μ'ν'=0 recast as a first-order symmetric hyperbolic (FOSH) generalized-harmonic system for u^α' = {g_μ'ν', Π_μ'ν', Φ_i'μ'ν'} (Lindblom et al. 2005).
Dependencies: eq:3+1_metric.

### eq:auto-1 (env 3, ~L221) — unlabeled
```latex
\tensor[]{e}{^{\hat{\alpha}^\prime}_{\mu^\prime}}s_{k^\prime}\tensor[]{A}{^{k^\prime\mu^\prime}_{\beta^\prime}}=v_{(\hat{\alpha}^\prime)}\tensor[]{e}{^{\hat{\alpha}^\prime}_{\beta^\prime}}
```
Definition of the left eigenvectors e^α̂'_β' and eigenvalues (characteristic speeds) v_(α̂') of the FOSH principal symbol contracted with the boundary normal s_k'; these define the characteristic fields u^α̂' = e^α̂'_β' u^β'.
Dependencies: FOSH, eq:def_s.

### eq:def_s (env 4, ~L225)
```latex
s^{t^\prime}=0,\qquad s^{k^\prime}=\frac{\gamma^{i^\prime k^\prime}\partial_{i^\prime}r^\prime}{\sqrt{\gamma^{i^\prime j^\prime}\partial_{i^\prime}r^\prime\partial_{j^\prime}r^\prime}}
```
Outward-directed spatial unit normal to the outer boundary of the Cauchy computational domain.
Dependencies: eq:3+1_metric.

### Bjorhus_bc (env 5, ~L229)
```latex
d_{t^\prime}u^{\hat{\alpha}^\prime}=D_{t^\prime}u^{\hat{\alpha}^\prime}+v_{(\hat{\alpha}^\prime)}\left(d_\perp u^{\hat{\alpha}^\prime}-\left. d_\perp u^{\hat{\alpha}^\prime}\right|_{\rm BC}\right)
```
Bjørhus method: on the boundary the time derivative of each incoming characteristic field is corrected by replacing its normal derivative d_⊥u^α̂' with the desired boundary value, leaving tangential derivatives unchanged.
Dependencies: eq:auto-1, eq:auto-2, eq:auto-3.

### eq:auto-2 (env 6, ~L233) — unlabeled
```latex
d_{t^\prime}u^{\hat{\alpha}^\prime}\equiv\tensor[]{e}{^{\hat{\alpha}^\prime}_{\beta^\prime}}\partial_{t^\prime}u^{\beta^\prime}, \qquad d_\perp u^{\hat{\alpha}^\prime}\equiv\tensor[]{e}{^{\hat{\alpha}^\prime}_{\beta^\prime}}s^{k^\prime}\partial_{k^\prime}u^{\beta^\prime}
```
Definitions of the characteristic-projected time derivative d_t' and normal derivative d_⊥ used in the Bjørhus condition.
Dependencies: eq:auto-1, eq:def_s.

### eq:auto-3 (env 7, ~L237) — unlabeled
```latex
D_{t^\prime}u^{\hat{\alpha}^\prime}\equiv\tensor[]{e}{^{\hat{\alpha}^\prime}_{\beta^\prime}}(-\tensor[]{A}{^{{k^\prime}\beta^\prime}_{\alpha^\prime}}\partial_{k^\prime}u^{\alpha^\prime}+F^{\beta^\prime})
```
Definition of D_t'u^α̂': the characteristic projection of the full right-hand side of the FOSH evolution equations (the value d_t'u^α̂' would take in the bulk).
Dependencies: FOSH, eq:auto-1.

### eq:bc_bjorhus (env 8, ~L347)
```latex
d_{t^\prime} u^{\hat{1}-}_{\mu^\prime \nu^\prime}=P^{{\rm P}\rho^\prime \tau^\prime}_{\mu^\prime \nu^\prime}\left[D_{t^\prime} u^{\hat{1}-}_{\rho^\prime \tau^\prime}-(\alpha+s_{j^\prime}\beta^{j^\prime})
\times(w_{\rho^\prime \tau^\prime}^{-}-\left.w_{\rho^\prime \tau^\prime}^{-}\right|_{\rm BC}-\gamma_2s^{ i^\prime}c_{i^\prime \rho^\prime \tau^\prime}^{3})\right]
```
The physical boundary condition on the incoming characteristic field u^1̂−: it drives the inward-propagating Weyl components w−_ρ'τ' toward their desired boundary value w−|_BC (supplied by CCM); c³ is a GH constraint field, γ₂ a constraint-damping parameter.
Dependencies: Bjorhus_bc, eq:auto-4, eq:wab_projection.

### eq:auto-4 (env 9, ~L352) — unlabeled
```latex
P^{{\rm P}\rho^\prime \tau^\prime}_{\mu^\prime \nu^\prime}\equiv\left(\tensor[]{P}{_{\mu^\prime}^{\rho^\prime}}\tensor[]{P}{_{\nu^\prime}^{ \tau^\prime}}-\frac{1}{2}P_{\mu^\prime \nu^\prime}P^{\rho^\prime \tau^\prime}\right)
```
Physical (transverse-traceless) projection operator built from the 2D projector P_μ'ν'.
Dependencies: eq:projection_operator.

### eq:projection_operator (env 10, ~L356)
```latex
P_{\mu^\prime \nu^\prime}=g_{\mu^\prime \nu^\prime}+n_{\mu^\prime}n_{\nu^\prime}-s_{\mu^\prime}s_{\nu^\prime}
```
Projector onto the 2D subspace orthogonal to both the hypersurface normal n_μ' and the boundary normal s_μ'.
Dependencies: eq:3+1_metric, eq:def_s.

### eq:wab_projection (env 11, ~L363)
```latex
w_{\rho^\prime \tau^\prime}^{-}=P^{{\rm P}\mu^\prime \nu^\prime}_{\rho^\prime \tau^\prime}(n^{\eta^\prime}+s^{\eta^\prime})(n^{ \alpha^\prime}+s^{\alpha^\prime})C_{\mu^\prime \eta^\prime \nu^\prime \alpha^\prime}
```
Inward-propagating components of the Weyl tensor C: the TT-projected double contraction with the outgoing null direction n+s; this is the physical incoming-radiation content at the boundary.
Dependencies: eq:auto-4, eq:projection_operator.

### eq:wab (env 12, ~L370)
```latex
w_{\rho^\prime \tau^\prime}^{-}=2(\psi_0^\prime\bar{m}_{\rho^\prime}\bar{m}_{\tau^\prime}+\bar{\psi}_0^\prime m_{\rho^\prime}m_{\tau^\prime})
```
w− rewritten in terms of the Cauchy Weyl scalar ψ₀′ and the angular tetrad vector m_μ'; this is the form into which the characteristic ψ₀ is inserted for CCM.
Dependencies: eq:wab_projection, eq:auto-5, eq:GH_psi0_def.

### eq:auto-5 (env 13, ~L374) — unlabeled
```latex
P_{\rho^\prime \tau^\prime}=m_{\rho^\prime}\bar{m}_{\tau^\prime}+m_{\tau^\prime}\bar{m}_{\rho^\prime}
```
Completeness identity expressing the 2D projector in terms of the complex null vector m and its conjugate.
Dependencies: eq:projection_operator, eq:GH_m.

### eq:GH_psi0_def (env 14, ~L378)
```latex
\psi_0^\prime=C_{\mu^\prime \nu^\prime \rho^\prime \tau^\prime}l^{\mu^\prime}m^{\nu^\prime}l^{\rho^\prime}m^{\tau^\prime}
```
Definition of the Newman-Penrose Weyl scalar ψ₀′ in the Cauchy (GH) tetrad — the ingoing transverse radiation component.
Dependencies: eq:gh_tetrad_l, eq:GH_m.

### eq:gh_tetrad_l, eq:gh_tetrad_k (env 15, ~L384; subequations parent label `eq:gh_tetrad`)
```latex
l^{\mu^\prime}=\frac{1}{\sqrt{2}}(n^{\mu^\prime}+s^{\mu^\prime})   % eq:gh_tetrad_l
k^{\mu^\prime}=\frac{1}{\sqrt{2}}(n^{\mu^\prime}-s^{\mu^\prime})   % eq:gh_tetrad_k
```
Unique Cauchy outgoing (l) and ingoing (k) null tetrad vectors fixed by the boundary-condition structure: built from the time normal n and boundary normal s.
Dependencies: eq:def_s.

### eq:GH_m (env 16, ~L390)
```latex
m^{\mu^\prime}l_{\mu^\prime}=0,\quad m^{\mu^\prime}k_{\mu^\prime}=0,\quad m^{\mu^\prime}\bar{m}_{\mu^\prime}=1
```
Conditions on the Cauchy angular tetrad vector m^μ'; m is fixed only up to a phase rotation m → m e^{iΘ}, on which w− does not depend.
Dependencies: eq:gh_tetrad_l.

---

## Section III — Summary of the characteristic evolution (`sec:CCE`)

### eq:BS_PFB (env 17, ~L411)
```latex
ds^2=-\left(e^{2\hat{\beta}}\frac{\hat{V}}{\hat{r}}-\hat{r}^2\hat{h}_{\hat{A}\hat{B}}\hat{U}^{\hat{A}}\hat{U}^{\hat{B}}\right)d\hat{u}^2-2e^{2\hat{\beta}}d\hat{u}d\hat{r}
-2\hat{r}^2\hat{h}_{\hat{A}\hat{B}}\hat{U}^{\hat{B}}d\hat{u}d\hat{x}^{\hat{A}}+\hat{r}^2\hat{h}_{\hat{A}\hat{B}}d\hat{x}^{\hat{A}}d\hat{x}^{\hat{B}}
```
Bondi-Sachs metric in partially flat Bondi-like (hatted) coordinates {r̂, x̂^Â, û} used by the SpECTRE characteristic system; degrees of freedom Ŵ, ĥ_ÂB̂, Û^Â, β̂.
Dependencies: none (root; Moxon et al. characteristic formulation).

### eq:BS_angular_gauge_condition (env 18, ~L416)
```latex
{\rm det}(\hat{h}_{\hat{A}\hat{B}})={\rm det}(q_{\hat{A}\hat{B}})=\sin^2\hat{\theta}
```
Bondi gauge condition: determinant of the angular metric equals that of the unit-sphere metric q_ÂB̂.
Dependencies: eq:BS_PFB.

### eq:auto-6 (env 19, ~L422) — unlabeled (subequations parent label `eq:falloff_partial_inertial`)
```latex
\lim_{\hat{r}\to\infty}\hat{W}=\mathcal{O}(\hat{r}^{-2}),\quad
\lim_{\hat{r}\to\infty}\hat{U}^{\hat{A}}=\mathcal{O}(\hat{r}^{-2}),\quad
\lim_{\hat{r}\to\infty}\hat{h}_{\hat{A}\hat{B}}=q_{\hat{A}\hat{B}}+\mathcal{O}(\hat{r}^{-1})
```
Required falloff rates of the hatted metric functions near future null infinity in partially flat Bondi-like coordinates.
Dependencies: eq:BS_PFB.

### eq:true_BS (env 20, ~L429)
```latex
\lim_{\hat{r}\to\infty}\hat{\beta}=\mathcal{O}(\hat{r}^{-1})
```
Additional condition that would promote partially flat Bondi-like coordinates to true Bondi-Sachs coordinates (asymptotically Minkowski, up to BMS); not imposed during evolution.
Dependencies: eq:auto-6.

### eq:null_wt (env 21, ~L443)
```latex
\left(\partial_{\underline{\lambda}}\right)^{\underline{a}}=\delta_{a^\prime}^{\underline{a}}\frac{n^{a^\prime}+s^{a^\prime}}{\alpha-\gamma_{i^\prime j^\prime}\beta^{i^\prime}s^{j^\prime}}
```
Outgoing null generator at the worldtube with affine parameter λ̲, built from the Cauchy normal vectors; converts the spacelike foliation into a null foliation.
Dependencies: eq:3+1_metric, eq:def_s.

### eq:BS_bondi_like (env 22, ~L450)
```latex
ds^2=-\left(e^{2\beta}\frac{V}{r}-r^2h_{AB}U^AU^B\right)du^2-2e^{2\beta}dudr
-2r^2h_{AB}U^Bdudx^A+r^2h_{AB}dx^Adx^B
```
Bondi-Sachs metric in (unhatted) Bondi-like coordinates {u, r, x^A}: same form as eq:BS_PFB but without asymptotic-flatness requirements on the metric functions.
Dependencies: eq:BS_angular_gauge_condition, eq:null_wt.

### eq:auto-7 (env 23, ~L456) — unlabeled
```latex
\lim_{r\to\infty}W=\mathcal{O}(r^{0}),\quad
\lim_{r\to\infty}U^A=\mathcal{O}(r^{0}),\quad
\lim_{r\to\infty}h_{AB}=\mathcal{O}(r^{0})
```
Relaxed falloff rates of the unhatted Bondi-like metric functions (finite limits allowed at null infinity).
Dependencies: eq:BS_bondi_like.

### eq:UA_asymptotic (env 24, ~L466)
```latex
U^A=U^{(0)A}+\mathcal{O}(r^{-1})
```
Definition of U^(0)A, the asymptotic (r→∞) value of the shift-like Bondi function U^A; removing it defines the partially flat coordinates.
Dependencies: eq:BS_bondi_like.

### eq:CCE_tetrad_m, eq:CCE_tetrad_k, eq:CCE_tetrad_l (env 25, ~L485; subequations parent label `eq:CCE_tetrad`)
```latex
m^{\mu}=-\frac{1}{\sqrt{2}r}\left(\sqrt{\frac{K+1}{2}}q^{\mu}-\frac{J}{\sqrt{2(1+K)}}\bar{q}^{\mu}\right)   % eq:CCE_tetrad_m
k^{\mu}=\sqrt{2}e^{-2\beta}\left[\delta^{\mu}_u-\frac{V}{2r}\delta^{\mu}_r+\frac{1}{2}\bar{U}q^{\mu}+\frac{1}{2}U\bar{q}^{\mu}\right]   % eq:CCE_tetrad_k
l^{\mu}=\frac{1}{\sqrt{2}}\delta^{\mu}_r   % eq:CCE_tetrad_l
```
Null tetrad of the characteristic (CCE) system in Bondi-like coordinates (applies equally with hatted variables in partially flat coordinates): m angular, k ingoing, l = ∂_r/√2 outgoing.
Dependencies: eq:BS_bondi_like, eq:CCE_all_variables:K, q_superA_expression.

### eq:CCE_all_variables:K (env 26, ~L495; subequations parent label `eq:CCE_all_variables`)
```latex
U\equiv U^{A}q_{A}, \quad J\equiv\frac{1}{2}q^{A}q^{B}h_{AB},
K\equiv\frac{1}{2}q^{A}\bar{q}^{B}h_{AB}=\sqrt{1+J\bar{J}}   % eq:CCE_all_variables:K
```
Spin-weighted Bondi scalars: U (shift), J (complex strain-like angular metric component), K (determinant-condition-fixed component, K=√(1+JJ̄)); corrects a typo in Moxon et al. Eq. (10e).
Dependencies: eq:BS_bondi_like, q_superA_expression.

### q_superA_expression, q_subA_expression (env 27, ~L503; subequations parent label `eq:q_hat_two_expressions`)
```latex
q^{A}\partial_{A}=-\partial_{\theta}-\frac{i}{\sin\theta}\partial_{\phi}   % q_superA_expression
q_{A}dx^{A}=-d\theta-i\sin\theta d\phi                                     % q_subA_expression
```
Complex dyad on the unit sphere in contravariant and covariant form.
Dependencies: none (root; dyad convention).

### eq:dyad_property (env 28, ~L509)
```latex
q^{A}q_{A}=0,\qquad q^{A}\bar{q}_{A}=2
```
Normalization/orthogonality identities of the complex dyad.
Dependencies: q_superA_expression.

### eq:psi0_CCE (env 29, ~L515)
```latex
\psi_0=\left(\frac{r\partial_{r}\beta-1}{4Kr}\right)\left[(1+K)\partial_{r}J-\frac{J^2\partial_{r}\bar{J}}{1+K}\right]+\frac{J(1+K^2)\partial_{r}J\partial_{r}\bar{J}}{8K^3}+\frac{1}{8K}\left[\frac{J^2\partial^2_{r}\bar{J}}{1+K}-(1+K)\partial_{r}^2J\right]-\frac{J\bar{J}^2(\partial_{r}J)^2+J^3(\partial_{r}\bar{J})^2}{16K^3}
```
Closed-form Bondi-like Weyl scalar ψ₀ in terms of radial derivatives of β, J, K; equally valid for ψ̂₀ with hatted variables. This is the quantity the characteristic system feeds back to the Cauchy boundary.
Dependencies: eq:CCE_tetrad_m, eq:CCE_all_variables:K.

---

## Section IV — Matching characteristic and Cauchy systems (`sec:CCM`)

### Sec. IV.A — Jacobians for CCM (`subsec:Jacobians_for_ccm`)

### eq:cauchy_null_radius_dependence (env 30, ~L702)
```latex
\begin{cases}
\underline{u}=t^\prime, &\\
\underline{x}^{\underline{A}}=\delta^{\underline{A}}_{A^\prime}x^{\prime\,A^{\prime}}, & \\
\underline{\lambda}=\underline{\lambda}(t^\prime,r^\prime).&
\end{cases}
```
Definition of null-radius (underlined) coordinates: time and angles inherited from Cauchy coordinates, radial coordinate replaced by the affine parameter λ̲ of the outgoing null rays.
Dependencies: eq:null_wt.

### eq:transformation_cauchy_null_radius (env 31, ~L711)
```latex
g_{\underline{\lambda}\underline{u}}=-1,\quad g_{\underline{\lambda}\underline{\lambda}}=0, \quad g_{\underline{\lambda}\underline{A}}=0, \quad g_{\underline{u}\underline{u}}=g_{t^{\prime}t^\prime},
g_{\underline{u}\underline{A}}=\delta_{\underline{A}}^{A^\prime} g_{t^{\prime}A^{\prime}}, \quad g_{\underline{A}\underline{B}}=\delta_{\underline{A}}^{A^\prime}\delta_{\underline{B}}^{B^\prime} g_{A^{\prime}B^{\prime}}
```
Metric components in null-radius coordinates: null-gauge conditions on the λ̲ components plus direct inheritance of the time/angular components from the Cauchy metric.
Dependencies: eq:cauchy_null_radius_dependence, eq:3+1_metric.

### eq:Jacobian_cauchy_null_radius (env 32, ~L716)
```latex
\frac{\partial(t^\prime,r^\prime,x^{\prime\, A^\prime})}{\partial(\underline{u},\underline{\lambda},\underline{x}^{\underline{A}})}=
\begin{pmatrix}
1& \partial_{\underline{\lambda}}t^\prime & 0 \\
0 & \partial_{\underline{\lambda}}r^\prime & 0 \\
0 & 0 & \delta_{\underline{A}}^{A^\prime}
\end{pmatrix}
```
Jacobian from null-radius to Cauchy coordinates.
Dependencies: eq:transformation_cauchy_null_radius.

### eq:auto-8 (env 33, ~L730) — unlabeled
```latex
r=\left[\frac{{\rm det}(g_{\underline{A}\underline{B}})}{{\rm det}(q_{\underline{A}\underline{B}})}\right]^{1/4}
```
Definition of the Bondi-like (areal) radius r from the determinant of the angular metric, enforcing the determinant gauge condition.
Dependencies: eq:BS_angular_gauge_condition, eq:transformation_cauchy_null_radius.

### eq:transformation_null_radius_bondi_like (env 34, ~L734)
```latex
\begin{cases}
u=\underline{u}, &\\
x^A=\delta^A_{\underline{A}}\underline{x}^{\underline{A}}, & \\
r=r(\underline{u},\underline{\lambda},\underline{x}^{\underline{A}}).&
\end{cases}
```
Null-radius → Bondi-like coordinate map: only the radial coordinate changes (λ̲ → areal r).
Dependencies: eq:auto-8, eq:cauchy_null_radius_dependence.

### eq:Jacobian_null_radius_bondi_like (env 35, ~L743)
```latex
\frac{\partial(u,r,x^{ A})}{\partial(\underline{u},\underline{\lambda},\underline{x}^{\underline{A}})}=
\begin{pmatrix}
1&0 & 0 \\
\partial_{\underline{u}}r & \partial_{\underline{\lambda}}r & \partial_{\underline{A}}r \\
0 & 0 & \delta_{\underline{A}}^{A}
\end{pmatrix}
```
Jacobian from null-radius to Bondi-like coordinates.
Dependencies: eq:transformation_null_radius_bondi_like.

### eq:Jacobian_bondi_like_null_radius (env 36, ~L753)
```latex
\frac{\partial(\underline{u},\underline{\lambda},\underline{x}^{\underline{A}})}{\partial(u,r,x^{ A})}=
\begin{pmatrix}
1&0 & 0 \\
-\partial_{\underline{u}}r/\partial_{\underline{\lambda}}r & (\partial_{\underline{\lambda}}r)^{-1} & -\delta^{\underline{A}}_{A}\partial_{\underline{A}}r/\partial_{\underline{\lambda}}r \\
0 & 0 & \delta^{\underline{A}}_{A}
\end{pmatrix}
```
Inverse Jacobian: Bondi-like to null-radius coordinates.
Dependencies: eq:Jacobian_null_radius_bondi_like.

### eq:du_xhat (env 37, ~L770)
```latex
\partial_u \hat{x}^{\hat{A}}=-\partial_A \hat{x}^{\hat{A}}U^{(0)A}
```
Evolution equation for the partially flat angular coordinates x̂^Â that removes the asymptotically constant part U^(0)A of U^A; solved numerically on the characteristic grid.
Dependencies: eq:UA_asymptotic.

### eq:transformation_bondi_like_inertial (env 38, ~L775)
```latex
\begin{cases}
\hat{u}=u, &\\
\hat{x}^{\hat{A}}=\hat{x}^{\hat{A}}(u,x^A), & \\
\hat{r}=r\hat{\omega}(u,x^A),&
\end{cases}
```
Bondi-like → partially flat Bondi-like coordinate map: angular relabeling plus conformal rescaling of the radius by ω̂ to restore the determinant condition.
Dependencies: eq:du_xhat, eq:omegahat, eq:BS_angular_gauge_condition.

### eq:omegahat (env 39, ~L784)
```latex
\hat{\omega}=\frac{1}{2}\sqrt{\hat{b}\bar{\hat{b}}-\hat{a}\bar{\hat{a}}}
```
Conformal factor of the angular map in terms of the spin-weighted Jacobian factors â, b̂.
Dependencies: eq:ahat.

### eq:ahat, eq:bhat (env 40, ~L788)
```latex
\hat{a}=\hat{q}^{\hat{A}}\partial_{\hat{A}}x^Aq_A   % eq:ahat  (spin-weight 2)
\hat{b}=\bar{\hat{q}}^{\hat{A}}\partial_{\hat{A}}x^Aq_A   % eq:bhat  (spin-weight 0)
```
Spin-weighted Jacobian factors of the angular map ∂_Â x^A contracted on both dyads.
Dependencies: q_superA_expression.

### eq:expand_angular_jacobian_nohat (env 41, ~L793)
```latex
\partial_{\hat{A}}x^A=
\frac{1}{4}
\begin{pmatrix} \hat{q}_{\hat{A}},\bar{\hat{q}}_{\hat{A}} \end{pmatrix}
\begin{pmatrix} \bar{\hat{a}} & \bar{\hat{b}} \\ \hat{b} & \hat{a} \end{pmatrix}
\begin{pmatrix} q^{A} \\ \bar{q}^{A} \end{pmatrix}
```
Dyad expansion of the angular Jacobian ∂_Â x^A; the middle-matrix determinant (with the 1/4) equals −ω̂².
Dependencies: eq:ahat, eq:omegahat, eq:dyad_property.

### eq:auto-9 (env 42, ~L810) — unlabeled
```latex
a=q^{A}\partial_{A}\hat{x}^{\hat{A}}\hat{q}_{\hat{A}}   % (spin-weight 2)
b=\bar{q}^{A}\partial_{A}\hat{x}^{\hat{A}}\hat{q}_{\hat{A}}   % (spin-weight 0)
```
Spin-weighted factors a, b of the inverse angular Jacobian ∂_A x̂^Â.
Dependencies: q_superA_expression.

### eq:oemga_b_a (env 43, ~L815)
```latex
\omega=\frac{1}{2}\sqrt{b\bar{b}-a\bar{a}}
```
Conformal factor ω(û, x̂^Â) associated with the inverse angular map.
Dependencies: eq:auto-9.

### eq:expand_angular_jacobian (env 44, ~L819)
```latex
\partial_{A}\hat{x}^{\hat{A}}=
\frac{1}{4}
\begin{pmatrix} q_{A},\bar{q}_{A} \end{pmatrix}
\begin{pmatrix} \bar{a} & \bar{b} \\ b & a \end{pmatrix}
\begin{pmatrix} \hat{q}^{\hat{A}} \\ \bar{\hat{q}}^{\hat{A}} \end{pmatrix}
```
Dyad expansion of the inverse angular Jacobian ∂_A x̂^Â.
Dependencies: eq:auto-9, eq:oemga_b_a, eq:dyad_property.

### eq:a_ahat_b_bhat (env 45, ~L836)
```latex
a=-\frac{\hat{a}}{\hat{\omega}^2}, \qquad b=\frac{\bar{\hat{b}}}{\hat{\omega}^2}
```
Relations between forward and inverse Jacobian factors implied by ∂_Â x^A ∂_A x̂^B̂ = δ_Â^B̂ at the same point.
Dependencies: eq:expand_angular_jacobian_nohat, eq:expand_angular_jacobian.

### eq:auto-10 (env 46, ~L840) — unlabeled
```latex
\omega\hat{\omega}=1
```
The two conformal factors are mutually inverse.
Dependencies: eq:a_ahat_b_bhat, eq:oemga_b_a, eq:omegahat.

### eq:Jacobian_bondi_like_inertial (env 47, ~L845)
```latex
\frac{\partial(\hat{r},\hat{x}^{\hat{A}},\hat{u})}{\partial(r,x^A,u)}=
\begin{pmatrix}
\hat{\omega} & r\partial_A\hat{\omega} & r\partial_u\hat{\omega} \\
0 & \partial_A \hat{x}^{\hat{A}} & \partial_u \hat{x}^{\hat{A}} \\
0 & 0 & 1
\end{pmatrix}
```
Jacobian from Bondi-like to partially flat Bondi-like coordinates.
Dependencies: eq:transformation_bondi_like_inertial.

### eq:jacobian_bondi_inertial (env 48, ~L855)
```latex
\frac{\partial(r,x^A,u)}{\partial(\hat{r},\hat{x}^{\hat{A}},\hat{u})}=
\begin{pmatrix}
\omega & r\delta^A_{\hat{A}}\partial_{A}\ln\omega & r\partial_u\ln\omega+rU^{(0)A}\partial_A\ln\omega \\
0 & \partial_{\hat{A}} x^{A} & U^{(0)A} \\
0 & 0 & 1
\end{pmatrix}
```
Inverse Jacobian (partially flat → Bondi-like), simplified using eq:du_xhat.
Dependencies: eq:Jacobian_bondi_like_inertial, eq:du_xhat, eq:auto-10.

### Sec. IV.B — Choice 1 (`sec:tetrad_transformation_Scenario_I`)

### eq:null_vector_CCE_rhat (env 49, ~L870)
```latex
l^{\hat{a}}=\frac{1}{\sqrt{2}}\left(\partial_{\hat{r}}\right)^{\hat{a}}=\frac{1}{\sqrt{2}}\left(\partial_{\underline{\lambda}}\hat{r}\right)^{-1}\delta_{\underline{a}}^{\hat{a}}\left(\partial_{\underline{\lambda}}\right)^{\underline{a}}
=\frac{1}{\sqrt{2}}e^{2\hat{\beta}}\delta_{\underline{a}}^{\hat{a}}\left(\partial_{\underline{\lambda}}\right)^{\underline{a}}
```
The characteristic outgoing null vector l^â expressed via the null-radius generator ∂_λ̲, using the Jacobians and β̂ = −½ln(∂_λ̲ r̂); shows l^â ∝ the Cauchy l^a'.
Dependencies: eq:CCE_tetrad_m, eq:Jacobian_bondi_like_null_radius, eq:jacobian_bondi_inertial, eq:auto-11.

### eq:auto-11 (env 50, ~L875) — unlabeled
```latex
\hat{\beta}=-\frac{1}{2}\ln(\partial_{\underline{\lambda}}\hat{r})
```
Relation between the Bondi lapse-like function β̂ and the affine rate of change of the areal radius [Eqs. (19a), (33a) of Moxon et al. 2020].
Dependencies: none (root; imported relation from Moxon:2020gha).

### eq:l_GH_and_l_CCE_transform (env 51, ~L889)
```latex
l^{a^\prime}=(\alpha-\gamma_{i^\prime j^\prime}\beta^{i^\prime}s^{j^\prime})e^{-2\hat{\beta}}l^{\hat{a}}\delta_{\hat{a}}^{a^\prime}
```
Explicit proportionality between the Cauchy and characteristic outgoing null vectors, fixing the boost factor between the two tetrads.
Dependencies: eq:null_vector_CCE_rhat, eq:null_wt, eq:gh_tetrad_l.

### eq:lorentz_transformation_I (env 52, ~L896)
```latex
\bm{l}\to \bm{l}, \quad \bm{k} \to \bm{k} + \bar{\kappa}\bm{m} + \kappa\bar{\bm{m}} + \kappa\bar{\kappa}\bm{l},
\bm{m}\to \bm{m} + \kappa\bm{l}, \quad \bar{\bm{m}} \to \bar{\bm{m}} + \bar{\kappa}\bm{l}
```
Type I (null rotation about l) Lorentz transformation of the tetrad; κ is a spin-weight-1 complex parameter.
Dependencies: none (root; standard NP tetrad freedom).

### eq:lorentz_transformation_II (env 53, ~L901)
```latex
\bm{l}\to A\bm{l}, \quad \bm{k}\to A^{-1}\bm{k},
\bm{m}\to e^{i\Theta}\bm{m}, \quad \bar{\bm{m}}\to e^{-i\Theta}\bar{\bm{m}}
```
Type II (boost + spin rotation) Lorentz transformation; A, Θ real.
Dependencies: none (root; standard NP tetrad freedom).

### eq:auto-12 (env 54, ~L909) — unlabeled
```latex
\hat{\psi}_0\to\hat{\psi}_0
```
ψ̂₀ is invariant under Type I transformations (no mixing with other Weyl scalars).
Dependencies: eq:lorentz_transformation_I, eq:GH_psi0_def.

### eq:lorentz_psi0_ii (env 55, ~L913)
```latex
\hat{\psi}_0\to A^{2}e^{2i\Theta}\hat{\psi}_0
```
ψ̂₀ rescaling under Type II transformations.
Dependencies: eq:lorentz_transformation_II, eq:GH_psi0_def.

### typeI_Ahat_to_Ap (env 56, ~L943)
```latex
\partial_{\hat{A}}=(\partial_{\hat{A}}x^A)\delta^{A^\prime}_A \partial_{A^\prime}+\partial_{\underline{\lambda}}\times\left[\frac{\partial_{\hat{A}}r}{\partial_{\underline{\lambda}}r}-\frac{\partial_{\underline{A}}r}{\partial_{\underline{\lambda}}r}(\partial_{\hat{A}}x^A) \delta_A^{\underline{A}}\right]
```
Angular basis vector ∂_Â expressed in Cauchy angular basis plus a term ∝ ∂_λ̲ (outgoing null), which is removable by a Type I transformation.
Dependencies: eq:Jacobian_cauchy_null_radius, eq:Jacobian_bondi_like_null_radius, eq:jacobian_bondi_inertial.

### eq:q_transformation_type_i (env 57, ~L950)
```latex
q^{\hat{\mu}}\approx \frac{1}{2}\hat{a}\delta_{\mu^\prime}^{\hat{\mu}}\bar{q}^{\mu^\prime}+\frac{1}{2}\bar{\hat{b}}\delta_{\mu^\prime}^{\hat{\mu}}q^{\mu^\prime}
```
Hatted dyad in terms of Cauchy dyad after dropping l-proportional terms via a Type I transformation (≈ denotes Type-I equivalence).
Dependencies: typeI_Ahat_to_Ap, eq:expand_angular_jacobian_nohat, eq:null_vector_CCE_rhat, eq:lorentz_transformation_I, q_superA_expression.

### eq:auto-13 (env 58, ~L956) — unlabeled
```latex
m^{\hat{\mu}}\approx-\frac{\delta^{\hat{\mu}}_{\mu^\prime}}{\sqrt{2}\hat{r}}\left[\left(\sqrt{\frac{\hat{K}+1}{2}}\frac{1}{2}\hat{a}-\frac{\hat{J}}{\sqrt{2(1+\hat{K})}}\frac{1}{2}\hat{b}\right)\bar{q}^{\mu^\prime}
+\left(\sqrt{\frac{\hat{K}+1}{2}}\frac{1}{2}\bar{\hat{b}}-\frac{\hat{J}}{\sqrt{2(1+\hat{K})}}\frac{1}{2}\bar{\hat{a}}\right)q^{\mu^\prime}\right]
```
Type-I-equivalent expression of the characteristic m^μ̂ in the Cauchy angular dyad basis.
Dependencies: eq:q_transformation_type_i, eq:CCE_tetrad_m.

### eq:new_m_cce_to_GH_abstract (env 59, ~L961)
```latex
m^{\hat{\mu}}\approx \delta^{\hat{\mu}}_{\mu^\prime}m^{\mu^{\prime}}
```
The transformed m differs from a purely Cauchy-angular vector m^μ' only by a Type I transformation.
Dependencies: eq:auto-13.

### eq:new_m_cce_to_GH (env 60, ~L965)
```latex
m^{a^{\prime}}= \hat{M}_{\theta^\prime}\left(\partial_{\theta^\prime}\right)^{a^\prime}+ \hat{M}_{\phi^\prime}\frac{i}{\sin\hat{\theta}(\theta^\prime)}\left(\partial_{\phi^\prime}\right)^{a^\prime}
```
Choice-1 Cauchy tetrad vector m^a' lying in the Cauchy angular subspace {θ', φ'}; satisfies eq:GH_m and is inserted (covariant form) into eq:wab.
Dependencies: eq:new_m_cce_to_GH_abstract, eq:auto-14.

### eq:auto-14 (env 61, ~L969) — unlabeled
```latex
4\hat{r} \hat{M}_{\theta^\prime}=(\hat{a}+\bar{\hat{b}})\sqrt{\hat{K}+1}-(\bar{\hat{a}}+\hat{b})\frac{\hat{J}}{\sqrt{(1+\hat{K})}},
4\hat{r}\hat{M}_{\phi^\prime}=(\bar{\hat{b}}-\hat{a})\sqrt{\hat{K}+1}-(\bar{\hat{a}}-\hat{b})\frac{\hat{J}}{\sqrt{(1+\hat{K})}}
```
Component coefficients M̂_θ', M̂_φ' of the Choice-1 m^a'.
Dependencies: eq:auto-13.

### eq:auto-15 (env 62, ~L976) — unlabeled
```latex
\left(\partial_{\theta^\prime}\right)^{a^\prime}=R^\prime_{\rm wt}\left(\cos\phi^\prime\cos\theta^\prime,\sin\phi^\prime\cos\theta^\prime,-\sin\theta^\prime\right),
\left(\partial_{\phi^\prime}\right)^{a^\prime}=R^\prime_{\rm wt}\sin\theta^\prime\left(-\sin\phi^\prime,\cos\phi^\prime,0\right)
```
Cartesian components of the Cauchy angular basis vectors at the worldtube radius R'_wt (needed because the Cauchy code uses Cartesian coordinates).
Dependencies: none (root; spherical-to-Cartesian basis).

### eq:auto-16 (env 63, ~L984) — unlabeled
```latex
\hat{A}=(\alpha-\gamma_{i^\prime j^\prime}\beta^{i^\prime}s^{j^\prime})e^{-2\hat{\beta}}
```
Type II Lorentz boost parameter Â relating l^â and l^a' for Choice 1.
Dependencies: eq:l_GH_and_l_CCE_transform, eq:lorentz_transformation_II.

### eq:psi0_prime_psi0_hat (env 64, ~L988)
```latex
\psi_0^{\prime}=\hat{A}^2\hat{\psi}_0
```
Choice-1 result: the Cauchy Weyl scalar is the hatted characteristic ψ̂₀ rescaled by Â²; the residual phase freedom Θ cancels in w−.
Dependencies: eq:auto-16, eq:lorentz_psi0_ii, eq:psi0_CCE.

### Sec. IV.C — Choice 2 (`sec:tetrad_transformation_Scenario_II`)

### eq:auto-17 (env 65, ~L1017) — unlabeled
```latex
\partial_A\approx\delta_A^{A^\prime}\partial_{A^\prime}
```
In Bondi-like coordinates the angular basis equals the Cauchy angular basis up to a Type-I-removable term.
Dependencies: eq:Jacobian_cauchy_null_radius, eq:Jacobian_bondi_like_null_radius, eq:lorentz_transformation_I.

### eq:new_m_cce_to_GH_B_abstract (env 66, ~L1021)
```latex
q^{\mu}\approx \delta^{\mu}_{\mu^\prime}q^{\mu^{\prime}}
```
The Bondi-like dyad is Type-I-equivalent to the Cauchy dyad.
Dependencies: eq:auto-17, q_superA_expression.

### eq:auto-18 (env 67, ~L1026) — unlabeled
```latex
m^\mu\approx\delta^{\mu}_{\mu^\prime}m^{\mu^{\prime}}
```
The Bondi-like m^μ is Type-I-equivalent to a Cauchy-angular vector m^μ'.
Dependencies: eq:new_m_cce_to_GH_B_abstract, eq:CCE_tetrad_m.

### eq:new_m_cce_to_GH_B (env 68, ~L1030)
```latex
m^{a^{\prime}}= M_{\theta^\prime}\left(\partial_{\theta^\prime}\right)^{a^\prime}+M_{\phi^\prime}\frac{i}{\sin\theta^\prime}\left(\partial_{\phi^\prime}\right)^{a^\prime}
```
Choice-2 Cauchy tetrad vector m^a' in the Cauchy angular subspace; satisfies eq:GH_m, used in eq:wab.
Dependencies: eq:auto-18, eq:auto-19.

### eq:auto-19 (env 69, ~L1034) — unlabeled
```latex
2r M_{\theta^\prime}=\sqrt{K+1}-\frac{J}{\sqrt{(1+K)}},
2r M_{\phi^\prime}=\sqrt{K+1}+\frac{J}{\sqrt{(1+K)}}
```
Component coefficients M_θ', M_φ' of the Choice-2 m^a'.
Dependencies: eq:auto-18.

### eq:J_BS_like, eq:K_BS_like (env 70, ~L1041)
```latex
J=\frac{\bar{b}^2\hat{J}+a^2\bar{\hat{J}}+2a\bar{b}\hat{K}}{4\omega^2}   % eq:J_BS_like
K=\sqrt{1+J\bar{J}}                                                      % eq:K_BS_like
```
Transformation of the evolved partially flat quantities Ĵ, K̂ back to the Bondi-like J, K (needed since the characteristic code evolves hatted variables).
Dependencies: eq:CCE_all_variables:K, eq:auto-9, eq:oemga_b_a.

### eq:auto-20 (env 71, ~L1052) — unlabeled
```latex
A=(\alpha-\gamma_{i^\prime j^\prime}\beta^{i^\prime}s^{j^\prime})e^{-2\beta}
```
Type II boost parameter A relating l^μ (Bondi-like) and l^μ' (Cauchy) for Choice 2.
Dependencies: eq:l_GH_and_l_CCE_transform, eq:lorentz_transformation_II.

### eq:auto-21 (env 72, ~L1056) — unlabeled
```latex
\psi_0^\prime=A^2\psi_0
```
Choice-2 result: Cauchy ψ₀′ from the Bondi-like ψ₀ via the Type II boost.
Dependencies: eq:auto-20, eq:lorentz_psi0_ii, eq:psi0_CCE, eq:J_BS_like.

### Sec. IV.D — Interpolating to the Cauchy coordinates (`sec:ccm_Cauchy_inertal`)

### eq:duhat_x (env 73, ~L1111)
```latex
\partial_{\hat{u}}x^A=U^{(0)A}
```
Evolution equation for the inverse angular map x^A(û, x̂^Â), read off the Jacobian eq:jacobian_bondi_inertial; evolved alongside eq:du_xhat to avoid expensive numerical inversion.
Dependencies: eq:jacobian_bondi_inertial.

### eq:auto-22 (env 74, ~L1115) — unlabeled
```latex
x^i=\left(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta\right)
```
Cartesian embedding of the Bondi-like angular coordinates on the unit sphere (spin-weight 0 quantities).
Dependencies: none (root; standard embedding).

### eq:auto-23 (env 75, ~L1119) — unlabeled
```latex
\eth x^i=q^BD_Bx^i, \qquad \bar{\eth}x^i=\bar{q}^BD_Bx^i
```
Spin-weighted (eth) derivatives of the Cartesian sphere coordinates; D_A is the covariant derivative of q_AB = (q_A q̄_B + q̄_A q_B)/2.
Dependencies: eq:auto-22, q_superA_expression.

### eq:duhat_x_cartesian (env 76, ~L1123)
```latex
\partial_{\hat{u}}x^i= \frac{1}{2}\hat{\mathcal{U}}^{(0)}\bar{\eth} x^i+\frac{1}{2}\bar{\hat{\mathcal{U}}}^{(0)}\eth x^i
```
Cartesian, spin-weight-friendly form of eq:duhat_x; the equation actually evolved numerically for the inverse angular map.
Dependencies: eq:duhat_x, eq:auto-23, eq:auto-24.

### eq:auto-24 (env 77, ~L1127) — unlabeled
```latex
\mathcal{U}^{(0)\hat{A}}\partial_{\hat{A}}x^B=U^{(0)B},\qquad
\mathcal{U}^{(0)}=\mathcal{U}^{(0)\hat{A}}q_{\hat{A}}
```
Definition of the auxiliary hatted-frame advection variable 𝒰^(0) corresponding to U^(0)A.
Dependencies: eq:duhat_x.

### eq:U_mathcal_U (env 78, ~L1132)
```latex
\mathcal{U}^{(0)}=\frac{1}{2\hat{\omega}^2}\left(\hat{\bar{b}}U^{(0)}-\hat{a}\bar{U}^{(0)}\right)
```
Explicit 𝒰^(0) in terms of U^(0) and the Jacobian factors â, b̂, ω̂.
Dependencies: eq:auto-24, eq:a_ahat_b_bhat.

### eq:auto-25 (env 79, ~L1136) — unlabeled
```latex
U^{(0)}=\frac{1}{2\omega^2}\left(\bar{b}\mathcal{U}^{(0)}-a\bar{\mathcal{U}}^{(0)}\right)
```
Inverse relation: U^(0) from 𝒰^(0) using the inverse-map factors a, b, ω.
Dependencies: eq:U_mathcal_U, eq:oemga_b_a.

### eq:du_xhat_cartesian (env 80, ~L1140)
```latex
\partial_u\hat{x}^{\hat{i}}=-\frac{1}{2}U^{(0)}\bar{\eth} \hat{x}^{\hat{i}}-\frac{1}{2}\bar{U}^{(0)}\eth  \hat{x}^{\hat{i}}
```
Cartesian version of eq:du_xhat for the forward angular map; evolved alongside eq:duhat_x_cartesian.
Dependencies: eq:du_xhat, eq:auto-23.

### eq:constraints_ca_cb (env 81, ~L1147) — [commented-out in tex: inside \mycomment{}]
```latex
C_a=a+\frac{\hat{a}}{\hat{\omega}^2}, \qquad C_b=\frac{\bar{\hat{b}}}{b\hat{\omega}^2}-1
```
Proposed CCM consistency constraints monitoring the mutual inverse relation eq:a_ahat_b_bhat between forward and inverse angular Jacobian factors. (Present only inside a \mycomment macro in the tex; not part of the published algorithm text.)
Dependencies: eq:a_ahat_b_bhat.

---

## Section V — Numerical Tests (`sec:tests`)

### eq:Teukolsky_metric (env 82, ~L1194)
```latex
ds^2=-dt^{\prime\,2}+(1+Af_{r^{\prime}r^{\prime}})dr^{\prime\,2}+2Bf_{r^{\prime}\theta^{\prime}}r^{\prime}dr^{\prime}d\theta^{\prime}
+2Bf_{r^{\prime}\phi^{\prime}}r^{\prime}\sin\theta^{\prime}dr^{\prime}d\phi^{\prime}+(1+Cf^{(1)}_{\theta^{\prime}\theta^{\prime}}+Af^{(2)}_{\theta^{\prime}\theta^{\prime}})r^{\prime\,2}d\theta^{\prime\,2}
+2(A-2C)f_{\theta^{\prime}\phi^{\prime}}r^{\prime\,2}\sin\theta^{\prime}d\theta^{\prime}d\phi^{\prime}
+(1+Cf^{(1)}_{\phi^{\prime}\phi^{\prime}}+Af^{(2)}_{\phi^{\prime}\phi^{\prime}})r^{\prime\,2}\sin^2\theta^{\prime}d\phi^{\prime\,2}
```
Perturbative Teukolsky-wave metric (Teukolsky 1982) on a flat background: linear quadrupolar (l=2, m=0) GW used as analytic reference in the small-amplitude test.
Dependencies: eq:auto-26, eq:Teukolsky_angular_basis, eq:auto-27.

### eq:auto-26 (env 83, ~L1203) — unlabeled (subequations parent label `eq:Teukolsky_ABC`)
```latex
A=3\left[\frac{F^{(2)}}{r^{\prime\,3}}+\frac{3F^{(1)}}{r^{\prime\,4}}+\frac{3F}{r^{\prime\,5}}\right],
B=-\left[\frac{F^{(3)}}{r^{\prime\,2}}+\frac{3F^{(2)}}{r^{\prime\,3}}+\frac{6F^{(1)}}{r^{\prime\,4}}+\frac{6F}{r^{\prime\,5}}\right],
C=\frac{1}{4}\left[\frac{F^{(4)}}{r^{\prime}}+\frac{2F^{(3)}}{r^{\prime\,2}}+\frac{9F^{(2)}}{r^{\prime\,3}}+\frac{21F^{(1)}}{r^{\prime\,4}}+\frac{21F}{r^{\prime\,5}}\right]
```
Radial amplitude functions A, B, C of the Teukolsky wave in terms of the free function F and its retarded-time derivatives.
Dependencies: eq:Gaussian_pulse_F, eq:Teukolsky_wave_outgoing.

### eq:Teukolsky_angular_basis (env 84, ~L1210)
```latex
f_{r^{\prime}r^{\prime}}=4\sqrt{\frac{\pi}{5}}Y_{20},\quad f_{r^{\prime}\theta^{\prime}}=2\sqrt{\frac{\pi}{5}}\partial_{\theta^{\prime}} Y_{20},
f_{r^{\prime}\phi^{\prime}}=0 ,\quad f^{(2)}_{\theta^{\prime}\theta^{\prime}}=-1,\quad f_{\theta^{\prime}\phi^{\prime}}=0,
f^{(1)}_{\theta^{\prime}\theta^{\prime}}=2\sqrt{\frac{\pi}{5}}\left(\partial_{\theta^{\prime}}^2-\cot\theta^\prime\partial_{\theta^{\prime}}-\frac{\partial^2_{\phi^{\prime}}}{\sin^2\theta^{\prime}}\right)Y_{20},
f^{(1)}_{\phi^{\prime}\phi^{\prime}}=-f^{(1)}_{\theta^{\prime}\theta^{\prime}},\quad f^{(2)}_{\phi^{\prime}\phi^{\prime}}=1-f_{r^{\prime}r^{\prime}}
```
Angular basis functions (l=2, m=0, even parity) of the Teukolsky-wave metric.
Dependencies: eq:auto-27.

### eq:auto-27 (env 85, ~L1220) — unlabeled
```latex
Y_{20}=\frac{1}{8}\sqrt{\frac{5}{\pi}}(1+3\cos2\theta^\prime)
```
The (l=2, m=0) spherical harmonic.
Dependencies: none (root; standard harmonic).

### eq:Gaussian_pulse_F (env 86, ~L1224)
```latex
F(u^\prime)=Xe^{-\frac{(u^\prime-r^\prime_c)^2}{\tau^2}}
```
Outgoing Gaussian pulse profile: retarded time u' = t'−r', center r'_c, width τ, amplitude X.
Dependencies: none (root; free-function choice).

### eq:Teukolsky_wave_outgoing (env 87, ~L1229)
```latex
F^{(n)}\equiv\left[\frac{d^nF(u^\prime)}{du^{\prime\,n}}\right]_{u^\prime=t^\prime-r^\prime}
```
Definition of the n-th retarded-time derivative of F entering eq:auto-26 (eq:Teukolsky_ABC).
Dependencies: eq:Gaussian_pulse_F.

### eq:gauge_constraint (env 88, ~L1252)
```latex
{\rm{Gauge~ constraint}}=\left\lVert \sqrt{\sum_{a=0}^3C_a^2}\right\rVert
```
Pointwise Euclidean L² norm (over the Cauchy grid) of the GH gauge constraint C_a [Eq. (40) of Lindblom et al. 2005]; numerical diagnostic.
Dependencies: none (root; imported GH constraint definition).

### eq:three_constraint (env 89, ~L1256)
```latex
{\rm{Three-Index~ constraint}}=\left\lVert \sqrt{\sum_{i=1}^3\sum_{a,b=0}^3C_{iab}^2}\right\rVert
```
L² norm of the GH three-index constraint C_iab [Eq. (26) of Lindblom et al. 2005]; numerical diagnostic.
Dependencies: none (root; imported GH constraint definition).

### eq:Teukolsky_h_20 (env 90, ~L1271)
```latex
[rh_{20}]_{\mathscr{I}^+}=\sqrt{\frac{6\pi}{5}}F^{(4)}
```
Perturbative prediction for the (l=2, m=0) strain harmonic at future null infinity; analytic reference for the X=1e-5 test.
Dependencies: eq:Teukolsky_metric, eq:Teukolsky_wave_outgoing.

### eq:bondi_violation_psi3 (env 91, ~L1323; subequations parent label `eq:bondi_violation`)
```latex
C_{\psi_4}\equiv\psi_4+\ddot{h}=0, \quad  C_{\psi_3}\equiv\psi_3-\frac{1}{\sqrt{2}}\eth \dot{h}=0,   % eq:bondi_violation_psi3
C_{\psi_s}\equiv\dot{\psi}_{s}+\frac{1}{\sqrt{2}}\eth \psi_{s+1}-\frac{3-s}{4}\bar{h}\psi_{s+2}=0 \quad (s=0,1,2)
```
Bondi-gauge constraints from the NP Bianchi identities relating Weyl scalars and the strain at scri; used as waveform-quality diagnostics.
Dependencies: none (root; imported from Moreschi 1986 / Iozzo et al. 2020).

### eq:im_psi2 (env 92, ~L1330)
```latex
C_{\rm{Im} \psi_2 }\equiv{\rm{Im}} \psi_2+{\rm{Im}} \left(\frac{1}{2}\eth^2h+\frac{1}{4}\bar{h}\dot{h}\right)=0
```
Reality condition on the Bondi mass aspect (Flanagan-Nichols); additional waveform-quality constraint.
Dependencies: none (root; imported constraint).

### eq:auto-28 (env 93, ~L1344) — unlabeled
```latex
F(v^\prime)=Xe^{-\frac{(v^\prime-r^\prime_c)^2}{\tau^2}}
```
Ingoing Gaussian profile in advanced time v' = t'+r' replacing F(u') for the Kerr-perturbation test, so the wave falls into the BH.
Dependencies: eq:Gaussian_pulse_F.

### eq:auto-29 (env 94, ~L1369) — unlabeled
```latex
\hat{J}(\hat{y},\hat{\theta},\hat{\phi})=
\begin{cases}
0, & \hat{y}\leq \hat{y}_{\min}, \\
\tensor[_{+2}]{Y}{_{20}}(\hat{\theta},\hat{\phi})\mathcal{J}(\hat{y}), &\hat{y}_{\min} \leq \hat{y} \leq \hat{y}_{\max}, \\
0, & \hat{y}\geq \hat{y}_{\max},
\end{cases}
```
Compactly supported spin-weight-2 initial data for Ĵ on the initial null slice (ŷ = 1 − 2R̂/r̂) for the characteristic-pulse-injection test.
Dependencies: eq:auto-30.

### eq:auto-30 (env 95, ~L1378) — unlabeled
```latex
\mathcal{J}(\hat{y})=4Z\frac{(\hat{y}_{\max}-\hat{y})(\hat{y}-\hat{y}_{\min})}{(\hat{y}_{\max}-\hat{y}_{\min})^2}e^{-\frac{(\hat{y}-\hat{y}_c)^2}{\tau^2}}
```
Radial profile of the injected pulse: polynomial envelope times Gaussian, amplitude Z, center ŷ_c, width τ.
Dependencies: none (root; free-data choice).

---

## Appendix material in tex `\begin{comment}` blocks (counted as equation environments by the source scan)

### eq:auto-31 (env 96, ~L1450) — unlabeled [commented-out in tex]
```latex
\Box^\prime x^{\prime\, \mu^\prime}=H^{\mu^\prime}
```
GH gauge condition: coordinates satisfy a wave equation with gauge source function H^μ'.
Dependencies: eq:3+1_metric.

### eq:auto-32 (env 97, ~L1457) — unlabeled [commented-out in tex]
```latex
d_{t^\prime}u^{\hat{\alpha}^\prime}+v_{(\hat{\alpha}^\prime)}d_\perp u^{\hat{\alpha}^\prime}=\tensor[]{e}{^{\hat{\alpha}^\prime}_{\beta^\prime}}(-\tensor[]{A}{^{i^\prime\beta^\prime}_{\alpha^\prime}}\tensor[]{P}{^{k^\prime}_{i^\prime}}\partial_{k^\prime}u^{\alpha^\prime}+F^{\beta^\prime})
```
FOSH system projected onto characteristic fields, split into normal and tangential parts (precursor of the Bjørhus replacement).
Dependencies: FOSH, eq:auto-1.

### eq:auto-33 (env 98, ~L1462) — unlabeled [commented-out in tex]
```latex
d_{t^\prime}u^{\hat{\alpha}^\prime}=D_{t^\prime}u^{\hat{\alpha}^\prime}
```
In the bulk (no boundary replacement) the Bjørhus condition reduces to the ordinary evolution.
Dependencies: Bjorhus_bc.

### eq:auto-34 (env 99, ~L1469) — unlabeled [commented-out in tex]
```latex
c_{i^\prime \rho^\prime \tau^\prime}^{3}=\partial_{i^\prime}g_{\rho^\prime \tau^\prime}-\Phi_{i^\prime \rho^\prime \tau^\prime}
```
The GH three-index constraint field appearing in eq:bc_bjorhus.
Dependencies: FOSH.

### eq:auto-35 (env 100, ~L1474) — unlabeled [commented-out in tex]
```latex
u^{\hat{1}-}_{\mu^\prime \nu^\prime}=\Pi_{\mu^\prime \nu^\prime}- s^{i^\prime}\Phi_{i^\prime \mu^\prime \nu^\prime}-\gamma_2g_{\mu^\prime \nu^\prime}
```
The incoming characteristic field u^1̂− of the GH system on which the physical boundary condition acts.
Dependencies: FOSH, eq:def_s.

### eq:cce_com:H (env 101, ~L1495; subequations parent label `eq:cce_com`) [commented-out in tex]
```latex
\hat{\beta}_{,\hat{r}}=S_{\beta}(\hat{J}),\quad
(\hat{r}^2\hat{Q})_{,\hat{r}}=S_Q(\hat{J},\hat{\beta}),\quad
\hat{U}_{,\hat{r}}=S_U(\hat{J},\hat{\beta},\hat{Q}),\quad
(\hat{r}^2\hat{W})_{,\hat{r}}=S_W(\hat{J},\hat{\beta},\hat{Q},\hat{U}),
(\hat{r}\hat{H})_{,\hat{r}}+L_H\hat{H}+L_{\bar{\hat{H}}}\bar{\hat{H}}=S_H(\hat{J},\hat{\beta},\hat{Q},\hat{U},\hat{W})   % eq:cce_com:H
```
The hierarchical radial hypersurface equations of the characteristic system, solved in order along each constant-û slice; the last (labeled eq:cce_com:H) determines Ĥ = ∂_û Ĵ for time stepping.
Dependencies: eq:BS_PFB, eq:CCE_all_variables:H.

### eq:CCE_all_variables:H (env 102, ~L1510; also contains a second occurrence of label eq:CCE_all_variables:K; subequations parent label `eq:CCE_all_variables`) [commented-out in tex]
```latex
\hat{U}\equiv\hat{U}^{\hat{A}}q_{\hat{A}},\quad
\hat{Q}\equiv\hat{r}^2e^{-2\hat{\beta}}q^{\hat{A}}\hat{h}_{\hat{A}\hat{B}}\partial_{\hat{r}}\hat{U}^{\hat{B}},\quad
\hat{J}\equiv\frac{1}{2}q^{\hat{A}}q^{\hat{B}}\hat{h}_{\hat{A}\hat{B}},
\hat{K}\equiv\frac{1}{2}q^{\hat{A}}\bar{q}^{\hat{B}}\hat{h}_{\hat{A}\hat{B}}=\sqrt{1+\hat{J}\bar{\hat{J}}},   % eq:CCE_all_variables:K (second occurrence)
\hat{H}=\partial_{\hat{u}}\hat{J}   % eq:CCE_all_variables:H
```
Hatted spin-weighted characteristic variables: Û, auxiliary Q̂, Ĵ, K̂, and the time derivative Ĥ = ∂_û Ĵ.
Dependencies: eq:BS_PFB, q_superA_expression.

### eq:auto-36 (env 103, ~L1520) — unlabeled [commented-out in tex]
```latex
q_{\hat{A}\hat{B}}=\frac{1}{2}(q_{\hat{A}}\bar{q}_{\hat{B}}+\bar{q}_{\hat{A}}q_{\hat{B}})
```
Unit-sphere metric expressed through the complex dyad.
Dependencies: q_superA_expression.

---

## Appendix A — Teukolsky wave in the perturbative limit (`app:1e-5`)

### eq:auto-37 (env 104, ~L1531) — unlabeled (subequations parent label `eq:Teukolsky_weyl_strain_news`)
```latex
\psi_0=-\sqrt{\frac{2\pi}{15}}\tensor[_{+2}]{Y}{_{20}}\left[(6\ddot{C}-3\ddot{A})+\frac{1}{2}r(3\dddot{B}+\dddot{A})\right],
\psi_1=\frac{1}{2}\sqrt{\frac{2\pi}{15}}\tensor[_{+1}]{Y}{_{20}}\left[r\dddot{A}+3\ddot{B}\right],
\psi_2=-\sqrt{\frac{\pi}{5}}Y_{20}\ddot{A},
\psi_3=\frac{1}{2}\sqrt{\frac{2\pi}{15}}\tensor[_{-1}]{Y}{_{20}}\left[r\dddot{A}-3\ddot{B}\right],
\psi_4=\sqrt{\frac{2\pi}{15}}\tensor[_{-2}]{Y}{_{20}}\left[(3\ddot{A}-6\ddot{C})+\frac{1}{2}r(3\dddot{B}+\dddot{A})\right],
N=-\sqrt{\frac{2\pi}{15}}\tensor[_{-2}]{Y}{_{20}}\left[(3\dot{A}-6\dot{C})+\frac{1}{2}r(3\ddot{B}+\ddot{A})\right],
h=-\sqrt{\frac{2\pi}{15}}\tensor[_{-2}]{Y}{_{20}}\left[(3A-6C)+\frac{1}{2}r(3\dot{B}+\dot{A})\right]
```
Analytic perturbative expressions for the Weyl scalars, News N, and strain h of the Teukolsky wave (from the Appendix of Teukolsky 1982).
Dependencies: eq:Teukolsky_metric, eq:auto-26, eq:auto-38.

### eq:auto-38 (env 105, ~L1542) — unlabeled
```latex
\tensor[_{-2}]{Y}{_{20}}=\tensor[_{+2}]{Y}{_{20}}=\frac{1}{4}\sqrt{\frac{15}{2\pi}}\sin^2\theta,
\tensor[_{-1}]{Y}{_{20}}=\tensor[_{+1}]{Y}{_{20}}=-\frac{1}{4}\sqrt{\frac{15}{2\pi}}\sin2\theta,
Y_{20}=\frac{1}{8}\sqrt{\frac{5}{\pi}}(1+3\cos2\theta)
```
Spin-weighted spherical harmonics ₛY₂₀ used in the analytic expressions.
Dependencies: none (root; standard harmonics).

### eq:Teukolsky_bulk_psi0 (env 106, ~L1550; subequations parent label `eq:Teukolsky_scri_weyl_strain_news`)
```latex
rh|_{\mathscr{I}^+}=\sqrt{\frac{6\pi}{5}}F^{(4)}\times\tensor[_{-2}]{Y}{_{20}},
rN|_{\mathscr{I}^+}=\sqrt{\frac{6\pi}{5}}F^{(5)}\times\tensor[_{-2}]{Y}{_{20}},
r\psi_4|_{\mathscr{I}^+}=-\sqrt{\frac{6\pi}{5}}F^{(6)}\times\tensor[_{-2}]{Y}{_{20}},
r^2\psi_3|_{\mathscr{I}^+}=\sqrt{\frac{6\pi}{5}}F^{(5)}\times\tensor[_{-1}]{Y}{_{20}},
r^3\psi_2|_{\mathscr{I}^+}=-\sqrt{\frac{9\pi}{5}}F^{(4)}\times Y_{20},
r^4\psi_1|_{\mathscr{I}^+}=\sqrt{\frac{27\pi}{10}}F^{(3)}\times \tensor[_{+1}]{Y}{_{20}},
r^5\psi_0=-\sqrt{\frac{27\pi}{10}}F^{(2)}\times\tensor[_{+2}]{Y}{_{20}}   % eq:Teukolsky_bulk_psi0
```
Asymptotic (scri) simplifications of the perturbative waveform quantities; the last line (eq:Teukolsky_bulk_psi0) holds throughout spacetime and shows ψ₀ ~ r⁻⁵, explaining the weakness of backscatter at large boundary radii.
Dependencies: eq:auto-37, eq:auto-26.

---

Registry totals: 106 equation environments; 74 tex-born labels; 38 auto labels (eq:auto-1 … eq:auto-38).
