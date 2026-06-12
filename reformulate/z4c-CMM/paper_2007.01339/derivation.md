# Derivations — arXiv:2007.01339 (Improved Cauchy-characteristic evolution system)

Every equation environment of `ref-paper/arxiv-2007.01339/src/characteristic_formulation.tex`, in document order. Labels are preserved verbatim from the tex; unlabeled environments carry sequential ids `eq:auto-1` … `eq:auto-55`. Coordinate-adornment macros (`\bs` = Bondi-Sachs ring, `\bl` = Bondi-like, `\ifc` = partially-flat hat, `\na` = numerical breve, `\rn` = null-radius underline, `\ca` = Cauchy prime) are kept as in the source. All entries are literature-grounded transcriptions (evidence type: citation/transcription), not independent derivations.

Subequations-group labels that wrap the environments below (not themselves equation-environment labels): `eq:BondiSachsFalloffs`, `eq:BondiLikeFalloff`, `eq:NullRadiusMetric`, `eq:BondiLikeScalarsFromBondiLikeMetric`, `eq:BondiHierarchy`, `eq:GeodesicEquations`, `eq:SpinWeightedJacobiansAB`, `eq:ifcSwScalars`, `eq:numerical-jacobian`, `eq:numerical-evolution-source-terms`, `eq:BetaNearScri`, `eq:WNearScri`, `eq:EvolutionNearScri`, `eq:NP_tetrads`, `eq:NP_orthonormal`, `eq:spin_coefficients`, `eq:weyl_scalar_identities`, `eq:PartialExpandPsi24`, `eq:FullExpandPsi01`, `eq:if_weyl_scalars`, `eq:bondi_weyl_scalars`. They are noted with the entry that contains them.

---

## Section II.A — Bondi-Sachs metric and Bondi-like coordinates (sec:bondi_sachs)

### eq:BondiSachsMetric
```latex
ds^2 =
  -\left(e^{2 \bs \beta} \frac{\bs V}{\bs r} - \bs r^2 \bs h_{\bs A \bs B} \bs U^{\bs A} \bs U^{\bs B} \right) d\bs u^2
  - 2 e^{2 \bs \beta} d\bs u d\bs r
  - 2 \bs r^2 \bs h_{\bs A \bs B} \bs U^{\bs B} d \bs u d\bs x^{\bs A}
  + \bs r^2 \bs h_{\bs A \bs B} d \bs x^{\bs A} d\bs x^{\bs B}.
```
Bondi-Sachs form of an asymptotically flat metric in spherical null coordinates {ů, r̊, x̊^Å}, with gauge conditions g_{r̊r̊}=0, g_{r̊Å}=0. Root definition (literature: Bondi 1962).
Dependencies: none.

### eq:auto-1
```latex
\det(h_{\bs A \bs B}) = \det(q_{\bs A \bs B}).
```
Determinant gauge condition: the conformal angular metric has unit-sphere determinant, making r areal.
Dependencies: eq:BondiSachsMetric.

### eq:auto-2 (subequations group `eq:BondiSachsFalloffs`)
```latex
\lim_{\bs r\rightarrow\infty} \bs \beta(\bs x^{\bs \alpha}) = \mathcal{O}(\bs r^{-1}),\\
\lim_{\bs r\rightarrow\infty} \bs V(\bs x^{\bs \alpha}) = \bs r + \mathcal{O}(\bs r^0),\\
\lim_{\bs r\rightarrow\infty} \bs U^{\bs A}(\bs x^{\bs \alpha}) = \mathcal{O}(\bs r^{-2}),\\
\lim_{\bs r\rightarrow\infty} \bs h_{\bs A \bs B}(\bs x^{\bs \alpha}) = q_{\bs A \bs B}(\bs x^{\bs A}) + \mathcal{O}(\bs r^{-1}),
```
Asymptotic Minkowski falloff rates demanded of each Bondi-Sachs metric component at scri+.
Dependencies: eq:BondiSachsMetric.

### eq:BondiLikeMetric
```latex
ds^2 =
  -\left(e^{2 \bl \beta} \frac{\bl V}{\bl r} - \bl r^2 \bl h_{\bl A \bl B} \bl U^{\bl A} \bl U^{\bl B} \right) d\bl u^2
  - 2 e^{2 \bl \beta} d\bl u d\bl r
  - 2 \bl r^2 \bl h_{\bl A \bl B} \bl U^{\bl B} d \bl u d\bl x^{\bl A}
  + \bl r^2 \bl h_{\bl A \bl B} d \bl x^{\bl A} d\bl x^{\bl B}.
```
Bondi-like metric used by CCE: same local form and gauge conditions as Bondi-Sachs but without the asymptotic falloff restrictions.
Dependencies: eq:BondiSachsMetric.

### eq:auto-3 (subequations group `eq:BondiLikeFalloff`)
```latex
\lim_{\bl{r}\rightarrow\infty} \bl \beta(\bl{x}^{\bl{\alpha}}) = \mathcal{O}(r^0),\\
\lim_{\bl{r}\rightarrow\infty} \bl W(\bl{x}^{\bl{\alpha}}) = \mathcal{O}(r^0),\\
\lim_{\bl{r}\rightarrow\infty} \bl U^{\bl{A}}(\bl{x}^{\bl{\alpha}}) = \mathcal{O}(r^0),\\
\lim_{\bl{r}\rightarrow\infty} \bl h_{\bl{A} \bl{B}}(\bl{x}^{\bl{\alpha}}) = \mathcal{O}(r^0).
```
Relaxed Bondi-like asymptotics: every metric component is merely finite at scri+.
Dependencies: eq:BondiLikeMetric.

### eq:auto-4
```latex
q_{A B} = \frac{1}{2} (q_{A} \bar{q}_{B} + \bar{q}_{A} q_{B}).
```
Decomposition of the unit-sphere metric in the complex dyad q^A with normalization q^A q̄_A = 2.
Dependencies: none (definition).

### eq:ComplexDyadChoice
```latex
q^A = \left\{-1, \frac{-i}{\sin\theta}\right\}.
```
Explicit dyad choice in standard spherical angles {θ, φ}.
Dependencies: eq:auto-4.

### eq:DyadRotation
```latex
q^A\rightarrow q^A e^{i \psi}.
```
In-plane dyad rotation symmetry of the unit-sphere metric, defining the spin rotation.
Dependencies: eq:auto-4.

### eq:auto-5
```latex
v \rightarrow v e^{i s \psi}
```
Definition of spin-weight s: transformation law of a scalar under the dyad rotation.
Dependencies: eq:DyadRotation.

### eq:QDef, WDef, eq:KDef (one align; also defines U and J)
```latex
U \equiv U^A q_A, \quad\text{(spin-weight 1)}\\
Q \equiv r^2 e^{-2 \beta} q^A h_{A B} \partial_r U^B, \quad\text{(spin-weight 1)} \label{eq:QDef}\\
r^2 W \equiv V - r, \quad\text{(spin-weight 0)} \label{WDef}\\
J \equiv \frac{1}{2} q^A q^B h_{A B}, \quad\text{(spin-weight 2)}\\
K \equiv \frac{1}{2} q^A \bar{q}^B h_{A B} = \sqrt{1 - J \bar{J}}, \quad\text{(spin-weight 0)} \label{eq:KDef}
```
Standard spin-weighted scalars of the Bondi-like metric: U (shift), Q (first-order auxiliary for ∂_r U), W (mass-aspect variable, r²W = V − r), J (the two gravitational degrees of freedom), K (fixed by the determinant condition; note the printed sign in eq:KDef — the determinant condition gives K = √(1+JJ̄) as used everywhere else in the paper).
Dependencies: eq:BondiLikeMetric, eq:auto-4 (eq:KDef also: eq:auto-1).

### eq:auto-6 (subequations group `eq:DefEthEthbar`)
```latex
\eth v = q_1^{A_1} \dots q_n^{A_n} q^B D_B v_{A_1 \dots A_n},\\
\bar{\eth} v = q_1^{A_1} \dots q_n^{A_n} \bar{q}^B D_B v_{A_1 \dots A_n}.
```
Definition of the spin-weighted angular derivatives ð, ð̄ via the unit-sphere covariant derivative D_A.
Dependencies: eq:auto-4.

### eq:auto-7
```latex
\eth v = -(\sin \theta)^s \left(\frac{\partial}{\partial \theta}
  + \frac{i}{\sin \theta} \frac{\partial}{\partial \phi}\right)
  \left[(\sin \theta )^{-s} v\right].
```
Coordinate form of ð for a spin-weight-s scalar with the dyad eq:ComplexDyadChoice.
Dependencies: eq:auto-6, eq:ComplexDyadChoice.

## Section II.B — Boundary transformations (sec:boundary_transforms)

### eq:auto-8
```latex
ds^2 = \left(-\alpha^2 + \beta^{\ca i} \beta^{\ca j} g_{\ca i \ca j} \right)d\ca t^2
  + 2\beta^{\ca i} g_{\ca i \ca j} d\ca x^{\ca j} d\ca t
  + g_{\ca i \ca j} d\ca x^{\ca i} d\ca x^{\ca j}.
```
ADM decomposition of the Cauchy metric on the extraction surface (α lapse, β^i' shift, g_i'j' spatial metric).
Dependencies: none (input data).

### eq:auto-9
```latex
s^{\ca t} = 0, \qquad
s^{\ca i} = \frac{g^{\ca i \ca j} \partial_{\ca i} \ca r}{\sqrt{g^{\ca i \ca j} \partial_{\ca i} \ca r \partial_{\ca j} \ca r}}.
```
Unit spatial normal to the constant-r' surface S_r', radial part of the candidate null vector.
Dependencies: eq:auto-8.

### eq:auto-10
```latex
n^{\ca t} = \frac{1}{\alpha}, \qquad
n^{\ca i} = \frac{-\beta^{\ca i}}{\alpha}.
```
Timelike unit hypersurface normal of the Cauchy slicing.
Dependencies: eq:auto-8.

### eq:auto-11
```latex
l^{\alpha^\prime} = \frac{n^{\alpha^\prime} + s^{\alpha^\prime}}
  {\alpha - g_{i^\prime j^\prime} \beta^{i^\prime} s^{j^\prime}}
```
Normalized outgoing null generator of the null cones at the worldtube (Bishop et al. 1998).
Dependencies: eq:auto-9, eq:auto-10.

### eq:auto-12 (subequations group `eq:NullRadiusMetric`)
```latex
g_{\rn{\lambda} \rn{u}} = l^{\ca \alpha} g_{\ca \alpha \ca t} = -1,\\
g_{\rn{\lambda} \rn{\lambda}} = l^{\ca \alpha} l^{\ca \beta} g_{\ca \alpha \ca \beta} = 0,\\
g_{\rn{\lambda} \rn{A}} = l^{\ca \alpha} \delta_{\rn{A}}{}^{\ca A}
  \frac{\partial x^{\ca i}}{\partial \ca x{}^{\ca A}} g_{\ca \alpha \ca i} = 0,\\
g_{\rn{u} \rn{u}} = g_{\ca t \ca t},\\
g_{\rn{u} \rn{A}} = \delta_{\rn{A}}{}^{\ca A}
  \frac{\partial \ca x^{\ca i}}{\partial \ca x{}^{\ca A}} g_{\ca t \ca i},\\
g_{\rn{A} \rn{B}} = \delta_{\rn{A}}{}^{\ca A} \delta_{\rn{B}}{}^{\ca B}
  \frac{\partial x^{\ca i}}{\partial \ca x{}^{\ca A}}
  \frac{\partial x^{\ca i}}{\partial \ca x{}^{\ca B}} g_{\ca i \ca j}.
```
Metric components in null-radius coordinates {u̲, λ̲, x̲^A̲} (affine λ along l), enforcing g_λλ = g_λA = 0.
Dependencies: eq:auto-8, eq:auto-11.

### eq:BondiLikeRadius
```latex
\bl{r} = \left[\frac{\det(g_{\rn{A} \rn{B}})}{\det(q_{\rn{A} \rn{B}})}\right]^{1/4},
```
Areal Bondi-like radius constructed from the angular metric determinant.
Dependencies: eq:auto-12.

### eq:auto-13 (subequations group `eq:BondiLikeScalarsFromBondiLikeMetric`)
```latex
\beta = -\frac{1}{2} \ln(- g^{\bl{u} \bl{r}}) = -\frac{1}{2} \ln(\partial_{\rn{\lambda}} r),\\
U = \frac{g^{\bl{u} \bl{A}}}{g^{\bl{u} \bl{r}}} q_{\bl{A}},\\
W = \frac{1}{r}\left(1 - \frac{g^{r r}}{g^{\bl u \bl r}}\right),\\
J = -\frac{1}{2} r^2 q_{\bl{A}} q_{\bl{B}} g^{\bl{A} \bl{B}},\\
K = \sqrt{1 + J\bar{J}},\\
Q = \bl{r}^2 (J \partial_{\rn{\lambda}} \bar{U} + K \partial_{\rn{\lambda}} U).
```
Bondi-like spin-weighted scalars evaluated from the up-index metric in null-radius coordinates at the worldtube.
Dependencies: eq:BondiLikeRadius, eq:auto-12, eq:QDef, WDef, eq:KDef.

### eq:auto-14
```latex
U_{,\rn{\lambda}} = -\left(g^{\rn{\lambda} \rn{A}}{}_{,\rn{\lambda}}
  + \frac{\bl{r}_{,\rn{\lambda} \rn{B}}}{\bl{r}_{,\rn{\lambda}}} g^{\rn{A} \rn{B}}
  + \frac{\bl{r}_{,\rn{B}}}{\bl{r}_{\rn{\lambda}}} g^{\rn{A} \rn{B}} \right)
  + 2 \beta_{,\rn{\lambda}} \left(U + g^{\rn{\lambda} \rn{A}} q_{\rn{A}}\right),\\
J_{,\rn{\lambda}} = -\frac{1}{2} \bl{r}^2 q_{\rn A} q_{\rn B} h^{\rn{A} \rn{B}}{}_{,\rn{\lambda}}
  - \frac{2 \bl{r}_{,\rn{\lambda}}}{r} J.
```
Identities (Barkett et al. 2019) showing λ-derivatives of U and J need only first derivatives of the input metric.
Dependencies: eq:auto-13.

### eq:auto-15
```latex
\beta_{,\rn{\lambda}} = \frac{\bl{r}}{8 \bl{r}_{,\rn{\lambda}}}
  \left(J_{,\rn{\lambda}} \bar{J}_{,\rn{\lambda}} + \left(K_{,\rn{\lambda}}\right)^2\right).
```
β_,λ from an Einstein-equation component, eliminating second derivatives r_,λλ from the worldtube computation.
Dependencies: eq:auto-13, eq:auto-14.

### eq:auto-16
```latex
H \equiv J_{,\bl{u}} = -\frac{1}{2} \bl{r}^2 q_{\rn A} q_{\rn B} h^{\rn{A} \rn{B}}{}_{,\rn{u}}
  - \frac{2 \bl{r}_{,\rn{u}}}{r} J.
```
Worldtube value of H ≡ ∂_u J, the time-derivative input of the evolution.
Dependencies: eq:auto-13.

## Section II.C — Hierarchical evolution system (sec:EquationHierarchy)

### eq:BondiHierarchyBeta (align; subequations group `eq:BondiHierarchy`)
```latex
\bl \beta_{, r} = S_{\bl \beta}(\bl J), \label{eq:BondiHierarchyBeta}\\
(\bl r^2 \bl Q)_{,\bl r} = S_{\bl Q}(\bl J,\bl \beta),\\
U_{,\bl r} = S_{\bl U}(\bl J,\bl \beta,\bl Q),\\
(\bl r^2 \bl W)_{,\bl r} = S_{\bl W}(\bl J,\bl \beta,\bl Q,\bl U),\\
(\bl r \bl H)_{,\bl r} + L_{\bl H}(\bl J,\bl \beta,\bl Q,\bl U,\bl W)\bl H
  + L_{\bl{\bar H}}(\bl J,\bl \beta,\bl Q,\bl U,\bl W) \bl{\bar{H}}
  = S_{\bl H}(\bl J,\bl \beta,\bl Q,\bl U,\bl W).
```
Schematic hierarchical hypersurface system: each radial ODE depends only on previously solved scalars; H then steps J in u.
Dependencies: eq:BondiLikeMetric, eq:QDef, WDef, eq:KDef.

## Section II.D — Observables at scri+ (sec:FoundationsScri)

### eq:BondiNews
```latex
N_{\bs A \bs B} = \lim_{\bs r\rightarrow \infty} \left(\bs r\partial_{\bs u} h_{\bs A \bs B} \right).
```
Bondi news tensor: leading time derivative of the angular metric in true Bondi-Sachs coordinates.
Dependencies: eq:BondiSachsMetric.

### eq:auto-17
```latex
N = \frac{1}{2}\bs{\bar{q}}^{\bs A} \bs{\bar{q}}^{\bs B} N_{\bs A \bs B}
  = \lim_{\bs r \rightarrow \infty} \left(\bs r \partial_{\bs u} \bs J \right).
```
Spin-weighted news scalar; gauge-invariant up to BMS.
Dependencies: eq:BondiNews, eq:auto-4.

### eq:auto-18 (subequations group `eq:GeodesicEquations`)
```latex
\partial_{\bl{u}} \bs x^{\bs A} = - U^{\bl{B}} \partial_{\bl{B}} \bs x^{\bs A},\\
\partial_{\bl{u}} \bs u = \left( - U^{\bl{B}} \partial_{\bl{B}} \bs u + 1\right)e^{2 \bl \beta},
```
Geodesic equations evolving inertial coordinates on scri+ in the standard CCE algorithm.
Dependencies: eq:BondiLikeMetric.

## Section III — Transformations from Bondi-like to Bondi coordinates (sec:bondi_transforms)

### eq:UpIndexBondiLikeMetric
```latex
g^{\mu \nu} = \left[
\begin{matrix}
  0 & l^2 e^{-2 \beta} & 0 \\
  l^2 e^{-2 \beta} & l^3 e^{-2 \beta} (W + l) & l^2 e^{-2\beta} U^A\\
  0 & l^2 e^{-2 \beta} U^B & l^2 h^{A B}
\end{matrix} \right].
```
Up-index Bondi-like metric in inverse radial coordinate l = 1/r; the working object for the gauge-transformation steps.
Dependencies: eq:BondiLikeMetric, WDef.

### eq:AsymptoticExpansionOfU
```latex
U^A = U^{(0) A}(u, x^A) + l U^{(R) A}(u, l, x^A).
```
Split of the shift into its asymptotically constant part U^(0)A and remainder U^(R)A.
Dependencies: eq:BondiLikeMetric.

### eq:XIfcEom
```latex
\partial_u \ifc{x}^{\ifc{A}} = - \partial_A \ifc{x}^{\ifc{A}} U^{(0) A},
```
Step 1: evolution equation for asymptotically inertial angular coordinates x̂^Â removing U^(0)A.
Dependencies: eq:auto-18, eq:AsymptoticExpansionOfU.

### eq:conformal_factor_determinant
```latex
\ifc{l} = l \sqrt{\det(\partial_A \ifc{x}^{\ifc{A}})}
  \left(\frac{\det(\ifc{q}_{\ifc{A} \ifc{B}})}{\det(q_{A B})}\right)^{1/4}
  \equiv \frac{l}{\ifc \omega(u, x^A)}.
```
Step 2: areal inverse radius l̂ for the new angles, defining the conformal factor ω̂.
Dependencies: eq:XIfcEom, eq:auto-1.

### eq:auto-19 (subequations group `eq:SpinWeightedJacobiansAB`)
```latex
\ifc a = \ifc q^{\ifc A} \partial_{\ifc A} x^{A} q_{A},\\
\ifc b = \ifc{\bar{q}}^{\ifc A} \partial_{\ifc A} x^{A} q_{A}.
```
Spin-weighted angular Jacobian factors â, b̂ of the angular transformation.
Dependencies: eq:XIfcEom, eq:auto-4.

### eq:ifc_omega
```latex
\ifc \omega = \frac{1}{2}\sqrt{\ifc b \ifc{\bar{b}} - \ifc a \ifc{\bar{a}}}.
```
Identity expressing the conformal factor through the spin-weighted Jacobians.
Dependencies: eq:conformal_factor_determinant, eq:auto-19.

### eq:ifcBeta, eq:ifcW (one align; subequations group `eq:ifcSwScalars`, also gives Ĵ and Û)
```latex
\ifc \beta = \beta - \frac{1}{2} \log \ifc \omega \label{eq:ifcBeta},\\
\ifc J = \frac{\ifc{\bar b}^2 J + \ifc a^2 \bar J + 2 \ifc a \ifc{\bar b} K}{4 \ifc \omega^2},\\
\ifc U = \frac{1}{2\ifc \omega^2} \left(\ifc{\bar b} (U - U^{(0)}) - \ifc a (\bar U - \bar U^{(0)})\right)
  - \frac{\ifc l e^{2 \ifc \beta}}{\ifc \omega} \left(\ifc \eth \ifc \omega \ifc K - \ifc{\bar{\eth}} \ifc \omega \ifc J\right),\\
\ifc W = W + (\ifc \omega - 1) \ifc l - \frac{2 \partial_{\ifc u} \ifc \omega}{\ifc \omega}
  - \frac{1}{\ifc \omega} \left(\ifc U \ifc {\bar \eth} \ifc \omega + \ifc{\bar U} \ifc \eth \ifc \omega\right)
  + \frac{e^{2 \ifc \beta} \ifc l}{2 \ifc \omega^2}\left(\ifc{\bar \eth} \ifc \omega^2 \ifc J
  + \ifc \eth \ifc \omega^2 \ifc{\bar{J}} - 2 \ifc \eth \ifc \omega \ifc{\bar \eth} \ifc \omega \ifc K\right) \label{eq:ifcW}.
```
Partially flat spin-weighted scalars β̂, Ĵ, Û, Ŵ after transformation steps 1–2; the new Û vanishes at scri+.
Dependencies: eq:conformal_factor_determinant, eq:auto-19, eq:ifc_omega, eq:UpIndexBondiLikeMetric, eq:AsymptoticExpansionOfU.

### eq:auto-20
```latex
\partial_{\ifc u} \ifc \omega = \frac{\ifc \omega}{4}(\ifc \eth\,\bar{\mathcal{U}}^{(0)} + \ifc{\bar{\eth}} \, \mathcal{U}^{(0)} )
  + \frac{1}{2} \left(\mathcal U^{(0)} \ifc{\bar{\eth}} \ifc \omega + \bar{\mathcal{U}}^{(0)} \ifc \eth \omega\right),\\
\mathcal{U}^{(0)} \equiv \frac{1}{2 \ifc \omega^2} \left(\ifc{\bar b} U^{(0)} - \ifc a \bar U^{(0)} \right).
```
Time derivative of the conformal factor ω̂ inferred from the determinant expansion, with the transformed asymptotic shift 𝒰^(0).
Dependencies: eq:ifc_omega, eq:auto-19, eq:AsymptoticExpansionOfU.

### eq:auto-21
```latex
\bs{u} = \bs u^{(0)}(\ifc u, \ifc{x}^{\ifc A}) + \ifc{l} \bs u^{(R)}(\ifc u, \ifc{l}, \ifc{x}^{\ifc{A}})
  = \int^{\ifc u} e^{2 \ifc \beta} + \ifc{l} \bs u^{(R)}(\ifc u, \ifc{l}, \ifc{x}^{\ifc{A}}).
```
Step 3: asymptotically inertial Bondi-Sachs time ů, leading part fixed by integrating e^{2β̂}, remainder ů^(R) to be determined.
Dependencies: eq:ifcBeta.

### eq:uRequation
```latex
0 = g^{\bs u \bs u} = 2 \left(e^{2 \ifc \beta} + \ifc l \partial_{\ifc u} \bs u^{(R)} \right) \partial_{\ifc l} \left(\ifc{l} \bs u^{(R)}\right) \ifc l^2 e^{-2 \ifc \beta}
  + \ifc l^3 \partial_{\ifc l}(\ifc l \bs u^{(R)})^2 (\ifc W + \ifc l) e^{-2 \ifc \beta}
  + 2 \ifc l^2 \left(\partial_{\ifc{A}} \bs u^{(0)} + \ifc l \partial_{\ifc A} \bs u^{(R)}\right)
    \partial_{\ifc l} \left(\ifc l \bs u^{(R)}\right) \ifc U^{\ifc A} e^{-2 \ifc \beta}
  + \ifc{l}^2 (\partial_{\ifc A} \bs u^{(0)} + \ifc l \partial_{\ifc A} \bs u^{(R)})
    (\partial_{\ifc B} \bs u^{(0)} + \ifc l \partial_{\ifc B} \bs u^{(R)}) \ifc h^{\ifc A \ifc B},
```
Null condition g^{ůů} = 0 fixing ů^(R) (elliptic equation; solved perturbatively in l̂ near scri+).
Dependencies: eq:auto-21, eq:UpIndexBondiLikeMetric, eq:ifcW.

### eq:auto-22
```latex
\bs{x}^{\bs{A}} = \delta^{\bs{A}}{}_{\ifc{A}} x^{\ifc{A}}
  + \ifc{l} \bs{x}^{(R) \bs{A}} (\ifc{u}, \ifc{l}, \ifc{x}^{\ifc{A}}).
```
Step 4 ansatz: subleading-in-l̂ angular coordinate change, untouched at scri+.
Dependencies: eq:uRequation.

### eq:bsxA
```latex
0 = g^{\bs u \bs A} = \partial_{\ifc u} \bs u\, \partial_{\ifc l}(\ifc l \bs x^{(R) \bs A}) \ifc l^2 e^{-2 \ifc \beta}
  + \partial_{\ifc l} \bs u \,\partial_{\ifc u} \bs x^{(R) \bs A} \ifc l^3 e^{-2\ifc \beta}
  + \partial_{\ifc l} \bs u\, \partial_{\ifc l} (\ifc l \bs x^{(R) \bs A}) \ifc l^3 (\ifc W + \ifc l) e^{-2 \ifc \beta}
  + \ifc l^2 (\delta^{\bs A}{}_{\ifc A} + \ifc l \partial_{\ifc A} \bs x^{(R) \bs A})\, \partial_{\ifc l} \bs u\,\ifc U^{\ifc A} e^{-2 \ifc \beta}
  + \ifc l^2 \partial_{\ifc l}(\ifc l \bs x^{\bs A}) \partial_{\ifc A} \bs u \,\ifc U^{\ifc A} e^{-2 \ifc \beta}
  + \ifc l^2 \partial_{\ifc B} \bs u (\delta_{\ifc A}{}^{\bs A} + \ifc l \partial_{\ifc A} \bs x^{(R) \bs A}) h^{\ifc A \ifc B}.
```
Bondi condition g^{ůÅ} = 0 constraining x̊^(R)Å order-by-order in l̂.
Dependencies: eq:auto-22, eq:auto-21, eq:UpIndexBondiLikeMetric.

### eq:Step5RadialCoordinate
```latex
\bs{l} = \ifc l \sqrt{\det(\partial_{\ifc A} \bs{x}^{\bs{A}})}
  \left(\frac{\det(\bs{q}_{\bs{A} \bs{B}})}{\det(q_{\ifc A \ifc B})}\right)^{1/4}
  \equiv \frac{\ifc l}{\bs \omega(\ifc u, \ifc l, \ifc x^{\ifc A})}
  \equiv \ifc{l} - \ifc{l}^2 \bs \omega^{(1)} + \mathcal{O}(\ifc{l}^3).
```
Step 5: areal radius l̊ restoring the determinant condition for x̊^Å, completing the Bondi-Sachs frame (unique up to BMS).
Dependencies: eq:bsxA, eq:conformal_factor_determinant.

### eq:auto-23
```latex
\bs \beta^{(0)} = \bs \beta^{(1)} = \bs J^{(0)} = 0,\\
\bs U^{(0)} = \bs U^{(1)} = \bs W^{(0)} = \bs W^{(1)} = 0,\\
\bs J^{(1)} = \ifc J^{(1)} + \ifc \eth^2 \bs u^{(0)}.
```
Compact asymptotic results for the Bondi-Sachs scalars after steps 1–5 (simplified using perturbative Einstein equations).
Dependencies: eq:Step5RadialCoordinate, eq:uRequation, eq:bsxA.

## Section III.B — Inference of the news (sec:news)

### eq:NewsDefinitionBondi
```latex
N = \lim_{\bs r \rightarrow \infty} \frac{1}{2 \bs r} \bs{\bar{q}}^{\bs{A}} \bs{\bar{q}}^{\bs{B}} \partial_{\bs{u}}\bs h_{\bs A \bs B}.
```
Gauge-invariant Bondi-Sachs news: defined in every gauge by transforming to Bondi-Sachs and evaluating this limit.
Dependencies: eq:BondiNews, eq:auto-17.

### eq:NewsDefinitionIncompletelyFlat
```latex
N = e^{-2 \ifc \beta^{(0)}} \left(\ifc{\bar H}^{(1)} + \ifc{\bar \eth} \ifc{\bar \eth} e^{2 \ifc \beta^{(0)}}\right).
```
News in partially flat quantities, obtained by expanding eq:NewsDefinitionBondi through transformation steps 3–5.
Dependencies: eq:NewsDefinitionBondi, eq:auto-23.

### eq:NewsDefinitionBondiLike
```latex
N = \frac{\ifc \omega^2 e^{-2 \beta^{(0)}}}{4} \bigg\{
  2 \left[\ifc{\bar b} \eth \bar U^{(0)} J^{(1)} + \ifc a \eth U^{(0)} \bar J^{(1)}
  + \left(\ifc{\bar b} \eth U^{(0)} + \ifc a \eth \bar U^{(0)}\right)\text{Re}\left(J^{(0)} \bar{J}^{(1)} \right)\right]
  + \left[\ifc{\bar b}^2 H^{(1)} + \ifc a^2 \bar H^{(1)}
  + \ifc{\bar b} \ifc a\left(2 \text{Re}\left(J^{(0)} \bar{H}^{(1)} + J^{(1)} \bar H^{(0)}\right)
  - \text{Re}\left(J^{(0)} \bar{H}^{(0)}\right) \text{Re}\left(J^{(0)} \bar J^{(1)}\right) \right)\right]
  + \frac{1}{2} \left[\ifc{\bar b}^2 \left(U^{(0)} \bar \eth + \bar U^{(0)} \eth\right) J^{(1)}
  + \ifc a^2 \left(U^{(0)} \bar \eth + \bar U^{(0)} \eth\right) \bar{J}^{(1)}
  + \ifc a \ifc{\bar b} \left(U^{(0)} \bar \eth + \bar U^{(0)} \eth\right) \text{Re}\left(J^{(0)} \bar J^{(1)}\right) \right]
  + 3\ifc \omega^2 \partial_u \ifc \omega \left[\ifc{\bar b} J^{(1)} + \ifc a^2 \bar{J}^{(1)}
  + \ifc a \ifc{\bar b} \text{Re}\left(J^{(0)} \bar J^{(1)}\right)\right] \bigg\}
  + \frac{\ifc \omega e^{-2 \beta^{(0)}}}{4} \left(\ifc b^2 \eth^2 + \ifc{\bar a}^2 \bar \eth^2
  + 2 \ifc b \ifc{\bar a} \eth \bar \eth \right)\left(\frac{e^{2 \beta^{(0)}}}{\ifc \omega}\right).
```
News expressed in an arbitrary Bondi-like gauge via the explicit transformations of β̂ (eq:ifcBeta) and Ĥ (eq:ifcH); simpler than previously published forms.
Dependencies: eq:NewsDefinitionIncompletelyFlat, eq:ifcBeta, eq:ifcH, eq:auto-19, eq:auto-20.

## Section IV — Compactified characteristic evolution equations (sec:compactified_evolution)

### eq:auto-24
```latex
\na u = \bl{u}, \qquad \na y = 1 - \frac{2 R}{\bl{r}},
  \qquad \na \theta = \bl{\theta}, \qquad \na \phi = \bl{\phi},
```
Definition of numerically adapted coordinates: compactified radial coordinate y̆ ∈ [−1, 1] (worldtube y̆ = −1).
Dependencies: eq:BondiLikeMetric, eq:auto-25.

### eq:auto-25
```latex
R(u, \theta, \phi) = \bl{r}|_{\Gamma}.
```
Worldtube Bondi-like radius function R, generally angle- and time-dependent.
Dependencies: eq:BondiLikeRadius.

### eq:numerical-jacobian-u (align; subequations group `eq:numerical-jacobian`)
```latex
\partial_{\bl r} = \frac{(1 - \na y)^2}{2 R} \partial_{\na y},\\
\partial_{\bl u} = \partial_{\na u} - (1 - \na y) \frac{\partial_{\na u} R}{R} \partial_{\na y}, \label{eq:numerical-jacobian-u}\\
\bl \eth = \na \eth - (1 - \na y) \frac{\na \eth R}{R} \partial_{\na y}.
```
Inverse Jacobian relating Bondi-like derivatives to numerical-coordinate derivatives.
Dependencies: eq:auto-24.

### eq:numerical-h-conversion
```latex
\na H = H + \partial_u R \partial_r J.
```
Conversion of the evolution quantity: H̆ = ∂_ŭ J̆ in numerical coordinates differs from H by the moving-worldtube term.
Dependencies: eq:numerical-jacobian-u, eq:auto-16.

### eq:form1
```latex
\partial_{\na y} F_1(\na x^{\na \mu}) = S_1(\na x^{\na \mu}),
```
Category-1 hypersurface equation (β and U): plain radial ODE, integrable by standard methods.
Dependencies: eq:BondiHierarchyBeta, eq:numerical-jacobian-u.

### eq:form2
```latex
(1- \na y) \partial_{\na y} F_2(\na x^{\na \mu}) + 2 F_2(\na x^{\na \mu})
  = S_2^{P}(\na x^{\na \mu}) + (1 - \na y) S_2^{R}(\na x^{\na{\mu}}).
```
Category-2 equation (Q and W): degenerate at scri+; quadratic-in-(1−y̆) source parts threaten pure-gauge logarithms.
Dependencies: eq:BondiHierarchyBeta, eq:numerical-jacobian-u.

### eq:form3
```latex
(1 - \na y) \partial_{\na y} F_3(\na{x}^{\na \mu})
  + \big[1 + (1 - \na y)L_3^G(\na x^{\na \mu}) L_3^J(\na x^{\na \mu})\big] F_3(\na x^{\na \mu})
  + (1-\na y)\bar L_3^G(\na x^{\na \mu}) L_3^J(\na x^{\na \mu}) \bar{F}_3(\na x^{\na \mu})
  = S_3^P(\na x^{\na \mu}) + (1-\na y)S_3^R(\na x^{\na \mu}).
```
Category-3 equation (H only): degenerate and coupled to the complex conjugate F̄₃, requiring real/imaginary decomposition.
Dependencies: eq:BondiHierarchyBeta, eq:numerical-jacobian-u, eq:numerical-h-conversion.

### eq:auto-26
```latex
\left(\left[\begin{array}{cc}(1-y)\partial_y + 1 & 0\\ 0 & (1-y)\partial_y + 1\end{array} \right]
  + (1 - y) \left[\begin{array}{cc} \text{Re}(L_3^J)\text{Re}(L_3^G) & \text{Re}(L_3^J) \text{Im}(L_3^G)\\
  \text{Im}(L_3^J) \text{Re}(L_3^G) & \text{Im}(L_3^J) \text{Im}(L_3^G)\end{array}\right]\right)
  \left[\begin{array}{c}\text{Re}(F_3)\\\text{Im}(F_3)\end{array}\right]
  = \left[\begin{array}{c}\text{Re}(S_3^P) + (1 - y) \text{Re}(S_3^R)\\ \text{Im}(S_3^P) + (1 - y) \text{Im}(S_3^R)\end{array}\right].
```
Real/imaginary matrix form of eq:form3 used for spectral solution (Handmer & Szilagyi 2015).
Dependencies: eq:form3.

### eq:HypersurfaceBeta
```latex
\partial_r \beta = - \frac{1}{8} r \left(\partial_r J \partial_r \bar{J}
  - \frac{\left(\partial_r (J \bar{J})\right)^2}{4 K^2}\right).
```
β hypersurface equation in Bondi-like coordinates: β on Σ_u from J alone.
Dependencies: eq:BondiHierarchyBeta, eq:KDef.

### eq:Betanumeric
```latex
\partial_{\na y} (\na \beta) = -\frac{1}{8} (1 - {\na y}) \left(\partial_{\na y} \na J \partial_{\na y} \na{\bar{J}}
  - \frac{\left(\partial_{\na y} (\na J \na{\bar{J}})\right)^2}{4 \na K^2}\right).
```
Compactified β equation; form eq:form1, integrable by traditional methods.
Dependencies: eq:HypersurfaceBeta, eq:numerical-jacobian-u, eq:form1.

### eq:HypersurfaceQ
```latex
\partial_r (Q r^2) = - r^2 \left(\Lambda_Q + \frac{\bar{\Lambda}_Q J}{K} + \frac{\partial_r \bar{\eth} J}{K}\right)
  + 2 r^4 \partial_r \left(\frac{\eth \beta}{r^2}\right),
```
Simplified Q hypersurface equation in Bondi-like coordinates (Q from J, β).
Dependencies: eq:BondiHierarchyBeta, eq:QDef, eq:auto-27.

### eq:auto-27
```latex
\Lambda_Q = - \tfrac{1}{2} \eth (\bar{J} \partial_r J) + \tfrac{1}{2} J \partial_r \eth \bar{J}
  - \tfrac{1}{2} \eth \bar{J} \partial_r J + \frac{\eth (J \bar{J}) \partial_r (J \bar{J})}{4 K^2}.
```
Nonlinear source bundle Λ_Q of the Q equation.
Dependencies: eq:KDef, eq:auto-6.

### eq:Qnumeric
```latex
2 \na Q + (1 - {\na y}) \partial_{\na y} \na Q = -4 \eth \na \beta
  - (1 - {\na y})\left(2 \Lambda_{\na Q} + \frac{2 \bar{\Lambda}_{\na Q} \na J}{\na K}
  - 2 \eth \partial_{\na y} \na \beta + \frac{\bar{\eth} \partial_{\na y} \na J}{\na K}
  - \frac{2 \eth R \partial_{\na y} \na \beta}{R} + \frac{\bar{\eth} R \partial_{\na y} \na J}{R \na K}\right),
```
Compactified Q equation, form eq:form2; first/second right-hand terms are S^P_Q̆ and S^R_Q̆.
Dependencies: eq:HypersurfaceQ, eq:numerical-jacobian-u, eq:form2, eq:auto-28.

### eq:auto-28
```latex
\Lambda_{\na Q} = - \tfrac{1}{4} \eth (\na{\bar{J}} \partial_{\na y} \na J)
  + \tfrac{1}{4} \na J \eth \partial_{\na y} \na{\bar{J}}
  - \tfrac{1}{4} \eth \na{\bar{J}} \partial_{\na y} \na J
  + \frac{\eth (\na J \na{\bar{J}}) \partial_{\na y} (\na J \na{\bar{J}})}{8 K^2}
  + \frac{\eth^\prime (R)\left(\na J \partial_{\na y} \na{\bar{J}} - \na{\bar{J}} \partial_{\na y} \na J\right)}{4 R}.
```
Numerical-coordinate source bundle Λ_Q̆ (includes worldtube-radius angular-derivative correction).
Dependencies: eq:auto-27, eq:numerical-jacobian-u.

### eq:HypersurfaceU
```latex
\partial_r U = \frac{e^{2 \beta}}{r^2} \left(K Q - J \bar{Q}\right).
```
U hypersurface equation (angular-shift analogue) from J, β, Q.
Dependencies: eq:BondiHierarchyBeta, eq:QDef.

### eq:Unumeric
```latex
\partial_{\na y} \na U = \frac{e^{2 \na \beta}}{2 R} \left(\na K \na Q - \na J \na{\bar{Q}}\right).
```
Compactified U equation, form eq:form1.
Dependencies: eq:HypersurfaceU, eq:numerical-jacobian-u, eq:form1.

### eq:HypersurfaceW
```latex
\partial_r (r^2 W) = 1 + \tfrac{1}{2} e^{2 \beta} (\Lambda_W + \Lambda_W) + (\eth \bar{U} + \bar{\eth} U) r
  + \tfrac{1}{4} (\partial_r \eth \bar{U} + \partial_r \bar{\eth} U) r^2,
```
W (mass aspect) hypersurface equation from J, β, Q, U.
Dependencies: eq:BondiHierarchyBeta, WDef, eq:auto-29.

### eq:auto-29
```latex
\Lambda_W = - \eth \beta \eth \bar{J} + \tfrac{1}{2} \bar{\eth} \bar{\eth} J + 2 \bar{\eth} \beta \bar{\eth} J
  + (\bar{\eth} \beta)^2 J + \bar{\eth} \bar{\eth} \beta J
  + \frac{\eth (J \bar{J}) \bar{\eth} (J \bar{J})}{8 K^3} + \frac{2 + J \bar J}{2 K}
  - \frac{\eth \bar{\eth} (J \bar{J})}{8 K} - \frac{\eth (J \bar{J}) \bar{\eth} \beta}{2 K}
  - \frac{\eth \bar{J} \bar{\eth} J}{4 K} - \frac{\bar J \eth \bar{\eth} J}{4 K}
  - K \eth \bar{\eth} \beta - K \eth \beta \bar{\eth} \beta
  + \tfrac{1}{4} (- K Q \bar{Q} + J \bar{Q}^2).
```
Nonlinear source bundle Λ_W of the W equation.
Dependencies: eq:KDef, eq:auto-6.

### eq:Wnumeric
```latex
2 \na W + (1 - \na y) \partial_{\na y} \na W
  = \left(\eth \na{\bar{U}} + \bar{\eth} \na U\right)
  + (1 - \na y)\left(\tfrac{1}{4} \eth \partial_{\na y} \na{\bar{U}} + \tfrac{1}{4} \bar{\eth} \partial_{\na y} \na U
  + \frac{1}{4} \partial_{\na y} \na U \frac{\bar{\eth} R}{R} + \frac{1}{4} \partial_{\na y} \na{\bar{U}} \frac{\eth R}{R}
  - \frac{1}{2 R} + \frac{e^{2 \na \beta} (\Lambda_{\na W} + \bar{\Lambda}_{\na W})}{4 R}\right).
```
Compactified W equation, form eq:form2 (Λ_W̆ obtained by breve substitution in Λ_W).
Dependencies: eq:HypersurfaceW, eq:numerical-jacobian-u, eq:form2.

### eq:HypersurfaceH
```latex
\partial_r (H r) + J r (\mathcal{L}_H H + \bar{\mathcal{L}}_H H )
  = - \frac{1}{2} \eth ((J + r \partial_r J ) \bar{U}) + (\mathcal{B}_H + \bar{\mathcal{B}}_H) J
  - J\left( \eth \bar U + \frac{1}{2} \bar \eth U\right) - K \eth U
  - \frac{1}{2} \partial_r (r \bar{\eth} J ) U + \frac{1}{2} \partial_r \partial_r J (r^2 W + r)
  + \frac{e^{2 \beta}}{r} \left[\frac{1}{2} \bar{\eth} J K \mathcal{C}_H
  + \frac{\bar{\mathcal{C}}_H J^2 \eth \bar J}{2 K}
  - \left(\mathcal{A}_H + \bar{\mathcal{A}_H} + \frac{\eth \bar{\mathcal{C}_H}}{2 K}\right) J
  + \eth \mathcal C_H - \frac{\eth (J \bar{\mathcal{C}}_H)}{2 K} + \mathcal{C}_H^2\right]
  + \partial_r (J) \left[1 + \frac{1}{2} r^2 \partial_r W + \frac{3}{2} r W
  - \frac{1}{2} \left(- \bar{\eth} \bar{U} J + \frac{\eth U \bar{J}}{K^2}\right) K r
  - \frac{1}{4} \mathcal D_H K^2 r\right]
  + \partial_r (\bar{J}) J^2 r \left[\frac{1}{4} \mathcal D_H + \frac{1}{2 K}\eth U \bar{J} \right],
```
Simplified H hypersurface equation in Bondi-like coordinates: determines ∂_u J on Σ_u from all previously solved scalars.
Dependencies: eq:BondiHierarchyBeta, eq:auto-16, eq:script-b.

### eq:script-b (align containing the definitions 𝒜_H, ℬ_H, 𝒞_H, 𝒟_H, ℒ_H)
```latex
\mathcal{A}_H = \frac{1}{4} \eth (\eth (\bar{J})) - \frac{\eth \bar{\eth} (J \bar{J}) - 2 \bar J \eth \bar{\eth} J }{16 K^3}
  - \frac{\bar J \eth \bar{\eth} J - 3}{4 K}
  + \frac{1}{2} \eth (\bar{J}) \left(\mathcal C_H + \frac{\bar{\eth} (J \bar{J}) J}{4 K^3}
  - \frac{\bar{\eth} J (2 J \bar J + 1)}{4 K^3} \right),\\
\mathcal{B}_H = \frac{1}{2 r} + W + \frac{1}{2} r \partial_r W + 2 \partial_r \beta (r W + 1)
  - \frac{\eth U \bar{J} \partial_r (J \bar{J}) r}{4 K}
  + \frac{r U}{4} \left(\bar{\eth} J \partial_r \bar{J} + \bar \eth (\bar J \partial_r J)
  - \bar J \partial_r \bar \eth J - \frac{\bar{\eth} (J \bar{J}) \partial_r (J \bar{J})}{2 K^2}\right) \label{eq:script-b},\\
\mathcal{C}_H = \eth \beta - \frac{1}{2} Q,\\
\mathcal{D}_H = \eth \bar U - \bar \eth U,\\
\mathcal{L}_H = \frac{1}{2} \left(- \partial_r \bar{J} + \frac{\bar{J} \partial_r (J \bar{J})}{2 K^2}\right).
```
Source/linear-operator bundles entering the H hypersurface equation.
Dependencies: eq:KDef, eq:auto-6.

### eq:Hnumeric
```latex
(1 - \na y) \partial_{\na y} \na H + \na H + (1 - \na y)\na J (\mathcal{L}_{\na H} \na H
  + \bar{\mathcal{L}}_{\na H} \na{\bar{H}}) = S_{\na H}^P + (1 - \na y) S_{\na H}^R,
```
Compactified H equation, form eq:form3, governing the numerical evolution quantity H̆.
Dependencies: eq:HypersurfaceH, eq:numerical-jacobian-u, eq:numerical-h-conversion, eq:form3, eq:auto-30.

### eq:auto-30 (subequations group `eq:numerical-evolution-source-terms`)
```latex
\mathcal{L}_{\na H} = \frac{1}{2} \left(- \partial_{\na y} \na{\bar J} + \frac{\na{\bar J} \partial_{\na y}(\na J \na{\bar J} )}{2 \na K^2}\right),\\
S_{\na H}^P = - \frac{1}{2} \left(\eth (\na J \na{\bar U}) + \bar \eth(\na J \na U) + 2 J \eth \bar U\right)
  + 2\na W \na J - \na K \eth \na U,\\
S_{\na H}^R = \na J (\mathcal{B}_{\na H} + \bar{\mathcal{B}}_{\na H})
  + \frac{1}{2} (1-\na y) \left(\na W + \frac{(1-\na y)}{2 R} + 2 \frac{\partial_{\na u} R}{R}\right) \partial_{\na y}^2 \na J
  - \frac{1}{2} \left(\na U \bar \eth \partial_{\na y} \na J + \na {\bar U} \eth \partial_{\na y} \na J\right)
  + \na J^2 \partial_{\na y} \na{\bar J} \left[\frac{\mathcal{D}_{\na H}}{4} + \frac{\eth \na U \na{\bar J}}{2 \na K}\right]
  + \frac{\partial_{\na y} \na J}{2} \left[- \frac{1}{2}(\eth \bar{\na U} + \bar \eth \na U)
  - \na U \frac{\bar \eth R}{R} - \na{\bar U} \frac{\eth R}{R} + \na W + (1- \na y) \partial_{\na y} \na W
  - \na K \left(\frac{\eth \na U \na{\bar J}}{\na K} - \bar \eth \na{\bar U} \na J\right)
  - \frac{1}{2}(\na K^2 + 1) \mathcal{D}_{\na H} - \eth \na{\bar U} \partial_{\na y} \na J\right]
  + \frac{e^{2 \na \beta}}{2 R} \left[\frac{J}{2 K^3} + \frac{1}{2} \bar{\eth}\na J \na K \mathcal{C}_{\na H}
  + \frac{\bar{\mathcal{C}}_{\na H} \na J^2 \eth \na{\bar J}}{2 \na K}
  - \left(\mathcal{A}_{\na H} + \bar{\mathcal{A}}_{\na H}\right)\na J + \eth \mathcal C_{\na H}
  - \frac{\bar{\mathcal{C}}_{\na H} \eth \na J}{2\na K} + \mathcal{C}_{\na H}^2\right],
```
Source terms ℒ_H̆, S^P_H̆, S^R_H̆ of the compactified H equation (𝒜, 𝒞, 𝒟 obtained by breve substitution).
Dependencies: eq:script-b, eq:numerical-jacobian-u.

### eq:auto-31
```latex
\mathcal{B}_{\na H} = \frac{1}{2} \bigg[\frac{1}{2 R} + \partial_{\na y} \na W
  + \left(\na W + \frac{1 - \na y}{2 R} + 2 \frac{\partial_{\na u} R}{R}\right) \partial_{\na y} \na \beta
  - \frac{\eth \na U \na{\bar J} \partial_{\na y} (\na J \na{\bar J})}{2 \na K}
  + \frac{\na U}{2} \left(\bar \eth \na J \partial_{\na y} \na{\bar J} + \bar \eth(\na{\bar J} \partial_{\na y} \na J)
  - \na{\bar J} \bar \eth \partial_{\na y} \na J
  - \frac{\bar \eth (\na J \na{\bar J}) \partial_{\na y} (\na J \na{\bar J})}{2 \na K^2}\right)\bigg]
```
Numerical-coordinate bundle ℬ_H̆ (the breve version of ℬ_H including moving-worldtube terms).
Dependencies: eq:script-b, eq:numerical-jacobian-u.

## Section V — Regularity-preserving CCE (sec:regularity_preservation)

### eq:QWabstract
```latex
\partial_l \left(\frac{F_2}{l^2}\right) = \frac{S_2^P}{l^3} + \frac{S^R_2}{l^2}
```
Form-2 equation re-expressed in l = 1/r, exposing the danger of l^{-1} right-hand-side scaling → l² ln(l) solutions.
Dependencies: eq:form2.

### eq:QWregularity
```latex
\frac{1}{2} \partial_{l}^2 S_2^P |_{l = 0} = - \partial_l S^R_2|_{l = 0}
```
Regularity (no-logarithm) condition on the sources of the Q and W equations.
Dependencies: eq:QWabstract.

### eq:Habstract
```latex
\partial_l \left(\frac{F_3}{l}\right) + \frac{1}{l}( L_{F_3} F_3 + L_{\bar F_3} \bar F_3)
  = \frac{S_3^P}{l^2} + \frac{S_3^R}{l}.
```
Form-3 (H) equation in l, with the analogous logarithm danger at linear order.
Dependencies: eq:form3.

### eq:Hregularity
```latex
\partial_l S_3^P |_{l = 0} = - S_3^R |_{l = 0}.
```
Regularity condition for the H equation.
Dependencies: eq:Habstract.

### eq:auto-32 (subequations group `eq:BetaNearScri`)
```latex
\ifc \beta^{(1)} = 0,\\
\ifc \beta^{(2)} = \frac{\ifc{\bar J}^{(1)} \ifc J^{(1)}}{16},
```
Asymptotic expansion of β̂ in the partially flat gauge (Ĵ^(0) = Û^(0) = 0); β̂^(0) fixed by boundary data; no regularity condition arises.
Dependencies: eq:HypersurfaceBeta.

### eq:SubleadingQRegularity (align containing Q̂^(0) and Q̂^(1))
```latex
\ifc Q^{(0)} = - 2 \ifc \eth \ifc \beta^{(0)},\\
\ifc Q^{(1)} = \ifc{\bar \eth} \ifc J^{(1)}. \label{eq:SubleadingQRegularity}
```
Leading and subleading parts of Q̂ from the Q hypersurface equation in the partially flat gauge.
Dependencies: eq:HypersurfaceQ, eq:auto-32, eq:QWregularity.

### eq:auto-33
```latex
\ifc J^{(2)} = 0,
```
Regularity condition from the O(l) part of the Q equation: the l² coefficient of Ĵ must vanish (key initial-data condition).
Dependencies: eq:SubleadingQRegularity, eq:QWregularity.

### eq:auto-34
```latex
\ifc U^{(1)} = 2 e^{2\ifc \beta^{(0)}} \ifc \eth \ifc \beta^{(0)},\\
\ifc U^{(2)} = -\frac{e^{2 \ifc \beta^{(0)}} (\ifc{\bar \eth} \ifc J^{(1)} + 2 \ifc{\bar \eth} \ifc \beta^{(0)} \ifc J^{(1)})}{2}.
```
Perturbative expansion of Û from its hypersurface equation; regular automatically.
Dependencies: eq:HypersurfaceU, eq:SubleadingQRegularity, eq:auto-32.

### eq:auto-35 (subequations group `eq:WNearScri`)
```latex
\ifc W^{(0)} = 0,\\
\ifc W^{(1)} = e^{2 \ifc \beta^{(0)}} - 1 + 2 e^{2 \ifc \beta^{(0)}} \ifc \eth \ifc{\bar \eth} \ifc \beta^{(0)}
  + 4 e^{2 \ifc \beta^{(0)}} \ifc \eth \ifc \beta^{(0)} \ifc{\bar \eth} \ifc \beta^{(0)}.
```
Order-by-order constraints from the Ŵ hypersurface equation; regularity directly satisfied.
Dependencies: eq:HypersurfaceW, eq:auto-34, eq:auto-32, eq:QWregularity.

### eq:auto-36 (subequations group `eq:EvolutionNearScri`)
```latex
\partial_{\ifc u}\ifc J^{(0)} = 0,\\
\partial_{\ifc u} \ifc J^{(2)} = 0,
```
Evolution stability of the regularity conditions: the H expansion preserves Ĵ^(0) = 0 and Ĵ^(2) = 0 in time.
Dependencies: eq:HypersurfaceH, eq:Hregularity, eq:auto-32, eq:auto-33, eq:auto-34, eq:auto-35.

### eq:UHat
```latex
\mathcal{U} = \frac{1}{2 \hat{\omega}^2} \left(\ifc{\bar b} U - \ifc a \bar{U} \right)
  - \frac{\ifc l e^{2 \ifc \beta}}{\ifc \omega} \left(\ifc \eth \ifc \omega \ifc K - \ifc{\bar \eth} \ifc \omega \ifc J\right).
```
Auxiliary shift 𝒰 with ∂_r𝒰 = ∂_r Û, enabling the U integration before U^(0) is known on the current hypersurface.
Dependencies: eq:auto-19, eq:ifc_omega, eq:Unumeric, eq:ifcBeta, eq:ifcW.

### eq:auto-37
```latex
U^{(0)} = \frac{1}{2 \ifc \omega^2} \left(\bar b\, \mathcal{U}^{(0)} - a\, \mathcal{\bar{U}}^{(0)}\right).
```
Inversion recovering the Bondi-like asymptotic shift U^(0) from the scri+ value of 𝒰.
Dependencies: eq:UHat.

### eq:auto-38
```latex
\partial_u \ifc{x}^{\ifc{A}} = - U^{(0) B} \partial_B \ifc{x}^{\ifc{A}},
```
Completed Jacobian evolution for the partially flat angular coordinates within the evolution loop.
Dependencies: eq:XIfcEom, eq:auto-37.

## Section VI — Newman-Penrose Weyl scalars (sec:weyl_scalars)

### eq:auto-39 (subequations group `eq:NP_tetrads`)
```latex
m^\mu = -\frac{1}{\sqrt{2} r}\left(\sqrt{\frac{K + 1}{2}}q^\mu
  - \sqrt{\frac{1}{2(1 + K)}} J \bar{q}^\mu\right),\\
n^\mu = \sqrt{2} e^{-2 \beta} \left(\delta^\mu{}_u - \frac{V}{2 r} \delta^\mu{}_r
  + \frac{1}{2} \bar{U} q^\mu + \frac{1}{2} U \bar{q}^\mu\right),\\
l^\mu = \frac{1}{\sqrt{2}} \delta^\mu{}_r.
```
NP null tetrad adapted to Bondi-like coordinates, asymptotically matching common numerical-relativity tetrads.
Dependencies: eq:BondiLikeMetric, eq:KDef.

### eq:auto-40 (subequations group `eq:NP_orthonormal`)
```latex
-l^\mu n_\mu = m^\mu \bar{m}_\mu = 1\\
l^\mu l_\mu = l^\mu m_\mu = n^\mu n_\mu = n^\mu m_\mu = m^\mu m_\mu = 0.
```
Standard normalization and orthogonality of the tetrad.
Dependencies: eq:auto-39.

### eq:auto-41
```latex
g_{\mu \nu} = - 2 l_{(\mu} n_{\nu)} + 2 m_{(\mu} \bar{m}_{\nu)}.
```
Metric reconstruction from the tetrad.
Dependencies: eq:auto-40.

### eq:NPCovariantDerivatives
```latex
D = l^\mu \nabla_\mu, \qquad \Delta = n^\mu \nabla_\mu, \qquad \delta = m^\mu \nabla_\mu.
```
NP directional derivatives along the tetrad legs.
Dependencies: eq:auto-39.

### eq:auto-42
```latex
\kappa = - m^\mu l^\nu \nabla_\nu l_\nu,\\
\rho = -m^\mu \bar{m}^\nu \nabla_\nu l_\mu,\\
\sigma = - m^\mu m^\nu \nabla_\nu l_\mu,\\
\tau = -m^\mu n^\nu \nabla_\nu l_\mu,\\
\nu = \bar{m}^\mu n^\nu \nabla_\nu n_\mu,\\
\mu = \bar{m}^\mu m^\nu \nabla_\nu n_\mu,\\
\lambda = \bar{m}^\mu \bar{m}^\nu \nabla_\nu n_\mu,\\
\pi = \bar{m}^\mu l^\nu \nabla_\nu n_\mu,\\
\epsilon = \frac{1}{2} (\bar{m}^\mu l^\nu \nabla_\nu m_\mu - n^\mu l^\nu \nabla_\nu l_\mu),\\
\beta_{\text{NP}} = \frac{1}{2} (\bar{m}^\mu m^\nu \nabla_\nu m_\mu - n^\mu m^\nu \nabla_\nu l_\mu),\\
\gamma = \frac{1}{2} (\bar{m}^\mu n^\nu \nabla_\nu m_\mu - n^\mu n^\nu \nabla_\nu l_\mu),\\
\alpha = \frac{1}{2} (\bar{m}^\mu \bar{m}^\nu \nabla_\nu m_\mu - n^\mu \bar{m}^\nu \nabla_\nu l_\mu).
```
Definitions of the twelve NP spin coefficients (β_NP disambiguated from Bondi β).
Dependencies: eq:auto-39, eq:NPCovariantDerivatives.

### eq:auto-43 (subequations group `eq:spin_coefficients`)
```latex
\kappa = 0,\qquad \rho = -\frac{1}{\sqrt{2} r},\\
\sigma = -\frac{(1 + K) \partial_r J}{4 \sqrt{2} K} + \frac{J^2 \partial_r \bar{J}}{4 \sqrt{2} K(1 + K)},\\
\tau = -\frac{(2 \bar{\eth} \beta - \bar Q) J }{4 r \sqrt{1 + K}} + \frac{(2\eth (\beta) - Q) \sqrt{1 + K}}{4 r},\\
\nu = \frac{\eth W \bar{J}}{2 e^{2 \beta} \sqrt{1 + K}} - \frac{\bar{\eth} W \sqrt{1 + K}}{2 e^{2 \beta}},\\
\mu = \frac{e^{-2 \beta}(\eth \bar U + \bar \eth U)}{2 \sqrt{2}} - \frac{e^{-2\beta} (r^2 W + r)}{\sqrt{2} r^2},\\
\lambda = \frac{e^{-2 \beta}}{\sqrt{2}}\bigg[\frac{1}{2} \bar J \left(\bar \eth U - \eth \bar U\right)
  - \frac{\bar J^2 \eth U}{2(1 + K)} + \frac{\bar \eth \bar U (1 + K)}{2}
  + \frac{r^2 W + r}{r} \left(\frac{\bar J^2 \partial_r J}{1 + K} - \partial_r \bar J (1 + K)\right)
  + \frac{U}{4 K} \left(\bar \eth \bar J (1+K) - \frac{\bar J^2 \bar \eth J}{1 + K}\right)
  + \frac{\bar U}{4 K} \left(\eth \bar J (1 + K) - \frac{\bar J^2 \eth J}{1 + K}\right)
  - \frac{\bar J^2 \partial_u J}{2 K (1 + K)} + \frac{(1 + K)\partial_u \bar J}{2K}\bigg],\\
\pi = \frac{(2 \eth \beta + Q)\bar J }{4 r \sqrt{1 + K}} - \frac{(2\bar \eth (\beta) +\bar Q) \sqrt{1 + K}}{4 r},\\
\epsilon = \frac{\partial_r (\beta)}{\sqrt 2} + \frac{J \partial_r (\bar{J}) - \bar{J} \partial_r (J)}{8 \sqrt 2 (1 + K)},\\
\beta_{\text{NP}} = \frac{1}{4 r \sqrt{1 + K}} \bigg[\bar \eth \beta J - \eth \beta (1 + K)
  - \frac{\eth (J \bar J)}{4 K} - \frac{J(\bar \eth (J \bar J) - \bar J \bar \eth J)}{4 K (1 + K)}
  + \frac{\bar \eth J (1 + 3 K)}{4 K} + \frac{1}{2} J \bar Q - \frac{1}{2} Q (1 + K)
  - \frac{1}{2} J \bar{\Theta} - \frac{1}{2}(1 + K) \Theta \bigg],\\
\alpha = \frac{1}{4 r \sqrt{1 + K}} \bigg[\eth \beta \bar J - \bar \eth \beta (1 + K)
  + \frac{\bar \eth (J \bar J)}{4 K} + \frac{\bar J(\eth (J \bar J) - J \eth \bar J)}{4 K (1 + K)}
  - \frac{\eth J (1 + 3 K)}{4 K} + \frac{1}{2} \bar J Q - \frac{1}{2} \bar Q (1 + K)
  + \frac{1}{2} \bar J \Theta + \frac{1}{2}(1 + K) \bar \Theta \bigg],\\
\gamma = \frac{e^{-2 \beta}}{\sqrt{2}}\bigg[\frac{1}{4} \left(J \bar \eth \bar U - \bar J \eth U\right)
  + \frac{K}{4}\left(\bar \eth U - \eth \bar U\right)
  + \frac{r^2 W + r}{8r(1 + K)} \left(\bar J \partial_r J - J \partial_r \bar J\right)
  + \frac{1}{2}(r \partial_r W + W)
  + \frac{U\left(\bar \eth(J \bar J) - 2 \bar J \bar \eth J\right)}{8 (1 +K)}
  - \frac{\bar U \left(\eth (J \bar J) - 2 J \eth \bar J\right)}{8(1 + K)}
  - \frac{\bar J \partial_u J}{4 (1 + K)} + \frac{J\partial_u \bar J}{4 (1 + K)}
  + \frac{1}{4} \bar U \Theta - \frac{1}{4} U \bar \Theta\bigg].
```
Closed-form spin coefficients in Bondi-like spin-weighted scalars (Θ = unit-sphere dyad connection).
Dependencies: eq:auto-42, eq:auto-39, eq:QDef, WDef, eq:KDef.

### eq:auto-44 (subequations group `eq:weyl_scalar_identities`)
```latex
\Psi_0 = D \sigma - \delta \kappa - (\rho + \bar \rho) \sigma - (3 \epsilon -\bar \epsilon) \sigma
  + (\tau - \bar \pi + \bar \alpha + 3 \beta_{\text{NP}}) \kappa,\\
\Psi_1 = D \beta_{\text{NP}} - \delta \epsilon - (\alpha + \pi) \sigma - (\bar \rho - \bar \epsilon) \beta_{\text{NP}}
  + (\mu + \gamma) \kappa + (\bar \alpha - \bar \pi) \epsilon,\\
\Psi_2 = D\mu - \delta \pi - (\bar \rho - \epsilon - \bar \epsilon) \mu - \sigma \lambda
  + (\bar \alpha - \beta_{\text{NP}} - \bar \pi) \pi + \nu \kappa,\\
\Psi_3 = D \nu - \Delta \pi - (\pi + \bar \tau) \mu - (\bar \pi + \tau) \lambda
  - (\gamma - \bar \gamma) \pi + (3 \epsilon + \bar \epsilon) \nu,\\
\Psi_4 = -\Delta \lambda + \bar{\delta}\nu - \lambda(\mu + \bar{\mu})
  - (3 \gamma - \bar{\gamma}) \lambda + (3 \alpha + \bar{\beta}_{\text{NP}} + \pi - \bar{\tau}) \nu,
```
Vacuum (R_μν = 0) NP identities expressing the Weyl scalars through spin coefficients.
Dependencies: eq:auto-42, eq:NPCovariantDerivatives.

### eq:auto-45
```latex
\Delta^{\text{SW}} = \sqrt{2} e^{-2 \beta} \left(\left(\delta^\mu{}_u - \frac{V}{2 r} \delta^\mu{}_r\right)\nabla_\mu
  + \frac{1}{2} \bar U \eth + \frac{1}{2} U \bar \eth\right),\\
\delta^{\text{SW}} = - \frac{1}{\sqrt{2} r} \left(\sqrt{\frac{K + 1}{2}}\eth - \frac{1}{2(1 + K)} J \bar\eth\right).
```
Spin-weighted generalizations of the NP scalar derivatives (coordinate-connection pieces absorbed into ð).
Dependencies: eq:NPCovariantDerivatives, eq:auto-39.

### eq:auto-46
```latex
\{\beta_{\text{NP}}, \alpha, \gamma, \Delta, \delta\}\rightarrow \{\beta_{\text{NP}}^{\text{SW}},
  \alpha^{\text{SW}}, \gamma^{\text{SW}}, \Delta^{\text{SW}}, \delta^{\text{SW}}\}.
```
Replacement rule: the Weyl identities are unchanged when Θ-dependent quantities are swapped for spin-weighted versions.
Dependencies: eq:auto-45, eq:auto-43.

### eq:auto-47 (subequations group `eq:PartialExpandPsi24`)
```latex
\Psi_2 = \partial_r \mu + \frac{1}{2 r} \left(\frac{J \bar \eth \pi}{\sqrt{1+K}} - \sqrt{1 + K} \eth \pi\right)
  + (\bar \epsilon + \epsilon - \bar \rho )\mu + (\bar \alpha^{\text{SW}} - \bar \pi - \beta_{\text{NP}}^{\text{SW}}) \pi - \lambda \sigma,\\
\Psi_3 = \frac{1}{\sqrt{2}} \partial_r \nu - \sqrt{2} e^{-2 \beta} \partial_u \pi
  + \frac{e^{-2 \beta}}{\sqrt 2}(r W + 1) \left(\partial_r \pi - U \bar \eth \pi - \bar U \eth \pi\right)
  - (\pi + \bar \tau) \mu - (\bar \pi + \tau) \lambda - (\gamma^{\text{SW}} - \bar{\gamma}^{\text{SW}}) \pi + (3 \epsilon + \bar{\epsilon}) \nu,\\
\Psi_4 = \frac{1}{2 r} \left(\sqrt{1 + K}\bar \eth \nu - \frac{\bar J \eth \nu}{\sqrt{1 + k}}\right)
  - e^{-2 \beta} \partial_u \lambda + \frac{e^{-2 \beta}}{2}(r W + 1) \left(\partial_r \lambda - U \bar \eth \lambda - \bar U \eth \lambda\right)
  - (\mu + \bar \mu + 3 \gamma^{\text{SW}} - \bar \gamma^{\text{SW}}) \lambda
  + (3 \alpha^{\text{SW}} + \bar \beta_{\text{NP}}^{\text{SW}} + \pi - \bar \tau)\nu
```
Bulk Ψ₂, Ψ₃, Ψ₄ after the spin-weighted replacements (full expansions in the companion notebook).
Dependencies: eq:auto-44, eq:auto-46, eq:auto-43.

### eq:auto-48 (subequations group `eq:FullExpandPsi01`)
```latex
\Psi_0 = \left(\frac{r \partial_r \beta - 1}{4 K r}\right)\left((1 + K) \partial_r J - \frac{J^2 \partial_r \bar J}{(1 + K)}\right)
  + \frac{J (1 + K^2) \partial_r J \partial_r \bar J}{8 K^3}
  + \frac{1}{8 K}\left(\frac{J^2 \partial_r^2 \bar J}{1 + K} - (1 + K) \partial_r^2 J \right)
  + \frac{-J \bar J^2 \partial_r J^2 - J^3 \partial_r \bar J^2}{16 K^3},\\
\Psi_1 = \frac{1}{4 \sqrt{2(K+1)} r} \left(- J \partial_r(2 \bar \eth \beta - \bar Q)
  + \bigg(\frac{J^2 \partial_r \bar J - (1 + K) \partial_r J}{4 K} + \frac{J}{r}\right) (2 \bar \eth \beta + \bar Q)
  - (1+K)\partial_r \left(2 \eth \beta - Q\right)
  - \left(\frac{1 + K}{r} + \frac{(1 + K)\bar J \partial_r J}{4 K} - \frac{J^2 \bar J \partial_r \bar{J}}{4K(1+K)}\right)\left(2 \eth \beta + Q\right)\bigg).
```
Concise closed bulk expressions for Ψ₀ and Ψ₁ (cancellations make these compact).
Dependencies: eq:auto-44, eq:auto-43.

### eq:auto-49
```latex
\Psi_n \sim r^{n-5}.
```
Peeling theorem falloff of the Weyl scalars.
Dependencies: none (literature).

### eq:auto-50 (subequations group `eq:if_weyl_scalars`)
```latex
\lim_{\hat{r} \rightarrow \infty} \hat{r}^5 \Psi_0^{\text{PF}}
  = \frac{3}{2}\left(\frac{1}{4}\ifc{\bar J}^{(1)} \ifc J^{(1)} {}^2 -\ifc J^{(3)}\right),\\
\lim_{\hat{r} \rightarrow \infty} \hat{r}^4 \Psi_1^{\text{PF}}
  = \frac{1}{8}\left(-12 \ifc \eth \ifc \beta^{(2)} +\ifc J^{(1)} \ifc{\bar{Q}}^{(1)} + 2\ifc Q^{(2)}\right),\\
\lim_{\hat{r} \rightarrow \infty} \hat{r}^3 \Psi_2^{\text{PF}}
  = -\frac{e^{-2 \ifc \beta^{(0)}}}{4}\left(e^{2 \ifc\beta^{(0)}} \eth \ifc{\bar Q}^{(1)} + \ifc \eth \ifc{\bar U}^{(2)}
  + \ifc{\bar \eth} \ifc U^{(2)} +\ifc J^{(1)} \ifc{\bar \eth} \ifc{\bar U}^{(1)}
  + \ifc J^{(1)} \ifc{\bar H}^{(1)} -2 \ifc W^{(2)}\right),\\
\lim_{\hat{r} \rightarrow \infty} \hat{r}^2 \Psi_3^{\text{PF}}
  = 2 \ifc{\bar \eth} \ifc \beta^{(0)} + 4 \ifc{\bar \eth} \ifc \beta^{(0)} \ifc{\bar \eth} \ifc \eth \ifc \beta^{(0)}
  + \ifc{\bar \eth} \ifc{\bar \eth} \ifc \eth \ifc \beta^{(0)}
  + \frac{e^{-2 \ifc \beta^{(0)}}}{2} \ifc \eth \ifc{\bar H}^{(1)}
  - e^{-2 \ifc \beta^{(0)}} \ifc \eth \ifc \beta^{(0)} \ifc{\bar H}^{(1)},\\
\lim_{\hat{r} \rightarrow \infty} \hat{r} \Psi_4^{\text{PF}}
  = -e^{-2 \ifc \beta^{(0)}} \partial_{\ifc u} \left[e^{-2 \ifc \beta^{(0)}} \left(\ifc{\bar{\eth}}\ifc{\bar{U}}^{(1)}
  + \ifc{\bar{H}}^{(1)}\right)\right].
```
Leading asymptotic Weyl scalars in the partially flat gauge: simple formulas in the expansion coefficients of the evolved scalars.
Dependencies: eq:auto-47, eq:auto-48, eq:auto-49, eq:auto-32, eq:auto-33, eq:auto-34, eq:auto-35.

### eq:auto-51 (subequations group `eq:bondi_weyl_scalars`)
```latex
\lim_{\bs r \rightarrow \infty} \bs r^5 \Psi_0^{\text{Bondi}}
  = \frac{3}{2}\left(\frac{1}{4}\bs{\bar J}^{(1)} \bs J^{(1)} {}^2 -\bs J^{(3)}\right)
  = \Psi_0^{\rm PF(5)} + 2 \ifc \eth \bs u \Psi_1^{\rm PF(4)} + \frac{3}{4} \left(\ifc \eth \bs u\right)^2 \Psi_2^{\rm PF(3)}
  + \frac{1}{2} \left(\ifc \eth \bs u\right)^3 \Psi_3^{\rm PF(2)} + \frac{1}{16} \left(\ifc \eth \bs u\right)^4 \Psi_4^{\rm PF(1)},\\
\lim_{\bs r \rightarrow \infty} \bs r^4 \Psi_1^{\text{Bondi}}
  = \frac{1}{8}\left(-12 \bs \eth \bs \beta^{(2)} +\bs J^{(1)} \bs{\bar{Q}}^{(1)} + 2\bs Q^{(2)}\right)
  = \Psi_1^{\rm PF(4)} + \frac{3}{2} \ifc \eth \bs u \Psi_2^{\rm PF(3)} + \frac{3}{4} \left(\ifc \eth \bs u\right)^2 \Psi_3^{\rm PF(3)}
  + \frac{1}{8} \left(\ifc \eth \bs u\right)^3 \Psi_4^{\rm PF(1)},\\
\lim_{\bs r \rightarrow \infty} \bs r^3 \Psi_2^{\text{Bondi}}
  = \frac{\bs W^{(2)}}{2} - \frac{1}{8} \bs \eth \bs \eth \bs{\bar J}^{(1)} + \frac{1}{8}\bs{\bar \eth} \bs{\bar \eth} \bs J^{(1)}
  - \frac{1}{4} \bs J^{(1)} \bs{\bar H}^{(1)}
  = \Psi_2^{\rm PF(3)} + \ifc \eth \bs u \Psi_3^{\rm PF(2)} + \frac{1}{4} \left(\ifc \eth \bs u\right)^2 \Psi_4^{\rm PF(1)},\\
\lim_{\bs r \rightarrow \infty} \bs r^2 \Psi_3^{\text{Bondi}}
  = \frac{1}{2} \bs \eth \bs{\bar H}^{(1)}
  = \Psi_3^{\rm PF (2)} + \frac{1}{2}\ifc \eth \bs u^{(0)} \Psi_4^{\rm PF (1)},\\
\lim_{\bs r \rightarrow \infty} \bs r \Psi_4^{\text{Bondi}}
  = -\partial_{\bs u} \bs{\bar{H}}^{(1)} = \Psi_4^{\rm PF (1)}.
```
Leading Weyl scalars in true Bondi-Sachs gauge and their relation to the partially flat scalars (Boyle 2016 transformations adapted to these conventions; all PF scalars evaluated at û(ů)).
Dependencies: eq:auto-50, eq:auto-23.

## Appendix B — Additional transformations to partially flat coordinates (sec:IF_extras)

### eq:ifcQ
```latex
\ifc{Q} = \ifc r^2 e^{-2 \ifc \beta} (\ifc K \partial_{\ifc r} \ifc U + \ifc J \partial_{\ifc r} \ifc{\bar U}),
```
Partially flat Q̂ from ∂_r̂Û via the Q definition in hatted coordinates (worldtube boundary computation).
Dependencies: eq:QDef, eq:auto-52, eq:ifcBeta.

### eq:auto-52
```latex
\partial_{\ifc r} \ifc U = \frac{1}{2 \ifc \omega^3}\left(\ifc{\bar b} \partial_r U - \ifc a \partial_r \bar U\right)
  + \frac{e^{2\ifc \beta}}{4 \ifc \omega} \left(\ifc J \ifc{\bar{\eth}} \ifc \omega - \ifc K \ifc \eth \ifc \omega\right)
    \left(\partial_{\ifc r} \ifc{\bar{J}} \partial_{\ifc r} \ifc J - \frac{\partial_{\ifc r}(\ifc J \ifc{\bar{J}})^2}{4 \ifc K^2}\right)
  + \frac{e^{2 \ifc \beta}}{\ifc \omega \ifc r} \left[\ifc{\bar{\eth}} \ifc \omega \left(\partial_{\ifc r} \ifc J-\frac{\ifc J}{\ifc r}\right)
  + \ifc{\bar{\eth}} \ifc \omega \left(\frac{\ifc K}{\ifc r} - \frac{\ifc J \partial_{\ifc r} \ifc{\bar J} + \ifc{\bar J} \partial_{\ifc r} \ifc J}{2 \ifc K}\right) \right].
```
Explicit ∂_r̂Û in terms of un-hatted derivatives and Jacobian factors, avoiding extra input-metric derivatives.
Dependencies: eq:auto-19, eq:ifc_omega, eq:ifcBeta, eq:ifcW.

### eq:ifcH
```latex
\ifc H = \frac{\partial_{\ifc u} \ifc \omega - \tfrac{1}{2} \left(\mathcal{U}^{(0)} \bar{\ifc \eth}\ifc \omega
  + \bar{\mathcal{U}}^{(0)} \ifc \eth \ifc \omega \right) }{\ifc \omega} \left(2 \ifc J - \ifc r \partial_{\ifc r} \ifc J\right)
  - \ifc J\ifc{\bar \eth} \mathcal{U}^{(0)} + \ifc K \ifc \eth \bar{\mathcal{U}}^{(0)}
  + \frac{1}{4 \ifc \omega} \left(\ifc{\bar b}^2 H + \ifc a^2 \bar H + \ifc{\bar b} a \frac{H \bar J + J \bar H}{K}\right)
  + \frac{1}{2} \left(\mathcal{U}^{(0)} \ifc{\bar \eth} \ifc J + \ifc{\eth}(\bar{\mathcal{U}}^{(0)} \ifc J) - \ifc J \ifc \eth \bar{\mathcal{U}}^{(0)}\right)
```
Partially flat Ĥ transformation for worldtube boundary data (critical for the regularity-preserving scheme).
Dependencies: eq:auto-16, eq:auto-19, eq:auto-20, eq:ifc_omega, eq:ifcW.

## Appendix C — Perturbative expansion of the transformations near scri+ (app:perturbative_transformations)

### eq:auto-53
```latex
2 \mathring u^{(1)} = -\eth \mathring u^{(0)} \bar \eth \mathring u^{(0)},\\
-4 \mathring u^{(2)} - 2 e^{- 2 \hat \beta} \partial_{\hat u} \mathring u^{(1)} \mathring u^{(1)}
  = 2 \partial_{\hat A} \mathring u^{(0)} \mathring u^{(1)} \hat U^{(1) \hat A} e^{-2 \hat \beta}
  - \left(\bar \eth \mathring u^{(0)}\right)^2 J^{(1)} - \left(\eth \mathring u^{(0)}\right)^2 \bar J^{(1)}.
```
First two orders ů^(1), ů^(2) of the perturbative solution of eq:uRequation.
Dependencies: eq:uRequation.

### eq:auto-54
```latex
-\bs x^{(1) \bs A} = \ifc q^{\ifc A \ifc B} \partial_{\ifc B} \bs u^{(0)},\\
-\bs x^{(2) \bs A} = \partial_{\ifc u} \bs u^{(1)} e^{-2 \ifc \beta} \bs x^{(1) \bs A}
  + \bs u^{(1)} \partial_{\ifc u} \bs x^{(1) \bs A} e^{-2 \ifc \beta}
  + \delta^{\bs A}{}_{\ifc A} \ifc U^{(1) \ifc A} \bs u^{(1)} e^{-2 \ifc \beta}
  + \bs x^{(1) \bs A} \partial_{\ifc A} \bs u^{(0)} \ifc U^{(1) \ifc A} e^{-2 \ifc \beta}
  + \partial_{\ifc B} \bs u^{(1)} \delta_{\ifc A}{}^{\bs A} \ifc q^{\ifc A \ifc B}
  + \partial_{\ifc B} \bs u^{(0)} \partial_{\ifc A} \bs x^{(1) \bs A} q^{\ifc A \ifc B}
  + \partial_{\ifc B} \bs u^{(0)} \delta{}_{\ifc A}{}^{\bs A} h^{(1) \ifc A \ifc B}.
```
First two orders of the subleading angular transformation x̊^(1)Å, x̊^(2)Å from eq:bsxA.
Dependencies: eq:bsxA, eq:auto-53.

### eq:auto-55
```latex
\bs \omega^{(1)} = \frac{1}{2} \ifc \eth \ifc{\bar \eth} \bs u^{(0)}.
```
Leading correction to the final conformal factor ω̊ = 1 + l̂ ω̊^(1) + O(l̂²).
Dependencies: eq:Step5RadialCoordinate, eq:auto-53, eq:auto-54.
