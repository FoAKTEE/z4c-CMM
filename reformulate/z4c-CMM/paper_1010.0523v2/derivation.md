# derivation.md — arXiv:1010.0523v2 (Ruiz, Hilditch, Bernuzzi)

Constraint preserving boundary conditions for the Z4c formulation of general relativity.

Every equation environment of `ref-paper/arxiv-1010.0523v2/src/Z4.tex` in document order.
Labeled equations keep their tex `\label{...}` verbatim as the node id; the 36 unlabeled
equation environments receive sequential ids `eq:auto-1` … `eq:auto-36` in document order.
All content is a literature-grounded transcription (evidence type: citation/transcription),
NOT an independently verified derivation. Equation environments in tex: 94; node ids: 116
(80 tex-born labels + 36 auto). Two `displaymath` blocks (the Laplace–Fourier form of the
Theta equation carrying the non-environment label `eq:theta-eq-LF`, and the gauge boundary
matrix `A`, plus the `A^{L+1}` widetext displaymath and the determinant condition) are not
equation environments and are described inline where relevant.

Conventions: see `convention.md`. `\p` = `\partial`, `\mbeta` = `\mathring\beta`,
`\hateq` = equality on the boundary surface, `\Lie` = Lie derivative.

---

## Section II — The Z4c formulation

### eq:Z4_ADM_1
```latex
\p_t\gamma_{ij} = -2\,\alpha\,K_{ij} + \Lie_\beta\gamma_{ij}
```
(Same `eqnarray` contains the unlabeled second row
`\p_t K_{ij} = -D_iD_j\alpha + \alpha[R_{ij} - 2K_{ik}K^k{}_j + K_{ij}K + 2\p_{(i}Z_{j)}]
+ 4\pi\alpha[\gamma_{ij}(S-\rho) - 2S_{ij}] + \Lie_\beta K_{ij}`.)
- Meaning: Z4 evolution equations for the 3-metric and extrinsic curvature; ADM equations augmented by the constraint field Z_i. ADM is recovered when Theta and Z_i vanish.
- Depends on: none (root; Z4 system imported from Bernuzzi–Hilditch 2009 [Bernuzzi:2009ex]).

### eq:Theta_dot
```latex
\p_t\Theta = \alpha\left[\tfrac{1}{2}H + \p_kZ^k\right] + \beta^i\,\Theta_{,i}
```
- Meaning: evolution of the scalar constraint field Theta, sourced by the Hamiltonian constraint and the divergence of Z.
- Depends on: eq:auto-1 (definition of H).

### eq:Z4_ADM_2
```latex
\p_tZ_{i} = \alpha\,M_i + \alpha\,\Theta_{,i} + \beta^j\,Z_{i,j}
```
- Meaning: evolution of the vector constraint field Z_i, sourced by the momentum constraint and gradients of Theta.
- Depends on: eq:auto-1 (definition of M^i).

### eq:auto-1
```latex
H = R - K_{ij}K^{ij} + K^2 - 16\pi\rho = 0,\qquad
M^i = D_j(K^{ij} - \gamma^{ij}K) - 8\pi S^i = 0
```
- Meaning: definitions of the Hamiltonian and momentum constraints of GR in ADM variables.
- Depends on: none (root definition).

### eq:ham-Z4
```latex
\p_0 H \simeq -2\,\p_iM^i
```
- Meaning: principal-part propagation of the Hamiltonian constraint under the Z4 evolution.
- Depends on: eq:Z4_ADM_1, eq:Theta_dot, eq:Z4_ADM_2, eq:auto-1, eq:auto-2.

### eq:mom-Z4
```latex
\p_0 M_i \simeq -\tfrac{1}{2}\p_iH + \p^j\p_jZ_i - \p_i\p^jZ_j
```
- Meaning: principal-part propagation of the momentum constraint under the Z4 evolution.
- Depends on: eq:Z4_ADM_1, eq:Theta_dot, eq:Z4_ADM_2, eq:auto-1, eq:auto-2.

### eq:auto-2
```latex
\p_0 = \frac{1}{\alpha}\left(\p_t - \beta^i\p_i\right)
```
- Meaning: definition of the derivative along the normal to the time slices (normal time derivative).
- Depends on: none (root definition).

### eq:auto-3
```latex
\tilde{\gamma}_{ij} = \gamma^{-\frac{1}{3}}\gamma_{ij},\quad
\hat{K} = \gamma^{ij}K_{ij} - 2\Theta,\quad
\chi = \gamma^{-\frac{1}{3}},\quad
\tilde{A}_{ij} = \gamma^{-\frac{1}{3}}(K_{ij} - \tfrac{1}{3}\gamma_{ij}K)
```
- Meaning: definitions of the Z4c conformal variables: conformal metric, modified trace of K, conformal factor, and trace-free conformal extrinsic curvature.
- Depends on: eq:Z4_ADM_1, eq:Theta_dot (variables gamma_ij, K_ij, Theta).

### eq:auto-4
```latex
\tilde{\Gamma}^{i} = 2\,\tilde{\gamma}^{ij}Z_j
+ \tilde{\gamma}^{ij}\tilde{\gamma}^{kl}\tilde{\gamma}_{jk,l}
```
- Meaning: definition of the Z4c conformal connection variable, absorbing the constraint Z_i.
- Depends on: eq:auto-3, eq:Z4_ADM_2 (variable Z_i).

### eq:auto-5
```latex
\tilde{\Gamma}_{\textrm{d}}{}^i
= \tilde{\gamma}^{ij}\tilde{\gamma}^{kl}\tilde{\gamma}_{jk,l}
= \gamma^{\frac{1}{3}}\gamma^{ij}\gamma^{kl}\left(\gamma_{jk,l}
- \tfrac{1}{3}\gamma_{kl,j}\right)
```
- Meaning: shorthand for the purely metric-derivative part of the conformal connection.
- Depends on: eq:auto-3.

### eq:Z4_decomp_first
```latex
\p_t \chi = \frac{2}{3}\chi[\alpha(\hat{K}+2\Theta) - D_i\beta^i]
```
(Same `align` contains the unlabeled rows
`\p_t\tilde{\gamma}_{ij} = -2\alpha\tilde{A}_{ij} + \beta^k\tilde{\gamma}_{ij,k}
+ 2\tilde{\gamma}_{(i|k}\beta^k_{,|j)} - \tfrac{2}{3}\tilde{\gamma}_{ij}\beta^k{}_{,k}` and
`\p_t\hat{K} = -D^iD_i\alpha + \alpha[\tilde{A}_{ij}\tilde{A}^{ij}
+ \tfrac{1}{3}(\hat{K}+2\Theta)^2] + 4\pi\alpha[S+\rho_{\rm ADM}] + \beta^iK_{,i}`.)
- Meaning: Z4c evolution equations in conformal variables for chi, the conformal metric and K-hat.
- Depends on: eq:Z4_ADM_1, eq:auto-3.

### eq:auto-6
```latex
\p_t \tilde{A}_{ij} = \chi[-D_iD_j\alpha + \alpha(R_{ij}-8\pi S_{ij})]^{\textrm{tf}}
+ \alpha[(\hat{K}+2\Theta)\tilde{A}_{ij} - 2\tilde{A}^k{}_i\tilde{A}_{kj}]
+ \beta^k\tilde{A}_{ij,k} + 2\tilde{A}_{(i|k}\beta^{k}{}_{,|j)}
- \tfrac{2}{3}\tilde{A}_{ij}\beta^{k}{}_{,k}
```
- Meaning: Z4c evolution of the trace-free conformal extrinsic curvature.
- Depends on: eq:Z4_ADM_1, eq:auto-3, eq:auto-8.

### eq:auto-7
```latex
\p_t \tilde{\Gamma}^{i} = -2\tilde{A}^{ij}\alpha_{,j}
+ 2\alpha[\tilde{\Gamma}^i_{jk}\tilde{A}^{jk}
- \tfrac{3}{2}\tilde{A}^{ij}\ln(\chi)_{,j}
- \tfrac{1}{3}\tilde{\gamma}^{ij}(2\hat{K}+\Theta)_{,j}
- 8\pi\tilde{\gamma}^{ij}S_j]
+ \tilde{\gamma}^{jk}\beta^i_{,jk}
+ \tfrac{1}{3}\tilde{\gamma}^{ij}\beta^k_{,kj}
+ \beta^j\tilde{\Gamma}^i_{,j}
- \tilde{\Gamma}_{\textrm{d}}{}^j\beta^i_{,j}
+ \tfrac{2}{3}\tilde{\Gamma}_{\textrm{d}}{}^i\beta^j_{,j}
```
- Meaning: Z4c evolution of the conformal connection variable.
- Depends on: eq:auto-4, eq:Z4_ADM_1, eq:Z4_ADM_2, eq:auto-3, eq:auto-5.

### eq:auto-8
```latex
R_{ij} = R^{\chi}{}_{ij} + \tilde{R}_{ij},\qquad
\tilde{R}^{\chi}{}_{ij} = \frac{1}{2\chi}\tilde{D}_i\tilde{D}_j\chi
+ \frac{1}{2\chi}\tilde{\gamma}_{ij}\tilde{D}^l\tilde{D}_l\chi
- \frac{1}{4\chi^2}\tilde{D}_i\chi\tilde{D}_j\chi
- \frac{3}{4\chi^2}\tilde{\gamma}_{ij}\tilde{D}^l\chi\tilde{D}_l\chi,\qquad
\tilde{R}_{ij} = -\tfrac{1}{2}\tilde{\gamma}^{lm}\tilde{\gamma}_{ij,lm}
+ \tilde{\gamma}_{k(i|}\tilde{\Gamma}^k_{|,j)}
+ \tilde{\Gamma}_{\textrm{d}}{}^k\tilde{\Gamma}_{(ij)k}
+ \tilde{\gamma}^{lm}\left(2\tilde{\Gamma}^k_{l(i}\tilde{\Gamma}_{j)km}
+ \tilde{\Gamma}^k_{im}\tilde{\Gamma}_{klj}\right)
```
- Meaning: split of the Ricci tensor into a conformal-factor part and a conformal-metric part, used in the A-tilde evolution.
- Depends on: eq:auto-3, eq:auto-4, eq:auto-5.

### eq:Conf_Constr_1
```latex
\Theta,\qquad 2Z_i = \tilde{\gamma}_{ij}\tilde{\Gamma}^{j}
- \tilde{\gamma}^{jk}\tilde{\gamma}_{ij,k}
```
- Meaning: the Z4 constraints expressed in conformal variables; Z_i recovered from the conformal connection.
- Depends on: eq:auto-4, eq:auto-3.

### eq:Conf_Constr_2
```latex
H = R - \tilde{A}^{ij}\tilde{A}_{ij} + \tfrac{2}{3}(\hat{K}+2\Theta)^2
- 16\pi\rho_{\rm ADM},\qquad
\tilde{M}^i = \p_j\tilde{A}^{ij} + \tilde{\Gamma}^i{}_{jk}\tilde{A}^{jk}
- \tfrac{2}{3}\tilde{\gamma}^{ij}\p_j(\hat{K}+2\Theta)
- \tfrac{3}{2}\tilde{A}^{ij}(\log\chi)_{,j},\qquad
D \equiv \ln(\det\tilde{\gamma}) = 0,\qquad
T \equiv \tilde{\gamma}^{ij}\tilde{A}_{ij} = 0
```
(The labels eq:Conf_Constr_1 and eq:Conf_Constr_2 bracket one `align` block containing the
full set of constraints in conformal variables, including the algebraic constraints D and T.)
- Meaning: Hamiltonian and momentum constraints in conformal variables plus the algebraic determinant and trace constraints, which are imposed continuously during evolution.
- Depends on: eq:auto-1, eq:auto-3.

### eq:punc_alpha
```latex
\p_t\alpha = \beta^i\alpha_{,i} - \mu_L\,\alpha^2\,\hat{K}
```
- Meaning: puncture-gauge slicing condition (generalized 1+log lapse), parametrized by mu_L.
- Depends on: eq:auto-3 (variable K-hat).

### eq:punc_beta
```latex
\p_t\beta^i = \beta^j\beta^i{}_{,j} + \mu_S\,\tilde{\Gamma}^i
- \eta\,\beta^i - \epsilon_\alpha\,\alpha\,\alpha^{,i}
+ \epsilon_\chi\,\tilde\gamma^{ij}\p_j\chi
```
- Meaning: generalized Gamma-driver shift with the new term proportional to the gradient of chi; (mu_S, eps_alpha, eps_chi) = (1,1,1/2) is the asymptotically harmonic shift.
- Depends on: eq:auto-4, eq:auto-3.

### eq:sec-lapse
```latex
\left(\p_0^2 - \mu_L\p_i\p^i\right)\alpha \simeq 0
```
- Meaning: principal part of the lapse: a wave equation with speed sqrt(mu_L).
- Depends on: eq:punc_alpha, eq:auto-2, eq:Z4_decomp_first.

### eq:sec-shift
```latex
\left(\p_0^2 - \gamma^{\frac{1}{3}}\frac{\mu_S}{\alpha^2}\p_j\p^j\right)\beta^i
\simeq \left(\frac{\gamma^{\frac{1}{3}}\mu_S}{\alpha^2\mu_L}
- \epsilon_\alpha\right)\p_0\p^i\alpha
+ \frac{1}{3\alpha}\left(\frac{\gamma^{\frac{1}{3}}\mu_S}{2}
- \epsilon_\chi\right)\gamma^{jk}\p_0\p^i\gamma_{jk}
```
- Meaning: principal part of the shift in fully second order form: a wave equation sourced by lapse and metric-trace derivatives.
- Depends on: eq:punc_beta, eq:punc_alpha, eq:Z4_ADM_1, eq:auto-2.

### eq:sec-gamma
```latex
\left(\p_0^2 - \p_l\p^l\right)\gamma_{ij} \simeq
\frac{1}{3}\gamma^{kl}\left(1 - \frac{2\epsilon_\chi}{\gamma^{\frac{1}{3}}\mu_S}\right)
\p_{i}\p_{j}\gamma_{kl}
+ \frac{2}{\alpha}\left(1 - \frac{\alpha^2\epsilon_\alpha}{\gamma^{\frac{1}{3}}\mu_S}\right)
\p_i\p_j\alpha
+ \frac{2}{\alpha}\left(1 - \frac{\alpha^2}{\gamma^{\frac{1}{3}}\mu_S}\right)
\gamma_{k(i}\p_{j)}\p_0\beta^k
```
- Meaning: principal part of the 3-metric: wave equation at speed of light with gauge sources.
- Depends on: eq:Z4_ADM_1, eq:punc_alpha, eq:punc_beta, eq:auto-2.

### eq:def_theta
```latex
2\Theta = \frac{1}{\alpha\mu_L}\p_0\alpha
- \frac{1}{2}\gamma^{ij}\p_0\gamma_{ij}
+ \frac{1}{\alpha}\p_i\beta^i
```
- Meaning: the constraint Theta viewed as defined by the gauge choice (lapse condition combined with metric trace).
- Depends on: eq:punc_alpha, eq:Z4_ADM_1, eq:auto-2.

### eq:def_Zi
```latex
2Z_i = \frac{1}{\mu_S\gamma^{1/3}}\Big(\alpha\gamma_{ij}\p_0\beta^j
+ \eta\beta_i + \epsilon_\alpha\alpha\alpha_{,i}
- \epsilon_\chi\gamma^{1/3}\p_i\chi\Big) - (\tilde\Gamma_{\textrm{d}})_i
```
- Meaning: the constraint Z_i viewed as defined by the shift condition.
- Depends on: eq:punc_beta, eq:auto-5, eq:auto-2.

### eq:sys-theta-Z
```latex
\Box\Theta \simeq 0,\qquad \Box Z_i \simeq 0
```
- Meaning: principal part of the Z4 constraint subsystem: Theta and Z_i obey wave equations and propagate at the speed of light.
- Depends on: eq:Theta_dot, eq:Z4_ADM_2.

### eq:sys-ham-mom
```latex
\Box H \simeq 0,\qquad \Box M_i \simeq 0
```
- Meaning: principal part of the Hamiltonian/momentum constraint propagation: also wave equations.
- Depends on: eq:ham-Z4, eq:mom-Z4.

### eq:auto-9
```latex
\gamma_{ss} = s^is^j\gamma_{ij},\quad \gamma_{qq} = q^{ij}\gamma_{ij},\quad
\gamma_{sA} = s^iq^j{}_A\gamma_{ij},\quad
\gamma_{AB}^{\textrm{TF}} = \left(q^i{}_Aq^j{}_B - \tfrac{1}{2}q_{AB}q^{ij}\right)\gamma_{ij},\quad
\beta_s = s^i\beta_i,\quad \beta_A = q^i{}_A\beta_i
```
- Meaning: 2+1 decomposition of the metric and shift against the spatial unit vector s^i with projector q^i_j = delta^i_j - s^i s_j.
- Depends on: none (root definition).

### eq:auto-10
```latex
\gamma_{ij} = \big(q^A{}_{(i}q^B{}_{j)} - \tfrac{1}{2}q^{AB}q_{ij}\big)\gamma^{\textrm{TF}}_{AB}
+ q^A{}_{(i}s_{j)}\gamma_{sA} + s_is_j\gamma_{ss} + q_{ij}\gamma_{qq},\qquad
\beta_{i} = s_i\beta_s + q^{A}_{i}\beta_A
```
- Meaning: reconstruction of metric and shift from the 2+1 decomposed scalar, vector and tensor pieces.
- Depends on: eq:auto-9.

### eq:auto-11
```latex
U^{\pm\sqrt{\mu_L}} = \p_0\alpha \pm \sqrt{\mu_L}\,\p_s\alpha,\qquad
U^{\pm\lambda} = \p_0\Lambda \pm \frac{\alpha^2\lambda}{\gamma^{1/3}\mu_S}\p_s\Lambda
- 2\frac{\alpha^2-\gamma^{1/3}\mu_S}{(\lambda^2\alpha^2-\gamma^{1/3}\mu_S)\alpha}
\left(\p_s\beta_s \pm \frac{\alpha^2\lambda}{\gamma^{1/3}\mu_S}\p_0\beta_s\right) + \dots,\qquad
U^{\pm 1} = \p_0\gamma_{qq} \pm \p_s\gamma_{qq},\qquad
U'^{\pm 1} = \p_0\beta_s \pm \frac{\gamma^{1/3}}{\alpha}\p_s\beta_s
+ \epsilon_\alpha\p_s\alpha \pm \frac{\gamma^{1/3}\mu_S}{\alpha^2\mu_L}\p_s\alpha
\mp \frac{\gamma^{1/3}\mu_S}{2\alpha}\p_0\Lambda
+ \frac{\epsilon_\chi-2\gamma^{1/3}\mu_S}{3\alpha}\p_s\Lambda
```
(Full lower-order coefficients of U^{±λ} as in tex, lines 505–537 of Z4.tex.)
- Meaning: scalar-sector characteristic variables of puncture-gauge Z4c with speeds ±(sqrt(mu_L), lambda, 1, 1).
- Depends on: eq:sec-lapse, eq:sec-shift, eq:sec-gamma, eq:auto-9, eq:auto-12.

### eq:auto-12
```latex
\Lambda = \gamma_{ss} + \gamma_{qq},\qquad
\lambda = \sqrt{\frac{2(2\gamma^{1/3}\mu_S - \epsilon_\chi)}{3\alpha^2}}
```
- Meaning: definitions of the metric trace variable Lambda and the gauge speed lambda; lambda = 0 makes the system only weakly hyperbolic.
- Depends on: eq:auto-9.

### eq:auto-13
```latex
U^{\pm\sqrt{\mu_S}}_{A} = \p_0\beta_A \pm \frac{\sqrt{\mu_S}\gamma^{1/3}}{\alpha}\p_s\beta_A,\qquad
U^{\pm 1}_{A} = \p_0\gamma_{sA} \pm \p_s\gamma_{sA}
- \frac{\alpha^2}{\gamma^{1/3}\mu_S}(\p_0\beta_A \pm \p_s\beta_A)
```
- Meaning: vector-sector characteristic variables with speeds ±(sqrt(mu_S), 1).
- Depends on: eq:sec-shift, eq:sec-gamma, eq:auto-9.

### eq:auto-14
```latex
U^{\pm 1}_{AB} = \p_0\gamma_{AB}^{\textrm{TF}} \pm \p_s\gamma_{AB}^{\textrm{TF}}
```
- Meaning: tensor-sector characteristic variables with speeds ±1.
- Depends on: eq:sec-gamma, eq:auto-9.

### eq:auto-15
```latex
U^\pm_{\Theta} = \p_0\Theta \pm \p_s\Theta,\quad
U^\pm_{s} = \p_0Z_s \pm \p_sZ_s,\quad
U^\pm_{A} = \p_0Z_A \pm \p_sZ_A,\quad
U^\pm_{H} = \p_0H \pm \p_sH,\quad
U'^\pm_{s} = \p_0M_s \pm \p_sM_s,\quad
U'^\pm_{A} = \p_0M_A \pm \p_sM_A
```
- Meaning: characteristic variables of the constraint subsystem, each with speeds ±1 (light speed).
- Depends on: eq:sys-theta-Z, eq:sys-ham-mom.

### eq:backg-metric
```latex
\textrm{d}\mathring{s}^2 = \mathring{g}_{ab}\textrm{d}x^a\textrm{d}x^b
= -\mathring{\alpha}^2\textrm{d}t^2
+ \mathring{\gamma}_{ij}(\textrm{d}x^i+\mathring{\beta}^i\textrm{d}t)
(\textrm{d}x^j+\mathring{\beta}^j\textrm{d}t)
```
- Meaning: background metric (alpha-ring, beta-ring, gamma-ring) against which boundary normal vectors are defined.
- Depends on: none (root definition).

### eq:auto-16
```latex
\mathring{\gamma}_{ij}\textrm{d}x^i\textrm{d}x^j
= \mathring{\psi}^4\left(\textrm{d}r^2 + r^2\textrm{d}\Omega^2\right)
```
- Meaning: assumed conformally flat form of the background 3-metric, defining the isotropic radius r.
- Depends on: eq:backg-metric.

### eq:auto-17
```latex
\mathring{g}_{ab}\mathring{n}^a\mathring{n}^b = -1,\quad
\mathring{g}_{ab}\mathring{s}^a\mathring{s}^b = 1,\quad
\mathring{g}_{ab}\mathring{n}^a\mathring{s}^b = 0
```
- Meaning: normalization of the background timelike normal n-ring and outgoing spatial normal s-ring at the boundary.
- Depends on: eq:backg-metric.

### eq:nullvector-k
```latex
\mathring{l}^a = \tfrac{1}{\sqrt{2}}(\mathring{n}^a + \mathring{s}^a),\qquad
\mathring{k}^a = \tfrac{1}{\sqrt{2}}(\mathring{n}^a + \sqrt{\nu_s}\,\mathring{s}^a)
```
(First two rows of the four-row `eqnarray`; the l-ring row is unlabeled and registered here.)
- Meaning: background outgoing characteristic vectors at unit speed (l-ring) and at the scalar-sector shift speed nu_s (k-ring).
- Depends on: eq:auto-17, eq:sec-shift.

### eq:nullvector-m
```latex
\mathring{j}^a = \tfrac{1}{\sqrt{2}}(\mathring{n}^a + \sqrt{\nu_T}\,\mathring{s}^a),\qquad
\mathring{m}^a = \tfrac{1}{\sqrt{2}}(\mathring{n}^a + \sqrt{\mu_L}\,\mathring{s}^a)
```
- Meaning: background outgoing characteristic vectors at the vector-sector shift speed nu_T (j-ring) and the lapse speed sqrt(mu_L) (m-ring).
- Depends on: eq:auto-17, eq:sec-shift, eq:sec-lapse.

### eq:general_CPBCs
```latex
\left(r^2\,\mathring{l}^a\p_a\right)^{L}\Theta \hateq 0,\qquad
\left(r^2\,\mathring{l}^a\p_a\right)^{L} Z_i \hateq 0
```
- Meaning: the high order (order L) absorbing constraint preserving boundary conditions on Theta and Z_i; L=0 recovers the Bona et al. conditions.
- Depends on: eq:sys-theta-Z, eq:nullvector-k, eq:auto-16.

### eq:BCs-alpha
```latex
\left(r^2\,\mathring{m}^a\p_a\right)^{L+1}\alpha \hateq h_\alpha
```
- Meaning: boundary condition on the lapse along the m-ring characteristic, with given data h_alpha.
- Depends on: eq:nullvector-m, eq:sec-lapse.

### eq:general_BCs_gauge_first
```latex
\left(r^2\,\mathring{k}^a\p_a\right)^{L+1}\beta_s \hateq h_s
```
- Meaning: boundary condition on the normal shift component along the k-ring characteristic.
- Depends on: eq:nullvector-k, eq:sec-shift, eq:auto-9.

### eq:BCs_lastII
```latex
\left(r^2\,\mathring{j}^a\p_a\right)^{L+1}\beta_A \hateq h_A,\qquad
\left(r^2\,\mathring{l}^a\p_a\right)^{L+1}\gamma_{AB}^{\textrm{TF}} \hateq h_{AB}^{\textrm{TF}}
```
(The beta_A row is unlabeled in tex and is registered with this node.)
- Meaning: boundary conditions on the tangential shift and the trace-free tangential metric, completing the ten conditions for the ten incoming characteristic variables.
- Depends on: eq:nullvector-k, eq:nullvector-m, eq:sec-gamma, eq:auto-9.

## Section III — Analytical setup

### eq:backg-metricP
```latex
\left.\textrm{d}\mathring{s}^2\right|_p = -\textrm{d}t^2
+ (\mathring\beta\,\textrm{d}t + \textrm{d}x)^2
+ \textrm{d}y^2 + \textrm{d}z^2
```
- Meaning: frozen-coefficient form of the background metric at a boundary point p after a slice-preserving coordinate transformation; beta-ring is the normal shift component at p, and the boundary becomes the plane x = 0.
- Depends on: eq:backg-metric.

### eq:auto-18
```latex
\mathring{n}^a\p_a = \p_t - \mathring\beta\,\p_x,\qquad
\mathring{s}^a\p_a = -\mathring\beta\,\p_x
```
- Meaning: explicit frozen-coefficient components of the background normal vectors. (Note: the tex prints s-ring as written; on the flat background the spatial normal points along -x.)
- Depends on: eq:backg-metricP, eq:auto-17.

### eq:full-sec-z4-alpha
```latex
\left(\p^2_0 - 2\,\p^l\p_l\right)\alpha = 0
```
- Meaning: frozen-coefficient lapse equation with 1+log slicing mu_L = 2/alpha at alpha=1: wave equation with speed sqrt(2).
- Depends on: eq:sec-lapse, eq:backg-metricP.

### eq:full-sec-z4-beta
```latex
\left(\p_0^2 - \p^l\p_l\right)\beta^i =
\left(\tfrac{1}{2} - \epsilon_\alpha\right)\p^i\p_0\alpha
+ \tfrac{1}{3}\left(\tfrac{1}{2} - \epsilon_\chi\right)\eta^{jk}\p_0\p^i\gamma_{jk}
```
- Meaning: frozen-coefficient shift equation with mu_S = 1; the asymptotically harmonic choice (eps_alpha, eps_chi) = (1, 1/2) decouples the gauge sector in the bulk (cascade property).
- Depends on: eq:sec-shift, eq:backg-metricP.

### eq:full-sec-z4-metric
```latex
\left(\p_0^2 - \p^l\p_l\right)\gamma_{ij} = \tfrac{1}{3}(1 - 2\epsilon_\chi)
\eta^{kl}\p_i\p_j\gamma_{kl} + 2(1 - \epsilon_\alpha)\p_i\p_j\alpha
```
- Meaning: frozen-coefficient metric equation; sources vanish for the asymptotically harmonic shift.
- Depends on: eq:sec-gamma, eq:backg-metricP.

### eq:2+1Z4_a
```latex
\left[\p_t^2 - 2\mathring\beta\,\p_t\p_x - (2-\mathring\beta^2)\p_x^2
- 2\p^A\p_A\right]\alpha = 0
```
- Meaning: 2+1 decomposed frozen-coefficient wave problem for the lapse.
- Depends on: eq:full-sec-z4-alpha, eq:auto-18, eq:auto-9.

### eq:2+1Z4_ac
```latex
\left[\p_t - (\sqrt{2}+\mathring{\beta})\p_x\right]^{L+1}\alpha \hateq h_\alpha
```
- Meaning: frozen-coefficient form of the lapse boundary condition.
- Depends on: eq:BCs-alpha, eq:backg-metricP.

### eq:Z4_betas_2+1
```latex
\Big[\p_t^2 - 2\mathring\beta\,\p_t\p_x - (1-\mathring\beta^2)\p_x^2
- \p^C\p_C\Big]\beta_{i} = -\tfrac{1}{2}\p_0\p_i\alpha
```
- Meaning: 2+1 decomposed wave problem for the shift components (project along s or transverse), coupled to the lapse only in the bulk.
- Depends on: eq:full-sec-z4-beta, eq:auto-9, eq:auto-18.

### eq:Z4_betas_2+1BC
```latex
\left[\p_t - (1+\mathring{\beta})\p_x\right]^{L+1}\beta_{i} \hateq h_{i}
```
- Meaning: frozen-coefficient boundary condition for the shift components.
- Depends on: eq:general_BCs_gauge_first, eq:BCs_lastII, eq:backg-metricP.

### eq:Z4_gammaAB_2+1
```latex
\Big[\p_t^2 - 2\mathring\beta\,\p_t\p_x - (1-\mathring\beta^2)\p_x^2
- \p^C\p_C\Big]\gamma^{\textrm{TF}}_{AB} = 0
```
- Meaning: decoupled wave problem for the trace-free tangential metric; well-posedness follows directly from the toy model of appendix B.
- Depends on: eq:full-sec-z4-metric, eq:auto-9, eq:auto-18.

### eq:Z4_gammaAB_2+1BC
```latex
\left[\p_t - (1+\mathring{\beta})\p_x\right]^{L+1}\gamma^{\textrm{TF}}_{AB}
\hateq h^{\textrm{TF}}_{AB}
```
- Meaning: frozen-coefficient boundary condition for gamma_AB^TF.
- Depends on: eq:BCs_lastII, eq:backg-metricP.

### eq:EOM-Lambda
```latex
\left[\p_t^2 - 2\mathring\beta\,\p_t\p_x - (1-\mathring\beta^2)\p_x^2
- \p^A\p_A\right]\Lambda = 0
```
- Meaning: wave problem for the metric trace variable Lambda, decoupled from the rest of the metric sector in the bulk.
- Depends on: eq:full-sec-z4-metric, eq:auto-12, eq:auto-18.

### eq:BC-Lambda
```latex
\left[\p_t - (1+\mathring{\beta})\p_x\right]^{L}
\Big[(\p_t - \mathring{\beta}\p_x)(\alpha - \Lambda)
- 2\p_x\beta_s + 2\p^A\beta_A\Big] \hateq 0
```
- Meaning: boundary condition for Lambda obtained from the CPBC on Theta using the gauge definition of Theta.
- Depends on: eq:general_CPBCs, eq:def_theta, eq:auto-12, eq:backg-metricP.

### eq:eomgammasA
```latex
\Big[\p_t^2 - 2\mathring\beta\,\p_t\p_x - (1-\mathring\beta^2)\p_x^2
- \p^C\p_C\Big]\gamma_{is} = 0
```
- Meaning: wave problem for the mixed metric components gamma_ss, gamma_sA (project along s or transverse); equations of motion decoupled but boundary conditions mutually coupled.
- Depends on: eq:full-sec-z4-metric, eq:auto-9, eq:auto-18.

### eq:Z4_gammasbc_2+1
```latex
\left[\p_t - (1+\mathring{\beta})\p_x\right]^L
\Big[(\p_t - \mathring{\beta}\p_x)\beta_s - \p^A\gamma_{sA}
- \p_x(\alpha - \gamma_{ss} + \Lambda/6)\Big] \hateq 0
```
- Meaning: boundary condition for gamma_ss obtained from the CPBC on Z_s using the gauge definition of Z_i.
- Depends on: eq:general_CPBCs, eq:def_Zi, eq:backg-metricP, eq:auto-12.

### eq:Z4_gammaAbc_2+1
```latex
\left[\p_t - (1+\mathring{\beta})\p_x\right]^L
\Big[(\p_t - \mathring{\beta}\p_x)\beta_A - \p^B\gamma^{\textrm{TF}}_{AB}
+ \p_x\gamma_{sA} + \p_A(\alpha - \Lambda/3 + \gamma_{ss}/2)\Big] \hateq 0
```
- Meaning: boundary condition for gamma_sA from the CPBC on Z_A; coupled to the gamma_ss boundary condition, so the two wave problems must be treated simultaneously.
- Depends on: eq:general_CPBCs, eq:def_Zi, eq:backg-metricP, eq:auto-12.

## Section IV — Well-posedness results

(The Laplace–Fourier transformed Theta wave problem appears in a `displaymath`
carrying the non-environment label `eq:theta-eq-LF`:
`[(s^2+\omega^2) - 2\mathring\beta s\p_x - (1-\mathring\beta^2)\p_x^2]\tilde\Theta = 0`,
`\tilde\Theta \hateq \tilde q`; it is used as input to eq:auto-19 and eq:MLF-Mat.)

### eq:auto-19
```latex
D\tilde\Theta = \frac{1}{\kappa}\left(\p_x\tilde\Theta
+ \gamma^2_\mu\,\mathring\beta\,s\,\tilde\Theta\right)
```
- Meaning: first order pseudo-differential reduction variable for the transformed Theta problem, with kappa = sqrt(|s|^2 + omega^2).
- Depends on: eq:sys-theta-Z, eq:backg-metricP, eq:generalform-U.

### eq:MLF-Mat
```latex
M(s,\omega) = \kappa\left(\begin{array}{cc}
-\gamma^2\,\mathring\beta\,s' & 1 \\
\gamma^4\,\lambda^2 & -\gamma^2\,\mathring\beta\,s'
\end{array}\right)
```
with W-tilde = (Theta-tilde, D Theta-tilde)^T, L(s,omega) = (1,0), gamma = 1/sqrt(1-beta-ring^2),
s' = s/kappa, omega' = omega/kappa, lambda^2 = s'^2 + gamma^{-2} omega'^2.
- Meaning: companion matrix of the first order reduction of the Theta constraint wave problem in Laplace–Fourier space.
- Depends on: eq:auto-19, eq:sys-theta-Z, eq:generalform-U, eq:backg-metricP.

### eq:auto-20
```latex
\tau_\pm = -\kappa\,\gamma^2\,(s'\,\mbeta \mp \lambda),\qquad
\hat{e}_\pm = (1, \pm\gamma^2\,\lambda)^T
```
- Meaning: eigenvalues and eigenvectors of M(s,omega) for the Theta constraint problem.
- Depends on: eq:MLF-Mat.

### eq:general_sol_theta
```latex
\tilde{W}(s,x,\omega) = \sigma\,\hat e^-\,\text{exp}(\tau^-\,x)
```
- Meaning: the L2 (decaying) solution of the transformed Theta problem; sigma is fixed by the boundary condition, sigma = q-tilde for L = 0.
- Depends on: eq:auto-20, eq:solution_lap-four.

### eq:fcpbcTheta
```latex
|\tilde W(s,0,\omega)| \leq C\,|\tilde q|
```
- Meaning: boundary-stability estimate for the Theta constraint with first order CPBCs: boundary values bounded by the boundary data.
- Depends on: eq:general_sol_theta, eq:general_CPBCs, eq:boundary_stable.

### eq:esTheta
```latex
\int_0^T\|W(\cdot,t)\|^2_\Sigma\,dt \leq C_T\int_0^T\|q\|^2_{\partial\Sigma}\,dt
```
- Meaning: physical-space L2 estimate from inverting the Laplace transform and Parseval; with q = 0 and constraint-satisfying initial data the constraints stay zero. Same holds for Z_i.
- Depends on: eq:fcpbcTheta, eq:estWoutS.

### eq:hoTheta-bc
```latex
L(s,\omega) = \frac{1}{2}\left(a_+^{L}+a_-^{L},
-\frac{a_+^{L}-a_-^{L}}{\lambda\,\gamma^{2}}\right)
```
- Meaning: boundary matrix for the high order (L >= 1) CPBCs on Theta, with a_pm = s' ± lambda, obtained via the appendix-B iteration procedure.
- Depends on: eq:h-bc, eq:auto-19, eq:general_CPBCs, eq:MLF-Mat.

### eq:auto-21
```latex
a_+^L\,\sigma = \tilde q
```
- Meaning: the integration constant for the high order Theta problem; since |a_+| > delta > 0 the system remains boundary stable, so the high order CPBC IBVP is well-posed and satisfies eq:esTheta.
- Depends on: eq:hoTheta-bc, eq:general_sol_theta, eq:bc_est-bc.

### eq:auto-22
```latex
D\tilde\alpha = \frac{1}{\kappa}\left(\p_x + \gamma^2_\alpha\,\mathring{\beta}\,s\right)\tilde\alpha,\qquad
D\tilde\beta_s = \frac{1}{\kappa}\left(\p_x + \gamma^2\,\mathring{\beta}\,s\right)\tilde\beta_s
```
- Meaning: first order reduction variables for the spherical reduction of the lapse/shift (gauge) subsystem.
- Depends on: eq:2+1Z4_a, eq:Z4_betas_2+1, eq:generalform-U.

### eq:auto-23
```latex
\tilde{W} = (\tilde\alpha, D\tilde\alpha, \tilde\beta_s, D\tilde\beta_s)^T,\qquad
M(s) = \kappa\left(\begin{array}{cccc}
-\gamma_\alpha^2\mathring\beta s' & 1 & 0 & 0 \\
\gamma^4_\alpha s'^2_\alpha & -\gamma^2_\alpha\mathring\beta s' & 0 & 0 \\
0 & 0 & -\gamma^2\mathring\beta s' & 1 \\
-2s'^2\mathring\beta\gamma^2\gamma^4_\alpha & s'(2+\mathring\beta^2)\gamma^2\gamma^2_\alpha/2
& \gamma^4 s'^2 & -\gamma^2\mathring\beta s'
\end{array}\right),\qquad
L(s) = \left(\begin{array}{cccc}
\sqrt{2}s' & -\gamma^{-2}_\alpha & 0 & 0 \\
0 & 0 & s' & -\gamma^{-2}
\end{array}\right),\qquad
\tilde{g} = \frac{1}{\kappa}\left(\begin{array}{c}
(\sqrt{2}-\mathring\beta)\tilde h_{\alpha} \\
(1-\mathring\beta)\tilde h_{s}
\end{array}\right)
```
- Meaning: first order form (state vector, companion matrix, boundary matrix, boundary data) of the spherically reduced gauge subsystem with first order BCs; here kappa = |s|.
- Depends on: eq:auto-22, eq:2+1Z4_a, eq:Z4_betas_2+1, eq:2+1Z4_ac, eq:Z4_betas_2+1BC, eq:generalform-U-BC.

### eq:auto-24
```latex
\gamma^{-2}_\alpha = 2 - \mathring\beta^2,\qquad \gamma^{-2} = 1 - \mathring\beta^2
```
- Meaning: definitions of the boost-type factors associated with the lapse speed sqrt(2) and the shift/metric speed 1.
- Depends on: eq:backg-metricP.

### eq:genralsol-gauge
```latex
\tilde{W}(s,x,\omega) = \sum_{i=1}^2 \sigma_i\,\hat e_i^-\,\textrm{exp}(\tau_i^-\,x)
```
- Meaning: L2 solution of the gauge subsystem, a superposition over the two decaying eigenmodes of M.
- Depends on: eq:auto-23, eq:solution_lap-four.

### eq:auto-25
```latex
\sigma_\alpha = \frac{\tilde g_\alpha}{2\sqrt{2}\,s'},\qquad
\sigma_s = \frac{\tilde g_s}{2\,s'} -
\frac{\tilde g_\alpha\left[(1-\mbeta)(2+\sqrt{2}+\mbeta(1+\sqrt{2})\right]\gamma^2_\alpha}{8\,s'}
```
- Meaning: integration constants for first order BCs; by the triangle inequality the gauge subsystem is boundary stable, so the estimate eq:esTheta holds for it.
- Depends on: eq:genralsol-gauge, eq:auto-23.

### eq:auto-26
```latex
\mathcal{L}_\mu = (\mu - \mbeta)\,s' - \frac{1}{\kappa\,\gamma^2_\mu}\,\p_x
```
- Meaning: pseudo-differential boundary operator at speed mu; the gauge boundary operator is calligraphic-L = (L_sqrt2, L_1)^T. Using the equations of motion the L=0 boundary becomes the algebraic form calligraphic-L W-tilde = A W-tilde (4x4 matrix A given in displaymath), and iteration gives calligraphic-L^{L+1} W-tilde = A^{L+1} W-tilde with eigenvalues a_+ = 2 sqrt(2) s', b_+ = 2 s' and shorthands F,G,H,J of beta-ring (A^{L+1} given in widetext displaymath).
- Depends on: eq:auto-23, eq:high-orderwave, eq:full-sec-z4-alpha, eq:full-sec-z4-beta.

### eq:auto-27
```latex
L(s,\omega) = \frac{1}{2}\left(\begin{array}{cccc}
a^{L+1}_+ & -a^{L+1}_+/(\sqrt{2}s'\gamma_\alpha^2) & 0 & 0 \\
-F(\mbeta)\gamma^4_\alpha a^{L+1}_+/\gamma^2 &
G(\mbeta)\gamma^2_\alpha a^{L+1}/(s'\gamma^2) & b^{L+1}_+ & -b^{L+1}_+/(s'\gamma^2)
\end{array}\right)
```
- Meaning: boundary matrix of the gauge subsystem with high order (order L+1) BCs in algebraic form.
- Depends on: eq:auto-26, eq:auto-23.

### eq:auto-28
```latex
\left(\begin{array}{c}\tilde g_\alpha \\ \tilde g_{s}\end{array}\right)
= \frac{1}{\kappa^{L+1}}\left(\begin{array}{c}
(\sqrt{2}-\mathring\beta)^{L+1}\tilde h_{\alpha} \\
(1-\mathring\beta)^{L+1}\tilde h_{s}
\end{array}\right)
```
- Meaning: rescaled boundary data for the high order gauge BCs.
- Depends on: eq:auto-27, eq:auto-23.

### eq:auto-29
```latex
\sigma_\alpha = \frac{\tilde g_\alpha}{2\,a^{L+1}},\qquad
\sigma_s = \tilde g_\alpha\left[
\frac{F(\mbeta)+\sqrt{2}\,G(\mbeta)\,\gamma_a^4}{4\,b^{L+1}\,\gamma^2}
- \frac{2(1+\sqrt{2})(1-\mbeta)}{8(2-\sqrt{2}\mbeta)\,a^{L+1}}\right]
+ \frac{\tilde g_s}{2\,b^{L+1}}
```
- Meaning: integration constants for the high order gauge BCs; since Re(s') > 0 and a^{L+1}, b^{L+1} are proportional to s', the gauge subsystem is boundary stable, hence well-posed; the metric wave problems follow similarly.
- Depends on: eq:auto-27, eq:auto-28, eq:genralsol-gauge.

## Section V — Numerical applications

### eq:Cmonitor
```latex
C \equiv \sqrt{H^2 + M^i M_i + \Theta^2 + Z^i Z_i}
```
- Meaning: pointwise constraint monitor combining all Z4c constraints.
- Depends on: eq:auto-1, eq:Conf_Constr_1, eq:Conf_Constr_2.

### eq:norm2
```latex
||C(\cdot,t)||_2 \equiv \sqrt{\int dr\, r^2\, C(r,t)^2}
```
- Meaning: radial 2-norm of the constraint monitor evaluated by the trapezium rule on the grid.
- Depends on: eq:Cmonitor.

### eq:theta _modes
```latex
U^\pm_\Theta = \partial_0 \Theta \pm \Theta_{,s}
```
- Meaning: incoming/outgoing characteristic fields of Theta used to quantify boundary absorption.
- Depends on: eq:auto-15.

### eq:refcoef
```latex
R \equiv \frac{|\tilde{U}^-_\Theta(k)|}{|\tilde{U}^+_\Theta(k)|}
```
- Meaning: experimental reflection coefficient: ratio of Fourier modes of incoming to outgoing Theta characteristics at the boundary.
- Depends on: eq:theta _modes.

### eq:auto-30
```latex
ds^2 = -\left(1-\frac{2M}{r}\right)dt^2 + \frac{4M}{r}\,dt\,dr
+ \left(1+\frac{2M}{r}\right)dr^2 + r^2\,d\Omega^2
```
- Meaning: Kerr-Schild (ingoing Eddington-Finkelstein) Schwarzschild initial data used in the excision black hole test.
- Depends on: none (root; standard exact solution).

## Appendix A — Well-posed problems (Kreiss theory)

### eq:def-system
```latex
\partial_t u = B\,\partial_{x^1} u + \sum_{A=2}^n C^A\,\partial_A u
```
- Meaning: generic constant-coefficient strongly hyperbolic first order system on the half-space x^1 >= 0.
- Depends on: none (root; Kreiss theory input).

### eq:formA
```latex
B = \left(\begin{array}{cc}-\Lambda^I & 0 \\ 0 & \Lambda^{II}\end{array}\right)
```
- Meaning: block-diagonal form of the non-singular boundary matrix B, separating incoming (Lambda^I, m of them) and outgoing modes.
- Depends on: eq:def-system.

### eq:BC
```latex
\left.L\,u^I(t,x)\right|_{x^1=0} \hateq g(t,x)
```
- Meaning: m boundary conditions imposed on the incoming part u^I with given data g.
- Depends on: eq:def-system, eq:formA.

### eq:lapl-four
```latex
\partial_x\tilde u = M(s,\omega)\,\tilde u,\qquad \textrm{on}\ x\in(0,\infty)
```
- Meaning: ODE in x satisfied by the Laplace–Fourier transform of u.
- Depends on: eq:def-system.

### eq:BC-lapl-four
```latex
L\tilde{u}^I \hateq \tilde{g},\qquad \textrm{at}\ x \hateq 0
```
- Meaning: transformed boundary condition.
- Depends on: eq:BC, eq:lapl-four.

### eq:auto-31
```latex
M(s,\omega) = B^{-1}\,(s\,\mathbb{I}_{n\times n} + i\,\omega_A\,C^A)
```
- Meaning: explicit symbol matrix of the transformed system.
- Depends on: eq:lapl-four, eq:def-system.

### eq:solution_lap-four
```latex
\tilde{u} = \sum_{i=1}^m \sigma_i\,e_i(s,\omega)\,\textrm{exp}(\lambda_i\,x^1)
```
- Meaning: general L2 solution as a sum over decaying eigenmodes of M; sigma_i fixed by the boundary conditions.
- Depends on: eq:lapl-four, eq:auto-31.

### eq:bc_general
```latex
\mathbb{D}(s,\omega)\,\sigma \hateq \tilde{g}
```
- Meaning: linear system for the integration constants obtained by inserting the mode solution into the boundary condition.
- Depends on: eq:solution_lap-four, eq:BC-lapl-four.

### sol:lf
```latex
u(t,x^i) = \tilde{u}(x^1)\,\textrm{exp}\left(s\,t + i\,\omega_A\,x^A\right)
```
- Meaning: nontrivial homogeneous solution that exists when Det D = 0 for some Re(s) > 0.
- Depends on: eq:bc_general, eq:solution_lap-four.

### sol:lf2
```latex
u(t,x^i) = \tilde{u}(\theta\,x^1)\,
\textrm{exp}\left(\zeta\,s\,t + i\,\zeta\,\omega_A\,x^A\right)
```
- Meaning: rescaled solution for any zeta > 0 by homogeneity; arbitrarily fast exponential growth, so Det D != 0 for Re(s) > 0 (determinant condition, stated in displaymath) is necessary for well-posedness.
- Depends on: sol:lf, eq:lapl-four.

### eq:boundary_stable
```latex
|\tilde{u}(s,0,\omega)| \leq C\,|\tilde{g}(s,\omega)|
```
- Meaning: boundary stability: the solution at the boundary is bounded by the boundary data, uniformly in s and omega.
- Depends on: eq:bc_general, eq:solution_lap-four.

### eq:estWoutS
```latex
\int_0^T\|u(t,\cdot)\|^2\,dt \leq \delta\int_0^T\left.\|g(t,\cdot)\|^2\right|_{x^1=0}\,dt
```
- Meaning: physical-space L2 estimate equivalent to boundary stability after inverting the Laplace–Fourier transform.
- Depends on: eq:boundary_stable.

### eq:gen_estimate
```latex
\int_0^T\|u(t,\cdot)\|^2\,dt + \int_0^T\left.\|u(t,\cdot)\|^2\right|_{x^1=0}\,dt \leq
\delta\left(\int_0^T\|F(t,\cdot)\|^2\,dt
+ \int_0^T\left.\|g(t,\cdot)\|^2\right|_{x^1=0}\,dt\right)
```
- Meaning: full Kreiss estimate including a bulk source F, obtained from the smooth symmetrizer H(s',omega') whose existence follows from boundary stability for strictly hyperbolic systems.
- Depends on: eq:boundary_stable, eq:estWoutS.

## Appendix B — Toy model (shifted wave equation with high order BCs)

### eq:eqMU
```latex
\left[\partial^2_0 - \mu^2\,\partial^l\partial_l\right]U(t,x^i) = F(t,x^i)
```
- Meaning: model wave equation at speed mu with source F on the half-space, trivial initial data.
- Depends on: eq:auto-2 (the p_0 operator); otherwise root.

### eq:bc_phi
```latex
\left[\partial_0 - \mu\,\partial_x\right]U(t,x^i) \hateq h
```
- Meaning: first order outgoing (Sommerfeld-type) boundary condition for the toy wave equation.
- Depends on: eq:eqMU.

### eq:LF-u
```latex
\left[(\mu^2-\mathring\beta^2)\p_x^2 + 2\mathring\beta\,s\,\p_x
- (s^2 + \mu^2\omega^2)\right]\tilde U = \tilde F
```
- Meaning: Laplace–Fourier transform of the toy wave equation in the frozen background eq:backg-metricP.
- Depends on: eq:eqMU, eq:backg-metricP.

### eq:BC-u
```latex
\left[s - (\mu + \mathring{\beta})\p_x\right]\tilde U \hateq \tilde h
```
- Meaning: transformed boundary condition.
- Depends on: eq:bc_phi, eq:backg-metricP.

### eq:auto-32
```latex
D\tilde U = \frac{1}{\kappa}\left(\p_x\tilde U
+ \gamma^2_\mu\,\mathring\beta\,s\,\tilde U\right)
```
- Meaning: first order reduction variable for the toy problem, with gamma_mu = 1/sqrt(mu^2 - beta-ring^2).
- Depends on: eq:LF-u.

### eq:generalform-U
```latex
\partial_x\tilde{W} = M(s,\omega)\,\tilde{W} + \tilde{f}
```
with W-tilde = (U-tilde, D U-tilde)^T and f-tilde = (gamma_mu^2/kappa)(0, F-tilde)^T.
- Meaning: canonical first order form of the transformed toy wave problem.
- Depends on: eq:LF-u, eq:auto-32.

### eq:generalform-U-BC
```latex
L(s,\omega)\tilde W = \tilde g
```
- Meaning: canonical algebraic form of the boundary condition for the first order reduction.
- Depends on: eq:BC-u, eq:auto-32.

### eq:LF-Mat
```latex
M(s,\omega) = \kappa\left(\begin{array}{cc}
-\gamma_\mu^2\,\mathring\beta\,s' & 1 \\
\mu^2\,\gamma^4_\mu\,\lambda^2 & -\gamma_\mu^2\,\mathring\beta\,s'
\end{array}\right),\qquad
L(s,\omega) = (\mu\,s', -\gamma^{-2}_\mu)
```
with lambda^2 = s'^2 + gamma_mu^{-2} omega'^2.
- Meaning: explicit companion and boundary matrices of the toy problem.
- Depends on: eq:generalform-U, eq:generalform-U-BC, eq:auto-32.

### eq:general_sol_phi
```latex
\tilde{W}(s,x,\omega) = \sigma\,e^{(\tau^- x)}\,e^-
```
- Meaning: L2 solution along the decaying eigenvector e^- of M (Re(tau^-) < 0).
- Depends on: eq:generalform-U, eq:LF-Mat, eq:solution_lap-four.

### eq:bc_est-phi
```latex
\mu\left(s' + \sqrt{s'^2 + \gamma^{-2}_\mu\,\omega'^2}\right)\sigma = \tilde g
```
- Meaning: boundary equation determining sigma for the first order BC.
- Depends on: eq:general_sol_phi, eq:LF-Mat.

### eq:bc_est-bc
```latex
\left|s' + \sqrt{s'^2 + \gamma^{-2}_\mu\,\omega'^2}\right| \geq \delta
```
- Meaning: lower bound on the boundary symbol, uniform for Re(s) > 0, the key non-degeneracy ensuring boundary stability.
- Depends on: eq:bc_est-phi.

### eq:Bon-stab-phi
```latex
|\tilde{W}(s,0,\omega)| \leq C\,|\tilde{g}|
```
- Meaning: boundary stability of the toy problem (with properly normalized eigenvector e^-).
- Depends on: eq:bc_est-bc, eq:general_sol_phi.

### eq:phi-estimate
```latex
\eta\int_0^\infty(|\kappa\tilde U|^2 + |\partial_x\tilde{U}|^2)\,dx
+ \left.(|\kappa\tilde{U}|^2 + |\partial_x\tilde{U}|^2)\right|_{x=0}
\leq C'\left[\frac{1}{\eta}\int_0^\infty|\tilde F|^2\,dx + |\tilde{h}|^2\right]
```
- Meaning: frequency-domain energy estimate for the sourced toy problem from a symmetrizer; inverting the transform and Parseval yields eq:gen_estimate.
- Depends on: eq:Bon-stab-phi, eq:gen_estimate.

### eq:bc_phi-ho
```latex
\left[\partial_0 - \mu\,\partial_x\right]^{m+1}U(t,x^i) \hateq h
```
- Meaning: high order (m >= 1) boundary conditions for the toy problem, known to reduce spurious reflections.
- Depends on: eq:bc_phi.

### eq:high-orderwave
```latex
\mathcal{L}^{m+1}\tilde U \hateq \left[\frac{\mu-\mathring\beta}{\kappa}\right]^{L+1} h
```
with calligraphic-L = (mu - beta-ring) s' - p_x/(kappa gamma_mu^2).
- Meaning: transformed high order boundary condition as a power of the pseudo-differential operator L.
- Depends on: eq:bc_phi-ho, eq:LF-u.

### eq:auto-33
```latex
\mathcal{L}\left(\begin{array}{c}\tilde U \\ D\tilde U\end{array}\right)
= A\left(\begin{array}{c}\tilde U \\ D\tilde U\end{array}\right)
- \frac{1}{\kappa^2}\left(\begin{array}{c}0 \\ \tilde F\end{array}\right)
```
- Meaning: action of the boundary operator reduced to an algebraic matrix A using the equation of motion.
- Depends on: eq:high-orderwave, eq:generalform-U.

### eq:A-BC
```latex
A = \left(\begin{array}{cc}
\mu\,s' & -\gamma_\mu^{-2} \\
-\mu^2\,\gamma^2_\mu\,\lambda^2 & \mu\,s'
\end{array}\right)
```
- Meaning: the algebraic boundary matrix whose eigenvalues are a_pm = mu(s' ± lambda).
- Depends on: eq:auto-33.

### eq:h-bc
```latex
L(s,\omega) = \frac{1}{2}\left(a_+^{m+1}+a_-^{m+1},
-\frac{a_+^{m+1}-a_-^{m+1}}{\mu\,\lambda\,\gamma_\mu^{2}}\right)
```
- Meaning: boundary matrix of the high order toy BCs by iteration of A; the integration constant satisfies a_+^{m+1} sigma = g-tilde and the system is boundary stable.
- Depends on: eq:A-BC, eq:high-orderwave.

### eq:full-estimate-W
```latex
\eta\int_0^\infty \sum_{j=0}^{m+1}|\kappa^{(m+1)-j}\p_x^j\tilde U|^2\,dx
+ \left.\sum_{j=0}^{m+1}|\kappa^{(m+1)-j}\p_x^j\tilde U|^2\right|_{x=0}
\leq C\left[\frac{1}{\eta}\int_0^\infty\sum_{j=0}^{m-1}|\kappa^{m-j}\partial_x^j\tilde F|^2\,dx
+ \left.\sum_{j=0}^{m-1}|\kappa^{(m-1)-j}\p_x^j\tilde F|^2\right|_{x=0}
+ |\tilde{h}|^2\right]
```
- Meaning: high order energy estimate: L2 control of higher derivatives of the solution in terms of given data.
- Depends on: eq:h-bc, eq:general_sol_phi, eq:Bon-stab-phi.

## Appendix C — Implementation of BCs in spherical symmetry

### eq:auto-34
```latex
\textrm{d}s^2 = \chi^{-1}\tilde{\gamma}_{rr}\textrm{d}r^2
+ \chi^{-1}\tilde{\gamma}_T r^2\textrm{d}\Omega^2
```
- Meaning: spherically symmetric line element in Z4c conformal variables (gamma-tilde_rr, gamma-tilde_T, chi); extrinsic curvature carried by (K-hat, A-tilde_rr, A-tilde_T).
- Depends on: eq:auto-3.

### eq:auto-35
```latex
D = \log\left(\tilde{\gamma}_{rr}\tilde{\gamma}_{T}^2\right) = 0,\qquad
T = \frac{\tilde{A}_{rr}}{\tilde{\gamma}_{rr}} + 2\frac{\tilde{A}_{T}}{\tilde{\gamma}_{T}} = 0
```
- Meaning: the algebraic constraints D and T specialized to spherical symmetry.
- Depends on: eq:Conf_Constr_2, eq:auto-34.

### eq:sph_bc_1st
```latex
\p_t\Theta \hateq -\p_r\Theta - \frac{1}{r}\Theta
```
- Meaning: first of the implemented second order CPBCs in spherical symmetry, linearized around flat space; Sommerfeld-like form replacing the Theta evolution equation at the boundary.
- Depends on: eq:general_CPBCs, eq:auto-34, eq:Theta_dot.

### eq:sph_bc_last
```latex
\p_t\tilde{\Gamma}^r \hateq -\frac{2}{\sqrt{3}}\p_r\tilde{\Gamma}^r
-\frac{2}{\sqrt{3}r}\tilde{\Gamma}^r - \frac{4}{3r^2}\beta^r
-\frac{1}{3}\p_r\Theta - \frac{2}{3}\p_r\hat{K},\qquad
\p_t\hat{K} \hateq -\sqrt{2}\,\p_r\hat{K} - \frac{\sqrt{2}}{r}\hat{K}
+ \frac{1}{r}\p_r\alpha,\qquad
\tilde{A}_{rr} \hateq -2\p_r\tilde{A}_{rr} - \frac{6}{r}\tilde{A}_{rr}
- \p_r\p_r\tilde{\gamma}_T - \frac{1}{2}\p_r\tilde{\Gamma}^r
- \frac{2}{3}(-2+\sqrt{2})\p_r\hat{K} + \frac{1}{3}\p_r\Theta
- \frac{2}{r^2}(\tilde{\gamma}_{rr}-\tilde{\gamma}_{T})
+ \frac{1}{r}\big(\frac{5}{2}\p_r\tilde{\gamma}_{rr} + 2\p_r\alpha
- \frac{1}{3}\Theta - \frac{2\sqrt{2}}{3}\hat{K} - 2\tilde{\Gamma}^r
- \p_r\tilde{\gamma}_T - \p_r\chi\big)
```
(The intermediate Gamma-tilde^r and K-hat rows of this `align` block are unlabeled in tex and
registered with this node; the label eq:sph_bc_last sits on the A-tilde_rr row.)
- Meaning: remaining implemented boundary conditions for (Gamma-tilde^r, K-hat, A-tilde_rr); A-tilde_T follows from the algebraic constraints. These replace the standard evolution equations at the outermost grid point.
- Depends on: eq:general_CPBCs, eq:BCs-alpha, eq:general_BCs_gauge_first, eq:auto-34, eq:auto-35, eq:sph_bc_1st.

### eq:auto-36
```latex
f_{N+i} = 6f_{N+i-1} - 15f_{N+i-2} + 20f_{N+i-3} - 15f_{N+i-4}
+ 6f_{N+i-5} - f_{N+i-6}
```
- Meaning: sixth order extrapolation used to populate ghostzones at the boundary point N for derivative and Kreiss–Oliger dissipation evaluation.
- Depends on: none (root; standard numerical stencil).
