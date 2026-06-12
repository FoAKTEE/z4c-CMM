# logic.md — equation dependency DAG for arXiv:2308.10361

One node per equation environment (node id = equation label; unlabeled
environments use eq:auto-N). `<- []` marks a root. Grouped by paper section.
This file is the source of truth mirrored by
`knowledge-database/paper_arxiv-2308.10361/nodes.jsonl`
(task_id eq-dag-stage1-2308.10361). All edges are literature-grounded
transcriptions of the dependency structure stated in the paper.

## Sec. II (sec:GH)

- `eq:3+1_metric` <- []
  - 3+1 (ADM) decomposition of the spacetime metric in Cauchy (primed) coordinates: lapse, shift, spatial metric
- `FOSH` <- ['eq:3+1_metric']
  - Vacuum Einstein equations as first-order symmetric hyperbolic generalized-harmonic system for {g,Pi,Phi}
- `eq:def_s` <- ['eq:3+1_metric']
  - Outward-directed spatial unit normal to the Cauchy outer boundary
- `eq:auto-1` <- ['FOSH', 'eq:def_s']
  - Left eigenvectors and characteristic speeds of the FOSH principal symbol along the boundary normal
- `eq:auto-2` <- ['eq:auto-1', 'eq:def_s']
  - Definitions of characteristic-projected time derivative d_t and normal derivative d_perp
- `eq:auto-3` <- ['FOSH', 'eq:auto-1']
  - Definition of D_t u: characteristic projection of the full FOSH right-hand side
- `Bjorhus_bc` <- ['eq:auto-1', 'eq:auto-2', 'eq:auto-3']
  - Bjorhus boundary method: replace normal derivative of incoming characteristic fields by desired boundary value
- `eq:projection_operator` <- ['eq:3+1_metric', 'eq:def_s']
  - 2D projector orthogonal to hypersurface normal n and boundary normal s
- `eq:auto-4` <- ['eq:projection_operator']
  - Physical (transverse-traceless) projection operator built from the 2D projector
- `eq:wab_projection` <- ['eq:auto-4', 'eq:projection_operator']
  - Inward-propagating Weyl components w-: TT-projected double contraction of Weyl tensor with n+s
- `eq:bc_bjorhus` <- ['Bjorhus_bc', 'eq:auto-4', 'eq:wab_projection']
  - Physical boundary condition on u^1- driving w- toward its CCM-supplied boundary value
- `eq:gh_tetrad_l` <- ['eq:def_s']  (labels: eq:gh_tetrad_l, eq:gh_tetrad_k)
  - Unique Cauchy outgoing/ingoing null tetrad vectors l,k from time normal n and boundary normal s
- `eq:GH_m` <- ['eq:gh_tetrad_l']
  - Orthonormality conditions on the Cauchy angular tetrad vector m (fixed up to phase rotation)
- `eq:GH_psi0_def` <- ['eq:gh_tetrad_l', 'eq:GH_m']
  - Definition of the Cauchy Weyl scalar psi0' in the GH tetrad (ingoing transverse radiation)
- `eq:auto-5` <- ['eq:projection_operator', 'eq:GH_m']
  - Completeness identity P = m mbar + mbar m for the 2D projector
- `eq:wab` <- ['eq:wab_projection', 'eq:auto-5', 'eq:GH_psi0_def']
  - w- rewritten in terms of psi0' and m; the slot where the characteristic psi0 is inserted for CCM
## Sec. III (sec:CCE)

- `eq:BS_PFB` <- []
  - Bondi-Sachs metric in partially flat Bondi-like (hatted) coordinates used by SpECTRE CCE
- `eq:BS_angular_gauge_condition` <- ['eq:BS_PFB']
  - Bondi determinant gauge condition: det(h_AB) equals unit-sphere determinant
- `eq:auto-6` <- ['eq:BS_PFB']
  - Falloff rates of hatted metric functions near future null infinity (eq:falloff_partial_inertial)
- `eq:true_BS` <- ['eq:auto-6']
  - Extra falloff on beta-hat promoting partially flat to true Bondi-Sachs coordinates
- `eq:null_wt` <- ['eq:3+1_metric', 'eq:def_s']
  - Outgoing null generator at the worldtube with affine parameter lambda, built from Cauchy normals
- `eq:BS_bondi_like` <- ['eq:BS_angular_gauge_condition', 'eq:null_wt']
  - Bondi-Sachs metric in (unhatted) Bondi-like coordinates without asymptotic-flatness requirements
- `eq:auto-7` <- ['eq:BS_bondi_like']
  - Relaxed O(r^0) falloff rates of the unhatted Bondi-like metric functions
- `eq:UA_asymptotic` <- ['eq:BS_bondi_like']
  - Definition of U^(0)A, the asymptotic value of the Bondi shift function U^A
- `q_superA_expression` <- []  (labels: q_superA_expression, q_subA_expression)
  - Complex dyad on the unit sphere, contravariant and covariant forms (eq:q_hat_two_expressions)
- `eq:dyad_property` <- ['q_superA_expression']
  - Dyad normalization identities q^A q_A = 0, q^A qbar_A = 2
- `eq:CCE_all_variables:K` <- ['eq:BS_bondi_like', 'q_superA_expression']
  - Spin-weighted Bondi scalars U, J, K with K = sqrt(1+J Jbar) (typo-corrected vs Moxon 2020 Eq. 10e)
- `eq:CCE_tetrad_m` <- ['eq:BS_bondi_like', 'eq:CCE_all_variables:K', 'q_superA_expression']  (labels: eq:CCE_tetrad_m, eq:CCE_tetrad_k, eq:CCE_tetrad_l)
  - CCE null tetrad (m, k, l) in Bondi-like coordinates; same form applies to hatted variables
- `eq:psi0_CCE` <- ['eq:CCE_tetrad_m', 'eq:CCE_all_variables:K']
  - Closed-form Bondi-like psi0 from radial derivatives of beta, J, K; fed back to the Cauchy boundary
## Sec. IV.A (subsec:Jacobians_for_ccm)

- `eq:cauchy_null_radius_dependence` <- ['eq:null_wt']
  - Definition of null-radius coordinates: Cauchy time/angles, radial coordinate replaced by affine parameter
- `eq:transformation_cauchy_null_radius` <- ['eq:cauchy_null_radius_dependence', 'eq:3+1_metric']
  - Metric components in null-radius coordinates: null gauge plus inherited Cauchy components
- `eq:Jacobian_cauchy_null_radius` <- ['eq:transformation_cauchy_null_radius']
  - Jacobian from null-radius to Cauchy coordinates
- `eq:auto-8` <- ['eq:BS_angular_gauge_condition', 'eq:transformation_cauchy_null_radius']
  - Bondi-like areal radius from the fourth root of the angular metric determinant ratio
- `eq:transformation_null_radius_bondi_like` <- ['eq:auto-8', 'eq:cauchy_null_radius_dependence']
  - Null-radius to Bondi-like coordinate map (only the radial coordinate changes)
- `eq:Jacobian_null_radius_bondi_like` <- ['eq:transformation_null_radius_bondi_like']
  - Jacobian from null-radius to Bondi-like coordinates
- `eq:Jacobian_bondi_like_null_radius` <- ['eq:Jacobian_null_radius_bondi_like']
  - Inverse Jacobian: Bondi-like to null-radius coordinates
- `eq:du_xhat` <- ['eq:UA_asymptotic']
  - Evolution equation for partially flat angular coordinates removing the asymptotic part of U^A
- `eq:ahat` <- ['q_superA_expression']  (labels: eq:ahat, eq:bhat)
  - Spin-weighted Jacobian factors a-hat, b-hat of the angular map contracted on the dyads
- `eq:omegahat` <- ['eq:ahat']
  - Conformal factor omega-hat of the angular map from a-hat, b-hat
- `eq:transformation_bondi_like_inertial` <- ['eq:du_xhat', 'eq:omegahat', 'eq:BS_angular_gauge_condition']
  - Bondi-like to partially flat Bondi-like map: angular relabeling plus conformal radius rescaling
- `eq:expand_angular_jacobian_nohat` <- ['eq:ahat', 'eq:omegahat', 'eq:dyad_property']
  - Dyad expansion of the angular Jacobian d x^A / d xhat^A; middle determinant equals -omegahat^2
- `eq:auto-9` <- ['q_superA_expression']
  - Spin-weighted factors a, b of the inverse angular Jacobian
- `eq:oemga_b_a` <- ['eq:auto-9']
  - Conformal factor omega associated with the inverse angular map
- `eq:expand_angular_jacobian` <- ['eq:auto-9', 'eq:oemga_b_a', 'eq:dyad_property']
  - Dyad expansion of the inverse angular Jacobian d xhat^A / d x^A
- `eq:a_ahat_b_bhat` <- ['eq:expand_angular_jacobian_nohat', 'eq:expand_angular_jacobian']
  - Relations a = -ahat/omegahat^2, b = bhatbar/omegahat^2 from Jacobian invertibility
- `eq:auto-10` <- ['eq:a_ahat_b_bhat', 'eq:oemga_b_a', 'eq:omegahat']
  - Identity omega * omegahat = 1: the conformal factors are mutually inverse
- `eq:Jacobian_bondi_like_inertial` <- ['eq:transformation_bondi_like_inertial']
  - Jacobian from Bondi-like to partially flat Bondi-like coordinates
- `eq:jacobian_bondi_inertial` <- ['eq:Jacobian_bondi_like_inertial', 'eq:du_xhat', 'eq:auto-10']
  - Inverse Jacobian (partially flat to Bondi-like), simplified using eq:du_xhat
## Sec. IV.B (sec:tetrad_transformation_Scenario_I)

- `eq:auto-11` <- []
  - Imported relation beta-hat = -1/2 ln(d rhat / d lambda) [Moxon 2020 Eqs. 19a, 33a]
- `eq:null_vector_CCE_rhat` <- ['eq:CCE_tetrad_m', 'eq:Jacobian_bondi_like_null_radius', 'eq:jacobian_bondi_inertial', 'eq:auto-11']
  - Characteristic outgoing null vector l-hat expressed via the null-radius generator; shows l-hat proportional to Cauchy l
- `eq:l_GH_and_l_CCE_transform` <- ['eq:null_vector_CCE_rhat', 'eq:null_wt', 'eq:gh_tetrad_l']
  - Explicit boost factor between Cauchy and characteristic outgoing null vectors
- `eq:lorentz_transformation_I` <- []
  - Type I Lorentz transformation (null rotation about l) with spin-weight-1 parameter kappa
- `eq:lorentz_transformation_II` <- []
  - Type II Lorentz transformation (boost A + spin rotation Theta)
- `eq:auto-12` <- ['eq:lorentz_transformation_I', 'eq:GH_psi0_def']
  - psi0-hat invariance under Type I transformations (no mixing of Weyl scalars)
- `eq:lorentz_psi0_ii` <- ['eq:lorentz_transformation_II', 'eq:GH_psi0_def']
  - psi0-hat rescaling A^2 e^{2i Theta} under Type II transformations
- `typeI_Ahat_to_Ap` <- ['eq:Jacobian_cauchy_null_radius', 'eq:Jacobian_bondi_like_null_radius', 'eq:jacobian_bondi_inertial']
  - Hatted angular basis in Cauchy angular basis plus an l-proportional, Type-I-removable term
- `eq:q_transformation_type_i` <- ['typeI_Ahat_to_Ap', 'eq:expand_angular_jacobian_nohat', 'eq:null_vector_CCE_rhat', 'eq:lorentz_transformation_I', 'q_superA_expression']
  - Hatted dyad in Cauchy dyad after dropping l-terms via Type I transformation
- `eq:auto-13` <- ['eq:q_transformation_type_i', 'eq:CCE_tetrad_m']
  - Type-I-equivalent expression of m-hat in the Cauchy angular dyad basis
- `eq:new_m_cce_to_GH_abstract` <- ['eq:auto-13']
  - m-hat equals a Cauchy-angular m' up to a Type I transformation
- `eq:auto-14` <- ['eq:auto-13']
  - Component coefficients M-hat_theta', M-hat_phi' of the Choice-1 m vector
- `eq:new_m_cce_to_GH` <- ['eq:new_m_cce_to_GH_abstract', 'eq:auto-14']
  - Choice-1 Cauchy tetrad vector m in the Cauchy angular subspace; satisfies eq:GH_m, used in eq:wab
- `eq:auto-15` <- []
  - Cartesian components of the Cauchy angular basis vectors at the worldtube radius
- `eq:auto-16` <- ['eq:l_GH_and_l_CCE_transform', 'eq:lorentz_transformation_II']
  - Type II boost parameter A-hat relating l-hat and Cauchy l for Choice 1
- `eq:psi0_prime_psi0_hat` <- ['eq:auto-16', 'eq:lorentz_psi0_ii', 'eq:psi0_CCE']
  - Choice-1 result: psi0' = A-hat^2 psi0-hat; phase freedom cancels in w-
## Sec. IV.C (sec:tetrad_transformation_Scenario_II)

- `eq:auto-17` <- ['eq:Jacobian_cauchy_null_radius', 'eq:Jacobian_bondi_like_null_radius', 'eq:lorentz_transformation_I']
  - Bondi-like angular basis equals Cauchy angular basis up to a Type-I-removable term
- `eq:new_m_cce_to_GH_B_abstract` <- ['eq:auto-17', 'q_superA_expression']
  - Bondi-like dyad Type-I-equivalent to Cauchy dyad
- `eq:auto-18` <- ['eq:new_m_cce_to_GH_B_abstract', 'eq:CCE_tetrad_m']
  - Bondi-like m Type-I-equivalent to a Cauchy-angular m'
- `eq:auto-19` <- ['eq:auto-18']
  - Component coefficients M_theta', M_phi' of the Choice-2 m vector
- `eq:new_m_cce_to_GH_B` <- ['eq:auto-18', 'eq:auto-19']
  - Choice-2 Cauchy tetrad vector m in the Cauchy angular subspace; satisfies eq:GH_m, used in eq:wab
- `eq:J_BS_like` <- ['eq:CCE_all_variables:K', 'eq:auto-9', 'eq:oemga_b_a']  (labels: eq:J_BS_like, eq:K_BS_like)
  - Transformation of evolved hatted J-hat, K-hat back to Bondi-like J, K
- `eq:auto-20` <- ['eq:l_GH_and_l_CCE_transform', 'eq:lorentz_transformation_II']
  - Type II boost parameter A relating Bondi-like l and Cauchy l for Choice 2
- `eq:auto-21` <- ['eq:auto-20', 'eq:lorentz_psi0_ii', 'eq:psi0_CCE', 'eq:J_BS_like']
  - Choice-2 result: psi0' = A^2 psi0
## Sec. IV.D (sec:ccm_Cauchy_inertal)

- `eq:duhat_x` <- ['eq:jacobian_bondi_inertial']
  - Evolution equation for the inverse angular map x^A(uhat, xhat) read off the inverse Jacobian
- `eq:auto-22` <- []
  - Cartesian unit-sphere embedding of the Bondi-like angular coordinates (spin-weight 0)
- `eq:auto-23` <- ['eq:auto-22', 'q_superA_expression']
  - Spin-weighted eth derivatives of the Cartesian sphere coordinates
- `eq:auto-24` <- ['eq:duhat_x']
  - Definition of the auxiliary hatted-frame advection variable U^(0) (mathcal U)
- `eq:duhat_x_cartesian` <- ['eq:duhat_x', 'eq:auto-23', 'eq:auto-24']
  - Cartesian spin-weighted form of eq:duhat_x; the equation evolved numerically for the inverse map
- `eq:U_mathcal_U` <- ['eq:auto-24', 'eq:a_ahat_b_bhat']
  - Explicit mathcal-U^(0) in terms of U^(0) and the Jacobian factors
- `eq:auto-25` <- ['eq:U_mathcal_U', 'eq:oemga_b_a']
  - Inverse relation: U^(0) from mathcal-U^(0) using inverse-map factors
- `eq:du_xhat_cartesian` <- ['eq:du_xhat', 'eq:auto-23']
  - Cartesian version of eq:du_xhat for the forward angular map
- `eq:constraints_ca_cb` <- ['eq:a_ahat_b_bhat']
  - Proposed CCM consistency constraints C_a, C_b monitoring the inverse-Jacobian relation (commented-out in tex)
## Sec. V (sec:tests)

- `eq:auto-27` <- []
  - The (l=2, m=0) spherical harmonic Y20
- `eq:Gaussian_pulse_F` <- []
  - Outgoing Gaussian pulse profile F(u') with amplitude X, center r_c, width tau
- `eq:Teukolsky_wave_outgoing` <- ['eq:Gaussian_pulse_F']
  - Definition of the n-th retarded-time derivative F^(n)
- `eq:auto-26` <- ['eq:Gaussian_pulse_F', 'eq:Teukolsky_wave_outgoing']
  - Radial amplitude functions A, B, C of the Teukolsky wave (eq:Teukolsky_ABC)
- `eq:Teukolsky_angular_basis` <- ['eq:auto-27']
  - Angular basis functions of the l=2, m=0 even-parity Teukolsky-wave metric
- `eq:Teukolsky_metric` <- ['eq:auto-26', 'eq:Teukolsky_angular_basis', 'eq:auto-27']
  - Perturbative Teukolsky-wave metric on a flat background (analytic test reference)
- `eq:gauge_constraint` <- []
  - L2 norm of the GH gauge constraint C_a (diagnostic; imported from Lindblom 2005)
- `eq:three_constraint` <- []
  - L2 norm of the GH three-index constraint C_iab (diagnostic; imported from Lindblom 2005)
- `eq:Teukolsky_h_20` <- ['eq:Teukolsky_metric', 'eq:Teukolsky_wave_outgoing']
  - Perturbative prediction for the (2,0) strain harmonic at scri
- `eq:bondi_violation_psi3` <- []
  - Bondi-gauge constraints from NP Bianchi identities (eq:bondi_violation; waveform-quality diagnostics)
- `eq:im_psi2` <- []
  - Reality condition on the Bondi mass aspect (Im psi2 constraint)
- `eq:auto-28` <- ['eq:Gaussian_pulse_F']
  - Ingoing Gaussian profile F(v') in advanced time for the Kerr-perturbation test
- `eq:auto-30` <- []
  - Radial profile mathcal-J of the injected characteristic pulse (envelope times Gaussian)
- `eq:auto-29` <- ['eq:auto-30']
  - Compactly supported spin-weight-2 initial data for J-hat on the initial null slice
## Appendix comment-block (tex \begin{comment})

- `eq:auto-31` <- ['eq:3+1_metric']
  - GH gauge condition: coordinates obey a wave equation with gauge source H (commented-out in tex)
- `eq:auto-32` <- ['FOSH', 'eq:auto-1']
  - FOSH system projected onto characteristic fields, normal/tangential split (commented-out in tex)
- `eq:auto-33` <- ['Bjorhus_bc']
  - Bulk reduction of the Bjorhus condition d_t u = D_t u (commented-out in tex)
- `eq:auto-34` <- ['FOSH']
  - GH three-index constraint field c^3 (commented-out in tex)
- `eq:auto-35` <- ['FOSH', 'eq:def_s']
  - Incoming characteristic field u^1- of the GH system (commented-out in tex)
- `eq:CCE_all_variables:H` <- ['eq:BS_PFB', 'q_superA_expression']  (labels: eq:CCE_all_variables:H, eq:CCE_all_variables:K)
  - Hatted spin-weighted characteristic variables U-hat, Q-hat, J-hat, K-hat, H-hat = du J-hat (commented-out in tex)
- `eq:cce_com:H` <- ['eq:BS_PFB', 'eq:CCE_all_variables:H']
  - Hierarchical radial hypersurface equations of the characteristic system (eq:cce_com; commented-out in tex)
- `eq:auto-36` <- ['q_superA_expression']
  - Unit-sphere metric expressed through the complex dyad (commented-out in tex)
## Appendix A (app:1e-5)

- `eq:auto-38` <- []
  - Spin-weighted spherical harmonics sY20 used in the analytic Teukolsky expressions
- `eq:auto-37` <- ['eq:Teukolsky_metric', 'eq:auto-26', 'eq:auto-38']
  - Analytic perturbative Weyl scalars, News, and strain of the Teukolsky wave (eq:Teukolsky_weyl_strain_news)
- `eq:Teukolsky_bulk_psi0` <- ['eq:auto-37', 'eq:auto-26']
  - Scri simplifications (eq:Teukolsky_scri_weyl_strain_news); psi0 ~ r^-5 valid in the bulk

Totals: 106 nodes; 169 edges; acyclic (verified by scripts/eqdag_check.py).