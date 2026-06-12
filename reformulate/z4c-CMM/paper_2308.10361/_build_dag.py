#!/usr/bin/env python3
"""Single source of truth for the arXiv:2308.10361 equation DAG.

Generates: logic.md and knowledge-database rows (via the validating append CLI).
Node tuple: (node_id, [equation_labels], [predecessors], section_anchor, summary)
Dependencies mirror reformulate/z4c-CMM/paper_2308.10361/derivation.md exactly.
"""
import json, subprocess, sys
from pathlib import Path

ROOT = Path("/data/haiyangw/claude/z4c-CMM")

N = []  # (node_id, labels, preds, anchor, summary)
def n(nid, labels, preds, anchor, summary):
    N.append((nid, labels, preds, anchor, summary))

S2 = "Sec. II (sec:GH)"
n("eq:3+1_metric", ["eq:3+1_metric"], [], S2, "3+1 (ADM) decomposition of the spacetime metric in Cauchy (primed) coordinates: lapse, shift, spatial metric")
n("FOSH", ["FOSH"], ["eq:3+1_metric"], S2, "Vacuum Einstein equations as first-order symmetric hyperbolic generalized-harmonic system for {g,Pi,Phi}")
n("eq:def_s", ["eq:def_s"], ["eq:3+1_metric"], S2, "Outward-directed spatial unit normal to the Cauchy outer boundary")
n("eq:auto-1", ["eq:auto-1"], ["FOSH", "eq:def_s"], S2, "Left eigenvectors and characteristic speeds of the FOSH principal symbol along the boundary normal")
n("eq:auto-2", ["eq:auto-2"], ["eq:auto-1", "eq:def_s"], S2, "Definitions of characteristic-projected time derivative d_t and normal derivative d_perp")
n("eq:auto-3", ["eq:auto-3"], ["FOSH", "eq:auto-1"], S2, "Definition of D_t u: characteristic projection of the full FOSH right-hand side")
n("Bjorhus_bc", ["Bjorhus_bc"], ["eq:auto-1", "eq:auto-2", "eq:auto-3"], S2, "Bjorhus boundary method: replace normal derivative of incoming characteristic fields by desired boundary value")
n("eq:projection_operator", ["eq:projection_operator"], ["eq:3+1_metric", "eq:def_s"], S2, "2D projector orthogonal to hypersurface normal n and boundary normal s")
n("eq:auto-4", ["eq:auto-4"], ["eq:projection_operator"], S2, "Physical (transverse-traceless) projection operator built from the 2D projector")
n("eq:wab_projection", ["eq:wab_projection"], ["eq:auto-4", "eq:projection_operator"], S2, "Inward-propagating Weyl components w-: TT-projected double contraction of Weyl tensor with n+s")
n("eq:bc_bjorhus", ["eq:bc_bjorhus"], ["Bjorhus_bc", "eq:auto-4", "eq:wab_projection"], S2, "Physical boundary condition on u^1- driving w- toward its CCM-supplied boundary value")
n("eq:gh_tetrad_l", ["eq:gh_tetrad_l", "eq:gh_tetrad_k"], ["eq:def_s"], S2, "Unique Cauchy outgoing/ingoing null tetrad vectors l,k from time normal n and boundary normal s")
n("eq:GH_m", ["eq:GH_m"], ["eq:gh_tetrad_l"], S2, "Orthonormality conditions on the Cauchy angular tetrad vector m (fixed up to phase rotation)")
n("eq:GH_psi0_def", ["eq:GH_psi0_def"], ["eq:gh_tetrad_l", "eq:GH_m"], S2, "Definition of the Cauchy Weyl scalar psi0' in the GH tetrad (ingoing transverse radiation)")
n("eq:auto-5", ["eq:auto-5"], ["eq:projection_operator", "eq:GH_m"], S2, "Completeness identity P = m mbar + mbar m for the 2D projector")
n("eq:wab", ["eq:wab"], ["eq:wab_projection", "eq:auto-5", "eq:GH_psi0_def"], S2, "w- rewritten in terms of psi0' and m; the slot where the characteristic psi0 is inserted for CCM")

S3 = "Sec. III (sec:CCE)"
n("eq:BS_PFB", ["eq:BS_PFB"], [], S3, "Bondi-Sachs metric in partially flat Bondi-like (hatted) coordinates used by SpECTRE CCE")
n("eq:BS_angular_gauge_condition", ["eq:BS_angular_gauge_condition"], ["eq:BS_PFB"], S3, "Bondi determinant gauge condition: det(h_AB) equals unit-sphere determinant")
n("eq:auto-6", ["eq:auto-6"], ["eq:BS_PFB"], S3, "Falloff rates of hatted metric functions near future null infinity (eq:falloff_partial_inertial)")
n("eq:true_BS", ["eq:true_BS"], ["eq:auto-6"], S3, "Extra falloff on beta-hat promoting partially flat to true Bondi-Sachs coordinates")
n("eq:null_wt", ["eq:null_wt"], ["eq:3+1_metric", "eq:def_s"], S3, "Outgoing null generator at the worldtube with affine parameter lambda, built from Cauchy normals")
n("eq:BS_bondi_like", ["eq:BS_bondi_like"], ["eq:BS_angular_gauge_condition", "eq:null_wt"], S3, "Bondi-Sachs metric in (unhatted) Bondi-like coordinates without asymptotic-flatness requirements")
n("eq:auto-7", ["eq:auto-7"], ["eq:BS_bondi_like"], S3, "Relaxed O(r^0) falloff rates of the unhatted Bondi-like metric functions")
n("eq:UA_asymptotic", ["eq:UA_asymptotic"], ["eq:BS_bondi_like"], S3, "Definition of U^(0)A, the asymptotic value of the Bondi shift function U^A")
n("q_superA_expression", ["q_superA_expression", "q_subA_expression"], [], S3, "Complex dyad on the unit sphere, contravariant and covariant forms (eq:q_hat_two_expressions)")
n("eq:dyad_property", ["eq:dyad_property"], ["q_superA_expression"], S3, "Dyad normalization identities q^A q_A = 0, q^A qbar_A = 2")
n("eq:CCE_all_variables:K", ["eq:CCE_all_variables:K"], ["eq:BS_bondi_like", "q_superA_expression"], S3, "Spin-weighted Bondi scalars U, J, K with K = sqrt(1+J Jbar) (typo-corrected vs Moxon 2020 Eq. 10e)")
n("eq:CCE_tetrad_m", ["eq:CCE_tetrad_m", "eq:CCE_tetrad_k", "eq:CCE_tetrad_l"], ["eq:BS_bondi_like", "eq:CCE_all_variables:K", "q_superA_expression"], S3, "CCE null tetrad (m, k, l) in Bondi-like coordinates; same form applies to hatted variables")
n("eq:psi0_CCE", ["eq:psi0_CCE"], ["eq:CCE_tetrad_m", "eq:CCE_all_variables:K"], S3, "Closed-form Bondi-like psi0 from radial derivatives of beta, J, K; fed back to the Cauchy boundary")

S4A = "Sec. IV.A (subsec:Jacobians_for_ccm)"
n("eq:cauchy_null_radius_dependence", ["eq:cauchy_null_radius_dependence"], ["eq:null_wt"], S4A, "Definition of null-radius coordinates: Cauchy time/angles, radial coordinate replaced by affine parameter")
n("eq:transformation_cauchy_null_radius", ["eq:transformation_cauchy_null_radius"], ["eq:cauchy_null_radius_dependence", "eq:3+1_metric"], S4A, "Metric components in null-radius coordinates: null gauge plus inherited Cauchy components")
n("eq:Jacobian_cauchy_null_radius", ["eq:Jacobian_cauchy_null_radius"], ["eq:transformation_cauchy_null_radius"], S4A, "Jacobian from null-radius to Cauchy coordinates")
n("eq:auto-8", ["eq:auto-8"], ["eq:BS_angular_gauge_condition", "eq:transformation_cauchy_null_radius"], S4A, "Bondi-like areal radius from the fourth root of the angular metric determinant ratio")
n("eq:transformation_null_radius_bondi_like", ["eq:transformation_null_radius_bondi_like"], ["eq:auto-8", "eq:cauchy_null_radius_dependence"], S4A, "Null-radius to Bondi-like coordinate map (only the radial coordinate changes)")
n("eq:Jacobian_null_radius_bondi_like", ["eq:Jacobian_null_radius_bondi_like"], ["eq:transformation_null_radius_bondi_like"], S4A, "Jacobian from null-radius to Bondi-like coordinates")
n("eq:Jacobian_bondi_like_null_radius", ["eq:Jacobian_bondi_like_null_radius"], ["eq:Jacobian_null_radius_bondi_like"], S4A, "Inverse Jacobian: Bondi-like to null-radius coordinates")
n("eq:du_xhat", ["eq:du_xhat"], ["eq:UA_asymptotic"], S4A, "Evolution equation for partially flat angular coordinates removing the asymptotic part of U^A")
n("eq:ahat", ["eq:ahat", "eq:bhat"], ["q_superA_expression"], S4A, "Spin-weighted Jacobian factors a-hat, b-hat of the angular map contracted on the dyads")
n("eq:omegahat", ["eq:omegahat"], ["eq:ahat"], S4A, "Conformal factor omega-hat of the angular map from a-hat, b-hat")
n("eq:transformation_bondi_like_inertial", ["eq:transformation_bondi_like_inertial"], ["eq:du_xhat", "eq:omegahat", "eq:BS_angular_gauge_condition"], S4A, "Bondi-like to partially flat Bondi-like map: angular relabeling plus conformal radius rescaling")
n("eq:expand_angular_jacobian_nohat", ["eq:expand_angular_jacobian_nohat"], ["eq:ahat", "eq:omegahat", "eq:dyad_property"], S4A, "Dyad expansion of the angular Jacobian d x^A / d xhat^A; middle determinant equals -omegahat^2")
n("eq:auto-9", ["eq:auto-9"], ["q_superA_expression"], S4A, "Spin-weighted factors a, b of the inverse angular Jacobian")
n("eq:oemga_b_a", ["eq:oemga_b_a"], ["eq:auto-9"], S4A, "Conformal factor omega associated with the inverse angular map")
n("eq:expand_angular_jacobian", ["eq:expand_angular_jacobian"], ["eq:auto-9", "eq:oemga_b_a", "eq:dyad_property"], S4A, "Dyad expansion of the inverse angular Jacobian d xhat^A / d x^A")
n("eq:a_ahat_b_bhat", ["eq:a_ahat_b_bhat"], ["eq:expand_angular_jacobian_nohat", "eq:expand_angular_jacobian"], S4A, "Relations a = -ahat/omegahat^2, b = bhatbar/omegahat^2 from Jacobian invertibility")
n("eq:auto-10", ["eq:auto-10"], ["eq:a_ahat_b_bhat", "eq:oemga_b_a", "eq:omegahat"], S4A, "Identity omega * omegahat = 1: the conformal factors are mutually inverse")
n("eq:Jacobian_bondi_like_inertial", ["eq:Jacobian_bondi_like_inertial"], ["eq:transformation_bondi_like_inertial"], S4A, "Jacobian from Bondi-like to partially flat Bondi-like coordinates")
n("eq:jacobian_bondi_inertial", ["eq:jacobian_bondi_inertial"], ["eq:Jacobian_bondi_like_inertial", "eq:du_xhat", "eq:auto-10"], S4A, "Inverse Jacobian (partially flat to Bondi-like), simplified using eq:du_xhat")

S4B = "Sec. IV.B (sec:tetrad_transformation_Scenario_I)"
n("eq:auto-11", ["eq:auto-11"], [], S4B, "Imported relation beta-hat = -1/2 ln(d rhat / d lambda) [Moxon 2020 Eqs. 19a, 33a]")
n("eq:null_vector_CCE_rhat", ["eq:null_vector_CCE_rhat"], ["eq:CCE_tetrad_m", "eq:Jacobian_bondi_like_null_radius", "eq:jacobian_bondi_inertial", "eq:auto-11"], S4B, "Characteristic outgoing null vector l-hat expressed via the null-radius generator; shows l-hat proportional to Cauchy l")
n("eq:l_GH_and_l_CCE_transform", ["eq:l_GH_and_l_CCE_transform"], ["eq:null_vector_CCE_rhat", "eq:null_wt", "eq:gh_tetrad_l"], S4B, "Explicit boost factor between Cauchy and characteristic outgoing null vectors")
n("eq:lorentz_transformation_I", ["eq:lorentz_transformation_I"], [], S4B, "Type I Lorentz transformation (null rotation about l) with spin-weight-1 parameter kappa")
n("eq:lorentz_transformation_II", ["eq:lorentz_transformation_II"], [], S4B, "Type II Lorentz transformation (boost A + spin rotation Theta)")
n("eq:auto-12", ["eq:auto-12"], ["eq:lorentz_transformation_I", "eq:GH_psi0_def"], S4B, "psi0-hat invariance under Type I transformations (no mixing of Weyl scalars)")
n("eq:lorentz_psi0_ii", ["eq:lorentz_psi0_ii"], ["eq:lorentz_transformation_II", "eq:GH_psi0_def"], S4B, "psi0-hat rescaling A^2 e^{2i Theta} under Type II transformations")
n("typeI_Ahat_to_Ap", ["typeI_Ahat_to_Ap"], ["eq:Jacobian_cauchy_null_radius", "eq:Jacobian_bondi_like_null_radius", "eq:jacobian_bondi_inertial"], S4B, "Hatted angular basis in Cauchy angular basis plus an l-proportional, Type-I-removable term")
n("eq:q_transformation_type_i", ["eq:q_transformation_type_i"], ["typeI_Ahat_to_Ap", "eq:expand_angular_jacobian_nohat", "eq:null_vector_CCE_rhat", "eq:lorentz_transformation_I", "q_superA_expression"], S4B, "Hatted dyad in Cauchy dyad after dropping l-terms via Type I transformation")
n("eq:auto-13", ["eq:auto-13"], ["eq:q_transformation_type_i", "eq:CCE_tetrad_m"], S4B, "Type-I-equivalent expression of m-hat in the Cauchy angular dyad basis")
n("eq:new_m_cce_to_GH_abstract", ["eq:new_m_cce_to_GH_abstract"], ["eq:auto-13"], S4B, "m-hat equals a Cauchy-angular m' up to a Type I transformation")
n("eq:auto-14", ["eq:auto-14"], ["eq:auto-13"], S4B, "Component coefficients M-hat_theta', M-hat_phi' of the Choice-1 m vector")
n("eq:new_m_cce_to_GH", ["eq:new_m_cce_to_GH"], ["eq:new_m_cce_to_GH_abstract", "eq:auto-14"], S4B, "Choice-1 Cauchy tetrad vector m in the Cauchy angular subspace; satisfies eq:GH_m, used in eq:wab")
n("eq:auto-15", ["eq:auto-15"], [], S4B, "Cartesian components of the Cauchy angular basis vectors at the worldtube radius")
n("eq:auto-16", ["eq:auto-16"], ["eq:l_GH_and_l_CCE_transform", "eq:lorentz_transformation_II"], S4B, "Type II boost parameter A-hat relating l-hat and Cauchy l for Choice 1")
n("eq:psi0_prime_psi0_hat", ["eq:psi0_prime_psi0_hat"], ["eq:auto-16", "eq:lorentz_psi0_ii", "eq:psi0_CCE"], S4B, "Choice-1 result: psi0' = A-hat^2 psi0-hat; phase freedom cancels in w-")

S4C = "Sec. IV.C (sec:tetrad_transformation_Scenario_II)"
n("eq:auto-17", ["eq:auto-17"], ["eq:Jacobian_cauchy_null_radius", "eq:Jacobian_bondi_like_null_radius", "eq:lorentz_transformation_I"], S4C, "Bondi-like angular basis equals Cauchy angular basis up to a Type-I-removable term")
n("eq:new_m_cce_to_GH_B_abstract", ["eq:new_m_cce_to_GH_B_abstract"], ["eq:auto-17", "q_superA_expression"], S4C, "Bondi-like dyad Type-I-equivalent to Cauchy dyad")
n("eq:auto-18", ["eq:auto-18"], ["eq:new_m_cce_to_GH_B_abstract", "eq:CCE_tetrad_m"], S4C, "Bondi-like m Type-I-equivalent to a Cauchy-angular m'")
n("eq:auto-19", ["eq:auto-19"], ["eq:auto-18"], S4C, "Component coefficients M_theta', M_phi' of the Choice-2 m vector")
n("eq:new_m_cce_to_GH_B", ["eq:new_m_cce_to_GH_B"], ["eq:auto-18", "eq:auto-19"], S4C, "Choice-2 Cauchy tetrad vector m in the Cauchy angular subspace; satisfies eq:GH_m, used in eq:wab")
n("eq:J_BS_like", ["eq:J_BS_like", "eq:K_BS_like"], ["eq:CCE_all_variables:K", "eq:auto-9", "eq:oemga_b_a"], S4C, "Transformation of evolved hatted J-hat, K-hat back to Bondi-like J, K")
n("eq:auto-20", ["eq:auto-20"], ["eq:l_GH_and_l_CCE_transform", "eq:lorentz_transformation_II"], S4C, "Type II boost parameter A relating Bondi-like l and Cauchy l for Choice 2")
n("eq:auto-21", ["eq:auto-21"], ["eq:auto-20", "eq:lorentz_psi0_ii", "eq:psi0_CCE", "eq:J_BS_like"], S4C, "Choice-2 result: psi0' = A^2 psi0")

S4D = "Sec. IV.D (sec:ccm_Cauchy_inertal)"
n("eq:duhat_x", ["eq:duhat_x"], ["eq:jacobian_bondi_inertial"], S4D, "Evolution equation for the inverse angular map x^A(uhat, xhat) read off the inverse Jacobian")
n("eq:auto-22", ["eq:auto-22"], [], S4D, "Cartesian unit-sphere embedding of the Bondi-like angular coordinates (spin-weight 0)")
n("eq:auto-23", ["eq:auto-23"], ["eq:auto-22", "q_superA_expression"], S4D, "Spin-weighted eth derivatives of the Cartesian sphere coordinates")
n("eq:auto-24", ["eq:auto-24"], ["eq:duhat_x"], S4D, "Definition of the auxiliary hatted-frame advection variable U^(0) (mathcal U)")
n("eq:duhat_x_cartesian", ["eq:duhat_x_cartesian"], ["eq:duhat_x", "eq:auto-23", "eq:auto-24"], S4D, "Cartesian spin-weighted form of eq:duhat_x; the equation evolved numerically for the inverse map")
n("eq:U_mathcal_U", ["eq:U_mathcal_U"], ["eq:auto-24", "eq:a_ahat_b_bhat"], S4D, "Explicit mathcal-U^(0) in terms of U^(0) and the Jacobian factors")
n("eq:auto-25", ["eq:auto-25"], ["eq:U_mathcal_U", "eq:oemga_b_a"], S4D, "Inverse relation: U^(0) from mathcal-U^(0) using inverse-map factors")
n("eq:du_xhat_cartesian", ["eq:du_xhat_cartesian"], ["eq:du_xhat", "eq:auto-23"], S4D, "Cartesian version of eq:du_xhat for the forward angular map")
n("eq:constraints_ca_cb", ["eq:constraints_ca_cb"], ["eq:a_ahat_b_bhat"], S4D, "Proposed CCM consistency constraints C_a, C_b monitoring the inverse-Jacobian relation (commented-out in tex)")

S5 = "Sec. V (sec:tests)"
n("eq:auto-27", ["eq:auto-27"], [], S5, "The (l=2, m=0) spherical harmonic Y20")
n("eq:Gaussian_pulse_F", ["eq:Gaussian_pulse_F"], [], S5, "Outgoing Gaussian pulse profile F(u') with amplitude X, center r_c, width tau")
n("eq:Teukolsky_wave_outgoing", ["eq:Teukolsky_wave_outgoing"], ["eq:Gaussian_pulse_F"], S5, "Definition of the n-th retarded-time derivative F^(n)")
n("eq:auto-26", ["eq:auto-26"], ["eq:Gaussian_pulse_F", "eq:Teukolsky_wave_outgoing"], S5, "Radial amplitude functions A, B, C of the Teukolsky wave (eq:Teukolsky_ABC)")
n("eq:Teukolsky_angular_basis", ["eq:Teukolsky_angular_basis"], ["eq:auto-27"], S5, "Angular basis functions of the l=2, m=0 even-parity Teukolsky-wave metric")
n("eq:Teukolsky_metric", ["eq:Teukolsky_metric"], ["eq:auto-26", "eq:Teukolsky_angular_basis", "eq:auto-27"], S5, "Perturbative Teukolsky-wave metric on a flat background (analytic test reference)")
n("eq:gauge_constraint", ["eq:gauge_constraint"], [], S5, "L2 norm of the GH gauge constraint C_a (diagnostic; imported from Lindblom 2005)")
n("eq:three_constraint", ["eq:three_constraint"], [], S5, "L2 norm of the GH three-index constraint C_iab (diagnostic; imported from Lindblom 2005)")
n("eq:Teukolsky_h_20", ["eq:Teukolsky_h_20"], ["eq:Teukolsky_metric", "eq:Teukolsky_wave_outgoing"], S5, "Perturbative prediction for the (2,0) strain harmonic at scri")
n("eq:bondi_violation_psi3", ["eq:bondi_violation_psi3"], [], S5, "Bondi-gauge constraints from NP Bianchi identities (eq:bondi_violation; waveform-quality diagnostics)")
n("eq:im_psi2", ["eq:im_psi2"], [], S5, "Reality condition on the Bondi mass aspect (Im psi2 constraint)")
n("eq:auto-28", ["eq:auto-28"], ["eq:Gaussian_pulse_F"], S5, "Ingoing Gaussian profile F(v') in advanced time for the Kerr-perturbation test")
n("eq:auto-30", ["eq:auto-30"], [], S5, "Radial profile mathcal-J of the injected characteristic pulse (envelope times Gaussian)")
n("eq:auto-29", ["eq:auto-29"], ["eq:auto-30"], S5, "Compactly supported spin-weight-2 initial data for J-hat on the initial null slice")

SC = "Appendix comment-block (tex \\begin{comment})"
n("eq:auto-31", ["eq:auto-31"], ["eq:3+1_metric"], SC, "GH gauge condition: coordinates obey a wave equation with gauge source H (commented-out in tex)")
n("eq:auto-32", ["eq:auto-32"], ["FOSH", "eq:auto-1"], SC, "FOSH system projected onto characteristic fields, normal/tangential split (commented-out in tex)")
n("eq:auto-33", ["eq:auto-33"], ["Bjorhus_bc"], SC, "Bulk reduction of the Bjorhus condition d_t u = D_t u (commented-out in tex)")
n("eq:auto-34", ["eq:auto-34"], ["FOSH"], SC, "GH three-index constraint field c^3 (commented-out in tex)")
n("eq:auto-35", ["eq:auto-35"], ["FOSH", "eq:def_s"], SC, "Incoming characteristic field u^1- of the GH system (commented-out in tex)")
n("eq:CCE_all_variables:H", ["eq:CCE_all_variables:H", "eq:CCE_all_variables:K"], ["eq:BS_PFB", "q_superA_expression"], SC, "Hatted spin-weighted characteristic variables U-hat, Q-hat, J-hat, K-hat, H-hat = du J-hat (commented-out in tex)")
n("eq:cce_com:H", ["eq:cce_com:H"], ["eq:BS_PFB", "eq:CCE_all_variables:H"], SC, "Hierarchical radial hypersurface equations of the characteristic system (eq:cce_com; commented-out in tex)")
n("eq:auto-36", ["eq:auto-36"], ["q_superA_expression"], SC, "Unit-sphere metric expressed through the complex dyad (commented-out in tex)")

SA = "Appendix A (app:1e-5)"
n("eq:auto-38", ["eq:auto-38"], [], SA, "Spin-weighted spherical harmonics sY20 used in the analytic Teukolsky expressions")
n("eq:auto-37", ["eq:auto-37"], ["eq:Teukolsky_metric", "eq:auto-26", "eq:auto-38"], SA, "Analytic perturbative Weyl scalars, News, and strain of the Teukolsky wave (eq:Teukolsky_weyl_strain_news)")
n("eq:Teukolsky_bulk_psi0", ["eq:Teukolsky_bulk_psi0"], ["eq:auto-37", "eq:auto-26"], SA, "Scri simplifications (eq:Teukolsky_scri_weyl_strain_news); psi0 ~ r^-5 valid in the bulk")


def check():
    ids = [x[0] for x in N]
    assert len(ids) == len(set(ids)), "duplicate node ids"
    idset = set(ids)
    for nid, labels, preds, _, _ in N:
        for p in preds:
            assert p in idset, f"dangling predecessor {p} of {nid}"
    # acyclicity (Kahn)
    preds = {nid: list(p) for nid, _, p, _, _ in N}
    indeg = {k: len(v) for k, v in preds.items()}
    succ = {k: [] for k in preds}
    for k, v in preds.items():
        for p in v:
            succ[p].append(k)
    q = [k for k, d in indeg.items() if d == 0]
    order = []
    while q:
        x = q.pop()
        order.append(x)
        for s in succ[x]:
            indeg[s] -= 1
            if indeg[s] == 0:
                q.append(s)
    assert len(order) == len(N), "cycle detected"
    labels = [l for _, ls, _, _, _ in N for l in ls]
    print(f"nodes={len(N)} labels(union)={len(set(labels))} "
          f"auto={len([l for l in set(labels) if l.startswith('eq:auto-')])} topo OK")
    return order


def gen_logic():
    out = ["# logic.md — equation dependency DAG for arXiv:2308.10361",
           "",
           "One node per equation environment (node id = equation label; unlabeled",
           "environments use eq:auto-N). `<- []` marks a root. Grouped by paper section.",
           "This file is the source of truth mirrored by",
           "`knowledge-database/paper_arxiv-2308.10361/nodes.jsonl`",
           "(task_id eq-dag-stage1-2308.10361). All edges are literature-grounded",
           "transcriptions of the dependency structure stated in the paper.",
           ""]
    cur = None
    for nid, labels, preds, anchor, summary in N:
        if anchor != cur:
            out += [f"## {anchor}", ""]
            cur = anchor
        lab = "" if labels == [nid] else f"  (labels: {', '.join(labels)})"
        out.append(f"- `{nid}` <- {preds if preds else '[]'}{lab}")
        out.append(f"  - {summary}")
    out.append("")
    out.append(f"Totals: {len(N)} nodes; "
               f"{sum(len(p) for _, _, p, _, _ in N)} edges; acyclic (verified by scripts/eqdag_check.py).")
    Path(ROOT / "reformulate/z4c-CMM/paper_2308.10361/logic.md").write_text("\n".join(out))
    print("wrote logic.md")


def append_db(order):
    by_id = {x[0]: x for x in N}
    cli = ROOT / "phys-agentic-loop/_common/knowledge_database.py"
    for nid in order:
        _, labels, preds, anchor, summary = by_id[nid]
        row = {
            "paper": "arxiv-2308.10361",
            "node_id": nid,
            "task_id": "eq-dag-stage1-2308.10361",
            "domain": "symbolic",
            "status": "solid",
            "summary": summary,
            "evidence": "ref-paper/arxiv-2308.10361/src/paper.tex + reformulate/z4c-CMM/paper_2308.10361/derivation.md",
            "equation_labels": labels,
            "predecessors": preds,
            "paper_anchor": anchor,
            "notes": "literature-grounded transcription; evidence type: citation/transcription, not independent derivation",
        }
        r = subprocess.run([sys.executable, str(cli), "append"],
                           input=json.dumps(row), capture_output=True, text=True, cwd=ROOT)
        if r.returncode != 0:
            print(f"FAIL {nid}\n{r.stdout}\n{r.stderr}")
            sys.exit(1)
    print(f"appended {len(order)} rows")


if __name__ == "__main__":
    order = check()
    gen_logic()
    if "--append" in sys.argv:
        append_db(order)
