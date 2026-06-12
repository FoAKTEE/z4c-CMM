# implementation_plan_python.md — arXiv:2308.10361

Code/work partition per logic.md node cluster. Target: Python (sympy for symbolic
checks, numpy/scipy + spherical/spin-weighted utilities for numerical checks).
Each task consumes its cluster's nodes from logic.md and emits verifier output for
claims.md/obligations.md. All tasks start from literature-grounded transcriptions;
their purpose is to upgrade evidence to symbolic derivation / numerical check.

## T1 — GH characteristic structure (cluster: eq:3+1_metric, FOSH, eq:def_s, eq:auto-1..3, Bjorhus_bc, eq:auto-31..35)
- sympy: build the 3+1 metric, define Π, Φ; verify the FOSH form reproduces
  ∂_t g = ... identities; construct s^k normalization; check u^1̂− = Π − s^iΦ_i − γ₂g
  is an eigenfield of s_k A^k with speed −(α + s_jβ^j)/... in a simple gauge.
- Deliverable: symbolic check that the characteristic decomposition is consistent
  for a generic diagonal test metric.

## T2 — Projectors and w− (cluster: eq:projection_operator, eq:auto-4, eq:auto-5, eq:wab_projection, eq:wab, eq:gh_tetrad_l, eq:GH_m, eq:GH_psi0_def, eq:bc_bjorhus)
- sympy: verify P² = P, trace properties of P^P; prove eq:auto-5 from eq:GH_m;
  derive eq:wab from eq:wab_projection + eq:GH_psi0_def symbolically (key check:
  the factor 2 and conjugate structure); verify w− is Θ-phase invariant.
- Deliverable: symbolic derivation closing the Sec. II chain.

## T3 — Bondi metrics and spin-weighted scalars (cluster: eq:BS_PFB, eq:BS_angular_gauge_condition, eq:auto-6/7, eq:true_BS, eq:BS_bondi_like, eq:UA_asymptotic, q_superA_expression, eq:dyad_property, eq:CCE_all_variables:K, eq:CCE_tetrad_m, eq:auto-36, eq:CCE_all_variables:H, eq:cce_com:H)
- sympy: instantiate the dyad on (θ, φ); verify eq:dyad_property; verify
  K = √(1+JJ̄) follows from det(h)=det(q); verify tetrad normalizations
  (l·k = −1, m·m̄ = 1, l·m = k·m = 0) against the Bondi-like metric.
- Deliverable: tetrad-consistency verifier for eq:CCE_tetrad.

## T4 — ψ₀ expression (cluster: eq:psi0_CCE)
- sympy: compute C_μνρτ l^μ m^ν l^ρ m^τ for the Bondi-like metric with generic
  J(r), β(r), K = √(1+JJ̄) and confirm eq:psi0_CCE term by term. This is the
  highest-value symbolic check (long expression; typo-prone). [OPEN] heavy CAS work.
- Numerical fallback: evaluate both sides on randomized smooth J(r), β(r) profiles.

## T5 — Coordinate maps and Jacobians (cluster: eq:cauchy_null_radius_dependence ... eq:jacobian_bondi_inertial, eq:auto-8/9/10, eq:ahat, eq:omegahat, eq:oemga_b_a, eq:expand_angular_jacobian(_nohat), eq:a_ahat_b_bhat)
- sympy: verify each Jacobian is the matrix inverse of its partner
  (eq:Jacobian_null_radius_bondi_like x eq:Jacobian_bondi_like_null_radius = 1;
  eq:Jacobian_bondi_like_inertial x eq:jacobian_bondi_inertial = 1 using
  eq:du_xhat and ωω̂ = 1); derive eq:a_ahat_b_bhat from the dyad expansions;
  verify det(middle matrix)/4 = −ω̂².
- Deliverable: complete symbolic closure of Sec. IV.A.

## T6 — Tetrad transformations, Choice 1 and Choice 2 (cluster: eq:null_vector_CCE_rhat, eq:auto-11, eq:l_GH_and_l_CCE_transform, eq:lorentz_transformation_I/II, eq:auto-12, eq:lorentz_psi0_ii, typeI_Ahat_to_Ap, eq:q_transformation_type_i, eq:auto-13..21, eq:new_m_cce_to_GH(_abstract/_B/_B_abstract), eq:psi0_prime_psi0_hat, eq:J_BS_like)
- sympy: verify Type I leaves ψ₀ invariant and Type II rescales by A²e^{2iΘ}
  (direct Weyl-scalar transformation); verify the boost factor in
  eq:l_GH_and_l_CCE_transform from eq:null_wt + eq:auto-11 + Jacobians; expand
  eq:q_transformation_type_i and the M̂/M coefficient equations from the dyad
  algebra; check m^a' from eq:new_m_cce_to_GH(_B) satisfies eq:GH_m.
- Deliverable: symbolic derivation of ψ₀' = Â²ψ̂₀ and ψ₀' = A²ψ₀.

## T7 — Angular interpolation maps (cluster: eq:duhat_x, eq:auto-22..25, eq:duhat_x_cartesian, eq:du_xhat_cartesian, eq:U_mathcal_U, eq:constraints_ca_cb)
- sympy + numpy: derive eq:U_mathcal_U and its inverse from the dyad expansions;
  numerically co-evolve forward and inverse angular maps for a toy U^(0)(u, x^A)
  and verify composition stays near identity (monitor eq:constraints_ca_cb).
- Deliverable: numerical demonstration of the cheap-inverse-map strategy.

## T8 — Teukolsky-wave analytics (cluster: eq:Teukolsky_metric, eq:auto-26/27/38, eq:Teukolsky_angular_basis, eq:Gaussian_pulse_F, eq:Teukolsky_wave_outgoing, eq:Teukolsky_h_20, eq:auto-37, eq:Teukolsky_bulk_psi0, eq:auto-28)
- sympy: linearized Einstein check R_μν = O(amplitude²) for eq:Teukolsky_metric;
  reproduce eq:auto-37 / eq:Teukolsky_bulk_psi0 from the metric with an NP tetrad;
  confirm the r^-5 law of ψ₀ and eq:Teukolsky_h_20 prefactor √(6π/5).
- Deliverable: independent analytic reference implementation for test comparisons.

## T9 — Diagnostics (cluster: eq:gauge_constraint, eq:three_constraint, eq:bondi_violation_psi3, eq:im_psi2, eq:auto-29/30)
- numpy/scri-style: implement the Bondi-constraint norms on mock waveform data;
  verify they vanish on exact perturbative Teukolsky data from T8 to truncation
  error.
- Deliverable: validation harness for stage-3.

Ordering: T1–T3 independent; T4 after T3; T5 independent; T6 after T2, T3, T5;
T7 after T5; T8 independent; T9 after T8.
