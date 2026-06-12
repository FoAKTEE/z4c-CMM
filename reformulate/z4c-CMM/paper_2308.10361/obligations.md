# obligations.md — arXiv:2308.10361

Open obligations with expected evidence type. Markers per markers.md.

1. [OPEN] Symbolic verification of eq:psi0_CCE (Weyl scalar from Bondi-like metric).
   Closing evidence: CAS computation of C l m l m for generic J(r), β(r) matching
   the published expression (implementation_plan_python.md T4). Owner: stage 2.
2. [OPEN] Symbolic verification of the Jacobian-inverse pairs
   (eq:Jacobian_null_radius_bondi_like vs eq:Jacobian_bondi_like_null_radius;
   eq:Jacobian_bondi_like_inertial vs eq:jacobian_bondi_inertial) and of
   eq:a_ahat_b_bhat, ωω̂ = 1. Closing evidence: sympy matrix products = identity (T5).
3. [OPEN] Symbolic verification of the Type I/II transformation laws of ψ₀
   (eq:auto-12, eq:lorentz_psi0_ii) directly from the Weyl-scalar definition.
   Closing evidence: sympy NP-tetrad computation (T6).
4. [OPEN] Symbolic check that m^a' of eq:new_m_cce_to_GH and eq:new_m_cce_to_GH_B
   satisfies eq:GH_m exactly (orthogonality + normalization), including the
   K = √(1+JJ̄) typo-correction claim C13. Closing evidence: sympy (T3, T6).
5. [OPEN] Verification of eq:Teukolsky_metric as a linearized vacuum solution and
   of the analytic waveform set eq:auto-37 / eq:Teukolsky_bulk_psi0 (prefactors,
   r^-5 law). Closing evidence: sympy linearized Ricci + NP scalars (T8).
6. [OPEN] Reproduction of at least one quantitative numerical result (e.g. the
   X = 2 deviation reduction, Fig. 6) from the public CCMData repository.
   Closing evidence: independent analysis of published data (stage 3).
7. [OPEN] Well-posedness of the matched CCM problem. The characteristic system is
   only weakly hyperbolic (Giannakopoulos et al.); no proof exists. Closing
   evidence: proof or counterexample — beyond this paper's scope. [FUTURE] for
   this pipeline; [BLOCKING] for any claim of continuum convergence (claim C5).
8. [FUTURE] Gauge-sector matching: extend CCM to the four incoming gauge modes
   (currently Sommerfeld). Stated as future work by the authors; required to
   eliminate the tertiary reflection (claim C10).
9. [FUTURE] Robust stability tests with high-frequency/non-smooth data and BBH
   mergers. Required before relying on CCM in production waveforms.
10. [OPEN] Regularity/invertibility of the angular maps (assumption 15): monitor
    eq:constraints_ca_cb (or equivalent) in any reimplementation. Closing
    evidence: numerical monitoring in T7.
11. [OPEN] Unit/frame audit: confirm spin-weight bookkeeping of â, b̂, a, b
    (sw 2/0) and of 𝒰^(0) against the dyad conventions, since the tex contains a
    sign-convention-sensitive dyad choice (q^A ∂_A = −∂_θ − i/sinθ ∂_φ). Closing
    evidence: sympy dyad algebra (T5/T7).
12. [OPEN] The commented-out tex equations (eq:constraints_ca_cb, eq:cce_com,
    eq:auto-31..36, eq:CCE_all_variables:H) are registered but not part of the
    published text; decide in stage 2 whether they are load-bearing for the
    implementation (the cce_com hierarchy IS load-bearing for any CCE
    reimplementation; its source terms live in Bishop:1997ik / Handmer:2014qha).
