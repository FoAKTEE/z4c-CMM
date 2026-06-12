# claims.md — arXiv:2308.10361

Schema per `_common/contracts/research_admission_contract.md`: each claim carries
working context, evidence type, dependencies, assumptions (by number from
assumptions.md), and a status marker. All claims below are transcribed from the
paper; none has been independently verified in this stage, so the strongest
admissible status here is literature-grounded [UNCHECKED] (or [ASSUMPTION] for
imported postulates).

## C1 — CCM construction is exact (no approximations)
- Working context: vacuum GR, GH Cauchy + SpECTRE characteristic system, outer
  boundary = worldtube; physical modes only.
- Claim: the mapping ψ₀(characteristic) → w−|_BC (Cauchy) via Jacobians + Type I/II
  Lorentz transformations (Choice 1: eq:psi0_prime_psi0_hat; Choice 2: eq:auto-21)
  is exact in the continuum — free of symmetry/perturbative approximations.
- Evidence type required: symbolic derivation. Current evidence: literature
  grounding (paper derivation, transcribed in derivation.md).
- Dependencies: eq:wab, eq:l_GH_and_l_CCE_transform, eq:psi0_CCE, eq:J_BS_like,
  Jacobian cluster (logic.md Sec. IV.A).
- Assumptions: 1, 2, 6, 9, 13, 14, 15.
- Status: [UNCHECKED] — closes via implementation_plan_python.md T5+T6.

## C2 — l^a' ∝ l^â at the worldtube
- Working context: worldtube surface, null-radius construction eq:null_wt.
- Claim: the characteristic outgoing null vector is proportional to the Cauchy one,
  with boost factor (α − γ_i'j'β^i's^j') e^{−2β̂} (eq:l_GH_and_l_CCE_transform).
- Evidence type required: symbolic derivation. Dependencies:
  eq:null_vector_CCE_rhat, eq:auto-11, eq:null_wt, eq:gh_tetrad_l.
- Assumptions: 6, 9. Status: [UNCHECKED] (T6).

## C3 — w− is invariant under the residual phase freedom Θ
- Claim: ψ₀'m̄m̄ + c.c. is invariant under m → m e^{iΘ} since ψ₀ → e^{2iΘ}ψ₀.
- Evidence type required: symbolic derivation (short). Dependencies: eq:wab,
  eq:lorentz_psi0_ii. Assumptions: 14. Status: [UNCHECKED] (T2).

## C4 — Constraint-preserving BCs need no matching
- Claim: the forty constraint-related incoming fields are correctly handled by
  Lindblom Eqs. (63)-(65) without approximation; CCM for them is unnecessary and
  not well-motivated (constraints are local).
- Evidence type: citation (Lindblom:2005qh) + argument. Dependencies: FOSH,
  Bjorhus_bc. Assumptions: 2, 3. Status: [ASSUMPTION] (imported).

## C5 — After matching physical+gauge BCs, CCM converges to the infinite-domain problem in the continuum limit
- Claim as stated in Sec. II; note the gauge part is NOT yet matched, so the claim
  is conditional/prospective.
- Evidence type required: proof (well-posedness + convergence). Status:
  [CONJECTURAL]/[OPEN] — blocked by weak hyperbolicity (assumption 8) and unmatched
  gauge modes (assumption 5).

## C6 — No numerical instabilities in the tested scenarios
- Working context: smooth axisymmetric data, three resolutions, >1000 code units;
  flat Teukolsky (X = 1e-5, 2), Kerr χ = 0.5, characteristic pulse injection.
- Evidence type: numerical simulation (paper's). Dependencies: eq:bc_bjorhus
  pipeline; diagnostics eq:gauge_constraint, eq:three_constraint.
- Assumptions: 8, 10, 16, 17. Status: [UNCHECKED] empirical (external simulation
  data at CCMData repo; not reproduced here). Does NOT establish well-posedness.

## C7 — CCM leaves GH constraint violations unchanged vs CCE, and they converge with resolution
- Evidence type: numerical simulation (paper Figs. 3, 5, 8, 12). Status:
  [UNCHECKED] empirical.

## C8 — In the perturbative test (X = 1e-5), CCE/CCM/analytic differences are below numerical error
- Claim: backscatter is negligible (ψ₀ ~ r^-5, eq:Teukolsky_bulk_psi0), so the
  matching term is unresolvable; agreement with eq:Teukolsky_h_20 cross-checks the
  code. Evidence type: numerical + controlled approximation. Assumptions: 11.
  Status: [UNCHECKED] empirical.

## C9 — In the nonlinear test (X = 2), CCM reduces waveform deviation from the reference by ~1 order of magnitude
- Also: CCM reduces Bondi-gauge constraint violations (eq:bondi_violation,
  eq:im_psi2) except C_ψ3, C_ψ4 where CCE ≈ CCM.
- Evidence type: numerical simulation vs SpEC reference. Assumptions: 12.
  Status: [UNCHECKED] empirical.

## C10 — A tertiary spurious reflection persists, attributed to unmatched gauge BCs
- Claim: pulse at t = 2r_out + (r_out − r_c) ≈ 103 arises from reflection at the
  boundary; attributed (not proven) to the Sommerfeld gauge subset.
- Evidence type: numerical observation + plausibility argument. Status:
  [HYPOTHESIS] (attribution), [UNCHECKED] empirical (timing).

## C11 — Kerr test: CCE-CCM difference below numerical error; QNM ringing matches perturbation theory
- Evidence type: numerical + external package (qnm). Assumptions: 18, 19.
  Status: [UNCHECKED] empirical.

## C12 — Characteristic-pulse injection: CCM transmits the pulse into the Cauchy domain; CCE reflects it
- Timing checks: entry at t₁, exit at t₂ = t₁ + 2R. Evidence type: numerical.
  Status: [UNCHECKED] empirical.

## C13 — K = √(1+JJ̄) (typo correction to Moxon:2020gha Eq. 10e)
- Claim: eq:CCE_all_variables:K is the correct expression. Evidence type: symbolic
  derivation from det(ĥ) = det(q). Status: [UNCHECKED] (T3; easy to close).

## C14 — eq:psi0_CCE is the correct Bondi-like ψ₀
- The long closed-form expression; imported from Moxon:2020gha. Evidence type
  required: symbolic derivation (CAS). Status: [UNCHECKED] (T4) — highest-value
  independent check.
