# Iteration 13 — N6 admitted: the composite Z4c-CCM boundary-condition set (keystone)

1. **Paper anchor** — DAG node N6; composite of N3+N4+N5 with P3 D5 cadence; table grounded in P1 eq:general_CPBCs / eq:BCs-alpha / eq:BCs_lastII + P3 eq:bc_bjorhus pattern.
2. **What shipped** — `reformulate/z4c-CMM/synthesis_z4c-CCM/z4c_ccm_boundary_conditions.md` (the explicit 10-mode BC table + Bjørhus form + cadence + obligations O-N6-1..4); `packages/zccm/zccm/boundary.py` (BCTargets, physical_target, bc_targets, physical_residual — chains N1+N2+N3 on GPU); `scripts/verify_n6_composite.py` + certified output; knowledge row N6 → solid; error rows: v1 threshold slip + v2 pass.
3. **Next-3 roadmap** — N7 (well-posedness sketch document: frozen-coefficient LF program, assumptions, what is/isn't proven), N8 (GPU test suite: 1D model evolution with CCM injection, R-coefficient measurement — the final-product gate), completion audit. No same-mode loop.
4. **Simplification flag** — boundary module reuses z4c_vars constructors (no duplication) ✓.
5. **Verifier output** (verbatim, certified at results/numerical/n6_composite_check.txt):
   [PASS] mode count: BC table supplies exactly 10 == 10 incoming-mode targets (2 phys + 4 constraint + 4 gauge)
   [PASS] physical target is symmetric trace-free (2 dof exactly)
   [PASS] P1 reduction limit: psi0_CCE = 0 gives all-homogeneous targets (max |target| = 0.0e+00)
   [PASS] wire-through N1+N2+N3 chain: rel residual 2.92e-16 < 1e-12 over 1,048,576 random states
   [PASS] transparency: outgoing zero residual; incoming residual equals analytic (2k)^2 a exactly
   OVERALL: PASS (22.2s, 4x A100)
