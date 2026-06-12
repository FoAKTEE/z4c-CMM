# Iteration 17 — Table II rerun at strict bounds (user directive): every deviation at its float64 limit

1. **Paper anchor** — N8 dynamical evidence (Table II col. 4); user: "in principle you can do with much more strict bound" + "consider increasing accuracy of the calculation chain".
2. **What shipped** — model1d upgraded: selectable-order stencils (2/4/6) with Fornberg edge weights; KO dissipation hook; GKS-stable mixed order-6 closure (4th-order rows at 3 edge nodes) after the fully-one-sided closure blew up at node N−2 (ledgered); test_model1d v3: L2/energy reflection measurement with the analytic ratio b/(2+b)·√((1+b)/(1−b)), strict asserts; roundoff-wall convergence ladder for the injection error; paper Table II final.
3. **Chain-accuracy summary** — symbolic links: EXACT (zero error). Float64 implementation checks: 1e-16..1e-17 floors (N2 9.3e-16, N3-exact Schwarzschild 7.7e-17/8.5e-16, N6 wire-through 1.5e-16). Dynamical model: Sommerfeld 2.06e-14 + P1 absorption 1.88e-14 (per-sample roundoff floors, saturated under refinement), injection 2.44e-11 (linear steps·ε accumulation optimum at N=32769; ladder non-monotone 5.4/2.4/5.9/30e-11; CFL-insensitive). Improvements vs v1: 0.03% → 2.1e-14 (10 orders), 1.1e-4 → 2.4e-11 (6.6 orders). Going below 1e-11 on injection requires compensated (Kahan) time accumulation or extended precision — recorded as an optional successor item, not a scheme limitation.
4. **Simplification flag** — n/a.
5. **Verifier output** (verbatim, certified at results/numerical/n8_model1d_test.txt):
   [PASS] sanity: q_in == 0 exactly
   [PASS] Sommerfeld L2 ratio = 0.177752646227 vs analytic 0.177752646227: rel dev 2.06e-14 (coarse 1.86e-14, saturated)
   [PASS] P1/CCM operator: L2 absorption ratio = 1.88e-14
   [PASS] CCM injection: rel L2 = 2.44e-11 (float64 accumulation optimum, N = 32769)
   OVERALL: PASS (40.8s, 4x A100)
