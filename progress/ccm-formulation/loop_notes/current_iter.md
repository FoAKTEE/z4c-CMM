# Iteration 10 — N1 admitted: Z4c worldtube variable maps + first zccm GPU module

1. **Paper anchor** — DAG node N1; P1 eq:Z4_decomp_first group → P2 Sec II.B worldtube input.
2. **What shipped** — `packages/zccm/` package seed (`zccm/z4c_vars.py`: Z4cState, to_adm/from_adm, four_metric(+inverse), boost_factor [N2], psi0_to_physical_datum [N3]); verifiers `scripts/verify_n1_varmap.py` (exact full-symbolic, 3.9 s) + `packages/zccm/tests/test_z4c_vars.py` (4× A100, 1,048,576 states, 22.3 s); knowledge row N1 → solid; error rows: verifier-side assertion-sign + jit-concretization fail, then pass; DAG doc N1 [SOLID]; figure refreshed.
3. **Next-3 roadmap** — N4 (CPBC retention statement + corner-compatibility obligation), N5 (gauge-sector replacement statement), N6 (composite scheme assembly: explicit 10-mode BC table + zccm boundary module). Distinct targets — no same-mode loop.
4. **Simplification flag** — n/a (modules small and single-purpose so far).
5. **Verifier output** (verbatim; certified at results/numerical/n1_varmap_check.txt and n1_varmap_gpu_test.txt):
   [PASS] chi recovery … residual det(gamma)*chi^3 - det(gt) == 0  (t=0.3s)
   [PASS] K recovery … vanishes given trace-free At  (t=1.7s)
   [PASS] At recovery … == -gt_ij/3 * gt^{kl}At_kl — vanishes given trace-free At  (t=2.4s)
   [PASS] 4-metric inverse: g4 * g4inv == identity (full symbols)  (t=2.6s)
   [PASS] KINEMATIC CLOSURE with D_i beta^i = beta^k_,k - (3/2) beta^k (ln chi)_,k … EXACT  (t=2.7s)
   [PASS] NAIVE reading D = beta^k_,k leaves residual = c * gamma_ij beta^k (ln chi)_,k with c = {1}  (t=3.8s)
   OVERALL: PASS (3.9s)
   JAX devices used (4): gpu:0..3
   [PASS] round trip … max rel residual 4.10e-16 < 1e-12 over 1,048,576 states  (t=22.3s)
   [PASS] g4 * g4inv == identity: 8.88e-16 < 1e-11
   [PASS] boost factor matches closed form (0.00e+00), positive for subluminal shift
   [PASS] physical datum TT with K = 4 (all residuals 0.0)
   OVERALL: PASS (22.3s)

Convention note (new, feeds the formulation document): P1's ∂_tχ equation
divergence is the FULL covariant D_iβ^i — proven by kinematic closure, with the
naive ∂_kβ^k reading leaving residual γ_ij β^k ∂_k ln χ (c = 1).
