# Iteration 15 — N8 admitted: the GPU final product (packages/zccm) with dynamical verification

1. **Paper anchor** — DAG node N8; P1 eq:eomgammasA model; P3 test-suite pattern; B7 reflection diagnostic realized numerically.
2. **What shipped** — `packages/zccm/zccm/model1d.py` (characteristic/Bjørhus BC enforcement, jit/scan FD2+RK4), `packages/zccm/tests/test_model1d.py` (4 configs fanned over 4 GPUs), `packages/zccm/README.md`; knowledge row N8 → solid; error rows: 2 structural fails (char labels swapped; orientation mirror — root causes ledgered) + pass.
3. **Next-3 roadmap** — completion audit: re-run all verifiers, structure check, RESEARCH_STATE marker scan, archive loop state, sweep commit, disarm, promise. N9 remains [FUTURE] documented obligation (not in exit conditions).
4. **Simplification flag** — model1d self-contained by design (model problem ≠ production boundary module); no duplication of z4c_vars.
5. **Verifier output** (verbatim, certified at results/numerical/n8_model1d_test.txt):
   [PASS] sanity: pure outgoing initial data has q_in == 0 exactly (max 0.0e+00)
   [PASS] Sommerfeld: measured q-ratio = 0.1305 vs exact b/(2+b) = 0.1304 (dev 0.03% < 3%); u-amplitude R = 0.2423 vs mirror formula 0.2422
   [PASS] P1/CCM operator: R = 3.24e-15 (fine) < q-ratio_Sommerfeld/100 = 1.30e-03 (absorbing at machine precision)
   [PASS] CCM injection: interior reproduces the exact incoming pulse, rel L2 = 1.09e-04 < 1% (two-way transparency)
   OVERALL: PASS (18.8s, 4x A100)

Physics note: the dynamical run REPRODUCES paper-3's bottleneck quantitatively
(naive Sommerfeld reflects at exactly the analytic coefficient) and the new
scheme's operator removes it to machine precision while letting CCE data
through — the formulation's three claims demonstrated in evolution, not just
algebra.
