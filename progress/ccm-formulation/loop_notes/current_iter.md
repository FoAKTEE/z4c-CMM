# Iteration 14 — N7 sketch admitted: LF/Kreiss program documented, model-level core verified

1. **Paper anchor** — DAG node N7; P1 eq:lapl-four / eq:bc_general / sol:lf / sol:lf2 machinery; P3 weak-hyperbolicity caveat.
2. **What shipped** — `reformulate/z4c-CMM/synthesis_z4c-CCM/n7_wellposedness_sketch.md` (assumptions 1–6 explicit; open items 1–4 as [FUTURE] obligations with human sign-off requested); `scripts/verify_n7_lf_sketch.py` + certified output; knowledge row N7 → solid (sketch scope); error row: pass (one cosmetic root-pairing fix, noted in trial row).
3. **Next-3 roadmap** — N8 (final product gate: 1D model evolution test with CCM injection on 4 GPUs measuring reflection coefficients + package README; then completion audit), N9 stays [FUTURE] (documented obligation, not in exit conditions). No same-mode loop.
4. **Simplification flag** — n/a.
5. **Verifier output** (verbatim, certified at results/numerical/n7_lf_sketch_check.txt):
   [PASS] LF interior roots are exactly {s/(1+b), -s/(1-b)} (set match)
   [PASS] admissibility: Re lam+ > 0 iff Re s > 0
   [PASS] L=0..3: D_L(s) = (2s)^(L+1) exactly; uniform ratio 2^(L+1); datum-independent (True)
   [PASS] Sommerfeld determinant D = s(2+b)/(1+b) != 0 (well-posed; defect is reflection per N5)
   [PASS] GPU sweep: |D_L(s)|/|s|^(L+1) == 2^(L+1) to 3.2e-14 over 4,194,304 samples, L = 0..3
   OVERALL: PASS (13.5s, 4x A100)
