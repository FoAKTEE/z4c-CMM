# Iteration 12 — N5 admitted: gauge sector — P1 BCs replace Sommerfeld, bottleneck removed at model level

1. **Paper anchor** — DAG node N5; P1 eq:BCs-alpha / eq:general_BCs_gauge_first / eq:eomgammasA; P3 Sommerfeld-gauge future-work bottleneck.
2. **What shipped** — `scripts/verify_n5_gauge.py` + certified output `results/numerical/n5_gauge_check.txt` (first-try PASS, 12.2 s); knowledge row N5 → solid (predecessor N4); DAG doc N5 [SOLID]; RESEARCH_STATE rows N4+N5; nodal note 10-iter window refreshed.
3. **Next-3 roadmap** — N6 (composite scheme: explicit 10-mode BC table document + zccm boundary module on GPU), N7 (well-posedness sketch document), N8 (GPU test suite: Teukolsky/pulse + R-coefficient). No same-mode loop.
4. **Simplification flag** — watch: N6 boundary module must reuse z4c_vars constructors.
5. **Verifier output** (verbatim, certified):
   [PASS] TT channel: gamma^{ij} d_t gamma_ij == 0 (1+log lapse source K == 0)
   [PASS] TT channel: Gammat^i = -d_j ht^{ij} == 0 (gamma-driver source)
   [PASS] => linearized 1+log + gamma-driver receive ZERO source from the CCM channel
   [PASS] dispersion: outgoing lam = w/c+ and incoming lam = -w/c- both solve P1 eq:eomgammasA
   [PASS] naive Sommerfeld reflection R = -b(1-b)/((1+b)(2-b)) != 0 for b != 0
   [PASS] paper-1 boundary-adapted operator reflection R = 0 IDENTICALLY
   [PASS] GPU sweep: residuals 2.2e-16 / 0.0; |R_naive| >= 0.023 for b in [0.05,0.6]; 4,194,304 samples (4x A100)
   OVERALL: PASS (12.2s)
