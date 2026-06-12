# Iteration 11 — N4 admitted: CPBCs retained verbatim; channel taxonomy refined

1. **Paper anchor** — DAG node N4; P1 eq:general_CPBCs / eq:theta _modes / eq:sys-theta-Z; P3 constraint-matching design argument.
2. **What shipped** — `scripts/verify_n4_cpbc_compat.py` + certified output `results/numerical/n4_cpbc_compat_check.txt`; knowledge row N4 → solid (predecessor N3); error rows: v1 negative-control taxonomy slip (h_xx is gauge, not constraint-sourcing) + v2 pass; DAG doc N4 [SOLID].
3. **Next-3 roadmap** — N5 (gauge sector: paper-1 gauge BCs replace P3 Sommerfeld; verifiable core = gauge-channel decoupling + the h_xx/lapse/shift channel bookkeeping from this iteration's taxonomy), N6 (composite 10-mode BC table + zccm boundary module), N7 (well-posedness sketch). No same-mode loop.
4. **Simplification flag** — n/a.
5. **Verifier output** (verbatim, certified):
   [PASS] TT channel: linearized H == 0 identically (off-shell)
   [PASS] TT channel: linearized M_i == 0 identically (off-shell)
   [PASS] => dot Theta, dot Z_i sources vanish; U-_Theta, U-_Z stay 0 under ANY CCM physical datum
   [PASS] NEGATIVE CONTROL (vector channel h_xy): M_y = -Derivative(G(t,x),t,x)/2 != 0
   [PASS] taxonomy: longitudinal h_xx(t,x) is ALSO constraint-free (gauge channel)
   [PASS] NEGATIVE CONTROL (tangential trace h_yy = h_zz = G): H = -2*Derivative(G,(x,2)) != 0
   [PASS] GPU sweep: TT residual 0.0 vs control amplitude 1.23e+01, 4,194,304 modes (4x A100)
   OVERALL: PASS (12.4s)

Physics note (feeds N5/N6): boundary-channel taxonomy — TT = physical (CCM-driven),
h_xx = gauge, h_xA = constraint-vector (sources M_A), tangential trace =
constraint-scalar (sources H). This is the linearized shadow of P1's
4 constraint + 4 gauge + 2 physical incoming-mode split.
