# Iteration 9 — N3 admitted: the physical-dictionary identity (linearized core)

1. **Paper anchor** — DAG node N3; P3 eq:wab / eq:bc_bjorhus pattern landed in P1 eq:BCs_lastII slot; P1 frame eq:nullvector-k.
2. **What shipped** — `scripts/verify_n3_dictionary.py` + certified output `results/numerical/n3_dictionary_check.txt`; knowledge rows: N3 → solid (predecessor N2); error rows: v1 null-leg convention fail (ledgered with diagnosis) + v2 pass; DAG doc updated (N3 [SOLID]); figure re-rendered.
3. **Next-3 roadmap** — N1 (Z4c → worldtube 4-metric round-trip; seeds `packages/zccm/`), N4 (CPBC retention + corner obligation note), N5 (gauge-sector replacement statement). Distinct targets — crash-triage clean (the N3 v1 fail was remediated by convention fix, not parameter tweak).
4. **Simplification flag** — n/a.
5. **Verifier output** (verbatim, certified at results/numerical/n3_dictionary_check.txt):
   [PASS] tetrad sanity (flat): l.l = 0, k.k = 0, l.k = -1, m.m = 0, m.mbar = 1  (t=0.6s)
   [PASS] psi0 = c * (d_t + d_s)^2 (F1 + i F2) OFF-SHELL, pure number c = 1/4  (t=5.2s)
   [PASS] w-_AB = K (psi0 mbar mbar + conj), component-exact for arbitrary F1, F2, pure number K = 4  (t=5.2s)
   JAX devices used (4): ['gpu:0', 'gpu:1', 'gpu:2', 'gpu:3']
   [PASS] GPU sweep: max relative residual 0.00e+00 < 1e-9 over 4,194,304 OFF-shell Fourier modes (independent omega, k)  (t=20.1s)
   OVERALL: PASS  (wall clock 20.1s)

Physics note: the incoming-mode extractor on the boundary is the OUTGOING null
direction operator l·∂ (it annihilates outgoing waves f(t−s)) — consistent with
paper-1's (r² l̊·∂)^{L+1} BC operator; v1's (∂_t − ∂_s) guess was the flipped leg.
