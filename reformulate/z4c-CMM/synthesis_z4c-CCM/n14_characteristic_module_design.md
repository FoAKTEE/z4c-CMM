# N14: in-process characteristic module for AthenaK (design, iter 21)

Goal (mission 3, user directive): Cauchy and characteristic evolutions run
SIMULTANEOUSLY inside AthenaK — in-process lockstep, no postprocessing.

## Equation source (all labels resolve in knowledge-database/paper_arxiv-2007.01339)

Hypersurface hierarchy on each u = const slice, given J (Bondi-like coords;
compactified numeric forms in the same source):

| step | equation | label | numeric form |
|---|---|---|---|
| 1 | dr beta = -(r/8)[drJ drJbar - (dr(J Jbar))^2/(4K^2)] | eq:HypersurfaceBeta | eq:Betanumeric |
| 2 | dr(Q r^2) = -r^2(Lam_Q + Lam_Qbar J/K + dr ethbar J/K) + 2 r^4 dr(eth beta/r^2) | eq:HypersurfaceQ | eq:Qnumeric |
| 3 | dr U = e^{2beta}(K Q - J Qbar)/r^2 | eq:HypersurfaceU | eq:Unumeric |
| 4 | W mass-aspect equation | eq:HypersurfaceW | eq:Wnumeric |
| 5 | evolution H = du J via the H hypersurface equation | eq:HypersurfaceH | eq:Hnumeric |
| out | psi0, psi1 closed forms | eq:auto-48 (eq:FullExpandPsi01) | — |

K = sqrt(1 + J Jbar). Compactified radius y = 1 - 2 r_wt/r maps
[r_wt, inf) -> [-1, 1]. Spin weights: J:2, U:1, Q:1, beta:0, W:0.

## Module shape (athenak/src/z4c/ccm/bondi.hpp/.cpp — host-side)

- Axisymmetric m = 0 first (covers every flat-background test of
  arXiv:2308.10361); fields are spin-s functions of (y, theta):
  Gauss-Legendre theta collocation (n_th ~ 24), spin-weighted associated
  Legendre basis for eth/ethbar (m = 0 ladder operators); radial: Chebyshev
  or 6th-order FD on y in [-1, 1], n_y ~ 32-48. Cost per step trivially
  host-sized (~1e4 points): no GPU needed; runs on rank 0 and broadcasts
  psi0(theta) (MPI).
- Stage A (iter 22): STANDALONE LINEAR solver — linearize the hierarchy
  symbolically (sympy transcription from the ledger labels above; verifier
  gates the linearization against the full forms at sample points), evolve
  Teukolsky worldtube data (analytic, from the N12-corrected solution), and
  gate: strain/psi0 at the worldtube and scri vs scripts/teuk_exact_waveform.py
  and the eq:teuk-psi0 datum chain. ALSO cross-oracle: the same worldtube
  data through SpECTRE CharacteristicExtract (ref-code/spectre-cce).
- Stage B (iter 23): worldtube Bondi data from the LIVE Cauchy state —
  reuse the cce/ sphere interpolation in-memory (no files): 4-metric +
  derivatives on the worldtube sphere -> Bondi-like scalars (worldtube
  transformation chain of 2007.01339 Sec. III, ledger labels eq:ifc*).
- Stage C (iter 24): lockstep coupling — characteristic du = Cauchy dt
  (subcycle if CFL differs); each cycle: Cauchy Task_End -> worldtube data
  -> hierarchy sweep -> H step -> psi0 at the worldtube -> Type-II boost
  (N2) -> Z4cCCMDatumPsi0 mode 4 (live) consumed by the next cycle's
  boundary task (one-step lag = P3's initialization ordering, O-N6-3).
- Nonlinear terms: added incrementally after the linear gates pass; the
  full forms are transcribed in the ledger (no new derivation needed).

## Verification ladder (each its own loop iteration)

1. linear hierarchy vs analytic Teukolsky (worldtube + scri), spectral/FD
   convergence measured;
2. SpECTRE CharacteristicExtract cross-oracle on identical worldtube data;
3. live-CCM Teukolsky battery vs the N12 analytic truth (the datum now
   carries the characteristic solution instead of the analytic formula —
   they must agree at linear order: a closed consistency loop);
4. injection test (paper Sec. V.C) driven from characteristic initial data;
5. X = 2 with the dissipation/ID config note (error-db 2026-06-12T17:29).

## Status

Design committed iter 21. Implementation begins iter 22 (stage A).
