# CCM-paper reproduction (mission 4) — Verified Research State

User directive (2026-06-13, verbatim): "proceed to apple to apple
numerical test reproduction and figure exact reproduction (within
resource limitation) of CCM paper" — arXiv:2308.10361.

Anchor: the authors' published simulation data + figure scripts,
ref-code/ccm-figures (github.com/Sizheng-Ma/CCM-Figures, cloned depth-1,
1.6 GB): X=1e-5/{277ccm,288ccm,299ccm,288cce,299cce}, X=2/{355ccm,366ccm,
377ccm,reference,...}, BH/{355ccm,366ccm,377ccm,366cce,377cce,reference},
pulse_on_characteristic_grid/{355ccm,366ccm,377ccm,cce} — Norms.dat /
Psi4.dat etc. + plot.ipynb per test.

## Iteration 1 (P1 figure reproduction): DONE
All four plot.ipynb generators ported headlessly (code cells extracted,
magics stripped, Agg backend, save-all-figures epilogue; ports live next
to the notebooks as plot_ported.py) and executed from their data dirs:
- X=1e-5: 2 figures; X=2: 8; BH: 2; pulse: 2 — 14 PDFs under
  results/figures/ccm_paper_repro/<test>/.
Deviations: serif font 'Georgia' unavailable (matplotlib fallback);
constrained_layout vs subplots_adjust warnings (their mpl version older);
content unaffected (data-driven curves).

## Iteration 2 (P1 audit + completion): DONE — M4-P1 [SOLID]
Port fix (error-db): jupyter per-cell figure semantics restored (the six
X=2 Bondi-violation cells had overdrawn the constraints figure). 20
figures now: X=1e-5: 2, X=2: 14 (7 quantities + con + ALL SIX
BondiViolationMin panels from the published violation.dat — no scri
needed), BH: 2, pulse: 2. Full mapping + coverage audit in
results/numerical/ccm_repro_p1_check.txt; spot-checks match the paper.
Not reproducible from the repo (data absent, documented): 12 appendix
panels (low_amp_{Psi0-4,News}, kerr_{Psi0-3,News,Strain}); 2 setup
illustrations out of scope.

## Iteration 3 (P2a X=1e-5 overlays): DONE — M4-P2a [SOLID]
Direct scri-to-scri overlay: ours (authors' CCE binary on exact analytic
worldtube data, X/radius matched) vs their published 299ccm/299cce —
1.92e-3 of peak, attributed by measurement to THEIR runs' default-
tolerance floor at X=1e-5 (tightening OUR tolerances: 25x improvement to
1.16e-4 vs the closed form). Time shift +1.000 = pure convention
(measured). Their CCM-vs-CCE = 3.6e-6 of peak (the matching-term effect,
quantified from their own data). Transcript ccm_repro_p2a_check.txt.

## Iteration 4 (P2c pulse overlays): DONE — M4-P2c [SOLID]
Their CCM InnerPsi0 trough vs our mirror x c2 (derived): amplitude ratio
0.9862, offset +2.33 (geometric expectation +1.3..+3.0). psi4 features
map 1:1 (+3..+4 offsets); fitted scale 2.73 = c2 within 5.4% (finite-r
residual documented). KEY: their live data matches our MIRROR, not our
live - localizes our O(dt)/geometry coupling artifacts (iter-45 reading
amended); actionable for the linear-in-time BC refinement. Transcript
ccm_repro_p2c_check.txt.

## Iteration 5 (P2b/P2d/P2e + exit): DONE — M4-P2bde [SOLID]
X=2: stability ACHIEVED (two t=80 batteries; the mission-3 NaN note
resolved); quantitative overlay blocked by the linear-ID constraint
violation (delay +6.6/+17.6 scaling with dissipation vs their +0.7;
XCTS out of scope) — limit documented, feasible slice delivered. BH:
out-of-resource documented (1000M excision = O(days); characteristic
side already validated on their binary). P2e: C_psi4 Bondi relation
verified on our oracle at 5.3e-4 of peak; C_psi3/psi2 eth conventions
queued. Transcript ccm_repro_p2bde_check.txt.

## MISSION 4 EXIT STATUS (iter 5): ALL CONDITIONS HOLD
1. P1 complete [SOLID iter 2]: 20/20 reproducible figures from their
   data via their scripts; absent-data panels + illustrations documented.
2. P2 overlays admitted: X=1e-5 [SOLID iter 3, 1.92e-3 of peak,
   attributed]; pulse [SOLID iter 4, 1.4% t1 amplitude]; X=2 limit
   documented with stable feasible slice [iter 5]; BH out-of-resource
   documented with the validated characteristic slice [iter 5].
3. Per-figure provenance + delta tables: ccm_repro_p1_check.txt,
   ccm_repro_p2a_check.txt, ccm_repro_p2c_check.txt,
   ccm_repro_p2bde_check.txt.
4. Ledger (M4-P1, M4-P2a, M4-P2c, M4-P2bde) + RESEARCH_STATE consistent;
   all commits pushed.

## Post-exit quality items (documented, not exit conditions)
- O(dt) live-coupling refinement with the P2c-localized O(30%) target.
- C_psi3/C_psi2/C_Im psi2 eth-convention bookkeeping for our outputs.
- psi4 finite-r-vs-scri 5.4% residual factor derivation.
- XCTS-class constraint-solved ID for quantitative X=2.

## Queue
- iter 3+: P2a X=1e-5 overlays (our mode3-5 batteries vs their ccm/cce
  Psi4 data; matched cadence; deltas + formulation differences).
- then: P2c pulse overlays (mode6 + mirror vs their ccm/cce); P2b X=2
  config attempt; P2d BH feasible slice + resource-limit documentation;
  P2e Bondi-violation norms (scri or ported norms).

# Mission 5 (production paper, USER DIRECTIVE 2026-06-13)

Full production-level formulation paper following the CCM paper's
organization with the complete Z4c-CCM equation chain + derivation DAG.
Loop state .claude/ralph-loop.local.md; plan: iter 1 skeleton + Sec II
(DONE: full Z4c system, 10-mode table, N-node content redistributed; 8
pages, 0 unresolved refs); iter 2 Sec III characteristic; iter 3 Sec IV
matching + DAG figure; iter 4 Sec V consolidation + I/VI + appendices;
iter 5 consistency pass + exit.
m5-iter 3 DONE: Sec IV matching complete (contract, six-row map,
retarded geometry + O-M5-1, cadence) + the TikZ derivation DAG fig:dag.
12 pages clean.

## X=2 FULL CAMPAIGN (user-authorized 8 GPUs / 2 days, 2026-06-13)
Launched detached (scripts/x2_full_campaign.sh): causally-disconnected
reference 768^3 box +-192 (boundary influence at r=36 only from t=156)
+ boundary-41 Sommerfeld/CCM(mode-3) pairs at h = 0.5/0.25/0.125; all at
the att2-stable X=2 config (diss=0.1, kappa1=0.1, nghost=4). Goal: the
paper-protocol X=2 error hierarchy (E[run] vs reference << E_h) with
shared linear ID (internally consistent; the XCTS difference documented).
STATUS: results/numerical/ccm_repro_x2_full/STATUS.log; hourly checks.
m5-iter 4 DONE: Sec V repro subsection (sec:repro + 3 overlay figs),
intro final, conclusions, appendices A/B/C; 15 pages, 0 unresolved.
X=2 campaign healthy (5/7 done; preliminary boundary error 1.1e-2 =
Somm=CCM, truncation 0.6; h=0.125 pair running).
