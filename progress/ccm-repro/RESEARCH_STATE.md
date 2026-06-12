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

## Queue
- iter 3+: P2a X=1e-5 overlays (our mode3-5 batteries vs their ccm/cce
  Psi4 data; matched cadence; deltas + formulation differences).
- then: P2c pulse overlays (mode6 + mirror vs their ccm/cce); P2b X=2
  config attempt; P2d BH feasible slice + resource-limit documentation;
  P2e Bondi-violation norms (scri or ported norms).
