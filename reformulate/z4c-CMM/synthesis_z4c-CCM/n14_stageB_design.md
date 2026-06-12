# N14 stage B: live worldtube map + embedding decision (iter 38)

## The worldtube-anchored gauge makes the BC map LOCAL

Anchor the linear Bondi-like gauge at the worldtube (zeta(r_wt) = 0 for
all u) instead of scri. Consequences (l=2 m=0, linear):
- J_wt  = (h_thth^ - h_phph^)/2 |_wt           (LOCAL Cauchy data)
- U_wt  = -h_{u theta}/r^2 |_wt                 (LOCAL; h_t-theta sector)
- beta_wt from h_rr|_wt and the LOCAL ODE RHS values of the gauge fields
- Q_wt, W_wt from first radial/time derivatives of h at r_wt (LOCAL)
- H_wt  = du J_wt                               (LOCAL time derivative)
- psi0 (the matching datum) is GAUGE-INVARIANT at linear order
  (background Weyl = 0); VERIFIED EXACTLY in the series algebra: psi0
  series is term-by-term invariant under k=0 (constant-in-r) gauge shifts
  (dr annihilates them) — transcript in the iter-38 commit.
- Cost: scri quantities (news/strain) carry the inertial-frame correction
  (standard CCE machinery) — NOT needed for CCM, where only psi0 at the
  worldtube feeds back. Ledgered as the scri-extraction note.

## Data contract (what AthenaK supplies per worldtube exchange, l=2 m=0)

On the sphere r = r_wt, projected on the l=2 m=0 harmonics (the cce/
sphere-interpolation machinery already computes metric components there):
  h_rr, h_{r theta}, h_thth^, h_phph^, h_{t theta};
  dr h_{...} (first radial derivatives), dt h_{...} (first time
  derivatives) for the J/Q/W rows. All linear-order; the nonlinear
  worldtube transformation (2007.01339 ifc chain) is the documented
  successor for X = 2 production.

## Embedding decision: Julia-in-process via libjulia

- CHOSEN: AthenaK's Z4c task list calls the ZccmJl solver in-process
  through the libjulia C API (jl_init once at startup; per-cycle calls
  pass raw double arrays; rank-0 only + MPI broadcast of psi0(theta)).
  Honors the Julia-only development constraint; genuine lockstep (no
  files); the solver is host-side small (n ~ 65 dense ops, microseconds
  in Float64).
- REJECTED: C++ port (violates the migration constraint for new
  development; revisit only if libjulia overhead proves prohibitive —
  measure first); file handshake (postprocessing — disqualified by the
  mission directive).
- Float64 suffices for production coupling (X = 2 datum is O(1)-relevant;
  the BigFloat path is the validation harness for the peeling-tiny linear
  datum — iter-37 scope note).

## SpECTRE cross-oracle (prep)

ref-code/spectre-cce holds CharacteristicExtract + PreprocessCceWorldtube
(v2026.06.09.01, runnable). Plan: write the Teukolsky worldtube modes in
the SpECTRE worldtube H5 format (PreprocessCceWorldtube ingests metric
modes at r_wt), run CharacteristicExtract, compare its psi0/news against
the ZccmJl solver on identical data. Execution: a later loop iteration
(stage C verification ladder, design-doc item 2).
