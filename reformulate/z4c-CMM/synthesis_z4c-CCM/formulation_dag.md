# Z4c-CCM Formulation DAG — synthesis of arXiv:1010.0523v2 + 2007.01339 + 2308.10361

**Goal (mission 2):** a NEW FORMULATION — Cauchy-characteristic matching (CCM) for the
**Z4c** evolution system. Organizing principle (per the mission hint): **physical
degrees of freedom**, **characteristic modes**, and **boundary conditions**.

**Thesis.** At a timelike worldtube boundary, Z4c has exactly ten incoming
characteristic fields, which paper 1 (Ruiz–Hilditch–Bernuzzi 2010) splits into
**4 constraint + 4 gauge + 2 physical** modes and supplies with ten boundary
conditions. Paper 3 (Ma et al. 2023) demonstrates, for the GH system, that CCM
need only act on the **2 physical** incoming modes: CCE's ψ₀, boosted into the
Cauchy frame by a Type-II transformation built purely from 3+1 quantities,
replaces the free physical boundary datum. Paper 2 (Moxon–Scheel–Teukolsky 2020)
provides the characteristic engine (worldtube data → hypersurface hierarchy →
ψ₀) that is agnostic to which Cauchy formulation feeds it. Therefore:

> **Z4c-CCM** = paper-1's boundary-condition set with its **physical slot**
> (the trace-free tangential-metric datum `h_AB^TF` of `eq:BCs_lastII`)
> driven by CCE's ψ₀ through paper-3's Choice-2 boost, while the **constraint
> slot** keeps Z4c CPBCs (`eq:general_CPBCs`) and the **gauge slot** keeps
> paper-1's gauge BCs — which simultaneously attacks paper-3's stated
> bottleneck (Sommerfeld gauge BCs causing spurious reflections).

Status markers follow `phys-agentic-loop/_common/contracts/markers.md`.
Source-grounded nodes cite verbatim equation labels resident in
`knowledge-database/paper_arxiv-<id>/nodes.jsonl` (P1 = 1010.0523v2,
P2 = 2007.01339, P3 = 2308.10361). N-nodes are the new-formulation work
program, mirrored in `knowledge-database/paper_z4c-CCM/nodes.jsonl`.

## The DAG

```mermaid
flowchart TD

subgraph SGA["Z4c Cauchy interior — P1"]
  A1["A1 Z4 covariant system\n(P1 eq:Z4_ADM_1, eq:Z4_ADM_2)"]
  A2["A2 Z4c conformal evolution,\ndamping kappa1, kappa2\n(P1 eq:Z4_decomp_first, eq:Theta_dot,\neq:def_theta, eq:def_Zi)"]
  A3["A3 puncture gauge + asymptotically\nharmonic shift (eps_alpha, eps_chi)\n(P1 eq:punc_alpha, eq:punc_beta,\neq:full-sec-z4-alpha, eq:full-sec-z4-beta)"]
  A4["A4 constraints H, M_i; Theta-Z subsystem;\nmonitor (P1 eq:ham-Z4, eq:mom-Z4,\neq:sys-ham-mom, eq:sys-theta-Z, eq:Cmonitor)"]
end

subgraph SGB["Z4c boundary characteristic structure — P1"]
  B1["B1 boundary frame, null vectors\nl-ring, k-ring, m-ring\n(P1 eq:backg-metric, eq:nullvector-k,\neq:nullvector-m)"]
  B2["B2 characteristic mode split:\n10 incoming = 4 constraint + 4 gauge\n+ 2 physical (P1 eq:2+1Z4_a,\neq:Z4_gammaAB_2+1, eq:theta _modes)"]
  B3["B3 constraint sector: order-L CPBCs\non Theta, Z_i (P1 eq:general_CPBCs,\neq:fcpbcTheta, eq:hoTheta-bc)"]
  B4["B4 gauge sector BCs: alpha, beta_perp,\nbeta_A (P1 eq:general_BCs_gauge_first,\neq:BCs-alpha, eq:BCs_lastII row 1)"]
  B5["B5 physical slot:\n(r^2 l-ring d)^(L+1) gamma_AB^TF = h_AB^TF\n(P1 eq:BCs_lastII row 2)"]
  B6["B6 Laplace-Fourier / Kreiss\nboundary-stability machinery\n(P1 eq:lapl-four, eq:bc_general,\nsol:lf, sol:lf2, eq:gen_estimate)"]
  B7["B7 reflection-coefficient diagnostic\n(P1 eq:refcoef)"]
end

subgraph SGC["Characteristic CCE engine — P2"]
  C1["C1 Bondi-Sachs / Bondi-like ansatz\n(P2 eq:BondiSachsMetric, eq:BondiLikeMetric)"]
  C2["C2 worldtube boundary data from\nCauchy metric (P2 Sec II.B: eq:auto-8,\neq:BondiLikeRadius)"]
  C3["C3 hypersurface hierarchy\nbeta->Q->U->W->H (P2 eq:HypersurfaceBeta,\neq:HypersurfaceQ, eq:HypersurfaceU,\neq:HypersurfaceW, eq:Habstract)"]
  C4["C4 partially flat gauge + regularity,\nno pure-gauge logs (P2 eq:UHat,\neq:Hregularity, eq:QWregularity)"]
  C5["C5 psi0 closed form on the\ncharacteristic grid (P2 eq:auto-48;\nP3 eq:psi0_CCE)"]
  C6["C6 scri+ outputs: news, Weyl scalars\n(P2 eq:BondiNews, eq:NewsDefinitionBondi,\neq:NewsDefinitionBondiLike)"]
end

subgraph SGD["Frame and psi0 matching chain — P3"]
  D1["D1 coordinate chain + Jacobians:\nCauchy -> null-radius -> Bondi-like -> PF\n(P3 eq:Jacobian_cauchy_null_radius,\neq:Jacobian_bondi_like_null_radius,\neq:ahat, eq:bhat, eq:omegahat)"]
  D2["D2 tetrad correspondence at worldtube:\nl-hat prop l-Cauchy (P3 eq:gh_tetrad_l,\neq:CCE_tetrad_l, eq:null_vector_CCE_rhat,\neq:l_GH_and_l_CCE_transform)"]
  D3["D3 Type I / Type II classification;\npsi0 scaling (P3 eq:lorentz_transformation_I,\neq:lorentz_transformation_II,\neq:lorentz_psi0_ii, typeI_Ahat_to_Ap)"]
  D4["D4 Choice-2 boost from 3+1 data:\nA = (alpha - gamma_ij beta^i s^j) e^(-2beta);\npsi0' = A^2 psi0 (P3 eq:auto-20, eq:auto-21)"]
  D5["D5 angular maps + time interpolation\nCCE grid <-> Cauchy grid\n(P3 eq:duhat_x, eq:du_xhat_cartesian)"]
end

subgraph SGE["GH-CCM reference implementation — P3"]
  E1["E1 GH FOSH + Bjorhus framework\n(P3 FOSH, Bjorhus_bc)"]
  E2["E2 physical projector + w- from psi0'\n(P3 eq:projection_operator,\neq:wab_projection, eq:wab, eq:GH_psi0_def)"]
  E3["E3 GH physical BC:\nw-|_BC <- CCM (P3 eq:bc_bjorhus)"]
  E4["E4 GH gauge BCs = Sommerfeld\n(P3 stated FUTURE bottleneck:\nspurious reflections); constraint BCs untouched"]
  E5["E5 diagnostics: Bondi-gauge violations,\nGH constraints (P3 eq:bondi_violation_psi3,\neq:im_psi2, eq:gauge_constraint,\neq:three_constraint)"]
end

subgraph SGN["NEW: Z4c-CCM formulation — mission 2 work program"]
  N1["N1 PRELIMINARY worldtube data from Z4c:\ng_munu(alpha, beta^i, chi, gamma~_ij)\n-> P2 Sec II.B input (variable map only)"]
  N2["N2 SOLID (verified iter 8) Z4c boundary tetrad\n<-> CCE tetrad: Choice-2 boost A holds verbatim;\nfactor sqrt(2)/2 exact; 4-GPU sweep residual 9e-16\n(scripts/verify_n2_boost.py)"]
  N3["N3 SOLID (verified iter 9) physical-mode\ninjection: (d_t + d_s)^2 gamma_AB^TF =\n4(psi0 mbar mbar + cc), psi0 = (1/4)(d_t+d_s)^2 h_mm,\nOFF-SHELL; CCE psi0' -> h_AB^TF datum of B5;\nL=0 Bjorhus form (scripts/verify_n3_dictionary.py)"]
  N4["N4 PRELIMINARY constraint sector:\nkeep order-L CPBCs (zero incoming\nTheta, Z modes); kappa1, kappa2 damping\nabsorbs residual injected violation"]
  N5["N5 HYPOTHESIS gauge sector: paper-1\ngauge BCs with asymptotically harmonic\nshift REPLACE paper-3 Sommerfeld\nplaceholder -> targets reflection bottleneck"]
  N6["N6 HYPOTHESIS Z4c-CCM composite scheme:\nfull 10-mode BC set + CCE coupling cadence\n(communication, interpolation, time stepping)"]
  N7["N7 OPEN well-posedness: LF/Kreiss boundary\nstability with inhomogeneous physical datum;\ncomposite system (weakly hyperbolic\nBondi-like sector caveat from P3)"]
  N8["N8 FUTURE GPU implementation (FINAL PRODUCT):\nJAX float64 code runnable on <=4 GPUs; tests\nTeukolsky wave, Kerr perturbation, pulse injection;\nR-coefficient + constraint monitor; every test <=10 min"]
  N9["N9 FUTURE damping-parameter interplay:\nkappa1, kappa2 vs injected-data error\n(Z4c analog of gamma2 c^3 term in eq:bc_bjorhus)"]
end

A1 --> A2
A2 --> A3
A2 --> A4
A2 --> B2
B1 --> B2
B2 --> B3
B2 --> B4
B2 --> B5
B6 --> B3
B6 --> B4
B3 --> B7

C1 --> C2
C1 --> C3
C2 --> C3
C3 --> C4
C3 --> C6
C4 --> C5
C4 --> C6

C1 --> D1
D1 --> D2
D2 --> D3
D3 --> D4
C5 --> D4

E1 --> E3
D4 --> E2
E2 --> E3
D5 --> E3
E3 --> E5
E4 --> E5

A2 --> N1
N1 --> C2
B1 --> N2
D4 --> N2
N2 --> N3
B5 --> N3
E2 --> N3
E3 --> N3
B3 --> N4
A4 --> N4
B4 --> N5
E4 --> N5
N3 --> N6
N4 --> N6
N5 --> N6
D5 --> N6
B6 --> N7
N6 --> N7
A2 --> N9
N6 --> N9
N7 --> N8
B7 --> N8
E5 --> N8
N9 --> N8

classDef p1 fill:#dbeafe,stroke:#1d4ed8,color:#111
classDef p2 fill:#dcfce7,stroke:#15803d,color:#111
classDef p3 fill:#fef3c7,stroke:#b45309,color:#111
classDef new fill:#fde2e2,stroke:#b91c1c,color:#111,stroke-width:2px
class A1,A2,A3,A4,B1,B2,B3,B4,B5,B6,B7 p1
class C1,C2,C3,C4,C5,C6 p2
class D1,D2,D3,D4,D5,E1,E2,E3,E4,E5 p3
class N1,N2,N3,N4,N5,N6,N7,N8,N9 new
```

## Node ledger (source grounding)

Every cited label below is a `node_id` in the corresponding
`knowledge-database/paper_arxiv-<id>/nodes.jsonl` (verified by
`scripts/zccm_dag_check.py`).

| node | paper | source equation labels (verbatim node_ids) | status |
|---|---|---|---|
| A1 | 1010.0523v2 | eq:Z4_ADM_1, eq:Z4_ADM_2 | [SOLID] transcription |
| A2 | 1010.0523v2 | eq:Z4_decomp_first, eq:Theta_dot, eq:def_theta, eq:def_Zi | [SOLID] transcription |
| A3 | 1010.0523v2 | eq:punc_alpha, eq:punc_beta, eq:full-sec-z4-alpha, eq:full-sec-z4-beta | [SOLID] transcription |
| A4 | 1010.0523v2 | eq:ham-Z4, eq:mom-Z4, eq:sys-ham-mom, eq:sys-theta-Z, eq:Cmonitor | [SOLID] transcription |
| B1 | 1010.0523v2 | eq:backg-metric, eq:nullvector-k, eq:nullvector-m | [SOLID] transcription |
| B2 | 1010.0523v2 | eq:2+1Z4_a, eq:Z4_gammaAB_2+1, eq:theta _modes | [SOLID] transcription |
| B3 | 1010.0523v2 | eq:general_CPBCs, eq:fcpbcTheta, eq:hoTheta-bc | [SOLID] transcription |
| B4 | 1010.0523v2 | eq:general_BCs_gauge_first, eq:BCs-alpha, eq:BCs_lastII | [SOLID] transcription |
| B5 | 1010.0523v2 | eq:BCs_lastII | [SOLID] transcription |
| B6 | 1010.0523v2 | eq:lapl-four, eq:bc_general, sol:lf, sol:lf2, eq:gen_estimate | [SOLID] transcription |
| B7 | 1010.0523v2 | eq:refcoef | [SOLID] transcription |
| C1 | 2007.01339 | eq:BondiSachsMetric, eq:BondiLikeMetric | [SOLID] transcription |
| C2 | 2007.01339 | eq:auto-8, eq:BondiLikeRadius | [SOLID] transcription |
| C3 | 2007.01339 | eq:HypersurfaceBeta, eq:HypersurfaceQ, eq:HypersurfaceU, eq:HypersurfaceW, eq:Habstract | [SOLID] transcription |
| C4 | 2007.01339 | eq:UHat, eq:Hregularity, eq:QWregularity | [SOLID] transcription |
| C5 | 2007.01339 + 2308.10361 | eq:auto-48 (P2); eq:psi0_CCE (P3) | [SOLID] transcription |
| C6 | 2007.01339 | eq:BondiNews, eq:NewsDefinitionBondi, eq:NewsDefinitionBondiLike | [SOLID] transcription |
| D1 | 2308.10361 | eq:Jacobian_cauchy_null_radius, eq:Jacobian_bondi_like_null_radius, eq:ahat, eq:bhat, eq:omegahat | [SOLID] transcription |
| D2 | 2308.10361 | eq:gh_tetrad_l, eq:CCE_tetrad_l, eq:null_vector_CCE_rhat, eq:l_GH_and_l_CCE_transform | [SOLID] transcription |
| D3 | 2308.10361 | eq:lorentz_transformation_I, eq:lorentz_transformation_II, eq:lorentz_psi0_ii, typeI_Ahat_to_Ap | [SOLID] transcription |
| D4 | 2308.10361 | eq:auto-20, eq:auto-21 | [SOLID] transcription |
| D5 | 2308.10361 | eq:duhat_x, eq:du_xhat_cartesian | [SOLID] transcription |
| E1 | 2308.10361 | FOSH, Bjorhus_bc | [SOLID] transcription |
| E2 | 2308.10361 | eq:projection_operator, eq:wab_projection, eq:wab, eq:GH_psi0_def | [SOLID] transcription |
| E3 | 2308.10361 | eq:bc_bjorhus | [SOLID] transcription |
| E4 | 2308.10361 | (design statement: gauge BCs stay Sommerfeld — P3 summary/obligations) | [SOLID] transcription |
| E5 | 2308.10361 | eq:bondi_violation_psi3, eq:im_psi2, eq:gauge_constraint, eq:three_constraint | [SOLID] transcription |
| N1 | z4c-CCM | new; consumes A2 variables, feeds C2 | [PRELIMINARY] |
| N2 | z4c-CCM | new; reuse of D4 with B1 frame; verifier scripts/verify_n2_boost.py → results/numerical/n2_boost_check.txt | [SOLID] verified |
| N3 | z4c-CCM | new; E2/E3 pattern into B5 slot; verifier scripts/verify_n3_dictionary.py → results/numerical/n3_dictionary_check.txt | [SOLID] verified |
| N4 | z4c-CCM | new; B3 + A4 reuse argument | [PRELIMINARY] |
| N5 | z4c-CCM | new; B4 replacing E4 | [HYPOTHESIS] |
| N6 | z4c-CCM | new; composite of N3, N4, N5 + D5 | [HYPOTHESIS] |
| N7 | z4c-CCM | new; B6 machinery, inhomogeneous data | [FUTURE] (analysis obligation) |
| N8 | z4c-CCM | new; P3 test suite + B7/E5/A4 diagnostics | [FUTURE] |
| N9 | z4c-CCM | new; A2 damping vs injected error | [FUTURE] |

## Why this formulation is new (delta over each paper)

1. **vs paper 1 (Z4c CPBC):** the physical datum `h_AB^TF` is no longer a local
   radiation-controlling guess (freezing-ψ₀ / hierarchical absorbing): it is the
   *actual* exterior solution's ψ₀ computed by CCE on the matched characteristic
   domain. Paper 1's boundary becomes transparent to nonlinear backscatter.
2. **vs paper 3 (GH-CCM):** the Cauchy formulation is Z4c (puncture-gauge NR
   codes: BAM/Einstein Toolkit lineage) instead of GH/SpECTRE, and both
   non-physical sectors improve: constraint sector gets Z4c's *damped*,
   light-speed constraint propagation with order-L CPBCs (vs untouched GH
   constraint BCs), and the gauge sector gets paper-1's well-posedness-analyzed
   gauge BCs (vs Sommerfeld — paper 3's own stated source of spurious
   reflections and its declared future work).
3. **vs paper 2 (CCE):** unchanged as an engine, but its worldtube input is now
   sourced from Z4c variables (N1) and its ψ₀ output is consumed at finite
   radius (via D4/N2) rather than only at scri+.

## Standing risks / obligations carried into the work program

- `[N7]` Composite well-posedness: paper 3 records the Bondi-like characteristic
  system as only weakly hyperbolic — CCM inherits this; paper 1's estimates are
  frozen-coefficient and (beyond the constraint subsystem) spherical-reduction
  only. The N7 analysis must state exactly what is and is not proven.
- `[N3]` Mode normalization: the identification of paper-1's incoming
  γ_AB^TF characteristic with w⁻ = 2(ψ₀′ m̄m̄ + ψ̄₀′ mm) needs an explicit
  second-order ↔ first-order dictionary (paper 1 is second order in space,
  paper 3's Bjørhus acts on FOSH fields); speeds and tetrad normalizations must
  be reconciled on the SAME background frame (B1).
- `[N2]` Paper-3's boost uses the characteristic β̂ on the worldtube; in the Z4c
  setup β̂ comes from C2/C3 with N1-supplied data — circular-dependency check
  needed at the initialization step (paper 3 resolves this for GH; verify order
  of operations survives the swap).
- `[N4]` CPBC compatibility: injected physical data must not source incoming
  constraint characteristics at corners/edges of the boundary (paper 1 discards
  tangential terms in 3D — its own flagged instability caveat).

## Verifier

`python3 scripts/zccm_dag_check.py` — checks (a) the mermaid graph above is
acyclic and every edge references a declared node; (b) every source label in
the node ledger exists as a `node_id` in the named paper's knowledge ledger;
(c) every N-node has a row in `knowledge-database/paper_z4c-CCM/nodes.jsonl`
whose status matches the ledger table. Output: `results/numerical/zccm_dag_check.txt`.
