# z4c-CMM — Cauchy-Characteristic Matching for the Z4c Evolution System

A verified-research repository: starting from three source papers, it derives,
**verifies, and implements a new formulation** — Cauchy-characteristic
matching (CCM) for the Z4c evolution system — under the
[phys-agentic-loop](phys-agentic-loop/) methodology (the agent proposes; the
verifier admits; every claim is ledgered with evidence).

## Source papers

| ID | Paper | Role |
|---|---|---|
| [arXiv:1010.0523](https://arxiv.org/abs/1010.0523) | Ruiz, Hilditch, Bernuzzi — *Constraint preserving boundary conditions for the Z4c formulation of GR* | Z4c boundary structure: 10-mode characteristic split, order-L CPBCs, gauge BCs, Laplace–Fourier/Kreiss machinery |
| [arXiv:2007.01339](https://arxiv.org/abs/2007.01339) | Moxon, Scheel, Teukolsky — *Improved Cauchy-characteristic evolution system* | CCE engine: worldtube data → hypersurface hierarchy → ψ₀ |
| [arXiv:2308.10361](https://arxiv.org/abs/2308.10361) | Ma et al. — *Fully relativistic 3D Cauchy-characteristic matching* | GH-CCM: ψ₀ → Type-II boost → physical-mode Bjørhus injection |

## The new formulation in one paragraph

At a timelike worldtube, Z4c has exactly **ten incoming characteristic
modes: 2 physical + 4 constraint + 4 gauge**. Z4c-CCM supplies them as
follows. The **2 physical** modes are driven by the *exact* incoming Weyl
object `U⁻_AB = TT₂[E + ε(s)B]_AB` toward the CCE-computed value
`−(ψ₀′ m̄_A m̄_B + c.c.)`, where `ψ₀′ = A²ψ₀^CCE` with the Type-II boost
`A = (α − γ_ij β^i s^j)e^{−2β̂}` (built from 3+1 quantities only — verified to
transplant from GH to Z4c verbatim), and `E_ij = R_ij + KK_ij − K_ikK^k_j`,
`B_ij = ε_i^{kl}D_kK_{lj}` realized from the Z4c state via Gauss–Codazzi
(verified nonperturbatively on Schwarzschild). The **4 constraint** modes
keep paper-1's order-L CPBCs unchanged — the CCM channel is provably
constraint-orthogonal. The **4 gauge** modes use paper-1's analyzed gauge BCs
in place of paper-3's Sommerfeld placeholder, removing its spurious-reflection
bottleneck (`R_Sommerfeld = −β̊(1−β̊)/((1+β̊)(2−β̊))` vs `R = 0` for the
boundary-adapted operator — proven analytically and reproduced dynamically to
float64 precision). With `ψ₀ → 0` the scheme reduces exactly to paper-1's
CPBC scheme: **Z4c-CCM ⊃ Z4c-CPBC**.

Full statement: [`reformulate/z4c-CMM/synthesis_z4c-CCM/z4c_ccm_boundary_conditions.md`](reformulate/z4c-CMM/synthesis_z4c-CCM/z4c_ccm_boundary_conditions.md) ·
roadmap DAG: [`formulation_dag.md`](reformulate/z4c-CMM/synthesis_z4c-CCM/formulation_dag.md) ·
well-posedness sketch: [`n7_wellposedness_sketch.md`](reformulate/z4c-CMM/synthesis_z4c-CCM/n7_wellposedness_sketch.md) ·
**method-only paper** (new content in dark blue): [`paper/z4c-CMM/zccm_formulation/`](paper/z4c-CMM/zccm_formulation/).

## Verification matrix

Every claim has a runnable verifier; every number below is reproduced by the
listed command (outputs certified under `results/numerical/`). All GPU runs:
JAX float64, ≤ 4× A100 (`CUDA_VISIBLE_DEVICES=0,1,2,3`), every run ≤ 60 s
(budget 10 min).

| Node | Claim | Verifier | Key numbers |
|---|---|---|---|
| Mission 1 | all 325 equations of the 3 papers as per-paper DAGs | `python3 scripts/eqdag_check.py` | acyclicity, edge closure, label coverage 80/80 + 48/48 + 74/74 |
| ONE DAG | single connected formulation DAG (36 nodes / 53 edges) | `python3 scripts/zccm_dag_check.py` | acyclic, connected, source-grounded, ledger-synced |
| N1 | Z4c ↔ ADM ↔ 4-metric maps; kinematic closure ⇒ `D_iβ^i` is the full covariant divergence | `python3 scripts/verify_n1_varmap.py` | exact (full-symbol); naive reading leaves residual `γ_ij β^k∂_k lnχ` (c = 1) |
| N2 | Type-II boost valid verbatim on the Z4c frame | `python3 scripts/verify_n2_boost.py` | exact at 12 rational pts; factor √2/2; GPU sweep 9.3e-16 over 4.2M samples |
| N3 | **exact dictionary** `(E+ε(s)B)(m,m) = −ψ₀` (partner `−ψ̄₄`); Gauss–Codazzi realization | `python3 scripts/verify_n3_exact.py` | generic-Weyl exact; Schwarzschild two-route 7.7e-17 / 3.9e-17, identity 8.5e-16 |
| N3 lim | linearized on-shell corollary `(∂_t+∂_s)²γ^TF = 4(ψ₀m̄m̄+cc)` | `python3 scripts/verify_n3_dictionary.py` | exact off-shell; GPU residual 0.0 |
| N4 | CCM channel constraint-orthogonal ⇒ CPBCs retained | `python3 scripts/verify_n4_cpbc_compat.py` | `H = M_i = 0` identically; negative controls pass |
| N5 | gauge sector decoupled; Sommerfeld reflects, P1 operator doesn't | `python3 scripts/verify_n5_gauge.py` | `R_Som = −b(1−b)/((1+b)(2−b))` exact; `R_P1 = 0` |
| N6 | composite 10-mode BC table; ψ₀→0 reduction; wire-through | `python3 scripts/verify_n6_composite.py` | reduction exact; chain residual 1.5e-16 over 1M states |
| N7 | Kreiss condition datum-independent: `D_L(s) = (2s)^{L+1}` | `python3 scripts/verify_n7_lf_sketch.py` | exact L = 0..3; uniform ratio 2^{L+1} |
| N8 | package + dynamical model at float64 limits | `python3 packages/zccm/tests/test_z4c_vars.py` · `test_model1d.py` | Sommerfeld L2 ratio dev **2.1e-14** (12 digits vs analytic); absorption **1.9e-14**; injection **2.4e-11** (roundoff-accumulation optimum) |

## Repository layout

```
phys-agentic-loop/      vendored methodology (submodule; pipelines, contracts, ledgers)
ref-paper/              imported LaTeX sources + sha256 PROVENANCE (gitignored, local-only)
reformulate/z4c-CMM/
  paper_<id>/           stage-1 decomposition: 10 artifacts per paper
                        (convention/derivation/logic/claims/obligations/...)
  synthesis_z4c-CCM/    THE FORMULATION: one DAG, BC table, well-posedness sketch
knowledge-database/     per-paper equation DAGs (325 nodes) + z4c-CCM program ledger
error-database/         every trial: 18 rows, each failure with root cause + fix
results/numerical/      certified verifier outputs        figures/  rendered DAG
packages/zccm/          FINAL PRODUCT: JAX float64 GPU implementation + tests
paper/z4c-CMM/zccm_formulation/   method-only PRD paper (REVTeX; dark blue = new)
scripts/                verifiers (one per DAG node) + DAG checkers + renderers
progress/               three-timescale mission notes (iteration/nodal/project)
```

## AthenaK production implementation (`athenak/`)

The formulation is implemented in a vendored, SHA-pinned copy of
[AthenaK](https://github.com/IAS-Astrophysics/athenak) (arXiv:2409.10383;
`athenak/PROVENANCE.md`), tracked in this repository by design:

- `athenak/src/z4c/z4c_ccm.hpp` — device-inline (Kokkos, CPU+GPU) CCM
  injection: the Bjorhus-style physical-mode datum
  `rhs(Ã_ab) += −(2χ/α)(ψ₀′ m̄_a m̄_b + c.c.)^TT`, with the Type-II boost
  `ψ₀′ = A²ψ₀`, `A = (α − γ_ij β^i s^j)e^{−2β̂}`, evaluated from the local Z4c
  state; applied inside the existing `Z4c_SomBC` task at all 12
  outer-boundary call sites (`z4c_Sbc.cpp`). `ccm_amp = 0` reduces exactly to
  the stock Sommerfeld scheme.
- Options (`<z4c>` input block): `ccm, ccm_amp, ccm_t0, ccm_sigma,
  ccm_betahat`. The analytic Gaussian ψ₀ provider exercises the full
  machinery; a CCE-coupled provider replaces one function
  (`Z4cCCMDatumPsi0`) — the `src/z4c/cce/` worldtube infrastructure is the
  documented coupling hook.
- Test battery (`scripts/test_athenak_ccm.sh` + `check_athenak_ccm.py`,
  input `athenak/inputs/z4c/ccm/z4c_ccm_test.athinput`): T1 reduction
  (Minkowski preserved with zero datum), T2 causal injection through the
  boundary, T3 linear response. Runs on CPU (OpenMP) and GPU (CUDA/A100)
  builds; outputs under `results/numerical/athenak_ccm/`.
- Build note: with CUDA 13 the vendored kokkos `nvcc_wrapper` default arch is
  patched to `sm_80`; the snapshot needs a CUDA-12.x toolkit
  (`PATH=/usr/local/cuda-12.9/bin:$PATH`) — both ledgered.

## The GPU package (`packages/zccm/`)

```python
import sys; sys.path.insert(0, "packages/zccm")
from zccm import Z4cState, four_metric, boost_factor, psi0_to_u_minus, u_minus_scalar
from zccm.boundary import bc_targets          # the 10-mode BC table
from zccm import model1d                      # dynamical model problem (order 2/4/6 FD)
```

- `z4c_vars` — variable maps (N1), Type-II boost (N2), exact `U⁻` datum
  constructor + `U⁻(m,m)` evaluator from 3+1 data (N3-exact).
- `boundary` — `BCTargets` for all 10 incoming modes; physical rows exact.
- `model1d` — boundary wave model: Fornberg stencils (order 6 interior with a
  GKS-stable mixed closure), Kreiss–Oliger hook, characteristic/Bjørhus BC
  enforcement; the reflection/injection testbed behind Table II of the paper.

## Reproduce everything

```bash
# environment (8x A100 host; uses at most 4 GPUs)
pip install --user --break-system-packages nvidia-cudnn-cu12   # enables JAX cuda12

# verification suite (each ≤ 60 s)
for v in scripts/verify_*.py scripts/eqdag_check.py scripts/zccm_dag_check.py; do
  CUDA_VISIBLE_DEVICES=0,1,2,3 timeout 600 python3 -u "$v" || echo "FAIL: $v"; done
CUDA_VISIBLE_DEVICES=0,1,2,3 timeout 600 python3 -u packages/zccm/tests/test_z4c_vars.py
CUDA_VISIBLE_DEVICES=0,1,2,3 timeout 600 python3 -u packages/zccm/tests/test_model1d.py

# paper
cd paper/z4c-CMM/zccm_formulation && latexmk -pdf main.tex
```

Note: `scripts/eqdag_check.py` needs the local `ref-paper/` imports
(stage-0 fast path: `curl https://arxiv.org/e-print/<id>`, unpack to
`ref-paper/arxiv-<id>/src/`, see any `PROVENANCE.md`).

## Provenance discipline

- **Knowledge ledger** (`knowledge-database/`): one row per admitted node —
  claim, evidence path, status, dependency edges. Equation nodes cite verbatim
  tex labels; DAG acyclicity is machine-checked.
- **Error ledger** (`error-database/paper_z4c-CCM/trials.jsonl`): all 18
  trials, including every failure with its root cause — among them: a
  kernel-timeout that forced the exact-sampling + GPU-sweep verification
  strategy; a null-leg convention flip caught by the checker; the
  characteristic-orientation mirror in the model problem; the GKS instability
  of one-sided 6th-order closures; and the user-review finding that the
  original dictionary was perturbative (fixed by the exact Weyl-object
  formulation in iter 16).
- **Commit history**: one commit per verified substage (`iter N [node]: ...`),
  each embedding its verbatim verifier transcript.

## Open obligations (successor work; human sign-off requested)

1. Coupled Cauchy↔characteristic well-posedness (Bondi-like sector is only
   weakly hyperbolic — open even for GH-CCM); oblique incidence + curvature;
   corners (`n7_wellposedness_sketch.md`, items 1–4).
2. N9: nonlinear constraint-damping addition to the physical Bjørhus form
   (κ₁/κ₂-weighted; Z4c analog of GH's γ₂c³ term).
3. Full-GR 3D tests (Teukolsky wave, Kerr perturbation, pulse injection) with
   a production Z4c code coupled to a SpECTRE-style CCE.
4. Optional: Kahan-compensated time accumulation to push the model-problem
   injection floor below the float64 steps×ε optimum (2.4e-11).

## Known repository caveat

`phys-agentic-loop/alignment.md` was replaced upstream (its commit `663804c`)
by an unrelated persona file, dangling the §15–§21 cross-references; the
binding scientific rules survive in `phys-agentic-loop/.claude/settings.json`
and `_common/contracts/`. Restorable from upstream commit `aa051c0`.
