# z4c-CCM formulation mission — Nodal Note

## 10-iter window

Iterations 7–12 (mission 2 so far). Knowledge ledger delta: S2-synthesis-dag,
N1, N2, N3, N4, N5 all → solid (6 admissions). Error ledger: 9 trial rows —
4 fail (N2 kernel_timeout; N3 null-leg flip; N1 assertion sign + jit
concretization; N4 negative-control taxonomy) each with root cause, 5 pass.
Node coverage: N1–N5 of the N1–N6 completion gate done; N6 (composite), N7
(sketch), N8 (GPU final product), N9 remain. Strategic redirects: user
constraints landed iter 8 (10-min tests, ≤4 GPUs, GPU code as final product,
commit-per-substage) — verification strategy switched to exact rational/point
sampling + JAX GPU sweeps; no parameter-tweak loops (every fail fixed by a
structural correction).

## Logic-DAG snapshot

ONE connected DAG, 36 nodes / 53 edges (zccm_dag_check: all PASS).
N-sector: N1 ✓ N2 ✓ N3 ✓ N4 ✓ N5 ✓ · N6 [HYPOTHESIS] next · N7/N8/N9 [FUTURE].
External deps: none beyond the three decomposed papers.

## Accepted-results snapshot

Mirrors RESEARCH_STATE Accepted Results Log: 7 rows solid (DAG, N2, N3, N1,
N4, N5 + mission-1 imports). Verifier outputs under results/numerical/
(n2_boost, n3_dictionary, n1_varmap{,_gpu}, n4_cpbc_compat, n5_gauge,
zccm_dag_check).

## Simplification cycle

None consumed; zccm package modules remain single-purpose. Watch: boundary
module (N6) should reuse z4c_vars constructors, not duplicate.

## Failure-mode drift

`branch_cut_ambiguity` reused for the N3 null-leg flip (closest existing tag;
arguably a convention-flip tag would fit better — flag for taxonomy review,
no enum extension yet).
