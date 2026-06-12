# zccm — GPU implementation of the Z4c-CCM boundary scheme

Mission-2 final product of the z4c-CMM project: a JAX float64 implementation
of the new Z4c Cauchy-characteristic-matching formulation
(`reformulate/z4c-CMM/synthesis_z4c-CCM/`), runnable on up to 4 GPUs.

## Modules (mapped to formulation-DAG nodes)

| module | nodes | content |
|---|---|---|
| `zccm.z4c_vars` | N1, N2, N3 | `Z4cState`, `to_adm`/`from_adm`, `four_metric(+inverse)` (worldtube data for CCE), `boost_factor` (Type-II boost, verified), `psi0_to_physical_datum` (w⁻ constructor, verified) |
| `zccm.boundary` | N6 | the 10-incoming-mode BC table: `bc_targets` (2 physical ← ψ₀, 4 constraint ← CPBC, 4 gauge ← P1), `physical_residual`, mode-count bookkeeping |
| `zccm.model1d` | N8 | dynamical model problem (P1 eq:eomgammasA): characteristic/Bjørhus BC enforcement, `char_decompose`, `evolve` (FD2 + RK4, jit/scan) |

## Tests (each ≤ 10 min, ≤ 4 GPUs; outputs in `results/numerical/`)

- `tests/test_z4c_vars.py` — exact maps, round trip (4e-16 over 1M states),
  boost and datum constructors.
- `tests/test_model1d.py` — dynamical reflection/injection: Sommerfeld
  reproduces the analytic q-ratio b/(2+b) to 0.03%; the P1/CCM operator
  absorbs at machine precision (R ~ 3e-15); CCM injection reproduces an
  exact incoming pulse to 1e-4 (two-way transparency).

## Usage

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 timeout 600 python3 packages/zccm/tests/test_z4c_vars.py
CUDA_VISIBLE_DEVICES=0,1,2,3 timeout 600 python3 packages/zccm/tests/test_model1d.py
```

```python
import sys; sys.path.insert(0, "packages/zccm")
from zccm import Z4cState, four_metric, boost_factor, psi0_to_physical_datum
from zccm.boundary import bc_targets
# state = Z4cState(...); psi0 from the CCE module of choice
# targets = bc_targets(state, psi0_cce, s_unit, beta_char)
```

Orientation note: `model1d` documents its characteristic-speed convention
explicitly (outgoing speed 1−b for the chosen sign of the shift); the
Sommerfeld reflection formula appears as the b → −b mirror of the N5 form —
both vanish iff b = 0, and the boundary-adapted operator absorbs in either
orientation.
