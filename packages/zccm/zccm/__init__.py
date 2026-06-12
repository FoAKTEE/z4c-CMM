"""zccm — GPU-runnable (JAX, float64) implementation of the Z4c-CCM scheme.

Mission-2 final product (z4c-CMM repo). Modules map to formulation-DAG nodes:
  z4c_vars  — N1 variable maps (Z4c state <-> ADM <-> 4-metric) +
              N2 Type-II boost factor + N3 physical-datum constructor.

Constraints (binding): float64 everywhere; runs on up to 4 GPUs
(CUDA_VISIBLE_DEVICES=0,1,2,3); every test <= 10 minutes.
"""
import jax

jax.config.update("jax_enable_x64", True)

from .z4c_vars import (  # noqa: E402,F401
    Z4cState, to_adm, from_adm, four_metric, four_metric_inverse,
    boost_factor, psi0_to_physical_datum,
)
