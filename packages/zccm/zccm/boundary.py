"""boundary — the Z4c-CCM composite boundary-condition table (DAG node N6).

Implements the 10-incoming-mode BC table of
reformulate/z4c-CMM/synthesis_z4c-CCM/z4c_ccm_boundary_conditions.md:
  rows 1-2  physical    (l.d)^2 gamma_AB^TF = 2(psi0' mbar mbar + cc)
                        with psi0' = A^2 psi0_CCE (N2 boost, N3 dictionary)
  rows 3-6  constraint  (r^2 l.d)^(L+1) {Theta, Z_s, Z_A} = 0   (P1 CPBCs, N4)
  rows 7-10 gauge       (r^2 X.d)^(L+1) {alpha, beta_s, beta_A} = h (P1, N5)

This module evaluates BC *targets* and *residuals* (the datum chain and the
characteristic bookkeeping); time integration belongs to the host evolution
code / the N8 model tests. Batched, float64, jit-safe.
"""
from typing import NamedTuple

import jax.numpy as jnp

from .z4c_vars import Z4cState, boost_factor, psi0_to_u_minus


class BCTargets(NamedTuple):
    """Targets for the 10 incoming modes (leading batch shape ...).

    physical : (..., 2, 2) TT matrix — EXACT target for the incoming Weyl
               object U^-_AB = TT2[E + eps(s)B]_AB (rows 1-2, CCM datum)
    theta    : (...,)      target for (r^2 l.d)^(L+1) Theta      (row 3)
    z_s      : (...,)      target for (r^2 l.d)^(L+1) Z_s        (row 4)
    z_A      : (..., 2)    target for (r^2 l.d)^(L+1) Z_A        (rows 5-6)
    lapse    : (...,)      datum h_alpha                         (row 7)
    shift_s  : (...,)      datum h_s                             (row 8)
    shift_A  : (..., 2)    datum h_A                             (rows 9-10)
    """
    physical: jnp.ndarray
    theta: jnp.ndarray
    z_s: jnp.ndarray
    z_A: jnp.ndarray
    lapse: jnp.ndarray
    shift_s: jnp.ndarray
    shift_A: jnp.ndarray


def n_incoming_modes(t: BCTargets) -> int:
    """Mode-count bookkeeping: physical TT carries 2 dof (symmetric traceless
    2x2), z_A and shift_A carry 2 each, scalars 1 each => must equal 10."""
    return 2 + 1 + 1 + t.z_A.shape[-1] + 1 + 1 + t.shift_A.shape[-1]


def physical_target(psi0_cce, state: Z4cState, s_unit, beta_char):
    """Rows 1-2 (EXACT): U^-_AB|_BC = -(psi0' mbar mbar + cc),
    psi0' = A^2 psi0_CCE, where U^-_AB = TT2[E + eps(s)B]_AB is the incoming
    Weyl object built from the Z4c state via Gauss-Codazzi
    (E = R + KK - KK, B = eps DK; scripts/verify_n3_exact.py — exact,
    nonperturbative). Chains N2 (boost_factor) and N3-exact
    (psi0_to_u_minus). psi0_cce -> 0 recovers paper-1's freezing-psi0
    absorbing scheme (verified limit).
    """
    A = boost_factor(state, s_unit, beta_char)
    return psi0_to_u_minus(A**2 * psi0_cce)


def bc_targets(state: Z4cState, psi0_cce, s_unit, beta_char,
               h_alpha=None, h_s=None, h_A=None) -> BCTargets:
    """Assemble the full 10-mode target table (CPBC rows are homogeneous;
    gauge data default to 0 per the formulation document)."""
    B = state.chi.shape
    zero = jnp.zeros(B, state.gamt.dtype)
    zero2 = jnp.zeros(B + (2,), state.gamt.dtype)
    return BCTargets(
        physical=physical_target(psi0_cce, state, s_unit, beta_char),
        theta=zero, z_s=zero, z_A=zero2,
        lapse=zero if h_alpha is None else h_alpha,
        shift_s=zero if h_s is None else h_s,
        shift_A=zero2 if h_A is None else h_A)


def physical_residual(u_minus_measured, targets: BCTargets):
    """Residual of rows 1-2 (EXACT): U^-_AB(E, B) - U^-_AB|_BC.

    u_minus_measured: (..., 2, 2) — the incoming Weyl object
    TT2[E + eps(s)B]_AB evaluated from the Z4c state and its spatial
    derivatives at the boundary (helpers: zccm.u_minus_scalar). Zero for any
    exact solution whose psi0 matches the injected datum (wire-through;
    nonperturbative). In the linearized flat limit this reduces to the
    superseded operator form (d_t+d_s)^2 gamma^TF = 4(psi0' mbm+cc) via
    U^- -> -(1/4)(d_t+d_s)^2 gamma^TF (on shell).
    """
    return u_minus_measured - targets.physical
