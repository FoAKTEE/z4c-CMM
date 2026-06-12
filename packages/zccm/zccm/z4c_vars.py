"""z4c_vars — Z4c variable maps for the worldtube (formulation-DAG node N1),
plus the N2 boost factor and the N3 physical-datum constructor.

Conventions (paper 1, arXiv:1010.0523v2, eq:Z4_decomp_first group):
  chi      = (det gamma)^(-1/3)          conformal factor, det(gamt) = 1
  gamt_ij  = chi * gamma_ij              conformal metric
  K        = Khat + 2*Theta              trace of extrinsic curvature
  At_ij    = chi * (K_ij - gamma_ij K/3) trace-free conformal extr. curvature
  g_00 = -alpha^2 + gamma_ij beta^i beta^j,  g_0i = gamma_ij beta^j,
  g_ij = gamma_ij.

All functions are batched over leading axes, float64, jit-safe. Index pairs
(i, j) are the trailing two axes of rank-2 fields.
"""
from typing import NamedTuple

import jax.numpy as jnp


class Z4cState(NamedTuple):
    """Z4c evolved variables at points (any leading batch shape)."""
    chi: jnp.ndarray      # (...,)
    gamt: jnp.ndarray     # (..., 3, 3)  det = 1
    At: jnp.ndarray       # (..., 3, 3)  gamt-trace-free
    Khat: jnp.ndarray     # (...,)
    Theta: jnp.ndarray    # (...,)
    alpha: jnp.ndarray    # (...,)
    beta: jnp.ndarray     # (..., 3)


def to_adm(s: Z4cState):
    """Z4c state -> (gamma_ij, K_ij, alpha, beta^i)."""
    gamma = s.gamt / s.chi[..., None, None]
    K = s.Khat + 2.0 * s.Theta
    K_ij = (s.At + s.gamt * (K / 3.0)[..., None, None]) / s.chi[..., None, None]
    return gamma, K_ij, s.alpha, s.beta


def from_adm(gamma, K_ij, alpha, beta, Theta=None):
    """(gamma_ij, K_ij, alpha, beta^i)[, Theta] -> Z4c state.

    Theta defaults to zero (constraint-satisfying data); Khat = K - 2 Theta.
    """
    if Theta is None:
        Theta = jnp.zeros(gamma.shape[:-2], gamma.dtype)
    chi = jnp.linalg.det(gamma) ** (-1.0 / 3.0)
    gamt = gamma * chi[..., None, None]
    gamma_inv = jnp.linalg.inv(gamma)
    K = jnp.einsum("...ij,...ij->...", gamma_inv, K_ij)
    At = (K_ij - gamma * (K / 3.0)[..., None, None]) * chi[..., None, None]
    return Z4cState(chi=chi, gamt=gamt, At=At, Khat=K - 2.0 * Theta,
                    Theta=Theta, alpha=alpha, beta=beta)


def four_metric(s: Z4cState):
    """Z4c state -> g_{mu nu} (..., 4, 4) on the worldtube (CCE input)."""
    gamma, _, alpha, beta = to_adm(s)
    beta_dn = jnp.einsum("...ij,...j->...i", gamma, beta)
    g = jnp.zeros(s.chi.shape + (4, 4), s.gamt.dtype)
    g = g.at[..., 0, 0].set(-alpha**2 + jnp.einsum("...i,...i->...", beta_dn, beta))
    g = g.at[..., 0, 1:].set(beta_dn)
    g = g.at[..., 1:, 0].set(beta_dn)
    g = g.at[..., 1:, 1:].set(gamma)
    return g


def four_metric_inverse(s: Z4cState):
    """Z4c state -> g^{mu nu} via the closed-form ADM inverse."""
    gamma, _, alpha, beta = to_adm(s)
    gamma_inv = jnp.linalg.inv(gamma)
    gi = jnp.zeros(s.chi.shape + (4, 4), s.gamt.dtype)
    a2 = alpha**2
    gi = gi.at[..., 0, 0].set(-1.0 / a2)
    gi = gi.at[..., 0, 1:].set(beta / a2[..., None])
    gi = gi.at[..., 1:, 0].set(beta / a2[..., None])
    gi = gi.at[..., 1:, 1:].set(gamma_inv - jnp.einsum("...i,...j->...ij", beta, beta)
                                / a2[..., None, None])
    return gi


def boost_factor(s: Z4cState, s_unit, beta_char):
    """N2 (verified): Type-II boost A = (alpha - gamma_ij beta^i s^j) e^{-2 beta_char}.

    s_unit: outward unit normal s^i (..., 3), unit w.r.t. gamma_ij.
    beta_char: characteristic beta-hat scalar on the worldtube (CCE variable).
    psi0_cauchy = A^2 * psi0_characteristic (paper 3 eq:auto-21; N2 verifier:
    scripts/verify_n2_boost.py, results/numerical/n2_boost_check.txt).
    """
    gamma, _, alpha, beta = to_adm(s)
    bs = jnp.einsum("...ij,...i,...j->...", gamma, beta, s_unit)
    return (alpha - bs) * jnp.exp(-2.0 * beta_char)


def psi0_to_physical_datum(psi0_prime):
    """N3 (verified): physical boundary datum from the boosted psi0'.

    Returns the TT 2x2 tangential matrix
        w-_AB = 4 (psi0' mbar_A mbar_B + conj(psi0') m_A m_B),
    with the tangential dyad m_A = (1, i)/sqrt(2) in the (e_y, e_z) basis of
    the boundary-adapted frame; equals (d_t + d_s)^2 gamma_AB^TF, the datum
    driving paper-1's eq:BCs_lastII physical slot (L=0 Bjorhus form).
    Verifier: scripts/verify_n3_dictionary.py (identity exact, off-shell).
    """
    # mbar⊗mbar = 0.5*[[1, -i], [-i, -1]] => (psi0' mbm + cc)_11 = Re psi0',
    # _12 = Im psi0'; with the verified K = 4:
    re, im = jnp.real(psi0_prime), jnp.imag(psi0_prime)
    w11 = 4.0 * re
    w12 = 4.0 * im
    return jnp.stack([jnp.stack([w11, w12], -1),
                      jnp.stack([w12, -w11], -1)], -2)
