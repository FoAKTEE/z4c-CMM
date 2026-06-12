"""model1d — dynamical model problem for the Z4c-CCM boundary scheme (N8).

The physical-sector boundary wave problem of paper 1 (eq:eomgammasA at normal
incidence): [d_t^2 - 2 b d_t d_x - (1 - b^2) d_x^2] u = 0 on x in [-X, 0],
boundary at x = 0, shift b. First-order reduction V = (u, Pi = d_t u,
Phi = d_x u):

    d_t u = Pi,   d_t Phi = d_x Pi,   d_t Pi = 2 b d_x Pi + (1 - b^2) d_x Phi.

The operator factors as (d_t + (1-b) d_x)(d_t - (1+b) d_x): with this sign
convention the OUTGOING direction (toward x = 0) has speed c_o = 1 - b and the
INCOMING direction speed c_i = 1 + b (left eigenvectors of the flux matrix
A = [[2b, 1-b^2], [1, 0]]: (1, 1-b) with speed 1+b moving LEFT = incoming;
(1, -(1+b)) with speed 1-b moving RIGHT = outgoing). Characteristic fields:

    q_in  = Pi + (1 - b) Phi   (vanishes on outgoing; the BC target)
    q_out = Pi - (1 + b) Phi   (vanishes on incoming; preserved by the BC)

Boundary conditions at x = 0 (characteristic/Bjorhus enforcement):

    sommerfeld : Pi + Phi = 0  <=>  q_in = q_out * b/(2 + b)   (reflects)
    p1         : q_in = 0                       (boundary-adapted; absorbs)
    ccm        : q_in = g(t)                    (CCM-injected datum; rows 1-2)

Analytic Sommerfeld reflection in THIS orientation is the b -> -b mirror of
the N5 formula: R = b(1+b)/((1-b)(2+b)) (both forms vanish iff b = 0).
For an exact incoming solution u = G(t + x/(1+b)) the CCM datum is
g(t) = q_in|_{x=0} = (1 + (1-b)/(1+b)) G'(t) = 2 G'(t)/(1+b).

Discretization: centered FD interior of selectable order (2/4/6; Fornberg
one-sided closures at the edges), RK4; BC enforced by characteristic
(Bjorhus-style) replacement at the boundary nodes after each substage.
float64, jit-safe, batched-free (single config per call — tests fan configs
across devices).
"""
from functools import partial

import numpy as _np
import jax
import jax.numpy as jnp


def _fornberg(z, x):
    """Fornberg weights for the 1st derivative at z from nodes x."""
    n = len(x)
    c = _np.zeros((n, 2))
    c1, c4 = 1.0, x[0] - z
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, 1)
        c2, c5 = 1.0, c4
        c4 = x[i] - z
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1]
                                    - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, 1]


_CENTERED = {2: _np.array([-0.5, 0.0, 0.5]),
             4: _np.array([1, -8, 0, 8, -1]) / 12.0,
             6: _np.array([-1, 9, -45, 0, 45, -9, 1]) / 60.0}


def _edge_weights(order):
    """Edge-row weights (rows i = 0..w-1, each over the first `win` nodes).

    Orders 2/4: one-sided/biased Fornberg rows on (order+1)-node windows.
    Order 6: GKS-stable MIXED closure — fully one-sided 6th-order rows are
    unstable under RK4 (ledgered), so rows 0..2 use 4th-order 5-node
    windows (one-sided, biased, centered), an O(h^4) closure localized at
    3 nodes (global accuracy ~ 5th order, Gustafsson).
    """
    if order < 6:
        w = order // 2
        nodes = _np.arange(order + 1, dtype=float)
        rows = [_fornberg(float(i), nodes) for i in range(w)]
        win = order + 1
    else:
        nodes = _np.arange(5, dtype=float)
        rows = [_fornberg(float(i), nodes) for i in range(3)]
        win = 5
    W = _np.zeros((len(rows), win))
    for i, r in enumerate(rows):
        W[i, :] = r
    return W, win


def _dx(f, h, order=2):
    w = order // 2
    cen = jnp.asarray(_CENTERED[order])
    interior = sum(cen[k + w] * f[w + k: f.shape[0] - w + k or None]
                   for k in range(-w, w + 1))
    d = jnp.zeros_like(f)
    d = d.at[w:-w].set(interior)
    ew_np, win = _edge_weights(order)
    ew = jnp.asarray(ew_np)
    d = d.at[:w].set(ew @ f[:win])
    d = d.at[-w:].set(-(ew @ f[-1: -(win + 1): -1])[::-1] * 1.0)
    return d / h


_KO = {6: _np.array([1, -6, 15, -20, 15, -6, 1]) / 64.0,
       8: _np.array([-1, 8, -28, 56, -70, 56, -28, 8, -1]) / 256.0}
# signs chosen so the operator DAMPS the highest mode: Fourier symbol of the
# applied term is -(sigma/h) sin^(2p)(kh/2) for both rows above.


def _ko(f, h, sigma, p2):
    """Kreiss-Oliger dissipation of order p2 (6 or 8): an O(h^(p2-1)) term,
    below the truncation of the matching interior stencil; applied only where
    the centered stencil fits (tapered to zero at the edges)."""
    if sigma == 0.0:
        return jnp.zeros_like(f)
    w = p2 // 2
    cko = jnp.asarray(_KO[p2])
    interior = sum(cko[k + w] * f[w + k: f.shape[0] - w + k or None]
                   for k in range(-w, w + 1))
    d = jnp.zeros_like(f)
    return d.at[w:-w].set(sigma / h * interior)


def _rhs(state, b, h, order, ko_sigma=0.0):
    u, Pi, Phi = state
    dPi, dPhi = _dx(Pi, h, order), _dx(Phi, h, order)
    p2 = 8 if order >= 6 else 6
    return (Pi + _ko(u, h, ko_sigma, p2),
            2 * b * dPi + (1 - b**2) * dPhi + _ko(Pi, h, ko_sigma, p2),
            dPi + _ko(Phi, h, ko_sigma, p2))


def _apply_bc(state, b, kind, gval):
    """Characteristic (Bjorhus-style) BC enforcement: PRESERVE the outgoing
    characteristic from the interior, SET the incoming one, reconstruct
    (Pi, Phi) from the pair. q_out = Pi - (1-b) Phi, q_in = Pi + (1+b) Phi
    => Phi = (q_in - q_out)/2, Pi = q_out + (1-b)(q_in - q_out)/2.

    Right boundary (index -1):
      p1 / ccm  : q_in = gval (0 for p1)            [rows 1-2 analog]
      sommerfeld: Pi + Phi = 0 => q_in = -q_out b/(2-b)  (reflective; its
                  u-amplitude ratio is exactly the N5 analytic R)
    Left edge (index 0): no inflow — q_out = 0, q_in preserved.
    """
    u, Pi, Phi = state
    qo = Pi[-1] - (1 + b) * Phi[-1]
    qi = jnp.where(kind == 0, qo * b / (2 + b),
                   jnp.where(kind == 2, gval, 0.0))
    # q_in - q_out = 2 Phi ; Pi = q_in - (1-b) Phi
    Phi = Phi.at[-1].set((qi - qo) / 2)
    Pi = Pi.at[-1].set(qi - (1 - b) * (qi - qo) / 2)
    # left edge: rightward movers (q_out) enter the domain -> no inflow,
    # keep the leftward q_in from the interior
    qi_l = Pi[0] + (1 - b) * Phi[0]
    Phi = Phi.at[0].set((qi_l - 0.0) / 2)
    Pi = Pi.at[0].set(qi_l - (1 - b) * (qi_l - 0.0) / 2)
    return (u, Pi, Phi)


@partial(jax.jit, static_argnames=("n_steps", "n_grid", "order", "ko_sigma"))
def evolve(b, kind, n_grid: int, n_steps: int, dt, h, x, t0, pulse_w,
           datum_amp, order: int = 2, ko_sigma: float = 0.0):
    """Evolve an outgoing Gaussian pulse (and, for kind==2, a CCM-injected
    incoming pulse). Returns (u, Pi, Phi, t_final, R_metrics).

    Outgoing pulse: u0 = exp(-((x - x0) / w)^2) moving with speed c+ = 1+b.
    Incoming/outgoing characteristic content at any time:
       q_out = Pi + (1 - b) Phi   (right-mover content, speed 1+b)
       q_in  = Pi - (1 + b) Phi   (left-mover content, speed 1-b)
    (left eigenvectors of the flux matrix; verified by the test's pre-hit
    diagnostics: q_in == 0 for the pure outgoing pulse.)
    R is measured as max|q_in| after reflection / max|q_out| before the hit.
    For kind==2 the datum injects G(t + x/(1-b)), G a Gaussian of amplitude
    datum_amp, and the test compares u to the exact incoming solution.
    """
    x0 = -0.5 * (-x[0])
    amp0 = jnp.where(kind == 2, 0.0, 1.0)   # injection runs start from vacuum
    u0 = amp0 * jnp.exp(-(((x - x0) / pulse_w) ** 2))
    # purely outgoing initial data: u(t, x) = f(x - (1-b) t) =>
    # Pi = -(1-b) f', Phi = f'  (so q_in = Pi + (1-b) Phi = 0)
    f1 = -2 * (x - x0) / pulse_w**2 * u0
    state = (u0, -(1 - b) * f1, f1)

    q_out0 = jnp.max(jnp.abs(state[1] - (1 + b) * state[2]))

    def datum(t):
        # exact incoming solution G(t + x/(1+b)) with G(s) = A exp(-((s-sc)/w)^2)
        # datum g(t) = q_in = 2 G'(t)/(1+b) at x = 0
        sc = t0
        G1 = (-2 * (t - sc) / pulse_w**2
              * datum_amp * jnp.exp(-(((t - sc) / pulse_w) ** 2)))
        return 2.0 * G1 / (1 + b)

    def step(carry, i):
        st, t = carry
        def f(s):
            return _rhs(s, b, h, order, ko_sigma)
        k1 = f(st)
        s2 = tuple(a + 0.5 * dt * k for a, k in zip(st, k1))
        s2 = _apply_bc(s2, b, kind, datum(t + 0.5 * dt))
        k2 = f(s2)
        s3 = tuple(a + 0.5 * dt * k for a, k in zip(st, k2))
        s3 = _apply_bc(s3, b, kind, datum(t + 0.5 * dt))
        k3 = f(s3)
        s4 = tuple(a + dt * k for a, k in zip(st, k3))
        s4 = _apply_bc(s4, b, kind, datum(t + dt))
        k4 = f(s4)
        new = tuple(a + dt / 6 * (p + 2 * q + 2 * r + s_)
                    for a, p, q, r, s_ in zip(st, k1, k2, k3, k4))
        new = _apply_bc(new, b, kind, datum(t + dt))
        return (new, t + dt), 0.0

    (state, tf), _ = jax.lax.scan(step, (state, 0.0), jnp.arange(n_steps))
    q_in = state[1] + (1 - b) * state[2]
    q_out = state[1] - (1 + b) * state[2]
    return state[0], state[1], state[2], tf, q_out0, jnp.max(jnp.abs(q_in)), q_out


def char_decompose(Pi, Phi, b):
    """Boundary-mode decomposition of the model system (the 1D analog of the
    10-mode table's bookkeeping): returns (q_out, q_in).

    For outgoing f(x - (1-b)t): Pi = -(1-b) f', Phi = f' => q_out = -2 f',
    q_in = 0. For incoming g(x + (1+b)t): Pi = (1+b) g', Phi = g' =>
    q_in = 2 g', q_out = 0. Hence R_u = max|q_in| / max|q_out| directly.
    """
    return Pi - (1 + b) * Phi, Pi + (1 - b) * Phi
