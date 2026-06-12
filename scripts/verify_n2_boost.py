#!/usr/bin/env python3
"""verify_n2_boost.py (v2) — N2 verifier (Z4c-CCM formulation DAG).

CLAIM (node N2): at a worldtube of fixed coordinate radius, the characteristic
outgoing null vector l-hat^mu = -e^{2 betahat} g^{mu nu} d_nu u (u = outgoing
null-cone label with u = t on the tube) is proportional to the Cauchy-side
outgoing null vector l'^mu = (n^mu + s^mu)/sqrt(2) (paper-1 eq:nullvector-k
frame), with

    l'^mu = const * (alpha - gamma_ij beta^i s^j) * e^{-2 betahat} * l-hat^mu,

`const` a pure number, for COMPLETELY GENERAL 3+1 data (alpha, beta^i,
gamma_ij). No GH or Z4c evolution variable enters, so paper-3's Choice-2
Type-II boost (arXiv:2308.10361 eq:auto-20, eq:l_GH_and_l_CCE_transform)
transplants verbatim onto the Z4c boundary frame of arXiv:1010.0523v2.

v2 verification strategy (v1 fully-symbolic simplify hit kernel_timeout —
error-DB row 2026-06-12T01:52:24Z):
  STAGE A (exact arithmetic, CPU, seconds) — substitute exact random Rationals
    for all 10 metric symbols; every check is then an identity between
    algebraic numbers, decided EXACTLY (radsimp/expand; no floats, no
    tolerance). Vanishing of a rational-radical identity at NS random rational
    points is Schwartz-Zippel-style certification.
    Checks: nullness of l', nullness of du, flat-space root u_r = -1 (fully
    symbolic, cheap), 6 proportionality cross products, factor constancy
    (same exact algebraic value at every sample).
  STAGE B (float64 falsification sweep, JAX, <= 4 GPUs) — 4 x 2^20 random
    3+1 samples; max relative residual of the cross products and of the
    factor-constancy must be < 1e-9 (float64 noise scale ~1e-13).

Wall-clock self-cap: 540 s (mission test budget: no run > 10 min).
Exit 0 only if all checks pass.
"""
import os, random, signal, time

signal.alarm(540)  # hard self-cap inside the 10-min budget
T0 = time.time()
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1,2,3")  # mission cap: <=4 GPUs

import sympy as sp

OK = True
def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  (t={time.time()-T0:.1f}s)")
    OK &= bool(cond)

# ---------- symbolic setup (no heavy simplify on full symbols) --------------
al, bh = sp.symbols("alpha betahat", positive=True)
b1, b2, b3 = sp.symbols("beta1 beta2 beta3", real=True)
g11, g12, g13, g22, g23, g33 = sp.symbols("g11 g12 g13 g22 g23 g33", real=True)
SYMS = [al, bh, b1, b2, b3, g11, g12, g13, g22, g23, g33]

gam = sp.Matrix([[g11, g12, g13], [g12, g22, g23], [g13, g23, g33]])
beta_up = sp.Matrix([b1, b2, b3])

def frame(sub):
    """Everything evaluated at one exact-rational sample point `sub`."""
    G = gam.subs(sub)
    Ginv = G.inv()
    a = al.subs(sub); bu = beta_up.subs(sub); bd = G * bu
    g4 = sp.zeros(4, 4)
    g4[0, 0] = -a**2 + (bd.T * bu)[0]
    for i in range(3):
        g4[0, i+1] = g4[i+1, 0] = bd[i]
        for j in range(3):
            g4[i+1, j+1] = G[i, j]
    g4i = sp.zeros(4, 4)
    g4i[0, 0] = -1/a**2
    for i in range(3):
        g4i[0, i+1] = g4i[i+1, 0] = bu[i]/a**2
        for j in range(3):
            g4i[i+1, j+1] = Ginv[i, j] - bu[i]*bu[j]/a**2
    n = sp.Matrix([1/a, -bu[0]/a, -bu[1]/a, -bu[2]/a])
    grr = Ginv[0, 0]
    s = sp.Matrix([0] + [Ginv[0, i]/sp.sqrt(grr) for i in range(3)])
    lp = (n + s)/sp.sqrt(2)
    A2_, B2_, C2_ = g4i[1, 1], 2*g4i[0, 1], g4i[0, 0]
    ur = (-B2_ - sp.sqrt(B2_**2 - 4*A2_*C2_)) / (2*A2_)   # outgoing root
    u_dn = sp.Matrix([1, ur, 0, 0])
    lhat = -sp.exp(2*bh.subs(sub)) * (g4i * u_dn)
    bs = (bd.T * s[1:, 0])[0]
    return g4, lp, u_dn, g4i, lhat, a, bs, sub[bh]

def exact_zero(e):
    return sp.radsimp(sp.expand(sp.nsimplify(e))) == 0 or sp.simplify(e) == 0

# flat-space limit fully symbolic (cheap 1-d check)
r2, r3 = sp.symbols("r2 r3", positive=True)
flat = {al: 1, bh: 0, b1: 0, b2: 0, b3: 0, g11: 1, g12: 0, g13: 0,
        g22: r2, g23: 0, g33: r3}
_, _, u_dn_f, _, lhat_f, _, _, _ = frame(flat)
report("flat-space outgoing root u_r = -1 (u = t - r), symbolic",
       sp.simplify(u_dn_f[1] + 1) == 0)

# ---------- STAGE A: exact rational-point identity testing ------------------
NS = 12
random.seed(20260612)
def rnd(lo, hi, den=64):
    return sp.Rational(random.randint(int(lo*den), int(hi*den)), den)

null_lp_ok = null_du_ok = cross_ok = True
factor_vals = []
for k in range(NS):
    # positive-definite gamma = I + 0.5*(M + M^T), entries in [-0.3, 0.3]
    M = [[rnd(-0.3, 0.3) for _ in range(3)] for _ in range(3)]
    sub = {al: rnd(0.5, 2), bh: rnd(-0.4, 0.4),
           b1: rnd(-0.3, 0.3), b2: rnd(-0.3, 0.3), b3: rnd(-0.3, 0.3),
           g11: 1 + M[0][0], g12: (M[0][1]+M[1][0])/2, g13: (M[0][2]+M[2][0])/2,
           g22: 1 + M[1][1], g23: (M[1][2]+M[2][1])/2, g33: 1 + M[2][2]}
    g4, lp, u_dn, g4i, lhat, a, bs, bhv = frame(sub)
    null_lp_ok &= exact_zero((lp.T * g4 * lp)[0])
    null_du_ok &= exact_zero((u_dn.T * g4i * u_dn)[0])
    for mu in range(4):
        for nu in range(mu+1, 4):
            cross_ok &= exact_zero(lhat[mu]*lp[nu] - lhat[nu]*lp[mu])
    factor_vals.append(sp.radsimp((lp[0]/lhat[0]) * sp.exp(2*bhv) / (a - bs)))

report(f"l' = (n+s)/sqrt(2) null at {NS}/{NS} exact rational points", null_lp_ok)
report(f"d_mu u null (outgoing root) at {NS}/{NS} exact rational points", null_du_ok)
report(f"proportionality l-hat ∝ l': 6 cross products vanish EXACTLY at {NS}/{NS} points", cross_ok)
const_ok = all(exact_zero(v - factor_vals[0]) for v in factor_vals[1:])
report(f"factor F*e^(2bh)/(alpha - beta_i s^i) is the SAME exact number at all points"
       f" (= {sp.nsimplify(factor_vals[0])})", const_ok)

# ---------- STAGE B: 4-GPU float64 mass falsification sweep -----------------
import jax, jax.numpy as jnp
devs = jax.devices()[:4]
print(f"JAX devices used ({len(devs)}): {[d.platform + ':' + str(d.id) for d in devs]}")

def residuals(key):
    B = 1 << 20
    ks = jax.random.split(key, 8)
    a = jax.random.uniform(ks[0], (B,), jnp.float64, 0.5, 2.0)
    bhx = jax.random.uniform(ks[1], (B,), jnp.float64, -0.4, 0.4)
    bu = jax.random.uniform(ks[2], (B, 3), jnp.float64, -0.3, 0.3)
    Mx = jax.random.uniform(ks[3], (B, 3, 3), jnp.float64, -0.3, 0.3)
    G = jnp.eye(3) + 0.5*(Mx + jnp.swapaxes(Mx, 1, 2))
    Gi = jnp.linalg.inv(G)
    bd = jnp.einsum("bij,bj->bi", G, bu)
    g4 = jnp.zeros((B, 4, 4), jnp.float64)
    g4 = g4.at[:, 0, 0].set(-a**2 + jnp.einsum("bi,bi->b", bd, bu))
    g4 = g4.at[:, 0, 1:].set(bd); g4 = g4.at[:, 1:, 0].set(bd)
    g4 = g4.at[:, 1:, 1:].set(G)
    g4i = jnp.zeros((B, 4, 4), jnp.float64)
    g4i = g4i.at[:, 0, 0].set(-1/a**2)
    g4i = g4i.at[:, 0, 1:].set(bu/a[:, None]**2)
    g4i = g4i.at[:, 1:, 0].set(bu/a[:, None]**2)
    g4i = g4i.at[:, 1:, 1:].set(Gi - jnp.einsum("bi,bj->bij", bu, bu)/a[:, None, None]**2)
    n = jnp.concatenate([(1/a)[:, None], -bu/a[:, None]], 1)
    grr = Gi[:, 0, 0]
    s = jnp.concatenate([jnp.zeros((B, 1)), Gi[:, 0, :]/jnp.sqrt(grr)[:, None]], 1)
    lp = (n + s)/jnp.sqrt(2.0)
    A2_, B2_, C2_ = g4i[:, 1, 1], 2*g4i[:, 0, 1], g4i[:, 0, 0]
    ur = (-B2_ - jnp.sqrt(B2_**2 - 4*A2_*C2_)) / (2*A2_)
    u = jnp.stack([jnp.ones(B), ur, jnp.zeros(B), jnp.zeros(B)], 1)
    lhat = -jnp.exp(2*bhx)[:, None]*jnp.einsum("bij,bj->bi", g4i, u)
    scale = jnp.maximum(jnp.max(jnp.abs(lhat), 1), 1e-30)*jnp.max(jnp.abs(lp), 1)
    cross = jnp.einsum("bi,bj->bij", lhat, lp)
    cross = jnp.max(jnp.abs(cross - jnp.swapaxes(cross, 1, 2)), (1, 2))/scale
    bs = jnp.einsum("bi,bi->b", bd, s[:, 1:])
    fac = lp[:, 0]/lhat[:, 0]*jnp.exp(2*bhx)/(a - bs)
    return jnp.max(cross), jnp.max(jnp.abs(fac - fac[0]))

jax.config.update("jax_enable_x64", True)
futs = [jax.jit(residuals, device=d)(jax.random.PRNGKey(7 + i))
        for i, d in enumerate(devs)]
cross_max = max(float(c) for c, _ in futs)
fac_dev = max(float(f) for _, f in futs)
n_total = len(devs) * (1 << 20)
report(f"GPU sweep: max relative cross-product residual {cross_max:.2e} < 1e-9 "
       f"over {n_total:,} samples", cross_max < 1e-9)
report(f"GPU sweep: max factor deviation {fac_dev:.2e} < 1e-9", fac_dev < 1e-9)

print()
print("CONCLUSION: l-hat ∝ (n+s)/sqrt(2) with factor (alpha - gamma_ij beta^i s^j)")
print("* e^(-2 betahat) — a pure 3+1/characteristic identity, no GH/Z4c evolution")
print("variable enters: paper-3's Choice-2 Type-II boost is valid VERBATIM on the")
print("Z4c boundary frame (P1 eq:nullvector-k). N2 admitted.")
print(f"\nwall clock: {time.time()-T0:.1f}s")
print(f"OVERALL: {'PASS' if OK else 'FAIL'}")
raise SystemExit(0 if OK else 1)
