#!/usr/bin/env python3
"""4-GPU test for packages/zccm z4c_vars (formulation-DAG nodes N1, N2, N3).

Per-device batches of random Z4c states (constraints det gamt = 1,
tr_gt At = 0 enforced by construction), float64:
  1. round trip from_adm(to_adm(state)) == state      (rtol 1e-12)
  2. four_metric * four_metric_inverse == identity    (atol 1e-11)
  3. boost_factor positivity + N2 closed form
  4. psi0_to_physical_datum: TT (symmetric, traceless) and K = 4 scaling

Budget: <= 10 min, <= 4 GPUs. Output -> results/numerical/ (run via bash with
timeout 600 | tee).
"""
import os, signal, sys, time

signal.alarm(540)
T0 = time.time()
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1,2,3")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import jax, jax.numpy as jnp
from zccm import (Z4cState, to_adm, from_adm, four_metric,
                  four_metric_inverse, boost_factor, psi0_to_physical_datum)

OK = True
def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  (t={time.time()-T0:.1f}s)")
    OK &= bool(cond)

devs = jax.devices()[:4]
print(f"JAX devices used ({len(devs)}): {[d.platform + ':' + str(d.id) for d in devs]}")

def random_state(key, B):
    ks = jax.random.split(key, 6)
    M = jax.random.uniform(ks[0], (B, 3, 3), jnp.float64, -0.3, 0.3)
    G = jnp.eye(3) + 0.5 * (M + jnp.swapaxes(M, 1, 2))
    gamt = G / jnp.linalg.det(G)[:, None, None] ** (1.0 / 3.0)   # det = 1
    A = jax.random.uniform(ks[1], (B, 3, 3), jnp.float64, -0.2, 0.2)
    A = 0.5 * (A + jnp.swapaxes(A, 1, 2))
    gti = jnp.linalg.inv(gamt)
    A = A - gamt * (jnp.einsum("bij,bij->b", gti, A) / 3.0)[:, None, None]
    return Z4cState(
        chi=jax.random.uniform(ks[2], (B,), jnp.float64, 0.5, 2.0),
        gamt=gamt, At=A,
        Khat=jax.random.uniform(ks[3], (B,), jnp.float64, -1.0, 1.0),
        Theta=jax.random.uniform(ks[4], (B,), jnp.float64, -0.1, 0.1),
        alpha=jax.random.uniform(ks[5], (B,), jnp.float64, 0.5, 2.0),
        beta=jax.random.uniform(ks[0], (B, 3), jnp.float64, -0.3, 0.3))

def device_checks(key):
    B = 1 << 18
    s = random_state(key, B)
    gamma, K_ij, alpha, beta = to_adm(s)
    s2 = from_adm(gamma, K_ij, alpha, beta, Theta=s.Theta)
    rt = jnp.max(jnp.stack([jnp.max(jnp.abs(a - b) / (1.0 + jnp.abs(a)))
                            for a, b in zip(s, s2)]))
    g4 = four_metric(s)
    g4i = four_metric_inverse(s)
    eye = jnp.broadcast_to(jnp.eye(4), g4.shape)
    inv = jnp.max(jnp.abs(jnp.einsum("bij,bjk->bik", g4, g4i) - eye))
    # boost: unit normal along gamma-normalized x-axis
    e1 = jnp.zeros((B, 3)).at[:, 0].set(1.0)
    nrm = jnp.sqrt(jnp.einsum("bij,bi,bj->b", gamma, e1, e1))
    s_unit = e1 / nrm[:, None]
    bchar = jax.random.uniform(key, (B,), jnp.float64, -0.4, 0.4)
    A2 = boost_factor(s, s_unit, bchar)
    bs = jnp.einsum("bij,bi,bj->b", gamma, beta, s_unit)
    bres = jnp.max(jnp.abs(A2 - (alpha - bs) * jnp.exp(-2 * bchar)))
    bpos = jnp.all(A2 > 0)
    # physical datum
    psi0 = (jax.random.uniform(key, (B,), jnp.float64, -1, 1)
            + 1j * jax.random.uniform(key, (B,), jnp.float64, -1, 1))
    w = psi0_to_physical_datum(psi0)
    tt = jnp.max(jnp.abs(w[:, 0, 0] + w[:, 1, 1]))           # traceless
    sym = jnp.max(jnp.abs(w[:, 0, 1] - w[:, 1, 0]))          # symmetric
    k4 = jnp.max(jnp.abs(w[:, 0, 0] - 4 * jnp.real(psi0)))   # K = 4
    return rt, inv, bres, bpos, tt, sym, k4

res = [jax.jit(device_checks, device=d)(jax.random.PRNGKey(31 + i))
       for i, d in enumerate(devs)]
B_tot = len(devs) * (1 << 18)
rt = max(float(r[0]) for r in res)
inv = max(float(r[1]) for r in res)
bres = max(float(r[2]) for r in res)
bpos = all(bool(r[3]) for r in res)
tt = max(float(r[4]) for r in res)
sym = max(float(r[5]) for r in res)
k4 = max(float(r[6]) for r in res)

report(f"round trip from_adm(to_adm(s)) == s: max rel residual {rt:.2e} < 1e-12 "
       f"over {B_tot:,} states", rt < 1e-12)
report(f"g4 * g4inv == identity: max abs residual {inv:.2e} < 1e-11", inv < 1e-11)
report(f"boost factor matches closed form (residual {bres:.2e}) and is positive "
       f"({bpos}) for subluminal shift", bres < 1e-14 and bpos)
report(f"physical datum TT (trace {tt:.1e}, asym {sym:.1e}) with K = 4 "
       f"(residual {k4:.1e})", tt < 1e-14 and sym < 1e-14 and k4 < 1e-14)

print(f"\nwall clock: {time.time()-T0:.1f}s")
print(f"OVERALL: {'PASS' if OK else 'FAIL'}")
raise SystemExit(0 if OK else 1)
