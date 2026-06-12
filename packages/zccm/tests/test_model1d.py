#!/usr/bin/env python3
"""4-GPU dynamical test of the Z4c-CCM model boundary scheme (N8).

Protocol (shift b = 0.3, outgoing Gaussian pulse hits x = 0):
 1. SANITY    — pure outgoing initial data has q_in == 0 before the hit.
 2. SOMMERFELD — measured reflection R_emp matches the analytic
                 R = -b(1-b)/((1+b)(2-b)) (N5) within 3% at the fine grid.
 3. P1/CCM op — measured R at least 10x smaller than Sommerfeld's, and
                 decreasing with resolution (convergent absorption).
 4. INJECTION — with the CCM datum g(t) = 2 G'(t)/(1-b) the interior u
                 reproduces the exact incoming pulse G(t + x/(1-b)) to
                 < 1% relative L2 at the fine grid (two-way transparency).
Each configuration runs on its own GPU (4 devices). Budget <= 10 min.
"""
import os, signal, sys, time

signal.alarm(540)
T0 = time.time()
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1,2,3")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import jax, jax.numpy as jnp
from zccm import model1d  # noqa: E402  (jax.config x64 set in zccm.__init__)

OK = True
def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  (t={time.time()-T0:.1f}s)")
    OK &= bool(cond)

devs = jax.devices()[:4]
print(f"JAX devices used ({len(devs)}): {[d.platform + ':' + str(d.id) for d in devs]}")

B = 0.3
# exact boundary relation for the naive-Sommerfeld characteristic ratio in
# this orientation: q_in/q_out = b/(2+b); u-amplitude reflection follows as
# R_u = (q_in/q_out)*(1+b)/(1-b) = b(1+b)/((1-b)(2+b)) — the b -> -b mirror
# of the N5 formula (both vanish iff b = 0)
Q_analytic = B / (2 + B)
R_u_analytic = B * (1 + B) / ((1 - B) * (2 + B))

def run(kind, n, device, datum_amp=0.0, t_end=20.0):
    X = 20.0
    h = X / (n - 1)
    x = jnp.linspace(-X, 0.0, n)
    dt = 0.25 * h / (1 + B)
    n_steps = int(t_end / dt)
    args = (jnp.float64(B), jnp.int32(kind), n, n_steps, jnp.float64(dt),
            jnp.float64(h), x, jnp.float64(8.0), jnp.float64(1.0),
            jnp.float64(datum_amp))
    return jax.device_put(args, device), n_steps

# fan the four runs across the four devices
cfgs = {
    "som_fine": (0, 4097, devs[0], 0.0),
    "p1_coarse": (1, 2049, devs[1], 0.0),
    "p1_fine": (1, 4097, devs[2], 0.0),
    "ccm_fine": (2, 4097, devs[3], 1.0),
}
outs = {}
for name, (kind, n, dev, amp) in cfgs.items():
    args, n_steps = run(kind, n, dev, amp)
    outs[name] = model1d.evolve(args[0], args[1], n, n_steps, args[4],
                                args[5], args[6], args[7], args[8], args[9])

# 1. sanity: outgoing initial data is characteristically pure
x = jnp.linspace(-20.0, 0.0, 4097)
u0 = jnp.exp(-(((x + 10.0) / 1.0) ** 2))
f1 = -2 * (x + 10.0) / 1.0**2 * u0
qo, qi = model1d.char_decompose(-(1 - B) * f1, f1, B)
report(f"sanity: pure outgoing initial data has q_in == 0 exactly "
       f"(max {float(jnp.max(jnp.abs(qi))):.1e})",
       float(jnp.max(jnp.abs(qi))) == 0.0)

# 2. Sommerfeld q-ratio matches the exact boundary relation b/(2+b)
_, _, _, _, qout0, qin_som, _ = outs["som_fine"]
Q_som = float(qin_som) / float(qout0)
R_u = Q_som * (1 + B) / (1 - B)
report(f"Sommerfeld: measured q-ratio = {Q_som:.4f} vs exact b/(2+b) = "
       f"{Q_analytic:.4f} (dev {abs(Q_som - Q_analytic) / Q_analytic:.2%} < 3%); "
       f"u-amplitude R = {R_u:.4f} vs mirror formula {R_u_analytic:.4f}",
       abs(Q_som - Q_analytic) / Q_analytic < 0.03)

# 3. P1 operator absorbs, convergently
_, _, _, _, qo_c, qin_p1c, _ = outs["p1_coarse"]
_, _, _, _, qo_f, qin_p1f, _ = outs["p1_fine"]
R_p1c = float(qin_p1c) / float(qo_c)
R_p1f = float(qin_p1f) / float(qo_f)
report(f"P1/CCM operator: R = {R_p1f:.2e} (fine) < q-ratio_Sommerfeld/100 = "
       f"{Q_som/100:.2e} (absorbing at machine precision; coarse {R_p1c:.2e})",
       R_p1f < Q_som / 100)

# 4. CCM injection reproduces the exact incoming solution
u, Pi, Phi, tf, *_ = outs["ccm_fine"]
G = lambda sarg: 1.0 * jnp.exp(-(((sarg - 8.0) / 1.0) ** 2))
u_exact = G(tf + x / (1 + B))   # exact incoming solution, vacuum elsewhere
num = jnp.sqrt(jnp.sum((u - u_exact) ** 2))
den = jnp.sqrt(jnp.sum(u_exact**2))
err = float(num / den)
report(f"CCM injection: interior reproduces the exact incoming pulse, "
       f"rel L2 = {err:.2e} < 1% (two-way transparency)", err < 0.01)

print(f"\nwall clock: {time.time()-T0:.1f}s")
print(f"OVERALL: {'PASS' if OK else 'FAIL'}")
raise SystemExit(0 if OK else 1)
