#!/usr/bin/env python3
"""verify_n3_exact.py — N3-EXACT verifier (Z4c-CCM formulation DAG).

Replaces the linearized dictionary (old eqs. 16-17) with the EXACT one.

CLAIMS:
 A. EXACT WEYL ALGEBRA (no perturbation, no field equations): for a GENERIC
    Weyl-symmetric tensor C_{abcd} (Riemann symmetries + cyclic + all traces
    zero; 10 independent components) in an orthonormal frame
    (e0 = n, e1 = s, e2, e3), with
       E_ij = C_{0i0j},
       B_ij = (1/2) eps_{ikl} C_{klj0},
       m^a  = (e2 + i e3)/sqrt(2),  l = (n+s)/sqrt(2),
       psi0 = -C(l, m, l, m),
    there exist a UNIQUE sign sigma and pure number kappa such that
       (E + sigma * eps(s)B)(m, m) = kappa * psi0,        [exact, all 10 dof]
    where (eps(s)B)_ij = eps_{ikl} s_k B_{lj} symmetrized. The opposite sign
    pairs with the k-leg scalar psi4-bar instead (verified). Since a
    tangentially-TT symmetric tensor U satisfies U_AB = U(mbar,mbar) m_A m_B
    + U(m,m) mbar_A mbar_B, this scalar identity is equivalent to
       U^-_AB := TT2[E + sigma eps(s)B]_AB
              = kappa (psi0 mbar_A mbar_B + conj(psi0) m_A m_B).
 B. LINEARIZED LIMIT (consistency with the superseded eqs. 16-17): for a TT
    perturbation h_AB(t,x) with E = R_ij[gamma] (linear), K = -h_dot/2,
    B = eps D K, ON SHELL (box h = 0), the same (sigma, kappa) give
       (E + sigma eps(s)B)(m,m) = (kappa/4)(d_t + d_s)^2 h_mm = kappa psi0_lin.
 C. NONPERTURBATIVE NUMERIC (4 GPUs, JAX autodiff, float64): Kerr-Schild
    Schwarzschild (fully nonlinear, constraints exact). Two routes:
      4D route: Riemann from autodiff -> E, B, psi0 directly;
      3+1 route: gamma_ij, K_ij = D_(i beta_j)/alpha (stationary) ->
                 E = R_ij[gamma] + K K_ij - K_ik K^k_j,
                 B = sigma_B eps_i^{kl} D_k K_{lj}  (Gauss-Codazzi, vacuum).
    Residuals: |E_3+1 - E_4D|, |B_3+1 - B_4D|, identity residual
    |(E + sigma eps(s)B)(m,m) - kappa psi0| < 1e-8 at random points
    r in [3M, 12M], random frames.

Budget: alarm 540 s, <= 4 GPUs. Exit 0 iff all checks pass.
"""
import os, signal, time
from itertools import product

signal.alarm(540)
T0 = time.time()
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1,2,3")

import sympy as sp

OK = True
def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  (t={time.time()-T0:.1f}s)", flush=True)
    OK &= bool(cond)

# ====================== A. exact Weyl algebra ================================
eta = sp.diag(-1, 1, 1, 1)
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
M = sp.Matrix(6, 6, lambda I, J: sp.Symbol(f"M{min(I,J)}{max(I,J)}"))

def Cc(a, b, c, d):
    """Riemann-symmetric tensor from the 6x6 symmetric pair matrix."""
    sgn = 1
    if a == b or c == d:
        return sp.Integer(0)
    if a > b:
        a, b, sgn = b, a, -sgn
    if c > d:
        c, d, sgn = d, c, -sgn
    return sgn * M[pairs.index((a, b)), pairs.index((c, d))]

# constraints: cyclic identity + vanishing traces (Weyl)
cons = [Cc(0, 1, 2, 3) + Cc(0, 2, 3, 1) + Cc(0, 3, 1, 2)]
for b in range(4):
    for d in range(b, 4):
        cons.append(sp.expand(sum(eta[a, a] * Cc(a, b, a, d) for a in range(4))))
syms = sorted(M.free_symbols, key=str)
sol = sp.solve([c for c in cons if c != 0], syms, dict=True)[0]
subC = {s: sol.get(s, s) for s in syms}
free = sorted({v for s, v in subC.items() for v in sp.sympify(v).free_symbols},
              key=str)
report(f"generic Weyl tensor constructed: {len(free)} independent components "
       "(cyclic + 10 trace conditions imposed)", len(free) == 10)

def C(a, b, c, d):
    return sp.expand(Cc(a, b, c, d).subs(subC))

eps3 = lambda i, j, k: sp.LeviCivita(i, j, k)
E = sp.Matrix(3, 3, lambda i, j: C(0, i + 1, 0, j + 1))
B = sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * sp.expand(
    sum(eps3(i + 1, k, l) * C(k, l, j + 1, 0)
        for k in range(1, 4) for l in range(1, 4))))

report("E_ij symmetric & traceless; B_ij symmetric & traceless (Weyl algebra)",
       all(sp.expand(E[i, j] - E[j, i]) == 0 for i in range(3) for j in range(3))
       and sp.expand(E.trace()) == 0
       and all(sp.expand(B[i, j] - B[j, i]) == 0 for i in range(3) for j in range(3))
       and sp.expand(B.trace()) == 0)

# frame vectors (orthonormal): n=e0, s=e1, m=(e2+ie3)/sqrt2; l,k null
lv = sp.Matrix([1, 1, 0, 0]) / sp.sqrt(2)
kv = sp.Matrix([1, -1, 0, 0]) / sp.sqrt(2)
mv = sp.Matrix([0, 0, 1, sp.I]) / sp.sqrt(2)

def con4(u, v, w, q):
    return sp.expand(sum(C(a, b, c, d) * u[a] * v[b] * w[c] * q[d]
                         for a in range(4) for b in range(4)
                         for c in range(4) for d in range(4)))

psi0 = sp.expand(-con4(lv, mv, lv, mv))
psi4 = sp.expand(-con4(kv, mv.conjugate(), kv, mv.conjugate()))

# (eps(s)B)_ij with s = e1: eps_{i1l} B_{lj}, symmetrized
sB = sp.Matrix(3, 3, lambda i, j: sum(eps3(i + 1, 1, l) * B[l - 1, j]
                                      for l in range(1, 4)))
sB = (sB + sB.T) / 2
m3 = sp.Matrix([0, 1, sp.I]) / sp.sqrt(2)     # spatial part of m

def mm(T):
    return sp.expand(sum(T[i, j] * m3[i] * m3[j] for i in range(3) for j in range(3)))

found = None
for sigma, kappa in product((1, -1), (1, -1, 2, -2, sp.Rational(1, 2), sp.Rational(-1, 2))):
    if sp.expand(mm(E + sigma * sB) - kappa * psi0) == 0:
        found = (sigma, kappa)
        break
report(f"EXACT identity (E + sigma*eps(s)B)(m,m) = kappa*psi0 with "
       f"(sigma, kappa) = {found}, generic Weyl (all 10 dof)", found is not None)
SIG, KAP = found if found else (1, 1)

psi4bar = sp.expand(-con4(kv, mv, kv, mv))   # = conj(psi4) for real components
found4 = None
for kappa in (1, -1, 2, -2):
    if sp.expand(mm(E - SIG * sB) - kappa * psi4bar) == 0:
        found4 = kappa
        break
report(f"opposite sign pairs with the outgoing scalar: "
       f"(E - sigma*eps(s)B)(m,m) = {found4}*conj(psi4)", found4 is not None)

# completeness of the dyad expansion for tangential-TT tensors
u1, u2 = sp.symbols("u1 u2", real=True)
U = sp.Matrix([[u1, u2], [u2, -u1]])
mA = sp.Matrix([1, sp.I]) / sp.sqrt(2)
Umm = sp.expand(sum(U[i, j] * mA[i] * mA[j] for i in range(2) for j in range(2)))
rec = sp.zeros(2, 2)
for A in range(2):
    for Bb in range(2):
        rec[A, Bb] = (Umm * mA.conjugate()[A] * mA.conjugate()[Bb]
                      + sp.conjugate(Umm) * mA[A] * mA[Bb])
report("dyad completeness: U_AB = U(m,m) mbar_A mbar_B + cc for any TT U "
       "(scalar identity <=> tensor identity)",
       all(sp.simplify(sp.expand(rec[i, j] - U[i, j])) == 0
           for i in range(2) for j in range(2)))

# ====================== B. linearized on-shell limit =========================
t, x = sp.symbols("t x", real=True)
F1 = sp.Function("F1", real=True)(t, x)
F2 = sp.Function("F2", real=True)(t, x)
h = sp.Matrix([[0, 0, 0], [0, F1, F2], [0, F2, -F1]])
# linear pieces: R_ij = -1/2 d_x^2 h_ij (TT, depends on x only),
# K_ij = -1/2 d_t h_ij, B_ij = SIGB eps_{ikl} d_k K_{lj}; here d_k = delta_kx d_x
Elin = -sp.diff(h, x, 2) / 2
K = -sp.diff(h, t) / 2
Blin = sp.Matrix(3, 3, lambda i, j: sum(eps3(i + 1, k + 1, 1) * 0 for k in range(3)))
Blin = sp.Matrix(3, 3, lambda i, j: sp.expand(
    sum(eps3(i + 1, 1, l) * sp.diff(K[l - 1, j], x) * (-1) for l in range(1, 4))))
# NOTE: B_ij = eps_i^{kl} D_k K_{lj} with D_k -> d_x delta_k1: eps_{i1l}(-d_x K_{lj})
# sign convention folded into SIGB below (determined, not assumed)
Blin_s = (Blin + Blin.T) / 2
hmm_ = F1 + sp.I * F2
target = sp.Rational(1, 4) * KAP * (sp.diff(hmm_, t) + sp.diff(hmm_, x)).diff(t) \
    + sp.Rational(1, 4) * KAP * (sp.diff(hmm_, t) + sp.diff(hmm_, x)).diff(x)
# target = (kappa/4)(d_t + d_x)^2 h_mm
onshell = {sp.Derivative(F1, (t, 2)): sp.Derivative(F1, (x, 2)),
           sp.Derivative(F2, (t, 2)): sp.Derivative(F2, (x, 2))}
SIGB_found = None
for SIGB in (1, -1):
    sBlin = sp.Matrix(3, 3, lambda i, j: sum(eps3(i + 1, 1, l) * SIGB * Blin_s[l - 1, j]
                                             for l in range(1, 4)))
    sBlin = (sBlin + sBlin.T) / 2
    lhs = mm(Elin + SIG * sBlin)
    if sp.expand((lhs - target).subs(onshell)) == 0:
        SIGB_found = SIGB
        break
report(f"linearized ON-SHELL limit reproduces superseded eqs. 16-17: "
       f"(E + sigma eps(s)B)(m,m) = (kappa/4)(d_t+d_s)^2 h_mm with Codazzi "
       f"sign sigma_B = {SIGB_found}", SIGB_found is not None)

# ====================== C. nonperturbative Schwarzschild =====================
import numpy as np
import jax, jax.numpy as jnp
jax.config.update("jax_enable_x64", True)
devs = jax.devices()[:4]
print(f"JAX devices used ({len(devs)}): {[d.platform + ':' + str(d.id) for d in devs]}",
      flush=True)

Mbh = 1.0
SIGBf = float(SIGB_found)
SIGf = float(SIG)
KAPf = float(KAP)

# constant Levi-Civita arrays (densitized later by sqrt(det))
def levi(n):
    e = np.zeros((n,) * n)
    from itertools import permutations
    for p_ in permutations(range(n)):
        # permutation parity
        pl, vis, sgn = list(p_), [False] * n, 1
        for i0 in range(n):
            if vis[i0]:
                continue
            j0, cyc = i0, 0
            while not vis[j0]:
                vis[j0] = True
                j0 = pl[j0]
                cyc += 1
            if cyc % 2 == 0:
                sgn = -sgn
        e[p_] = sgn
    return jnp.asarray(e)

LC4, LC3 = levi(4), levi(3)

def gfun(X):
    r = jnp.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2)
    lmu = jnp.array([1.0, X[1] / r, X[2] / r, X[3] / r])
    return jnp.diag(jnp.array([-1.0, 1.0, 1.0, 1.0])) + 2 * (Mbh / r) * jnp.outer(lmu, lmu)

def christoffel(gfunc, X):
    gi = jnp.linalg.inv(gfunc(X))
    d = jax.jacfwd(gfunc)(X)                       # d[a,b,c] = d_c g_ab
    return 0.5 * (jnp.einsum("ad,dcb->abc", gi, d)
                  + jnp.einsum("ad,dbc->abc", gi, d)
                  - jnp.einsum("ad,bcd->abc", gi, d))

def riemann_up(gfunc, X):
    G = christoffel(gfunc, X)
    dG = jax.jacfwd(lambda Y: christoffel(gfunc, Y))(X)  # dG[a,b,c,e] = d_e G^a_bc
    return (jnp.transpose(dG, (0, 1, 3, 2)) - dG
            + jnp.einsum("ace,ebd->abcd", G, G)
            - jnp.einsum("ade,ebc->abcd", G, G))

def routes(X, v1, v2):
    g = gfun(X)
    gi = jnp.linalg.inv(g)
    Rl = jnp.einsum("ea,ebcd->abcd", g, riemann_up(gfun, X))
    # ADM split
    gam = g[1:, 1:]
    gami = jnp.linalg.inv(gam)
    beta_dn = g[0, 1:]
    beta_up = gami @ beta_dn
    alpha = jnp.sqrt(beta_dn @ beta_up - g[0, 0])
    n_up = jnp.concatenate([jnp.array([1.0]), -beta_up]) / alpha
    n_dn = g @ n_up
    P = jnp.eye(4) + jnp.outer(n_up, n_dn)         # spatial projector P^mu_nu
    # 4D route
    Enn = jnp.einsum("acbd,c,d->ab", Rl, n_up, n_up)
    E4 = jnp.einsum("ab,ai,bj->ij", Enn, P[:, 1:], P[:, 1:])
    sq4 = jnp.sqrt(-jnp.linalg.det(g))
    eps4 = sq4 * LC4
    epsud = jnp.einsum("abef,ec,fd->abcd", eps4, gi, gi)
    dualC = 0.5 * jnp.einsum("abef,efcd->abcd", epsud, Rl)
    Bnn = jnp.einsum("acbd,c,d->ab", dualC, n_up, n_up)
    B4 = jnp.einsum("ab,ai,bj->ij", Bnn, P[:, 1:], P[:, 1:])
    B4 = (B4 + B4.T) / 2
    # 3+1 route (stationary: K_ij = D_(i beta_j) / alpha)
    t0c = X[:1]
    gamf = lambda xs: gfun(jnp.concatenate([t0c, xs]))[1:, 1:]
    betf = lambda xs: gfun(jnp.concatenate([t0c, xs]))[0, 1:]
    xs0 = X[1:]
    G3 = christoffel(gamf, xs0)
    dbeta = jax.jacfwd(betf)(xs0)                  # dbeta[j,i] = d_i beta_j
    Db = dbeta.T - jnp.einsum("mij,m->ij", G3, beta_dn)
    def Kfun(xs):
        gp = gamf(xs)
        gip = jnp.linalg.inv(gp)
        G3p = christoffel(gamf, xs)
        dbp = jax.jacfwd(betf)(xs)
        bdp = betf(xs)
        Dbp = dbp.T - jnp.einsum("mij,m->ij", G3p, bdp)
        gfull = gfun(jnp.concatenate([t0c, xs]))
        bup = gip @ bdp
        al = jnp.sqrt(bdp @ bup - gfull[0, 0])
        return (Dbp + Dbp.T) / (2 * al)
    Kij = Kfun(xs0)
    R3up = riemann_up(gamf, xs0)
    Ric3 = jnp.einsum("abad->bd", R3up)
    trK = jnp.einsum("ij,ij->", gami, Kij)
    E3 = Ric3 + trK * Kij - Kij @ gami @ Kij
    dK = jax.jacfwd(Kfun)(xs0)                     # dK[l,j,k] = d_k K_lj
    DK = (jnp.transpose(dK, (2, 0, 1))
          - jnp.einsum("mkl,mj->klj", G3, Kij)
          - jnp.einsum("mkj,lm->klj", G3, Kij))
    sq3 = jnp.sqrt(jnp.linalg.det(gam))
    eps3lo = sq3 * LC3
    epsmix = jnp.einsum("ikl,km,ln->imn", eps3lo, gami, gami)  # eps_i^{mn}
    B3 = SIGBf * jnp.einsum("ikl,klj->ij", epsmix, DK)
    B3 = (B3 + B3.T) / 2
    # frame: s radial unit (gamma-norm), dyad from Gram-Schmidt on v1, v2
    rad = xs0 / jnp.linalg.norm(xs0)
    s_up = rad / jnp.sqrt(rad @ gam @ rad)
    w1 = v1 - (v1 @ gam @ s_up) * s_up
    e2 = w1 / jnp.sqrt(w1 @ gam @ w1)
    w2 = v2 - (v2 @ gam @ s_up) * s_up - (v2 @ gam @ e2) * e2
    e3v = w2 / jnp.sqrt(w2 @ gam @ w2)
    m3c = (e2 + 1j * e3v) / jnp.sqrt(2.0)
    s4 = jnp.concatenate([jnp.array([0.0]), s_up])
    l4 = ((n_up + s4) / jnp.sqrt(2.0)).astype(jnp.complex128)
    m4 = jnp.concatenate([jnp.array([0.0 + 0.0j]), m3c])
    psi0_4d = -jnp.einsum("abcd,a,b,c,d->", Rl.astype(jnp.complex128), l4, m4, l4, m4)
    s_dn = gam @ s_up
    sxB = jnp.einsum("ikl,k,lj->ij", epsmix, s_dn, B3)
    sxB = (sxB + sxB.T) / 2
    lhs = jnp.einsum("ij,i,j->", (E3 + SIGf * sxB).astype(jnp.complex128), m3c, m3c)
    resE = jnp.max(jnp.abs(E3 - E4)) / (1 + jnp.max(jnp.abs(E4)))
    resB = jnp.max(jnp.abs(B3 - B4)) / (1 + jnp.max(jnp.abs(B4)))
    resI = jnp.abs(lhs - KAPf * psi0_4d) / (1 + jnp.abs(psi0_4d))
    return resE, resB, resI

routes_batch = jax.jit(jax.vmap(routes))

def device_sweep(key, dev):
    k1, k2, k3, k4 = jax.random.split(key, 4)
    B = 256
    r = jax.random.uniform(k1, (B,), jnp.float64, 3.0, 12.0)
    v = jax.random.normal(k2, (B, 3), jnp.float64)
    xs = r[:, None] * v / jnp.linalg.norm(v, axis=1)[:, None]
    X = jnp.concatenate([jnp.zeros((B, 1)), xs], axis=1)
    v1 = jax.random.uniform(k3, (B, 3), jnp.float64, -1, 1)
    v2 = jax.random.uniform(k4, (B, 3), jnp.float64, -1, 1)
    rE, rB, rI = routes_batch(jax.device_put(X, dev), jax.device_put(v1, dev),
                              jax.device_put(v2, dev))
    return float(jnp.max(rE)), float(jnp.max(rB)), float(jnp.max(rI))

out = [device_sweep(jax.random.PRNGKey(81 + i), d) for i, d in enumerate(devs)]
rE = max(o[0] for o in out)
rB = max(o[1] for o in out)
rI = max(o[2] for o in out)
n_pts = len(devs) * 256
report(f"Schwarzschild (Kerr-Schild, M=1, r in [3,12]M, random frames, "
       f"{n_pts} pts): |E_3+1 - E_4D| rel {rE:.2e} < 1e-8", rE < 1e-8)
report(f"Schwarzschild: |B_3+1(Codazzi, sigma_B) - B_4D| rel {rB:.2e} < 1e-8",
       rB < 1e-8)
report(f"Schwarzschild NONPERTURBATIVE identity: "
       f"|(E + sigma eps(s)B)(m,m) - kappa psi0| rel {rI:.2e} < 1e-8", rI < 1e-8)

print(flush=True)
print(f"CONCLUSION: the EXACT dictionary is (E + ({SIG})*eps(s)B)(m,m) = "
      f"({KAP})*psi0 as pure", flush=True)
print("Weyl algebra (generic 10-dof tensor), with E, B realized from the Z4c", flush=True)
print("state via Gauss-Codazzi (on-shell vacuum, Codazzi sign sigma_B="
      f"{SIGB_found}) and", flush=True)
print("verified nonperturbatively on Schwarzschild; the superseded linearized", flush=True)
print("eqs. 16-17 are its on-shell flat limit. Physical BC rows become", flush=True)
print("conditions on U^-_AB = TT2[E + sigma eps(s)B]_AB driven to", flush=True)
print("kappa (psi0' mbar mbar + cc).", flush=True)
print(f"\nwall clock: {time.time()-T0:.1f}s", flush=True)
print(f"OVERALL: {'PASS' if OK else 'FAIL'}", flush=True)
raise SystemExit(0 if OK else 1)
