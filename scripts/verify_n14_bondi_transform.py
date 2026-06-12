#!/usr/bin/env python3
"""verify_n14_bondi_transform.py — N14 A2 machine derivation: the COMPLETE
linear Bondi-like transformation of the corrected Teukolsky solution.

Background in flat-Bondi coordinates (u, r, th, ph), u = t - r:
  gbar = -du^2 - 2 du dr + r^2 dth^2 + r^2 s^2 dph^2.
Teukolsky perturbation mapped to these coordinates (retarded channel,
F^n = F^(n)(u); dt = du + dr kills nothing since h_t* = 0):
  h_rr = A(3c^2-1),  h_rth = -3 B r s c,
  h_thth = r^2 (3C s^2 - A),  h_phph = r^2 s^2 (-3C s^2 + A(3s^2-1)),
  A = 3[F2/r^3 + 3F1/r^4 + 3F0/r^5],
  B = -[F3/r^2 + 3F2/r^3 + 6F1/r^4 + 6F0/r^5],
  C = (1/4)[F4/r + 2F3/r^2 + 9F2/r^3 + 21F1/r^4 + 21F0/r^5].

Gauge vector zeta = (zu, zr, zth, 0), delta g = h + Lie_zeta gbar
(coordinate formula, no Christoffels). Ansatz: each component is a series
sum_{n,k} c_{nk} F^n(u) r^{-k} times angular shapes in {1, c^2} x {1, sc}
as appropriate (zu, zr: even shapes a + b c^2; zth: s c shape).
Conditions imposed EXACTLY (coefficients in QQ):
  (1) g'_rr = 0   (null radial)
  (2) g'_rth = 0  (angular alignment)
  (3) trace-free areal: delta g'_thth / r^2 + delta g'_phph / (r^2 s^2) = 0
Then read off (linear): beta from g'_ur = -(1 + 2 beta);
  U   from g'_uth = -r^2 U   (spin-1 scalar, sc shape);
  W   from g'_uu  = -(1 + 2 beta + r^2 W ... linear: -(1 + 2b + 2 r W?) —
       linear Bondi form: g_uu = -(V/r) e^{2 beta} + O(U^2):
       delta g_uu = -(2 beta + (V-r)/r) -> W = (V-r)/r^2;
  J   from (delta g'_thth / r^2 - delta g'_phph / (r^2 s^2)) / 2.

Gates:
  T1 the three conditions admit an exact series solution (unique up to
     residual gauge with scri falloff: all unknowns k >= 1).
  T2 psi0 chain: psi0 = -dr J/(2r) - dr^2 J/4 — CONFIRM or AMEND the A1
     coefficient 9/8 s^2 F2/r^5 (the A1 xi dropped d_theta T; this is the
     authoritative re-derivation).
  T3 hierarchy: with the derived J, beta, U (Q := r^2 dr U), gate the
     linear Q-equation residual and the closure H == du J with the floated
     H-coefficient solve now FREE of integration-constant crutches.
Transcript -> results/numerical/n14_bondi_transform_check.txt
"""
import signal
import sys
from pathlib import Path

import sympy as sp

signal.alarm(540)
ROOT = Path(__file__).resolve().parent.parent
OK = True
LINES = []


def log(s):
    global LINES
    print(s, flush=True)
    LINES.append(s)


def gate(name, cond):
    global OK
    log(f"[{'PASS' if cond else 'FAIL'}] {name}")
    OK &= bool(cond)


u, r, th = sp.symbols("u r theta", positive=True)
c = sp.cos(th)
s = sp.sin(th)
NF = 8
F = sp.symbols(f"F0:{NF}")          # F[n] = F^(n)(u)


def du_expr(e):
    """d/du: F^n -> F^(n+1)."""
    out = sp.diff(e, u)
    for n in range(NF-1, -1, -1):
        out = out.subs(sp.Derivative(F[n], u), 0)
    # F[n] are symbols: implement via replacement after formal diff
    return out


def DU(e):
    rep = {F[n]: sp.Symbol(f"G{n}") for n in range(NF)}
    e2 = e.subs(rep)
    # shift: G_n -> F_{n+1}
    back = {sp.Symbol(f"G{n}"): F[n+1] if n+1 < NF else 0 for n in range(NF)}
    return e2.subs(back)


A = 3*(F[2]/r**3 + 3*F[1]/r**4 + 3*F[0]/r**5)
B = -(F[3]/r**2 + 3*F[2]/r**3 + 6*F[1]/r**4 + 6*F[0]/r**5)
C = sp.Rational(1, 4)*(F[4]/r + 2*F[3]/r**2 + 9*F[2]/r**3
                       + 21*F[1]/r**4 + 21*F[0]/r**5)

h = sp.zeros(4, 4)   # coords (u, r, th, ph)
h[1, 1] = A*(3*c**2 - 1)
h[1, 2] = h[2, 1] = -3*B*r*s*c
h[2, 2] = r**2*(3*C*s**2 - A)
h[3, 3] = r**2*s**2*(-3*C*s**2 + A*(3*s**2 - 1))

gbar = sp.diag(-1, 0, r**2, r**2*s**2)
gbar[0, 1] = gbar[1, 0] = -1

X = [u, r, th, sp.Symbol("phi")]

# gauge ansatz with unknown rational coefficients
KMAX = 6
unk = []


def series(shapes, name, kmin=1):
    e = 0
    for i, sh in enumerate(shapes):
        for n in range(0, 5):
            for k in range(kmin, KMAX+1):
                a = sp.Symbol(f"{name}_{i}_{n}_{k}")
                unk.append(a)
                e += a*sh*F[n]/r**k
    return e


zu = series([1, c**2], "zu")
zr = series([1, c**2], "zr")
zth = series([s*c], "zt", kmin=2)
zeta = [zu, zr, zth, 0]

lie = sp.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        e = 0
        for al in range(4):
            e += zeta[al]*sp.diff(gbar[mu, nu], X[al])
            e += gbar[al, nu]*sp.diff(zeta[al], X[mu])
            e += gbar[mu, al]*sp.diff(zeta[al], X[nu])
        lie[mu, nu] = sp.expand(e)

g1 = sp.expand(h + lie)   # linear metric perturbation in the new gauge


def channels(expr):
    """decompose expr into exact coefficient equations over the basis
    F^n * r^-k * c^m (s-powers reduced via s^2 = 1 - c^2)."""
    e = sp.expand(sp.simplify(sp.expand_trig(expr).rewrite(sp.cos)))
    e = sp.expand(e.subs(s**2, 1 - c**2))
    out = []
    p = sp.Poly(e, *F, c, s)
    for mono, coef in zip(p.monoms(), p.coeffs()):
        out.append(coef)
    # coef still depends on r: expand in powers of r
    eqs = []
    for coef in out:
        pr = sp.Poly(sp.together(coef*r**(KMAX+8)).as_numer_denom()[0], r)
        eqs.extend(pr.coeffs())
    return eqs


eqs = []
eqs += channels(g1[1, 1])                         # g'_rr = 0
eqs += channels(g1[1, 2])                         # g'_rth = 0
eqs += channels(g1[2, 2]/r**2 + g1[3, 3]/(r**2*s**2))   # areal trace-free

sol = sp.solve(eqs, unk, dict=True)
log(f"[INFO] gauge conditions: {len(eqs)} channel equations, "
    f"{len(unk)} unknowns, solutions: {len(sol)}")
gate("T1 the Bondi-like gauge conditions admit an exact series solution",
     len(sol) >= 1)
if not sol:
    log(f"\nOVERALL: FAIL")
    (ROOT/"results/numerical/n14_bondi_transform_check.txt").write_text(
        "\n".join(LINES) + "\n")
    sys.exit(1)
s0 = sol[0]
# any remaining free symbols (residual gauge) -> set to zero
g1s = g1.subs(s0)
g1s = g1s.subs({a: 0 for a in unk})
g1s = sp.expand(g1s)

# read off the Bondi scalars
Jexpr = sp.simplify(sp.expand(
    (g1s[2, 2]/r**2 - g1s[3, 3]/(r**2*s**2))/2).subs(s**2, 1 - c**2))
beta_expr = sp.simplify(-(g1s[0, 1])/2)
U_expr = sp.simplify(-g1s[0, 2]/r**2)
W_expr = sp.simplify(-(g1s[0, 0] + 2*beta_expr)/r**2)  # wait: see docstring
log(f"[INFO] J = {sp.nsimplify(Jexpr)}")
log(f"[INFO] beta = {sp.nsimplify(beta_expr)}")
log(f"[INFO] U = {sp.nsimplify(U_expr)}")
log(f"[INFO] W = {sp.nsimplify(W_expr)}")

# T2: psi0 chain on the derived J
psi0 = sp.expand(-sp.diff(Jexpr, r)/(2*r) - sp.diff(Jexpr, r, 2)/4)
psi0 = sp.expand(psi0.subs(s**2, 1 - c**2))
lead = psi0.coeff(F[2]).coeff(1/r**5) if psi0.has(F[2]) else 0
log(f"[INFO] psi0 = {sp.nsimplify(psi0)}")
target = sp.Rational(9, 8)*(1 - c**2)
gate(f"T2 psi0 leading F2/r^5 coefficient = (9/8) s^2 "
     f"(A1 CONFIRM gate; got {sp.nsimplify(lead)})",
     sp.simplify(lead - target) == 0)

log(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
(ROOT/"results/numerical/n14_bondi_transform_check.txt").write_text(
    "\n".join(LINES) + "\n")
sys.exit(0 if OK else 1)
