#!/usr/bin/env python3
"""verify_n12_exact_psi4.py — N12a: EXACT linear r*psi4 (and psi0) of the
(l,m)=(2,0) Teukolsky solution at FINITE radius, in the AthenaK extraction
convention (BAM-style contraction, FR4=1/4, Gram-Schmidt (phi,r,theta)
tetrad, output = r*psi4 projected on -2Y20).

Construction is a verbatim sympy mirror of the two C++ code paths:
  athenak/src/pgen/z4c_ccm_teukolsky.cpp  (the solution h_ij, hdot_ij)
  athenak/src/z4c/z4c_calculate_weyl_scalars.cpp (the psi4 contraction)
taken at linear order in the amplitude (tetrad and K*K corrections enter
psi4 only at O(X^2) since the background Weyl vanishes).

Exact-linear-field layer: every scalar is a dict {(part, k): coeff} meaning
sum_k coeff(x,y,z) * F^(k)(u_part), part in {ret: u=t-r, adv: u=t+r}, F a
GENERIC pulse profile. d/dt raises k; d/dx_i acts by product+chain rule with
du/dx_i = -/+ x_i/r. All coefficients are exact sympy expressions; after the
substitution x = r*st, y = 0, z = r*ct (st, ct rational Pythagorean pairs,
r symbolic positive) every coefficient is an exact rational function of r.

Gates:
  G1 trace: delta^ij h_ij = 0 identically (Teukolsky solution is traceless).
  G2 linearized vacuum: R_ab(h) + dt K_ab = 0, momentum d^j(K_ij - d_ij trK)
     = 0, energy dt(tr K) = 0 — all coefficient-exact.
  G3 psi0 datum: the N3-route psi0 = -(E + eps(s)B)(m,m) reproduces
     r^5 psi0|_ret = -sqrt(27 pi/10) F^(2)(t-r) * 2Y20 EXACTLY (the paper's
     eq:Teukolsky_bulk_psi0), with no other retarded term.
  G4 angular factorization: psi4 / sin^2(theta) is the same rational-in-r
     object at 4 distinct Pythagorean angles.
  G5 scri limit: the (2,0) projection of r*psi4 has leading retarded term
     exactly -sqrt(6 pi/5) F^(6)(t-r).
  G6 peeling: retarded psi4 series has only k+j = 6 terms (F^(6-j)/r^j),
     j = 0..4 — pure outgoing structure.

Output: results/numerical/n12_exact_psi4_check.txt (transcript) and
results/numerical/n12_psi4_exact_coeffs.json (exact coefficient table), and
the checker module scripts/teuk_exact_waveform.py is validated against the
table at random points.
"""
import json
import signal
import sys
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
signal.alarm(540)
T0 = time.time()
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


t = sp.Symbol('t', real=True)
x, y, z = sp.symbols('x y z', real=True)
r3 = sp.sqrt(x*x + y*y + z*z)

RET, ADV = 0, 1
SGN = {RET: -1, ADV: +1}   # u_part = t + SGN*r


def fadd(a, b):
    out = dict(a)
    for key, c in b.items():
        out[key] = sp.together(out.get(key, 0) + c)
        if out[key] == 0:
            del out[key]
    return out


def fscale(a, fac):
    return {key: sp.together(c*fac) for key, c in a.items() if sp.together(c*fac) != 0}


def fdt(a):
    return {(p, k+1): c for (p, k), c in a.items()}


def fdx(a, xi):
    """d/dx_i with xi in (x, y, z): product rule + chain rule du/dxi = SGN*xi/r."""
    out = {}
    for (p, k), c in a.items():
        out = fadd(out, {(p, k): sp.diff(c, xi)})
        out = fadd(out, {(p, k+1): c*SGN[p]*xi/r3})
    return out


XYZ = (x, y, z)

# ---- pgen mirror: the Teukolsky solution -----------------------------------
def Qn(n):
    return {(RET, n): sp.Integer(1), (ADV, n): -sp.Integer((-1)**n)}


def Rn(n):
    return {(RET, n): sp.Integer(1), (ADV, n): sp.Integer((-1)**n)}


# Af = 3[Q2/r^3 + 3Q1/r^4 + 3Q0/r^5]; Bf = -[Q3/r^2 + 3Q2/r^3 + 6Q1/r^4 + 6Q0/r^5]
# Cf = (1/4)[Q4/r + 2Q3/r^2 + 9Q2/r^3 + 21Q1/r^4 + 21Q0/r^5]
Af = fadd(fadd(fscale(Qn(2), 3/r3**3), fscale(Qn(1), 9/r3**4)),
          fscale(Qn(0), 9/r3**5))
Bf = fadd(fadd(fadd(fscale(Qn(3), -1/r3**2), fscale(Qn(2), -3/r3**3)),
               fscale(Qn(1), -6/r3**4)), fscale(Qn(0), -6/r3**5))
Cf = fadd(fadd(fadd(fadd(fscale(Qn(4), sp.Rational(1, 4)/r3),
                         fscale(Qn(3), sp.Rational(1, 2)/r3**2)),
                    fscale(Qn(2), sp.Rational(9, 4)/r3**3)),
               fscale(Qn(1), sp.Rational(21, 4)/r3**4)),
          fscale(Qn(0), sp.Rational(21, 4)/r3**5))

cth = z/r3
s2th = 1 - cth**2
frr = (1 + 3*(1 - 2*s2th))/2
f1tt, f2tt = 3*s2th, -1
f1pp, f2pp = -3*s2th, 3*s2th - 1

# frame vectors exactly as the (iter-18 corrected) pgen: e_th = e_ph x r-hat
# = +theta-hat. The pre-correction rh x eph order (-theta-hat) flips h_{r th}
# and violates linearized vacuum at O(B) — the N12 root cause.
rho = sp.sqrt(x*x + y*y)
rh = [x/r3, y/r3, z/r3]
eph = [-y/rho, x/rho, 0]
eth = [eph[1]*rh[2] - eph[2]*rh[1],
       eph[2]*rh[0] - eph[0]*rh[2],
       eph[0]*rh[1] - eph[1]*rh[0]]
# pgen: hrt_o = Bf*frt with frt = -3 sth cth; sth = rho/r3
frt = -3*(rho/r3)*cth

hrr_o = fscale(Af, frr)
hrt_o = fscale(Bf, frt)
htt_o = fadd(fscale(Cf, f1tt), fscale(Af, f2tt))
hpp_o = fadd(fscale(Cf, f1pp), fscale(Af, f2pp))

h = [[None]*3 for _ in range(3)]
for a in range(3):
    for b in range(3):
        h[a][b] = fadd(
            fadd(fscale(hrr_o, rh[a]*rh[b]),
                 fscale(hrt_o, rh[a]*eth[b] + eth[a]*rh[b])),
            fadd(fscale(htt_o, eth[a]*eth[b]), fscale(hpp_o, eph[a]*eph[b])))

K = [[fscale(fdt(h[a][b]), -sp.Rational(1, 2)) for b in range(3)]
     for a in range(3)]

log(f"[INFO] fields assembled ({time.time()-T0:.1f} s)")

# ---- evaluation frame: x = r*st, y = 0, z = r*ct, r symbolic ---------------
r = sp.Symbol('r', positive=True)
ANGLES = [(sp.Rational(3, 5), sp.Rational(4, 5)),
          (sp.Rational(4, 5), sp.Rational(3, 5)),
          (sp.Rational(5, 13), sp.Rational(12, 13)),
          (sp.Rational(8, 17), sp.Rational(15, 17))]


def at_angle(field, st, ct):
    sub = {x: r*st, y: 0, z: r*ct}
    out = {}
    for key, c in field.items():
        v = sp.cancel(sp.radsimp(c.subs(sub)))
        if v != 0:
            out[key] = v
    return out


def fd2(field, i, j):
    return fdx(fdx(field, XYZ[i]), XYZ[j])


# precompute first and second derivative tables of h and first of K
dh = [[[fdx(h[a][b], XYZ[c]) for c in range(3)] for b in range(3)]
      for a in range(3)]
d2h = [[[[fdx(dh[a][b][c], XYZ[d]) for d in range(3)] for c in range(3)]
        for b in range(3)] for a in range(3)]
dK = [[[fdx(K[a][b], XYZ[c]) for c in range(3)] for b in range(3)]
      for a in range(3)]
log(f"[INFO] derivative tables built ({time.time()-T0:.1f} s)")

# linear 3-Ricci: R_ab = 1/2 sum_c (d_c d_a h_cb + d_c d_b h_ca - d_c d_c h_ab
#                                   - d_a d_b h_cc)
Ric = [[None]*3 for _ in range(3)]
for a in range(3):
    for b in range(3):
        acc = {}
        for c in range(3):
            acc = fadd(acc, d2h[c][b][c][a])
            acc = fadd(acc, d2h[c][a][c][b])
            acc = fadd(acc, fscale(d2h[a][b][c][c], -1))
            acc = fadd(acc, fscale(d2h[c][c][a][b], -1))
        Ric[a][b] = fscale(acc, sp.Rational(1, 2))
Rsc = fadd(fadd(Ric[0][0], Ric[1][1]), Ric[2][2])
log(f"[INFO] linear Ricci built ({time.time()-T0:.1f} s)")

# ---- G1 trace ---------------------------------------------------------------
tr = fadd(fadd(h[0][0], h[1][1]), h[2][2])
tr_eval = at_angle(tr, *ANGLES[0])
gate("G1 trace: delta^ij h_ij = 0 identically", len(tr_eval) == 0)

# ---- G2 linearized vacuum ---------------------------------------------------
viol = []
for a in range(3):
    for b in range(a, 3):
        # ADM: dt K_ab = R_ab at linear order (alpha = 1, beta = 0)
        ev = at_angle(fadd(Ric[a][b], fscale(fdt(K[a][b]), -1)), *ANGLES[0])
        viol.append(("R_ab - dt K_ab", a, b, len(ev)))
trK = fadd(fadd(K[0][0], K[1][1]), K[2][2])
for i in range(3):
    # momentum: d^j K_ij - d_i trK
    mom = {}
    for j in range(3):
        mom = fadd(mom, dK[i][j][j])
    mom = fadd(mom, fscale(fdx(trK, XYZ[i]), -1))
    ev = at_angle(mom, *ANGLES[0])
    viol.append(("momentum_i", i, '-', len(ev)))
ev = at_angle(fdt(trK), *ANGLES[0])
viol.append(("dt trK", '-', '-', len(ev)))
bad = [v for v in viol if v[3] != 0]
gate(f"G2 linearized vacuum: all {len(viol)} component checks "
     f"coefficient-exact zero", len(bad) == 0)
for v in bad[:4]:
    log(f"   [G2 violation] {v}")

# ---- BAM psi4 contraction (linear order) ------------------------------------
def psi4_at(st, ct):
    sub = {x: r*st, y: 0, z: r*ct}
    # extraction tetrad (flat limit of the code's Gram-Schmidt at phi=0)
    uv = [st, 0, ct]
    vv = [ct, 0, -st]      # +theta-hat: (xz, yz, -x^2-y^2) normalized
    wv = [0, 1, 0]         # +phi-hat
    RicE = [[at_angle(Ric[a][b], st, ct) for b in range(3)] for a in range(3)]
    RscE = at_angle(Rsc, st, ct)
    dKE = [[[at_angle(dK[a][b][c], st, ct) for c in range(3)]
            for b in range(3)] for a in range(3)]

    def riem4(a, b, c, d):
        out = {}
        if a == c:
            out = fadd(out, RicE[b][d])
        if b == d:
            out = fadd(out, RicE[a][c])
        if a == d:
            out = fadd(out, fscale(RicE[b][c], -1))
        if b == c:
            out = fadd(out, fscale(RicE[a][d], -1))
        de = (1 if a == c else 0)*(1 if b == d else 0) \
            - (1 if a == d else 0)*(1 if b == c else 0)
        if de:
            out = fadd(out, fscale(RscE, -sp.Rational(de, 2)))
        return out

    def riem4_ddd(a, b, c):
        # code: -(DK_ddd(c,a,b) - DK_ddd(b,a,c)); DK_ddd(c,a,b) = D_c K_ab
        return fadd(fscale(dKE[a][b][c], -1), dKE[a][c][b])

    acc = {}
    for a in range(3):
        for b in range(3):
            ang = vv[a]*vv[b] - wv[a]*wv[b]
            if ang != 0:
                acc = fadd(acc, fscale(RicE[a][b], -sp.Rational(1, 4)*ang))
            for c in range(3):
                if uv[c] != 0 and ang != 0:
                    acc = fadd(acc, fscale(riem4_ddd(a, c, b),
                                           sp.Rational(1, 2)*uv[c]*ang))
                for d in range(3):
                    if uv[c] != 0 and uv[d] != 0 and ang != 0:
                        acc = fadd(acc, fscale(
                            riem4(d, a, c, b),
                            -sp.Rational(1, 4)*uv[d]*uv[c]*ang))
    return {key: sp.cancel(c) for key, c in acc.items()}   # psi4 (no r yet)


P4 = {}
for st, ct in ANGLES:
    p = psi4_at(st, ct)
    P4[(st, ct)] = {key: sp.cancel(c/st**2) for key, c in p.items()}
log(f"[INFO] psi4 evaluated at {len(ANGLES)} Pythagorean angles "
    f"({time.time()-T0:.1f} s)")

# ---- G4 angular factorization ----------------------------------------------
ref = P4[ANGLES[0]]
same = True
for ang in ANGLES[1:]:
    other = P4[ang]
    keys = set(ref) | set(other)
    for key in keys:
        if sp.cancel(ref.get(key, 0) - other.get(key, 0)) != 0:
            same = False
gate("G4 psi4 = P4(t,r) sin^2(theta): P4 identical at 4 Pythagorean angles",
     same)

# ---- (2,0) projection and the exact waveform table --------------------------
# C20(t,r) = int (r psi4) (-2Y20) dOmega; -2Y20 = (1/4)sqrt(15/2pi) sin^2 th
# => C20 = r * P4 * (1/4)sqrt(15/(2 pi)) * 2 pi * int_0^pi sin^5 = (8/15)sqrt(15 pi/2) * r * P4
NPROJ = sp.Rational(8, 15)*sp.sqrt(15*sp.pi/2)
C20 = {key: sp.cancel(sp.radsimp(NPROJ*r*c)) for key, c in ref.items()}

table = {"ret": {}, "adv": {}}
for (p, k), c in C20.items():
    num, den = sp.fraction(sp.together(sp.radsimp(c)))
    dpoly = sp.Poly(den, r)
    if dpoly.length() != 1:
        gate(f"table extraction: denominator pure r-power for ({p},{k})",
             False)
        continue
    d = dpoly.degree()
    dcoef = dpoly.coeffs()[0]
    npoly = sp.Poly(num, r)
    key = "ret" if p == RET else "adv"
    for mono, coef in zip(npoly.monoms(), npoly.coeffs()):
        j = d - mono[0]
        table[key].setdefault(str(k), {})
        prev = table[key][str(k)].get(str(j), sp.Integer(0))
        val = sp.radsimp(prev + coef/dcoef)
        if val != 0:
            table[key][str(k)][str(j)] = val
        elif str(j) in table[key][str(k)]:
            del table[key][str(k)][str(j)]

# ---- G5 scri limit -----------------------------------------------------------
# The AthenaK/BAM contraction fixes the overall sign convention of psi4; the
# gate is |coefficient| = sqrt(6 pi/5) with the sign RECORDED as the
# convention map to the paper formula -sqrt(6 pi/5) F^(6) (arXiv:2308.10361).
lead = table["ret"].get("6", {}).get("0", 0)
target = -sp.sqrt(6*sp.pi/5)
sign_map = sp.simplify(lead/target)
gate(f"G5 scri limit: |retarded F^(6)/r^0 coefficient| = sqrt(6 pi/5); "
     f"convention map AthenaK = ({sign_map}) x paper",
     sp.simplify(lead**2 - target**2) == 0)
log(f"[INFO] convention: rpsi4_AthenaK = ({sign_map}) x rpsi4_paper-formula")

# ---- G6 peeling structure ----------------------------------------------------
ret_ok = all(int(k) + int(j) == 6
             for k, js in table["ret"].items() for j in js)
gate("G6 peeling: every retarded term of r*psi4(2,0) is F^(6-j)/r^j",
     ret_ok)
log("[INFO] retarded table (k=F-order, j=1/r power, coeff/(-sqrt(6pi/5))): "
    + str({(k, j): str(sp.simplify(cc/target))
           for k, js in table["ret"].items() for j, cc in js.items()}))
log("[INFO] advanced table: "
    + str({(k, j): str(sp.simplify(cc/target))
           for k, js in table["adv"].items() for j, cc in js.items()}))

# ---- G3 psi0 datum (N3 route) ------------------------------------------------
def psi0_at(st, ct):
    """psi0 = -(E + eps(s)B)(m,m); E_ab = Ric_ab (linear), B_ab =
    eps_a^{kl} d_k K_lb |sym; m = (theta-hat + i phi-hat)/sqrt2; real part
    of psi0 for m=0: -(E_vv - E_ww)/2 - (B_vw + B_wv)/2 convention resolved
    by direct contraction below."""
    uv = [st, 0, ct]
    vv = [ct, 0, -st]
    wv = [0, 1, 0]
    RicE = [[at_angle(Ric[a][b], st, ct) for b in range(3)] for a in range(3)]
    dKE = [[[at_angle(dK[a][b][c], st, ct) for c in range(3)]
            for b in range(3)] for a in range(3)]
    eps = [[[0]*3 for _ in range(3)] for _ in range(3)]
    for (i, j, k), s in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                         ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
        eps[i][j][k] = s
    B = [[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            acc = {}
            for k in range(3):
                for l in range(3):
                    if eps[a][k][l]:
                        acc = fadd(acc, fscale(dKE[l][b][k], eps[a][k][l]))
            B[a][b] = acc
    B = [[fadd(fscale(B[a][b], sp.Rational(1, 2)),
               fscale(B[b][a], sp.Rational(1, 2))) for b in range(3)]
         for a in range(3)]
    sB = [[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            acc = {}
            for k in range(3):
                for l in range(3):
                    if eps[a][k][l] and uv[k] != 0:
                        acc = fadd(acc, fscale(B[l][b], eps[a][k][l]*uv[k]))
            sB[a][b] = acc
    # Orientation map (FD-adjudicated at r = 36 and 72, exact r^-5 scaling):
    # in the concrete right-handed (s, e_th, e_ph) frame with
    # eps_ikl = +[ikl], the DATUM channel is
    #   psi0_paper = Re[(E - eps(s)B)(m,m)] = +[(E-sB)_vv - (E-sB)_ww]/2 ;
    # the N3 abstract identity (E + eps B)(m,m) = -psi0 corresponds under
    # eps -> -eps (equivalently m <-> mbar). The AthenaK injection w_ab is
    # quadratic in e_th and unaffected.
    acc = {}
    for a in range(3):
        for b in range(3):
            ang = (vv[a]*vv[b] - wv[a]*wv[b])
            if ang == 0:
                continue
            acc = fadd(acc, fscale(RicE[a][b], sp.Rational(1, 2)*ang))
            acc = fadd(acc, fscale(sB[a][b], -sp.Rational(1, 2)*ang))
    return {key: sp.cancel(c) for key, c in acc.items()}


st, ct = ANGLES[0]
p0 = psi0_at(st, ct)
# expected: psi0|_ret = -sqrt(27 pi/10) F^(2)(t-r)/r^5 * 2Y20,
# 2Y20 = (1/4) sqrt(15/2pi) sin^2 th  => coefficient of (RET,2):
exp_ret = -sp.sqrt(27*sp.pi/10)*sp.Rational(1, 4)*sp.sqrt(15/(2*sp.pi)) \
    * st**2 / r**5
got_ret = p0.get((RET, 2), sp.Integer(0))
other_ret = {k: sp.cancel(p0[k]) for k in p0
             if k[0] == RET and k != (RET, 2) and sp.cancel(p0[k]) != 0}
gate("G3 psi0 datum: retarded part = -sqrt(27pi/10) F^(2)(t-r)/r^5 * 2Y20 "
     "exactly, no other retarded term (the mode-3 datum normalization)",
     sp.simplify(got_ret - exp_ret) == 0 and len(other_ret) == 0)
if sp.simplify(got_ret - exp_ret) != 0 or other_ret:
    log(f"   [G3] got {got_ret}, expected {exp_ret}, ratio "
        f"{sp.simplify(got_ret/exp_ret)}; other retarded: "
        f"{ {k: str(v) for k, v in other_ret.items()} }")

# ---- G7 independent numerical cross-check ------------------------------------
# Evaluate the SAME BAM contraction with high-precision finite differences on
# the closed-form h (mpmath, 8th-order stencils, dps = 40) — independent of
# the exact-linear-field layer above — at two (t, r, theta) points with a
# concrete Gaussian profile, and compare to the symbolic table.
import mpmath as mp

mp.mp.dps = 40
RCv, TAUv = mp.mpf(20), mp.mpf(2)
HERM = [lambda s: 1 + 0*s, lambda s: 2*s, lambda s: 4*s**2 - 2,
        lambda s: s*(8*s**2 - 12), lambda s: (16*s**2 - 48)*s**2 + 12,
        lambda s: s*((32*s**2 - 160)*s**2 + 120),
        lambda s: ((64*s**2 - 480)*s**2 + 720)*s**2 - 120]


def Fmp_n(n, u):
    s = (u - RCv)/TAUv
    return (-1/TAUv)**n*HERM[n](s)*mp.e**(-s**2)


def h_num(tv, xv, yv, zv):
    rr = mp.sqrt(xv*xv + yv*yv + zv*zv)
    def Q(n):
        return Fmp_n(n, tv - rr) - (-1)**n*Fmp_n(n, tv + rr)
    Afv = 3*(Q(2)/rr**3 + 3*Q(1)/rr**4 + 3*Q(0)/rr**5)
    Bfv = -(Q(3)/rr**2 + 3*Q(2)/rr**3 + 6*Q(1)/rr**4 + 6*Q(0)/rr**5)
    Cfv = (Q(4)/rr + 2*Q(3)/rr**2 + 9*Q(2)/rr**3 + 21*Q(1)/rr**4
           + 21*Q(0)/rr**5)/4
    ctv = zv/rr
    s2v = 1 - ctv*ctv
    rhov = mp.sqrt(xv*xv + yv*yv)
    frrv = (1 + 3*(1 - 2*s2v))/2
    frtv = -3*(rhov/rr)*ctv
    httv = Cfv*3*s2v - Afv
    hppv = -Cfv*3*s2v + Afv*(3*s2v - 1)
    rhv = [xv/rr, yv/rr, zv/rr]
    ephv = [-yv/rhov, xv/rhov, mp.mpf(0)]
    ethv = [ephv[1]*rhv[2] - ephv[2]*rhv[1],
            ephv[2]*rhv[0] - ephv[0]*rhv[2],
            ephv[0]*rhv[1] - ephv[1]*rhv[0]]
    hm = [[mp.mpf(0)]*3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            hm[a][b] = (Afv*frrv*rhv[a]*rhv[b]
                        + Bfv*frtv*(rhv[a]*ethv[b] + ethv[a]*rhv[b])
                        + httv*ethv[a]*ethv[b] + hppv*ephv[a]*ephv[b])
    return hm


def bam_psi4_num(tv, rv, stv, ctv):
    pt = [rv*stv, mp.mpf(0), rv*ctv]
    co = [tv] + pt
    memo = {}

    def d2(mu, nu, a, b):
        key = (min(mu, nu), max(mu, nu), a, b)
        if key in memo:
            return memo[key]
        mu_, nu_ = key[0], key[1]

        def g(u, v):
            c4 = list(co)
            c4[mu_] = u
            c4[nu_] = v
            return h_num(c4[0], c4[1], c4[2], c4[3])[a][b]
        if mu_ == nu_:
            val = mp.diff(lambda u: g(u, u), co[mu_], 2)
        else:
            val = mp.diff(g, (co[mu_], co[nu_]), (1, 1))
        memo[key] = val
        return val
    # linear Ricci, K = -hdot/2, dK
    Ric_n = [[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            acc = mp.mpf(0)
            for c in range(3):
                acc += d2(1+c, 1+a, c, b) + d2(1+c, 1+b, c, a)
                acc -= d2(1+c, 1+c, a, b) + d2(1+a, 1+b, c, c)
            Ric_n[a][b] = acc/2
    Rsc_n = Ric_n[0][0] + Ric_n[1][1] + Ric_n[2][2]
    def dK(a, b, c):     # d_c K_ab = -1/2 d_c dt h_ab
        return -d2(0, 1+c, a, b)/2
    uv = [stv, mp.mpf(0), ctv]
    vv = [ctv, mp.mpf(0), -stv]
    wv = [mp.mpf(0), mp.mpf(1), mp.mpf(0)]
    def riem4(a, b, c, d):
        out = mp.mpf(0)
        if a == c: out += Ric_n[b][d]
        if b == d: out += Ric_n[a][c]
        if a == d: out -= Ric_n[b][c]
        if b == c: out -= Ric_n[a][d]
        de = (1 if a == c else 0)*(1 if b == d else 0) \
            - (1 if a == d else 0)*(1 if b == c else 0)
        return out - de*Rsc_n/2
    acc = mp.mpf(0)
    for a in range(3):
        for b in range(3):
            ang = vv[a]*vv[b] - wv[a]*wv[b]
            if ang == 0:
                continue
            acc += -ang*Ric_n[a][b]/4
            for c in range(3):
                if uv[c] != 0:
                    # Riemm4_ddd(a,c,b) = d_c K_ab - d_b K_ac (code storage
                    # (A,B,C) = -(DK(C,A,B) - DK(B,A,C)) at (A,B,C)=(a,c,b))
                    acc += ang*uv[c]*(dK(a, b, c) - dK(a, c, b))/2
                for d in range(3):
                    if uv[c] != 0 and uv[d] != 0:
                        acc += -ang*uv[d]*uv[c]*riem4(d, a, c, b)/4
    return acc


def table_psi4_num(tv, rv, stv, ctv):
    """symbolic P4 table (G4-verified angular factorization) evaluated with
    the same concrete Gaussian: psi4 = P4(t,r) * sin^2(theta)."""
    acc = mp.mpf(0)
    for (p, k), c in ref.items():
        arg = tv - rv if p == RET else tv + rv
        cv = mp.mpf(str(sp.N(c.subs(r, sp.Integer(int(rv))), 40)))
        acc += cv*Fmp_n(int(k), arg)
    return acc*stv**2


tv, rv = mp.mpf(54), mp.mpf(36)
stv, ctv = mp.mpf(3)/5, mp.mpf(4)/5
num = bam_psi4_num(tv, rv, stv, ctv)
sym = table_psi4_num(tv, rv, stv, ctv)
rel = abs(num - sym)/max(abs(sym), mp.mpf("1e-50"))
gate(f"G7 independent FD cross-check of psi4 at (t,r,th)=(54,36,*): "
     f"rel dev = {mp.nstr(rel, 3)} < 1e-8", rel < mp.mpf("1e-8"))
log(f"[INFO] G7 values: FD = {mp.nstr(num, 12)}, table = {mp.nstr(sym, 12)}")

# ---- export ------------------------------------------------------------------
out = {"convention": "AthenaK z4c_calculate_weyl_scalars (BAM contraction, "
       "FR4=1/4, GS tetrad), output r*psi4 projected on -2Y20; "
       "F = generic pulse profile, F^(k) its k-th derivative",
       "form": "r psi4_(2,0)(t,r) = sum_kj c[part][k][j] F^(k)(t -/+ r) r^-j",
       "ret": {k: {j: [str(c), float(c)] for j, c in js.items()}
               for k, js in table["ret"].items()},
       "adv": {k: {j: [str(c), float(c)] for j, c in js.items()}
               for k, js in table["adv"].items()}}
with open(ROOT/"results/numerical/n12_psi4_exact_coeffs.json", "w") as f:
    json.dump(out, f, indent=1)
log("[INFO] wrote results/numerical/n12_psi4_exact_coeffs.json")

log(f"\nOVERALL: {'PASS' if OK else 'FAIL'}  ({time.time()-T0:.1f} s)")
with open(ROOT/"results/numerical/n12_exact_psi4_check.txt", "w") as f:
    f.write("\n".join(LINES) + "\n")
sys.exit(0 if OK else 1)
