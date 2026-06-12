#!/usr/bin/env python3
"""verify_n14_linear_psi0_chain.py — N14 stage A1: the linear characteristic
psi0 chain closes against the N12-exact datum.

Chain (all linear order, outgoing Teukolsky (2,0), corrected N12 solution):
 1. Bondi-like J from the Cauchy metric on null slices u = t - r:
      J_naive = (h_thth - h_phph)/2  [orthonormal components]
              = (3/4) sin^2(th) [F4/r + 2F3/r^2 + 3F2/r^3 + 3F1/r^4 + 3F0/r^5]
    PLUS the angular gauge shift removing g_{r theta} (Bondi condition):
      dr xi^th = -h_{r theta}/r^2,  xi(scri) = 0,  h_{r th} = -3 B r s c,
      B = -[F3/r^2 + 3F2/r^3 + 6F1/r^4 + 6F0/r^5]
      => xi^th = 3 s c [F3/(2r^2) + F2/r^3 + (3/2)F1/r^4 + (6/5)F0/r^5]
      delta J = d_th xi^th - cot(th) xi^th = -3 s^2 R_xi
 2. J_full = (3/4) s^2 [F4/r - F2/r^3 - 3F1/r^4 - (9/5)F0/r^5]
    (the F3/r^2 term cancels EXACTLY — peeling restored).
 3. Linearized eq:psi0cce (method paper; K = 1, beta = 0):
      psi0_lin = -dr J/(2r) - (1/4) dr^2 J   [dr at fixed u]
 4. Result: psi0_lin = (3/4) s^2 [(3/2)F2/r^5 + 9F1/r^6 + 9F0/r^7];
    leading term (9/8) s^2 F2(u)/r^5 = -1 x the N12 datum
      -sqrt(27 pi/10) F2(t-r)/r^5 * 2Y20,  2Y20 = (1/4)sqrt(15/2pi) s^2,
    i.e. |coeff| = 9/8 exactly (convention map = the known psi0 sign).

Gates:
  C1 gauge cancellation: J_full has NO F3/r^2 term (exact rational算).
  C2 leading coefficient: |9/8| matches |sqrt(27pi/10) * (1/4)sqrt(15/2pi)|
     exactly (sympy).
  C3 numeric closure: psi0_lin evaluated by finite differences of J_full
     (mpmath, off the closed form) matches the term-by-term series at
     (u, r) samples to 1e-25.
  C4 subleading structure reported (F1/r^6, F0/r^7 tails are the
     finite-radius CCE-frame content; the datum is the leading peeling
     coefficient).
"""
import signal
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
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


# ---- exact series algebra: fields = dict {(n, k): Fraction} meaning
# sum c_{nk} F^(n)(u) / r^k  (retarded outgoing part; advanced dropped)
def dr(f):
    """d/dr at fixed u: only 1/r^k differentiates."""
    return {(n, k+1): c*(-k) for (n, k), c in f.items()}


J_naive = {(4, 1): Fraction(3, 4), (3, 2): Fraction(3, 2),
           (2, 3): Fraction(9, 4), (1, 4): Fraction(9, 4),
           (0, 5): Fraction(9, 4)}
# gauge shift contribution: delta J = -3 s^2 R_xi, with the s^2 already
# factored out of all entries here (every term carries sin^2 theta)
R_xi = {(3, 2): Fraction(1, 2), (2, 3): Fraction(1, 1),
        (1, 4): Fraction(3, 2), (0, 5): Fraction(6, 5)}
J_full = dict(J_naive)
for key, c in R_xi.items():
    J_full[key] = J_full.get(key, Fraction(0)) - 3*c
J_full = {k: c for k, c in J_full.items() if c != 0}
log(f"[INFO] J_full series (n=F-order, k=1/r power, coeff x sin^2): "
    f"{ {k: str(c) for k, c in sorted(J_full.items())} }")

gate("C1 gauge cancellation: no F3/r^2 term in J_full",
     (3, 2) not in J_full)

# psi0_lin = -dr J/(2r) - (1/4) dr^2 J
dJ = dr(J_full)
d2J = dr(dJ)
psi0 = {}
for (n, k), c in dJ.items():
    psi0[(n, k+1)] = psi0.get((n, k+1), Fraction(0)) - Fraction(1, 2)*c
for (n, k), c in d2J.items():
    psi0[(n, k)] = psi0.get((n, k), Fraction(0)) - Fraction(1, 4)*c
psi0 = {k: c for k, c in psi0.items() if c != 0}
log(f"[INFO] psi0_lin series: { {k: str(c) for k, c in sorted(psi0.items())} }")

lead = psi0.get((2, 5), Fraction(0))
datum_mag = sp.sqrt(27*sp.pi/10)*sp.Rational(1, 4)*sp.sqrt(15/(2*sp.pi))
gate(f"C2 leading coefficient: |{lead}| == sqrt(27pi/10)*(1/4)sqrt(15/2pi) "
     f"= {sp.nsimplify(datum_mag)} exactly",
     sp.simplify(sp.Rational(lead.numerator, lead.denominator)
                 - datum_mag) == 0)
no_lower = all(k >= 5 for (n, k) in psi0)
gate("C4 peeling: no psi0 term below 1/r^5", no_lower)
log("[INFO] C4 subleading CCE-frame tails: "
    + str({k: str(c) for k, c in sorted(psi0.items()) if k != (2, 5)}))

# ---- C3 numeric closure (mpmath FD of the closed-form J_full) ------------
mp.mp.dps = 35
RC, TAU = mp.mpf(20), mp.mpf(2)
HERM = [lambda s: 1+0*s, lambda s: 2*s, lambda s: 4*s**2-2,
        lambda s: s*(8*s**2-12), lambda s: (16*s**2-48)*s**2+12,
        lambda s: s*((32*s**2-160)*s**2+120),
        lambda s: ((64*s**2-480)*s**2+720)*s**2-120]


def Fn(n, u):
    s = (u - RC)/TAU
    return (-1/TAU)**n*HERM[n](s)*mp.e**(-s**2)


def Jnum(u, r):
    return sum(mp.mpf(c.numerator)/c.denominator*Fn(n, u)/r**k
               for (n, k), c in J_full.items())


def psi0_series(u, r):
    return sum(mp.mpf(c.numerator)/c.denominator*Fn(n, u)/r**k
               for (n, k), c in psi0.items())


worst = mp.mpf(0)
for u, r in ((mp.mpf(18), mp.mpf(41)), (mp.mpf(22), mp.mpf(60)),
             (mp.mpf(20), mp.mpf(36))):
    dJn = mp.diff(lambda rr: Jnum(u, rr), r)
    d2Jn = mp.diff(lambda rr: Jnum(u, rr), r, 2)
    lhs = -dJn/(2*r) - d2Jn/4
    rhs = psi0_series(u, r)
    rel = abs(lhs - rhs)/max(abs(rhs), mp.mpf("1e-60"))
    worst = max(worst, rel)
gate(f"C3 numeric closure of psi0_lin = -J'/(2r) - J''/4 vs series: "
     f"worst rel dev = {mp.nstr(worst, 3)} < 1e-25", worst < mp.mpf("1e-25"))

log("[INFO] convention map: psi0_lin(eq:psi0cce frame) = (-1) x the N12 "
    "datum -sqrt(27pi/10)F2(t-r)/r^5*2Y20 at leading peeling order — same "
    "sign-map class as the N12 psi4 result")
log("[INFO] stage-B consequence: the live worldtube J needs h_thth, h_phph "
    "AND the radial integral of h_{r theta}/r^2 over the characteristic "
    "domain (the ifc worldtube-transformation chain) — the angular metric "
    "alone violates peeling by an F3/r^4 psi0 term")

log(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
with open(ROOT/"results/numerical/n14_linear_check.txt", "w") as f:
    f.write("\n".join(LINES) + "\n")
sys.exit(0 if OK else 1)
