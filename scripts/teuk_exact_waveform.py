#!/usr/bin/env python3
"""teuk_exact_waveform.py — exact linear r*psi4 (2,0) of the Teukolsky
solution at FINITE radius, in the AthenaK extraction convention.

Coefficients are loaded from results/numerical/n12_psi4_exact_coeffs.json,
which is emitted (with exact symbolic provenance and gates) by
scripts/verify_n12_exact_psi4.py [ledger node N12].

    rpsi4_20_exact(t, r, X, rc, tau)
      = sum_kj c_ret[k][j] F^(k)(t-r) r^-j + c_adv[k][j] F^(k)(t+r) r^-j,
    F(u) = X exp(-((u-rc)/tau)^2),  F^(k) via physicists' Hermite H_k.
"""
import json
from pathlib import Path

import numpy as np

_COEF = json.loads(
    (Path(__file__).resolve().parent.parent /
     "results/numerical/n12_psi4_exact_coeffs.json").read_text())


def F_n(n, u, X, rc, tau):
    s = (u - rc)/tau
    H = [lambda s: 1.0 + 0*s,
         lambda s: 2*s,
         lambda s: 4*s**2 - 2,
         lambda s: s*(8*s**2 - 12),
         lambda s: (16*s**2 - 48)*s**2 + 12,
         lambda s: s*((32*s**2 - 160)*s**2 + 120),
         lambda s: ((64*s**2 - 480)*s**2 + 720)*s**2 - 120,
         lambda s: s*(((128*s**2 - 1344)*s**2 + 3360)*s**2 - 1680),
         lambda s: (((256*s**2 - 3584)*s**2 + 13440)*s**2 - 13440)*s**2
                   + 1680][n](s)
    return X*(-1.0/tau)**n*H*np.exp(-s**2)


def rpsi4_20_exact(t, r, X, rc, tau):
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    for part, sign in (("ret", -1.0), ("adv", +1.0)):
        for k, js in _COEF[part].items():
            fk = F_n(int(k), t + sign*r, X, rc, tau)
            for j, (expr, val) in js.items():
                out += val*fk/r**int(j)
    return out
