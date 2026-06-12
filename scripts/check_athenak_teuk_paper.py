#!/usr/bin/env python3
"""check_athenak_teuk_paper.py — evaluate the EXACT-parameter Teukolsky
battery (arXiv:2308.10361 Sec. V.A figure protocol): r_c = 20, tau = 2,
boundary 41, extraction r = 36, reference boundary 124.

Analytic scri formula (ledger eq:Teukolsky_bulk_psi0 chain, leading 1/r):
  r psi4 (2,0) = -sqrt(6 pi/5) F^(6)(t - r),  F = X exp(-((u-r_c)/tau)^2).

Checks:
  W1 (INFO + order-of-magnitude PASS) peak |r psi4(2,0)| vs the scri formula
     at r = 36 (r/tau = 18; near-zone terms O(tau/r) ~ 6%).
  W2 m = 0 reality: max |Im(2,0)| << real peak.
  W3 boundary error vs causally disconnected reference:
     E(run) = ||psi4_run - psi4_ref||_2 / ||psi4_ref||_2 on a common time
     grid. X <= 1e-4: PASS iff both < 0.05 and |E(ccm)-E(somm)| < 0.1 max(E)
     (the paper's perturbative indistinguishability). X = 2: INFO.
  W5 boundary error below the numerical-error estimate (the paper's bottom-
     panel conclusion): ||ccm - ref|| < ||ccm - ccm_mid|| (two-resolution
     difference at the same boundary treatment). X <= 1e-4 only.
  W4 constraints bounded on the ccm run: final L2-H < 10x initial.
"""
import sys
from pathlib import Path

import numpy as np

base = Path(sys.argv[1]); AMP = float(sys.argv[2])
RC, TAU, REXT = 20.0, 2.0, 36.0
OK = True

def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    OK &= bool(cond)

def F_n(n, u):
    s = (u - RC) / TAU
    H = [1, 2*s, 4*s**2 - 2, 8*s**3 - 12*s, 16*s**4 - 48*s**2 + 12,
         32*s**5 - 160*s**3 + 120*s,
         64*s**6 - 480*s**4 + 720*s**2 - 120][n]
    return AMP * (-1.0/TAU)**n * H * np.exp(-s**2)

def wf(run):
    fr = np.loadtxt(base/run/"waveforms"/f"rpsi4_real_{int(REXT):04d}.txt")
    fi = np.loadtxt(base/run/"waveforms"/f"rpsi4_imag_{int(REXT):04d}.txt")
    return fr[:, 0], fr[:, 3], fi[:, 3]   # t, Re(2,0), Im(2,0)

# common time grid = ref grid clipped to all runs' coverage
t_ref, ref_re, _ = wf("ref")
runs = {}
for name in ("somm", "somm_mid", "ccm", "ccm_mid", "ccm_low"):
    tr, re_, im_ = wf(name)
    runs[name] = (tr, re_, im_)
tmax = min([t_ref[-1]] + [runs[n][0][-1] for n in runs])
t = t_ref[t_ref <= tmax]
ref = np.interp(t, t_ref, ref_re)
re = {n: np.interp(t, runs[n][0], runs[n][1]) for n in runs}

ana = -np.sqrt(6*np.pi/5) * F_n(6, t - REXT)
win = (t > RC + REXT - 4*TAU) & (t < RC + REXT + 4*TAU)

for name, y in (("somm", re["somm"]), ("ccm", re["ccm"]), ("ref", ref)):
    ratio = np.max(np.abs(y[win])) / np.max(np.abs(ana[win]))
    print(f"[INFO] W1 {name}: peak |r psi4(2,0)| / scri formula = {ratio:.3f}"
          f" (extraction r = {REXT:g}, tau/r = {TAU/REXT:.3f})")
if AMP <= 1e-4:
    # the deficit from 1 is FD truncation (ratio climbs monotonically with
    # resolution: 0.29 / 0.47 / 0.61 at 104/132/164^3, 2nd-order stencils at
    # ~4 points per tau), not the finite extraction radius
    report("W1 peak amplitudes within an order of magnitude of the scri "
           "formula", all(0.1 < np.max(np.abs(y[win]))/np.max(np.abs(ana[win]))
                          < 10 for y in (re["somm"], re["ccm"], ref)))
else:
    print("[INFO] W1 gate skipped: the linear scri formula does not apply "
          f"at X = {AMP:g} (paper compares against a nonlinear reference)")

im_peak = max(np.max(np.abs(runs["somm"][2])), np.max(np.abs(runs["ccm"][2])))
im_ratio = im_peak / np.max(np.abs(ref))
report(f"W2 m=0 reality: max |Im(2,0)| / peak = {im_ratio:.2e} < 1e-3",
       im_ratio < 1e-3)

den = np.linalg.norm(ref)
Es = np.linalg.norm(re["somm"] - ref) / den
Ec = np.linalg.norm(re["ccm"] - ref) / den
gap = abs(Ec - Es) / max(Ec, Es, 1e-300)
if AMP <= 1e-4:
    report(f"W3 boundary error vs reference: E(somm) = {Es:.3e}, "
           f"E(ccm) = {Ec:.3e} — both < 0.05, gap {gap:.1%} < 10% "
           f"(perturbative indistinguishability)",
           Es < 0.05 and Ec < 0.05 and gap < 0.10)
else:
    print(f"[INFO] W3 nonlinear regime: E(somm) = {Es:.3e}, E(ccm) = {Ec:.3e}"
          f" (gap {gap:.1%})")

res_est_c = np.linalg.norm(re["ccm"] - re["ccm_mid"]) / den
res_est_s = np.linalg.norm(re["somm"] - re["somm_mid"]) / den
print(f"[INFO] numerical-error estimates: ccm |164-132| = {res_est_c:.3e}, "
      f"somm |164-132| = {res_est_s:.3e}")
if AMP <= 1e-4:
    report(f"W5 boundary error below numerical error: E(ccm) = {Ec:.3e} < "
           f"|ccm164 - ccm132| = {res_est_c:.3e}", Ec < res_est_c)

def hnorm(run):
    d = np.loadtxt(base/run/"teuk.user.hst")
    return d[0, 4], d[-1, 4]
h0, h1 = hnorm("ccm")
report(f"W4 constraints bounded (ccm run): L2-H {h0:.2e} -> {h1:.2e} "
       f"(< 10x initial)", h1 < 10*h0)

print(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
sys.exit(0 if OK else 1)
