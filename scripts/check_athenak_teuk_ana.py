#!/usr/bin/env python3
"""check_athenak_teuk_ana.py — N12b gate: the AthenaK waveform converges to
the EXACT analytic finite-radius r*psi4(2,0) (scripts/teuk_exact_waveform.py,
derived and gated by verify_n12_exact_psi4.py).

Window t in [40, 58] (pulse passage, boundary-clean; see battery header).
Each rung is compared against the exact waveform AT ITS OWN sample times —
never interpolated between rungs (linear interpolation across the snapped
output cadence, ~0.33 at L2, injects O(dt^2) error far above the 6th-order
simulation error and was the dominant term in the first checker version).

  A1 exact match at the finest rung: E(L3) < 0.005 (measured 2.0e-3; the
     2nd-order baseline at h = 1/2 had a 39% peak deficit).
  A2 measured convergence order on both rung pairs >= 3.5 (6th-order
     stencils; measured 5.4 / 5.6).
  A3 finest-rung fidelity: best-fit amplitude within 1 +- 0.002 and
     best-fit time shift |dt*| <= dt_step(L3)/2 — the waveform IS the
     analytic one, not a rescaled or shifted cousin.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from teuk_exact_waveform import rpsi4_20_exact  # noqa: E402

base = Path(sys.argv[1])
X, RC, TAU, REXT = 1.0e-5, 20.0, 2.0, 36.0
H = {"L1": 0.5, "L2": 82.0/248.0, "L3": 82.0/296.0}
OK = True


def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    OK &= bool(cond)


E = {}
fit = {}
for run in ("L1", "L2", "L3"):
    d = np.loadtxt(base/run/"waveforms"/"rpsi4_real_0036.txt")
    t, y = d[:, 0], d[:, 3]
    win = (t >= 40.0) & (t <= 58.0)
    t, y = t[win], y[win]
    ex = rpsi4_20_exact(t, REXT, X, RC, TAU)
    den = np.linalg.norm(ex)
    E[run] = np.linalg.norm(y - ex)/den
    best = (np.inf, 0.0, 1.0)
    for dt in np.arange(-0.5, 0.5, 0.002):
        exs = rpsi4_20_exact(t - dt, REXT, X, RC, TAU)
        a = np.dot(y, exs)/np.dot(exs, exs)
        e = np.linalg.norm(y - a*exs)/den
        if e < best[0]:
            best = (e, dt, a)
    fit[run] = best
    print(f"[INFO] {run} (h = {H[run]:.4f}): E = {E[run]:.4e}, best-fit "
          f"shift = {best[1]:+.3f}, amplitude = {best[2]:.4f}, "
          f"E_after-fit = {best[0]:.2e}")

report(f"A1 exact match at the finest rung: E(L3) = {E['L3']:.3e} < 0.005",
       E["L3"] < 0.005)

p12 = np.log(E["L1"]/E["L2"])/np.log(H["L1"]/H["L2"])
p23 = np.log(E["L2"]/E["L3"])/np.log(H["L2"]/H["L3"])
report(f"A2 measured convergence order to the ANALYTIC waveform: "
       f"p(L1,L2) = {p12:.2f}, p(L2,L3) = {p23:.2f}, both >= 3.5",
       p12 >= 3.5 and p23 >= 3.5)

dt3 = 0.25*H["L3"]
e3, s3, a3 = fit["L3"]
report(f"A3 finest-rung fidelity: |amp - 1| = {abs(a3-1):.2e} <= 0.002 and "
       f"|shift| = {abs(s3):.3f} <= dt_step/2 = {dt3/2:.3f}",
       abs(a3 - 1) <= 0.002 and abs(s3) <= dt3/2)

print(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
sys.exit(0 if OK else 1)
