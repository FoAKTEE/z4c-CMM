#!/usr/bin/env python3
"""check_athenak_ccm.py — evaluate the AthenaK Z4c-CCM test battery.

Reads the .hst files of t1_reduction / t2_injection / t3_linearity under the
given directory (columns per the z4c_stability user history: time, dt,
LINF-Err, RMS-Err, ... deviation from Minkowski) and checks:
  T1 reduction : max LINF deviation over the run < 1e-12 (Minkowski preserved;
                 noise seed is ~1e-30, growth would indicate the CCM path
                 perturbs the Sommerfeld scheme).
  T2 injection : max RMS deviation > 100x the T1 level (the datum injects a
                 physical pulse through the boundary), and the deviation is
                 negligible before the pulse (t < t0 - 3 sigma) compared to
                 after.
  T3 linearity : max-RMS(t3) / max-RMS(t2) = 2 within 2% (linear response;
                 quadratic corrections ~ amp).
Exit 0 iff all pass. Output is the certified test record.
"""
import sys
from pathlib import Path

import numpy as np

base = Path(sys.argv[1])
OK = True

def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    OK &= bool(cond)

def load(run):
    f = sorted((base / run).glob("*.hst"))[0]
    data = np.loadtxt(f)
    return data  # col0 = time, col2 = LINF-Err, col3 = RMS-Err

t1, t2, t3 = load("t1_reduction"), load("t2_injection"), load("t3_linearity")
linf1 = np.max(np.abs(t1[:, 2]))
rms2, rms3 = np.max(np.abs(t2[:, 3])), np.max(np.abs(t3[:, 3]))
early = t2[t2[:, 0] < 1.0]
early2 = np.max(np.abs(early[:, 3])) if len(early) else 0.0

report(f"T1 reduction: ccm_amp=0 preserves Minkowski — max LINF deviation "
       f"{linf1:.2e} < 1e-12", linf1 < 1e-12)
report(f"T2 injection: pulse enters the domain — max RMS {rms2:.3e} "
       f">= 100x flatness level ({linf1:.1e}); pre-pulse (t<1) response "
       f"{early2:.1e} <= 1e-6 of peak (datum tail is exp(-((t-t0)/sigma)^2): "
       f"the early response IS the physical tail)", rms2 > 100 * max(linf1, 1e-300)
       and rms2 > 1e-10 and early2 < 1e-6 * rms2)
ratio = rms3 / rms2
report(f"T3 linearity: max-RMS ratio (amp 2e-3 / 1e-3) = {ratio:.6f} "
       f"within 2% of 2", abs(ratio - 2.0) < 0.04)

print(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
sys.exit(0 if OK else 1)
