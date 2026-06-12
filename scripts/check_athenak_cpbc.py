#!/usr/bin/env python3
"""check_athenak_cpbc.py — N13 gates for the constraint-preserving
Gamma-tilde boundary row (formulation eq:bc-cpbc; arXiv:1010.0523v2
eq:general_CPBCs, linearized realization).

  B1 Minkowski + CPBC preserved: final LINF deviation < 1e-12 (amp = 0).
  B2 constraint absorption: after the pulse reaches the boundary (t >= 62),
     the CPBC run's L2 Hamiltonian and momentum norms stay BELOW the
     Sommerfeld run's: time-integrated ratio < 1.0, with the measured value
     reported (the paper's CPBC tests show absorption, not reflection).
  B3 waveform unharmed: E vs the exact analytic waveform (N12 module,
     window t in [40, 58]) within 1.2x of the Sommerfeld run's.
  B4 convergence maintained under CPBC: order(164^3, 248^3) >= 3.5 vs exact.
  B5 full Z4c-CCM operator (CPBC + CCM datum): stable, constraints bounded
     (< 10x initial), waveform E within 1.2x of the CPBC run.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from teuk_exact_waveform import rpsi4_20_exact  # noqa: E402

base = Path(sys.argv[1])
X, RC, TAU, REXT = 1.0e-5, 20.0, 2.0, 36.0
OK = True


def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    OK &= bool(cond)


def hst(run):
    d = np.loadtxt(base/run/"teuk.user.hst")
    return d  # cols: time dt LINF-dev RMS-dev L2-H L2-M


def Ewf(run):
    d = np.loadtxt(base/run/"waveforms"/"rpsi4_real_0036.txt")
    t, y = d[:, 0], d[:, 3]
    win = (t >= 40.0) & (t <= 58.0)
    ex = rpsi4_20_exact(t[win], REXT, X, RC, TAU)
    return np.linalg.norm(y[win] - ex)/np.linalg.norm(ex)


# B1
fl = hst("flat")
report(f"B1 Minkowski + CPBC: final LINF-dev = {fl[-1, 2]:.2e} < 1e-12",
       fl[-1, 2] < 1e-12)

# B2
hs, hc = hst("somm"), hst("cpbc")
n = min(len(hs), len(hc))
hs, hc = hs[:n], hc[:n]
late = hs[:, 0] >= 62.0
ratH = np.trapz(hc[late, 4], hs[late, 0])/np.trapz(hs[late, 4], hs[late, 0])
ratM = np.trapz(hc[late, 5], hs[late, 0])/np.trapz(hs[late, 5], hs[late, 0])
report(f"B2 boundary constraint absorption (t >= 62): "
       f"int L2-H ratio cpbc/somm = {ratH:.3f} < 1.0, "
       f"int L2-M ratio = {ratM:.3f} < 1.0", ratH < 1.0 and ratM < 1.0)

# B3
Es, Ec = Ewf("somm"), Ewf("cpbc")
report(f"B3 waveform unharmed: E(cpbc) = {Ec:.3e} <= 1.2 x E(somm) = "
       f"{1.2*Es:.3e}", Ec <= 1.2*Es)

# B4
Ec2 = Ewf("cpbc2")
p = np.log(Ec/Ec2)/np.log(0.5/(82.0/248.0))
report(f"B4 convergence with CPBC: E(164) = {Ec:.3e}, E(248) = {Ec2:.3e}, "
       f"order p = {p:.2f} >= 3.5", p >= 3.5)

# B5
hcc = hst("ccmcp")
Ecc = Ewf("ccmcp")
report(f"B5 CPBC + CCM datum: constraints bounded "
       f"(L2-H {hcc[0, 4]:.2e} -> {hcc[-1, 4]:.2e} < 10x) and "
       f"E = {Ecc:.3e} <= 1.2 x E(cpbc)",
       hcc[-1, 4] < 10*hcc[0, 4] and Ecc <= 1.2*Ec)

print(f"\nOVERALL: {'PASS' if OK else 'FAIL'}")
sys.exit(0 if OK else 1)
