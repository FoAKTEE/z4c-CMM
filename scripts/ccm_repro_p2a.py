#!/usr/bin/env python3
"""ccm_repro_p2a.py — mission 4 P2a: apples-to-apples X=1e-5 Teukolsky
strain overlay at scri: the authors' published CCM/CCE curves
(ref-code/ccm-figures/X=1e-5) vs OUR SpECTRE-CharacteristicExtract run on
the exact analytic worldtube data (results/numerical/spectre_oracle,
X=1e-3, worldtube 41 = THEIR worldtube radius), rescaled linearly to
X=1e-5, with the measured (not assumed) time shift.

Shared truth: h_pert(t) = sqrt(6pi/5) F''''(t - c) (their notebook uses
c = 21 = r_out - r_c, the scri arrival of the outgoing pulse; our oracle
file was built with tube phase u = t, so ours peaks at u = 20 -> expect
a +1 shift, measured below)."""
import math, pathlib, sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

ROOT = pathlib.Path("/data/haiyangw/claude/z4c-CMM")
THEIRS = ROOT/"ref-code/ccm-figures/X=1e-5"
X_THEIRS, X_OURS = 1.0e-5, 1.0e-5
RC, TAU = 20.0, 2.0

def F4(t, X, rc, tau):
    x = np.asarray(t) - rc
    e = X*np.exp(-x**2/tau**2)
    return e*(16*x**4 - 48.0*x**2*tau**2 + 12*tau**4)/tau**8

# ---- their curves ----
ccm = np.loadtxt(THEIRS/"299ccm/Strain.dat")
cce = np.loadtxt(THEIRS/"299cce/Strain.dat")

# ---- our oracle strain (2,0) ----
f = h5py.File(ROOT/"results/numerical/spectre_oracle_x1e5/TeukolskyOracleReductionTight.h5", "r")
d = np.asarray(f["SpectreR0041.cce/Strain"])
t_o, h_o = d[:, 0], d[:, 13]*(X_THEIRS/X_OURS)        # linear rescale

# ---- measured time shift: align our oracle to THEIR time axis by
# maximizing overlap with the shared perturbative envelope ----
pert = lambda t, c: math.sqrt(6*math.pi/5)*F4(t, X_THEIRS, RC + (c - RC), TAU)
# fit shift of OUR curve against sqrt(6pi/5) F4(t - c): scan c
tgrid = np.arange(5.0, 40.0, 0.05)
ho_i = interp1d(t_o, h_o, bounds_error=False, fill_value=0.0)(tgrid)
best = (1e30, 0.0)
for c in np.arange(19.0, 23.0, 0.0005):
    model = math.sqrt(6*math.pi/5)*F4(tgrid - (c - RC), X_THEIRS, RC, TAU)
    r = np.max(np.abs(ho_i - model))
    if r < best[0]:
        best = (r, c)
res_o, c_ours = best
best = (1e30, 0.0)
ccm_i = interp1d(ccm[:, 0], ccm[:, 1], bounds_error=False, fill_value=0.0)(tgrid)
for c in np.arange(19.0, 23.0, 0.0005):
    model = math.sqrt(6*math.pi/5)*F4(tgrid - (c - RC), X_THEIRS, RC, TAU)
    r = np.max(np.abs(ccm_i - model))
    if r < best[0]:
        best = (r, c)
res_t, c_theirs = best
shift = c_theirs - c_ours
pk = np.max(np.abs(ccm_i))
print(f"measured pulse-center: ours c = {c_ours:.3f}, theirs c = {c_theirs:.3f}, shift = {shift:.3f}")
print(f"shape residual vs perturbative (of their peak {pk:.3e}):")
print(f"  their CCM : {res_t/pk:.3e}")
print(f"  our oracle: {res_o/pk:.3e}")

# ---- direct scri-to-scri overlay after the measured shift ----
ho_shift = interp1d(t_o + shift, h_o, bounds_error=False, fill_value=0.0)(tgrid)
d_ccm = np.max(np.abs(ho_shift - ccm_i))
cce_i = interp1d(cce[:, 0], cce[:, 1], bounds_error=False, fill_value=0.0)(tgrid)
d_cce = np.max(np.abs(ho_shift - cce_i))
print(f"direct overlay max|ours - their CCM| = {d_ccm:.3e} = {d_ccm/pk:.3e} of peak")
print(f"direct overlay max|ours - their CCE| = {d_cce:.3e} = {d_cce/pk:.3e} of peak")

# ---- figure ----
fig, ax = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                       constrained_layout=True)
model = math.sqrt(6*math.pi/5)*F4(tgrid - (c_theirs - RC), X_THEIRS, RC, TAU)
ax[0].semilogy(tgrid, np.abs(model), "b", lw=1.5, label="Perturbative (shared truth)")
ax[0].semilogy(tgrid, np.abs(cce_i), "k", lw=1.5, label="their CCE (299cce)")
ax[0].semilogy(tgrid, np.abs(ccm_i), "r", lw=1.5, label="their CCM (299ccm)")
ax[0].semilogy(tgrid, np.abs(ho_shift), "g--", lw=1.2,
               label=f"OUR SpECTRE-CCE on exact worldtube data\n(X = 1e-5 matched, tight tolerances; shift {shift:+.2f})")
ax[0].set_ylim(1e-13, 5e-5)
ax[0].set_ylabel(r"$|h_{2,0}|$")
ax[0].set_title(r"Teukolsky wave $X=10^{-5}$, worldtube 41: scri strain, theirs vs ours")
ax[0].legend(fontsize=8)
ax[1].semilogy(tgrid, np.abs(ccm_i - model), "r", lw=1.2,
               label=f"|their CCM - pert| (max {res_t/pk:.2e} of peak)")
ax[1].semilogy(tgrid, np.abs(ho_shift - model), "g--", lw=1.2,
               label=f"|ours - pert| (max {res_o/pk:.2e})")
ax[1].semilogy(tgrid, np.abs(ho_shift - ccm_i), "k", lw=1.2,
               label=(f"|ours - their CCM| (max {d_ccm/pk:.2e})\n"
                      "NOTE: ours = CCE transport of EXACT worldtube data,"
                      "\nso this curve is dominated by THEIR deviation"
                      "\n(triangle bound: differs from |their CCM - pert|"
                      f"\nby at most |ours - pert| = {res_o/pk:.1e} pointwise)"))
ax[1].set_ylim(1e-13, 3e-7)
ax[1].set_xlabel("Time (their axis)")
ax[1].set_ylabel("Differences")
ax[1].legend(fontsize=8)
out = ROOT/"results/figures/ccm_paper_repro_p2/x1e-5"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out/"strain_overlay.pdf", bbox_inches="tight")
print("wrote", out/"strain_overlay.pdf")
