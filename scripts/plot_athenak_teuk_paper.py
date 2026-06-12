#!/usr/bin/env python3
"""plot_athenak_teuk_paper.py — EXACT-parameter Teukolsky figures in the
panel structure of arXiv:2308.10361 (Sec. V.A), all curves from THIS
repository's AthenaK runs (no paper graphics reused):

  waveform figure (analog of fig:low_amp_strain_only / large_amp_strain_only)
    top: (l=2, m=0) r psi4 at r = 36 — reference, Sommerfeld, CCM, and the
         analytic scri formula -sqrt(6 pi/5) F^(6)(t - r);
    bottom (log): |CCM - ref|, |Somm - ref|, |CCM - Somm|, |CCM - analytic|,
         and the two-resolution numerical-error estimates (dashed), the
         paper's "error est." curves.

  constraint figure (analog of fig:low_amp_con / gauge_constraint_X2)
    two panels: L2 Hamiltonian (top) and momentum (bottom) constraints vs t
    for the CCM run at three resolutions (104^3, 132^3, 164^3).

Usage: plot_athenak_teuk_paper.py <battery-dir> <amp> <wave.png> <con.png>
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = Path(sys.argv[1]); AMP = float(sys.argv[2])
out_wave, out_con = sys.argv[3], sys.argv[4]
RC, TAU, REXT = 20.0, 2.0, 36.0

def F_n(n, u):
    s = (u - RC)/TAU
    H = [1, 2*s, 4*s**2 - 2, 8*s**3 - 12*s, 16*s**4 - 48*s**2 + 12,
         32*s**5 - 160*s**3 + 120*s,
         64*s**6 - 480*s**4 + 720*s**2 - 120][n]
    return AMP*(-1.0/TAU)**n*H*np.exp(-s**2)

def wf(run):
    d = np.loadtxt(base/run/"waveforms"/f"rpsi4_real_{int(REXT):04d}.txt")
    return d[:, 0], d[:, 3]

t_ref, ref_raw = wf("ref")
series = {n: wf(n) for n in ("somm", "somm_mid", "ccm", "ccm_mid", "ccm_low")}
tmax = min([t_ref[-1]] + [series[n][0][-1] for n in series])
t = t_ref[t_ref <= tmax]
ref = np.interp(t, t_ref, ref_raw)
y = {n: np.interp(t, *series[n]) for n in series}
ana = -np.sqrt(6*np.pi/5)*F_n(6, t - REXT)

# ---- waveform figure ----
fig, (axw, axd) = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True,
                               gridspec_kw={"height_ratios": [1.1, 1]})
axw.plot(t, ref, "k-", lw=1.8, label=r"reference (boundary at $104$)")
axw.plot(t, y["somm"], "C0--", lw=1.4, label="Sommerfeld")
axw.plot(t, y["ccm"], "C3:", lw=2.0, label="Z4c-CCM (Teukolsky datum)")
if AMP <= 1e-4:   # the linear scri formula has no meaning at X = 2
    axw.plot(t, ana, color="0.6", lw=0.9,
             label=r"$-\sqrt{6\pi/5}\,F^{(6)}(t-r)$ (scri formula)")
axw.set_ylabel(r"$r\,\psi_4^{(2,0)}\ (r=36)$")
axw.set_title(rf"Teukolsky wave, $X = {AMP:g}$:  $r_c=20$, $\tau=2$, "
              rf"boundary $41$, $h=1/2$")
axw.legend(fontsize=8, loc="upper left")
eps = 1e-300
axd.semilogy(t, np.abs(y["ccm"] - ref) + eps, "C3-", lw=1.6,
             label="|CCM $-$ reference|")
axd.semilogy(t, np.abs(y["somm"] - ref) + eps, "C0-", lw=1.2,
             label="|Sommerfeld $-$ reference|")
axd.semilogy(t, np.abs(y["ccm"] - y["somm"]) + eps, "C2-", lw=1.0,
             label="|CCM $-$ Sommerfeld|")
if AMP <= 1e-4:
    axd.semilogy(t, np.abs(y["ccm"] - ana) + eps, color="0.6", lw=0.9,
                 label="|CCM $-$ scri formula|")
axd.semilogy(t, np.abs(y["ccm"] - y["ccm_mid"]) + eps, "C3--", lw=1.0,
             label=r"CCM error est. ($164^3 - 132^3$)")
axd.semilogy(t, np.abs(y["somm"] - y["somm_mid"]) + eps, "C0--", lw=1.0,
             label=r"Somm. error est. ($164^3 - 132^3$)")
pk = np.max(np.abs(ref))
axd.set_ylim(pk*1e-10, pk*3)
axd.set_ylabel(r"$|\Delta\, r\psi_4^{(2,0)}|$")
axd.set_xlabel(r"$t\ [M]$")
axd.legend(fontsize=7, loc="lower right", ncol=2)
fig.tight_layout()
fig.savefig(out_wave, dpi=180)
print("wrote", out_wave)

# ---- constraint figure ----
def hst(run):
    d = np.loadtxt(base/run/"teuk.user.hst")
    return d[:, 0], d[:, 4], d[:, 5]   # t, L2-H, L2-M

fig2, (axh, axm) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
for run, sty, lab in (("ccm_low", "C0-", r"CCM $104^3\ (h=0.788)$"),
                      ("ccm_mid", "C2-", r"CCM $132^3\ (h=0.621)$"),
                      ("ccm", "C3-", r"CCM $164^3\ (h=0.5)$"),
                      ("ref", "k--", r"reference $416^3\ (h=0.5)$")):
    th, hh, mm = hst(run)
    axh.semilogy(th, hh, sty, lw=1.4, label=lab)
    axm.semilogy(th, mm, sty, lw=1.4, label=lab)
axh.set_ylabel(r"$\|H\|_{L_2}$")
axh.set_title(rf"Constraint violation, $X = {AMP:g}$ (CCM run, three resolutions;"
              "\nnorms volume-averaged over each run's own domain)", fontsize=10)
axh.legend(fontsize=8)
axm.set_ylabel(r"$\|M\|_{L_2}$")
axm.set_xlabel(r"$t\ [M]$")
fig2.tight_layout()
fig2.savefig(out_con, dpi=180)
print("wrote", out_con)
