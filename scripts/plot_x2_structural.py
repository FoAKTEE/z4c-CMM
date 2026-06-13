#!/usr/bin/env python3
"""plot_x2_structural.py — the X=2 figure from OUR full campaign
(results/numerical/ccm_repro_x2_full): a causally disconnected
far-boundary reference (768^3, box +-192) and boundary-41
Sommerfeld/CCM pairs on a resolution ladder, all at the stable X=2
configuration (diss 0.1, kappa1 0.1, nghost 4, CFL 0.25, t=80).
Overwrites figures/zccm_x2_structural.pdf with current data."""
import pathlib
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path("/data/haiyangw/claude/z4c-CMM")
BASE = ROOT/"results/numerical/ccm_repro_x2_full"
OUT = ROOT/"paper/z4c-CMM/zccm_formulation/figures/zccm_x2_structural.pdf"

def wf(d):
    a = np.loadtxt(BASE/d/"waveforms/rpsi4_real_0036.txt", comments="#")
    return a[:, 0], a[:, 3]

# common grid (the ladder rungs carry slightly different dt time labels)
tg = np.arange(1.0, 79.0, 0.25)
def gi(d):
    t, w = wf(d)
    return interp1d(t, w, bounds_error=False, fill_value=0.0)(tg)

# runs available now (h=0.125 still computing -> use the complete rungs)
runs = {d: gi(d) for d in ("reference_768", "somm_h050", "ccm_h050",
                           "somm_h025", "ccm_h025")}
have0125 = (BASE/"somm_h0125/waveforms/rpsi4_real_0036.txt").exists() and \
           (BASE/"ccm_h0125/waveforms/rpsi4_real_0036.txt").exists()
R = runs["reference_768"]; nR = np.linalg.norm(R); pk = np.max(np.abs(R))

eb = np.linalg.norm(runs["somm_h050"] - R)/nR          # boundary error h=0.5
match5 = np.linalg.norm(runs["ccm_h050"] - runs["somm_h050"])/nR
Eh = np.linalg.norm(runs["somm_h050"] - runs["somm_h025"])/nR  # truncation

plt.rcParams.update({"font.size": 9, "legend.fontsize": 7.6,
                     "lines.linewidth": 1.0})
fig, ax = plt.subplots(2, 1, figsize=(7, 7.4), sharex=True,
                       constrained_layout=True)

# Top: waveforms (normalized by reference peak)
ax[0].plot(tg, R/pk, "k-", lw=1.6,
           label=r"far-boundary reference ($768^3$, box $\pm192$)")
ax[0].plot(tg, runs["somm_h050"]/pk, "r--", lw=1.0,
           label=r"Sommerfeld, boundary 41, $h{=}0.5$")
ax[0].plot(tg, runs["ccm_h050"]/pk, "b:", lw=1.3,
           label=r"Z4c-CCM, boundary 41, $h{=}0.5$")
ax[0].plot(tg, runs["somm_h025"]/pk, color="0.6", lw=0.8,
           label=r"Sommerfeld, $h{=}0.25$")
ax[0].set_ylabel(r"$r\psi_4^{(2,0)}/\max_{\rm ref}$ at $r=36$")
ax[0].set_title(r"Teukolsky $X=2$: boundary-41 runs vs.\ the "
                r"far-boundary reference", fontsize=10)
ax[0].legend(loc="upper left")

# Bottom: the error hierarchy
ax[1].semilogy(tg, np.abs(runs["somm_h050"] - R)/pk, "r-", lw=1.1,
               label=rf"$|{{\rm Somm}}_{{0.5}}-{{\rm ref}}|$: boundary error "
                     rf"$E={eb:.2e}$")
ax[1].semilogy(tg, np.abs(runs["ccm_h050"] - R)/pk, "b--", lw=1.1,
               label=r"$|{\rm CCM}_{0.5}-{\rm ref}|$ (coincident: peeling)")
ax[1].semilogy(tg, np.abs(runs["somm_h050"] - runs["somm_h025"])/pk,
               color="0.5", lw=1.0,
               label=rf"$|h_{{0.5}}-h_{{0.25}}|$: truncation $E_h={Eh:.2f}$")
ax[1].semilogy(tg, np.abs(runs["ccm_h050"] - runs["somm_h050"])/pk, "g-",
               lw=1.1,
               label=rf"$|{{\rm CCM}}-{{\rm Somm}}|$ (matching) $={match5:.1e}$")
ax[1].set_xlim(0, 80)
ax[1].set_ylim(1e-12, 5)
ax[1].set_xlabel("Time")
ax[1].set_ylabel(r"$|\Delta\,r\psi_4^{(2,0)}|/\max_{\rm ref}$")
ax[1].legend(loc="lower right", ncol=1)
fig.savefig(OUT, bbox_inches="tight")
print(f"wrote {OUT}")
print(f"boundary error E[Somm]=E[CCM]={eb:.4e}; matching |CCM-Somm|={match5:.3e}; "
      f"truncation E_h={Eh:.4e}; h=0.125 complete: {have0125}")
