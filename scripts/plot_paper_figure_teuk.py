#!/usr/bin/env python3
"""plot_paper_figure_teuk.py — publication figure for the method paper:
exact-parameter Teukolsky test (arXiv:2308.10361 Sec. V.A protocol), two
columns (X = 1e-5, X = 2), waveform + log-difference rows, vector PDF.
All curves from results/numerical/athenak_teuk_paper/ (commit 0eac498 runs).

Usage: plot_paper_figure_teuk.py <out.pdf>
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
    "legend.fontsize": 7.2, "mathtext.fontset": "cm", "font.family": "serif",
    "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT/"results/numerical/athenak_teuk_paper"
RC, TAU, REXT = 20.0, 2.0, 36.0

def F6(u, X):
    s = (u - RC)/TAU
    H = ((64*s**2 - 480)*s**2 + 720)*s**2 - 120
    return X*(-1.0/TAU)**6*H*np.exp(-s**2)

def wf(tag, run):
    d = np.loadtxt(BASE/tag/run/"waveforms"/"rpsi4_real_0036.txt")
    return d[:, 0], d[:, 3]

fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True,
                         gridspec_kw={"height_ratios": [1, 1.15],
                                      "hspace": 0.08, "wspace": 0.26})

for col, (tag, X, title) in enumerate(
        ((("paper_X1e-5"), 1e-5, r"$X = 10^{-5}$"),
         (("paper_X2"), 2.0, r"$X = 2$"))):
    t_ref, ref_raw = wf(tag, "ref")
    y = {}
    runs = {n: wf(tag, n) for n in ("somm", "somm_mid", "ccm", "ccm_mid")}
    tmax = min([t_ref[-1]] + [runs[n][0][-1] for n in runs])
    t = t_ref[t_ref <= tmax]
    ref = np.interp(t, t_ref, ref_raw)
    y = {n: np.interp(t, *runs[n]) for n in runs}

    axw, axd = axes[0][col], axes[1][col]
    axw.plot(t, ref/X, "k-", lw=1.3, label="reference (boundary $104$)")
    axw.plot(t, y["somm"]/X, "--", color="#1f77b4", lw=1.0, label="Sommerfeld")
    axw.plot(t, y["ccm"]/X, ":", color="#d62728", lw=1.6,
             label=r"Z4c-CCM, datum $\psi_0^{\rm Teuk}$")
    if X <= 1e-4:
        axw.plot(t, -np.sqrt(6*np.pi/5)*F6(t - REXT, X)/X, "-", color="0.65",
                 lw=0.7, label=r"$-\sqrt{6\pi/5}\,F^{(6)}(t-r)$")
    axw.set_title(title)
    if col == 0:
        axw.set_ylabel(r"$r\,\psi_4^{(2,0)}/X$")
        axw.legend(loc="upper left", frameon=False, handlelength=1.7)

    eps = 1e-300
    axd.semilogy(t, np.abs(y["ccm"] - ref)/X + eps, "-", color="#d62728",
                 lw=1.1, label=r"$|{\rm CCM}-{\rm ref}|$")
    axd.semilogy(t, np.abs(y["somm"] - ref)/X + eps, "-", color="#1f77b4",
                 lw=0.8, label=r"$|{\rm Somm.}-{\rm ref}|$")
    axd.semilogy(t, np.abs(y["ccm"] - y["somm"])/X + eps, "-",
                 color="#2ca02c", lw=0.8, label=r"$|{\rm CCM}-{\rm Somm.}|$")
    axd.semilogy(t, np.abs(y["ccm"] - y["ccm_mid"])/X + eps, "--",
                 color="#d62728", lw=0.8,
                 label=r"CCM err. est. $(164^3{-}132^3)$")
    axd.semilogy(t, np.abs(y["somm"] - y["somm_mid"])/X + eps, "--",
                 color="#1f77b4", lw=0.8,
                 label=r"Somm. err. est. $(164^3{-}132^3)$")
    pk = np.max(np.abs(ref))/X
    axd.set_ylim(pk*1e-10, pk*30)
    axd.set_xlabel(r"$t$")
    if col == 0:
        axd.set_ylabel(r"$|\Delta\, r\psi_4^{(2,0)}|/X$")
        axd.legend(loc="upper left", frameon=False, ncol=2, handlelength=1.7,
                   columnspacing=1.0)
    axd.set_xlim(0, 80)

out = sys.argv[1]
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
