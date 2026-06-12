#!/usr/bin/env python3
"""plot_athenak_teuk.py — figure in the style of the CCM 2024 paper's
flat-background Teukolsky test (arXiv:2308.10361 Sec. V.A): waveform panel,
deviation-from-reference panel (log), and constraint panel — all from THIS
repository's AthenaK runs.

Usage: plot_athenak_teuk.py <run-dir> <amp> <out.png> [<run-dir2> <amp2>]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RC, TAU, REXT = 4.0, 1.0, 8.0

def F_n(n, u, amp):
    s = (u - RC)/TAU
    H = [1, 2*s, 4*s**2 - 2, 8*s**3 - 12*s, 16*s**4 - 48*s**2 + 12,
         32*s**5 - 160*s**3 + 120*s,
         64*s**6 - 480*s**4 + 720*s**2 - 120][n]
    return amp*(-1.0/TAU)**n*H*np.exp(-s**2)

def wf(base, run):
    d = np.loadtxt(Path(base)/run/"waveforms"/"rpsi4_real_0008.txt")
    return d[:, 0], d[:, 3]

def hst(base, run):
    import glob
    d = np.loadtxt(sorted(Path(base, run).glob("*.hst"))[0])
    return d[:, 0], d[:, 4]

def panel_set(axw, axd, axc, base, amp, title):
    t, s = wf(base, "somm"); _, c = wf(base, "ccm"); _, r = wf(base, "ref")
    n = min(len(s), len(c), len(r)); t, s, c, r = t[:n], s[:n], c[:n], r[:n]
    ana = -np.sqrt(6*np.pi/5)*F_n(6, t - REXT, amp)
    axw.plot(t, r, "k-", lw=1.8, label="reference (double domain)")
    axw.plot(t, s, "C0--", lw=1.4, label="Sommerfeld")
    axw.plot(t, c, "C3:", lw=1.8, label="Z4c-CCM (Teukolsky self-datum)")
    axw.plot(t, ana*np.max(np.abs(r))/np.max(np.abs(ana)), color="0.7", lw=0.9,
             label=r"scri formula $\propto F^{(6)}(t-r)$ (scaled)")
    axw.set_ylabel(r"$r\,\psi_4^{(2,0)}$")
    axw.set_title(title)
    axw.legend(fontsize=7, loc="upper left")
    eps = 1e-300
    axd.semilogy(t, np.abs(s - r) + eps, "C0--", lw=1.4, label="|Sommerfeld $-$ ref|")
    axd.semilogy(t, np.abs(c - r) + eps, "C3:", lw=1.8, label="|CCM $-$ ref|")
    axd.semilogy(t, np.abs(c - s) + eps, "C2-", lw=1.0, label="|CCM $-$ Sommerfeld|")
    axd.set_ylabel(r"$|\Delta\, r\psi_4^{(2,0)}|$")
    axd.set_ylim(np.max(np.abs(r))*1e-9, np.max(np.abs(r)))
    axd.legend(fontsize=7, loc="upper left")
    for run, sty, lab in (("ref", "k-", "reference"), ("somm", "C0--", "Sommerfeld"),
                          ("ccm", "C3:", "CCM")):
        th, hh = hst(base, run)
        axc.semilogy(th, hh, sty, lw=1.4, label=lab)
    axc.set_ylabel(r"$\|H\|_{L_2}$")
    axc.set_xlabel(r"$t\ [M]$")
    axc.legend(fontsize=7)

args = sys.argv[1:]
sets = [(args[0], float(args[1]))] + ([(args[3], float(args[4]))] if len(args) > 4 else [])
out = args[2]
ncol = len(sets)
fig, axes = plt.subplots(3, ncol, figsize=(5.2*ncol, 7.5), sharex="col",
                         squeeze=False)
titles = [rf"$X = {amp:g}$" + (" (perturbative)" if amp <= 1e-4 else " (nonlinear)")
          for _, amp in sets]
for cidx, ((base, amp), ti) in enumerate(zip(sets, titles)):
    panel_set(axes[0][cidx], axes[1][cidx], axes[2][cidx], base, amp, ti)
fig.suptitle("AthenaK Z4c-CCM: flat-background Teukolsky-wave test "
             "(protocol of arXiv:2308.10361 Sec. V.A)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(out, dpi=180)
print("wrote", out)
