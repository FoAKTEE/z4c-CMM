#!/usr/bin/env python3
"""plot_paper_figures_v2.py — iter 48: method-paper figures from the
CORRECTED batteries (post pgen-frame + time-label fixes; ledger N12/N14).

Fig 1  zccm_teuk_exact.pdf   exact-reference Teukolsky ladder (N12b):
       gpu6 L1/L2/L3 (6th-order, nghost=4) vs the EXACT analytic
       rpsi4(2,0) at r=36 (the corrected program replaces the
       causally-disconnected reference run with the analytic truth).
       Data: results/numerical/athenak_teuk_ana/gpu6/L*/waveforms/,
       exact: scripts/teuk_exact_waveform.py (N12a coefficients).

Fig 2  zccm_live_ccm.pdf     live-CCM results (N14):
       (a) Sec V.C pulse: live psi0 datum (mode-6 diag) vs the ZccmJl
           mirror reference (transparent-vs-mirror separation = the CCM
           physics; data results/numerical/n14_pulse_ref.csv +
           n14_native/mode6/run.log);
       (b) the transmitted burst rpsi4(2,0) at r=36 (mode 6 vs the
           no-datum control);
       (c) SpECTRE cross-oracle: Psi4 and Psi0 (2,0) at scri vs the
           admitted closed forms -sqrt(6pi/5) F6 and -(9/8)c2 F2
           (results/numerical/spectre_oracle/TeukolskyOracleReduction.h5).
"""
import sys, math, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path("/data/haiyangw/claude/z4c-CMM")
sys.path.insert(0, str(ROOT/"scripts"))
from teuk_exact_waveform import rpsi4_20_exact, F_n

OUT = ROOT/"paper/z4c-CMM/zccm_formulation/figures"
X, RC, TAU = 1.0e-5, 20.0, 2.0

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "lines.linewidth": 1.0})


def wf(path):
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 3]                  # time, (2,0)


# ---------------- Fig 1: exact-reference ladder ----------------
fig, ax = plt.subplots(2, 1, figsize=(6.6, 5.2), sharex=True,
                       gridspec_kw={"height_ratios": [2.0, 1.3]})
cols = {"L1": "C0", "L2": "C1", "L3": "C2"}
hs = {"L1": 0.5, "L2": 0.3306, "L3": 0.2770}
Es = {"L1": 4.9864e-2, "L2": 5.3607e-3, "L3": 2.0189e-3}
tref = np.linspace(36.0, 80.0, 1200)
ax[0].plot(tref, [rpsi4_20_exact(t, 36.0, X, RC, TAU)/X for t in tref],
           "k-", lw=1.6, label="exact analytic (N12)")
for L in ("L1", "L2", "L3"):
    t, w = wf(ROOT/f"results/numerical/athenak_teuk_ana/gpu6/{L}/waveforms"
                   "/rpsi4_real_0036.txt")
    ax[0].plot(t, w/X, color=cols[L], ls="--", lw=0.9,
               label=f"{L}: h = {hs[L]:.3f}")
    ex = np.array([rpsi4_20_exact(tt, 36.0, X, RC, TAU) for tt in t])
    pk = np.max(np.abs(ex))
    ax[1].semilogy(t, np.abs(w - ex)/pk, color=cols[L], lw=0.9,
                   label=f"{L}: E = {Es[L]:.2e}")
ax[0].set_ylabel(r"$r\psi_4^{(2,0)}/X$ at $r=36$")
ax[0].legend(loc="upper right", ncol=2)
ax[0].set_xlim(36, 80)
ax[1].set_ylabel(r"$|{\rm run}-{\rm exact}|/{\rm peak}$")
ax[1].set_xlabel(r"$t$")
ax[1].legend(loc="upper right", title="order 5.39 (L1,L2), 5.52 (L2,L3)")
fig.tight_layout()
fig.savefig(OUT/"zccm_teuk_exact.pdf")
print("wrote", OUT/"zccm_teuk_exact.pdf")

# ---------------- Fig 2: live CCM ----------------
fig, ax = plt.subplots(3, 1, figsize=(6.6, 7.4))
# (a) pulse datum: live vs mirror
ref = np.loadtxt(ROOT/"results/numerical/n14_pulse_ref.csv",
                 delimiter=",", comments="#")
tsl, psl = [], []
for line in open(ROOT/"results/numerical/n14_native/mode6/run.log"):
    if line.startswith("ccm5-diag"):
        kv = dict(tok.split("=", 1) for tok in line.split()[1:])
        tsl.append(float(kv["t"])); psl.append(float(kv["psi0"]))
ax[0].plot(ref[:, 1], ref[:, 2]/1e-6, "k-", lw=1.2,
           label="ZccmJl mirror (CCE) reference")
ax[0].plot(tsl, np.array(psl)/1e-6, "C3o", ms=2.5,
           label="live datum (mode 6, CCM)")
for tmark, lab in ((52.5, r"$t_1$"), (161.0, r"$t_2=t_1+2R$")):
    ax[0].axvline(tmark, color="0.6", ls=":", lw=0.8)
    ax[0].text(tmark + 1, ax[0].get_ylim()[1]*0.8, lab, fontsize=8)
ax[0].set_ylabel(r"$\psi_0^t(r_B)\ [10^{-6}]$")
ax[0].set_xlabel(r"$t$")
ax[0].legend(loc="lower left")
ax[0].set_title("(a) Sec. V.C pulse: transparent (live CCM) vs mirror "
                "(quiescent-BC) worldtube", fontsize=8.5)
# (b) transmitted burst
t6, w6 = wf(ROOT/"results/numerical/n14_native/mode6/waveforms"
                 "/rpsi4_real_0036.txt")
tc, wc = wf(ROOT/"results/numerical/n14_native/mode6_ctrl/waveforms"
                 "/rpsi4_real_0036.txt")
ax[1].plot(t6, w6/1e-4, "C0-", lw=1.0, label="mode 6 (pulse injected)")
ax[1].plot(tc, wc/1e-4, "0.5", lw=0.8, ls="--",
           label=r"control (no datum): $2.9\times10^{-10}$ of peak")
ax[1].axvline(161.0, color="0.6", ls=":", lw=0.8)
ax[1].set_ylabel(r"$r\psi_4^{(2,0)}\ [10^{-4}]$ at $r=36$")
ax[1].set_xlabel(r"$t$")
ax[1].legend(loc="upper left")
ax[1].set_title("(b) the burst is transmitted through the Cauchy domain "
                "and exits at $t_2$", fontsize=8.5)
# (c) SpECTRE oracle
import h5py
f = h5py.File(ROOT/"results/numerical/spectre_oracle"
                   "/TeukolskyOracleReduction.h5", "r")
XO = 1.0e-3
c2 = 4*math.sqrt(2*math.pi/15)
for nm, n, A, col, lab in (
        ("Psi4", 6, -math.sqrt(6*math.pi/5), "C0",
         r"$\Psi_4$: $-\sqrt{6\pi/5}\,F^{(6)}$"),
        ("Psi0", 2, -9/8*c2, "C3",
         r"$\Psi_0$: $-(9/8)\,c_2\,F^{(2)}$")):
    d = np.asarray(f[f"SpectreR0041.cce/{nm}"])
    t, v = d[:, 0], d[:, 13]
    ax[2].plot(t, v/XO, color=col, lw=2.2, alpha=0.35,
               label=f"SpECTRE {nm} (2,0)")
    ax[2].plot(t, [A*F_n(n, tt, XO, RC, TAU)/XO for tt in t], color=col,
               lw=0.9, ls="--", label=lab)
ax[2].set_xlim(5, 40)
ax[2].set_ylabel(r"scri output $/Z$")
ax[2].set_xlabel(r"$u$")
ax[2].legend(loc="upper right", fontsize=7)
ax[2].set_title("(c) SpECTRE CharacteristicExtract on the exact worldtube "
                "data vs the admitted closed forms", fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT/"zccm_live_ccm.pdf")
print("wrote", OUT/"zccm_live_ccm.pdf")
