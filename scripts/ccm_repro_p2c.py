#!/usr/bin/env python3
"""ccm_repro_p2c.py — mission 4 P2c: Sec V.C pulse-injection overlays vs
the authors' published pulse_on_characteristic_grid data.

psi0 panel: their 377ccm/InnerPsi0.dat (modal (2,0) at their boundary 41)
vs OUR native-solver curves converted by the DERIVED s^2->modal factor
c2 = 4 sqrt(2pi/15): the ZccmJl mirror reference and the mode-6 live
series. psi4 panel: their scri Psi4 (ccm + cce) vs our r=36 rpsi4(2,0)
scaled by a FITTED constant (reported on the figure; the derived modal
factor c2 accounts for most of it)."""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path("/data/haiyangw/claude/z4c-CMM")
TH = ROOT/"ref-code/ccm-figures/pulse_on_characteristic_grid"
C2 = 4*np.sqrt(2*np.pi/15)

ip0 = np.loadtxt(TH/"377ccm/InnerPsi0.dat")
p4c = np.loadtxt(TH/"377ccm/Psi4.dat")
p4e = np.loadtxt(TH/"cce/Psi4.dat")
ref = np.loadtxt(ROOT/"results/numerical/n14_pulse_ref.csv", delimiter=",",
                 comments="#")
ts, ps = [], []
for line in open(ROOT/"results/numerical/n14_native/mode6/run.log"):
    if line.startswith("ccm5-diag"):
        kv = dict(tok.split("=", 1) for tok in line.split()[1:])
        ts.append(float(kv["t"])); ps.append(float(kv["psi0"]))
ts, ps = np.array(ts), np.array(ps)
our4 = np.loadtxt(ROOT/"results/numerical/n14_native/mode6/waveforms"
                       "/rpsi4_real_0036.txt", comments="#")

# measured numbers
i_t = np.argmin(ip0[:, 1])
mir_t = ref[np.argmin(ref[:, 2]), 1]
mir_v = np.min(ref[:, 2])*C2
print(f"psi0 trough: theirs {ip0[i_t,1]:.4e} at t={ip0[i_t,0]:.2f} | "
      f"our mirror x c2 {mir_v:.4e} at t={mir_t:.2f} | "
      f"amp ratio {mir_v/ip0[i_t,1]:.4f} | dt = {ip0[i_t,0]-mir_t:+.2f} "
      f"(geometric expectation +1.3..+3.0: tube 41 vs 40)")
scale = np.min(our4[:, 3])/np.min(p4c[:, 1])
print(f"psi4 fitted scale ours/theirs = {scale:.4f} (c2 = {C2:.4f}; "
      f"residual factor {scale/C2:.4f} = finite-r 36 vs scri + conv)")

fig, ax = plt.subplots(2, 1, figsize=(7.2, 8.4), constrained_layout=True)
ax[0].plot(ip0[:, 0], ip0[:, 1]/1e-6, "k-", lw=1.4,
           label="THEIRS: CCM InnerPsi0 (boundary 41, full CCE matching)")
ax[0].plot(ref[:, 1], ref[:, 2]*C2/1e-6, "b--", lw=1.2,
           label=r"OURS mirror solver $\times\,c_2$ (quiescent tube; CCE analog)")
ax[0].plot(ts, ps*C2/1e-6, "r.", ms=3,
           label=r"OURS live mode-6 $\times\,c_2$ (two-way coupled)")
ax[0].set_xlim(0, 200)
ax[0].set_xlabel("Time (each code's own axis)")
ax[0].set_ylabel(r"$\psi_0^{(2,0)}$ at the boundary $[10^{-6}]$")
ax[0].legend(fontsize=8, loc="upper right")
ax[0].set_title("Sec. V.C pulse: boundary $\\psi_0$ — theirs vs ours "
                "(derived modal factor, no fit)", fontsize=10)
ax[1].plot(p4c[:, 0], p4c[:, 1]/1e-4, "r-", lw=1.4, label="THEIRS: CCM Psi4 (scri)")
ax[1].plot(p4e[:, 0], p4e[:, 1]/1e-4, "k-", lw=1.0, label="THEIRS: CCE Psi4 (mirror)")
ax[1].plot(our4[:, 0], our4[:, 3]/scale/1e-4, "b--", lw=1.2,
           label=f"OURS: rpsi4(2,0) at r=36 / {scale:.2f} (fitted; c2 within 5%)")
ax[1].set_xlim(0, 250)
ax[1].set_xlabel("Time (each code's own axis)")
ax[1].set_ylabel(r"$\psi_4^{(2,0)}\ [10^{-4}]$")
ax[1].legend(fontsize=8, loc="upper left")
ax[1].set_title("transmitted burst: their scri waveform vs our finite-radius "
                "extraction (feature offsets +2..+4 = tube 41 vs 40 geometry)",
                fontsize=10)
out = ROOT/"results/figures/ccm_paper_repro_p2/pulse"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out/"pulse_overlay.pdf", bbox_inches="tight")
print("wrote", out/"pulse_overlay.pdf")
