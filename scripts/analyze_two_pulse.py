#!/usr/bin/env python3
"""analyze_two_pulse.py — the S4 GATE analysis for the arXiv:2308.10361 X=2
Teukolsky two-pulse feature.

The terminal mission gate: with XCTS-solved constraint-satisfying initial data
on a large Cauchy domain (boundary r >= 100), the X=2 rPsi4 waveform extracted
at a finite radius must show the paper's TWO-PULSE structure —
  1. the PRIMARY physical pulse, and
  2. for a NON-transparent outer boundary (Sommerfeld/CCE), a SECONDARY spurious
     pulse from the boundary reflection,
whereas the Z4c-CCM transparent worldtube SUPPRESSES that secondary reflection.
Reproducing the pulse structure + arrival times (and the CCM-vs-Sommerfeld
contrast) proves the CCM implementation is bug-free.

This reads the AthenaK rPsi4 waveform of the S3 CCM run (and, if given, a paired
Sommerfeld run), finds the pulse peaks, reports their arrival times, and compares
to the light-travel-time predictions. Arrival-time model (paper Sec V.A.2,
generalized to an outer boundary / worldtube radius r_b and extraction radius
r_ext, pulse center r_c):
  - primary (outgoing pulse passing r_ext):            t_1 ~ r_ext - r_c
  - secondary (reflection off the boundary, back in):  t_2 ~ 2*r_b - r_ext - r_c  (+ 2(r_b-r_ext) after t_1)
  - tertiary (second boundary bounce):                 t_3 ~ t_2 + 2*r_b
These are the geodesic estimates; the script reports measured-vs-predicted and
the CCM/Sommerfeld secondary-pulse amplitude ratio (the transparency metric).

Usage:
  analyze_two_pulse.py --ccm <ccm_run_dir> [--somm <somm_run_dir>] \
     --r-ext 100 --r-c 20 --r-b 128 [--wavefile rpsi4_real_0100.txt] [--out fig.pdf]
Set the radii to the actual S3 configuration. The waveform column convention
matches plot_x2_structural.py (col 0 = time, col 3 = the (2,0) real part).
"""
import argparse
import pathlib
import sys

import numpy as np

try:
    from scipy.signal import find_peaks
except Exception:  # scipy optional; fall back to a simple local-max finder
    find_peaks = None


def load_waveform(run_dir, wavefile):
    p = pathlib.Path(run_dir) / "waveforms" / wavefile
    if not p.is_file():
        # tolerate a flat layout
        p = pathlib.Path(run_dir) / wavefile
    if not p.is_file():
        raise FileNotFoundError(f"no waveform at {p}")
    a = np.loadtxt(p, comments="#")
    return a[:, 0], a[:, 3]  # time, rPsi4^(2,0) real


def peaks(t, w, prominence_frac=0.05):
    """Return (times, amps) of |w| local maxima above prominence_frac*max|w|."""
    aw = np.abs(w)
    thr = prominence_frac * aw.max() if aw.max() > 0 else 0.0
    if find_peaks is not None:
        idx, _ = find_peaks(aw, height=thr, prominence=thr)
    else:
        idx = [i for i in range(1, len(aw) - 1)
               if aw[i] > aw[i - 1] and aw[i] >= aw[i + 1] and aw[i] >= thr]
        idx = np.asarray(idx, dtype=int)
    return t[idx], w[idx], aw[idx]


def predicted_times(r_ext, r_c, r_b):
    t1 = r_ext - r_c
    t2 = t1 + 2.0 * (r_b - r_ext)          # reflection off the outer boundary
    t3 = t2 + 2.0 * r_b                      # second bounce
    return {"primary": t1, "secondary": t2, "tertiary": t3}


def nearest(peaktimes, t_pred, tol):
    if len(peaktimes) == 0:
        return None
    d = np.abs(np.asarray(peaktimes) - t_pred)
    i = int(np.argmin(d))
    return (float(peaktimes[i]), float(d[i])) if d[i] <= tol else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccm", required=True, help="CCM run dir")
    ap.add_argument("--somm", default=None, help="paired Sommerfeld run dir")
    ap.add_argument("--r-ext", type=float, required=True)
    ap.add_argument("--r-c", type=float, default=20.0)
    ap.add_argument("--r-b", type=float, required=True, help="outer boundary / worldtube radius")
    ap.add_argument("--wavefile", default=None, help="rpsi4 file name (default rpsi4_real_<rext>.txt)")
    ap.add_argument("--tol", type=float, default=5.0, help="arrival-time match tolerance")
    ap.add_argument("--out", default="/data/haiyangw/claude/z4c-CMM/results/numerical/x2_xcts/two_pulse.pdf")
    args = ap.parse_args()

    wavefile = args.wavefile or f"rpsi4_real_{int(round(args.r_ext)):04d}.txt"
    pred = predicted_times(args.r_ext, args.r_c, args.r_b)
    print(f"predicted arrival times (r_ext={args.r_ext}, r_c={args.r_c}, r_b={args.r_b}): {pred}")

    t_ccm, w_ccm = load_waveform(args.ccm, wavefile)
    pk_t, pk_w, pk_a = peaks(t_ccm, w_ccm)
    print(f"CCM peaks (t, amp): {list(zip(np.round(pk_t,2), np.round(pk_w,4)))}")
    for name, tp in pred.items():
        m = nearest(pk_t, tp, args.tol)
        print(f"  {name}: predicted t={tp:.2f} -> "
              + (f"measured t={m[0]:.2f} (|dt|={m[1]:.2f})" if m else "NO peak within tol"))

    somm_secondary = ccm_secondary = None
    t_somm = w_somm = None
    if args.somm:
        t_somm, w_somm = load_waveform(args.somm, wavefile)
        ps_t, ps_w, ps_a = peaks(t_somm, w_somm)
        ms = nearest(ps_t, pred["secondary"], args.tol)
        mc = nearest(pk_t, pred["secondary"], args.tol)
        # secondary amplitudes (transparency metric: CCM should suppress it)
        if ms is not None:
            somm_secondary = float(np.interp(ms[0], t_somm, np.abs(w_somm)))
        ccm_secondary = float(np.interp(pred["secondary"], t_ccm, np.abs(w_ccm)))
        print(f"SECONDARY (spurious reflection) amplitude: "
              f"Sommerfeld~{somm_secondary}, CCM~{ccm_secondary}")
        if somm_secondary and ccm_secondary is not None:
            print(f"  CCM/Somm secondary ratio = {ccm_secondary/somm_secondary:.3e} "
                  f"(transparency: << 1 means CCM suppressed the reflection)")

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pathlib.Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
        if t_somm is not None:
            ax.plot(t_somm, w_somm, color="0.6", lw=1.0, label="Sommerfeld (reflective)")
        ax.plot(t_ccm, w_ccm, "b-", lw=1.2, label="Z4c-CCM (transparent)")
        for name, tp in pred.items():
            ax.axvline(tp, ls=":", color="r", lw=0.8)
            ax.text(tp, ax.get_ylim()[1], name[:4], fontsize=7, color="r",
                    rotation=90, va="top")
        ax.set_xlabel("Time"); ax.set_ylabel(r"$r\psi_4^{(2,0)}$ at $r_{\rm ext}$")
        ax.set_title(f"X=2 Teukolsky two-pulse gate (XCTS ID, r_b={args.r_b})", fontsize=10)
        ax.legend(fontsize=8)
        fig.savefig(args.out)
        print(f"wrote {args.out}")
    except Exception as e:
        print(f"(figure skipped: {e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
