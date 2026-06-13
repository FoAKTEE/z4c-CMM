#!/usr/bin/env python3
"""plot_xcts_convergence.py — human-readable XCTS elliptic-solver convergence,
in the style of the SpECTRE BBH initial-data convergence plot.

Reads the SolveXctsVacuum reductions H5 (default
/data/haiyangw/nr/spectre/x2_xcts_run/XctsTeukolskyX2Reductions*.h5) and renders:

  Top    : nonlinear (Newton-Raphson) residual and linear-solver (GMRES /
           multigrid) residual vs. cumulative linear-solver iteration, log-y —
           the canonical SpECTRE elliptic-convergence figure.
  Bottom : the Hamiltonian and Momentum constraint L2 norms vs. iteration
           (the physically meaningful XCTS gate: they fall from O(X^2)=O(4) of
           the raw free data to the discretization floor of the solved ID).

SpECTRE writes each reduction as a 2-D dataset carrying a `Legend` attribute
naming the columns; this script discovers the relevant datasets/columns by
their legends rather than hard-coding indices, so it adapts to the actual file.

Usage: plot_xcts_convergence.py [run_dir] [out.pdf]
"""
import glob
import sys
from pathlib import Path

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = sys.argv[1] if len(sys.argv) > 1 else "/data/haiyangw/nr/spectre/x2_xcts_run"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(
    "/data/haiyangw/claude/z4c-CMM/results/numerical/x2_xcts/xcts_convergence.pdf")


def legend(dset):
    """Return the lower-cased column names from a dataset's Legend attribute."""
    leg = dset.attrs.get("Legend")
    if leg is None:
        return []
    out = []
    for c in leg:
        out.append((c.decode() if isinstance(c, bytes) else str(c)))
    return out


def collect(h5):
    """Walk the file; return {path: (array, legend)} for every 2-D dataset."""
    found = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and obj.ndim == 2 and obj.shape[0] >= 1:
            found[name] = (np.asarray(obj), legend(obj))
    h5.visititems(visit)
    return found


def main():
    files = sorted(glob.glob(f"{RUN}/XctsTeukolskyX2Reductions*.h5"))
    if not files:
        print(f"NO reductions file under {RUN}")
        return 1
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(files[0], "r") as h5:
        data = collect(h5)
    print(f"reductions file: {files[0]}")
    for name, (arr, leg) in data.items():
        print(f"  {name}: shape={arr.shape} legend={leg}")

    # Identify residual datasets (nonlinear + linear) and constraint-norm sets.
    resid_sets = {}
    for name, (arr, leg) in data.items():
        low = [c.lower() for c in leg]
        ri = next((i for i, c in enumerate(low) if "residual" in c), None)
        ii = next((i for i, c in enumerate(low)
                   if "iteration" in c or "cumulative" in c), None)
        if ri is not None:
            resid_sets[name] = (arr, ii, ri)

    ham = {}
    for name, (arr, leg) in data.items():
        low = [c.lower() for c in leg]
        ii = next((i for i, c in enumerate(low) if "iteration" in c), 0)
        for i, c in enumerate(low):
            if "hamiltonian" in c:
                ham.setdefault("Hamiltonian", (name, arr, ii, i))
            if "momentum" in c:
                ham.setdefault("Momentum", (name, arr, ii, i))

    fig, ax = plt.subplots(2, 1, figsize=(7, 7.6), constrained_layout=True)

    # Top: solver residuals.
    if resid_sets:
        for name, (arr, ii, ri) in sorted(resid_sets.items()):
            x = arr[:, ii] if ii is not None else np.arange(arr.shape[0])
            lab = name.split("/")[-1].replace(".dat", "")
            ax[0].semilogy(x, np.abs(arr[:, ri]), marker=".", ms=3, lw=1.0,
                           label=lab)
        ax[0].legend(fontsize=8)
    else:
        ax[0].text(0.5, 0.5, "no residual dataset found", ha="center",
                   transform=ax[0].transAxes)
    ax[0].set_xlabel("cumulative linear-solver iteration")
    ax[0].set_ylabel("elliptic-solver residual")
    ax[0].set_title("XCTS solve for the Teukolsky $X{=}2$ free data: "
                    "elliptic convergence", fontsize=10)
    ax[0].grid(True, which="both", alpha=0.3)

    # Bottom: constraint norms.
    if ham:
        for key, (name, arr, ii, ci) in ham.items():
            x = arr[:, ii] if arr.shape[1] > ii else np.arange(arr.shape[0])
            ax[1].semilogy(x, np.abs(arr[:, ci]), marker=".", ms=3, lw=1.0,
                           label=f"{key} constraint $L_2$")
        ax[1].legend(fontsize=8)
    else:
        ax[1].text(0.5, 0.5, "no constraint-norm dataset found", ha="center",
                   transform=ax[1].transAxes)
    ax[1].set_xlabel("nonlinear (Newton) step")
    ax[1].set_ylabel(r"constraint $L_2$ norm")
    ax[1].grid(True, which="both", alpha=0.3)

    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
