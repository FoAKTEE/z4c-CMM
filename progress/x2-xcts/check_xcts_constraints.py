#!/usr/bin/env python3
# S2 gate: read the XCTS reductions H5, report Newton convergence + the
# Hamiltonian/Momentum constraint L2-norm history. PASS if converged and
# final H/M norms <= 1e-8 (discretization-level, not O(X^2)=O(4)).
import sys, glob
import numpy as np, h5py
run = sys.argv[1] if len(sys.argv) > 1 else "/data/haiyangw/nr/spectre/x2_xcts_run"
red = glob.glob(f"{run}/XctsTeukolskyX2Reductions*.h5")
if not red:
    print("NO reductions file in", run); sys.exit(1)
f = h5py.File(red[0], "r")
def find(name):
    hits = []
    f.visititems(lambda n, o: hits.append(n) if (isinstance(o, h5py.Dataset) and name.lower() in n.lower()) else None)
    return hits
print("datasets:", [n for n in find("")][:30])
for key in ("Norms", "HamiltonianConstraint", "MomentumConstraint"):
    for ds in find(key):
        d = np.asarray(f[ds])
        if d.ndim == 2 and d.shape[0] >= 1:
            print(f"{ds}: shape {d.shape}, last row {d[-1]}")
