#!/usr/bin/env bash
# N12b: convergence of the AthenaK runs to the EXACT analytic finite-radius
# waveform (verify_n12_exact_psi4.py table), 6th-order stencils (nghost=4).
# Ladder on the paper box (boundary 41): h = 82/164, 82/248, 82/296.
# Comparison window t in [40, 58]: the (2,0) pulse passes the extraction
# sphere r = 36 at t ~ 50-60 and the earliest boundary-reflection
# contamination of r = 36 arrives at t ~ 58 (outgoing tail reaches 41 at
# ~53, returns 36 at ~58). tlim = 60 keeps every rung under the 10-min cap
# (328^3 measured-projected to ~670 s at 0.0875 s/cycle-per-164^3-unit and
# was replaced by 296^3, ~435 s projected).
# Usage: scripts/test_athenak_teuk_ana.sh <athena-binary> <tag> [launcher...]
set -eu
BIN=$(realpath "$1"); TAG="$2"; shift 2; LAUNCH=("$@")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT="$ROOT/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput"
OUT="$ROOT/results/numerical/athenak_teuk_ana/$TAG"
mkdir -p "$OUT"

run() { # name nx mb
  local d="$OUT/$1"; mkdir -p "$d"
  local t0=$SECONDS
  ( cd "$d" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" \
      problem/amp=1.0e-5 z4c/ccm_amp=1.0e-5 z4c/ccm=false \
      mesh/nghost=4 time/tlim=60.0 \
      mesh/nx1="$2" mesh/nx2="$2" mesh/nx3="$2" \
      meshblock/nx1="$3" meshblock/nx2="$3" meshblock/nx3="$3" \
      > run.log 2>&1 ) \
    || { echo "RUN FAILED: $d"; tail -5 "$d/run.log"; exit 1; }
  echo "done: $d  ($((SECONDS - t0)) s)"
}

run L1 164  82
run L2 248 124
run L3 296 148

python3 "$ROOT/scripts/check_athenak_teuk_ana.py" "$OUT"
