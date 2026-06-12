#!/usr/bin/env bash
# N13 battery: constraint-preserving Gamma-tilde row (z4c/cpbc=true) vs the
# stock all-Sommerfeld boundary, on the CORRECTED Teukolsky solution (N12),
# 6th-order stencils, paper box (boundary 41), X = 1e-5.
#   flat  : Minkowski (amp = 0) + CPBC, 64^3, t = 10  — exact preservation
#   somm  : stock boundary,   164^3, t = 80
#   cpbc  : CPBC row,         164^3, t = 80
#   cpbc2 : CPBC row,         248^3, t = 80  — convergence partner
#   ccmcp : CPBC + CCM datum, 164^3, t = 80  — full Z4c-CCM boundary operator
# Usage: scripts/test_athenak_cpbc.sh <athena-binary> <tag> [launcher...]
set -eu
BIN=$(realpath "$1"); TAG="$2"; shift 2; LAUNCH=("$@")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT="$ROOT/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput"
OUT="$ROOT/results/numerical/athenak_cpbc/$TAG"
mkdir -p "$OUT"

run() { # name extra-args...
  local d="$OUT/$1"; shift; mkdir -p "$d"
  local t0=$SECONDS
  ( cd "$d" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" \
      mesh/nghost=4 problem/amp=1.0e-5 z4c/ccm_amp=1.0e-5 "$@" \
      > run.log 2>&1 ) \
    || { echo "RUN FAILED: $d"; tail -5 "$d/run.log"; exit 1; }
  echo "done: $d  ($((SECONDS - t0)) s)"
}

run flat  z4c/cpbc=true z4c/ccm=false problem/amp=0.0 z4c/ccm_amp=0.0 \
    time/tlim=10.0 mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 \
    meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=32
run somm  z4c/cpbc=false z4c/ccm=false
run cpbc  z4c/cpbc=true  z4c/ccm=false
run cpbc2 z4c/cpbc=true  z4c/ccm=false \
    mesh/nx1=248 mesh/nx2=248 mesh/nx3=248 \
    meshblock/nx1=124 meshblock/nx2=124 meshblock/nx3=124
run ccmcp z4c/cpbc=true  z4c/ccm=true z4c/ccm_mode=3

python3 "$ROOT/scripts/check_athenak_cpbc.py" "$OUT"
