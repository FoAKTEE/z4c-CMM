#!/usr/bin/env bash
# Teukolsky-wave reproduction battery (arXiv:2308.10361 Sec. V.A protocol):
#   somm : Sommerfeld boundary (ccm off),      box 12, 48^3
#   ccm  : CCM with Teukolsky self-datum,      box 12, 48^3
#   ref  : causally-disconnected reference,    box 24, 96^3 (same h)
# Usage: scripts/test_athenak_teuk.sh <athena-binary> <tag> <amp> [launcher...]
set -eu
BIN=$(realpath "$1"); TAG="$2"; AMP="$3"; shift 3; LAUNCH=("$@")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT="$ROOT/athenak/inputs/z4c/ccm/z4c_teuk.athinput"
OUT="$ROOT/results/numerical/athenak_teuk/$TAG"
mkdir -p "$OUT"

run() { # name extra-args...
  local d="$OUT/$1"; shift; mkdir -p "$d"
  ( cd "$d" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" problem/amp="$AMP" \
      z4c/ccm_amp="$AMP" "$@" > run.log 2>&1 ) \
    || { echo "RUN FAILED: $d"; tail -5 "$d/run.log"; exit 1; }
  echo "done: $d"
}

run somm z4c/ccm=false
run ccm  z4c/ccm=true z4c/ccm_mode=3
run ref  z4c/ccm=false \
    mesh/nx1=96 mesh/nx2=96 mesh/nx3=96 \
    mesh/x1min=-24 mesh/x1max=24 mesh/x2min=-24 mesh/x2max=24 \
    mesh/x3min=-24 mesh/x3max=24

python3 "$ROOT/scripts/check_athenak_teuk.py" "$OUT" "$AMP"
