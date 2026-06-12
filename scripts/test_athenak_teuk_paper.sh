#!/usr/bin/env bash
# EXACT-PARAMETER Teukolsky battery (arXiv:2308.10361 Sec. V.A figures):
# pulse r_c = 20, tau = 2; Cauchy outer boundary (worldtube analog) at 41;
# psi4 extraction at r = 36; causally-disconnected reference at 104:
# boundary-sourced contamination needs t_b < 80 - (104 - 36) = 12 to reach
# the extraction sphere by tlim = 80, but at t = 12 the outgoing front
# (r ~ 20 + t + 4 tau) is at r ~ 40, so the field at 104 is zero to
# exp(-((104-40)/tau)^2) = e^-1024 — exact zero in float64. (The 124-box
# variant with the unconditional 124 - 36 > tlim margin projects to ~580 s,
# breaching the 10-min per-run cap; measured 0.90 s/cycle at 496^3 on 4 A100s.)
# Runs (each within the 10-min cap):
#   somm     : Sommerfeld,             164^3 (h = 0.5)
#   somm_mid : Sommerfeld,             132^3 (resolution-error estimate)
#   ccm      : CCM Teukolsky datum,    164^3
#   ccm_mid  : CCM,                    132^3 (resolution-error estimate)
#   ccm_low  : CCM,                    104^3 (third resolution, constraints)
#   ref      : reference,  box 104,    416^3 (same h = 0.5)
# Usage: scripts/test_athenak_teuk_paper.sh <athena-binary> <tag> <amp> [launcher...]
set -eu
BIN=$(realpath "$1"); TAG="$2"; AMP="$3"; shift 3; LAUNCH=("$@")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT="$ROOT/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput"
OUT="$ROOT/results/numerical/athenak_teuk_paper/$TAG"
mkdir -p "$OUT"

run() { # name extra-args...
  local d="$OUT/$1"; shift; mkdir -p "$d"
  local t0=$SECONDS
  ( cd "$d" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" problem/amp="$AMP" \
      z4c/ccm_amp="$AMP" "$@" > run.log 2>&1 ) \
    || { echo "RUN FAILED: $d"; tail -5 "$d/run.log"; exit 1; }
  echo "done: $d  ($((SECONDS - t0)) s)"
}

run somm     z4c/ccm=false
run somm_mid z4c/ccm=false \
    mesh/nx1=132 mesh/nx2=132 mesh/nx3=132 \
    meshblock/nx1=66 meshblock/nx2=66 meshblock/nx3=66
run ccm      z4c/ccm=true z4c/ccm_mode=3
run ccm_mid  z4c/ccm=true z4c/ccm_mode=3 \
    mesh/nx1=132 mesh/nx2=132 mesh/nx3=132 \
    meshblock/nx1=66 meshblock/nx2=66 meshblock/nx3=66
run ccm_low  z4c/ccm=true z4c/ccm_mode=3 \
    mesh/nx1=104 mesh/nx2=104 mesh/nx3=104 \
    meshblock/nx1=52 meshblock/nx2=52 meshblock/nx3=52
run ref      z4c/ccm=false \
    mesh/nx1=416 mesh/nx2=416 mesh/nx3=416 \
    mesh/x1min=-104 mesh/x1max=104 mesh/x2min=-104 mesh/x2max=104 \
    mesh/x3min=-104 mesh/x3max=104 \
    meshblock/nx1=208 meshblock/nx2=208 meshblock/nx3=208

python3 "$ROOT/scripts/check_athenak_teuk_paper.py" "$OUT" "$AMP"
