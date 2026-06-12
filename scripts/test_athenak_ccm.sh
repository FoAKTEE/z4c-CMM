#!/usr/bin/env bash
# Z4c-CCM test battery for AthenaK (CPU or GPU binary).
# Usage: scripts/test_athenak_ccm.sh <athena-binary> <tag> [launcher...]
# e.g. 4-GPU: scripts/test_athenak_ccm.sh build_gpu_mpi/src/athena gpu4 \
#        mpirun -np 4
# Runs (each well under the 10-min budget):
#   T1 reduction : ccm_amp = 0     -> Minkowski stays flat to roundoff
#   T2 injection : ccm_amp = 1e-3  -> pulse enters; deviation O(amp)
#   T3 linearity : ccm_amp = 2e-3  -> response scales linearly
# Emits per-run .hst files and a summary via scripts/check_athenak_ccm.py.
set -eu
BIN=$(realpath "$1"); TAG="$2"; shift 2; LAUNCH=("$@")
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT="$ROOT/athenak/inputs/z4c/ccm/z4c_ccm_test.athinput"
OUT="$ROOT/results/numerical/athenak_ccm/$TAG"
mkdir -p "$OUT"

run() { # name amp
  local d="$OUT/$1"; mkdir -p "$d"
  ( cd "$d" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" z4c/ccm_amp="$2" \
      > run.log 2>&1 ) || { echo "RUN FAILED: $1"; tail -5 "$OUT/$1/run.log"; exit 1; }
  echo "done: $1 (amp=$2)"
}

run t1_reduction 0.0
run t2_injection 1.0e-3
run t3_linearity 2.0e-3

python3 "$ROOT/scripts/check_athenak_ccm.py" "$OUT"
