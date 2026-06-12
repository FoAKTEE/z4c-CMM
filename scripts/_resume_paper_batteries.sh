#!/usr/bin/env bash
set -eu
ROOT=/data/haiyangw/claude/z4c-CMM
BIN=$ROOT/athenak/build_gpu_teuk/src/athena
INPUT=$ROOT/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput
LAUNCH=(mpirun -np 4 --mca pml ucx --mca osc ucx)
D=$ROOT/results/numerical/athenak_teuk_paper/paper_X1e-5/ref
rm -rf "$D"; mkdir -p "$D"
t0=$SECONDS
( cd "$D" && timeout 600 "${LAUNCH[@]}" "$BIN" -i "$INPUT" problem/amp=1.0e-5 \
    z4c/ccm_amp=1.0e-5 z4c/ccm=false \
    mesh/nx1=416 mesh/nx2=416 mesh/nx3=416 \
    mesh/x1min=-104 mesh/x1max=104 mesh/x2min=-104 mesh/x2max=104 \
    mesh/x3min=-104 mesh/x3max=104 \
    meshblock/nx1=208 meshblock/nx2=208 meshblock/nx3=208 > run.log 2>&1 ) \
  || { echo "RUN FAILED: $D"; tail -5 "$D/run.log"; exit 1; }
echo "done: $D  ($((SECONDS - t0)) s)"
python3 $ROOT/scripts/check_athenak_teuk_paper.py \
  $ROOT/results/numerical/athenak_teuk_paper/paper_X1e-5 1.0e-5 \
  | tee $ROOT/results/numerical/athenak_teuk_paper/paper_X1e-5_summary.txt
echo "=== X=2 battery ==="
$ROOT/scripts/test_athenak_teuk_paper.sh "$BIN" paper_X2 2.0 "${LAUNCH[@]}" \
  | tee $ROOT/results/numerical/athenak_teuk_paper/paper_X2_summary.txt
