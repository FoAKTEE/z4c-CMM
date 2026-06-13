#!/usr/bin/env bash
# x2_full_campaign.sh — USER-AUTHORIZED long campaign (2026-06-13):
# X=2 Teukolsky on 8 GPUs, up to 2 days. Reference (causally
# disconnected, 768^3 box +-192) + boundary-41 Somm/CCM pairs at
# h = 0.5 / 0.25 / 0.125. All: diss=0.1 kappa1=0.1 nghost=4 CFL 0.25
# tlim=80 extraction r=36 (the att2-stable X=2 configuration).
set -u
ROOT=/data/haiyangw/claude/z4c-CMM
EXE=$ROOT/athenak/build_gpu_teuk/src/athena
INP=$ROOT/athenak/inputs/z4c/ccm/z4c_teuk_paper.athinput
BASE=$ROOT/results/numerical/ccm_repro_x2_full
STATUS=$BASE/STATUS.log
mkdir -p $BASE
echo "$(date -u +%FT%TZ) CAMPAIGN START (pid $$)" >> $STATUS

run () {
  local name=$1; shift
  local dir=$BASE/$name
  mkdir -p $dir && cd $dir
  echo "$(date -u +%FT%TZ) START $name" >> $STATUS
  mpirun -np 8 --mca pml ucx --mca osc ucx $EXE -i $INP \
    problem/amp=2.0 z4c/diss=0.1 z4c/damp_kappa1=0.1 mesh/nghost=4 \
    time/tlim=80 "$@" > run.log 2>&1
  local rc=$?
  local last=$(grep -E "^time=|cycle=" run.log | tail -1 | head -c 100)
  local nan=$(grep -ci nan run.log || true)
  echo "$(date -u +%FT%TZ) END $name rc=$rc nan=$nan last=[$last]" >> $STATUS
}

# B1/B2: h=0.5 boundary 41 (fast smoke of the campaign config)
run somm_h050 z4c/ccm=false \
  mesh/nx1=164 mesh/nx2=164 mesh/nx3=164 \
  meshblock/nx1=82 meshblock/nx2=82 meshblock/nx3=82
run ccm_h050 z4c/ccm=true z4c/ccm_mode=3 z4c/ccm_amp=2.0 \
  mesh/nx1=164 mesh/nx2=164 mesh/nx3=164 \
  meshblock/nx1=82 meshblock/nx2=82 meshblock/nx3=82

# R1: causally disconnected reference, 768^3 box +-192, h=0.5
run reference_768 z4c/ccm=false \
  mesh/nx1=768 mesh/nx2=768 mesh/nx3=768 \
  mesh/x1min=-192 mesh/x1max=192 mesh/x2min=-192 mesh/x2max=192 \
  mesh/x3min=-192 mesh/x3max=192 \
  meshblock/nx1=384 meshblock/nx2=384 meshblock/nx3=384

# B3/B4: h=0.25 boundary 41
run somm_h025 z4c/ccm=false \
  mesh/nx1=328 mesh/nx2=328 mesh/nx3=328 \
  meshblock/nx1=164 meshblock/nx2=164 meshblock/nx3=164
run ccm_h025 z4c/ccm=true z4c/ccm_mode=3 z4c/ccm_amp=2.0 \
  mesh/nx1=328 mesh/nx2=328 mesh/nx3=328 \
  meshblock/nx1=164 meshblock/nx2=164 meshblock/nx3=164

# B5/B6: h=0.125 boundary 41
run somm_h0125 z4c/ccm=false \
  mesh/nx1=656 mesh/nx2=656 mesh/nx3=656 \
  meshblock/nx1=328 meshblock/nx2=328 meshblock/nx3=328
run ccm_h0125 z4c/ccm=true z4c/ccm_mode=3 z4c/ccm_amp=2.0 \
  mesh/nx1=656 mesh/nx2=656 mesh/nx3=656 \
  meshblock/nx1=328 meshblock/nx2=328 meshblock/nx3=328

echo "$(date -u +%FT%TZ) CAMPAIGN COMPLETE" >> $STATUS
