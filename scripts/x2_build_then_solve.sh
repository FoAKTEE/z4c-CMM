#!/usr/bin/env bash
# Mission 6 orchestrator: wait for the SolveXcts rebuild to finish, then
# (on success) run the X=2 Teukolsky XCTS solve. One tracked operation so
# the loop need not spin through the ~20-min compile + the solve.
set -u
SRC=/data/haiyangw/nr/spectre/src/PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.cpp
BIN=/data/haiyangw/nr/spectre/build/bin/SolveXcts
LOG=/tmp/spectre_build_s1.log
echo "$(date -u +%FT%TZ) orchestrator: waiting for SolveXcts build"
for i in $(seq 1 90); do            # up to 90 min
  if grep -q "BUILD_EXIT=0" "$LOG" 2>/dev/null && [ "$BIN" -nt "$SRC" ]; then
    echo "$(date -u +%FT%TZ) build OK (binary fresh); launching XCTS solve"
    bash /data/haiyangw/claude/z4c-CMM/scripts/run_xcts_teukolsky_x2.sh 48
    echo "$(date -u +%FT%TZ) orchestrator DONE"
    exit 0
  fi
  if grep -qE "BUILD_EXIT=[1-9]|Error [0-9]" "$LOG" 2>/dev/null; then
    echo "$(date -u +%FT%TZ) BUILD FAILED — see $LOG"; tail -20 "$LOG"; exit 2
  fi
  sleep 60
done
echo "$(date -u +%FT%TZ) TIMEOUT waiting for build"; exit 3
