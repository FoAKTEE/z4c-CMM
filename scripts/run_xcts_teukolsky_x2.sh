#!/usr/bin/env bash
# S2: run the X=2 Teukolsky XCTS solve. Usage: bash run_xcts_teukolsky_x2.sh [ncores]
set -u
source /data/haiyangw/claude/z4c-CMM/scripts/spectre_xcts_env.sh
NCORE=${1:-96}
OUT=/data/haiyangw/nr/spectre/x2_xcts_run
mkdir -p "$OUT" && cd "$OUT"
cp /data/haiyangw/nr/spectre/SolveXctsTeukolskyX2.yaml .
echo "$(date -u +%FT%TZ) XCTS solve START (ncore=$NCORE)"
# Vacuum XCTS executable (no hydro observation; registers TeukolskyWave).
# mpi-linux-x86_64 Charm++: +pN must go via charmrun (-> setarch -R mpirun -np N).
"$SPECTRE_BUILD/bin/charmrun" +p"$NCORE" "$SPECTRE_BUILD/bin/SolveXctsVacuum" \
  --input-file SolveXctsTeukolskyX2.yaml > xcts_solve.log 2>&1
echo "XCTS_EXIT=$? $(date -u +%FT%TZ)"
ls -la XctsTeukolskyX2Volume*.h5 XctsTeukolskyX2Reductions.h5 2>/dev/null
tail -5 xcts_solve.log
