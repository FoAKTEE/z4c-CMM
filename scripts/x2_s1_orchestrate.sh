#!/usr/bin/env bash
# M6 S1 orchestrator: wait for the SolveXctsVacuum build to link, validate the
# input file, run the X=2 Teukolsky XCTS solve, and post-process (constraint
# gate + SpECTRE-style convergence figure). One autonomous chain; writes a
# concise running log so the parent can evaluate the S1 gate on completion.
set -u
BUILDLOG=/tmp/spectre_build_vacuum_s2b.log
LOG=/tmp/x2_s1_orchestrate.log
: > "$LOG"
say(){ echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

NCORE=${1:-32}   # XCTS solve PEs. 1600 elements: +p32 (~50 elem/rank) is the strong-scaling
                 # sweet spot -- +p96 (~17 elem/rank) was communication-bound (no speedup), +p48
                 # also under-utilized. 32 maximizes per-core efficiency and frees cores.
REPO=/data/haiyangw/claude/z4c-CMM
BIN=/data/haiyangw/nr/spectre/build/bin/SolveXctsVacuum
OUT=/data/haiyangw/nr/spectre/x2_xcts_run

# 1. wait for the build to finish (BUILD_EXIT line)
say "waiting for SolveXctsVacuum build link"
while ! grep -q "BUILD_EXIT=" "$BUILDLOG" 2>/dev/null; do sleep 15; done
BX=$(grep -o "BUILD_EXIT=[0-9]*" "$BUILDLOG" | tail -1 | cut -d= -f2)
say "build BUILD_EXIT=$BX"
if [ "$BX" != "0" ] || [ ! -x "$BIN" ]; then
    say "BUILD FAILED or binary missing — stopping; see $BUILDLOG"
    grep -E " error: |\] Error " "$BUILDLOG" | head -5 | tee -a "$LOG"
    exit 1
fi
ls -la --time-style=full-iso "$BIN" | tee -a "$LOG"

source "$REPO/scripts/spectre_xcts_env.sh"
mkdir -p "$OUT" && cd "$OUT"
cp /data/haiyangw/nr/spectre/SolveXctsTeukolskyX2.yaml .

# 2. validate the input file (fast). Informational: a non-zero exit here does
#    not abort — the solve re-parses and will fail fast on a real yaml error.
say "check-options ..."
"$BIN" --check-options --input-file SolveXctsTeukolskyX2.yaml > check_options.log 2>&1
say "check-options exit=$? (see check_options.log)"
tail -5 check_options.log | tee -a "$LOG"

# 3. run the XCTS solve. This Charm++ is the mpi-linux-x86_64 build, so +pN must
#    go through charmrun (which invokes `setarch -R mpirun -np N <bin>`), NOT
#    directly on the binary (that prints the charmrun usage and exits).
say "XCTS solve START (charmrun +p$NCORE)"
timeout 36000 "$(dirname "$BIN")/charmrun" +p"$NCORE" "$BIN" \
    --input-file SolveXctsTeukolskyX2.yaml > xcts_solve.log 2>&1
SX=$?
say "XCTS_EXIT=$SX"
tail -12 xcts_solve.log | tee -a "$LOG"
ls -la XctsTeukolskyX2*.h5 2>/dev/null | tee -a "$LOG"

# 4. post-process: constraint gate + convergence figure
say "--- constraint gate (check_xcts_constraints.py) ---"
python3 "$REPO/scripts/check_xcts_constraints.py" "$OUT" 2>&1 | tee -a "$LOG"
say "--- convergence figure (plot_xcts_convergence.py) ---"
python3 "$REPO/scripts/plot_xcts_convergence.py" "$OUT" \
    "$REPO/results/numerical/x2_xcts/xcts_convergence.pdf" 2>&1 | tee -a "$LOG"
say "S1 ORCHESTRATOR DONE"
