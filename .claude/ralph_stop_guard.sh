#!/usr/bin/env bash
# ralph_stop_guard.sh — progress-aware Stop-hook guard.
#
# Consults _common/loop/loop_gate.py, which decides whether the Ralph loop
# should keep going. The guard blocks the stop (forces continuation) ONLY while
# the gate verdict is `continue` — i.e. while the loop is making verified
# progress and is within its iteration / wall-clock budget. On any `halt:*`
# verdict it stays silent and lets the agent stop gracefully; loop_gate has
# already written .claude/HUMAN_REVIEW_REQUIRED.md for the halts that need
# human attention (no_progress / stuck_counter / time_budget / max_iterations).
#
# PROTOCOL SELF-INJECTION: when the gate decision carries a due `self_inject`
# block (every `self_inject_interval` iterations, default 5), the guard re-emits
# the listed protocol files (loop state + alignment + admission contract) inside
# the block reason, so the loop session re-reads its operating contract on a
# fixed cadence. Escalation directives from the gate (`escalation_required` /
# `escalation_suggested`) are surfaced the same way.
#
# If loop_gate.py is unavailable (missing file or no python3), the guard falls
# back to the legacy `iteration < max_iterations` rule, so a tooling error never
# silently kills the loop and never silently lets it run away unbounded.

set -u

# Non-loop sessions must never be held open by the guard: the independent
# pipeline-feedback observer exports CLAUDE_OBSERVER=1; any other one-shot
# tool/headless session may export RALPH_GUARD_BYPASS=1. Bail out for both.
[ -n "${CLAUDE_OBSERVER:-}" ] && exit 0
[ -n "${RALPH_GUARD_BYPASS:-}" ] && exit 0

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || printf '%s' "$(cd "$(dirname "$0")/.." && pwd)")"
STATE="$REPO_ROOT/.claude/ralph-loop.local.md"

[ -f "$STATE" ] || exit 0

# Terminal freeze (cheapest point): when the loop is disarmed, do not even
# evaluate the gate — a mutating `decide` on every post-disarm session stop is
# how disarmed terminal records used to drift (F-020). Sessions stop freely.
ACTIVE=$(awk '/^---$/{n++; next} n==1 && /^active:/ {print $2}' "$STATE" | head -1)
[ "$ACTIVE" = "true" ] || exit 0

emit_block() {
    # $1 = reason text; JSON-encode via python3 when available.
    local reason="$1" encoded
    if command -v python3 >/dev/null 2>&1; then
        encoded="$(printf '%s' "$reason" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')"
    else
        encoded="\"${reason//\"/\'}\""
    fi
    printf '{"decision":"block","reason":%s}\n' "$encoded"
}

legacy_guard() {
    # Legacy behavior: block while active and iteration < max_iterations.
    local iter max
    iter=$(awk '/^---$/{n++; next} n==1 && /^iteration:/ {print $2}' "$STATE" | head -1)
    max=$(awk '/^---$/{n++; next} n==1 && /^max_iterations:/ {print $2}' "$STATE" | head -1)
    [ -n "$iter" ] && [ -n "$max" ] || exit 0
    case "$iter$max" in *[!0-9]*) exit 0 ;; esac
    if [ "$iter" -lt "$max" ]; then
        emit_block "Ralph loop active — iteration ${iter}/${max} (legacy guard; loop_gate.py unavailable). Advance the counter in .claude/ralph-loop.local.md, ship the next verified step, then continue. If stalled >30 min, run pipelines/6-escalation/spec.md with a 30-min cap."
    fi
    exit 0
}

# Locate loop_gate.py: an explicit override, the consumer's submodule, this
# repo's _common/ (canonical _common/loop/ path first, then the compatibility
# shim), or the methodology checkout that ships the guard.
GATE=""
for cand in \
    "${PHYS_AGENTIC_LOOP:-}/_common/loop/loop_gate.py" \
    "${PHYS_AGENTIC_LOOP:-}/_common/loop_gate.py" \
    "$REPO_ROOT/phys-agentic-loop/_common/loop/loop_gate.py" \
    "$REPO_ROOT/phys-agentic-loop/_common/loop_gate.py" \
    "$REPO_ROOT/_common/loop/loop_gate.py" \
    "$REPO_ROOT/_common/loop_gate.py" \
    "$(cd "$(dirname "$0")/.." 2>/dev/null && pwd)/_common/loop/loop_gate.py" \
    "$(cd "$(dirname "$0")/.." 2>/dev/null && pwd)/_common/loop_gate.py"; do
    if [ -n "$cand" ] && [ -f "$cand" ]; then GATE="$cand"; break; fi
done

if [ -z "$GATE" ] || ! command -v python3 >/dev/null 2>&1; then
    legacy_guard
fi

out="$(python3 "$GATE" decide --repo-root "$REPO_ROOT" --state-file "$STATE" --write-gate 2>/dev/null)"
decision="$(printf '%s' "$out" | python3 -c 'import json,sys
try: print(json.load(sys.stdin)["decision"])
except Exception: pass' 2>/dev/null)"

# No usable verdict (loop_gate errored) -> legacy fallback rather than a silent stop.
[ -n "$decision" ] || legacy_guard

case "$decision" in
    continue)
        # Build the block reason in python: base reason + escalation directives
        # + the self-injection emit (file contents) when the gate says it is due.
        printf '%s' "$out" | python3 -c '
import json, sys
d = json.load(sys.stdin)
parts = []
parts.append(
    "Ralph loop progress gate: %s. Advance the counter in "
    ".claude/ralph-loop.local.md, ship the next verified step (append a "
    "result/knowledge-ledger row; use error rows for trial activity), then "
    "continue. The gate halts the loop automatically if progress stalls."
    % d.get("reason", ""))
if "escalation_required" in d:
    parts.insert(0, "*** ESCALATION REQUIRED (forced): %s ***" % d["escalation_required"])
if "escalation_suggested" in d:
    parts.append("ESCALATION SUGGESTED: %s" % d["escalation_suggested"])
si = d.get("self_inject") or {}
if si.get("due"):
    if si.get("ok"):
        emit = ["", "===== PROTOCOL SELF-INJECTION (every %s iterations; binding — re-read and apply) =====" % si.get("interval")]
        for f in si.get("files", []):
            try:
                body = open(f, encoding="utf-8", errors="replace").read()
            except OSError as e:
                body = "(unreadable: %s)" % e
            emit.append("----- BEGIN %s -----\n%s\n----- END %s -----" % (f, body, f))
        emit.append("===== END PROTOCOL SELF-INJECTION =====")
        parts.append("\n".join(emit))
    else:
        parts.append(
            "WARNING: protocol self-injection FAILED this window "
            "(missing: %s) — failure %s/%s. Restore the missing protocol files; "
            "at the limit the gate forces a pipelines/6-escalation run."
            % (", ".join(si.get("missing", [])), si.get("failures"), si.get("fail_limit")))
print(json.dumps({"decision": "block", "reason": "\n\n".join(parts)}))
' 2>/dev/null || emit_block "Ralph loop progress gate: continue. Advance the counter in .claude/ralph-loop.local.md, ship the next verified step, then continue."
        ;;
    halt:*)
        # Graceful stop. loop_gate wrote .claude/HUMAN_REVIEW_REQUIRED.md where a
        # human needs to look. Emit nothing so the normal stop proceeds.
        : ;;
esac

exit 0
