#!/usr/bin/env bash
# inject_infra.sh — cross-client reference implementation
#
# Reads the phys-agentic-loop methodology (alignment + _common + state templates)
# from disk and emits a <session-start-briefing> envelope. Consumer repos should
# copy this file into their own .claude/ and extend the tail with their local
# task context.
#
# METHODOLOGY-DIR INDIRECTION: consumers vendor the methodology as a subdir or
# git submodule (conventionally ./phys-agentic-loop), overridable via
# $PHYS_AGENTIC_LOOP. All emits resolve against $PAL_DIR, never against the
# consumer root blindly (F-002).
#
# Wired into .claude/settings.json as the SessionStart hook, and also usable as a
# manual bootstrap for clients that do not support hooks.

set -u

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || printf '%s' "$(cd "$(dirname "$0")/.." && pwd)")"

# Resolve the methodology root: explicit override, consumer submodule, or this
# repo itself (when the methodology IS the checkout).
PAL_DIR="${PHYS_AGENTIC_LOOP:-}"
if [ -z "$PAL_DIR" ] || [ ! -f "$PAL_DIR/alignment.md" ]; then
    if [ -f "$REPO_ROOT/phys-agentic-loop/alignment.md" ]; then
        PAL_DIR="$REPO_ROOT/phys-agentic-loop"
    else
        PAL_DIR="$REPO_ROOT"
    fi
fi

emit_file() {
    local label="$1" path="$2"
    printf '\n----- BEGIN %s (%s) -----\n' "$label" "$path"
    if [ -f "$path" ]; then
        cat "$path"
    else
        printf '(missing: %s)\n' "$path"
    fi
    printf '\n----- END %s -----\n' "$label"
}

printf '<session-start-briefing enforcement="MANDATORY">\n'
printf '\n=== PHYS-AGENTIC-LOOP INFRA (methodology: %s; consumer: %s) ===\n' "$PAL_DIR" "$REPO_ROOT"

emit_file "INDEX.md"                         "$PAL_DIR/INDEX.md"
emit_file "alignment.md"                     "$PAL_DIR/alignment.md"
emit_file "_common/contracts/markers.md"               "$PAL_DIR/_common/contracts/markers.md"
emit_file "_common/contracts/note_discipline.md"       "$PAL_DIR/_common/contracts/note_discipline.md"
emit_file "_common/contracts/progress_principles.md"   "$PAL_DIR/_common/contracts/progress_principles.md"
emit_file "_common/contracts/research_admission_contract.md" "$PAL_DIR/_common/contracts/research_admission_contract.md"
emit_file "_common/contracts/project_structure.md"     "$PAL_DIR/_common/contracts/project_structure.md"
emit_file "notes/research_state_template.md" "$PAL_DIR/notes/research_state_template.md"
emit_file "notes/multi_timescale_tracking_template.md" "$PAL_DIR/notes/multi_timescale_tracking_template.md"
emit_file "notes/research_note_template.md"  "$PAL_DIR/notes/research_note_template.md"
emit_file "notes/ralph_loop_local_template.md" "$PAL_DIR/notes/ralph_loop_local_template.md"

# --- project-structure check (contract: _common/contracts/project_structure.md) ---
printf '\n=== PROJECT STRUCTURE CHECK ===\n'
if command -v python3 >/dev/null 2>&1 && [ -f "$PAL_DIR/_common/structure_check.py" ]; then
    python3 "$PAL_DIR/_common/structure_check.py" --repo-root "$REPO_ROOT" 2>/dev/null \
        || printf '(structure_check failed to run)\n'
else
    printf '(structure_check unavailable)\n'
fi

printf '\n=== END PHYS-AGENTIC-LOOP INFRA ===\n'
printf '</session-start-briefing>\n'
