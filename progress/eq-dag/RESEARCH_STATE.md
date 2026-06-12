# z4c-CMM — Equation-DAG Mission — Verified Research State

Mission: import arXiv:1010.0523v2, arXiv:2007.01339, arXiv:2308.10361 (stage 0),
decompose all three (stage 1), and store **all equations** as a directed acyclic
graph in the knowledge database (`knowledge-database/paper_<id>/nodes.jsonl`).
Phase: mission complete (stage 0 + stage 1 + equation DAG admitted).
Branch: `main`. Loop state: `.claude/ralph-loop.local.md` (session `z4c-cmm-eq-dag-20260612`).

## Source Library

| ID | Source | Kind | Status | Notes |
|---|---|---|---|---|
| 1010.0523v2 | ref-paper/arxiv-1010.0523v2/ | paper | [SOLID] imported, tex source | Z4c constraint-preserving boundary conditions (Ruiz, Hilditch, Bernuzzi) |
| 2007.01339 | ref-paper/arxiv-2007.01339/ | paper | [SOLID] imported, tex source | improved CCE system (Moxon, Scheel, Teukolsky) |
| 2308.10361 | ref-paper/arxiv-2308.10361/ | paper | [SOLID] imported, tex source | 3D Cauchy-characteristic matching (Ma et al.) |

## Working Context

Per-paper symbol tables live in `reformulate/z4c-CMM/paper_<id>/convention.md`.

## Active Claims

(none — mission claims discharged; stage-2 work seeds live in `reformulate/z4c-CMM/paper_<id>/{claims,obligations,result_seed}.md`)

## Accepted Results Log

| Claim | Evidence type | Evidence / verifier | Assumptions / deps | Status | Open obligations |
|---|---|---|---|---|---|
| Three papers imported with tex source + checksummed provenance | artifact + checksum | ref-paper/arxiv-<id>/PROVENANCE.md; structure_check OK | arXiv e-print availability | [SOLID] | none |
| Stage-1 decomposition complete (10 artifacts × 3 papers) | literature-grounded transcription | reformulate/z4c-CMM/paper_<id>/ | transcription fidelity, not independent derivation | [SOLID] | stage-2 items in per-paper obligations.md |
| All equations of all three papers stored as an acyclic dependency graph in the knowledge DB (325 equation nodes: 116 + 103 + 106; 202 tex labels covered 80/48/74) | automated check | `python3 scripts/eqdag_check.py` → OVERALL: PASS; certified output: results/numerical/eqdag_check_all.txt | label-coverage gate is against tex \label{} + derivation.md registry; edges mirror logic.md | [SOLID] | none |

## Next Work Steps

- `[FUTURE]` stage-2 implementation per `reformulate/z4c-CMM/paper_<id>/implementation_plan_python.md` (T-clusters), starting with the per-paper obligations files (suspected tex typos in 1010.0523v2 eq:sph_bc_last and the 2007.01339 eq:KDef sign discrepancy are recorded there).
- `[FUTURE]` stage-3 adversarial validation of transcription fidelity on a sampled equation subset.

## Appendix

### Abandoned Methods

(none)

### Audit references

- knowledge DB: `knowledge-database/paper_arxiv-<id>/{nodes.jsonl,summary.csv}`; HTML views under `human-read/knowledge-database/` (gitignored).
- mission verifier: `scripts/eqdag_check.py`; certified run: `results/numerical/eqdag_check_all.txt`.
- alignment.md anomaly: upstream commit `663804c` in `phys-agentic-loop/` replaced the §15–§21 protocol with an unrelated "PUA" persona skill; §-references in specs/templates/hooks now dangle. Substantive binding rules survive in the UserPromptSubmit wrapper of `phys-agentic-loop/.claude/settings.json`. Flagged to owner 2026-06-12; methodology repo left unmodified.
