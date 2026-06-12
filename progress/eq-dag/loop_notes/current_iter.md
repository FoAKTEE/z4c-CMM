# Iterations 2–4 — stage-1 decomposition + equation DAG (three parallel sub-agents); iteration 5 — consolidation + completion

1. **Paper anchor** — S1-decompose-{1010.0523v2, 2007.01339, 2308.10361}: all equation environments of each paper.
2. **What shipped** — 10 stage-1 artifacts per paper under `reformulate/z4c-CMM/paper_<id>/`; 325 equation rows appended to the knowledge DB via the validating CLI (116 + 103 + 106), predecessors mirroring each `logic.md`; summary.csv regenerated; HTML views rendered.
3. **Next-3 roadmap** — mission exit (this iteration); stage-2 and stage-3 work deferred as [FUTURE] in RESEARCH_STATE.
4. **Simplification flag** — n/a (no metric-bearing code iterations).
5. **Verifier output** (orchestrator re-run, certified copy at results/numerical/eqdag_check_all.txt):
   - arxiv-1010.0523v2: acyclicity PASS (117 nodes) · edge closure PASS · tex-label coverage PASS 80/80 · derivation registry PASS · 116 equation nodes (80 tex-born + 36 auto)
   - arxiv-2007.01339: acyclicity PASS (104 nodes) · edge closure PASS · coverage PASS 48/48 · registry PASS · 103 equation nodes (48 + 55 auto)
   - arxiv-2308.10361: acyclicity PASS (107 nodes) · edge closure PASS · coverage PASS 74/74 · registry PASS · 106 equation nodes (74 + 38 auto)
   - OVERALL: PASS (exit=0)
