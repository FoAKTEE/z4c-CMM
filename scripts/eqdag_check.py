#!/usr/bin/env python3
"""eqdag_check.py — mission verifier for the z4c-CMM equation-DAG.

Per paper, checks:
  1. ACYCLICITY  — predecessors edges in knowledge-database/paper_<P>/nodes.jsonl
                   (latest non-amended row per node_id) form a DAG.
  2. EDGE CLOSURE — every predecessor named by a node exists as a node.
  3. LABEL COVERAGE — every \\label{...} attached to an equation environment in
                   the imported tex (ref-paper/arxiv-<id>/src/*.tex) appears in
                   the union of equation_labels over DB nodes, and every such
                   label appears verbatim in reformulate/<proj>/paper_<id>/derivation.md.
  4. COUNT REPORT — number of equation environments in tex vs labels covered
                   (unlabeled equations must be covered via eq:auto-* labels
                   registered in derivation.md; reported, and gated by the
                   derivation.md registry).

Exit 0 only if every hard check passes for every requested paper.
Usage: python3 scripts/eqdag_check.py [--paper ID ...] [--project z4c-CMM]
"""
import argparse, json, re, sys
from pathlib import Path

EQ_ENVS = r"equation|align|eqnarray|gather|multline|alignat|flalign"
BEGIN_RE = re.compile(r"\\begin\{(" + EQ_ENVS + r")\*?\}")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
COMMENT_RE = re.compile(r"(?<!\\)%.*")

def tex_equation_info(srcdir: Path):
    """Return (n_environments, labels_in_eq_envs) across all .tex files."""
    n_env, labels = 0, set()
    for tex in sorted(srcdir.glob("*.tex")):
        text = "\n".join(COMMENT_RE.sub("", ln) for ln in
                         tex.read_text(encoding="utf-8", errors="replace").splitlines())
        for m in re.finditer(r"\\begin\{(" + EQ_ENVS + r")(\*?)\}(.*?)\\end\{\1\2\}",
                             text, re.S):
            n_env += 1
            labels.update(LABEL_RE.findall(m.group(3)))
    return n_env, labels

def latest_nodes(jsonl: Path):
    rows = {}
    if jsonl.exists():
        for line in jsonl.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("status") != "amended":
                rows[r["node_id"]] = r
    return rows

def acyclic(nodes):
    """Kahn's algorithm over predecessors edges; returns (ok, cycle_members)."""
    preds = {nid: [p for p in (r.get("predecessors") or []) if p in nodes]
             for nid, r in nodes.items()}
    indeg = {nid: len(ps) for nid, ps in preds.items()}
    succ = {nid: [] for nid in nodes}
    for nid, ps in preds.items():
        for p in ps:
            succ[p].append(nid)
    queue = [n for n, d in indeg.items() if d == 0]
    seen = 0
    while queue:
        n = queue.pop()
        seen += 1
        for s in succ[n]:
            indeg[s] -= 1
            if indeg[s] == 0:
                queue.append(s)
    return seen == len(nodes), sorted(n for n, d in indeg.items() if d > 0)

def check_paper(root: Path, pid: str, project: str) -> bool:
    print(f"\n===== paper arxiv-{pid} =====")
    ok = True
    srcdir = root / f"ref-paper/arxiv-{pid}/src"
    deriv = root / f"reformulate/{project}/paper_{pid}/derivation.md"
    db = root / f"knowledge-database/paper_arxiv-{pid}/nodes.jsonl"

    n_env, tex_labels = tex_equation_info(srcdir)
    nodes = latest_nodes(db)
    eq_nodes = {nid: r for nid, r in nodes.items() if r.get("equation_labels")}
    db_labels = set(l for r in eq_nodes.values() for l in r["equation_labels"])
    deriv_text = deriv.read_text(encoding="utf-8", errors="replace") if deriv.exists() else ""

    # 1-2. acyclicity + edge closure
    dag_ok, cyc = acyclic(nodes)
    print(f"[{'PASS' if dag_ok else 'FAIL'}] acyclicity over {len(nodes)} nodes" +
          ("" if dag_ok else f" — cycle members: {cyc}"))
    dangling = sorted({p for r in nodes.values() for p in (r.get("predecessors") or [])
                       if p not in nodes})
    closure_ok = not dangling
    print(f"[{'PASS' if closure_ok else 'FAIL'}] edge closure" +
          ("" if closure_ok else f" — dangling predecessors: {dangling}"))

    # 3. label coverage
    miss_db = sorted(tex_labels - db_labels)
    miss_deriv = sorted(l for l in tex_labels if l not in deriv_text)
    cov_db_ok, cov_deriv_ok = not miss_db, not miss_deriv
    print(f"[{'PASS' if cov_db_ok else 'FAIL'}] tex labels covered by DB equation_labels "
          f"({len(tex_labels - set(miss_db))}/{len(tex_labels)})" +
          ("" if cov_db_ok else f" — missing: {miss_db}"))
    print(f"[{'PASS' if cov_deriv_ok else 'FAIL'}] tex labels present in derivation.md" +
          ("" if cov_deriv_ok else f" — missing: {miss_deriv}"))

    # 4. registry / count report
    auto = sorted(l for l in db_labels if l.startswith("eq:auto-"))
    reg_miss = sorted(l for l in db_labels if l not in deriv_text)
    reg_ok = not reg_miss
    print(f"[{'PASS' if reg_ok else 'FAIL'}] every DB label registered in derivation.md" +
          ("" if reg_ok else f" — unregistered: {reg_miss}"))
    print(f"[INFO] tex equation environments: {n_env}; DB labels: {len(db_labels)} "
          f"(tex-born {len(db_labels) - len(auto)}, auto-assigned {len(auto)}); "
          f"equation nodes: {len(eq_nodes)}")
    return ok and dag_ok and closure_ok and cov_db_ok and cov_deriv_ok and reg_ok

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", action="append",
                    default=None, help="arxiv id (repeatable)")
    ap.add_argument("--project", default="z4c-CMM")
    ap.add_argument("--repo-root", default=".")
    a = ap.parse_args()
    papers = a.paper or ["1010.0523v2", "2007.01339", "2308.10361"]
    root = Path(a.repo_root)
    results = [check_paper(root, p, a.project) for p in papers]
    print(f"\nOVERALL: {'PASS' if all(results) else 'FAIL'}")
    sys.exit(0 if all(results) else 1)

if __name__ == "__main__":
    main()
