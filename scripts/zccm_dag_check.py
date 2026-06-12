#!/usr/bin/env python3
"""zccm_dag_check.py — verifier for the Z4c-CCM synthesis DAG.

Checks, against reformulate/z4c-CMM/synthesis_z4c-CCM/formulation_dag.md:
  1. MERMAID GRAPH — every `X --> Y` edge references a declared node id and the
     graph is acyclic (Kahn).
  2. SOURCE GROUNDING — every equation label cited in the "Node ledger" table
     exists as a node_id in the named paper's knowledge ledger
     (knowledge-database/paper_arxiv-<id>/nodes.jsonl).
  3. N-NODE LEDGER SYNC — every N-node row in the table has a row in
     knowledge-database/paper_z4c-CCM/nodes.jsonl whose status agrees with the
     table's marker ([PRELIMINARY]→preliminary, [HYPOTHESIS]→hypothesis,
     [FUTURE]→future, [SOLID]→solid).

Exit 0 only if all checks pass.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "reformulate/z4c-CMM/synthesis_z4c-CCM/formulation_dag.md"
PAPER_IDS = {"1010.0523v2", "2007.01339", "2308.10361"}
MARKER2STATUS = {"SOLID": "solid", "PRELIMINARY": "preliminary",
                 "HYPOTHESIS": "hypothesis", "FUTURE": "future"}

def latest_nodes(jsonl: Path):
    rows = {}
    if jsonl.exists():
        for line in jsonl.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                if r.get("status") != "amended":
                    rows[r["node_id"]] = r
    return rows

def label_index(rows):
    """All resolvable labels: node_ids plus every equation_labels entry."""
    idx = set(rows)
    for r in rows.values():
        idx.update(r.get("equation_labels") or [])
    return idx

def main():
    ok = True
    text = DOC.read_text(encoding="utf-8")
    mermaid = re.search(r"```mermaid\n(.*?)```", text, re.S)
    if not mermaid:
        print("[FAIL] no mermaid block found"); sys.exit(1)
    block = mermaid.group(1)

    declared = set(re.findall(r'^\s*([A-Z]\d+)\["', block, re.M))
    edges = re.findall(r"^\s*([A-Z]\d+)\s*-->\s*([A-Z]\d+)\s*$", block, re.M)
    undeclared = sorted({n for e in edges for n in e if n not in declared})
    print(f"[{'PASS' if not undeclared else 'FAIL'}] edge endpoints declared "
          f"({len(declared)} nodes, {len(edges)} edges)" +
          ("" if not undeclared else f" — undeclared: {undeclared}"))
    ok &= not undeclared

    indeg = {n: 0 for n in declared}
    succ = {n: [] for n in declared}
    for a, b in edges:
        if a in declared and b in declared:
            indeg[b] += 1
            succ[a].append(b)
    queue = [n for n, d in indeg.items() if d == 0]
    seen = 0
    while queue:
        n = queue.pop()
        seen += 1
        for s in succ[n]:
            indeg[s] -= 1
            if indeg[s] == 0:
                queue.append(s)
    acyclic = seen == len(declared)
    print(f"[{'PASS' if acyclic else 'FAIL'}] mermaid graph acyclic" +
          ("" if acyclic else f" — cycle members: {sorted(n for n, d in indeg.items() if d > 0)}"))
    ok &= acyclic

    # ONE DAG: the graph must be a single weakly-connected component.
    adj = {n: set() for n in declared}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    comp, stack = set(), [next(iter(declared))]
    while stack:
        n = stack.pop()
        if n not in comp:
            comp.add(n)
            stack.extend(adj[n] - comp)
    one_dag = comp == declared
    print(f"[{'PASS' if one_dag else 'FAIL'}] single connected DAG (ONE graph)" +
          ("" if one_dag else f" — disconnected: {sorted(declared - comp)}"))
    ok &= one_dag

    ledgers = {pid: label_index(latest_nodes(ROOT / f"knowledge-database/paper_arxiv-{pid}/nodes.jsonl"))
               for pid in PAPER_IDS}
    nrows = latest_nodes(ROOT / "knowledge-database/paper_z4c-CCM/nodes.jsonl")

    table = re.search(r"\| node \| paper \|.*?\n((?:\|.*\n)+)", text)
    rows = [[c.strip() for c in ln.strip().strip("|").split("|")]
            for ln in table.group(1).splitlines() if ln.startswith("|")]
    rows = [r for r in rows if len(r) == 4 and not set(r[0]) <= {"-"}]

    miss_labels, miss_n, status_bad = [], [], []
    for node, paper, labels, status in rows:
        if node.startswith("N"):
            m = re.search(r"\[(\w+)\]", status)
            want = MARKER2STATUS.get(m.group(1), "?") if m else "?"
            row = nrows.get(node)
            if row is None:
                miss_n.append(node)
            elif row.get("status") != want:
                status_bad.append(f"{node}: ledger={row.get('status')} table={want}")
            continue
        pids = [p for p in PAPER_IDS if p in paper]
        if not pids or labels.startswith("("):
            continue
        for lab in re.split(r"[;,]", labels):
            lab = re.sub(r"\((P\d|[^)]*)\)", "", lab).strip()
            if not lab or lab.startswith("new"):
                continue
            if not any(lab in ledgers[p] for p in pids):
                miss_labels.append(f"{node}: {lab}")
    print(f"[{'PASS' if not miss_labels else 'FAIL'}] source labels resolve in paper ledgers" +
          ("" if not miss_labels else f" — missing: {miss_labels}"))
    print(f"[{'PASS' if not miss_n else 'FAIL'}] N-nodes present in z4c-CCM ledger" +
          ("" if not miss_n else f" — missing: {miss_n}"))
    print(f"[{'PASS' if not status_bad else 'FAIL'}] N-node statuses consistent" +
          ("" if not status_bad else f" — mismatched: {status_bad}"))
    ok &= not miss_labels and not miss_n and not status_bad

    print(f"\nOVERALL: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
