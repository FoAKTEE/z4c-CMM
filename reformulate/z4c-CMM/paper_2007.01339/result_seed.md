# Result seed — arXiv:2007.01339 (stage-1 initial result-log entries)

Markers per `phys-agentic-loop/_common/contracts/markers.md`; statuses per the admission contract. Evidence for every entry below is *transcription from the source tex* — literature grounding, never promoted to checked derivation.

## RS-1 — Equation corpus transcribed
- name: full equation decomposition of arXiv:2007.01339
- working_context: stage-1 decomposition; tex `ref-paper/arxiv-2007.01339/src/characteristic_formulation.tex`
- claim: all 100 equation environments (48 tex labels + 55 auto labels) are transcribed in `derivation.md` with dependencies mirrored in `logic.md` and the knowledge DB (task `eq-dag-stage1-2007.01339`).
- evidence_type: literature grounding / transcription; verifier: `scripts/eqdag_check.py --paper 2007.01339` (output pasted in stage-1 report).
- status: `checked` (as a transcription artifact, not as physics) — `[SOLID]` for coverage, conditional on verifier PASS.

## RS-2 — Five-step gauge transformation (claim C-1)
- claim: explicit Bondi-like → Bondi-Sachs transformation with partially flat intermediate gauge (eq:XIfcEom … eq:auto-23).
- evidence_type: literature grounding.
- dependencies: C-1 nodes; assumptions 1–5, 8.
- status: `unchecked` — `[PRELIMINARY]` `[UNCHECKED]`; closes via O-2.

## RS-3 — Logarithm-free evolution in partially flat gauge (claim C-2)
- claim: partially flat gauge + Ĵ^(2)=0 initial data ⇒ regular (polynomial in 1/r) evolution (eq:QWregularity, eq:Hregularity, eq:auto-32 … eq:auto-36).
- evidence_type: literature grounding (paper's perturbative proof sketch).
- status: `conditional` on assumptions 5–7; `[PRELIMINARY]` `[UNCHECKED]`; closes via O-3.

## RS-4 — Compactified hypersurface equations (claim C-3)
- claim: eq:Betanumeric / eq:Qnumeric / eq:Unumeric / eq:Wnumeric / eq:Hnumeric equivalent to the classic characteristic system.
- evidence_type: literature grounding; companion Mathematica package cited as the paper's own rederivation (unobtained — O-7).
- status: `unchecked` — `[UNCHECKED]` `[BLOCKING]` for downstream numerical use until O-1 closes.

## RS-5 — News formulas (claim C-4)
- claim: eq:NewsDefinitionIncompletelyFlat (partially flat) and eq:NewsDefinitionBondiLike (arbitrary Bondi-like).
- evidence_type: literature grounding.
- status: `unchecked` — `[PRELIMINARY]` `[UNCHECKED]`; closes via O-4.

## RS-6 — Spin coefficients and asymptotic Weyl scalars (claim C-5)
- claim: eq:auto-43 (spin coefficients), eq:auto-48 (bulk Ψ0, Ψ1), eq:auto-50/eq:auto-51 (asymptotic Weyl scalars, PF and Bondi gauges).
- evidence_type: literature grounding.
- status: `unchecked` — `[PRELIMINARY]` `[UNCHECKED]`; closes via O-5.

## RS-7 — Implementable evolution roadmap (claim C-6)
- claim: 𝒰-trick (eq:UHat, eq:auto-37, eq:auto-38) preserves the hierarchical structure while imposing partially flat gauge.
- evidence_type: literature grounding (algorithm construction); numerical demonstration `[FUTURE]` (companion paper).
- status: `unchecked` — `[PRELIMINARY]`.

## RS-8 — Source-text inconsistency (claim C-7)
- claim: eq:KDef prints K = √(1−JJ̄), inconsistent with K = √(1+JJ̄) used in eq:auto-13 and throughout; v1 equation errors acknowledged by the authors.
- evidence_type: counterexample (internal textual inconsistency).
- status: `checked` (observation) — `[SOLID]`; consequence: `[OPEN]` O-6 reconciliation before implementation.

## Open-obligation visibility
All `[OPEN]`/`[BLOCKING]`/`[FUTURE]` items live in `obligations.md` (O-1 … O-11) and must remain visible in any downstream report rendered from this log.
