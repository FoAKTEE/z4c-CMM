# Logic DAG — arXiv:2007.01339

One node per equation label (tex labels preserved verbatim; unlabeled environments carry `eq:auto-*` ids assigned in document order in `derivation.md`). Each node lists its direct predecessors (equation labels within this paper). Roots have `[]`. This file is the source of truth mirrored by `knowledge-database/paper_arxiv-2007.01339/nodes.jsonl` (task `eq-dag-stage1-2007.01339`). 103 nodes total: 48 tex-born labels + 55 auto labels. The graph is acyclic; cross-section edges (e.g. eq:NewsDefinitionBondiLike → eq:ifcH in Appendix B) follow content, not page order.

## Section II.A — Bondi-Sachs and Bondi-like metrics (sec:bondi_sachs)

| node | predecessors |
|---|---|
| eq:BondiSachsMetric | [] |
| eq:auto-1 | [eq:BondiSachsMetric] |
| eq:auto-2 | [eq:BondiSachsMetric] |
| eq:BondiLikeMetric | [eq:BondiSachsMetric] |
| eq:auto-3 | [eq:BondiLikeMetric] |
| eq:auto-4 | [] |
| eq:ComplexDyadChoice | [eq:auto-4] |
| eq:DyadRotation | [eq:auto-4] |
| eq:auto-5 | [eq:DyadRotation] |
| eq:QDef | [eq:BondiLikeMetric, eq:auto-4] |
| WDef | [eq:BondiLikeMetric] |
| eq:KDef | [eq:BondiLikeMetric, eq:auto-4, eq:auto-1] |
| eq:auto-6 | [eq:auto-4] |
| eq:auto-7 | [eq:auto-6, eq:ComplexDyadChoice] |

## Section II.B — Boundary transformations (sec:boundary_transforms)

| node | predecessors |
|---|---|
| eq:auto-8 | [] |
| eq:auto-9 | [eq:auto-8] |
| eq:auto-10 | [eq:auto-8] |
| eq:auto-11 | [eq:auto-9, eq:auto-10] |
| eq:auto-12 | [eq:auto-8, eq:auto-11] |
| eq:BondiLikeRadius | [eq:auto-12] |
| eq:auto-13 | [eq:BondiLikeRadius, eq:auto-12, eq:QDef, WDef, eq:KDef] |
| eq:auto-14 | [eq:auto-13] |
| eq:auto-15 | [eq:auto-13, eq:auto-14] |
| eq:auto-16 | [eq:auto-13] |

## Section II.C — Hierarchical evolution system (sec:EquationHierarchy)

| node | predecessors |
|---|---|
| eq:BondiHierarchyBeta | [eq:BondiLikeMetric, eq:QDef, WDef, eq:KDef] |

## Section II.D — Observables at scri+ (sec:FoundationsScri)

| node | predecessors |
|---|---|
| eq:BondiNews | [eq:BondiSachsMetric] |
| eq:auto-17 | [eq:BondiNews, eq:auto-4] |
| eq:auto-18 | [eq:BondiLikeMetric] |

## Section III — Bondi-like → Bondi transformations (sec:bondi_transforms)

| node | predecessors |
|---|---|
| eq:UpIndexBondiLikeMetric | [eq:BondiLikeMetric, WDef] |
| eq:AsymptoticExpansionOfU | [eq:BondiLikeMetric] |
| eq:XIfcEom | [eq:auto-18, eq:AsymptoticExpansionOfU] |
| eq:conformal_factor_determinant | [eq:XIfcEom, eq:auto-1] |
| eq:auto-19 | [eq:XIfcEom, eq:auto-4] |
| eq:ifc_omega | [eq:conformal_factor_determinant, eq:auto-19] |
| eq:ifcBeta | [eq:conformal_factor_determinant, eq:auto-19, eq:ifc_omega, eq:UpIndexBondiLikeMetric, eq:AsymptoticExpansionOfU] |
| eq:ifcW | [eq:conformal_factor_determinant, eq:auto-19, eq:ifc_omega, eq:UpIndexBondiLikeMetric, eq:AsymptoticExpansionOfU] |
| eq:auto-20 | [eq:ifc_omega, eq:auto-19, eq:AsymptoticExpansionOfU] |
| eq:auto-21 | [eq:ifcBeta] |
| eq:uRequation | [eq:auto-21, eq:UpIndexBondiLikeMetric, eq:ifcW] |
| eq:auto-22 | [eq:uRequation] |
| eq:bsxA | [eq:auto-22, eq:auto-21, eq:UpIndexBondiLikeMetric] |
| eq:Step5RadialCoordinate | [eq:bsxA, eq:conformal_factor_determinant] |
| eq:auto-23 | [eq:Step5RadialCoordinate, eq:uRequation, eq:bsxA] |

## Section III.B — News inference (sec:news)

| node | predecessors |
|---|---|
| eq:NewsDefinitionBondi | [eq:BondiNews, eq:auto-17] |
| eq:NewsDefinitionIncompletelyFlat | [eq:NewsDefinitionBondi, eq:auto-23] |
| eq:NewsDefinitionBondiLike | [eq:NewsDefinitionIncompletelyFlat, eq:ifcBeta, eq:ifcH, eq:auto-19, eq:auto-20] |

## Section IV — Compactified evolution equations (sec:compactified_evolution)

| node | predecessors |
|---|---|
| eq:auto-25 | [eq:BondiLikeRadius] |
| eq:auto-24 | [eq:BondiLikeMetric, eq:auto-25] |
| eq:numerical-jacobian-u | [eq:auto-24] |
| eq:numerical-h-conversion | [eq:numerical-jacobian-u, eq:auto-16] |
| eq:form1 | [eq:BondiHierarchyBeta, eq:numerical-jacobian-u] |
| eq:form2 | [eq:BondiHierarchyBeta, eq:numerical-jacobian-u] |
| eq:form3 | [eq:BondiHierarchyBeta, eq:numerical-jacobian-u, eq:numerical-h-conversion] |
| eq:auto-26 | [eq:form3] |
| eq:HypersurfaceBeta | [eq:BondiHierarchyBeta, eq:KDef] |
| eq:Betanumeric | [eq:HypersurfaceBeta, eq:numerical-jacobian-u, eq:form1] |
| eq:auto-27 | [eq:KDef, eq:auto-6] |
| eq:HypersurfaceQ | [eq:BondiHierarchyBeta, eq:QDef, eq:auto-27] |
| eq:auto-28 | [eq:auto-27, eq:numerical-jacobian-u] |
| eq:Qnumeric | [eq:HypersurfaceQ, eq:numerical-jacobian-u, eq:form2, eq:auto-28] |
| eq:HypersurfaceU | [eq:BondiHierarchyBeta, eq:QDef] |
| eq:Unumeric | [eq:HypersurfaceU, eq:numerical-jacobian-u, eq:form1] |
| eq:auto-29 | [eq:KDef, eq:auto-6] |
| eq:HypersurfaceW | [eq:BondiHierarchyBeta, WDef, eq:auto-29] |
| eq:Wnumeric | [eq:HypersurfaceW, eq:numerical-jacobian-u, eq:form2] |
| eq:script-b | [eq:KDef, eq:auto-6] |
| eq:HypersurfaceH | [eq:BondiHierarchyBeta, eq:auto-16, eq:script-b] |
| eq:auto-30 | [eq:script-b, eq:numerical-jacobian-u] |
| eq:auto-31 | [eq:script-b, eq:numerical-jacobian-u] |
| eq:Hnumeric | [eq:HypersurfaceH, eq:numerical-jacobian-u, eq:numerical-h-conversion, eq:form3, eq:auto-30, eq:auto-31] |

## Section V — Regularity-preserving CCE (sec:regularity_preservation)

| node | predecessors |
|---|---|
| eq:QWabstract | [eq:form2] |
| eq:QWregularity | [eq:QWabstract] |
| eq:Habstract | [eq:form3] |
| eq:Hregularity | [eq:Habstract] |
| eq:auto-32 | [eq:HypersurfaceBeta] |
| eq:SubleadingQRegularity | [eq:HypersurfaceQ, eq:auto-32, eq:QWregularity] |
| eq:auto-33 | [eq:SubleadingQRegularity, eq:QWregularity] |
| eq:auto-34 | [eq:HypersurfaceU, eq:SubleadingQRegularity, eq:auto-32] |
| eq:auto-35 | [eq:HypersurfaceW, eq:auto-34, eq:auto-32, eq:QWregularity] |
| eq:auto-36 | [eq:HypersurfaceH, eq:Hregularity, eq:auto-32, eq:auto-33, eq:auto-34, eq:auto-35] |
| eq:UHat | [eq:auto-19, eq:ifc_omega, eq:Unumeric, eq:ifcBeta, eq:ifcW] |
| eq:auto-37 | [eq:UHat] |
| eq:auto-38 | [eq:XIfcEom, eq:auto-37] |

## Section VI — Newman-Penrose Weyl scalars (sec:weyl_scalars)

| node | predecessors |
|---|---|
| eq:auto-39 | [eq:BondiLikeMetric, eq:KDef] |
| eq:auto-40 | [eq:auto-39] |
| eq:auto-41 | [eq:auto-40] |
| eq:NPCovariantDerivatives | [eq:auto-39] |
| eq:auto-42 | [eq:auto-39, eq:NPCovariantDerivatives] |
| eq:auto-43 | [eq:auto-42, eq:auto-39, eq:QDef, WDef, eq:KDef] |
| eq:auto-44 | [eq:auto-42, eq:NPCovariantDerivatives] |
| eq:auto-45 | [eq:NPCovariantDerivatives, eq:auto-39] |
| eq:auto-46 | [eq:auto-45, eq:auto-43] |
| eq:auto-47 | [eq:auto-44, eq:auto-46, eq:auto-43] |
| eq:auto-48 | [eq:auto-44, eq:auto-43] |
| eq:auto-49 | [] |
| eq:auto-50 | [eq:auto-47, eq:auto-48, eq:auto-49, eq:auto-32, eq:auto-33, eq:auto-34, eq:auto-35] |
| eq:auto-51 | [eq:auto-50, eq:auto-23] |

## Appendix B — Extra partially flat transformations (sec:IF_extras)

| node | predecessors |
|---|---|
| eq:auto-52 | [eq:auto-19, eq:ifc_omega, eq:ifcBeta, eq:ifcW] |
| eq:ifcQ | [eq:QDef, eq:auto-52, eq:ifcBeta] |
| eq:ifcH | [eq:auto-16, eq:auto-19, eq:auto-20, eq:ifc_omega, eq:ifcW] |

## Appendix C — Perturbative expansions near scri+ (app:perturbative_transformations)

| node | predecessors |
|---|---|
| eq:auto-53 | [eq:uRequation] |
| eq:auto-54 | [eq:bsxA, eq:auto-53] |
| eq:auto-55 | [eq:Step5RadialCoordinate, eq:auto-53, eq:auto-54] |

## Implementation clusters (consumed by implementation_plan_python.md)

- C1 foundations/dyad: eq:auto-4 … eq:auto-7, eq:QDef, WDef, eq:KDef.
- C2 worldtube boundary: eq:auto-8 … eq:auto-16, eq:BondiLikeRadius.
- C3 gauge transformation: eq:UpIndexBondiLikeMetric … eq:auto-23, eq:ifcQ, eq:auto-52, eq:ifcH, eq:auto-53 … eq:auto-55.
- C4 compactified hypersurface solvers: eq:auto-24 … eq:auto-31, eq:form1-3, eq:Hypersurface*/\*numeric.
- C5 regularity/evolution loop: eq:QWabstract … eq:auto-38.
- C6 news and Weyl outputs: eq:BondiNews, eq:auto-17, news definitions, eq:auto-39 … eq:auto-51.
