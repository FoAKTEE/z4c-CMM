# Claims — arXiv:2007.01339

Schema per research_admission_contract.md: working context, claim, evidence type, dependencies, status. All statuses here reflect *our* verification state, not the paper's confidence: transcribed literature content is literature-grounded, NOT a verified derivation.

## C-1 Five-step Bondi-like → Bondi-Sachs transformation
- Working context: asymptotically flat vacuum exterior; Bondi-like metric eq:BondiLikeMetric; assumptions 1–5, 8.
- Claim: the sequence eq:XIfcEom → eq:conformal_factor_determinant → eq:auto-21/eq:uRequation → eq:auto-22/eq:bsxA → eq:Step5RadialCoordinate maps any Bondi-like metric to a true Bondi-Sachs frame, unique up to BMS, with asymptotic results eq:auto-23.
- Evidence type: literature grounding (paper derivation + companion Mathematica notebook, not independently rerun here).
- Dependencies: eq:UpIndexBondiLikeMetric, eq:ifcBeta, eq:ifcW, eq:auto-19, eq:ifc_omega.
- Status: `[PRELIMINARY]` `[UNCHECKED]` — symbolic re-derivation obligation O-2.

## C-2 Partially flat gauge ⇒ no pure-gauge logarithms
- Working context: partially flat gauge Ĵ^(0)=Û^(0)=0; regular initial Ĵ with Ĵ^(2)=0; assumptions 5–7.
- Claim: the regularity conditions eq:QWregularity and eq:Hregularity are satisfied order-by-order for all hypersurface equations, and are preserved in time (eq:auto-36), so all scalars remain polynomial in 1/r for the entire evolution.
- Evidence type: literature grounding (perturbative expansion argument in Section V).
- Dependencies: eq:QWabstract, eq:Habstract, eq:auto-32 … eq:auto-36, eq:SubleadingQRegularity, eq:auto-33.
- Status: `[PRELIMINARY]` `[UNCHECKED]` — order-by-order expansion re-derivable symbolically (O-3).

## C-3 Compactified hypersurface equations are equivalent to the classic system
- Working context: any Bondi-like coordinates; compactified y̆ = 1 − 2R/r; assumptions 3, 12, 13.
- Claim: eq:Betanumeric, eq:Qnumeric, eq:Unumeric, eq:Wnumeric, eq:Hnumeric are a coordinate transformation + simplification of the Bishop et al. (1997) equations, falling into the three normal forms eq:form1–eq:form3.
- Evidence type: literature grounding; paper cites the companion Mathematica package as the rederivation.
- Dependencies: eq:BondiHierarchyBeta, eq:numerical-jacobian-u, eq:numerical-h-conversion, eq:auto-27 … eq:auto-31.
- Status: `[PRELIMINARY]` `[UNCHECKED]` — highest-value verification target (O-1); long source bundles are error-prone.

## C-4 News formulas
- Working context: partially flat gauge at scri+ (C-1, C-2 hold); BMS-frame caveat (assumption 10).
- Claim: the Bondi news equals eq:NewsDefinitionIncompletelyFlat in partially flat quantities and eq:NewsDefinitionBondiLike in an arbitrary Bondi-like gauge, compatible with Bishop et al. but simpler.
- Evidence type: literature grounding.
- Dependencies: eq:NewsDefinitionBondi, eq:auto-23, eq:ifcBeta, eq:ifcH.
- Status: `[PRELIMINARY]` `[UNCHECKED]` (O-4).

## C-5 Closed-form spin coefficients and asymptotic Weyl scalars
- Working context: NP tetrad eq:auto-39; vacuum (assumption 2); peeling (assumption 14); tetrad-convention caveat.
- Claim: spin coefficients take the closed forms eq:auto-43; bulk Ψ0, Ψ1 reduce to eq:auto-48; leading asymptotic Weyl scalars are eq:auto-50 (partially flat) and eq:auto-51 (Bondi-Sachs, related by Boyle 2016 BMS transformations).
- Evidence type: literature grounding (Θ-cancellation check done in the companion package per the paper).
- Dependencies: eq:auto-42, eq:auto-44, eq:auto-46, eq:auto-32 … eq:auto-35, eq:auto-23.
- Status: `[PRELIMINARY]` `[UNCHECKED]` (O-5).

## C-6 The evolution loop is implementable without breaking the hierarchy
- Working context: computational roadmap, Section V.C; assumptions 7, 11, 15.
- Claim: 𝒰 (eq:UHat) satisfies ∂_r𝒰 = ∂_r Û, so U^(0) (eq:auto-37) and the Jacobian evolution (eq:auto-38) can be computed within the hierarchical pass; together with the worldtube transformations eq:ifcBeta, eq:ifcQ, eq:ifcW, eq:ifcH this defines a complete regularity-preserving algorithm.
- Evidence type: literature grounding (algorithmic construction).
- Dependencies: eq:UHat, eq:auto-37, eq:auto-38, eq:ifcQ, eq:auto-52, eq:ifcH.
- Status: `[PRELIMINARY]` `[UNCHECKED]` — numerical demonstration deferred to companion paper `[FUTURE]`.

## C-7 Known transcription hazards in the source
- Working context: comparison of repeated expressions within the tex.
- Claim: eq:KDef prints K = √(1 − JJ̄) while eq:auto-13 and all uses give K = √(1 + JJ̄); v1 equation errors were corrected after publication (acknowledgments). Therefore implementations must validate against the companion package, not the printed text alone.
- Evidence type: counterexample (internal inconsistency observed in transcription).
- Dependencies: eq:KDef, eq:auto-13.
- Status: `[SOLID]` as an observation of the text; consequence handled by O-1.
