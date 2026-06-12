# logic.md — equation dependency DAG, arXiv:1010.0523v2

Source of truth mirrored by knowledge-database/paper_arxiv-1010.0523v2/nodes.jsonl
(task_id eq-dag-stage1-1010.0523v2). One node per equation; node id = equation label
(tex labels verbatim, incl. `sol:lf`, `sol:lf2`, `eq:theta _modes`; 36 auto ids
eq:auto-1..eq:auto-36 for unlabeled environments). `<- []` marks root nodes.
Graph is acyclic (verified by scripts/eqdag_check.py: OVERALL PASS).
All edges are literature-grounded structural dependencies (which earlier equations or
definitions each equation is derived from or uses), not independently re-derived.

## Section II — The Z4c formulation

- `eq:Z4_ADM_1` <- []
- `eq:Theta_dot` <- ['eq:auto-1']
- `eq:Z4_ADM_2` <- ['eq:auto-1']
- `eq:auto-1` <- []
- `eq:ham-Z4` <- ['eq:Theta_dot', 'eq:Z4_ADM_1', 'eq:Z4_ADM_2', 'eq:auto-1', 'eq:auto-2']
- `eq:mom-Z4` <- ['eq:Theta_dot', 'eq:Z4_ADM_1', 'eq:Z4_ADM_2', 'eq:auto-1', 'eq:auto-2']
- `eq:auto-2` <- []
- `eq:auto-3` <- ['eq:Theta_dot', 'eq:Z4_ADM_1']
- `eq:auto-4` <- ['eq:Z4_ADM_2', 'eq:auto-3']
- `eq:auto-5` <- ['eq:auto-3']
- `eq:Z4_decomp_first` <- ['eq:Z4_ADM_1', 'eq:auto-3']
- `eq:auto-6` <- ['eq:Z4_ADM_1', 'eq:auto-3', 'eq:auto-8']
- `eq:auto-7` <- ['eq:Z4_ADM_1', 'eq:Z4_ADM_2', 'eq:auto-3', 'eq:auto-4', 'eq:auto-5']
- `eq:auto-8` <- ['eq:auto-3', 'eq:auto-4', 'eq:auto-5']
- `eq:Conf_Constr_1` <- ['eq:auto-3', 'eq:auto-4']
- `eq:Conf_Constr_2` <- ['eq:auto-1', 'eq:auto-3']
- `eq:punc_alpha` <- ['eq:auto-3']
- `eq:punc_beta` <- ['eq:auto-3', 'eq:auto-4']
- `eq:sec-lapse` <- ['eq:Z4_decomp_first', 'eq:punc_alpha', 'eq:auto-2']
- `eq:sec-shift` <- ['eq:punc_alpha', 'eq:punc_beta', 'eq:Z4_ADM_1', 'eq:auto-2']
- `eq:sec-gamma` <- ['eq:punc_alpha', 'eq:punc_beta', 'eq:Z4_ADM_1', 'eq:auto-2']
- `eq:def_theta` <- ['eq:punc_alpha', 'eq:Z4_ADM_1', 'eq:auto-2']
- `eq:def_Zi` <- ['eq:punc_beta', 'eq:auto-2', 'eq:auto-5']
- `eq:sys-theta-Z` <- ['eq:Theta_dot', 'eq:Z4_ADM_2']
- `eq:sys-ham-mom` <- ['eq:ham-Z4', 'eq:mom-Z4']
- `eq:auto-9` <- []
- `eq:auto-10` <- ['eq:auto-9']
- `eq:auto-11` <- ['eq:sec-lapse', 'eq:sec-shift', 'eq:sec-gamma', 'eq:auto-12', 'eq:auto-9']
- `eq:auto-12` <- ['eq:auto-9']
- `eq:auto-13` <- ['eq:sec-shift', 'eq:sec-gamma', 'eq:auto-9']
- `eq:auto-14` <- ['eq:sec-gamma', 'eq:auto-9']
- `eq:auto-15` <- ['eq:sys-theta-Z', 'eq:sys-ham-mom']
- `eq:backg-metric` <- []
- `eq:auto-16` <- ['eq:backg-metric']
- `eq:auto-17` <- ['eq:backg-metric']
- `eq:nullvector-k` <- ['eq:sec-shift', 'eq:auto-17']
- `eq:nullvector-m` <- ['eq:sec-lapse', 'eq:sec-shift', 'eq:auto-17']
- `eq:general_CPBCs` <- ['eq:nullvector-k', 'eq:sys-theta-Z', 'eq:auto-16']
- `eq:BCs-alpha` <- ['eq:nullvector-m', 'eq:sec-lapse']
- `eq:general_BCs_gauge_first` <- ['eq:nullvector-k', 'eq:sec-shift', 'eq:auto-9']
- `eq:BCs_lastII` <- ['eq:nullvector-k', 'eq:nullvector-m', 'eq:sec-gamma', 'eq:auto-9']

## Section III — Analytical setup

- `eq:backg-metricP` <- ['eq:backg-metric']
- `eq:auto-18` <- ['eq:backg-metricP', 'eq:auto-17']
- `eq:full-sec-z4-alpha` <- ['eq:backg-metricP', 'eq:sec-lapse']
- `eq:full-sec-z4-beta` <- ['eq:backg-metricP', 'eq:sec-shift']
- `eq:full-sec-z4-metric` <- ['eq:backg-metricP', 'eq:sec-gamma']
- `eq:2+1Z4_a` <- ['eq:full-sec-z4-alpha', 'eq:auto-18', 'eq:auto-9']
- `eq:2+1Z4_ac` <- ['eq:backg-metricP', 'eq:BCs-alpha']
- `eq:Z4_betas_2+1` <- ['eq:full-sec-z4-beta', 'eq:auto-18', 'eq:auto-9']
- `eq:Z4_betas_2+1BC` <- ['eq:general_BCs_gauge_first', 'eq:backg-metricP', 'eq:BCs_lastII']
- `eq:Z4_gammaAB_2+1` <- ['eq:full-sec-z4-metric', 'eq:auto-18', 'eq:auto-9']
- `eq:Z4_gammaAB_2+1BC` <- ['eq:backg-metricP', 'eq:BCs_lastII']
- `eq:EOM-Lambda` <- ['eq:full-sec-z4-metric', 'eq:auto-12', 'eq:auto-18']
- `eq:BC-Lambda` <- ['eq:general_CPBCs', 'eq:backg-metricP', 'eq:def_theta', 'eq:auto-12']
- `eq:eomgammasA` <- ['eq:full-sec-z4-metric', 'eq:auto-18', 'eq:auto-9']
- `eq:Z4_gammasbc_2+1` <- ['eq:general_CPBCs', 'eq:backg-metricP', 'eq:auto-12', 'eq:def_Zi']
- `eq:Z4_gammaAbc_2+1` <- ['eq:general_CPBCs', 'eq:backg-metricP', 'eq:auto-12', 'eq:def_Zi']

## Section IV — Well-posedness results

- `eq:auto-19` <- ['eq:backg-metricP', 'eq:generalform-U', 'eq:sys-theta-Z']
- `eq:MLF-Mat` <- ['eq:backg-metricP', 'eq:generalform-U', 'eq:sys-theta-Z', 'eq:auto-19']
- `eq:auto-20` <- ['eq:MLF-Mat']
- `eq:general_sol_theta` <- ['eq:solution_lap-four', 'eq:auto-20']
- `eq:fcpbcTheta` <- ['eq:general_sol_theta', 'eq:boundary_stable', 'eq:general_CPBCs']
- `eq:esTheta` <- ['eq:fcpbcTheta', 'eq:estWoutS']
- `eq:hoTheta-bc` <- ['eq:general_CPBCs', 'eq:auto-19', 'eq:MLF-Mat', 'eq:h-bc']
- `eq:auto-21` <- ['eq:general_sol_theta', 'eq:hoTheta-bc', 'eq:bc_est-bc']
- `eq:auto-22` <- ['eq:generalform-U', 'eq:Z4_betas_2+1', 'eq:2+1Z4_a']
- `eq:auto-23` <- ['eq:generalform-U-BC', 'eq:Z4_betas_2+1BC', 'eq:Z4_betas_2+1', 'eq:2+1Z4_ac', 'eq:2+1Z4_a', 'eq:auto-22']
- `eq:auto-24` <- ['eq:backg-metricP']
- `eq:genralsol-gauge` <- ['eq:solution_lap-four', 'eq:auto-23']
- `eq:auto-25` <- ['eq:genralsol-gauge', 'eq:auto-23']
- `eq:auto-26` <- ['eq:full-sec-z4-alpha', 'eq:full-sec-z4-beta', 'eq:high-orderwave', 'eq:auto-23']
- `eq:auto-27` <- ['eq:auto-23', 'eq:auto-26']
- `eq:auto-28` <- ['eq:auto-23', 'eq:auto-27']
- `eq:auto-29` <- ['eq:genralsol-gauge', 'eq:auto-27', 'eq:auto-28']

## Section V — Numerical applications

- `eq:Cmonitor` <- ['eq:Conf_Constr_1', 'eq:Conf_Constr_2', 'eq:auto-1']
- `eq:norm2` <- ['eq:Cmonitor']
- `eq:theta _modes` <- ['eq:auto-15']
- `eq:refcoef` <- ['eq:theta _modes']
- `eq:auto-30` <- []

## Appendix A — Well-posed problems (Kreiss theory)

- `eq:def-system` <- []
- `eq:formA` <- ['eq:def-system']
- `eq:BC` <- ['eq:def-system', 'eq:formA']
- `eq:lapl-four` <- ['eq:def-system']
- `eq:BC-lapl-four` <- ['eq:lapl-four', 'eq:BC']
- `eq:auto-31` <- ['eq:def-system', 'eq:lapl-four']
- `eq:solution_lap-four` <- ['eq:lapl-four', 'eq:auto-31']
- `eq:bc_general` <- ['eq:solution_lap-four', 'eq:BC-lapl-four']
- `sol:lf` <- ['eq:solution_lap-four', 'eq:bc_general']
- `sol:lf2` <- ['eq:lapl-four', 'sol:lf']
- `eq:boundary_stable` <- ['eq:solution_lap-four', 'eq:bc_general']
- `eq:estWoutS` <- ['eq:boundary_stable']
- `eq:gen_estimate` <- ['eq:boundary_stable', 'eq:estWoutS']

## Appendix B — Toy model (shifted wave equation with high order BCs)

- `eq:eqMU` <- ['eq:auto-2']
- `eq:bc_phi` <- ['eq:eqMU']
- `eq:LF-u` <- ['eq:backg-metricP', 'eq:eqMU']
- `eq:BC-u` <- ['eq:backg-metricP', 'eq:bc_phi']
- `eq:auto-32` <- ['eq:LF-u']
- `eq:generalform-U` <- ['eq:auto-32', 'eq:LF-u']
- `eq:generalform-U-BC` <- ['eq:auto-32', 'eq:BC-u']
- `eq:LF-Mat` <- ['eq:generalform-U-BC', 'eq:generalform-U', 'eq:auto-32']
- `eq:general_sol_phi` <- ['eq:solution_lap-four', 'eq:generalform-U', 'eq:LF-Mat']
- `eq:bc_est-phi` <- ['eq:general_sol_phi', 'eq:LF-Mat']
- `eq:bc_est-bc` <- ['eq:bc_est-phi']
- `eq:Bon-stab-phi` <- ['eq:general_sol_phi', 'eq:bc_est-bc']
- `eq:phi-estimate` <- ['eq:gen_estimate', 'eq:Bon-stab-phi']
- `eq:bc_phi-ho` <- ['eq:bc_phi']
- `eq:high-orderwave` <- ['eq:bc_phi-ho', 'eq:LF-u']
- `eq:auto-33` <- ['eq:high-orderwave', 'eq:generalform-U']
- `eq:A-BC` <- ['eq:auto-33']
- `eq:h-bc` <- ['eq:high-orderwave', 'eq:A-BC']
- `eq:full-estimate-W` <- ['eq:general_sol_phi', 'eq:Bon-stab-phi', 'eq:h-bc']

## Appendix C — Implementation of BCs in spherical symmetry

- `eq:auto-34` <- ['eq:auto-3']
- `eq:auto-35` <- ['eq:Conf_Constr_2', 'eq:auto-34']
- `eq:sph_bc_1st` <- ['eq:general_CPBCs', 'eq:Theta_dot', 'eq:auto-34']
- `eq:sph_bc_last` <- ['eq:general_BCs_gauge_first', 'eq:general_CPBCs', 'eq:sph_bc_1st', 'eq:BCs-alpha', 'eq:auto-34', 'eq:auto-35']
- `eq:auto-36` <- []

## Node clusters (implementation partition, see implementation_plan_python.md)

- C1 Z4/Z4c system & constraints: eq:Z4_ADM_1, eq:Theta_dot, eq:Z4_ADM_2, eq:auto-1, eq:ham-Z4, eq:mom-Z4, eq:auto-2, eq:auto-3, eq:auto-4, eq:auto-5, eq:Z4_decomp_first, eq:auto-6, eq:auto-7, eq:auto-8, eq:Conf_Constr_1, eq:Conf_Constr_2
- C2 gauge & second order principal part: eq:punc_alpha, eq:punc_beta, eq:sec-lapse, eq:sec-shift, eq:sec-gamma, eq:def_theta, eq:def_Zi, eq:sys-theta-Z, eq:sys-ham-mom
- C3 2+1 decomposition & characteristic variables: eq:auto-9, eq:auto-10, eq:auto-11, eq:auto-12, eq:auto-13, eq:auto-14, eq:auto-15
- C4 boundary geometry & CPBCs: eq:backg-metric, eq:auto-16, eq:auto-17, eq:nullvector-k, eq:nullvector-m, eq:general_CPBCs, eq:BCs-alpha, eq:general_BCs_gauge_first, eq:BCs_lastII
- C5 frozen-coefficient reduction: eq:backg-metricP, eq:auto-18, eq:full-sec-z4-alpha, eq:full-sec-z4-beta, eq:full-sec-z4-metric, eq:2+1Z4_a, eq:2+1Z4_ac, eq:Z4_betas_2+1, eq:Z4_betas_2+1BC, eq:Z4_gammaAB_2+1, eq:Z4_gammaAB_2+1BC, eq:EOM-Lambda, eq:BC-Lambda, eq:eomgammasA, eq:Z4_gammasbc_2+1, eq:Z4_gammaAbc_2+1
- C6 Kreiss theory (Appendix A): eq:def-system, eq:formA, eq:BC, eq:lapl-four, eq:BC-lapl-four, eq:auto-31, eq:solution_lap-four, eq:bc_general, sol:lf, sol:lf2, eq:boundary_stable, eq:estWoutS, eq:gen_estimate
- C7 toy model (Appendix B): eq:eqMU, eq:bc_phi, eq:LF-u, eq:BC-u, eq:auto-32, eq:generalform-U, eq:generalform-U-BC, eq:LF-Mat, eq:general_sol_phi, eq:bc_est-phi, eq:bc_est-bc, eq:Bon-stab-phi, eq:phi-estimate, eq:bc_phi-ho, eq:high-orderwave, eq:auto-33, eq:A-BC, eq:h-bc, eq:full-estimate-W
- C8 constraint-subsystem well-posedness: eq:auto-19, eq:MLF-Mat, eq:auto-20, eq:general_sol_theta, eq:fcpbcTheta, eq:esTheta, eq:hoTheta-bc, eq:auto-21
- C9 gauge-subsystem well-posedness: eq:auto-22, eq:auto-23, eq:auto-24, eq:genralsol-gauge, eq:auto-25, eq:auto-26, eq:auto-27, eq:auto-28, eq:auto-29
- C10 numerical diagnostics & spherical implementation: eq:Cmonitor, eq:norm2, eq:theta _modes, eq:refcoef, eq:auto-30, eq:auto-34, eq:auto-35, eq:sph_bc_1st, eq:sph_bc_last, eq:auto-36
