# implementation_plan_python.md — arXiv:1010.0523v2

Code/work partition per logic.md node cluster (C1–C10). Target: Python (sympy + numpy/scipy).
Each task consumes the node ids listed in logic.md and produces a verifier-checkable artifact.
All tasks start from literature-grounded transcriptions; passing the stated check promotes the
corresponding claims in claims.md.

## T1 (C1) Z4c symbolic system — sympy
Implement ADM/Z4 variables and the conformal map (eq:auto-3, eq:auto-4, eq:auto-5).
Checks: (a) conformal evolution equations (eq:Z4_decomp_first, eq:auto-6, eq:auto-7) follow
from eq:Z4_ADM_1–eq:Z4_ADM_2 by direct substitution; (b) Ricci split eq:auto-8 reproduces
R_ij of gamma_ij = chi^{-1} gamma-tilde_ij; (c) constraint propagation eq:ham-Z4, eq:mom-Z4
in the principal part; (d) eq:Conf_Constr_1/eq:Conf_Constr_2 equal eq:auto-1 under the map.

## T2 (C2) Principal symbol & second order form — sympy
From eq:punc_alpha, eq:punc_beta and the C1 system derive the fully second order principal
part; verify eq:sec-lapse, eq:sec-shift, eq:sec-gamma and gauge definitions eq:def_theta,
eq:def_Zi; verify Box Theta ≃ 0, Box Z_i ≃ 0, Box H ≃ 0, Box M_i ≃ 0 (eq:sys-theta-Z,
eq:sys-ham-mom).

## T3 (C3) Characteristic analysis — sympy linear algebra
Build the 2+1 projected principal symbol (eq:auto-9, eq:auto-10); verify the characteristic
variables and speeds in scalar (eq:auto-11, eq:auto-12), vector (eq:auto-13), tensor
(eq:auto-14) and constraint (eq:auto-15) sectors by diagonalizing the symbol; confirm weak
hyperbolicity when lambda = 0 and count incoming modes (assumption 8).

## T4 (C4) Boundary geometry & CPBC operators — sympy
Symbolic background metric (eq:backg-metric, eq:auto-16, eq:auto-17), characteristic vectors
(eq:nullvector-k, eq:nullvector-m); implement the boundary operators of eq:general_CPBCs and
eq:BCs-alpha–eq:BCs_lastII as composable differential operators; check L = 0 reduces to the
Bona et al. conditions.

## T5 (C5) Frozen coefficient reduction — sympy
Apply eq:backg-metricP, mu_L = 2/alpha, mu_S = 1 to T2 output; verify
eq:full-sec-z4-alpha/-beta/-metric, the 2+1 wave problems eq:2+1Z4_a–eq:Z4_gammaAbc_2+1, the
cascade property at (eps_alpha, eps_chi) = (1, 1/2), and the bulk coupling at (0,0).
Also check the sign/consistency of s-ring in eq:auto-18 ([OPEN] O3 in obligations.md).

## T6 (C6+C7) Kreiss machinery & toy model — sympy + mpmath
Generic first order Laplace–Fourier pipeline (eq:def-system–eq:gen_estimate): companion
matrix, decaying eigenmodes, boundary determinant D(s,omega), determinant-condition scanner
over Re(s) > 0. Toy model: verify eq:LF-Mat eigenstructure, the bound eq:bc_est-bc by
numerical minimization over (s', omega'), the iteration eq:auto-33 → eq:A-BC → eq:h-bc, and
boundary stability eq:Bon-stab-phi for m = 0,1,2.

## T7 (C8) Constraint subsystem boundary stability — sympy + numpy
Verify M (eq:MLF-Mat) eigenpairs eq:auto-20, the solution eq:general_sol_theta, sigma = q-tilde
(L=0) and a_+^L sigma = q-tilde (eq:auto-21); numerically certify |a_+| ≥ delta > 0 on a
(s', omega') grid with Re(s') > 0 ⇒ estimates eq:fcpbcTheta, eq:esTheta for L = 0,1,2.

## T8 (C9) Gauge subsystem boundary stability — sympy + numpy
Verify the 4x4 system eq:auto-23 (re-derive M, L from T5 wave problems), eigen-decomposition,
integration constants eq:auto-25 / eq:auto-29 including reconstruction of the unspecified
shorthands F, G, H, J of beta-ring ([OPEN] O4); certify boundary stability over a beta-ring,
(s', omega') grid.

## T9 (C10a) Spherical Z4c evolution code — numpy
1D spherical reduction (eq:auto-34, eq:auto-35), method-of-lines RK4, 4th order FD,
Kreiss–Oliger dissipation, ghostzone extrapolation eq:auto-36; implement Sommerfeld, 1st and
2nd order CPBCs as boundary replacement equations (eq:sph_bc_1st, eq:sph_bc_last); also
re-derive eq:sph_bc_1st–eq:sph_bc_last from T4 operators linearized on flat space (check).

## T10 (C10b) Diagnostics & reproduction runs — numpy/scipy
Constraint monitor eq:Cmonitor, norm eq:norm2, characteristic fields eq:theta _modes,
reflection coefficient eq:refcoef (FFT at boundary). Reproduce the perturbed-flat test
(convergence factor, R(k) curves), star and black hole tests qualitatively (puncture initial
data; Kerr–Schild eq:auto-30 with excision optional, [FUTURE]).

Dependencies: T1→T2→{T3,T5}; T4→T5; T6→{T7,T8}; T5→{T7,T8}; {T4,T9}→T10.
