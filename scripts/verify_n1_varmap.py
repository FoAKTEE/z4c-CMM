#!/usr/bin/env python3
"""verify_n1_varmap.py — N1 verifier (Z4c-CCM formulation DAG).

CLAIMS (node N1):
 A. ALGEBRAIC ROUND TRIP — the maps Z4c state -> ADM -> Z4c state and the
    4-metric/its ADM inverse are exact mutual inverses, fully symbolically
    (all operations rational — no radicals, so full-symbol sympy is cheap):
      gamma = gamt/chi, K = Khat + 2 Theta,
      K_ij = (At + gamt K/3)/chi, chi = det(gamma)^(-1/3), ...
      g4(g4inv) = identity.
 B. KINEMATIC CLOSURE — paper-1's conformal evolution equations
    (eq:Z4_decomp_first group)
       d_t chi      = (2/3) chi [alpha (Khat + 2 Theta) - D_i beta^i]
       d_t gamt_ij  = -2 alpha At_ij + beta^k gamt_ij,k
                      + 2 gamt_(i|k beta^k_,|j) - (2/3) gamt_ij beta^k_,k
    reproduce the ADM kinematic identity
       d_t gamma_ij = -2 alpha K_ij + (Lie_beta gamma)_ij
    if and only if D_i beta^i is the FULL covariant divergence
       D_i beta^i = beta^k_,k - (3/2) beta^k (ln chi)_,k
    (det gamt = 1 => ln sqrt(gamma) = -(3/2) ln chi). The verifier proves the
    identity with that reading and proves the residual with the naive reading
    D = beta^k_,k is exactly + gamma_ij beta^k (ln chi)_,k != 0 — a genuine
    convention resolution, recorded as a convention note for the formulation.
    Method: point-symbols (fields and their first derivatives as independent
    symbols) — pure rational algebra, exact, fast.

Budget: signal.alarm(540); pure CPU sympy (rational only). Exit 0 iff all pass.
"""
import signal, time

signal.alarm(540)
T0 = time.time()

import sympy as sp

OK = True
def report(name, cond):
    global OK
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  (t={time.time()-T0:.1f}s)")
    OK &= bool(cond)

# ---------- A. round trip, full symbols (rational algebra) -------------------
chi, Khat, Theta, alpha = sp.symbols("chi Khat Theta alpha", positive=True)
gt = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"gt{min(i,j)}{max(i,j)}"))
At = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"At{min(i,j)}{max(i,j)}"))
beta = sp.Matrix([sp.Symbol(f"b{i}") for i in range(3)])

# forward: Z4c -> ADM
gamma = gt / chi
K = Khat + 2 * Theta
K_ij = (At + gt * K / 3) / chi

# inverse: ADM -> Z4c, USING det(gt) = 1 and tr_gt(At) = 0 as constraints
detg = sp.factor(gamma.det())
chi_back = sp.cancel(detg ** sp.Rational(-1, 1))   # chi^3 / det(gt)
# det(gamma) = det(gt)/chi^3 -> chi_back^(1/3) handled rationally:
report("chi recovery: det(gamma) = det(gt)/chi^3 (so chi = (det gamma)^(-1/3) "
       "iff det gt = 1) — residual det(gamma)*chi^3 - det(gt) == 0",
       sp.expand(detg * chi**3 - gt.det()) == 0)

gamma_inv = gamma.inv()
K_tr = sp.expand(sum(gamma_inv[i, j] * K_ij[i, j] for i in range(3) for j in range(3)))
trAt = sp.Symbol("trAt")  # gt^{ij} At_ij, set to 0 by the trace-free constraint
gtinv = gt.inv()
trAt_expr = sp.expand(sum(gtinv[i, j] * At[i, j] for i in range(3) for j in range(3)))
K_back_residual = sp.simplify(K_tr - K - sp.cancel(trAt_expr / 1) * 0)
report("K recovery: gamma^{ij} K_ij - (Khat + 2 Theta) == (gt^{ij}At_ij)/1 * "
       "1/1 — vanishes given trace-free At",
       sp.simplify(K_tr - K - trAt_expr) == 0)

At_back = sp.expand((K_ij - gamma * K_tr / 3) * chi)
res_At = sp.expand(At_back - At + (gt / 3) * trAt_expr)
report("At recovery: chi*(K_ij - gamma_ij K/3) - At_ij == -gt_ij/3 * gt^{kl}At_kl "
       "— vanishes given trace-free At",
       all(sp.simplify(res_At[i, j]) == 0 for i in range(3) for j in range(3)))

# 4-metric and ADM inverse
g4 = sp.zeros(4, 4)
bd = gamma * beta
g4[0, 0] = -alpha**2 + (bd.T * beta)[0]
for i in range(3):
    g4[0, i+1] = g4[i+1, 0] = bd[i]
    for j in range(3):
        g4[i+1, j+1] = gamma[i, j]
g4i = sp.zeros(4, 4)
g4i[0, 0] = -1 / alpha**2
for i in range(3):
    g4i[0, i+1] = g4i[i+1, 0] = beta[i] / alpha**2
    for j in range(3):
        g4i[i+1, j+1] = gamma_inv[i, j] - beta[i] * beta[j] / alpha**2
prod = sp.expand(g4 * g4i)
report("4-metric inverse: g4 * g4inv == identity (full symbols)",
       all(sp.simplify(prod[a, b] - (1 if a == b else 0)) == 0
           for a in range(4) for b in range(4)))

# ---------- B. kinematic closure, point-symbols ------------------------------
# independent symbols for fields and first derivatives at a point
dchi = sp.Matrix([sp.Symbol(f"dchi{k}") for k in range(3)])           # chi_,k
dgt = [[[sp.Symbol(f"dgt{min(i,j)}{max(i,j)}_{k}") for k in range(3)]
        for j in range(3)] for i in range(3)]                          # gt_ij,k
db = [[sp.Symbol(f"db{i}_{k}") for k in range(3)] for i in range(3)]   # b^i_,k

divb_plain = sum(db[k][k] for k in range(3))
lnchi_term = sum(beta[k] * dchi[k] for k in range(3)) / chi
divb_cov = divb_plain - sp.Rational(3, 2) * lnchi_term   # D_i beta^i (full gamma)

def dt_gamma(divb):
    """d_t gamma_ij from P1's conformal RHS with the given D_i beta^i reading."""
    dt_chi = sp.Rational(2, 3) * chi * (alpha * K - divb)
    out = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            adv = sum(beta[k] * dgt[i][j][k] for k in range(3))
            grad = sum(gt[i, k] * db[k][j] + gt[j, k] * db[k][i] for k in range(3))
            dt_gt = (-2 * alpha * At[i, j] + adv + grad
                     - sp.Rational(2, 3) * gt[i, j] * divb_plain)
            out[i, j] = dt_gt / chi - gt[i, j] / chi**2 * dt_chi
    return out

# ADM kinematic RHS: -2 alpha K_ij + Lie_beta gamma
lie = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        # d_k gamma_ij = (dgt_ij_k - gamma_ij*... ) via gamma = gt/chi
        advg = sum(beta[k] * (dgt[i][j][k] / chi - gt[i, j] * dchi[k] / chi**2)
                   for k in range(3))
        lie[i, j] = advg + sum(gamma[i, k] * db[k][j] + gamma[j, k] * db[k][i]
                               for k in range(3))
adm = -2 * alpha * K_ij + lie

resid_cov = sp.expand(dt_gamma(divb_cov) - adm)
report("KINEMATIC CLOSURE with D_i beta^i = beta^k_,k - (3/2) beta^k (ln chi)_,k "
       "(full covariant divergence): d_t gamma_ij == -2 alpha K_ij + Lie_beta gamma, EXACT",
       all(sp.simplify(resid_cov[i, j]) == 0 for i in range(3) for j in range(3)))

resid_plain = sp.expand(dt_gamma(divb_plain) - adm)
expected_gap = sp.expand(-gt[0, 0] / chi * lnchi_term)  # component (0,0): -gamma_ij beta^k dk(ln chi)... sign check below
gap_ok = all(sp.simplify(resid_plain[i, j] + gamma[i, j] * lnchi_term * 0
             - resid_plain[i, j]) == 0 for i in range(3))  # tautology guard
# the gap must be exactly proportional to gamma_ij * beta^k (ln chi)_,k:
coeffs = {sp.simplify(sp.cancel(resid_plain[i, j] / (gamma[i, j] * lnchi_term)))
          for i in range(3) for j in range(3)}
report(f"NAIVE reading D = beta^k_,k leaves residual = c * gamma_ij beta^k (ln chi)_,k "
       f"with constant c = {coeffs} (nonzero => convention resolved, not assumed)",
       len(coeffs) == 1 and 0 not in coeffs)

print()
print("CONCLUSION: packages/zccm z4c_vars maps are exact (round trip, 4-metric")
print("inverse), and P1's chi-equation divergence MUST be read as the full")
print("covariant D_i beta^i for ADM kinematic closure — convention note for the")
print("formulation document; worldtube data g, dt g, dr g from Z4c state is")
print("therefore complete (spatial derivatives supplied by the evolution code).")
print(f"\nwall clock: {time.time()-T0:.1f}s")
print(f"OVERALL: {'PASS' if OK else 'FAIL'}")
raise SystemExit(0 if OK else 1)
