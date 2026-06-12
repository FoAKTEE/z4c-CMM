# psi0chain.jl — port of scripts/verify_n14_linear_psi0_chain.py (stage A1).
# NOTE (honest-admission): the A1 J_full gauge shift is UNDER RE-EXAMINATION
# (iter-27 error-DB row: the d_theta T term was dropped; the machine
# re-derivation verify_n14_bondi_transform.py adjudicates). This port
# reproduces the A1 ledgered numbers as the migration baseline; if A1 is
# amended, the J_FULL constant here is amended in the same commit as the
# ledger row.

export J_FULL_A1, psi0_lin_series, psi0_chain_gates

# J_full (A1): coefficients x sin^2(theta), basis F^(n)(u)/r^k
const J_FULL_A1 = LinSeries(
    (4, 1) => big(3)//4, (2, 3) => big(-3)//4,
    (1, 4) => big(-9)//4, (0, 5) => big(-27)//20)

"""psi0_lin = -J'/(2r) - J''/4 on a retarded series (exact)."""
function psi0_lin_series(J::LinSeries)
    dJ = dr_series(J)
    d2J = dr_series(dJ)
    series_add(series_scale(mul_rpow(dJ, -1), -1//2),
               series_scale(d2J, -1//4))
end

"""The A1 gates as exact predicates; returns (pass::Bool, report::String)."""
function psi0_chain_gates()
    J = J_FULL_A1
    c1 = !haskey(J, (3, 2))                       # gauge cancellation
    p0 = psi0_lin_series(J)
    c2 = get(p0, (2, 5), big(0)//1) == 9//8             # leading = datum magnitude
    c4 = all(k >= 5 for ((n, k), c) in p0)         # peeling
    rep = "C1=$(c1) C2=$(c2) (lead $(get(p0,(2,5),0//1))) C4=$(c4); " *
          "psi0 = $(sort(collect(p0)))"
    (c1 && c2 && c4, rep)
end
