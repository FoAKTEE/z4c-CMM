# derive_beta_hierarchy.jl — iter 44: machine-derive the beta-corrected
# linear hierarchy on the worldtube-anchored series. The anchored gauge has
# beta != 0 (r-independent at linear order) and non-decaying tails
# (k = 0, 1 monomials) that the beta = 0 Bondi-gauge hierarchy does not
# support: feeding anchored BCs into the ringed solver fails at O(10^3) in
# psi0 (measured iter 44). The anchored series satisfies the full linear
# Bondi-Sachs system BY CONSTRUCTION (it is a gauge transform of an exact
# solution), so each equation's residual against the beta = 0 hierarchy
# must be EXACTLY a beta-linear source term — fit it with exact rationals.
using ZccmJl
const Z = ZccmJl

J, U, W, B = J_ANCH, U_ANCH, W_ANCH, BETA_ANCH

# linear U-eq: dr U = Q/r^2 exactly (e^{2beta} Q with Q0 = 0 background) —
# Q is DEFINED by inversion, no independent fit possible at this order.
Q = Z.ang_rpow(Z.ang_dr(U), 2)
H = Z.ang_du(J)

# shapes of the beta sources (beta is P-shaped, spin 0):
#   eth beta        = -6 beta_t * sc      (spin 1, Q channel)
#   ethbar eth beta = -6 beta_t * P       (spin 0, W channel)
#   eth eth beta    = +6 beta_t * s^2     (spin 2, H channel)
eB = Z.eth_spin(B, 0)
ebeB = Z.ethbar_spin(eB, 1)
eeB = Z.eth_spin(eB, 1)

function fit_source(res::Z.AngSeries, base::Z.AngSeries, powers, name)
    basis = Dict(p => Z.ang_rpow(base, p) for p in powers)
    chans = Set{NTuple{4,Int}}(keys(res))
    for b in values(basis)
        union!(chans, keys(b))
    end
    chans = sort(collect(chans))
    syms = sort(collect(keys(basis)))
    A = Rational{BigInt}[get(basis[sy], ch, big(0)//1)
                         for ch in chans, sy in syms]
    rhs = Rational{BigInt}[get(res, ch, big(0)//1) for ch in chans]
    x, ok = Z.exact_solve(A, rhs)
    if ok
        terms = ["($(x[i])) * r^$(syms[i])" for i in eachindex(syms) if x[i] != 0]
        println("$name source = [", join(terms, " + "), "] * (shape)")
    else
        println("$name: NOT representable in the given basis; residual channels:")
        for ch in chans
            v = get(res, ch, big(0)//1)
            v != 0 && println("   $ch => $v")
        end
    end
    (x, ok, syms)
end

# (1) Q hypersurface: dr(r^2 Q) + r^2 dr(ethbar J) = S_Q[beta]
lhsQ = Z.ang_dr(Z.ang_rpow(Q, 2))
rhsQ = Z.ascale(Z.ang_rpow(Z.ang_dr(Z.ethbar_spin(J, 2)), 2), -1)
resQ = Z.aadd(lhsQ, Z.ascale(rhsQ, -1))
fit_source(resQ, eB, -2:2, "Q-eq")

# (2) W hypersurface: dr(r^2 W) - [(1/2) ethb ethb J + 2 r ethbar U
#     + (r^2/2) dr(ethbar U)] = S_W[beta]
ebbJ = Z.ethbar_spin(Z.ethbar_spin(J, 2), 1)
eUb = Z.ethbar_spin(U, 1)
sumU = Z.ascale(eUb, 2)
lhsW = Z.ang_dr(Z.ang_rpow(W, 2))
rhsW = Z.aadd(Z.ascale(ebbJ, 1//2),
              Z.aadd(Z.ang_rpow(sumU, 1),
                     Z.ascale(Z.ang_rpow(Z.ang_dr(sumU), 2), 1//4)))
resW = Z.aadd(lhsW, Z.ascale(rhsW, -1))
fit_source(resW, ebeB, -2:2, "W-eq")

# (3) H evolution: dr(r H) - [aH eth U + bH r dr^2 J + cH dr J
#     + (eH eth Q + cJ J)/r] = S_H[beta], with the iter-30 pinned
#     coefficients aH = -1/8*... — read them from hierarchy_closure
hc = hierarchy_closure()
println("plain-hierarchy H coefficients: ", hc.h_coeffs)
co = hc.h_coeffs
lhsH = Z.ang_dr(Z.ang_rpow(H, 1))
rhsH = Z.aadd(Z.ascale(Z.eth_spin(U, 1), co[:aH]),
       Z.aadd(Z.ascale(Z.ang_rpow(Z.ang_dr(Z.ang_dr(J)), 1), co[:bH]),
       Z.aadd(Z.ascale(Z.ang_dr(J), co[:cH]),
       Z.aadd(Z.ascale(Z.ang_rpow(Z.eth_spin(Q, 1), -1), co[:eH]),
              Z.ascale(Z.ang_rpow(J, -1), co[:cJ])))))
resH = Z.aadd(lhsH, Z.ascale(rhsH, -1))
fit_source(resH, eeB, -2:2, "H-eq")

# The residual is NOT pure ð²β: the plain-family fit of the H coefficients
# was degenerate (the plain solution manifold cannot distinguish all basis
# elements). JOINT refit: fit dr(rH) against the full basis ON THE ANCHORED
# series (richer family), then gate the same equation on the plain series.
function full_basis(Jx, Ux, Qx, Bx)
    duBx = Z.ang_du(Bx)
    b = Dict{Tuple{Symbol,Int},Z.AngSeries}()
    add!(nm, ser, powers) = for p in powers
        b[(nm, p)] = Z.ang_rpow(ser, p)
    end
    add!(:ethU, Z.eth_spin(Ux, 1), -1:1)
    add!(:d2J, Z.ang_dr(Z.ang_dr(Jx)), 0:1)
    add!(:dJ, Z.ang_dr(Jx), -1:0)
    add!(:J, Jx, -2:-1)
    add!(:ethQ, isempty(Qx) ? Z.AngSeries() : Z.eth_spin(Qx, 1), -1:0)
    if !isempty(Bx)
        eeBx = Z.eth_spin(Z.eth_spin(Bx, 0), 1)
        add!(:eeB, eeBx, -1:1)
        if !isempty(duBx)
            add!(:eeduB, Z.eth_spin(Z.eth_spin(duBx, 0), 1), 0:1)
        end
    end
    b
end
bA = full_basis(J, U, Q, B)
lhsHA = Z.ang_dr(Z.ang_rpow(Z.ang_du(J), 1))
bt = bondi_transform_teukolsky()
QP = Z.ang_rpow(Z.ang_dr(bt.U), 2)
bP = full_basis(bt.J, bt.U, QP, Z.AngSeries())
lhsHP = Z.ang_dr(Z.ang_rpow(Z.ang_du(bt.J), 1))
# joint channel set: anchored channels tagged :A, plain tagged :P
syms = sort(collect(union(keys(bA), keys(bP))))
chansA = Set{NTuple{4,Int}}(keys(lhsHA))
for s in syms
    haskey(bA, s) && union!(chansA, keys(bA[s]))
end
chansP = Set{NTuple{4,Int}}(keys(lhsHP))
for s in syms
    haskey(bP, s) && union!(chansP, keys(bP[s]))
end
rows = vcat([(:A, ch) for ch in sort(collect(chansA))],
            [(:P, ch) for ch in sort(collect(chansP))])
A = Rational{BigInt}[
    (tag == :A ? get(get(bA, sy, Z.AngSeries()), ch, big(0)//1)
               : get(get(bP, sy, Z.AngSeries()), ch, big(0)//1))
    for (tag, ch) in rows, sy in syms]
rhs = Rational{BigInt}[
    (tag == :A ? get(lhsHA, ch, big(0)//1) : get(lhsHP, ch, big(0)//1))
    for (tag, ch) in rows]
x, ok = Z.exact_solve(A, rhs)
println("H-eq JOINT fit (anchored + plain): consistent = $ok")
ok && for i in eachindex(syms)
    x[i] != 0 && println("   $(syms[i]) => $(x[i])")
end
