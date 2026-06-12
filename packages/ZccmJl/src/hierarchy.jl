# hierarchy.jl — linear l=2 (m=0) hypersurface-hierarchy closure on the
# machine-derived Bondi fields (iter 30: J, beta = 0, U, W). Axisymmetric
# spin-weighted eth in the 2007.01339 q^A qbar_A = 2 convention:
#   eth    f_s = d_theta f - s cot(theta) f      (m = 0)
#   ethbar f_s = d_theta f + s cot(theta) f
# cot on q=0 terms is polynomial ONLY when the term carries an s^2 factor;
# spin-2 fields are stored s^2-reduced as (p=0, p=2) pairs with opposite
# sign — ang_cot_s2 reconstructs the factor and returns c*s * (the cofactor).

export eth_spin, ethbar_spin, hierarchy_closure

"""cot * f for an s^2-divisible q=0 series (stored as (v at p) + (-v at p+2)
pairs, i.e. f = s^2 * sum v c^p). Returns q=1 series c^{p+1} s * v."""
function ang_cot_s2(a::AngSeries)
    # reconstruct cofactor: f = s^2 * g, g[(n,k,p)] = v from pairs
    g = Dict{NTuple{3,Int},Rational{BigInt}}()
    rem = copy(a)
    while !isempty(rem)
        # process in ascending p so the s^2-cofactor peels correctly
        key = sort(collect(keys(rem)), by = k -> (k[1], k[2], k[3]))[1]
        (n, k, p, q) = key
        v = rem[key]
        q == 0 || error("cot_s2 on q=1 term")
        g[(n, k, p)] = get(g, (n, k, p), big(0)//1) + v
        delete!(rem, key)
        partner = (n, k, p+2, 0)
        pv = get(rem, partner, big(0)//1) + v
        pv == 0 ? delete!(rem, partner) : (rem[partner] = pv)
    end
    # validity: reconstruct and compare
    chk = AngSeries()
    for ((n, k, p), v) in g
        aadd!(chk, (n, k, p, 0), v); aadd!(chk, (n, k, p+2, 0), -v)
    end
    chk == a || error("series is not s^2-divisible")
    out = AngSeries()
    for ((n, k, p), v) in g
        aadd!(out, (n, k, p+1, 1), v)   # cot * s^2 c^p = c^{p+1} s
    end
    out
end

"""spin-aware eth/ethbar for the shapes appearing at l=2 m=0."""
function eth_spin(a::AngSeries, s::Int)
    isempty(a) && return AngSeries()
    cot_part = s == 0 ? AngSeries() :
        (all(k[4] == 1 for k in keys(a)) ? ang_cot_mul(a) : ang_cot_s2(a))
    aadd(ang_dtheta(a), ascale(cot_part, -s))
end
function ethbar_spin(a::AngSeries, s::Int)
    isempty(a) && return AngSeries()
    cot_part = s == 0 ? AngSeries() :
        (all(k[4] == 1 for k in keys(a)) ? ang_cot_mul(a) : ang_cot_s2(a))
    aadd(ang_dtheta(a), ascale(cot_part, s))
end

"""Gate the linear hierarchy on the derived fields. Returns a NamedTuple
(q_residual, w_residual, h_report)."""
function hierarchy_closure()
    bt = bondi_transform_teukolsky()
    J, U, W = bt.J, bt.U, bt.W
    @assert isempty(bt.beta)

    # Q from the U equation: dr U = Q / r^2  =>  Q = r^2 dr U  (spin 1)
    Q = ang_rpow(ang_dr(U), 2)

    # (Q1) Q hypersurface equation, linear, beta = 0:
    #   dr(r^2 Q) = -r^2 dr(ethbar J)
    lhsQ = ang_dr(ang_rpow(Q, 2))
    rhsQ = ascale(ang_rpow(ang_dr(ethbar_spin(J, 2)), 2), -1)
    q_res = aadd(lhsQ, ascale(rhsQ, -1))

    # (Q2) W hypersurface equation, linear, beta = 0 (background subtracted):
    #   dr(r^2 W) = (1/4)(ethbar ethbar J + eth eth Jbar)
    #               + r (eth Ubar + ethbar U) + (r^2/4) dr(eth Ubar + ethbar U)
    # m = 0 reality: Jbar has spin -2 with the same shape; eth_{-2} = ethbar_2
    # on the shape level, so eth eth Jbar = ethbar ethbar J.
    ebbJ = ethbar_spin(ethbar_spin(J, 2), 1)
    # eth Ubar (spin -1) = ethbar-type on the shape: d_th + cot
    eUbar = ethbar_spin(U, 1)
    sumU = ascale(eUbar, 2)             # eth Ubar + ethbar U = 2 ethbar U
    lhsW = ang_dr(ang_rpow(W, 2))
    rhsW = aadd(ascale(ebbJ, 1//2),     # (1/4)(2 ebbJ)
                aadd(ang_rpow(sumU, 1),
                     ascale(ang_rpow(ang_dr(sumU), 2), 1//4)))
    w_res = aadd(lhsW, ascale(rhsW, -1))

    # (Q3) H := du J; pin the linear H equation coefficients:
    #   dr(r H) = aH eth U + bH r dr^2 J + cH dr J + (eH eth Q + cJ J)/r
    H = ang_du(J)
    lhsH = ang_dr(ang_rpow(H, 1))
    basis = Dict(
        :aH => eth_spin(U, 1),
        :bH => ang_rpow(ang_dr(ang_dr(J)), 1),
        :cH => ang_dr(J),
        :eH => ang_rpow(eth_spin(Q, 1), -1),
        :cJ => ang_rpow(J, -1))
    # collect channels: lhsH - sum coeff * basis = 0 per (n,k,p,q)
    chans = Set{NTuple{4,Int}}(keys(lhsH))
    for b in values(basis)
        union!(chans, keys(b))
    end
    # build the linear system rows: [basis columns] x = lhsH
    syms = collect(keys(basis))
    A = Rational{BigInt}[get(basis[sy], ch, big(0)//1)
                         for ch in sort(collect(chans)), sy in syms]
    rhs = Rational{BigInt}[get(lhsH, ch, big(0)//1)
                           for ch in sort(collect(chans))]
    # exact least-norm attempt: solve via rationals (A may be tall)
    x, consistent = exact_solve(A, rhs)
    rep = consistent ?
        "H-equation closure CONSISTENT; pinned " *
        join(["$(syms[i]) = $(x[i])" for i in eachindex(syms)], ", ") :
        "H-equation closure INCONSISTENT (residual channels nonzero)"
    (q_res=q_res, w_res=w_res, h_consistent=consistent, h_report=rep,
     h_coeffs=Dict(zip(syms, x)))
end

"""Exact Gaussian elimination for tall systems; returns (x, consistent)."""
function exact_solve(A::Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    m, n = size(A)
    M = hcat(copy(A), copy(b))
    row = 1
    pivots = Int[]
    for col in 1:n
        p = findfirst(i -> M[i, col] != 0, row:m)
        p === nothing && continue
        p += row - 1
        M[[row, p], :] = M[[p, row], :]
        M[row, :] ./= M[row, col]
        for i in 1:m
            i != row && M[i, col] != 0 && (M[i, :] .-= M[i, col]*M[row, :])
        end
        push!(pivots, col)
        row += 1
        row > m && break
    end
    x = zeros(Rational{BigInt}, n)
    for (i, col) in enumerate(pivots)
        x[col] = M[i, end]
    end
    consistent = all(iszero, M[row:end, end])
    (x, consistent)
end
