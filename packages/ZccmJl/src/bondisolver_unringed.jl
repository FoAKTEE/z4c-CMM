# bondisolver_unringed.jl — iter 44: UNRINGED spectral scheme for the
# beta-corrected linear hierarchy (anchored-gauge capable).
#
# Why: the worldtube-anchored gauge (the live-data gauge of the six-row
# map) has beta != 0 and non-decaying tails (k = 0, 1 monomials in J, U,
# W). The ringed variables jr = r J, ur = r^2 U, ... DIVERGE at scri on
# such data — the ringed solver cannot represent the anchored solution at
# all (measured failure: psi0 off by O(10^3), iter 44). In the plain
# t-scalars (Jt, Qt, Ut, Wt, Ht) every equation is regular on y in
# [-1, 1] and the anchored solution is finite everywhere.
#
# Machine-derived system (test/derive_beta_hierarchy.jl; exact-rational
# fits on the anchored + plain families jointly; beta r-independent):
#   Q:  dr(r^2 Q) = -r^2 dr(ethbar J) - 4 r eth(beta)
#   U:  dr U = Q / r^2                  (no beta source at linear order)
#   W:  dr(r^2 W) = (1/2) ethb ethb J + 2 r ethbar U
#                   + (r^2/2) dr(ethbar U) - (4/3) ethbar eth beta
#   H:  dr(r H) = (1/2) r dr^2 J + dr J + eth eth beta / r
#                 - (1/2) eth Q / r - eth U
# t-scalar forms (l = 2, m = 0; J = Jt s^2, Q = Qt sc, U = Ut sc,
# W = Wt P, beta = bt P; eth beta = -6 bt sc, ethb eth beta = -6 bt P,
# eth eth beta = 6 bt s^2, eth Q = -Qt s^2, eth U = -Ut s^2):
#   2 Qt + (1-y) Qt' = -4 (1-y) Jt' + 24 bt
#   Ut' = Qt / (2 rwt)
#   2 Wt + (1-y) Wt' = (1-y) Jt / rwt + 2 Ut + (1-y) Ut'/2 + 4 bt (1-y)/rwt
#   Ht + (1-y) Ht' = ((1-y)/(4 rwt)) [(1-y)^2 Jt'' - 2 (1-y) Jt']
#                    + ((1-y)^2/(2 rwt)) Jt' + 3 bt (1-y)/rwt
#                    + Qt (1-y)/(4 rwt) + Ut
# scri rows are AUTOMATIC (regular limits; no L'Hopital, no IBP).
# psi0t = -dr Jt/(2 r) - dr^2 Jt/4 (J-only; gauge-invariant — verified
# symbolically: psi0_from_J(J_ANCH) == psi0_from_J(J_bondi)).

using LinearAlgebra

export UnringedOps, ut_sweep, ut_evolve_sat, ut_psi0_at, ut_pulse_id,
       ut_evolve_sat_probe

struct UnringedOps{T<:Real}
    FQ; FU; FW; FH
    g::ChebRegGrid{T}
end

function UnringedOps(g::ChebRegGrid{T}) where {T}
    n = length(g.y)
    one_y = 1 .- g.y
    AQ = 2*Matrix{T}(I, n, n) + Diagonal(one_y)*g.D
    AQ[1, :] .= 0; AQ[1, 1] = 1
    AU = copy(g.D); AU[1, :] .= 0; AU[1, 1] = 1
    AW = 2*Matrix{T}(I, n, n) + Diagonal(one_y)*g.D
    AW[1, :] .= 0; AW[1, 1] = 1
    AH = Matrix{T}(I, n, n) + Diagonal(one_y)*g.D
    AH[1, :] .= 0; AH[1, 1] = 1
    UnringedOps{T}(lu(AQ), lu(AU), lu(AW), lu(AH), g)
end

"""beta-corrected hierarchy sweep on unringed t-scalars. bc is the tube
tuple (Qt, Ut, Wt, Ht); bt the (r-independent) beta t-scalar."""
function ut_sweep(ops::UnringedOps{T}, Jt::Vector{T},
                  bc::NTuple{4,T}, bt::T) where {T}
    g = ops.g
    n = length(g.y)
    one_y = 1 .- g.y
    djt = g.D*Jt
    d2jt = g.D*djt
    bQ = -4 .* one_y .* djt .+ 24bt
    bQ[1] = bc[1]
    Qt = ops.FQ\bQ
    bU = Qt ./ (2g.rwt); bU[1] = bc[2]
    Ut = ops.FU\bU
    dut = g.D*Ut
    bW = one_y .* Jt ./ g.rwt .+ 2 .* Ut .+ one_y .* dut ./ 2 .+
         4bt .* one_y ./ g.rwt
    bW[1] = bc[3]
    Wt = ops.FW\bW
    bH = (one_y ./ (4g.rwt)) .* (one_y.^2 .* d2jt .- 2 .* one_y .* djt) .+
         (one_y.^2 ./ (2g.rwt)) .* djt .+ 3bt .* one_y ./ g.rwt .+
         Qt .* one_y ./ (4g.rwt) .+ Ut
    bH[1] = bc[4]
    Ht = ops.FH\bH
    (Qt, Ut, Wt, Ht)
end

"""SAT RK4 evolution of unringed Jt with a worldtube-BC closure:
bcfun(u) -> (Jt, Qt, Ut, Wt, Ht, beta_t) tube scalars (e.g. the
worldtube_map outputs scaled by 1/rwt-powers)."""
function ut_evolve_sat(g::ChebRegGrid{T}, u0::T, u1::T, du::T,
                       Jt0::Vector{T}, bcfun; sigma::T = T(1)) where {T}
    ops = UnringedOps(g)
    N = length(g.y) - 1
    taupen = sigma*(1/(2g.rwt))/(T(1)/N^2)
    Jt = copy(Jt0)
    u = u0
    rhs(uu, Jx) = begin
        b = bcfun(uu)
        (_, _, _, Ht) = ut_sweep(ops, Jx, (b[2], b[3], b[4], b[5]), b[6])
        Ht[1] += taupen*(b[1] - Jx[1])
        Ht
    end
    while u < u1 - du/2
        k1 = rhs(u, Jt); k2 = rhs(u + du/2, Jt .+ du/2 .* k1)
        k3 = rhs(u + du/2, Jt .+ du/2 .* k2); k4 = rhs(u + du, Jt .+ du .* k3)
        Jt .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
    end
    Jt
end

"""2308.10361 Sec V.C characteristic-pulse initial data, t-scalar nodal
values: J = +2Y20 * Jcal(y) with +2Y20 = (1/4) sqrt(15/2pi) s^2, so
Jt = (1/4) sqrt(15/2pi) * Jcal;  Jcal = 4 Z (ymax-y)(y-ymin)/(ymax-ymin)^2
* exp(-((y-yc)/tauy)^2) on [ymin, ymax], 0 outside (the support-edge kink
is at the e^-28 level of peak — spectrally smooth). Paper parameters:
yc = 0, ymin/max = -+0.8, tauy = 0.15, Z = 1e-3."""
function ut_pulse_id(g::ChebRegGrid{T}, Z::T; ymin::T = T(-0.8),
                     ymax::T = T(0.8), yc::T = T(0),
                     tauy::T = T(0.15)) where {T}
    pref = sqrt(T(15)/(2*T(pi)))/4
    map(g.y) do y
        (y <= ymin || y >= ymax) ? zero(T) :
            pref*4Z*(ymax - y)*(y - ymin)/(ymax - ymin)^2 *
            exp(-((y - yc)/tauy)^2)
    end
end

"""ut_evolve_sat + a psi0 probe time series: returns (Jt_final, us, vals)
where vals[i] = psi0t(rprobe) on the cone us[i] (every `stride` substeps)."""
function ut_evolve_sat_probe(g::ChebRegGrid{T}, u0::T, u1::T, du::T,
                             Jt0::Vector{T}, bcfun, rprobe::T;
                             sigma::T = T(1), stride::Int = 32) where {T}
    ops = UnringedOps(g)
    N = length(g.y) - 1
    taupen = sigma*(1/(2g.rwt))/(T(1)/N^2)
    Jt = copy(Jt0)
    u = u0
    us = T[u]; vals = T[ut_psi0_at(g, Jt, rprobe)]
    rhs(uu, Jx) = begin
        b = bcfun(uu)
        (_, _, _, Ht) = ut_sweep(ops, Jx, (b[2], b[3], b[4], b[5]), b[6])
        Ht[1] += taupen*(b[1] - Jx[1])
        Ht
    end
    step = 0
    while u < u1 - du/2
        k1 = rhs(u, Jt); k2 = rhs(u + du/2, Jt .+ du/2 .* k1)
        k3 = rhs(u + du/2, Jt .+ du/2 .* k2); k4 = rhs(u + du, Jt .+ du .* k3)
        Jt .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
        step += 1
        if step % stride == 0
            push!(us, u); push!(vals, ut_psi0_at(g, Jt, rprobe))
        end
    end
    (Jt, us, vals)
end

"""psi0t at radius rstar >= rwt from the unringed Jt state (barycentric
interpolation + spectral y-derivatives; J-only, gauge-invariant)."""
function ut_psi0_at(g::ChebRegGrid{T}, Jt::Vector{T}, rstar::T) where {T}
    n = length(g.y)
    ystar = 1 - 2g.rwt/rstar
    wb = T[(-1)^k for k in 0:n-1]; wb[1] /= 2; wb[end] /= 2
    hit = findfirst(y -> y == ystar, g.y)
    row = zeros(T, n)
    if hit === nothing
        d = wb ./ (ystar .- g.y)
        row .= d ./ sum(d)
    else
        row[hit] = 1
    end
    djt = g.D*Jt; d2jt = g.D*djt
    j1 = sum(row .* djt); j2 = sum(row .* d2jt)
    oy = 1 - ystar
    drJ  = oy^2/(2g.rwt)*j1
    d2rJ = oy^2/(4g.rwt^2)*(oy^2*j2 - 2oy*j1)
    -drJ/(2rstar) - d2rJ/4
end
