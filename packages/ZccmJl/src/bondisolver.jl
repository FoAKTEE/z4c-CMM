# bondisolver.jl — numeric linear l=2 (m=0) characteristic solver on the
# compactified grid y = 1 - 2 r_wt / r in [-1, 1] (r_wt at y = -1, scri at
# y = +1). Scalarized system (shape constants validated by hierarchy.jl):
#   dr(r^2 Qt) = -4 r^2 dr Jt
#   dr Ut      = Qt / r^2
#   dr(r^2 Wt) = 2 Jt + 2 r Ut + (r^2/2) dr Ut
#   dr(r Ht)   = -(1/2) Ut + (1/4) r dr^2 Jt - (3/2) Jt / r,   du Jt = Ht
#   psi0t      = -dr Jt/(2r) - dr^2 Jt/4
# (J = s^2 Jt, Q = sc Qt, U = sc Ut, W = P Wt, H = s^2 Ht, psi0 = s^2 psi0t,
#  P = 3c^2 - 1.) Radial scheme: FD differentiation matrix (4th order) +
# cumulative integration from the worldtube; u-evolution: classical RK4.
# Worldtube BCs supplied by a closure (for the standalone gate: the exact
# Teukolsky closed forms; in live CCM: the worldtube transformation).

export BondiGrid, hierarchy_sweep!, evolve_J, psi0_worldtube,
       teuk_Jt, teuk_Ht, teuk_Qt, teuk_Ut, teuk_Wt

struct BondiGrid{T<:Real}
    rwt::T
    y::Vector{T}      # uniform in [-1, 1]
    r::Vector{T}
    D::Matrix{T}      # d/dy, 4th-order centered + one-sided edges
end

function BondiGrid(rwt::T, ny::Int) where {T<:Real}
    y = collect(range(T(-1), T(1), length=ny))
    h = y[2] - y[1]
    D = zeros(T, ny, ny)
    for i in 3:ny-2
        D[i, i-2:i+2] .= T[1, -8, 0, 8, -1]./(12h)
    end
    # 4th-order one-sided rows
    D[1, 1:5] .= T[-25, 48, -36, 16, -3]./(12h)
    D[2, 1:5] .= T[-3, -10, 18, -6, 1]./(12h)
    D[ny-1, ny-4:ny] .= -T[1, -6, 18, -10, -3]./(12h)
    D[ny, ny-4:ny] .= -T[-3, 16, -36, 48, -25]./(12h)
    # r(y): scri point y = 1 maps to r = Inf — use a large finite cap via
    # y capped slightly below 1 to keep the linear gates finite
    yc = copy(y)
    yc[end] = y[end] - h/1000
    r = @. 2rwt/(1 - yc)
    BondiGrid{T}(rwt, y, r, D)
end

"""d/dr on the grid: dr = ((1-y)^2 / (2 rwt)) d/dy."""
function ddr(g::BondiGrid{T}, f::Vector{T}) where {T}
    fac = @. (1 - g.y)^2/(2g.rwt)
    fac .* (g.D*f)
end

"""cumulative integral from y[1] (worldtube), 4th-order and fully LOCAL
(marching with the integral of the cubic through 4 neighboring samples):
inner values never couple to the scri-endpoint divergence of raw-r
integrands — the global matrix variant did and diverged (ledgered)."""
function cumint(g::BondiGrid{T}, f::Vector{T}) where {T}
    n = length(f)
    h = g.y[2] - g.y[1]
    out = zeros(T, n)
    for i in 1:n-1
        # integrate over [y_i, y_{i+1}] the cubic through samples j0..j0+3
        j0 = clamp(i - 1, 1, n - 3)
        a = i - j0   # local offset of the interval start (0, 1 or 2)
        # exact weights: integral of Lagrange cubics on [a, a+1] (unit h)
        w = ntuple(m -> begin
            # L_m(x) over nodes 0..3; integrate on [a, a+1]
            nodes = (0, 1, 2, 3)
            num(x) = prod(x - nodes[k] for k in 1:4 if k != m)
            den = prod(nodes[m] - nodes[k] for k in 1:4 if k != m)
            # Simpson-exact integration of the cubic num(x): 2-pt Gauss is
            # exact for cubics
            x1 = a + (1 - 1/sqrt(T(3)))/2
            x2 = a + (1 + 1/sqrt(T(3)))/2
            (num(x1) + num(x2))/(2den)
        end, 4)
        out[i+1] = out[i] + h*sum(w[m]*f[j0+m-1] for m in 1:4)
    end
    out
end

"""One hierarchy sweep on a slice: given Jt(y) and worldtube BCs, return
(Qt, Ut, Wt, Ht). bc = (Qt_wt, Ut_wt, Wt_wt, Ht_wt)."""
function hierarchy_sweep!(g::BondiGrid{T}, Jt::Vector{T},
                          bc::NTuple{4,T}) where {T}
    dJ = ddr(g, Jt)
    d2J = ddr(g, dJ)
    r = g.r
    # Q: r^2 Qt = (rwt^2 Qt_wt) + cumint_r(-4 r^2 dr Jt)
    SQ = @. -4r^2*dJ
    # integrate in r: ∫ S dr = ∫ S (dr/dy) dy, dr/dy = 2rwt/(1-y)^2
    drdy = @. 2g.rwt/(1 - g.y)^2
    drdy[end] = drdy[end-1]    # capped scri point
    r2Q = g.rwt^2*bc[1] .+ cumint(g, SQ .* drdy)
    Qt = r2Q ./ r.^2
    # U
    SU = Qt ./ r.^2
    Ut = bc[2] .+ cumint(g, SU .* drdy)
    dU = ddr(g, Ut)
    # W
    SW = @. 2Jt + 2r*Ut + (r^2/2)*dU
    r2W = g.rwt^2*bc[3] .+ cumint(g, SW .* drdy)
    Wt = r2W ./ r.^2
    # H
    SH = @. -Ut/2 + (r/4)*d2J - (3//2)*Jt/r
    rH = g.rwt*bc[4] .+ cumint(g, SH .* drdy)
    Ht = rH ./ r
    (Qt, Ut, Wt, Ht)
end

"""psi0t at the worldtube from the slice."""
function psi0_worldtube(g::BondiGrid{T}, Jt::Vector{T}) where {T}
    dJ = ddr(g, Jt)
    d2J = ddr(g, dJ)
    -dJ[1]/(2g.rwt) - d2J[1]/4
end

# ---- exact Teukolsky closed forms (scalarized; X, rc, tau via teuk_F) ----
teuk_Jt(u::T, r::T, X, rc, tau) where {T} =
    T(3)/4*(teuk_F(4, u, T(X), T(rc), T(tau))/r
            + teuk_F(2, u, T(X), T(rc), T(tau))/r^3)
teuk_Ht(u::T, r::T, X, rc, tau) where {T} =
    T(3)/4*(teuk_F(5, u, T(X), T(rc), T(tau))/r
            + teuk_F(3, u, T(X), T(rc), T(tau))/r^3)
teuk_Ut(u::T, r::T, X, rc, tau) where {T} =
    (-T(3)/2*teuk_F(4, u, T(X), T(rc), T(tau))/r^2
     + 3teuk_F(3, u, T(X), T(rc), T(tau))/r^3
     + T(9)/4*teuk_F(2, u, T(X), T(rc), T(tau))/r^4)
teuk_Qt(u::T, r::T, X, rc, tau) where {T} =
    (3teuk_F(4, u, T(X), T(rc), T(tau))/r
     - 9teuk_F(3, u, T(X), T(rc), T(tau))/r^2
     - 9teuk_F(2, u, T(X), T(rc), T(tau))/r^3)
teuk_Wt(u::T, r::T, X, rc, tau) where {T} =
    (-T(3)/2*teuk_F(4, u, T(X), T(rc), T(tau))/r^2
     - T(3)/2*teuk_F(3, u, T(X), T(rc), T(tau))/r^3
     - T(3)/4*teuk_F(2, u, T(X), T(rc), T(tau))/r^4)

"""Evolve Jt from u0 to u1 with RK4 (du steps), worldtube data from the
closed forms; returns the final slice."""
function evolve_J(g::BondiGrid{T}, u0::T, u1::T, du::T,
                  X, rc, tau) where {T}
    Jt = T[teuk_Jt(u0, r, X, rc, tau) for r in g.r]
    u = u0
    bcs(uu) = (teuk_Qt(uu, g.rwt, X, rc, tau),
               teuk_Ut(uu, g.rwt, X, rc, tau),
               teuk_Wt(uu, g.rwt, X, rc, tau),
               teuk_Ht(uu, g.rwt, X, rc, tau))
    rhs(uu, J) = begin
        Jb = copy(J)
        Jb[1] = teuk_Jt(uu, g.rwt, X, rc, tau)   # worldtube Dirichlet
        (_, _, _, Ht) = hierarchy_sweep!(g, Jb, bcs(uu))
        Ht
    end
    while u < u1 - du/2
        k1 = rhs(u, Jt)
        k2 = rhs(u + du/2, Jt .+ du/2 .* k1)
        k3 = rhs(u + du/2, Jt .+ du/2 .* k2)
        k4 = rhs(u + du, Jt .+ du .* k3)
        Jt .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
        Jt[1] = teuk_Jt(u, g.rwt, X, rc, tau)
    end
    Jt
end
