# bondisolver_reg.jl — O-N14-1: the (1-y)-REGULARIZED linear l=2 solver
# (the paper's eq:*numeric philosophy applied to the verified linear system
# of iter 31). Ringed variables, all O(1) and finite at scri:
#   jr = r Jt,  qr = r Qt,  ur = r^2 Ut,  wr = r^2 Wt,  hr = r Ht = du jr.
# Exact transformed ODEs (derived from the verified r-forms with
# r = 2 rwt/(1-y), dr = ((1-y)^2/(2 rwt)) dy; scri-limit checks in tests):
#   Q: qr + (1-y) qr' = 4 jr - 4 (1-y) jr'
#   U: (1-y) ur' - 2 ur = qr
#   W: (1-y) wr' = 2 jr + 2 ur + qr/2
#   H: hr' = [ -ur - 2 jr + ((1-y)^2/2) jr'' - 2 (1-y) jr' ] / (4 rwt)
# psi0t at the worldtube from jr (regular variable, bounded weights):
#   Jt = jr/r;  psi0t = -dr Jt/(2 r) - dr^2 Jt / 4.
# Pole ODEs integrated by implicit trapezoid (2nd order, A-stable) marching
# from the worldtube; the scri endpoint set by the regular limit.

export RegGrid, rofy, reg_sweep, reg_evolve_J, reg_psi0_worldtube,
       teuk_jr, teuk_qr, teuk_ur, teuk_wr, teuk_hr

struct RegGrid{T<:Real}
    rwt::T
    y::Vector{T}
    D::Matrix{T}     # d/dy, 4th order (reuse pattern from BondiGrid)
end

function RegGrid(rwt::T, ny::Int) where {T<:Real}
    g = BondiGrid(rwt, ny)
    D = copy(g.D)
    h = g.y[2] - g.y[1]
    # GKS-stable mixed closure (mission-2 model1d lesson, proven there):
    # 4th-order interior + LOW-order edge rows; the 4th-order one-sided
    # rows with RK4 carry a refinement-growing edge mode (rate ~ 1/h;
    # ledgered iter 33).
    D[1, :] .= 0;  D[1, 1:3] .= T[-3, 4, -1]./(2h)
    D[2, :] .= 0;  D[2, 1:3] .= T[-1, 0, 1]./(2h)
    D[ny-1, :] .= 0;  D[ny-1, ny-2:ny] .= T[-1, 0, 1]./(2h)
    D[ny, :] .= 0;  D[ny, ny-2:ny] .= T[1, -4, 3]./(2h)
    RegGrid{T}(rwt, g.y, D)
end

rofy(g::RegGrid{T}, y::T) where {T} = 2g.rwt/(1 - y)

# closed forms in ringed variables, EXPANDED in 1/r so the scri node
# (r = Inf, 1/r = 0) is exact — the r*(F/r) form produced Inf*0 = NaN.
function teuk_jr(u::T, r::T, X, rc, tau) where {T}
    ir = 1/r
    T(3)/4*(teuk_F(4, u, T(X), T(rc), T(tau))
            + teuk_F(2, u, T(X), T(rc), T(tau))*ir^2)
end
function teuk_qr(u::T, r::T, X, rc, tau) where {T}
    ir = 1/r
    (3teuk_F(4, u, T(X), T(rc), T(tau))
     - 9teuk_F(3, u, T(X), T(rc), T(tau))*ir
     - 9teuk_F(2, u, T(X), T(rc), T(tau))*ir^2)
end
function teuk_ur(u::T, r::T, X, rc, tau) where {T}
    ir = 1/r
    (-T(3)/2*teuk_F(4, u, T(X), T(rc), T(tau))
     + 3teuk_F(3, u, T(X), T(rc), T(tau))*ir
     + T(9)/4*teuk_F(2, u, T(X), T(rc), T(tau))*ir^2)
end
function teuk_wr(u::T, r::T, X, rc, tau) where {T}
    ir = 1/r
    (-T(3)/2*teuk_F(4, u, T(X), T(rc), T(tau))
     - T(3)/2*teuk_F(3, u, T(X), T(rc), T(tau))*ir
     - T(3)/4*teuk_F(2, u, T(X), T(rc), T(tau))*ir^2)
end
function teuk_hr(u::T, r::T, X, rc, tau) where {T}
    ir = 1/r
    T(3)/4*(teuk_F(5, u, T(X), T(rc), T(tau))
            + teuk_F(3, u, T(X), T(rc), T(tau))*ir^2)
end

"""implicit-trapezoid march of a' = (S(y) - a*alpha(y))/(1-y) type pole ODEs,
written generically as a' = (c1(y)*a + src(y))/(1-y); the last node takes
the regular limit value lim_val."""
function pole_march(y::Vector{T}, a1::T, c1::Vector{T}, src::Vector{T},
                    lim_val::T) where {T}
    n = length(y)
    a = zeros(T, n)
    a[1] = a1
    for i in 1:n-2
        h = y[i+1] - y[i]
        fi = (c1[i]*a[i] + src[i])/(1 - y[i])
        # implicit: a_{i+1} = a_i + h/2 (f_i + f_{i+1}),
        # f_{i+1} = (c1*a + src)_{i+1}/(1-y_{i+1})
        d = 1 - y[i+1]
        a[i+1] = (a[i] + h/2*(fi + src[i+1]/d))/(1 - h/2*c1[i+1]/d)
    end
    a[n] = lim_val
    a
end

"""hierarchy sweep in regularized variables. Worldtube BCs (qr, ur, wr, hr)
at y = -1; scri limits supplied analytically by the caller via jr."""
function reg_sweep(g::RegGrid{T}, jr::Vector{T},
                   bc::NTuple{4,T}) where {T}
    djr = g.D*jr
    d2jr = g.D*djr
    n = length(g.y)
    # Q: qr' = (-qr + 4 jr - 4 (1-y) jr')/(1-y): c1 = -1,
    # src = 4 jr - 4 (1-y) jr'; scri limit qr = 4 jr (src at y=1)
    srcQ = T[4jr[i] - 4*(1 - g.y[i])*djr[i] for i in 1:n]
    qr = pole_march(g.y, bc[1], fill(T(-1), n), srcQ, 4jr[end])
    # U: ur' = (2 ur + qr)/(1-y). The homogeneous mode ~ (1-y)^-2 GROWS
    # toward scri: outward marching is ill-posed (implicit-trapezoid
    # denominator crosses zero — ledgered). March INWARD from scri (stable
    # direction) starting at the regularity value ur(scri) = -qr(scri)/2;
    # the worldtube BC becomes a consistency CHECK, returned to the caller.
    ur = zeros(T, n)
    ur[n] = -qr[n]/2
    for i in n-1:-1:1
        h = g.y[i] - g.y[i+1]          # negative
        d1 = 1 - g.y[i+1]
        d0 = 1 - g.y[i]
        # at scri use the L'Hopital limit: 3 ur' = -qr'  (no ur values
        # needed — the earlier form read the not-yet-filled array; ledgered)
        f1 = d1 == 0 ? -(g.D*qr)[n]/T(3) :
                        (2ur[i+1] + qr[i+1])/d1
        # implicit at i: ur_i = ur_{i+1} + h/2 (f1 + (2 ur_i + qr_i)/d0)
        ur[i] = (ur[i+1] + h/2*(f1 + qr[i]/d0))/(1 - h*1/d0)
    end
    # W: wr' = (2jr + 2ur + qr/2)/(1-y), no homogeneous growth (c1 = 0);
    # forward march; the last node by 4th-order extrapolation (0/0 there)
    srcW = T[2jr[i] + 2ur[i] + qr[i]/2 for i in 1:n]
    wr = pole_march(g.y, bc[3], zeros(T, n), srcW, T(0))
    wr[end] = 4wr[end-1] - 6wr[end-2] + 4wr[end-3] - wr[end-4]
    # H by ANALYTIC integration by parts (removes the D^2-composition edge
    # modes that drove the refinement-growing layer — ledgered):
    #   integral of (1/4) r' d2r Jt dr' = (1/4)[r dr Jt - Jt] exactly, so
    #   hr(y) = hr_wt + (1/4)[G(y) - G(-1)]
    #           - (1/(4 rwt)) int ur dy' - (3/(4 rwt)) int jr dy',
    #   G = ((1-y)^2/(2 rwt)) dy jr - (1-y) jr / rwt.
    # Only FIRST derivatives of jr appear.
    G = T[(1 - g.y[i])^2/(2g.rwt)*djr[i] - (1 - g.y[i])*jr[i]/g.rwt
          for i in 1:n]
    cumU = zeros(T, n); cumJ = zeros(T, n)
    for i in 1:n-1
        h = g.y[i+1] - g.y[i]
        cumU[i+1] = cumU[i] + h/2*(ur[i] + ur[i+1])
        cumJ[i+1] = cumJ[i] + h/2*(jr[i] + jr[i+1])
    end
    hr = T[bc[4] + (G[i] - G[1])/4 - cumU[i]/(4g.rwt) - 3cumJ[i]/(4g.rwt)
           for i in 1:n]
    (qr, ur, wr, hr)
end

"""psi0t at the worldtube from jr (one-sided 4th-order y-derivatives)."""
function reg_psi0_worldtube(g::RegGrid{T}, jr::Vector{T}) where {T}
    djr = (g.D*jr)[1]
    d2jr = (g.D*(g.D*jr))[1]
    rwt = g.rwt
    y1 = g.y[1]                      # -1
    w = (1 - y1)^2/(2rwt)            # dr -> w * dy at the worldtube
    # Jt = jr/r: dJt/dr = (w djr)/r - jr/r^2; r = rwt at y1
    dJ = w*djr/rwt - jr[1]/rwt^2
    # d2Jt/dr2 = (w d/dy)(dJt/dr): compute via chain on the expression
    # d2(jr/r) = w*d/dy[w djr/r - jr/r^2]
    #          = w*(dw/dy*djr/r + w*d2jr/r - w*djr/r^2*drdy_inv ... ) —
    # safer: use dr(1/r) = -1/r^2 and dr(jr) = w djr:
    # d2Jt = dr(dJ) = dr(w djr)/r - (w djr)/r^2 - (w djr)/r^2 + 2 jr/r^3
    # dr(w djr) = w * d/dy(w djr) = w*(dwdy*djr + w*d2jr)
    dwdy = -2*(1 - y1)/(2rwt)        # d/dy of (1-y)^2/(2rwt)
    drwdjr = w*(dwdy*djr + w*d2jr)
    d2J = drwdjr/rwt - 2w*djr/rwt^2 + 2jr[1]/rwt^3
    -dJ/(2rwt) - d2J/4
end

"""RK4 evolution of jr with du jr = hr; worldtube data from closed forms."""
function reg_evolve_J(g::RegGrid{T}, u0::T, u1::T, du::T,
                      X, rc, tau) where {T}
    rwt = g.rwt
    jr = T[teuk_jr(u0, rofy(g, y), X, rc, tau) for y in g.y]
    u = u0
    bcs(uu) = (teuk_qr(uu, rwt, X, rc, tau), teuk_ur(uu, rwt, X, rc, tau),
               teuk_wr(uu, rwt, X, rc, tau), teuk_hr(uu, rwt, X, rc, tau))
    # BC imposition (iter 34 ledger): SAT penalty FAILED both regimes —
    # O(1) strength is stiff for explicit RK4 (du*tau/h >> 2.8: NaN) and the
    # v-scaled SBP strength under-tracks the time-dependent inflow BC
    # (boundary lag; psi0 rel 49-63, WORSE than pinning's 1.6-5.0).
    # Default = Dirichlet pinning (best known); O-N14-1 continues via the
    # hierarchy-based psi0 (no edge d2 of the evolved field) next.
    rhs(uu, J) = begin
        Jb = copy(J)
        Jb[1] = teuk_jr(uu, rwt, X, rc, tau)
        (_, _, _, hr) = reg_sweep(g, Jb, bcs(uu))
        hr
    end
    while u < u1 - du/2
        k1 = rhs(u, jr)
        k2 = rhs(u + du/2, jr .+ du/2 .* k1)
        k3 = rhs(u + du/2, jr .+ du/2 .* k2)
        k4 = rhs(u + du, jr .+ du .* k3)
        jr .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
        jr[1] = teuk_jr(u, rwt, X, rc, tau)
    end
    jr
end
