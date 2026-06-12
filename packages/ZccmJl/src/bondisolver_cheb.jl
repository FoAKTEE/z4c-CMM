# bondisolver_cheb.jl — O-N14-1 (iter 36): Chebyshev spectral radial scheme.
# The iter-35 conditioning analysis: psi0_wt/J ~ 1e-8 (peeling), so the
# field needs ~1e-11 relative accuracy — spectral or bust (FD-4 ceiling
# 3.7e-7). Chebyshev-Gauss-Lobatto collocation in y; pole equations as
# GLOBAL dense solves (the scri row enforces regularity automatically since
# (1-y) = 0 there); integrals via the exact Chebyshev antiderivative.

using LinearAlgebra

export ChebRegGrid, cheb_sweep, cheb_evolve_J, cheb_psi0_hierarchy,
       ChebSweepOps, cheb_sweep_fact, cheb_evolve_J_sat,
       cheb_psi0_at, cheb_evolve_J_sat_bc

struct ChebRegGrid{T<:Real}
    rwt::T
    y::Vector{T}
    D::Matrix{T}
    Aint::Matrix{T}    # cumulative-integral operator from y = -1
end

"""Trefethen Chebyshev differentiation matrix on ascending CGL nodes."""
function cheb_D(n::Int, ::Type{T}) where {T}
    N = n - 1
    x = T[-cospi(T(k)/N) for k in 0:N]      # ascending: -1 ... +1
    c = ones(T, n); c[1] = 2; c[end] = 2
    c .*= T[(-1)^k for k in 0:N]
    D = zeros(T, n, n)
    for i in 1:n, j in 1:n
        if i != j
            D[i, j] = (c[i]/c[j])/(x[i] - x[j])
        end
    end
    for i in 1:n
        D[i, i] = -sum(D[i, j] for j in 1:n if j != i)
    end
    (x, D)
end

"""cumulative integral operator from y = -1: spectrally exact via the
differentiation matrix (solve D F = f with F(-1) = 0; valid here because
the regularized integrands are SMOOTH and finite up to scri — the iter-32
global-solve failure was a divergent raw-r integrand, not this)."""
function cheb_Aint(y::Vector{T}, D::Matrix{T}) where {T}
    n = length(y)
    A = copy(D); A[1, :] .= 0; A[1, 1] = 1
    P = Matrix{T}(I, n, n); P[1, 1] = 0
    A\P
end

function ChebRegGrid(rwt::T, n::Int) where {T<:Real}
    y, D = cheb_D(n, T)
    Aint = cheb_Aint(y, D)
    ChebRegGrid{T}(rwt, y, D, Aint)
end

"""spectral hierarchy sweep (regularized scalars; global solves)."""
function cheb_sweep(g::ChebRegGrid{T}, jr::Vector{T},
                    bc::NTuple{4,T}) where {T}
    n = length(g.y)
    one_y = 1 .- g.y
    djr = g.D*jr
    # Q: qr + (1-y) qr' = 4 jr - 4 (1-y) jr'   [global solve; BC row 1]
    A = Diagonal(ones(T, n)) + Diagonal(one_y)*g.D
    b = 4jr .- 4one_y.*djr
    A = Matrix(A); A[1, :] .= 0; A[1, 1] = 1; b[1] = bc[1]
    qr = A\b
    # U: -2 ur + (1-y) ur' = qr  [global; scri row automatic; BC row 1 =
    # consistency — impose it and let the system pick the regular solution]
    A2 = -2*Matrix{T}(I, n, n) + Diagonal(one_y)*g.D
    b2 = copy(qr); A2[1, :] .= 0; A2[1, 1] = 1; b2[1] = bc[2]
    ur = A2\b2
    # W: (1-y) wr' = 2 jr + 2 ur + qr/2  [global; row n by L'Hopital:
    # -wr'(1) = d/dy(src)(1)]
    A3 = Diagonal(one_y)*g.D
    b3 = 2jr .+ 2ur .+ qr./2
    A3 = Matrix(A3); A3[1, :] .= 0; A3[1, 1] = 1
    srcp = g.D*b3
    A3[n, :] .= -g.D[n, :]; b3w = copy(b3); b3w[1] = bc[3]; b3w[n] = srcp[n]
    wr = A3\b3w
    # H via IBP + spectral antiderivatives:
    G = one_y.^2 ./ (2g.rwt).*djr .- one_y.*jr./g.rwt
    hr = bc[4] .+ (G .- G[1])./4 .- (g.Aint*ur)./(4g.rwt) .-
         3 .*(g.Aint*jr)./(4g.rwt)
    (qr, ur, wr, hr)
end

"""psi0t at the worldtube via the hierarchy (Q-eq jr' + H-identity d2J)."""
function cheb_psi0_hierarchy(g::ChebRegGrid{T}, jr::Vector{T},
                             bc::NTuple{4,T}) where {T}
    qr, ur, wr, hr = cheb_sweep(g, jr, bc)
    rwt = g.rwt
    y1 = g.y[1]
    w = (1 - y1)^2/(2rwt)
    djr1 = (4jr[1] - qr[1] - (1 - y1)*(g.D*qr)[1])/(4*(1 - y1))
    dJt = w*djr1/rwt - jr[1]/rwt^2
    drH = w*(g.D*hr)[1]
    d2Jt = (drH + ur[1]/rwt^2/2 + T(3)/2*(jr[1]/rwt)/rwt)*4/rwt
    -dJt/(2rwt) - d2Jt/4
end

"""RK4 evolution with Dirichlet pinning (worldtube data from closed forms)."""
function cheb_evolve_J(g::ChebRegGrid{T}, u0::T, u1::T, du::T,
                       X, rc, tau) where {T}
    rwt = g.rwt
    jr = T[teuk_jr(u0, rofy_c(g, y), X, rc, tau) for y in g.y]
    u = u0
    bcs(uu) = (teuk_qr(uu, rwt, X, rc, tau), teuk_ur(uu, rwt, X, rc, tau),
               teuk_wr(uu, rwt, X, rc, tau), teuk_hr(uu, rwt, X, rc, tau))
    rhs(uu, J) = begin
        Jb = copy(J)
        Jb[1] = teuk_jr(uu, rwt, X, rc, tau)
        (_, _, _, hr) = cheb_sweep(g, Jb, bcs(uu))
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

rofy_c(g::ChebRegGrid{T}, y::T) where {T} = 2g.rwt/(1 - y)
export rofy_c

"""u-independent factorized sweep operators (LU once, reuse every stage —
BigFloat-feasible; the per-stage dense solves were the runtime wall)."""
struct ChebSweepOps{T<:Real}
    FQ; FU; FW
    g::ChebRegGrid{T}
end

function ChebSweepOps(g::ChebRegGrid{T}) where {T}
    n = length(g.y)
    one_y = 1 .- g.y
    A = Matrix(Diagonal(ones(T, n)) + Diagonal(one_y)*g.D)
    A[1, :] .= 0; A[1, 1] = 1
    A2 = -2*Matrix{T}(I, n, n) + Diagonal(one_y)*g.D
    A2[1, :] .= 0; A2[1, 1] = 1
    A3 = Matrix(Diagonal(one_y)*g.D)
    A3[1, :] .= 0; A3[1, 1] = 1
    A3[n, :] .= -g.D[n, :]
    ChebSweepOps{T}(lu(A), lu(A2), lu(A3), g)
end

function cheb_sweep_fact(ops::ChebSweepOps{T}, jr::Vector{T},
                         bc::NTuple{4,T}) where {T}
    g = ops.g
    n = length(g.y)
    one_y = 1 .- g.y
    djr = g.D*jr
    b = 4jr .- 4one_y.*djr; b[1] = bc[1]
    qr = ops.FQ\b
    b2 = copy(qr); b2[1] = bc[2]
    ur = ops.FU\b2
    b3 = 2jr .+ 2ur .+ qr./2
    srcp = g.D*b3
    b3w = copy(b3); b3w[1] = bc[3]; b3w[n] = srcp[n]
    wr = ops.FW\b3w
    G = one_y.^2 ./ (2g.rwt).*djr .- one_y.*jr./g.rwt
    hr = bc[4] .+ (G .- G[1])./4 .- (g.Aint*ur)./(4g.rwt) .-
         3 .*(g.Aint*jr)./(4g.rwt)
    (qr, ur, wr, hr)
end

"""psi0t at an arbitrary radius rstar >= rwt (interior point of the
characteristic domain) by barycentric interpolation of jr and its spectral
y-derivatives. Needs ONLY the jr state (no sweep, no BCs): the linear
hierarchy gives psi0t = -dr Jt/(2r) - dr^2 Jt/4 with Jt = jr/r. The native
CCM coupling needs this because the Cauchy-boundary datum at radius r_B
lives on the cone u = t - r_B evaluated at r_B > rwt (causality lag)."""
function cheb_psi0_at(g::ChebRegGrid{T}, jr::Vector{T}, rstar::T) where {T}
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
    djr = g.D*jr; d2jr = g.D*djr
    j0 = sum(row .* jr); j1 = sum(row .* djr); j2 = sum(row .* d2jr)
    r = rstar
    dydr = 2g.rwt/r^2; d2ydr2 = -4g.rwt/r^3
    drJ  = j1*dydr/r - j0/r^2
    d2rJ = (j2*dydr^2 + j1*d2ydr2)/r - 2j1*dydr/r^2 + 2j0/r^3
    -drJ/(2r) - d2rJ/4
end

"""SAT evolution with a general worldtube-BC closure: bcfun(u) returns the
ringed tube values (jr, qr, ur, wr, hr) — e.g. the worldtube_map output on
live Cauchy scalars (the production CCM path). jr0 is the initial slice
(start quiescent: the anchored-gauge series has a + b/r tails whose ringed
jr = r J diverges at scri, so a series init is NOT representable; with
F(u0) ~ 0 the zero slice is exact)."""
function cheb_evolve_J_sat_bc(g::ChebRegGrid{T}, u0::T, u1::T, du::T,
                              jr0::Vector{T}, bcfun;
                              sigma::T = T(1)) where {T}
    ops = ChebSweepOps(g)
    N = length(g.y) - 1
    taupen = sigma*(1/(2g.rwt))/(T(1)/N^2)
    jr = copy(jr0)
    u = u0
    rhs(uu, J) = begin
        b = bcfun(uu)
        (_, _, _, hr) = cheb_sweep_fact(ops, J, (b[2], b[3], b[4], b[5]))
        hr[1] += taupen*(b[1] - J[1])
        hr
    end
    while u < u1 - du/2
        k1 = rhs(u, jr); k2 = rhs(u + du/2, jr .+ du/2 .* k1)
        k3 = rhs(u + du/2, jr .+ du/2 .* k2); k4 = rhs(u + du, jr .+ du .* k3)
        jr .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
    end
    jr
end

"""SAT evolution (sigma = 1; O-N14-1 resolution) using factorized ops."""
function cheb_evolve_J_sat(g::ChebRegGrid{T}, u0::T, u1::T, du::T,
                           X, rc, tau; sigma::T = T(1)) where {T}
    ops = ChebSweepOps(g)
    rwt = g.rwt
    N = length(g.y) - 1
    taupen = sigma*(1/(2rwt))/(T(1)/N^2)
    jr = T[teuk_jr(u0, rofy_c(g, y), X, rc, tau) for y in g.y]
    u = u0
    bcs(uu) = (teuk_qr(uu, rwt, X, rc, tau), teuk_ur(uu, rwt, X, rc, tau),
               teuk_wr(uu, rwt, X, rc, tau), teuk_hr(uu, rwt, X, rc, tau))
    rhs(uu, J) = begin
        (_, _, _, hr) = cheb_sweep_fact(ops, J, bcs(uu))
        hr[1] += taupen*(teuk_jr(uu, rwt, X, rc, tau) - J[1])
        hr
    end
    while u < u1 - du/2
        k1 = rhs(u, jr); k2 = rhs(u + du/2, jr .+ du/2 .* k1)
        k3 = rhs(u + du/2, jr .+ du/2 .* k2); k4 = rhs(u + du, jr .+ du .* k3)
        jr .+= du/6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        u += du
    end
    jr
end
