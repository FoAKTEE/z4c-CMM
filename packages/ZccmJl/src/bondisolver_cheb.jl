# bondisolver_cheb.jl — O-N14-1 (iter 36): Chebyshev spectral radial scheme.
# The iter-35 conditioning analysis: psi0_wt/J ~ 1e-8 (peeling), so the
# field needs ~1e-11 relative accuracy — spectral or bust (FD-4 ceiling
# 3.7e-7). Chebyshev-Gauss-Lobatto collocation in y; pole equations as
# GLOBAL dense solves (the scri row enforces regularity automatically since
# (1-y) = 0 there); integrals via the exact Chebyshev antiderivative.

using LinearAlgebra

export ChebRegGrid, cheb_sweep, cheb_evolve_J, cheb_psi0_hierarchy

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
