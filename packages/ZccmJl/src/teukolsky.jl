# teukolsky.jl — port of scripts/teuk_exact_waveform.py (N12, ledger-solid).
# F(u) = X exp(-((u-rc)/tau)^2); F^(n) via physicists' Hermite polynomials.
# Exact finite-radius waveform (AthenaK extraction convention):
#   rpsi4_20(t, r) = sqrt(6pi/5) [ F6 + 2F5/r + 3F4/r^2 + 3F3/r^3
#                    + (3/2)F2/r^4 ](t-r)  -  sqrt(6pi/5)(3/2)F2(t+r)/r^4
# (the N12 verifier's exact table: ret coeffs/(-sqrt(6pi/5)) =
#  {6:-1, 5:-2, 4:-3, 3:-3, 2:-3/2}; adv {2: +3/2}).

export teuk_F, rpsi4_20_exact

const HERMITE = (
    s -> one(s),
    s -> 2s,
    s -> 4s^2 - 2,
    s -> s*(8s^2 - 12),
    s -> (16s^2 - 48)*s^2 + 12,
    s -> s*((32s^2 - 160)*s^2 + 120),
    s -> ((64s^2 - 480)*s^2 + 720)*s^2 - 120,
    s -> s*(((128s^2 - 1344)*s^2 + 3360)*s^2 - 1680),
    s -> (((256s^2 - 3584)*s^2 + 13440)*s^2 - 13440)*s^2 + 1680,
)

"""F^(n)(u) for the Gaussian pulse; works for Float64 and BigFloat."""
function teuk_F(n::Integer, u::T, X::T, rc::T, tau::T) where {T<:Real}
    s = (u - rc)/tau
    return X*(-1/tau)^n*HERMITE[n+1](s)*exp(-s^2)
end

# exact rational table (N12; coefficients of -sqrt(6pi/5))
const RET_TABLE = ((6, 0, -1//1), (5, 1, -2//1), (4, 2, -3//1),
                   (3, 3, -3//1), (2, 4, -3//2))
const ADV_TABLE = ((2, 4, 3//2),)

"""Exact linear r*psi4(2,0) at finite radius, AthenaK convention (N12)."""
function rpsi4_20_exact(t::T, r::T, X::T, rc::T, tau::T) where {T<:Real}
    pref = -sqrt(T(6)*T(pi)/5)
    out = zero(T)
    for (n, k, c) in RET_TABLE
        out += pref*T(c)*teuk_F(n, t - r, X, rc, tau)/r^k
    end
    for (n, k, c) in ADV_TABLE
        out += pref*T(c)*teuk_F(n, t + r, X, rc, tau)/r^k
    end
    return out
end
