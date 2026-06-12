# bonditransform.jl — the COMPLETE linear Bondi-like transformation of the
# corrected Teukolsky solution (N14 stage A2; adjudicates the iter-27 A1
# re-examination). Exact angular-shape-indexed series algebra: keys
# (n, k, p, q) = F^(n)(u) * r^(-k) * c^p * s^q with c = cos(th), s = sin(th),
# s^2 reduced to 1 - c^2 so q in {0, 1}.
#
# Flat-Bondi background (u, r, th, ph): gbar = -du^2 - 2 du dr + r^2 dOmega^2.
# Gauge delta g_mn = h_mn + zeta^a d_a gbar_mn + gbar_an d_m zeta^a
#                    + gbar_ma d_n zeta^a, zeta = (zu, zr, zth, 0):
#   dg_rr  = h_rr - 2 dr zu                          = 0   -> zu
#   dg_rth = h_rth + r^2 dr zth - dth zu             = 0   -> zth
#   dg_thth/r^2 + dg_phph/(r^2 s^2)
#          = trace(h_ang) + 4 zr/r + 2 dth zth + 2 (c/s) zth = 0 -> zr
# (all explicit; scri falloff = integration from infinity, no constants).
# Read-off:
#   J    = (dg_thth/r^2 - dg_phph/(r^2 s^2))/2
#        = (h_thth/r^2 - h_phph/(r^2 s^2))/2 + dth zth - (c/s) zth
#   beta = (du zu + dr zu + dr zr)/2          [dg_ur = -2 beta]
#   U    = -dg_uth/r^2 = -(r^2 du zth - dth zu - dth zr)/r^2
#   W    = -(dg_uu + 2 beta)/r,  dg_uu = -2 du zu - 2 du zr

export AngSeries, bondi_transform_teukolsky, ang_dtheta, ang_dr, ang_du,
       ang_cot_mul, psi0_from_J

const AngSeries = Dict{NTuple{4,Int},Rational{BigInt}}

aadd!(a::AngSeries, key, v) = begin
    vv = get(a, key, big(0)//1) + v
    vv == 0 ? delete!(a, key) : (a[key] = vv)
    a
end
aadd(a::AngSeries, b::AngSeries) = begin
    out = copy(a)
    for (k, v) in b
        aadd!(out, k, v)
    end
    out
end
ascale(a::AngSeries, f) = begin
    ff = Rational{BigInt}(f)
    AngSeries(k => v*ff for (k, v) in a if v*ff != 0)
end

"""reduce s^q with q >= 2 via s^2 = 1 - c^2 (maintains q in {0,1})."""
function s2reduce(key::NTuple{4,Int}, v::Rational{BigInt})
    n, k, p, q = key
    q < 2 && return AngSeries(key => v)
    out = AngSeries()
    for (kk, vv) in s2reduce((n, k, p, q-2), v)
        nn, kk2, pp, qq = kk
        aadd!(out, (nn, kk2, pp, qq), vv)
        aadd!(out, (nn, kk2, pp+2, qq), -vv)
    end
    out
end

ang_du(a::AngSeries) = AngSeries((n+1, k, p, q) => v for ((n, k, p, q), v) in a)
ang_dr(a::AngSeries) =
    AngSeries((n, k+1, p, q) => -k*v for ((n, k, p, q), v) in a if k != 0)
ang_rpow(a::AngSeries, m::Int) =
    AngSeries((n, k - m, p, q) => v for ((n, k, p, q), v) in a)

"""d/dtheta on c^p s^q (q in {0,1}); returns reduced series."""
function ang_dtheta(a::AngSeries)
    out = AngSeries()
    for ((n, k, p, q), v) in a
        # d(c^p s^q) = -p c^(p-1) s^(q+1) + q c^(p+1) s^(q-1)
        if p > 0
            for (kk, vv) in s2reduce((n, k, p-1, q+1), -p*v)
                aadd!(out, kk, vv)
            end
        end
        if q > 0
            aadd!(out, (n, k, p+1, q-1), q*v)
        end
    end
    out
end

"""(c/s) * series — requires every term to carry q = 1 (else singular)."""
function ang_cot_mul(a::AngSeries)
    out = AngSeries()
    for ((n, k, p, q), v) in a
        q == 1 || error("cot multiplication on a q=$q term — singular")
        aadd!(out, (n, k, p+1, 0), v)
    end
    out
end

"""integrate in r with scri falloff (no constants); log channel forbidden."""
function ang_integrate_r(a::AngSeries, label::String)
    out = AngSeries()
    for ((n, k, p, q), v) in a
        k == 1 && error("log channel in $label at n=$n p=$p q=$q")
        out[(n, k-1, p, q)] = v//(1 - k)
    end
    out
end

# ---- the corrected Teukolsky solution in flat-Bondi coordinates ----------
const _A = ((2, 3, 3), (1, 4, 9), (0, 5, 9))            # (n, k, coeff): A
const _BR = ((3, 1, -1), (2, 2, -3), (1, 3, -6), (0, 4, -6))  # B*r
const _C = ((4, 1, 1//4), (3, 2, 1//2), (2, 3, 9//4),
            (1, 4, 21//4), (0, 5, 21//4))


function shape_mul(tab, p::Int, q::Int, fac)::AngSeries
    out = AngSeries()
    for (n, k, cf) in tab
        for (kk, vv) in s2reduce((n, k, p, q), Rational{BigInt}(cf)*Rational{BigInt}(fac))
            aadd!(out, kk, vv)
        end
    end
    out
end


"""Run the full transformation; returns a NamedTuple with J, beta, U, W,
zu, zth, zr and the verification residuals."""
function bondi_transform_teukolsky()
    # h_rr = A (3c^2 - 1)
    h_rr = aadd(shape_mul(_A, 2, 0, 3), shape_mul(_A, 0, 0, -1))
    # h_rth = -3 (B r) s c
    h_rth = shape_mul(_BR, 1, 1, -3)
    # h_thth / r^2 = 3C(1 - c^2) - A
    h_T = aadd(aadd(shape_mul(_C, 0, 0, 3), shape_mul(_C, 2, 0, -3)),
               shape_mul(_A, 0, 0, -1))
    # h_phph / (r^2 s^2) = -3C(1 - c^2) + A(2 - 3c^2)
    h_P = aadd(aadd(shape_mul(_C, 0, 0, -3), shape_mul(_C, 2, 0, 3)),
               aadd(shape_mul(_A, 0, 0, 2), shape_mul(_A, 2, 0, -3)))

    # zu: dr zu = h_rr / 2
    zu = ang_integrate_r(ascale(h_rr, 1//2), "zu")
    # zth: r^2 dr zth = dth zu - h_rth
    zth = ang_integrate_r(
        ang_rpow(aadd(ang_dtheta(zu), ascale(h_rth, -1)), -2), "zth")
    # zr: 4 zr / r = -(trace + 2 dth zth + 2 cot zth)
    trace = aadd(h_T, h_P)
    zr = ang_rpow(ascale(aadd(trace,
                              aadd(ascale(ang_dtheta(zth), 2),
                                   ascale(ang_cot_mul(zth), 2))), -1//4), 1)

    # read-off
    J = aadd(ascale(aadd(h_T, ascale(h_P, -1)), 1//2),
             aadd(ang_dtheta(zth), ascale(ang_cot_mul(zth), -1)))
    beta = ascale(aadd(ang_du(zu), aadd(ang_dr(zu), ang_dr(zr))), 1//2)
    U = ascale(aadd(ang_du(zth),
                    ang_rpow(ascale(aadd(ang_dtheta(zu), ang_dtheta(zr)),
                                    -1), -2)), -1)
    dg_uu = ascale(aadd(ang_du(zu), ang_du(zr)), -2)
    W = ang_rpow(ascale(aadd(dg_uu, ascale(beta, 2)), -1), -1)

    # residual checks (the conditions themselves, recomputed)
    res_rr = aadd(h_rr, ascale(ang_dr(zu), -2))
    res_rth = aadd(h_rth, aadd(ang_rpow(ang_dr(zth), 2),
                               ascale(ang_dtheta(zu), -1)))
    res_areal = aadd(trace, aadd(ascale(ang_rpow(zr, -1), 4),
                                 aadd(ascale(ang_dtheta(zth), 2),
                                      ascale(ang_cot_mul(zth), 2))))
    (J=J, beta=beta, U=U, W=W, zu=zu, zth=zth, zr=zr,
     res_rr=res_rr, res_rth=res_rth, res_areal=res_areal)
end

"""psi0 = -dr J/(2r) - dr^2 J / 4 on an angular series."""
function psi0_from_J(J::AngSeries)
    dJ = ang_dr(J)
    aadd(ascale(ang_rpow(dJ, -1), -1//2), ascale(ang_dr(dJ), -1//4))
end
