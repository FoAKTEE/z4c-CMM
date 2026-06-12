# worldtubemap.jl — the LIVE worldtube map (N14 stage B, iter 41-42):
# exact local closed forms in the consistent worldtube-anchored gauge
# (derivation: test/derive_wtmap_v2.jl; all coefficients verified exact on
# the Teukolsky case at r_wt = 41; the r_wt powers follow the analytic
# structure — beta fully hand-verified channel-by-channel).
# Inputs: l=2 m=0 scalarized Cauchy data at r = r_wt (the data contract):
#   hTT      = (h_thth^ - h_phph^)/(2 s^2)  [s^2-scalar]
#   dt_hTT   = du of hTT
#   hrr      = h_rr P-scalar (P = 3c^2 - 1)
#   trace    = (h_thth^/r^2 + h_phph^/(r^2 s^2)) P-scalar
#   dr_trace, dth_trace_sc (= sc-scalar of d_theta trace = -6x trace-scalar),
#   dth_dr_trace_sc, hrth (sc-scalar of coordinate h_{r theta}),
#   dt_hrth, ethb_hrth_P (P-scalar of (d_theta + cot) h_rth = trace-free part)
# Outputs: ringed worldtube BC scalars (jr, qr, ur, wr, hr) + beta.
export worldtube_map

function worldtube_map(rwt::T; hTT::T, dt_hTT::T, hrr::T, trace::T,
                       dr_trace::T, dth_trace_sc::T, dth_dr_trace_sc::T,
                       hrth::T, dt_hrth::T, ethb_hrth_P::T) where {T<:Real}
    Jw = hTT
    Hw = dt_hTT
    Uw = -dth_trace_sc/(4rwt)
    Qw = -dth_trace_sc/4 - rwt*dth_dr_trace_sc/4 - 3hrth/rwt + dt_hrth
    Bw = hrr/4 - trace/8 - rwt*dr_trace/8 + ethb_hrth_P/(4rwt)
    Ww = 3hrr/(4rwt) - dr_trace/4
    (jr = rwt*Jw, qr = rwt*Qw, ur = rwt^2*Uw, wr = rwt^2*Ww,
     hr = rwt*Hw, beta = Bw)
end

"evaluate the s^2/sc/P-scalar of an AngSeries at (u-channel F^n, r)"
function series_scalar(a::AngSeries, kind::Symbol, u::T, r::T,
                       X, rc, tau) where {T<:Real}
    out = zero(T)
    for ((n, k, p, q), v) in a
        vv = v
        if kind == :s2 && q == 0 && p == 0
        elseif kind == :sc && q == 1 && p == 1
        elseif kind == :P && q == 0 && p == 2
            vv = v//3
        else
            continue
        end
        out += T(vv)*teuk_F(n, u, T(X), T(rc), T(tau))/r^k
    end
    out
end
export series_scalar
