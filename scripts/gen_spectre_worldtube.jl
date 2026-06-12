# gen_spectre_worldtube.jl — iter 46: generate a SpECTRE-CCE Bondi-Sachs
# worldtube H5 file from the EXACT linear Teukolsky solution in the Bondi
# gauge (ZccmJl machine-derived series, N14 stage A).
#
# Format (reverse-engineered from CceExecutables/Tests/BondiSachsCceR0200.h5
# legends): datasets Beta, DuR, R, W real (cols: Time, Re(0,0), Re(1,0),
# Re(1,1), Im(1,1), Re(2,0), ...; m >= 0 compression); J, DrJ, H, Q, U
# complex spin-weighted (cols: Time, then Re/Im pairs over (l,m), m from
# -l to l, l = 0..16). LMax = 16 -> 290/579 columns.
#
# The data is prescribed DIRECTLY on a constant-Bondi-radius worldtube
# r_B = 41 (R modal (0,0) = 41 sqrt(4pi), DuR = 0): our series ARE the
# Bondi-coordinate fields, so the dataset is an exact linear solution —
# no Cauchy-coordinate gauge shift enters (zr would only appear when
# converting data given on a Cauchy-coordinate sphere).
#
# Mode conversions (t-scalars -> SWSH modal coefficients, m = 0):
#   spin 2 (J, DrJ, H): field = ft s^2,  s^2 = 4 sqrt(2pi/15) 2Y20
#       -> modal (2,0) = 4 sqrt(2pi/15) ft
#   spin 1 (Q, U): field = ft sc,  sc = -sqrt(8pi/15) 1Y20
#       -> modal (2,0) = -sqrt(8pi/15) ft
#   spin 0 const (R): modal (0,0) = R sqrt(4pi)
# t-scalars at the tube (phase u = file time; pulse F(u) peaks at u = rc):
#   Jt = teuk_jr/r, Qt = teuk_qr/r, Ut = teuk_ur/r^2, Wt = teuk_wr/r^2,
#   Ht = teuk_hr/r; drJt = (3/4)(-F5/r - F4/r^2 - F3/r^3 - 3 F2/r^4).
using ZccmJl
using HDF5
using Printf

const rB = 41.0
const X, rc, tau = 1.0e-3, 20.0, 2.0
const LMAX = 16
const NREAL = (LMAX + 1)^2 + 1            # 290
const NCPLX = 2*(LMAX + 1)^2 + 1          # 579
const times = collect(0.0:0.05:45.0)
const nt = length(times)

# real-compressed (2,0) column (0-based after Time: (0,0)=1,(1,0)=2,
# (1,1)=3,4 -> (2,0)=5); complex (2,0): flat index l^2+l+m = 6 -> cols
# 1+2*6 = 13 (Re), 14 (Im) [0-based]
const RE20_REAL = 5
const RE20_CPLX = 13

c2 = 4*sqrt(2pi/15)          # spin-2 modal factor
c1 = -sqrt(8pi/15)           # spin-1 modal factor

function drJt(u, r)
    F(n) = teuk_F(n, u, X, rc, tau)
    3/4*(-F(5)/r - F(4)/r^2 - F(3)/r^3 - 3*F(2)/r^4)
end

real_ds = Dict("Beta" => zeros(nt, NREAL), "DuR" => zeros(nt, NREAL),
               "R" => zeros(nt, NREAL), "W" => zeros(nt, NREAL))
cplx_ds = Dict("J" => zeros(nt, NCPLX), "DrJ" => zeros(nt, NCPLX),
               "H" => zeros(nt, NCPLX), "Q" => zeros(nt, NCPLX),
               "U" => zeros(nt, NCPLX))
for (i, t) in enumerate(times)
    u = t
    for d in values(real_ds); d[i, 1] = t; end
    for d in values(cplx_ds); d[i, 1] = t; end
    real_ds["R"][i, 1 + 1] = rB*sqrt(4pi)        # Re(0,0)
    cplx_ds["J"][i, 1 + RE20_CPLX] = c2*teuk_jr(u, rB, X, rc, tau)/rB
    cplx_ds["DrJ"][i, 1 + RE20_CPLX] = c2*drJt(u, rB)
    cplx_ds["H"][i, 1 + RE20_CPLX] = c2*teuk_hr(u, rB, X, rc, tau)/rB
    cplx_ds["Q"][i, 1 + RE20_CPLX] = c1*teuk_qr(u, rB, X, rc, tau)/rB
    cplx_ds["U"][i, 1 + RE20_CPLX] = c1*teuk_ur(u, rB, X, rc, tau)/rB^2
end
# W is spin-0 but P-shaped: W = Wt * P, P = 3c^2 - 1 = sqrt(16pi/5) Y20
# -> modal (2,0) = sqrt(16pi/5) * Wt
c0P = sqrt(16pi/5)
for (i, t) in enumerate(times)
    real_ds["W"][i, 1 + RE20_REAL] = c0P*teuk_wr(t, rB, X, rc, tau)/rB^2
end

legend_real = ["Time", "Re(0,0)"]
for l in 1:LMAX
    push!(legend_real, "Re($l,0)")
    for m in 1:l
        push!(legend_real, "Re($l,$m)"); push!(legend_real, "Im($l,$m)")
    end
end
legend_cplx = ["Time"]
for l in 0:LMAX, m in -l:l
    push!(legend_cplx, "Re($l,$m)"); push!(legend_cplx, "Im($l,$m)")
end
@assert length(legend_real) == NREAL && length(legend_cplx) == NCPLX

out = "/data/haiyangw/claude/z4c-CMM/results/numerical/spectre_oracle/TeukolskyBondiCceR0041.h5"
mkpath(dirname(out))
rm(out; force=true)
# HDF5.jl reverses dimension order vs C-order readers (Julia is
# column-major): write the transpose so spectre/h5py see (ntimes, ncols)
h5open(out, "w") do f
    for (nm, d) in real_ds
        f[nm*".dat"] = permutedims(d)
        attributes(f[nm*".dat"])["Legend"] = legend_real
    end
    for (nm, d) in cplx_ds
        f[nm*".dat"] = permutedims(d)
        attributes(f[nm*".dat"])["Legend"] = legend_cplx
    end
end
@printf("wrote %s (%d rows; X=%.0e rc=%g tau=%g)\n", out, nt, X, rc, tau)
@printf("peak |J(2,0)| modal = %.6e\n", maximum(abs.(cplx_ds["J"][:, 1+RE20_CPLX])))
