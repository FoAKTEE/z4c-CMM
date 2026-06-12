# check_spectre_oracle.jl — iter 46: gate the SpECTRE CharacteristicExtract
# cross-oracle output (TeukolskyOracleReduction.h5) against the exact
# linear Teukolsky scri content (ZccmJl/N12 machine-derived).
#
# The worldtube input (gen_spectre_worldtube.jl) is the EXACT Bondi-gauge
# linear solution at r_B = 41, X = 1e-3, rc = 20, tau = 2, only (2,0).
# Exact retarded-time shapes at scri (t-scalars):
#   lim r J     = (3/4) F4(u)            -> Strain ~ F4
#   News        ~ F5
#   lim r psi4  ~ F6   (N12: paper-convention amplitude sqrt(6pi/5))
#   lim r^5 psi0 = -(9/8) F2(u)          (the datum channel itself!)
# SpECTRE's inertial retarded time carries a constant offset for linear
# data: fit (A, c) in  data(t) ~ A * Fn(t - c)  per quantity; gates:
#   - shape residual after fit < 1e-2 of peak (News/Psi4/Strain),
#     < 5e-2 for Psi0 (spectre's own tolerance class for Psi0 is 1e4x);
#   - fitted time shifts agree across quantities (< 0.2);
#   - fitted amplitudes reported as convention ratios vs the naive
#     predictions (documented, not gated: spectre normalization factors).
using ZccmJl
using HDF5
using Printf

const base = "/data/haiyangw/claude/z4c-CMM/results/numerical/spectre_oracle"
const X, rc, tau = 1.0e-3, 20.0, 2.0
const c2 = 4*sqrt(2pi/15)              # s^2 -> (2,0) modal factor

f = h5open(joinpath(base, "TeukolskyOracleReduction.h5"), "r")
grp = f["SpectreR0041.cce"]
println("oracle output subfiles: ", keys(grp))

function modal20(name)
    d = read(grp[name])                # HDF5.jl: (ncols, ntimes) col-major view
    dm = size(d, 1) < size(d, 2) ? d : permutedims(d)   # -> (ncols, ntimes)
    t = dm[1, :]
    re20 = dm[14, :]                   # 1-based: col 14 = Re(2,0) (flat idx 6)
    (t, re20)
end

function fitshape(t, v, n)
    pk = maximum(abs.(v))
    best = (1e30, 0.0, 0.0)
    for c in -10.0:0.01:15.0
        F = [teuk_F(n, tt - c, X, rc, tau) for tt in t]
        den = sum(F .* F)
        den == 0 && continue
        A = sum(F .* v)/den
        r = maximum(abs.(v .- A .* F))/pk
        r < best[1] && (best = (r, A, c))
    end
    best
end

quants = (("Strain", 4, c2*3/4), ("News", 5, c2*3/8),
          ("Psi4", 6, -c2*0.75), ("Psi0", 2, -c2*9/8))
shifts = Float64[]
for (nm, n, Apred) in quants
    t, v = modal20(nm)
    pk = maximum(abs.(v))
    res, A, c = fitshape(t, v, n)
    push!(shifts, c)
    gate = nm == "Psi0" ? 5e-2 : 1e-2
    @printf("%-6s peak %.3e | fit F%d: residual %.3e of peak [gate < %.0e] %s | shift c = %+.2f | A = %+.6e (A/naive = %+.3f)\n",
            nm, pk, n, res, gate, res < gate ? "PASS" : "FAIL", c, A, A/Apred)
end
ds = maximum(shifts) - minimum(shifts)
@printf("shift consistency: max spread = %.3f [gate < 0.2] %s\n", ds,
        ds < 0.2 ? "PASS" : "FAIL")
close(f)
