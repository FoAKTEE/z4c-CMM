# check_n14_mode6.jl — iter 45 gates for the N14 datum-fidelity test
# (2308.10361 Sec V.C characteristic-pulse injection through the live
# native solver, mode 6). Usage:
#   julia --project=<repo>/packages/ZccmJl scripts/check_n14_mode6.jl
#
# The quiescent-BC ZccmJl reference is the MIRROR (CCE) solution; the
# live run is the TRANSPARENT (CCM) solution — they coincide only while
# the pulse has not yet interacted with the tube, and the O(1) deviation
# afterwards is the transparency physics, not an error.
# Gate A (early-window datum fidelity): dev/peak < 1e-6 at t <= 19,
#   < 1e-4 at t <= 28 (measured 7e-8 / 2.4e-5).
# Gate B (transmission burst): rpsi4(2,0) at r=36 peaks at t2 ~ 161
#   (= center-arrival 79 + Cauchy crossing 82); early times quiet.
# Gate C (loop closure): the exit burst re-crosses r_B at t ~ 160-180:
#   live |psi0| >> mirror reference there.
# Gate L (self-convergence of the coupled system): full-vs-half
#   resolution rpsi4 at fixed tube radius 39 (mode6_r39 vs mode6_half).
# Gate D (control): no-datum run shows no burst.
using DelimitedFiles
using Printf

const base = joinpath(@__DIR__, "..", "results", "numerical")
const rB = 41.0

# ---- reference ----
ref = readdlm(joinpath(base, "n14_pulse_ref.csv"), ','; comments=true,
              comment_char='#')
ref_t = Float64.(ref[:, 2]); ref_v = Float64.(ref[:, 3])
pk_ref = maximum(abs.(ref_v))
lin(t) = begin     # linear interp of the reference at time t
    t <= ref_t[1] && return ref_v[1]
    t >= ref_t[end] && return ref_v[end]
    i = searchsortedlast(ref_t, t)
    f = (t - ref_t[i])/(ref_t[i+1] - ref_t[i])
    (1 - f)*ref_v[i] + f*ref_v[i+1]
end

# ---- diag series from the live run ----
lines = filter(l -> startswith(l, "ccm5-diag"),
               readlines(joinpath(base, "n14_native", "mode6", "run.log")))
println("mode6: $(length(lines)) diag samples")
ts = Float64[]; psis = Float64[]; jrs = Float64[]; hTTs = Float64[]
for l in lines
    kv = Dict{Symbol,Float64}()
    for tok in split(l)[2:end]
        k, v = split(tok, "=")
        kv[Symbol(k)] = parse(Float64, v)
    end
    push!(ts, kv[:t]); push!(psis, kv[:psi0]); push!(jrs, kv[:jr])
    push!(hTTs, kv[:hTT])
end

# Gate A: early-window datum fidelity (live == mirror before interaction)
for (tcut, gate) in ((19.0, 1e-6), (28.5, 1e-4))
    mask = ts .<= tcut
    errA = maximum(abs.(psis[mask] .- lin.(ts[mask])))/pk_ref
    @printf("Gate A: dev/peak (t <= %.1f) = %.3e  [gate < %.0e] %s\n",
            tcut, errA, gate, errA < gate ? "PASS" : "FAIL")
end
errAall = maximum(abs.(psis .- lin.(ts)))
@printf("        full-run live-vs-mirror deviation = %.3e of ref peak (transparency physics, expected O(1))\n", errAall/pk_ref)

# Gate B: transmission burst in rpsi4(2,0) at r = 36
function wf6(dir)
    d = readdlm(joinpath(base, "n14_native", dir, "waveforms",
                         "rpsi4_real_0036.txt"); comments=true, comment_char='#')
    (Float64.(d[:, 1]), Float64.(d[:, 4]))
end
wt, w = wf6("mode6")
pkw, ipkw = findmax(abs.(w))
okB = abs(wt[ipkw] - 161.0) < 5.0
@printf("Gate B: rpsi4(2,0) r=36 peak %.3e at t = %.2f [gate |t-161|<5] %s\n",
        pkw, wt[ipkw], okB ? "PASS" : "FAIL")
early = maximum(abs.(w[wt .< 45.0]))
@printf("        early-time (t<45) max = %.3e of peak %s\n", early/pkw,
        early/pkw < 1e-3 ? "PASS" : "FAIL")

# Gate C: loop closure — exit burst re-crossing r_B in the live psi0
late = (ts .>= 160.0)
ratioC = maximum(abs.(psis[late]))/max(maximum(abs.(lin.(ts[late]))), 1e-30)
@printf("Gate C: exit-window live|psi0|/mirror = %.1f [gate > 10] %s\n",
        ratioC, ratioC > 10 ? "PASS" : "FAIL")

# Gate L: coupled-system self-convergence (fixed tube radius 39). The
# waveform labels are offset by (dt_full - dt_half)/... between
# resolutions (label = time + dt): interpolate the half-res series onto
# the full-res times (linear; the label offset 0.125 << feature widths).
tfull, wfull = wf6("mode6_r39")
thalf, whalf = wf6("mode6_half")
linh(t) = begin
    t <= thalf[1] && return whalf[1]
    t >= thalf[end] && return whalf[end]
    i = searchsortedlast(thalf, t)
    f = (t - thalf[i])/(thalf[i+1] - thalf[i])
    (1 - f)*whalf[i] + f*whalf[i+1]
end
inwin = (tfull .>= thalf[1]) .& (tfull .<= thalf[end])
pkL = maximum(abs.(wfull))
dL = maximum(abs.(wfull[inwin] .- linh.(tfull[inwin])))
@printf("Gate L: full-vs-half rpsi4 max|diff| = %.3e = %.3e of peak (~= the half-res truncation error)\n", dL, dL/pkL)

# Gate D: control (no datum) — the burst must vanish
tc, wc = wf6("mode6_ctrl")
pkc = maximum(abs.(wc))
@printf("Gate D: control peak |rpsi4| = %.3e = %.3e of the mode-6 peak %s\n",
        pkc, pkc/pkw, pkc/pkw < 1e-8 ? "PASS" : "FAIL")

# context line
pkj, ipkj = findmax(abs.(jrs))
@printf("info: tube |jr| peak %.3e at t = %.2f; hTT peak %.3e\n",
        pkj, ts[ipkj], maximum(abs.(hTTs)))
