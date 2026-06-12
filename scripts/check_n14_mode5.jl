# check_n14_mode5.jl — iter 44 gates for N14 stage C part 2 (live CCM).
# Usage: julia --project=packages/ZccmJl scripts/check_n14_mode5.jl
#
# (G3) waveform gates: mode4_v2 (phase-fixed analytic stub) and mode5
#      (live worldtube data) vs mode3 (exact self-datum), rel to peak.
# (G4) datum-level gates from the ccm5-diag lines: sampled contract
#      scalars, mapped BCs, and the served psi0 datum vs the analytic
#      Teukolsky references (rwt = 40, r_B = 41; retarded phase t - r).
using ZccmJl
using DelimitedFiles
using Printf

const base = joinpath(@__DIR__, "..", "results", "numerical", "n14_native")
const X, rc, tau = 1.0e-5, 20.0, 2.0
const rwt, rB = 40.0, 41.0

# ---- G3: waveforms -------------------------------------------------------
function wf(dir)
    d = readdlm(joinpath(base, dir, "waveforms", "rpsi4_real_0036.txt");
                comments=true, comment_char='#')
    (d[:, 1], d[:, 4])           # column 4 = the (l=2, m=0) mode
end
t3, w3 = wf("mode3")
peak = maximum(abs.(w3))
for dir in ("mode4_v2", "mode5")
    td, wd = wf(dir)
    @assert td == t3
    md = maximum(abs.(wd .- w3))
    @printf("G3 %-9s vs mode3: max|diff| = %.3e  rel-to-peak = %.3e\n",
            dir, md, md/peak)
end
t4, w4 = wf("mode4_v2"); t5, w5 = wf("mode5")
m45 = maximum(abs.(w5 .- w4))
@printf("G3 mode5 vs mode4_v2: max|diff| = %.3e  rel-to-peak = %.3e\n",
        m45, m45/peak)

# ---- G4: datum-level diagnostics ----------------------------------------
# analytic live scalars at the tube (Cauchy time t; retarded phase t - rwt)
function scalars_exact(t)
    u = t - rwt
    A(n) = teuk_F(n, u, X, rc, tau)
    Aser = 3*(A(2)/rwt^3 + 3A(1)/rwt^4 + 3A(0)/rwt^5)
    drA = 3*(-3A(2)/rwt^4 - 12A(1)/rwt^5 - 15A(0)/rwt^6)
    Bser = -(A(3)/rwt^2 + 3A(2)/rwt^3 + 6A(1)/rwt^4 + 6A(0)/rwt^5)
    Cser = (A(4)/rwt + 2A(3)/rwt^2 + 9A(2)/rwt^3 + 21A(1)/rwt^4 + 21A(0)/rwt^5)/4
    dtB = -(A(4)/rwt^2 + 3A(3)/rwt^3 + 6A(2)/rwt^4 + 6A(1)/rwt^5)
    dtC = (A(5)/rwt + 2A(4)/rwt^2 + 9A(3)/rwt^3 + 21A(2)/rwt^4 + 21A(1)/rwt^5)/4
    dtA = 3*(A(3)/rwt^3 + 3A(2)/rwt^4 + 3A(1)/rwt^5)
    (hTT=(3/2)*(2Cser - Aser), dt_hTT=(3/2)*(2dtC - dtA), hrr=Aser,
     trace=-Aser, dr_trace=-drA, hrth=-3*Bser*rwt, dt_hrth=-3*dtB*rwt)
end
function bc_exact(t)
    s = scalars_exact(t)
    worldtube_map(rwt; hTT=s.hTT, dt_hTT=s.dt_hTT, hrr=s.hrr, trace=s.trace,
        dr_trace=s.dr_trace, dth_trace_sc=-6*s.trace,
        dth_dr_trace_sc=-6*s.dr_trace, hrth=s.hrth, dt_hrth=s.dt_hrth,
        ethb_hrth_P=s.hrth)
end
psi0_exact(t) = -9/8*teuk_F(2, t - rB, X, rc, tau)/rB^5

lines = filter(l -> startswith(l, "ccm5-diag"),
               readlines(joinpath(base, "mode5", "run.log")))
println("G4: $(length(lines)) diag samples")
names = (:hTT, :dthTT, :hrr, :tr, :drtr, :hrth, :dthrth,
         :jr, :qr, :ur, :wr, :hr, :beta, :psi0)
function refs_at(t)
    s = scalars_exact(t)
    m = bc_exact(t)
    Dict(:hTT => s.hTT, :dthTT => s.dt_hTT, :hrr => s.hrr,
         :tr => s.trace, :drtr => s.dr_trace, :hrth => s.hrth,
         :dthrth => s.dt_hrth, :jr => m.jr, :qr => m.qr,
         :ur => m.ur, :wr => m.wr, :hr => m.hr, :beta => m.beta,
         :psi0 => psi0_exact(t))
end
samples = map(lines) do l
    kv = Dict{Symbol,Float64}()
    for tok in split(l)[2:end]
        k, v = split(tok, "=")
        kv[Symbol(k)] = parse(Float64, v)
    end
    kv
end
# per-channel scale = the analytic pulse peak of THAT channel over the run
scale = Dict(nm => maximum(abs(refs_at(kv[:t])[nm]) for kv in samples)
             for nm in names)
worst = Dict(nm => 0.0 for nm in names)
worst_t = Dict(nm => 0.0 for nm in names)
for kv in samples
    r = refs_at(kv[:t])
    for nm in names
        scale[nm] == 0.0 && continue
        rel = abs(kv[nm] - r[nm])/scale[nm]
        if rel > worst[nm]
            worst[nm] = rel; worst_t[nm] = kv[:t]
        end
    end
end
for nm in names
    @printf("G4 %-7s worst err/peak = %.3e (t = %.1f, peak scale %.3e)\n",
            nm, worst[nm], worst_t[nm], scale[nm])
end
