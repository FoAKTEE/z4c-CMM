# Honest unit tests: every ported component must reproduce the admitted
# Python-verified numbers (provenance: scripts/*.py + results/numerical/).
using Test
using ZccmJl

@testset "teukolsky: F^(n) vs Python reference values" begin
    # reference values computed by scripts/teuk_exact_waveform.py F_n
    # (numpy float64) at (n, u, X=1e-5, rc=20, tau=2) — regenerate via
    # python3 -c "...">> see test/python_reference.jl for provenance.
    include("python_reference.jl")
    for (n, u, want) in FREF
        got = teuk_F(n, u, 1.0e-5, 20.0, 2.0)
        @test isapprox(got, want; rtol=1e-13, atol=1e-300)
    end
end

@testset "teukolsky: exact waveform vs Python reference" begin
    include("python_reference.jl")
    for (t, r, want) in WREF
        got = rpsi4_20_exact(t, r, 1.0e-5, 20.0, 2.0)
        @test isapprox(got, want; rtol=1e-13)
    end
end

@testset "waveform BigFloat consistency (arbitrary precision)" begin
    setprecision(BigFloat, 256) do
        t, r = big"54.0", big"36.0"
        hi = rpsi4_20_exact(t, r, big"1.0e-5", big"20.0", big"2.0")
        lo = rpsi4_20_exact(54.0, 36.0, 1.0e-5, 20.0, 2.0)
        @test isapprox(Float64(hi), lo; rtol=1e-13)
    end
end

@testset "psi0 chain (A1 gates, exact rationals)" begin
    pass, rep = psi0_chain_gates()
    @info rep
    @test pass
    # the exact subleading tails from the A1 transcript
    p0 = psi0_lin_series(J_FULL_A1)
    @test get(p0, (1, 6), 0//1) == 27//4
    @test get(p0, (0, 7), 0//1) == 27//4
end

@testset "linseries algebra identities" begin
    a = LinSeries((2, 3) => 5//2)
    @test dr_series(a) == LinSeries((2, 4) => -15//2)
    @test du_series(a) == LinSeries((3, 3) => 5//2)
    @test integrate_r(dr_series(a), "t") == a  # d/dr then integrate = id
    @test_throws ErrorException integrate_r(LinSeries((0, 1) => 1//1), "t")
end

@testset "bondi transform (A2 machine derivation, exact)" begin
    bt = bondi_transform_teukolsky()
    # T1: the three gauge conditions hold identically
    @test isempty(bt.res_rr)
    @test isempty(bt.res_rth)
    @test isempty(bt.res_areal)
    # J must be pure s^2-shape: every key with q = 0 must pair as (1 - c^2)
    # — check J has only shapes {(p=0,q=0), (p=2,q=0)} with opposite signs
    @test all(q == 0 && (p == 0 || p == 2) for ((n, k, p, q), v) in bt.J)
    for ((n, k, p, q), v) in bt.J
        p == 0 && @test get(bt.J, (n, k, 2, 0), big(0)//1) == -v
    end
    # T2: psi0 chain — the A1 adjudication
    p0 = psi0_from_J(bt.J)
    lead00 = get(p0, (2, 5, 0, 0), big(0)//1)
    lead20 = get(p0, (2, 5, 2, 0), big(0)//1)
    @info "psi0 F2/r^5 s^2-coefficient = $(lead00) (c^2 partner $(lead20))"
    @info "J = $(sort(collect(bt.J)))"
    @info "psi0 = $(sort(collect(p0)))"
    @info "beta = $(sort(collect(bt.beta)))"
    @info "U = $(sort(collect(bt.U)))"
    @info "W = $(sort(collect(bt.W)))"
    # A1 AMENDED (iter 30): the complete transformation yields psi0 =
    # -(9/8) s^2 F2(u)/r^5 EXACTLY — equal to the paper datum
    # -sqrt(27pi/10)F2/r^5*2Y20 including sign, with NO subleading tails.
    # A1's +9/8 lead, 27/4 tails, and the "-1 convention map" were artifacts
    # of the incomplete gauge shift (d_theta T dropped).
    @test lead00 == -9//8 && lead20 == 9//8
    @test length(p0) == 2                       # pure datum: no other terms
    # the exact Bondi J of the outgoing Teukolsky solution:
    @test bt.J == AngSeries((4, 1, 0, 0) => big(3)//4, (4, 1, 2, 0) => big(-3)//4,
                            (2, 3, 0, 0) => big(3)//4, (2, 3, 2, 0) => big(-3)//4)
    @test isempty(bt.beta)                      # beta = 0 at linear order
    # T3: falloff sanity — no growing powers anywhere
    for f in (bt.J, bt.beta, bt.U, bt.W)
        @test all(k >= 0 for ((n, k, p, q), v) in f)
    end
end

@testset "hierarchy closure on derived fields (Q1/Q2/Q3)" begin
    hc = hierarchy_closure()
    @info "Q-eq residual: $(sort(collect(hc.q_res)))"
    @info "W-eq residual: $(sort(collect(hc.w_res)))"
    @info hc.h_report
    @test isempty(hc.q_res)            # Q1
    @test isempty(hc.w_res)            # Q2
    @test hc.h_consistent              # Q3
end

@testset "bondi solver: S1 hierarchy sweep convergence" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    u0 = 22.0
    errs = Float64[]
    for ny in (65, 129, 257)
        g = BondiGrid(41.0, ny)
        Jt = [teuk_Jt(u0, r, X, rc, tau) for r in g.r]
        bc = (teuk_Qt(u0, 41.0, X, rc, tau), teuk_Ut(u0, 41.0, X, rc, tau),
              teuk_Wt(u0, 41.0, X, rc, tau), teuk_Ht(u0, 41.0, X, rc, tau))
        Qt, Ut, Wt, Ht = hierarchy_sweep!(g, Jt, bc)
        inner = g.y .< 0.0   # inner half (worldtube side); scri-regularized
                             # variables are the production successor
        eq = maximum(abs.(Qt[inner] .- [teuk_Qt(u0, r, X, rc, tau) for r in g.r[inner]]))
        eu = maximum(abs.(Ut[inner] .- [teuk_Ut(u0, r, X, rc, tau) for r in g.r[inner]]))
        eh = maximum(abs.(Ht[inner] .- [teuk_Ht(u0, r, X, rc, tau) for r in g.r[inner]]))
        push!(errs, max(eq/abs(bc[1]), eu/abs(bc[2]), eh/abs(bc[4])))
    end
    p12 = log2(errs[1]/errs[2]); p23 = log2(errs[2]/errs[3])
    @info "S1 sweep rel errors: $(errs); orders $(p12), $(p23)"
    @test errs[3] < 1e-5
    @test p23 > 3.0
end

@testset "bondi solver: S2 evolution + S3 psi0 gates" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    g = BondiGrid(41.0, 129)
    u0, u1 = 14.0, 26.0
    errs = Float64[]
    for du in (0.1, 0.05)
        Jt = evolve_J(g, u0, u1, du, X, rc, tau)
        ref = [teuk_Jt(u1, r, X, rc, tau) for r in g.r]
        inner = g.y .< 0.0
        push!(errs, maximum(abs.(Jt[inner] .- ref[inner]))/maximum(abs.(ref[inner])))
    end
    @info "S2 evolution rel errors (du=0.1, 0.05): $(errs)"
    @test errs[2] < 2e-4
    # S3a: the psi0 extraction operator on EXACT slices (u = 20, F2 peak)
    g2 = BondiGrid(41.0, 129)
    Jt20 = [teuk_Jt(20.0, r, X, rc, tau) for r in g2.r]
    p0x = psi0_worldtube(g2, Jt20)
    want20 = -9/8*teuk_F(2, 20.0, X, rc, tau)/41.0^5
    relx = abs(p0x - want20)/abs(want20)
    @info "S3a extraction on exact slice: rel = $(relx)"
    @test relx < 1e-7
    # S3b KNOWN-BROKEN (O-N14-1): psi0 from the EVOLVED field is limited by
    # an edge boundary-layer mode of the marching scheme (error GROWS with
    # refinement: rel 4.1 at ny=257 -> 15.0 at ny=513; invisible in inner
    # norms — S2 passes — fatal for edge second derivatives). Candidate
    # fixes ledgered: edge dissipation, (1-y)-regularized variables
    # (eq:*numeric forms), hierarchy-based psi0 evaluation.
    Jte = evolve_J(g, u0, 20.0, 0.05, X, rc, tau)
    p0e = psi0_worldtube(g, Jte)
    @test_broken abs(p0e - want20)/abs(want20) < 5e-3
end
