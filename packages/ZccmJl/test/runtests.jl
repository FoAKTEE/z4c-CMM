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
