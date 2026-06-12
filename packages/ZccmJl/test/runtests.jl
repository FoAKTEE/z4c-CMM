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
