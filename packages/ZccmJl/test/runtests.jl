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

@testset "regularized solver (O-N14-1): sweep + evolution + psi0" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    u0 = 22.0
    # R1: regularized sweep vs closed forms (FULL grid now, scri included)
    errs = Float64[]
    for ny in (129, 257, 513)
        g = RegGrid(41.0, ny)
        jr = [teuk_jr(u0, rofy(g, y), X, rc, tau) for y in g.y]
        bc = (teuk_qr(u0, 41.0, X, rc, tau), teuk_ur(u0, 41.0, X, rc, tau),
              teuk_wr(u0, 41.0, X, rc, tau), teuk_hr(u0, 41.0, X, rc, tau))
        qr, ur, wr, hr = reg_sweep(g, jr, bc)
        scale = maximum(abs.(jr))
        e = max(maximum(abs.(qr .- [teuk_qr(u0, rofy(g, y), X, rc, tau) for y in g.y])),
                maximum(abs.(ur .- [teuk_ur(u0, rofy(g, y), X, rc, tau) for y in g.y])),
                maximum(abs.(hr .- [teuk_hr(u0, rofy(g, y), X, rc, tau) for y in g.y])))/scale
        push!(errs, e)
    end
    @info "R1 regularized sweep rel errors: $(errs); orders $(log2(errs[1]/errs[2])), $(log2(errs[2]/errs[3]))"
    @test errs[3] < 1e-4
    @test log2(errs[2]/errs[3]) > 1.6     # trapezoid pole march: 2nd order
    # R2/R3: evolution to the F2 peak and psi0 at the worldtube
    want20 = -9/8*teuk_F(2, 20.0, X, rc, tau)/41.0^5
    rels = Float64[]
    for (ny, du) in ((257, 0.05), (513, 0.025))
        g = RegGrid(41.0, ny)
        jr = reg_evolve_J(g, 14.0, 20.0, du, X, rc, tau)
        p0 = reg_psi0_worldtube(g, jr)
        push!(rels, abs(p0 - want20)/abs(want20))
    end
    @info "R3 psi0_wt rel errors (refining): $(rels)"
    # O-N14-1 REMAINS OPEN (iter 33): the edge mode survives the
    # regularized variables, the integration-by-parts H (no D^2
    # composition), the inward-marched U, and the mixed low-order edge
    # closure (reduced: 11.8 -> 5.0 at the fine point, still O(1) and
    # refinement-growing). Sharpened root cause: Dirichlet NODE-PINNING of
    # the worldtube BC is not an SBP-SAT imposition — the discrete energy
    # estimate fails at the boundary row. Next: SAT penalty (iter 34).
    @test_broken rels[1] < 5e-3
    @test_broken rels[2] < rels[1]
end

@testset "spectral solver (O-N14-1 DISCHARGED, iter 37)" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    want20 = -9/8*teuk_F(2, 20.0, X, rc, tau)/41.0^5
    # exact-slice sweep: machine-exact (ringed Teukolsky fields are
    # low-degree polynomials in y)
    g = ChebRegGrid(41.0, 33)
    jr = [teuk_jr(20.0, rofy_c(g, y), X, rc, tau) for y in g.y]
    bc = (teuk_qr(20.0, 41.0, X, rc, tau), teuk_ur(20.0, 41.0, X, rc, tau),
          teuk_wr(20.0, 41.0, X, rc, tau), teuk_hr(20.0, 41.0, X, rc, tau))
    qr, ur, wr, hr = cheb_sweep(g, jr, bc)
    @test maximum(abs.(qr .- [teuk_qr(20.0, rofy_c(g, y), X, rc, tau) for y in g.y])) < 1e-18
    @test abs(cheb_psi0_hierarchy(g, jr, bc) - want20)/abs(want20) < 1e-6
    # Float64 evolved psi0: documented conditioning floor (~0.05; the datum
    # is r^-5-peeling-tiny: psi0/J ~ 1e-8 — iter-35 analysis)
    g65 = ChebRegGrid(41.0, 65)
    jre = cheb_evolve_J_sat(g65, 16.0, 20.0, 0.00078125, X, rc, tau)
    bc20 = (teuk_qr(20.0, 41.0, X, rc, tau), teuk_ur(20.0, 41.0, X, rc, tau),
            teuk_wr(20.0, 41.0, X, rc, tau), teuk_hr(20.0, 41.0, X, rc, tau))
    relF = abs(cheb_psi0_hierarchy(g65, jre, bc20) - want20)/abs(want20)
    @info "Float64 evolved psi0 rel = $(relF) (conditioning floor ~0.05)"
    @test relF < 0.2
    # BigFloat: the O-N14-1 gate (rel < 5e-3, measured 3.3e-4 at du=3.9e-4;
    # 5.3e-3 at du=7.8e-4 — use the passing point; ~80 s)
    setprecision(BigFloat, 128) do
        Xb, rcb, taub = big"1.0e-5", big"20.0", big"2.0"
        wantb = -big(9)/8*teuk_F(2, big"20.0", Xb, rcb, taub)/big"41.0"^5
        gb = ChebRegGrid(big"41.0", 49)
        jrb = cheb_evolve_J_sat(gb, big"16.0", big"20.0", big"0.000390625",
                                Xb, rcb, taub)
        bcb = (teuk_qr(big"20.0", big"41.0", Xb, rcb, taub),
               teuk_ur(big"20.0", big"41.0", Xb, rcb, taub),
               teuk_wr(big"20.0", big"41.0", Xb, rcb, taub),
               teuk_hr(big"20.0", big"41.0", Xb, rcb, taub))
        relB = Float64(abs(cheb_psi0_hierarchy(gb, jrb, bcb) - wantb)/abs(wantb))
        @info "BigFloat evolved psi0 rel = $(relB) (gate 5e-3)"
        @test relB < 5e-3
    end
end

@testset "beta-corrected hierarchy (machine-derived, exact)" begin
    # derivation: test/derive_beta_hierarchy.jl (iter 44). The anchored
    # series satisfies the full linear system; the residuals against the
    # beta = 0 hierarchy are EXACT beta-linear sources.
    ZZ = ZccmJl
    J, U, W, B = J_ANCH, U_ANCH, W_ANCH, BETA_ANCH
    Q = ZZ.ang_rpow(ang_dr(U), 2)
    H = ang_du(J)
    eB = eth_spin(B, 0)
    # Q: dr(r^2 Q) + r^2 dr(ethbar J) = -4 r eth(beta)
    resQ = ZZ.aadd(ZZ.ang_dr(ZZ.ang_rpow(Q, 2)),
                   ZZ.ang_rpow(ZZ.ang_dr(ethbar_spin(J, 2)), 2))
    @test resQ == ZZ.ascale(ZZ.ang_rpow(eB, 1), -4)
    # W: dr(r^2 W) - [(1/2) ethb ethb J + 2 r ethbar U + (r^2/2) dr(ethbar U)]
    #    = -(4/3) ethbar eth beta
    ebbJ = ethbar_spin(ethbar_spin(J, 2), 1)
    sumU = ZZ.ascale(ethbar_spin(U, 1), 2)
    rhsW = ZZ.aadd(ZZ.ascale(ebbJ, 1//2),
                   ZZ.aadd(ZZ.ang_rpow(sumU, 1),
                           ZZ.ascale(ZZ.ang_rpow(ZZ.ang_dr(sumU), 2), 1//4)))
    resW = ZZ.aadd(ZZ.ang_dr(ZZ.ang_rpow(W, 2)), ZZ.ascale(rhsW, -1))
    @test resW == ZZ.ascale(ethbar_spin(eB, 1), -4//3)
    # H (joint refit, valid on plain AND anchored families):
    # dr(rH) = (1/2) r dr^2 J + dr J + eth eth beta / r - (1/2) eth Q / r - eth U
    rhsH = ZZ.aadd(ZZ.ascale(ZZ.ang_rpow(ZZ.ang_dr(ZZ.ang_dr(J)), 1), 1//2),
           ZZ.aadd(ZZ.ang_dr(J),
           ZZ.aadd(ZZ.ang_rpow(eth_spin(eB, 1), -1),
           ZZ.aadd(ZZ.ascale(ZZ.ang_rpow(eth_spin(Q, 1), -1), -1//2),
                   ZZ.ascale(eth_spin(U, 1), -1)))))
    @test ZZ.ang_dr(ZZ.ang_rpow(H, 1)) == rhsH
    # the same H equation on the PLAIN family (beta = 0)
    bt = bondi_transform_teukolsky()
    QP = ZZ.ang_rpow(ZZ.ang_dr(bt.U), 2)
    rhsHP = ZZ.aadd(ZZ.ascale(ZZ.ang_rpow(ZZ.ang_dr(ZZ.ang_dr(bt.J)), 1), 1//2),
            ZZ.aadd(ZZ.ang_dr(bt.J),
            ZZ.aadd(ZZ.ascale(ZZ.ang_rpow(eth_spin(QP, 1), -1), -1//2),
                    ZZ.ascale(eth_spin(bt.U, 1), -1))))
    @test ZZ.ang_dr(ZZ.ang_rpow(ang_du(bt.J), 1)) == rhsHP
end

@testset "interior-point psi0 + anchored end-to-end (stage C part 2)" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    rwt = 41.0
    want(rs) = -9/8*teuk_F(2, 20.0, X, rc, tau)/rs^5
    # (a) symbolic gauge invariance of the J-only psi0 formula: the anchored
    # gauge adds a(u) + b(u)/r tails to J (J_ANCH keys k = 0, 1) which are
    # EXACTLY psi0-silent: psi0_from_J(J_ANCH) == psi0_from_J(J_bondi)
    bt = bondi_transform_teukolsky()
    @test psi0_from_J(J_ANCH) == psi0_from_J(bt.J)
    # (b) exact-slice interior-point psi0 (machine-exact jr): direct
    # spectral D^2 + barycentric rows at production-like radii
    g = ChebRegGrid(rwt, 33)
    jr = [teuk_jr(20.0, rofy_c(g, y), X, rc, tau) for y in g.y]
    for rs in (41.0, 42.025, 60.0, 82.0)
        rel = abs(cheb_psi0_at(g, jr, rs) - want(rs))/abs(want(rs))
        @info "exact-slice psi0_at($rs) rel = $rel"
        @test rel < 1e-6
    end
    # (c) evolved plain-BC interior probe (isolates probe machinery from
    # gauge effects; same configuration as the admitted evolved-psi0 gate)
    g65 = ChebRegGrid(rwt, 65)
    jre = cheb_evolve_J_sat(g65, 16.0, 20.0, 0.00078125, X, rc, tau)
    relp = abs(cheb_psi0_at(g65, jre, 42.025) - want(42.025))/abs(want(42.025))
    @info "evolved plain-BC psi0_at(42.025) rel = $relp (floor class ~0.05)"
    @test relp < 0.2
    # (d) the live (anchored-gauge) BCs are NOT representable in the ringed
    # scheme: jr = r J diverges at scri on the anchored tails (k = 0, 1).
    # Pin the documented failure (iter 44): psi0 off by O(10^3).
    function map_bc_T(::Type{T}, u, rw) where {T}
        A(n) = teuk_F(n, u, T(X), T(rc), T(tau))
        Aser = 3*(A(2)/rw^3 + 3A(1)/rw^4 + 3A(0)/rw^5)
        drA = 3*(-3A(2)/rw^4 - 12A(1)/rw^5 - 15A(0)/rw^6)
        Bser = -(A(3)/rw^2 + 3A(2)/rw^3 + 6A(1)/rw^4 + 6A(0)/rw^5)
        Cser = (A(4)/rw + 2A(3)/rw^2 + 9A(2)/rw^3 + 21A(1)/rw^4 + 21A(0)/rw^5)/4
        dtB = -(A(4)/rw^2 + 3A(3)/rw^3 + 6A(2)/rw^4 + 6A(1)/rw^5)
        dtC = (A(5)/rw + 2A(4)/rw^2 + 9A(3)/rw^3 + 21A(2)/rw^4 + 21A(1)/rw^5)/4
        dtA = 3*(A(3)/rw^3 + 3A(2)/rw^4 + 3A(1)/rw^5)
        worldtube_map(rw; hTT=T(3)/2*(2Cser - Aser), dt_hTT=T(3)/2*(2dtC - dtA),
            hrr=Aser, trace=-Aser, dr_trace=-drA, dth_trace_sc=6Aser,
            dth_dr_trace_sc=6drA, hrth=-3*Bser*rw, dt_hrth=-3*dtB*rw,
            ethb_hrth_P=-3*Bser*rw)
    end
    ringed_bc(u) = (m = map_bc_T(Float64, u, rwt); (m.jr, m.qr, m.ur, m.wr, m.hr))
    jra = cheb_evolve_J_sat_bc(g65, 8.0, 20.0, 0.00078125,
                               zeros(Float64, 65), ringed_bc)
    relring = abs(cheb_psi0_at(g65, jra, 42.025) - want(42.025))/abs(want(42.025))
    @info "ringed solver on anchored BCs: probe rel = $relring (DOCUMENTED FAILURE)"
    @test relring > 10.0    # the ringed scheme must NOT be used with live BCs
    # (e) UNRINGED beta-corrected solver (bondisolver_unringed.jl):
    # machine-exact sweeps on BOTH gauges + the anchored end-to-end gates
    g33 = ChebRegGrid(rwt, 33)
    uops = UnringedOps(g33)
    rof(y) = y == 1.0 ? Inf : 2rwt/(1 - y)
    QA = ZccmJl.ang_rpow(ang_dr(U_ANCH), 2)
    HA = ang_du(J_ANCH)
    sv(ser, kind, r) = series_scalar(ser, kind, 20.0, r, X, rc, tau)
    JtA = [sv(J_ANCH, :s2, rof(y)) for y in g33.y]
    btA = sv(BETA_ANCH, :P, rwt)
    bcA = (sv(QA, :sc, rwt), sv(U_ANCH, :sc, rwt), sv(W_ANCH, :P, rwt),
           sv(HA, :s2, rwt))
    QtA, UtA, WtA, HtA = ut_sweep(uops, JtA, bcA, btA)
    @test maximum(abs.(QtA .- [sv(QA, :sc, rof(y)) for y in g33.y])) < 1e-18
    @test maximum(abs.(UtA .- [sv(U_ANCH, :sc, rof(y)) for y in g33.y])) < 1e-18
    @test maximum(abs.(WtA .- [sv(W_ANCH, :P, rof(y)) for y in g33.y])) < 1e-18
    @test maximum(abs.(HtA .- [sv(HA, :s2, rof(y)) for y in g33.y])) < 1e-18
    for rs in (42.025, 60.0, 82.0)
        @test abs(ut_psi0_at(g33, JtA, rs) - want(rs))/abs(want(rs)) < 1e-6
    end
    # anchored end-to-end (production path) Float64 + du refinement:
    # measured (iter 44): probe 8.8e-5 at du = 7.8e-4, 7.1e-6 at halved du
    unr_bc(u) = (m = map_bc_T(Float64, u, rwt);
                 (m.jr/rwt, m.qr/rwt, m.ur/rwt^2, m.wr/rwt^2, m.hr/rwt, m.beta))
    Jta = ut_evolve_sat(g65, 8.0, 20.0, 0.00078125, zeros(Float64, 65), unr_bc)
    relwa = abs(ut_psi0_at(g65, Jta, rwt) - want(rwt))/abs(want(rwt))
    relpa = abs(ut_psi0_at(g65, Jta, 42.025) - want(42.025))/abs(want(42.025))
    @info "UNRINGED anchored end-to-end psi0: tube rel = $relwa, probe rel = $relpa"
    @test relwa < 0.2
    @test relpa < 5e-4
    Jta2 = ut_evolve_sat(g65, 8.0, 20.0, 0.000390625, zeros(Float64, 65), unr_bc)
    relpa2 = abs(ut_psi0_at(g65, Jta2, 42.025) - want(42.025))/abs(want(42.025))
    @info "UNRINGED anchored probe at du/2: rel = $relpa2 (truncation-dominated)"
    @test relpa2 < 1.5e-5
    # (f) BigFloat-128 anchored end-to-end (no precision pathology;
    # n = 49 spectral truncation limits the probe here, measured 1.3e-3)
    setprecision(BigFloat, 128) do
        rwb = big"41.0"
        gb = ChebRegGrid(rwb, 49)
        unr_bcb(u) = (m = map_bc_T(BigFloat, u, rwb);
                      (m.jr/rwb, m.qr/rwb, m.ur/rwb^2, m.wr/rwb^2,
                       m.hr/rwb, m.beta))
        Jtb = ut_evolve_sat(gb, big"10.0", big"20.0", big"0.00078125",
                            zeros(BigFloat, 49), unr_bcb)
        wantb(rs) = -big(9)/8*teuk_F(2, big"20.0", big"1.0e-5", big"20.0",
                                     big"2.0")/rs^5
        relwb = Float64(abs(ut_psi0_at(gb, Jtb, rwb) - wantb(rwb))/abs(wantb(rwb)))
        relpb = Float64(abs(ut_psi0_at(gb, Jtb, big"42.025") -
                    wantb(big"42.025"))/abs(wantb(big"42.025")))
        @info "BigFloat-128 unringed anchored psi0: tube rel = $relwb, probe rel = $relpb"
        @test relwb < 2e-2
        @test relpb < 1e-2
    end
end

@testset "characteristic pulse injection (2308.10361 Sec V.C, N14 fidelity)" begin
    # ingoing l=2 m=0 J-pulse on the initial null slice; quiescent tube
    # BCs; psi0 probe at the Cauchy boundary radius. Production geometry
    # (mode 6): rwt = 40, r_B = 41, Z = 1e-3, paper pulse shape.
    rwt, rB, Z = 40.0, 41.0, 1.0e-3
    g = ChebRegGrid(rwt, 65)
    Jt0 = ut_pulse_id(g, Z)
    # normalization: Jt = (1/4) sqrt(15/2pi) * Jcal; Jcal(yc) = Z exactly,
    # and y = 0 is a CGL node for even N
    @test any(y -> y == 0.0, g.y)
    @test isapprox(maximum(Jt0), sqrt(15/(2pi))/4*Z; rtol=1e-12)
    @test Jt0[1] == 0.0 && Jt0[end] == 0.0          # compact support
    qbc(u) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # short-window evolution (u: -40 -> 20, i.e. t at r_B up to 61 —
    # includes the first trough) + radial-refinement ladder
    series = Dict{Int,Vector{Float64}}()
    local ts
    for n in (49, 65, 97)
        gn = ChebRegGrid(rwt, n)
        (_, us, vals) = ut_evolve_sat_probe(gn, -rwt, 20.0, 0.00078125,
                                            ut_pulse_id(gn, Z), qbc, rB)
        series[n] = vals
        n == 65 && (ts = us .+ rB)
    end
    pk = maximum(abs.(series[65]))
    e49 = maximum(abs.(series[49] .- series[65]))/pk
    e97 = maximum(abs.(series[97] .- series[65]))/pk
    @info "pulse psi0 ladder: n49-vs-n65 = $e49, n97-vs-n65 = $e97 of window peak"
    @test e49 < 2e-3        # measured 3.9e-4 over the full run
    @test e97 < 1e-5        # measured 2.1e-7 of peak (n=65 converged)
    # conditioning: psi0 is O(1e-3) of the J data (vs 1e-8 for the purely
    # outgoing Teukolsky test — the fidelity-test design point)
    @test pk/maximum(abs.(Jt0)) > 1e-3   # measured 7.0e-3 full-run
    # regression pins from the admitted full-run reference
    # (results/numerical/n14_pulse_ref.csv): first trough at t = 52.50,
    # psi0(trough) = -1.718782e-6
    itr = argmin(series[65])
    @test abs(ts[itr] - 52.50) < 0.5
    @test isapprox(series[65][itr], -1.718782e-6; rtol=1e-3)
end

@testset "worldtube map (stage B): identity + anchored evolution + psi0" begin
    X, rc, tau = 1.0e-5, 20.0, 2.0
    rwt = 41.0
    # identity gate: map outputs == anchored-gauge series values at wt
    for u in (18.0, 20.0, 23.0)
        # analytic 'live' Cauchy data scalars at the worldtube
        A(n) = teuk_F(n, u, X, rc, tau)
        Aser = 3*(A(2)/rwt^3 + 3A(1)/rwt^4 + 3A(0)/rwt^5)
        drA = 3*(-3A(2)/rwt^4 - 12A(1)/rwt^5 - 15A(0)/rwt^6)
        Bser = -(A(3)/rwt^2 + 3A(2)/rwt^3 + 6A(1)/rwt^4 + 6A(0)/rwt^5)
        Cser = (A(4)/rwt + 2A(3)/rwt^2 + 9A(2)/rwt^3 + 21A(1)/rwt^4 + 21A(0)/rwt^5)/4
        dtB = -(A(4)/rwt^2 + 3A(3)/rwt^3 + 6A(2)/rwt^4 + 6A(1)/rwt^5)
        hTT = 3*(2Cser - Aser)/4 * 2     # (3/2)(2C-A): s^2-scalar
        hTT = (3/2)*(2Cser - Aser)
        dtC = (A(5)/rwt + 2A(4)/rwt^2 + 9A(3)/rwt^3 + 21A(2)/rwt^4 + 21A(1)/rwt^5)/4
        dtA = 3*(A(3)/rwt^3 + 3A(2)/rwt^4 + 3A(1)/rwt^5)
        dt_hTT = (3/2)*(2dtC - dtA)
        hrr = Aser                        # P-scalar of A(3c^2-1)
        trace = -Aser
        dr_trace = -drA
        dth_trace_sc = 6Aser              # d_theta(-A P) -> sc: +6A
        dth_dr_trace_sc = 6drA
        hrth = -3*Bser*rwt / (-3)         # h_rth = -3 B r sc: sc-scalar = -3 B r
        hrth = -3*Bser*rwt
        dt_hrth = -3*dtB*rwt
        ethb_hrth_P = -3*Bser*rwt         # ethbar(sc X) = P X: scalar = -3 B r
        m = worldtube_map(rwt; hTT=hTT, dt_hTT=dt_hTT, hrr=hrr, trace=trace,
            dr_trace=dr_trace, dth_trace_sc=dth_trace_sc,
            dth_dr_trace_sc=dth_dr_trace_sc, hrth=hrth, dt_hrth=dt_hrth,
            ethb_hrth_P=ethb_hrth_P)
        # gate vs the anchored series at wt
        # the identity is exact; Float64 evaluation differs only by
        # cancellation roundoff (W ~ 1e-11 from 1e-9-scale terms) — gate
        # at the cancellation-conditioned tolerance and EXACTLY in BigFloat
        for (got, ser, kind, w) in ((m.jr, J_ANCH, :s2, rwt),
                                    (m.ur, U_ANCH, :sc, rwt^2),
                                    (m.wr, W_ANCH, :P, rwt^2),
                                    (m.beta, BETA_ANCH, :P, 1.0))
            want = w*series_scalar(ser, kind, u, rwt, X, rc, tau)
            @test isapprox(got, want; rtol=1e-9, atol=1e-30)
        end
    end
end

@testset "worldtube map identity in BigFloat (exact)" begin
    setprecision(BigFloat, 192) do
        X, rc, tau = big"1.0e-5", big"20.0", big"2.0"
        rwt = big"41.0"
        u = big"20.0"
        A(n) = teuk_F(n, u, X, rc, tau)
        Aser = 3*(A(2)/rwt^3 + 3A(1)/rwt^4 + 3A(0)/rwt^5)
        drA = 3*(-3A(2)/rwt^4 - 12A(1)/rwt^5 - 15A(0)/rwt^6)
        Bser = -(A(3)/rwt^2 + 3A(2)/rwt^3 + 6A(1)/rwt^4 + 6A(0)/rwt^5)
        Cser = (A(4)/rwt + 2A(3)/rwt^2 + 9A(2)/rwt^3 + 21A(1)/rwt^4 + 21A(0)/rwt^5)/4
        dtB = -(A(4)/rwt^2 + 3A(3)/rwt^3 + 6A(2)/rwt^4 + 6A(1)/rwt^5)
        dtC = (A(5)/rwt + 2A(4)/rwt^2 + 9A(3)/rwt^3 + 21A(2)/rwt^4 + 21A(1)/rwt^5)/4
        dtA = 3*(A(3)/rwt^3 + 3A(2)/rwt^4 + 3A(1)/rwt^5)
        m = worldtube_map(rwt; hTT=big(3)/2*(2Cser-Aser), dt_hTT=big(3)/2*(2dtC-dtA),
            hrr=Aser, trace=-Aser, dr_trace=-drA, dth_trace_sc=6Aser,
            dth_dr_trace_sc=6drA, hrth=-3*Bser*rwt, dt_hrth=-3*dtB*rwt,
            ethb_hrth_P=-3*Bser*rwt)
        for (got, ser, kind, w) in ((m.jr, J_ANCH, :s2, rwt),
                                    (m.ur, U_ANCH, :sc, rwt^2),
                                    (m.wr, W_ANCH, :P, rwt^2),
                                    (m.beta, BETA_ANCH, :P, big"1.0"))
            want = w*series_scalar(ser, kind, u, rwt, X, rc, tau)
            @test Float64(abs(got - want)/max(abs(want), big"1e-60")) < 1e-35
        end
    end
end
