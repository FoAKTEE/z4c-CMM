# smoke_unringed.jl — iter 44 verification ladder for the unringed
# beta-corrected solver (before runtests integration).
using ZccmJl
const Z = ZccmJl

X, rc, tau = 1.0e-5, 20.0, 2.0
rwt = 41.0
u0 = 20.0
want(rs) = -9/8*teuk_F(2, u0, X, rc, tau)/rs^5

# ---- V1: exact-slice sweep, PLAIN gauge (beta = 0) ----
g = ChebRegGrid(rwt, 33)
ops = UnringedOps(g)
rof(y) = y == 1.0 ? Inf : 2rwt/(1 - y)
Jt = [teuk_jr(u0, rof(y), X, rc, tau)/rof(y) for y in g.y]
Jt[end] = 0.0   # Jt = jr/r -> 0 at scri (plain gauge)
bc = (teuk_qr(u0, rwt, X, rc, tau)/rwt, teuk_ur(u0, rwt, X, rc, tau)/rwt^2,
      teuk_wr(u0, rwt, X, rc, tau)/rwt^2, teuk_hr(u0, rwt, X, rc, tau)/rwt)
Qt, Ut, Wt, Ht = ut_sweep(ops, Jt, bc, 0.0)
eQ = maximum(abs.(Qt .- [y == 1.0 ? 0.0 : teuk_qr(u0, rof(y), X, rc, tau)/rof(y) for y in g.y]))
eU = maximum(abs.(Ut .- [y == 1.0 ? 0.0 : teuk_ur(u0, rof(y), X, rc, tau)/rof(y)^2 for y in g.y]))
eW = maximum(abs.(Wt .- [y == 1.0 ? 0.0 : teuk_wr(u0, rof(y), X, rc, tau)/rof(y)^2 for y in g.y]))
eH = maximum(abs.(Ht .- [y == 1.0 ? 0.0 : teuk_hr(u0, rof(y), X, rc, tau)/rof(y) for y in g.y]))
println("V1 plain exact-slice sweep errors: Q=$eQ U=$eU W=$eW H=$eH")

# ---- V2: exact-slice sweep, ANCHORED gauge (beta != 0) ----
QA = Z.ang_rpow(Z.ang_dr(U_ANCH), 2)
HA = Z.ang_du(J_ANCH)
sv(ser, kind, r) = series_scalar(ser, kind, u0, r, X, rc, tau)
JtA = [sv(J_ANCH, :s2, rof(y)) for y in g.y]
btA = sv(BETA_ANCH, :P, rwt)
bcA = (sv(QA, :sc, rwt), sv(U_ANCH, :sc, rwt), sv(W_ANCH, :P, rwt),
       sv(HA, :s2, rwt))
QtA, UtA, WtA, HtA = ut_sweep(ops, JtA, bcA, btA)
eQA = maximum(abs.(QtA .- [sv(QA, :sc, rof(y)) for y in g.y]))
eUA = maximum(abs.(UtA .- [sv(U_ANCH, :sc, rof(y)) for y in g.y]))
eWA = maximum(abs.(WtA .- [sv(W_ANCH, :P, rof(y)) for y in g.y]))
eHA = maximum(abs.(HtA .- [sv(HA, :s2, rof(y)) for y in g.y]))
println("V2 anchored exact-slice sweep errors: Q=$eQA U=$eUA W=$eWA H=$eHA")

# ---- V3: psi0 probe on the anchored exact slice ----
for rs in (41.0, 42.025, 60.0, 82.0)
    rel = abs(ut_psi0_at(g, JtA, rs) - want(rs))/abs(want(rs))
    println("V3 anchored exact-slice psi0_at($rs) rel = $rel")
end

# ---- V4: evolution with PLAIN BCs (regression vs ringed numbers) ----
g65 = ChebRegGrid(rwt, 65)
plain_bc(u) = (teuk_jr(u, rwt, X, rc, tau)/rwt,
               teuk_qr(u, rwt, X, rc, tau)/rwt,
               teuk_ur(u, rwt, X, rc, tau)/rwt^2,
               teuk_wr(u, rwt, X, rc, tau)/rwt^2,
               teuk_hr(u, rwt, X, rc, tau)/rwt, 0.0)
Jte = ut_evolve_sat(g65, 16.0, 20.0, 0.00078125,
                    [teuk_jr(16.0, rof(y), X, rc, tau)/max(rof(y), rwt) *
                     (y == 1.0 ? 0.0 : 1.0) for y in g65.y], plain_bc)
relp = abs(ut_psi0_at(g65, Jte, 42.025) - want(42.025))/abs(want(42.025))
relw = abs(ut_psi0_at(g65, Jte, rwt) - want(rwt))/abs(want(rwt))
println("V4 plain-BC evolved psi0: tube rel = $relw  probe(42.025) rel = $relp")

# ---- V5: THE gate — anchored end-to-end via worldtube_map ----
function map_bc(u)
    A(n) = teuk_F(n, u, X, rc, tau)
    Aser = 3*(A(2)/rwt^3 + 3A(1)/rwt^4 + 3A(0)/rwt^5)
    drA = 3*(-3A(2)/rwt^4 - 12A(1)/rwt^5 - 15A(0)/rwt^6)
    Bser = -(A(3)/rwt^2 + 3A(2)/rwt^3 + 6A(1)/rwt^4 + 6A(0)/rwt^5)
    Cser = (A(4)/rwt + 2A(3)/rwt^2 + 9A(2)/rwt^3 + 21A(1)/rwt^4 + 21A(0)/rwt^5)/4
    dtB = -(A(4)/rwt^2 + 3A(3)/rwt^3 + 6A(2)/rwt^4 + 6A(1)/rwt^5)
    dtC = (A(5)/rwt + 2A(4)/rwt^2 + 9A(3)/rwt^3 + 21A(2)/rwt^4 + 21A(1)/rwt^5)/4
    dtA = 3*(A(3)/rwt^3 + 3A(2)/rwt^4 + 3A(1)/rwt^5)
    m = worldtube_map(rwt; hTT=(3/2)*(2Cser - Aser), dt_hTT=(3/2)*(2dtC - dtA),
        hrr=Aser, trace=-Aser, dr_trace=-drA, dth_trace_sc=6Aser,
        dth_dr_trace_sc=6drA, hrth=-3*Bser*rwt, dt_hrth=-3*dtB*rwt,
        ethb_hrth_P=-3*Bser*rwt)
    (m.jr/rwt, m.qr/rwt, m.ur/rwt^2, m.wr/rwt^2, m.hr/rwt, m.beta)
end
Jta = ut_evolve_sat(g65, 8.0, 20.0, 0.00078125, zeros(Float64, 65), map_bc)
relwa = abs(ut_psi0_at(g65, Jta, rwt) - want(rwt))/abs(want(rwt))
relpa = abs(ut_psi0_at(g65, Jta, 42.025) - want(42.025))/abs(want(42.025))
println("V5 ANCHORED end-to-end psi0: tube rel = $relwa  probe(42.025) rel = $relpa")
