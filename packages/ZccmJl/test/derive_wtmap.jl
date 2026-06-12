using ZccmJl, LinearAlgebra
const Z = ZccmJl
RWT = 41

"evaluate an AngSeries at r = RWT: dict (n,p,q) => exact rational"
function at_wt(a)
    out = Dict{NTuple{3,Int},Rational{BigInt}}()
    for ((n, k, p, q), v) in a
        key = (n, p, q)
        out[key] = get(out, key, big(0)//1) + v//big(RWT)^k
    end
    Dict(k => v for (k, v) in out if v != 0)
end

"subtract the r_wt value as k=0 entries (worldtube anchor)"
function anchor(a)
    out = copy(a)
    for ((n, p, q), v) in at_wt(a)
        key = (n, 0, p, q)
        nv = get(out, key, big(0)//1) - v
        nv == 0 ? delete!(out, key) : (out[key] = nv)
    end
    out
end

bt = bondi_transform_teukolsky()
zu_a = anchor(bt.zu)
zth_a = anchor(bt.zth)

# rebuild h pieces (mirror bonditransform internals via its tables)
A_ = ((2,3,3),(1,4,9),(0,5,9)); BR_ = ((3,1,-1),(2,2,-3),(1,3,-6),(0,4,-6))
C_ = ((4,1,1//4),(3,2,1//2),(2,3,9//4),(1,4,21//4),(0,5,21//4))
sm(tab,p,q,f) = Z.shape_mul(tab,p,q,f)
h_rr = Z.aadd(sm(A_,2,0,3), sm(A_,0,0,-1))
h_rth = sm(BR_,1,1,-3)
h_T = Z.aadd(Z.aadd(sm(C_,0,0,3), sm(C_,2,0,-3)), sm(A_,0,0,-1))
h_P = Z.aadd(Z.aadd(sm(C_,0,0,-3), sm(C_,2,0,3)),
             Z.aadd(sm(A_,0,0,2), sm(A_,2,0,-3)))
trace = Z.aadd(h_T, h_P)

# anchored zr (algebraic) and read-offs
zr_a = Z.ang_rpow(Z.ascale(Z.aadd(trace,
        Z.aadd(Z.ascale(Z.ang_dtheta(zth_a), 2),
               Z.ascale(Z.ang_cot_mul(zth_a), 2))), -1//4), 1)
J_a = Z.aadd(Z.ascale(Z.aadd(h_T, Z.ascale(h_P, -1)), 1//2),
             Z.aadd(Z.ang_dtheta(zth_a), Z.ascale(Z.ang_cot_mul(zth_a), -1)))
beta_a = Z.ascale(Z.aadd(Z.ang_du(zu_a), Z.aadd(Z.ang_dr(zu_a), Z.ang_dr(zr_a))), 1//2)
U_a = Z.ascale(Z.aadd(Z.ang_du(zth_a),
        Z.ang_rpow(Z.ascale(Z.aadd(Z.ang_dtheta(zu_a), Z.ang_dtheta(zr_a)), -1), -2)), -1)
dg_uu = Z.ascale(Z.aadd(Z.ang_du(zu_a), Z.ang_du(zr_a)), -2)
W_a = Z.ang_rpow(Z.ascale(Z.aadd(dg_uu, Z.ascale(beta_a, 2)), -1), -1)

# psi0 invariance gate on the ANCHORED J (gauge invariance end-to-end)
p0a = psi0_from_J(J_a)
p0o = psi0_from_J(bt.J)
println("G0 psi0 identical in anchored gauge: ", p0a == p0o)

# scalarize per shape: J,s^2 pairs (p,q)=(0,0)/(2,0); U: (1,1); W,beta: P-shape
function scal(d, kind)
    out = Dict{Int,Rational{BigInt}}()  # n => coeff
    for ((n, p, q), v) in d
        if kind == :s2 && q == 0 && p == 0
            out[n] = get(out, n, big(0)//1) + v
        elseif kind == :sc && q == 1 && p == 1
            out[n] = get(out, n, big(0)//1) + v
        elseif kind == :P && q == 0 && p == 2
            out[n] = get(out, n, big(0)//1) + v//3   # P = 3c^2 - 1
        end
    end
    Dict(k => v for (k, v) in out if v != 0)
end

# local data candidates at wt (scalarized)
cand = Dict(
 :hTT => scal(at_wt(Z.ascale(Z.aadd(h_T, Z.ascale(h_P, -1)), 1//2)), :s2),
 :hrr => scal(at_wt(h_rr), :P),
 :trace => scal(at_wt(trace), :P),
 :hrth => scal(at_wt(h_rth), :sc),
 :dr_hTT => scal(at_wt(Z.ang_dr(Z.ascale(Z.aadd(h_T, Z.ascale(h_P, -1)), 1//2))), :s2),
 :dr_trace => scal(at_wt(Z.ang_dr(trace)), :P),
 :dr_hrth => scal(at_wt(Z.ang_dr(h_rth)), :sc),
 :dt_hTT => scal(at_wt(Z.ang_du(Z.ascale(Z.aadd(h_T, Z.ascale(h_P, -1)), 1//2))), :s2),
 :dt_hrth => scal(at_wt(Z.ang_du(h_rth)), :sc),
 :dt_trace => scal(at_wt(Z.ang_du(trace)), :P),
)
targets = Dict(:J => scal(at_wt(J_a), :s2), :H => scal(at_wt(Z.ang_du(J_a)), :s2),
               :U => scal(at_wt(U_a), :sc), :beta => scal(at_wt(beta_a), :P),
               :W => scal(at_wt(W_a), :P))

"solve target = sum c_i cand_i over the F^n channels (exact); allow r_wt powers"
function fitmap(tname, names)
    chans = Set{Int}()
    for nm in names; union!(chans, keys(cand[nm])); end
    union!(chans, keys(targets[tname]))
    chans = sort(collect(chans))
    A = Rational{BigInt}[get(cand[nm], n, big(0)//1) for n in chans, nm in names]
    b = Rational{BigInt}[get(targets[tname], n, big(0)//1) for n in chans]
    x, ok = Z.exact_solve(A, b)
    println("$tname: ", ok ? "LOCAL MAP " * join(["$(names[i])*($(x[i]))" for i in eachindex(names) if x[i] != 0], " + ") : "NOT representable with $(names)")
    ok
end
println("-- map fits (coefficients exact rationals; r_wt = 41 absorbed) --")
fitmap(:J, [:hTT])
fitmap(:H, [:dt_hTT])
fitmap(:U, [:hrth, :dr_hrth, :dt_hrth])
fitmap(:beta, [:hrr, :trace, :dr_trace, :hrth])
fitmap(:W, [:hrr, :trace, :dr_trace, :dt_trace, :hrth, :dt_hrth])
