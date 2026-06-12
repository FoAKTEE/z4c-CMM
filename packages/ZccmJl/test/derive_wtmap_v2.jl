using ZccmJl
const Z = ZccmJl
RWT = 41
function at_wt(a)
    out = Dict{NTuple{3,Int},Rational{BigInt}}()
    for ((n,k,p,q),v) in a
        key=(n,p,q); out[key]=get(out,key,big(0)//1)+v//big(RWT)^k
    end
    Dict(k=>v for (k,v) in out if v!=0)
end
function anchor(a)
    out=copy(a)
    for ((n,p,q),v) in at_wt(a)
        key=(n,0,p,q); nv=get(out,key,big(0)//1)-v
        nv==0 ? delete!(out,key) : (out[key]=nv)
    end
    out
end
bt = bondi_transform_teukolsky()
zu_a = anchor(bt.zu); zth_a = anchor(bt.zth)
A_=((2,3,3),(1,4,9),(0,5,9)); BR_=((3,1,-1),(2,2,-3),(1,3,-6),(0,4,-6))
C_=((4,1,1//4),(3,2,1//2),(2,3,9//4),(1,4,21//4),(0,5,21//4))
sm(t,p,q,f)=Z.shape_mul(t,p,q,f)
h_rr=Z.aadd(sm(A_,2,0,3),sm(A_,0,0,-1))
h_rth=sm(BR_,1,1,-3)
h_T=Z.aadd(Z.aadd(sm(C_,0,0,3),sm(C_,2,0,-3)),sm(A_,0,0,-1))
h_P=Z.aadd(Z.aadd(sm(C_,0,0,-3),sm(C_,2,0,3)),Z.aadd(sm(A_,0,0,2),sm(A_,2,0,-3)))
trace=Z.aadd(h_T,h_P)
S = Z.aadd(trace, Z.aadd(Z.ascale(Z.ang_dtheta(zth_a),2),
                         Z.ascale(Z.ang_cot_mul(zth_a),2)))
zr_a = Z.ang_rpow(Z.ascale(S,-1//4),1)
beta_a = Z.ascale(Z.aadd(Z.ang_du(zu_a),
                  Z.aadd(Z.ang_dr(zu_a),Z.ang_dr(zr_a))),1//2)
# hand formula as a SERIES (valid only AT wt, so compare at_wt):
ethb_hrth = Z.aadd(Z.ang_dtheta(h_rth), Z.ang_cot_mul(h_rth))
hand = Z.aadd(Z.aadd(Z.ascale(h_rr,1//4), Z.ascale(trace,-1//8)),
       Z.aadd(Z.ascale(Z.ang_rpow(Z.ang_dr(trace),1),-1//8),
              Z.ascale(Z.ang_rpow(ethb_hrth,-1),1//4)))
d_beta = at_wt(beta_a); d_hand = at_wt(hand)
ks = union(keys(d_beta), keys(d_hand))
println("beta_a vs hand formula at wt (channel: script | hand):")
for k in sort(collect(ks))
    println("  $k: ", get(d_beta,k,big(0)//1), " | ", get(d_hand,k,big(0)//1))
end
# purity check of beta_a at wt: P-shape means (n,0,0) = -(n,2,0)/3
println("purity of beta_a at wt:")
for ((n,p,q),v) in sort(collect(d_beta))
    p==2 && println("  n=$n: p0=", get(d_beta,(n,0,0),big(0)//1), " vs -(p2)/3=", -v//3)
end

# CORRECTED anchored gauge: the zu anchor feeds the xi ODE
# xi_a = xi_old + dtheta(zu_wt)/r - (const so xi_a(wt) = 0)
zuwt = Dict{NTuple{4,Int},Rational{BigInt}}()   # k=0 series of zu at wt
for ((n,p,q),v) in at_wt(bt.zu)
    zuwt[(n,0,p,q)] = v
end
corr = Z.ang_rpow(Z.ang_dtheta(zuwt), -1)       # dtheta(zuwt)/r
zth_c = anchor(Z.aadd(bt.zth, corr))            # consistent + anchored
# verify the g_rtheta condition for the corrected gauge:
res = Z.aadd(h_rth, Z.aadd(Z.ang_rpow(Z.ang_dr(zth_c), 2),
                           Z.ascale(Z.ang_dtheta(anchor(bt.zu)), -1)))
println("corrected gauge g_rtheta residual: ", isempty(Dict(k=>v for (k,v) in res if v!=0)) ? "ZERO" : sort(collect(res)))
# rebuild S, zr, beta, U with the corrected xi
zu_c = anchor(bt.zu)
S_c = Z.aadd(trace, Z.aadd(Z.ascale(Z.ang_dtheta(zth_c),2),
                           Z.ascale(Z.ang_cot_mul(zth_c),2)))
zr_c = Z.ang_rpow(Z.ascale(S_c,-1//4),1)
beta_c = Z.ascale(Z.aadd(Z.ang_du(zu_c),
                  Z.aadd(Z.ang_dr(zu_c),Z.ang_dr(zr_c))),1//2)
U_c = Z.ascale(Z.aadd(Z.ang_du(zth_c),
      Z.ang_rpow(Z.ascale(Z.aadd(Z.ang_dtheta(zu_c),Z.ang_dtheta(zr_c)),-1),-2)),-1)
J_c = Z.aadd(Z.ascale(Z.aadd(h_T,Z.ascale(h_P,-1)),1//2),
             Z.aadd(Z.ang_dtheta(zth_c),Z.ascale(Z.ang_cot_mul(zth_c),-1)))
# gates: psi0 invariance on the CONSISTENT anchored J; beta vs hand formula
println("psi0 invariant (consistent anchored J): ",
        psi0_from_J(J_c) == psi0_from_J(bt.J))
d_beta_c = at_wt(beta_c)
println("beta_c at wt vs hand:")
for k in sort(collect(union(keys(d_beta_c), keys(d_hand))))
    println("  $k: ", get(d_beta_c,k,big(0)//1), " | ", get(d_hand,k,big(0)//1))
end
println("U_c at wt: ", sort(collect(at_wt(U_c))))
println("J_c at wt == naive hTT at wt: ",
        at_wt(J_c) == at_wt(Z.ascale(Z.aadd(h_T,Z.ascale(h_P,-1)),1//2)))

# refit U_c and W_c against local candidates (corrected gauge)
function scal2(d, kind)
    out = Dict{Int,Rational{BigInt}}()
    for ((n,p,q),v) in d
        if kind == :sc && q==1 && p==1
            out[n]=get(out,n,big(0)//1)+v
        elseif kind == :P && q==0 && p==2
            out[n]=get(out,n,big(0)//1)+v//3
        end
    end
    Dict(k=>v for (k,v) in out if v!=0)
end
candU = Dict(
 :dth_trace => scal2(at_wt(Z.ang_dtheta(trace)), :sc),
 :hrth => scal2(at_wt(h_rth), :sc),
 :dt_hrth => scal2(at_wt(Z.ang_du(h_rth)), :sc),
 :dr_hrth => scal2(at_wt(Z.ang_dr(h_rth)), :sc),
 :dt_dth_hrr => scal2(at_wt(Z.ang_du(Z.ang_dtheta(h_rr))), :sc),
)
tU = scal2(at_wt(U_c), :sc)
function fit2(target, cands, names)
    chans = Set{Int}(keys(target))
    for nm in names; union!(chans, keys(cands[nm])); end
    chans = sort(collect(chans))
    A = Rational{BigInt}[get(cands[nm],n,big(0)//1) for n in chans, nm in names]
    b = Rational{BigInt}[get(target,n,big(0)//1) for n in chans]
    x, ok = Z.exact_solve(A, b)
    ok ? join(["$(names[i])*($(x[i]))" for i in eachindex(names) if x[i]!=0], " + ") : "NOT representable"
end
println("U_c fit: ", fit2(tU, candU, [:dth_trace, :hrth, :dt_hrth, :dr_hrth, :dt_dth_hrr]))
# W with the corrected gauge
dg_uu_c = Z.ascale(Z.aadd(Z.ang_du(zu_c), Z.ang_du(zr_c)), -2)
W_c = Z.ang_rpow(Z.ascale(Z.aadd(dg_uu_c, Z.ascale(beta_c,2)),-1),-1)
tW = scal2(at_wt(W_c), :P)
candW = Dict(
 :hrr => scal2(at_wt(h_rr), :P),
 :dr_trace => scal2(at_wt(Z.ang_dr(trace)), :P),
 :dt_trace => scal2(at_wt(Z.ang_du(trace)), :P),
 :ethb_hrth => scal2(at_wt(Z.aadd(Z.ang_dtheta(h_rth),Z.ang_cot_mul(h_rth))), :P),
 :dt_ethb_hrth => scal2(at_wt(Z.ang_du(Z.aadd(Z.ang_dtheta(h_rth),Z.ang_cot_mul(h_rth)))), :P),
 :dr_dt_trace => scal2(at_wt(Z.ang_dr(Z.ang_du(trace))), :P),
)
println("W_c fit: ", fit2(tW, candW, [:hrr, :dr_trace, :dt_trace, :ethb_hrth, :dt_ethb_hrth, :dr_dt_trace]))
