import numpy as np
# Reimplement the VERIFIED TeukolskyMetric (z4c_ccm_teukolsky.cpp) in Python:
# F-ladder (physicists' Hermite), Q/R combos, A/B/C radial, l=2 m=0 angular,
# frame eth = eph x rh. Then test: (1) h trace-free, (2) high-order FD of h
# accurate, (3) det(delta+h) ~ 1 at X=2 (conformal-metric sanity).
def Fn(n,u,X,rc,tau):
    s=(u-rc)/tau; s2=s*s
    H=[1.0,2*s,4*s2-2,s*(8*s2-12),(16*s2-48)*s2+12,s*((32*s2-160)*s2+120),
       ((64*s2-480)*s2+720)*s2-120][n]
    pref=(-1.0/tau)**n
    return X*pref*H*np.exp(-s2)
def TeukQ(n,t,r,X,rc,tau):
    sgn=1.0 if n%2==0 else -1.0
    return Fn(n,t-r,X,rc,tau)-sgn*Fn(n,t+r,X,rc,tau)
def hmetric(t,x,y,z,X,rc,tau):
    r=max((x*x+y*y+z*z)**0.5,1e-12); ir=1/r
    Q=[TeukQ(n,t,r,X,rc,tau) for n in range(5)]
    ir2,ir3,ir4,ir5=ir*ir,ir**3,ir**4,ir**5
    Af=3*(Q[2]*ir3+3*Q[1]*ir4+3*Q[0]*ir5)
    Bf=-(Q[3]*ir2+3*Q[2]*ir3+6*Q[1]*ir4+6*Q[0]*ir5)
    Cf=0.25*(Q[4]*ir+2*Q[3]*ir2+9*Q[2]*ir3+21*Q[1]*ir4+21*Q[0]*ir5)
    cth=z*ir; s2=max(1-cth*cth,0.0); sth=s2**0.5
    frr=0.5*(1+3*(1-2*s2)); frt=-3*sth*cth
    f1tt=3*s2; f2tt=-1.0; f1pp=-3*s2; f2pp=3*s2-1
    hrr=Af*frr; hrt=Bf*frt; htt=Cf*f1tt+Af*f2tt; hpp=Cf*f1pp+Af*f2pp
    rh=np.array([x*ir,y*ir,z*ir])
    rho=max((x*x+y*y)**0.5,1e-300)
    eph=np.array([-y/rho,x/rho,0.0]) if rho>1e-12*r else np.array([0.0,1.0,0.0])
    eth=np.cross(eph,rh)
    h=np.zeros((3,3))
    for a in range(3):
        for b in range(3):
            h[a,b]=hrr*rh[a]*rh[b]+hrt*(rh[a]*eth[b]+eth[a]*rh[b])+htt*eth[a]*eth[b]+hpp*eph[a]*eph[b]
    return h
X,rc,tau=2.0,20.0,2.0
pts=[(12.,5.,8.),(20.,0.,3.),(-7.,15.,22.),(30.,-10.,5.)]
print("X=2 Teukolsky free-data checks (t=0):")
for (x,y,z) in pts:
    h=hmetric(0,x,y,z,X,rc,tau)
    tr=np.trace(h)
    g=np.eye(3)+h; det=np.linalg.det(g)
    # 4th-order FD of h vs analytic-ish (compare 4th vs 6th order to gauge FD floor)
    eps=1e-3; dh4=np.zeros((3,3,3))
    for k,dx in enumerate([(eps,0,0),(0,eps,0),(0,0,eps)]):
        hp1=hmetric(0,x+dx[0],y+dx[1],z+dx[2],X,rc,tau)
        hm1=hmetric(0,x-dx[0],y-dx[1],z-dx[2],X,rc,tau)
        hp2=hmetric(0,x+2*dx[0],y+2*dx[1],z+2*dx[2],X,rc,tau)
        hm2=hmetric(0,x-2*dx[0],y-2*dx[1],z-2*dx[2],X,rc,tau)
        dh4[k]=(-hp2+8*hp1-8*hm1+hm2)/(12*eps)
    # compare to 2nd-order to estimate convergence
    dh2=np.zeros((3,3,3))
    for k,dx in enumerate([(eps,0,0),(0,eps,0),(0,0,eps)]):
        hp1=hmetric(0,x+dx[0],y+dx[1],z+dx[2],X,rc,tau)
        hm1=hmetric(0,x-dx[0],y-dx[1],z-dx[2],X,rc,tau)
        dh2[k]=(hp1-hm1)/(2*eps)
    fd_spread=np.max(np.abs(dh4-dh2))
    print(f"  ({x:5.0f},{y:4.0f},{z:4.0f}): |h|max={np.max(np.abs(h)):.3e} tr(h)={tr:+.2e} det(d+h)-1={det-1:+.3e} |dh4-dh2|={fd_spread:.1e}")
