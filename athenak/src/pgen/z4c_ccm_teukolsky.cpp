//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file z4c_ccm_teukolsky.cpp
//! \brief Teukolsky-wave initial data (even parity, l=2, m=0) on flat space, for the
//! Z4c-CCM reproduction of the CCM 2024 paper tests (arXiv:2308.10361 Sec. V):
//! metric per its eq:Teukolsky_metric with radial functions eq:Teukolsky_ABC,
//! angular basis eq:Teukolsky_angular_basis, Gaussian profile eq:Gaussian_pulse_F
//! (ledger: reformulate/z4c-CMM/paper_2308.10361/derivation.md).
//!
//! Regularity at the origin uses the standard combination
//!     Q^(n)(t,r) = F^(n)(t-r) - (-1)^n F^(n)(t+r),
//! for which A, B, C are finite as r -> 0 (e.g. A -> -2 F^(5)(t)/5 ... verified
//! numerically by the battery: finite fields + converging constraints).
//! K_ij = -1/2 d_t h_ij (alpha = 1, beta = 0 at this order) uses
//!     d_t Q^(n) = R^(n+1),   R^(n)(t,r) = F^(n)(t-r) + (-1)^n F^(n)(t+r).
//!
//! <problem> parameters: amp (X), r_c (pulse center), tau (width).
//! History output: Minkowski-deviation norms (as z4c_stability).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/adm.hpp"
#include "mesh/mesh.hpp"
#include "z4c/z4c.hpp"
#include "z4c/z4c_ccm.hpp"   // TeukolskyF (shared with the CCM datum provider)
#include "driver/driver.hpp"
#include "pgen/pgen.hpp"

void Z4cCCMTeukolskyErrors(HistoryData *pdata, Mesh *pm);

namespace {

//----------------------------------------------------------------------------------------
//! \fn TeukQ / TeukR
//! \brief regular (Q) and time-derivative (R) retarded/advanced combinations
KOKKOS_INLINE_FUNCTION
static Real TeukQ(const int n, const Real t, const Real r,
                  const Real X, const Real rc, const Real tau) {
  const Real sgn = (n % 2 == 0) ? 1.0 : -1.0;
  return z4c::TeukolskyF(n, t - r, X, rc, tau)
       - sgn * z4c::TeukolskyF(n, t + r, X, rc, tau);
}

KOKKOS_INLINE_FUNCTION
static Real TeukR(const int n, const Real t, const Real r,
                  const Real X, const Real rc, const Real tau) {
  const Real sgn = (n % 2 == 0) ? 1.0 : -1.0;
  return z4c::TeukolskyF(n, t - r, X, rc, tau)
       + sgn * z4c::TeukolskyF(n, t + r, X, rc, tau);
}

//----------------------------------------------------------------------------------------
//! \fn TeukolskyMetric
//! \brief h_ij (Cartesian, orthonormal-frame transform) and dot h_ij at (t, x, y, z)
KOKKOS_INLINE_FUNCTION
static void TeukolskyMetric(const Real t, const Real x, const Real y, const Real z,
                            const Real X, const Real rc, const Real tau,
                            Real h[3][3], Real hdot[3][3]) {
  const Real r = fmax(sqrt(x*x + y*y + z*z), 1e-12);
  const Real ir = 1.0/r;
  // radial functions: A = 3[Q2/r^3 + 3 Q1/r^4 + 3 Q0/r^5], etc.
  // (eq:Teukolsky_ABC with the regular Q combinations)
  const Real Q0 = TeukQ(0,t,r,X,rc,tau), Q1 = TeukQ(1,t,r,X,rc,tau);
  const Real Q2 = TeukQ(2,t,r,X,rc,tau), Q3 = TeukQ(3,t,r,X,rc,tau);
  const Real Q4 = TeukQ(4,t,r,X,rc,tau);
  const Real R1 = TeukR(1,t,r,X,rc,tau), R2 = TeukR(2,t,r,X,rc,tau);
  const Real R3 = TeukR(3,t,r,X,rc,tau), R4 = TeukR(4,t,r,X,rc,tau);
  const Real R5 = TeukR(5,t,r,X,rc,tau);
  const Real ir2 = ir*ir, ir3 = ir2*ir, ir4 = ir3*ir, ir5 = ir4*ir;

  const Real Af = 3.0*(Q2*ir3 + 3.0*Q1*ir4 + 3.0*Q0*ir5);
  const Real Bf = -(Q3*ir2 + 3.0*Q2*ir3 + 6.0*Q1*ir4 + 6.0*Q0*ir5);
  const Real Cf = 0.25*(Q4*ir + 2.0*Q3*ir2 + 9.0*Q2*ir3
                        + 21.0*Q1*ir4 + 21.0*Q0*ir5);
  const Real Ad = 3.0*(R3*ir3 + 3.0*R2*ir4 + 3.0*R1*ir5);
  const Real Bd = -(R4*ir2 + 3.0*R3*ir3 + 6.0*R2*ir4 + 6.0*R1*ir5);
  const Real Cd = 0.25*(R5*ir + 2.0*R4*ir2 + 9.0*R3*ir3
                        + 21.0*R2*ir4 + 21.0*R1*ir5);

  // angular functions (l=2, m=0 even parity; eq:Teukolsky_angular_basis):
  // f_rr = (1+3cos2th)/2, f_rth = -3 sin th cos th, f1_thth = 3 sin^2 th,
  // f2_thth = -1, f1_phph = -f1_thth, f2_phph = 1 - f_rr = 3 sin^2 th - 1
  const Real cth = z*ir;
  const Real s2 = fmax(1.0 - cth*cth, 0.0);          // sin^2 theta
  const Real sth = sqrt(s2);
  const Real frr = 0.5*(1.0 + 3.0*(1.0 - 2.0*s2));   // (1+3cos2th)/2
  const Real frt = -3.0*sth*cth;
  const Real f1tt = 3.0*s2,  f2tt = -1.0;
  const Real f1pp = -3.0*s2, f2pp = 3.0*s2 - 1.0;

  // orthonormal spherical frame components
  const Real hrr_o  = Af*frr,            hrr_d  = Ad*frr;
  const Real hrt_o  = Bf*frt,            hrt_d  = Bd*frt;
  const Real htt_o  = Cf*f1tt + Af*f2tt, htt_d  = Cd*f1tt + Ad*f2tt;
  const Real hpp_o  = Cf*f1pp + Af*f2pp, hpp_d  = Cd*f1pp + Ad*f2pp;

  // frame vectors (Cartesian); axis guard: at sth -> 0 pick e_ph = y-hat.
  // e_th = e_ph x r-hat = +theta-hat: the rh x eph order used before iter-18
  // gives -theta-hat, which flips the h_{r theta} cross term and violates the
  // linearized vacuum equations at O(B) (5% of |Ricci|; root-caused by the
  // exact-reference program, ledger N12: vacuum residual 1.8e-4 -> 1.9e-34
  // under the flip at 30-digit precision).
  Real rh[3] = {x*ir, y*ir, z*ir};
  Real eph[3], eth[3];
  const Real rho = sqrt(fmax(x*x + y*y, 1e-300));
  if (rho > 1e-12*r) {
    eph[0] = -y/rho; eph[1] = x/rho; eph[2] = 0.0;
  } else {
    eph[0] = 0.0; eph[1] = 1.0; eph[2] = 0.0;
  }
  eth[0] = eph[1]*rh[2] - eph[2]*rh[1];
  eth[1] = eph[2]*rh[0] - eph[0]*rh[2];
  eth[2] = eph[0]*rh[1] - eph[1]*rh[0];

  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      h[a][b] = hrr_o*rh[a]*rh[b]
              + hrt_o*(rh[a]*eth[b] + eth[a]*rh[b])
              + htt_o*eth[a]*eth[b] + hpp_o*eph[a]*eph[b];
      hdot[a][b] = hrr_d*rh[a]*rh[b]
                 + hrt_d*(rh[a]*eth[b] + eth[a]*rh[b])
                 + htt_d*eth[a]*eth[b] + hpp_d*eph[a]*eph[b];
    }
  }
}

} // anonymous namespace

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem()
//! \brief Teukolsky-wave initial data on flat space
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_hist_func = &Z4cCCMTeukolskyErrors;
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pz4c == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "z4c_ccm_teukolsky requires a <z4c> block" << std::endl;
    exit(EXIT_FAILURE);
  }

  const Real X  = pin->GetOrAddReal("problem", "amp", 1e-5);
  const Real rc = pin->GetOrAddReal("problem", "r_c", 4.0);
  const Real tau = pin->GetOrAddReal("problem", "tau", 1.0);

  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int ksg = indcs.ks - indcs.ng, keg = indcs.ke + indcs.ng;
  int jsg = indcs.js - indcs.ng, jeg = indcs.je + indcs.ng;
  int isg = indcs.is - indcs.ng, ieg = indcs.ie + indcs.ng;
  int nmb = pmbp->nmb_thispack;

  auto &adm = pmbp->padm->adm;
  auto &z4c = pmbp->pz4c->z4c;
  auto &u_adm_indcs = indcs;

  par_for("pgen_teukolsky", DevExeSpace(), 0, nmb-1, ksg, keg, jsg, jeg, isg, ieg,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x1v = CellCenterX(i-u_adm_indcs.is, u_adm_indcs.nx1, x1min, x1max);
    Real x2v = CellCenterX(j-u_adm_indcs.js, u_adm_indcs.nx2, x2min, x2max);
    Real x3v = CellCenterX(k-u_adm_indcs.ks, u_adm_indcs.nx3, x3min, x3max);

    Real h[3][3], hdot[3][3];
    TeukolskyMetric(0.0, x1v, x2v, x3v, X, rc, tau, h, hdot);

    for (int a = 0; a < 3; ++a) {
      for (int b = a; b < 3; ++b) {
        adm.g_dd(m,a,b,k,j,i) = (a == b ? 1.0 : 0.0) + h[a][b];
        adm.vK_dd(m,a,b,k,j,i) = -0.5*hdot[a][b];
      }
    }
    adm.alpha(m,k,j,i) = 1.0;
    z4c.alpha(m,k,j,i) = 1.0;
    z4c.vTheta(m,k,j,i) = 0.0;
  });

  switch (indcs.ng) {
    case 2: pmbp->pz4c->ADMToZ4c<2>(pmbp, pin); break;
    case 3: pmbp->pz4c->ADMToZ4c<3>(pmbp, pin); break;
    case 4: pmbp->pz4c->ADMToZ4c<4>(pmbp, pin); break;
  }
  pmbp->pz4c->Z4cToADM(pmbp);
  switch (indcs.ng) {
    case 2: pmbp->pz4c->ADMConstraints<2>(pmbp); break;
    case 3: pmbp->pz4c->ADMConstraints<3>(pmbp); break;
    case 4: pmbp->pz4c->ADMConstraints<4>(pmbp); break;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Z4cCCMTeukolskyErrors
//! \brief history: Minkowski-deviation LINF/RMS of the z4c state + constraint norms
void Z4cCCMTeukolskyErrors(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 4;
  pdata->label[0] = "LINF-dev";
  pdata->label[1] = "RMS-dev";
  pdata->label[2] = "L2-H";
  pdata->label[3] = "L2-M";

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;

  auto &z4c = pmbp->pz4c->z4c;
  auto &con = pmbp->pz4c->con;

  Real linf = 0.0, rms = 0.0, hnorm = 0.0, mnorm = 0.0;
  Kokkos::parallel_reduce("teuk_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &mx, Real &sum2, Real &hh, Real &mm_) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    Real dev = 0.0;
    for (int a = 0; a < 3; ++a) {
      for (int b = a; b < 3; ++b) {
        Real d = z4c.g_dd(m,a,b,k,j,i) - (a == b ? 1.0 : 0.0);
        dev = fmax(dev, fabs(d));
        sum2 += d*d;
        d = z4c.vA_dd(m,a,b,k,j,i);
        dev = fmax(dev, fabs(d));
        sum2 += d*d;
      }
    }
    mx = fmax(mx, dev);
    hh += SQR(con.H(m,k,j,i));
    mm_ += con.M(m,k,j,i);
  }, Kokkos::Max<Real>(linf), Kokkos::Sum<Real>(rms),
     Kokkos::Sum<Real>(hnorm), Kokkos::Sum<Real>(mnorm));

  pdata->hdata[0] = linf;
  pdata->hdata[1] = sqrt(rms/(12.0*nmkji));
  pdata->hdata[2] = sqrt(hnorm/nmkji);
  pdata->hdata[3] = sqrt(fmax(mnorm,0.0)/nmkji);
}
