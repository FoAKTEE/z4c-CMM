#ifndef Z4C_Z4C_CCM_HPP_
#define Z4C_Z4C_CCM_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file z4c_ccm.hpp
//! \brief Cauchy-characteristic matching (CCM) physical-mode injection for Z4c.
//!
//! Implements the physical rows (1-2) of the Z4c-CCM boundary-condition set
//! (z4c-CMM repository, reformulate/z4c-CMM/synthesis_z4c-CCM/
//! z4c_ccm_boundary_conditions.md): the two incoming physical modes at the
//! outer boundary are driven by a psi0 datum supplied by a characteristic
//! (CCE) evolution, transported into the Cauchy frame by the Type-II boost
//!     A = (alpha - gamma_ij beta^i s^j) e^{-2 betahat},   psi0' = A^2 psi0.
//! The injection enters the boundary RHS of the trace-free tangential part of
//! Atilde as the Bjorhus-style inhomogeneous datum (verified linearized
//! dictionary, scripts/verify_n3_dictionary.py + exact form verify_n3_exact.py):
//!     (d_t + d_s)^2 gamma^TF_AB = 4 ( psi0' mbar_A mbar_B + cc ),
//! with Atilde^TT ~ -(chi/2 alpha) d_t gamma^TT this adds
//!     rhs(Atilde_ab) += -(2 chi / alpha) * w_ab,
//!     w_ab = Re(psi0') (e_th e_th - e_ph e_ph)_ab
//!          + Im(psi0') (e_th e_ph + e_ph e_th)_ab,
//! on top of the homogeneous Sommerfeld advection (which kills the incoming
//! content; datum -> 0 reduces EXACTLY to the standard Sommerfeld scheme).
//! The constraint sector (Theta, Gam) and gauge sector keep their existing
//! boundary conditions: the injected channel is constraint-orthogonal at
//! linear order (verified, scripts/verify_n4_cpbc_compat.py).
//!
//! The psi0 datum here is an analytic provider (Gaussian pulse, uniform over
//! the boundary 2-sphere in the local frame) exercising the full machinery;
//! a CCE-coupled provider replaces Z4cCCMDatumPsi0 (the cce/ worldtube
//! infrastructure supplies the outgoing direction; reading psi0 back is the
//! documented successor hook).

#include <math.h>
#include "athena.hpp"
#include "athena_tensor.hpp"
#include "coordinates/cell_locations.hpp"
#include "z4c/z4c.hpp"

namespace z4c {

//----------------------------------------------------------------------------------------
//! \fn TeukolskyF
//! \brief n-th derivative of the Gaussian profile F(u) = X exp(-((u-rc)/tau)^2)
//! (arXiv:2308.10361 eq:Gaussian_pulse_F / eq:Teukolsky_wave_outgoing), via
//! physicists' Hermite polynomials: F^(n)(u) = X (-1/tau)^n H_n(s) e^{-s^2}.
KOKKOS_INLINE_FUNCTION
static Real TeukolskyF(const int n, const Real u, const Real X,
                       const Real rc, const Real tau) {
  const Real sarg = (u - rc)/tau;
  const Real s2 = sarg*sarg;
  Real H = 1.0;
  switch (n) {
    case 0: H = 1.0; break;
    case 1: H = 2.0*sarg; break;
    case 2: H = 4.0*s2 - 2.0; break;
    case 3: H = sarg*(8.0*s2 - 12.0); break;
    case 4: H = (16.0*s2 - 48.0)*s2 + 12.0; break;
    case 5: H = sarg*((32.0*s2 - 160.0)*s2 + 120.0); break;
    case 6: H = ((64.0*s2 - 480.0)*s2 + 720.0)*s2 - 120.0; break;
    default: H = 0.0; break;
  }
  Real pref = 1.0;
  for (int q = 0; q < n; ++q) { pref *= -1.0/tau; }
  return X*pref*H*exp(-s2);
}

//----------------------------------------------------------------------------------------
//! \fn Z4cCCMDatumPsi0
//! \brief analytic psi0 datum (CCE-frame psi0 BEFORE the Type-II boost).
//! Modes (ccm_mode):
//!   1: uniform Gaussian pulse over the sphere (machinery test);
//!   2: (l=2, m=0) pulse  psi0 = amp g(t) * +2Y20,  +2Y20 = (1/4)sqrt(15/2pi) sin^2(th)
//!      — the characteristic-pulse-injection test of arXiv:2308.10361 Sec. V.C
//!      (its eq:auto-29 prescribes a +2Y20 J-pulse on the initial null slice);
//!   3: Teukolsky self-datum — the EXACT linear-order incoming Weyl scalar of the
//!      outgoing Teukolsky wave (its eq:Teukolsky_bulk_psi0):
//!        psi0(t, r, th) = -sqrt(27 pi/10) F^(2)(t - r) * +2Y20 / r^5,
//!      with (amp, t0, sigma) = (X, r_c, tau) of the z4c_ccm_teukolsky pgen —
//!      the self-consistent CCM datum for the flat-background Teukolsky tests
//!      (Sec. V.A) at linear order, no external CCE required.
KOKKOS_INLINE_FUNCTION
static void Z4cCCMDatumPsi0(const int mode, const Real t, const Real x1,
                            const Real x2, const Real x3, const Real amp,
                            const Real t0, const Real sigma,
                            Real *re, Real *im) {
  *re = 0.0; *im = 0.0;
  const Real r = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
  const Real s2th = (r > 1e-12) ? fmax(1.0 - SQR(x3/r), 0.0) : 0.0;
  // +2Y20 = (1/4) sqrt(15/(2 pi)) sin^2 theta  (eq:auto-38)
  const Real y22 = 0.25*sqrt(15.0/(2.0*M_PI))*s2th;
  switch (mode) {
    case 1: {
      const Real arg = (t - t0)/sigma;
      *re = amp*exp(-arg*arg);
      break;
    }
    case 2: {
      const Real arg = (t - t0)/sigma;
      *re = amp*exp(-arg*arg)*y22;
      break;
    }
    case 3: {
      const Real F2 = TeukolskyF(2, t - r, amp, t0, sigma);
      *re = -sqrt(27.0*M_PI/10.0)*F2*y22/(r*r*r*r*r);
      break;
    }
    default: break;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Z4cCCMInjection
//! \brief add the CCM physical-mode datum to the boundary RHS of Atilde at (m,k,j,i).
//! Call AFTER Z4cSommerfeld at the same point. No-op when the datum vanishes.
KOKKOS_INLINE_FUNCTION
static void Z4cCCMInjection(const Z4c::Z4c_vars& z4c, const Z4c::Z4c_vars& rhs,
    const RegionIndcs &indcs, const DualArray1D<RegionSize> &size,
    const int mode, const Real tcur, const Real amp, const Real t0,
    const Real sigma, const Real betahat, const Real chi_div_floor,
    const int m, const int k, const int j, const int i) {
  // ---- boundary point coordinates and coordinate-radial direction ----------
  Real &x1min = size.d_view(m).x1min;
  Real &x1max = size.d_view(m).x1max;
  Real &x2min = size.d_view(m).x2min;
  Real &x2max = size.d_view(m).x2max;
  Real &x3min = size.d_view(m).x3min;
  Real &x3max = size.d_view(m).x3max;
  Real x1v = CellCenterX(i-indcs.is, indcs.nx1, x1min, x1max);
  Real x2v = CellCenterX(j-indcs.js, indcs.nx2, x2min, x2max);
  Real x3v = CellCenterX(k-indcs.ks, indcs.nx3, x3min, x3max);
  Real r = sqrt(SQR(x1v) + SQR(x2v) + SQR(x3v));
  Real sc[3] = {x1v/r, x2v/r, x3v/r};   // coordinate radial unit (Sommerfeld frame)

  // ---- physical 3-metric gamma_ij = gtilde_ij / chi at the point -----------
  const Real chi = fmax(z4c.chi(m,k,j,i), chi_div_floor);
  Real gam[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      gam[a][b] = z4c.g_dd(m,a,b,k,j,i) / chi;
    }
  }

  // ---- gamma-normalized outward normal s^a ---------------------------------
  Real snorm = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      snorm += gam[a][b]*sc[a]*sc[b];
    }
  }
  snorm = sqrt(fmax(snorm, 1e-30));
  Real s_u[3] = {sc[0]/snorm, sc[1]/snorm, sc[2]/snorm};

  // ---- tangential dyad (flat-frame angular vectors with axis guard) --------
  // e_ph ~ z-hat x s, e_th = s x e_ph (coordinate frame; consistent with the
  // pseudo-radial Sommerfeld frame; exact gamma-orthonormalization is part of
  // the full-GR successor item).
  Real eph[3], eth[3];
  const Real rho2 = SQR(sc[0]) + SQR(sc[1]);
  if (rho2 > 1e-12) {
    const Real irho = 1.0/sqrt(rho2);
    eph[0] = -sc[1]*irho; eph[1] = sc[0]*irho; eph[2] = 0.0;
  } else {  // on the z-axis: any tangential pair works
    eph[0] = 0.0; eph[1] = 1.0; eph[2] = 0.0;
  }
  // e_th = e_ph x s = +theta-hat (the s x eph order gives -theta-hat; the
  // real part of the injection w_ab is quadratic in e_th and unaffected,
  // but Im psi0' would flip — fixed alongside the pgen frame correction,
  // ledger N12)
  eth[0] = eph[1]*sc[2] - eph[2]*sc[1];
  eth[1] = eph[2]*sc[0] - eph[0]*sc[2];
  eth[2] = eph[0]*sc[1] - eph[1]*sc[0];

  // ---- Type-II boost: A = (alpha - gamma_ij beta^i s^j) e^{-2 betahat} -----
  const Real alpha = fmax(z4c.alpha(m,k,j,i), 1e-12);
  Real bs = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      bs += gam[a][b]*z4c.beta_u(m,a,k,j,i)*s_u[b];
    }
  }
  const Real Aboost = (alpha - bs) * exp(-2.0*betahat);

  // ---- datum and boosted psi0' ---------------------------------------------
  Real p0re, p0im;
  Z4cCCMDatumPsi0(mode, tcur, x1v, x2v, x3v, amp, t0, sigma, &p0re, &p0im);
  const Real A2 = Aboost*Aboost;
  const Real pre = A2*p0re;
  const Real pim = A2*p0im;
  if (pre == 0.0 && pim == 0.0) {
    return;  // exact Sommerfeld reduction
  }

  // ---- w_ab = (psi0' mbar mbar + cc)_ab  (TT in the tangent plane) ---------
  // mbar = (e_th - i e_ph)/sqrt(2):
  // (psi0' mbar_a mbar_b + cc) = Re psi0' (eth eth - eph eph)
  //                            + Im psi0' (eth eph + eph eth)
  const Real coef = -2.0*chi/alpha;   // rhs(Atilde) += coef * w_ab
  for (int a = 0; a < 3; ++a) {
    for (int b = a; b < 3; ++b) {
      const Real w_ab = pre*(eth[a]*eth[b] - eph[a]*eph[b])
                      + pim*(eth[a]*eph[b] + eph[a]*eth[b]);
      rhs.vA_dd(m,a,b,k,j,i) += coef*w_ab;
    }
  }
}

} // namespace z4c
#endif // Z4C_Z4C_CCM_HPP_
