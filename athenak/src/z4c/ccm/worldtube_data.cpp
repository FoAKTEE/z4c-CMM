//========================================================================================
// worldtube_data.cpp — see worldtube_data.hpp.
//========================================================================================
#include <cmath>
#include <memory>
#include <vector>

#include "athena.hpp"
#include "geodesic-grid/spherical_grid.hpp"
#include "coordinates/adm.hpp"
#include "mesh/mesh.hpp"
#include "z4c/ccm/worldtube_data.hpp"

namespace z4c_ccm {

namespace {

// shape inner products on the sphere (m = 0):
//   <s^2, s^2> dOmega = 2pi * 32/15 ;  <P, P> = 2pi * 16/15 ... computed
//   numerically per-grid from the quadrature for consistency.

struct SphereSet {
  std::vector<std::unique_ptr<SphericalGrid>> grids;  // rs-dr, rs, rs+dr
  double rs, dr;
};

SphereSet* sset = nullptr;

// project the four contract scalars from interpolated ADM data on grid g.
// vals layout: interp_vals(angle, var) for the I_ADM_* range interpolated.
// dt h = -2 K at linear order (alpha = 1 + O(h); the ADM alpha slot of
// u_adm is NOT filled by Z4cToADM — reading it was the iter-44 NaN bug).
// Angles not owned by this rank carry zero-filled interp values (gamma = 0
// would alias to h = -delta): skip them in the data sums — detected by
// gxx == 0 exactly (a real metric value is ~1; interp_indcs is private).
// The shape norms run over ALL angles (identical on every rank), so the
// MPI_Allreduce SUM assembles the full projection.
void project(SphericalGrid* g, double radius,
             const DualArray2D<Real>& vals, int off_g, int off_k,
             double* hTT, double* hrr, double* trace,
             double* hrth, double* dt_hTT, double* dt_hrth) {
  double n_s2 = 0.0, n_P = 0.0, n_sc = 0.0;
  double aTT = 0.0, arr = 0.0, atr = 0.0, art = 0.0, dTT = 0.0, drt = 0.0;
  const int na = g->nangles;
  for (int ip = 0; ip < na; ++ip) {
    const Real th = g->polar_pos.h_view(ip, 0);
    const Real wgt = g->solid_angles.h_view(ip);
    const Real st = std::sin(th), ct = std::cos(th);
    // shapes + norms (global on every rank)
    const Real s2 = st*st;
    const Real P = 3.0*ct*ct - 1.0;
    const Real sc = st*ct;
    n_s2 += wgt*s2*s2; n_P += wgt*P*P; n_sc += wgt*sc*sc;
    if (vals.h_view(ip, off_g+0) == 0.0) continue;       // off-rank angle
    const Real ph = g->polar_pos.h_view(ip, 1);
    const Real sp = std::sin(ph), cp = std::cos(ph);
    // unit frame vectors
    const Real rh[3] = {st*cp, st*sp, ct};
    const Real thh[3] = {ct*cp, ct*sp, -st};
    const Real phh[3] = {-sp, cp, 0.0};
    // gamma - delta and K in Cartesian
    Real hc[3][3], kc[3][3];
    const Real gxx = vals.h_view(ip, off_g+0), gxy = vals.h_view(ip, off_g+1);
    const Real gxz = vals.h_view(ip, off_g+2), gyy = vals.h_view(ip, off_g+3);
    const Real gyz = vals.h_view(ip, off_g+4), gzz = vals.h_view(ip, off_g+5);
    const Real kxx = vals.h_view(ip, off_k+0), kxy = vals.h_view(ip, off_k+1);
    const Real kxz = vals.h_view(ip, off_k+2), kyy = vals.h_view(ip, off_k+3);
    const Real kyz = vals.h_view(ip, off_k+4), kzz = vals.h_view(ip, off_k+5);
    hc[0][0] = gxx - 1.0; hc[0][1] = hc[1][0] = gxy; hc[0][2] = hc[2][0] = gxz;
    hc[1][1] = gyy - 1.0; hc[1][2] = hc[2][1] = gyz; hc[2][2] = gzz - 1.0;
    kc[0][0] = kxx; kc[0][1] = kc[1][0] = kxy; kc[0][2] = kc[2][0] = kxz;
    kc[1][1] = kyy; kc[1][2] = kc[2][1] = kyz; kc[2][2] = kzz;
    auto proj2 = [&](const Real a[3], const Real b[3], Real m[3][3]) {
      Real acc = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) acc += a[i]*m[i][j]*b[j];
      return acc;
    };
    const Real h_tt = proj2(thh, thh, hc), h_pp = proj2(phh, phh, hc);
    const Real h_rt = proj2(rh, thh, hc);
    const Real h_rr = proj2(rh, rh, hc);
    // dt h = -2 K (linear order)
    const Real dth_tt = -2.0*proj2(thh, thh, kc);
    const Real dth_pp = -2.0*proj2(phh, phh, kc);
    const Real dth_rt = -2.0*proj2(rh, thh, kc);
    aTT += wgt*0.5*(h_tt - h_pp)*s2;
    dTT += wgt*0.5*(dth_tt - dth_pp)*s2;
    arr += wgt*h_rr*P;
    atr += wgt*(h_tt + h_pp)*P;
    art += wgt*(radius*h_rt)*sc;        // coordinate h_{r theta} = r h_rt^
    drt += wgt*(radius*dth_rt)*sc;
  }
  *hTT = aTT/n_s2; *dt_hTT = dTT/n_s2;
  *hrr = arr/n_P; *trace = atr/n_P;
  *hrth = art/n_sc; *dt_hrth = drt/n_sc;
}

}  // namespace

WtScalars sample_worldtube(void* pmbp_, double rs, double dr) {
  auto* pmbp = static_cast<MeshBlockPack*>(pmbp_);
  if (sset == nullptr || sset->rs != rs || sset->dr != dr) {
    delete sset;
    sset = new SphereSet();
    sset->rs = rs; sset->dr = dr;
    for (int k = -1; k <= 1; ++k) {
      sset->grids.emplace_back(
          std::make_unique<SphericalGrid>(pmbp, 8, rs + k*dr));
    }
  }
  const int IG = pmbp->padm->I_ADM_GXX;     // 6 metric comps from here
  const int IK = pmbp->padm->I_ADM_KXX;     // 6 K comps
  WtScalars out{};
  double tr3[3];
  for (int k = 0; k < 3; ++k) {
    auto& g = sset->grids[k];
    g->InterpolateToSphere(IK + 6, pmbp->padm->u_adm);   // g_ij + K_ij
    double hTT, hrr, trc, hrth, dTT, drt;
    project(g.get(), sset->rs + (k-1)*sset->dr, g->interp_vals,
            IG, IK, &hTT, &hrr, &trc, &hrth, &dTT, &drt);
    tr3[k] = trc;
    if (k == 1) {
      out.hTT = hTT; out.hrr = hrr; out.trace = trc;
      out.hrth = hrth; out.dt_hTT = dTT; out.dt_hrth = drt;
    }
  }
  out.dr_trace = (tr3[2] - tr3[0])/(2.0*sset->dr);
#if MPI_PARALLEL_ENABLED
  // interp_vals are already globally assembled by InterpolateToSphere?
  // The geodesic-grid interpolator fills only local points; reduce the
  // projected scalars (off-rank points contribute zero).
  double buf[7] = {out.hTT, out.dt_hTT, out.hrr, out.trace, out.dr_trace,
                   out.hrth, out.dt_hrth};
  MPI_Allreduce(MPI_IN_PLACE, buf, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  out = {buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6]};
#endif
  return out;
}

BondiSolver::WtBC worldtube_map(double rwt, const WtScalars& s) {
  // the verified six-row local map (iters 39-42); theta-derivative entries
  // analytic at l=2: dth_trace_sc = -6 trace, dth_dr_trace_sc = -6 dr_trace,
  // ethb_hrth_P = hrth
  const double dth_trace_sc = -6.0*s.trace;
  const double dth_dr_trace_sc = -6.0*s.dr_trace;
  const double ethb_hrth_P = s.hrth;
  BondiSolver::WtBC b;
  b.jr = rwt*s.hTT;
  b.hr = rwt*s.dt_hTT;
  b.ur = rwt*rwt*(-dth_trace_sc/(4.0*rwt));
  b.qr = rwt*(-dth_trace_sc/4.0 - rwt*dth_dr_trace_sc/4.0
              - 3.0*s.hrth/rwt + s.dt_hrth);
  b.wr = rwt*rwt*(3.0*s.hrr/(4.0*rwt) - s.dr_trace/4.0);
  // beta t-scalar (anchored gauge; feeds the beta-corrected sweep sources)
  b.beta = s.hrr/4.0 - s.trace/8.0 - rwt*s.dr_trace/8.0
           + ethb_hrth_P/(4.0*rwt);
  return b;
}

}  // namespace z4c_ccm
