//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file z4c_Sbc.cpp
//! \brief placeholder for Sommerfeld boundary condition

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "z4c/z4c.hpp"
#include "z4c/z4c_ccm.hpp"
#include <functional>
#include "z4c/ccm/bondi_solver.hpp"
#include "coordinates/cell_locations.hpp"

namespace z4c {

//----------------------------------------------------------------------------------------
//! \fn void Z4c::Z4cSommerfeld
//! \brief apply Sommerfeld BCs to the given set of points
KOKKOS_INLINE_FUNCTION
static void Z4cSommerfeld(const Z4c::Z4c_vars& z4c, const Z4c::Z4c_vars& rhs,
    const RegionIndcs &indcs, const DualArray1D<RegionSize> &size,
    const int cpbc, const int face,
    const int m, const int k, const int j, const int i) {
  // -------------------------------------------------------------------------------------
  // Scratch data
  //

  // First derivatives
  // Scalars
  AthenaPointTensor<Real, TensorSymm::NONE, 3, 1> dKhat_d;
  AthenaPointTensor<Real, TensorSymm::NONE, 3, 1> dTheta_d;

  // Vectors
  AthenaPointTensor<Real, TensorSymm::NONE, 3, 2> dGam_du;

  // Tensors
  AthenaPointTensor<Real, TensorSymm::SYM2, 3, 3> dA_ddd;


  // Psuedoradial vector
  AthenaPointTensor<Real, TensorSymm::NONE, 3, 1> s_u;

  Real idx[] = {1./size.d_view(m).dx1, 1./size.d_view(m).dx2, 1./size.d_view(m).dx3};

  // -------------------------------------------------------------------------------------
  // First derivatives
  // We force all derivatives to be calculated at second-order, as this was found to
  // be necessary for stability in Athena++.
  //
  for (int a = 0; a < 3; a++) {
    dKhat_d(a) = Dx<2>(a, idx, z4c.vKhat, m, k, j, i);
    dTheta_d(a) = Dx<2>(a, idx, z4c.vTheta, m, k, j, i);
  }
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      dGam_du(b,a) = Dx<2>(b, idx, z4c.vGam_u, m, a, k, j, i);
    }
  }
  for (int a = 0; a < 3; a++) {
    for (int b = a; b < 3; b++) {
      for (int c = 0; c < 3; c++) {
        dA_ddd(c, a, b) = Dx<2>(c, idx, z4c.vA_dd, m, a, b, k, j, i);
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // Compute psuedo-radial vector
  //
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
  s_u(0) = x1v/r;
  s_u(1) = x2v/r;
  s_u(2) = x3v/r;

  // -------------------------------------------------------------------------------------
  // Boundary RHS for scalars
  //
  rhs.vTheta(m,k,j,i) = - z4c.vTheta(m,k,j,i)/r;
  rhs.vKhat(m,k,j,i) = - sqrt(2.) * z4c.vKhat(m,k,j,i)/r;
  for (int a = 0; a < 3; a++) {
    rhs.vTheta(m,k,j,i) -= s_u(a) * dTheta_d(a);
    rhs.vKhat(m,k,j,i) -= sqrt(2.) * s_u(a) * dKhat_d(a);
  }

  // -------------------------------------------------------------------------------------
  // Boundary RHS for Gamma
  //
  for (int a = 0; a < 3; a++) {
    rhs.vGam_u(m,a,k,j,i) = - z4c.vGam_u(m, a, k, j, i)/r;
    for (int b = 0; b < 3; b++) {
      rhs.vGam_u(m,a,k,j,i) -= s_u(b) * dGam_du(b,a);
    }
  }

  // -------------------------------------------------------------------------------------
  // Boundary RHS for A_ab
  //
  for (int a = 0; a < 3; a++) {
    for (int b = a; b < 3; b++) {
      rhs.vA_dd(m,a,b,k,j,i) = - z4c.vA_dd(m,a,b,k,j,i)/r;
      for (int c = 0; c < 3; c++) {
        rhs.vA_dd(m,a,b,k,j,i) -= s_u(c) * dA_ddd(c,a,b);
      }
    }
  }

  // -------------------------------------------------------------------------------------
  // N13 constraint-preserving Gamma-tilde row (Z4c-CCM formulation
  // eq:bc-cpbc; arXiv:1010.0523v2 eq:general_CPBCs, linearized realization
  // around flat at the boundary — the regime of the paper's derivation):
  // the constraint vector Z_i = (Gamma_i - Gamma_i[gt])/2 is the advected
  // object (radiative row, same family as the other rows), while the metric
  // part dt Gamma[gt]_i = d_j(dt gt_ij) = -2 d_j(alpha At_ij) keeps its
  // evolution value. Plain Sommerfeld on Gamma reflects constraint
  // violations; this row lets them leave. All inputs are STATE derivatives
  // (race-free; second derivatives via Dxx/Dxy at order 2, like the rest of
  // this kernel).
  if (cpbc == 1) {
    // v1 constraint-damped row: stock advection + causal-rate sink on
    // Z_a = (Gamma_a - d_j gt_aj)/2 (zeroth order; stable but absorption
    // unmeasurable — ledger iter 19).
    for (int a = 0; a < 3; a++) {
      Real divg = 0.0;
      for (int b = 0; b < 3; b++) {
        divg += Dx<2>(b, idx, z4c.g_dd, m, a, b, k, j, i);
      }
      const Real Z = 0.5*(z4c.vGam_u(m,a,k,j,i) - divg);
      rhs.vGam_u(m,a,k,j,i) -= 2.0*Z/r;
    }
  } else if (cpbc == 2) {
    // O-N13-1 v2 — LEDGERED UNSTABLE (error-db iter 20: NaN t = 19 on the
    // X = 1e-5 battery): one-sided inward stencils do not rescue the
    // closure; the instability is in the row structure (no damping on the
    // Gamma-At face loop). Kept compilable as the documented attempt; do
    // not use in production. The admissible closure is to be derived from
    // the frozen-coefficient LF/GKS analysis (O-N13-1, redirected).
    // Original design: characteristic-exact advected Z row, ALL face-normal
    // derivative reads one-sided INWARD; tangential edges centered (O-N6-4).
    //   Z_a   = (Gamma_a - d_j gt_aj)/2
    //   rhs(Gamma_a) = -2 d_j(alpha At_aj) + 2 (-s.dZ_a - 2 Z_a/r)
    // (falloff power 2: Z ~ r^-2 for outgoing constraint content).
    const int axis = face/2;
    const int dirsgn = (face % 2) ? -1 : +1;   // inward index step
    const Real idn = idx[axis];
    const int dk = (axis == 2) ? dirsgn : 0;
    const int dj = (axis == 1) ? dirsgn : 0;
    const int di = (axis == 0) ? dirsgn : 0;

    // one-sided first derivative along the face axis (2nd order, inward)
#define DN1(Q, ...) (dirsgn*idn*( \
      -1.5*Q(m, ##__VA_ARGS__, k,      j,      i) \
      +2.0*Q(m, ##__VA_ARGS__, k+dk,   j+dj,   i+di) \
      -0.5*Q(m, ##__VA_ARGS__, k+2*dk, j+2*dj, i+2*di)))
    // one-sided second derivative along the face axis (2nd order)
#define DN2(Q, ...) (idn*idn*( \
       2.0*Q(m, ##__VA_ARGS__, k,      j,      i) \
      -5.0*Q(m, ##__VA_ARGS__, k+dk,   j+dj,   i+di) \
      +4.0*Q(m, ##__VA_ARGS__, k+2*dk, j+2*dj, i+2*di) \
      -1.0*Q(m, ##__VA_ARGS__, k+3*dk, j+3*dj, i+3*di)))
    // one-sided-normal derivative of a centered tangential derivative
#define DNT(Q, tdir, ...) (dirsgn*idn*( \
      -1.5*Dx<2>(tdir, idx, Q, m, ##__VA_ARGS__, k,      j,      i) \
      +2.0*Dx<2>(tdir, idx, Q, m, ##__VA_ARGS__, k+dk,   j+dj,   i+di) \
      -0.5*Dx<2>(tdir, idx, Q, m, ##__VA_ARGS__, k+2*dk, j+2*dj, i+2*di)))

    const Real alp = z4c.alpha(m,k,j,i);
    Real dal[3];
    for (int b = 0; b < 3; b++) {
      dal[b] = (b == axis) ? DN1(z4c.alpha)
                           : Dx<2>(b, idx, z4c.alpha, m, k, j, i);
    }
    for (int a = 0; a < 3; a++) {
      // first derivatives of gt_aj and At_aj (normal one-sided)
      Real dg[3][3], dA[3][3], dGam[3];
      for (int c = 0; c < 3; c++) {
        for (int b = 0; b < 3; b++) {
          dg[c][b] = (c == axis) ? DN1(z4c.g_dd, a, b)
                                 : Dx<2>(c, idx, z4c.g_dd, m, a, b, k, j, i);
          dA[c][b] = (c == axis) ? DN1(z4c.vA_dd, a, b)
                                 : Dx<2>(c, idx, z4c.vA_dd, m, a, b, k, j, i);
        }
        dGam[c] = (c == axis) ? DN1(z4c.vGam_u, a)
                              : Dx<2>(c, idx, z4c.vGam_u, m, a, k, j, i);
      }
      Real divg = 0.0;
      for (int b = 0; b < 3; b++) divg += dg[b][b];
      const Real Z = 0.5*(z4c.vGam_u(m,a,k,j,i) - divg);
      // d_c (d_j gt_aj): second derivatives with the normal rule
      Real sdZ = 0.0;
      for (int c = 0; c < 3; c++) {
        Real ddivg_c = 0.0;
        for (int b = 0; b < 3; b++) {
          if (c == axis && b == axis) {
            ddivg_c += DN2(z4c.g_dd, a, b);
          } else if (c == axis) {
            ddivg_c += DNT(z4c.g_dd, b, a, b);
          } else if (b == axis) {
            ddivg_c += DNT(z4c.g_dd, c, a, b);
          } else if (b == c) {
            ddivg_c += Dxx<2>(b, idx, z4c.g_dd, m, a, b, k, j, i);
          } else {
            ddivg_c += Dxy<2>(c, b, idx, z4c.g_dd, m, a, b, k, j, i);
          }
        }
        sdZ += s_u(c)*0.5*(dGam[c] - ddivg_c);
      }
      Real dtGamdef = 0.0;
      for (int b = 0; b < 3; b++) {
        dtGamdef += -2.0*(alp*dA[b][b]
                          + z4c.vA_dd(m,a,b,k,j,i)*dal[b]);
      }
      rhs.vGam_u(m,a,k,j,i) = dtGamdef + 2.0*(-sdZ - 2.0*Z/r);
    }
#undef DN1
#undef DN2
#undef DNT
  }
}


//---------------------------------------------------------------------------------------
//! \fn TaskStatus Z4c::Z4cBoundaryRHS
//! \brief placeholder for the Sommerfield Boundary conditions for z4c
TaskStatus Z4c::Z4cBoundaryRHS(Driver *pdriver, int stage) {
  auto &pm = pmy_pack->pmesh;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;

  int nmb = pmy_pack->nmb_thispack;
  int is = indcs.is;
  int ie = indcs.ie;
  int js = indcs.js;
  int je = indcs.je;
  int ks = indcs.ks;
  int ke = indcs.ke;

  auto &z4c_ = z4c;
  auto &rhs_ = rhs;
  bool &user_Sbc = opt.user_Sbc;

  const int cpbc_on = opt.cpbc;

  // CCM physical-mode injection (z4c_ccm.hpp); no-op unless <z4c> ccm = true
  const bool ccm_on = opt.ccm;
  const int ccm_mode = opt.ccm_mode;
  Real ccm_amp = opt.ccm_amp;
  // N14 stage C (USER DIRECTIVE: native in-AthenaK implementation): live
  // psi0 datum from the in-process characteristic solver. Rank 0 advances
  // the solver in lockstep with the Cauchy time and broadcasts the scalar;
  // mode 4 kernels read it through the amp slot (psi0 = scalar * sin^2 th).
  // The analytic Teukolsky worldtube stub gates the machinery end-to-end;
  // the live sphere-projected data path replaces it next stage.
  if (opt.ccm && opt.ccm_mode == 4) {
    Real psi0_live = 0.0;
    int rank = 0;
#if MPI_PARALLEL_ENABLED
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == 0) {
      static z4c_ccm::BondiSolver* solver = nullptr;
      const Real rwt = pm->mesh_size.x1max;
      if (solver == nullptr) {
        solver = new z4c_ccm::BondiSolver(rwt, 65, pm->time, 0.0078125);
        solver->init_teukolsky(opt.ccm_amp, opt.ccm_t0, opt.ccm_sigma);
      }
      const Real X = opt.ccm_amp, rc = opt.ccm_t0, sg = opt.ccm_sigma;
      std::function<z4c_ccm::BondiSolver::WtBC(double)> bc =
        [&](double uu) {
          return z4c_ccm::BondiSolver::teuk_bc(uu, rwt, X, rc, sg);
        };
      solver->advance(pm->time, bc);
      psi0_live = solver->psi0_worldtube(
          z4c_ccm::BondiSolver::teuk_bc(pm->time, rwt, X, rc, sg));
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&psi0_live, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif
    ccm_amp = psi0_live;
  }
  const Real ccm_t0 = opt.ccm_t0;
  const Real ccm_sigma = opt.ccm_sigma;
  const Real ccm_betahat = opt.ccm_betahat;
  const Real ccm_chifloor = opt.chi_min_floor;
  const Real tcur = pm->time;

  // We only need to apply this condition for outflow boundaries
  if (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::user
      || pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::user) {
    par_for("z4crhs_bc_x1", DevExeSpace(), 0, (nmb-1), ks, ke, js, je,
    KOKKOS_LAMBDA(int m, int k, int j) {
      // Inner boundary
      switch(mb_bcs.d_view(m,BoundaryFace::inner_x1)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 0, m, k, j, is);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, j, is); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 0, m, k, j, is);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, j, is); }
            }
          break;
        default:
          break;
      }
      // Outer boundary
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x1)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 1, m, k, j, ie);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, j, ie); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 1, m, k, j, ie);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, j, ie); }
            }
          break;
        default:
          break;
      }
    });
  }
  if (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::user
      || pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::user) {
    par_for("z4crhs_bc_x2", DevExeSpace(), 0, (nmb-1), ks, ke, is, ie,
    KOKKOS_LAMBDA(int m, int k, int i) {
      // Inner boundary
      switch(mb_bcs.d_view(m,BoundaryFace::inner_x2)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 2, m, k, js, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, js, i); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 2, m, k, js, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, js, i); }
            }
          break;
        default:
          break;
      }
      // Outer boundary
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x2)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 3, m, k, je, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, je, i); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 3, m, k, je, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, k, je, i); }
            }
          break;
        default:
          break;
      }
    });
  }
  if (pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::user
      || pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::outflow
      || pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::diode
      || pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::vacuum
      || pm->mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::user) {
    par_for("z4crhs_bc_x3", DevExeSpace(), 0, (nmb-1), js, je, is, ie,
    KOKKOS_LAMBDA(int m, int j, int i) {
      // Inner boundary
      switch(mb_bcs.d_view(m,BoundaryFace::inner_x3)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 4, m, ks, j, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, ks, j, i); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 4, m, ks, j, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, ks, j, i); }
            }
          break;
        default:
          break;
      }
      // Outer boundary
      switch (mb_bcs.d_view(m,BoundaryFace::outer_x3)) {
        case BoundaryFlag::diode:
        case BoundaryFlag::vacuum:
        case BoundaryFlag::outflow:
            Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 5, m, ke, j, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, ke, j, i); }
          break;
        case BoundaryFlag::user:
            if (user_Sbc) {
              Z4cSommerfeld(z4c_, rhs_, indcs, size, cpbc_on, 5, m, ke, j, i);
            if (ccm_on) { Z4cCCMInjection(z4c_, rhs_, indcs, size, ccm_mode, tcur, ccm_amp, ccm_t0, ccm_sigma, ccm_betahat, ccm_chifloor, m, ke, j, i); }
            }
          break;
        default:
          break;
      }
    });
  }


  return TaskStatus::complete;
}

} // end namespace z4c
