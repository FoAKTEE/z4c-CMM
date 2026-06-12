#ifndef Z4C_CCM_WORLDTUBE_DATA_HPP_
#define Z4C_CCM_WORLDTUBE_DATA_HPP_
//========================================================================================
// worldtube_data — live l=2 m=0 worldtube data for the native CCM solver
// (N14 stage C part 2). Samples the Cauchy ADM state on three spheres
// (rs - dr, rs, rs + dr) with the geodesic-grid interpolator, projects the
// data-contract scalars (n14_stageB_design.md):
//   hTT  = (h_thth^ - h_phph^)/(2 s^2)      [s^2 projection]
//   hrr, trace                               [P = 3c^2-1 projections]
//   hrth = r h_{r theta}^                    [s c projection]
//   dt_* from K (linear: dt h = -2 alpha K), dr_trace by centered diff.
// The theta-derivative entries of the contract are ANALYTIC at l=2:
//   dth_trace_sc = -6 trace, dth_dr_trace_sc = -6 dr_trace,
//   ethb_hrth_P = hrth  (shape algebra, verified iters 39-42).
// worldtube_map (C++ port of the ZccmJl-verified six-row map) produces the
// ringed BCs. Host-only; collective (all ranks call; reduction inside).
//========================================================================================
#include "z4c/ccm/bondi_solver.hpp"

namespace z4c {  // fwd
class MeshBlockPack;
}

namespace z4c_ccm {

struct WtScalars {
  double hTT, dt_hTT, hrr, trace, dr_trace, hrth, dt_hrth;
};

// sample + project the live scalars at radius rs (centered dr stencil).
// Collective over MPI; the returned values are valid on every rank.
WtScalars sample_worldtube(void* pmbp_, double rs, double dr);

// the verified six-row local map -> ringed BCs (+ jr Dirichlet value)
BondiSolver::WtBC worldtube_map(double rwt, const WtScalars& s);

}  // namespace z4c_ccm
#endif  // Z4C_CCM_WORLDTUBE_DATA_HPP_
