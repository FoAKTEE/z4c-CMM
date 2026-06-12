#ifndef Z4C_CCM_BONDI_SOLVER_HPP_
#define Z4C_CCM_BONDI_SOLVER_HPP_
//========================================================================================
// bondi_solver — native in-process characteristic solver (N14 stage C).
// USER DIRECTIVE (2026-06-12): implement directly in AthenaK, no embedding
// bridge. This is a line-by-line port of the ZccmJl-VERIFIED algorithm
// (packages/ZccmJl: Chebyshev-CGL spectral radial scheme, LU-factorized
// hierarchy sweep, SAT RK4 evolution, hierarchy-based psi0) — the Julia
// package remains the verification harness whose admitted numbers gate
// this port (ledger N14 iters 36-42); solver R&D stays in Julia per the
// migration constraint, production lives here.
// Linear l=2 m=0 system (verified iter 31):
//   dr(r^2 Qt) = -4 r^2 dr Jt ;  dr Ut = Qt/r^2 ;
//   dr(r^2 Wt) = 2 Jt + 2 r Ut + (r^2/2) dr Ut ;
//   dr(r Ht) = -(1/2)Ut + (1/4) r dr^2 Jt - (3/2) Jt/r ;  du Jt = Ht ;
//   psi0t = -dr Jt/(2r) - dr^2 Jt/4   (hierarchy-evaluated at the tube).
// Host-only, rank-0 + broadcast; never called from device code.
//========================================================================================
#include <vector>

namespace z4c_ccm {

class BondiSolver {
 public:
  // worldtube BC scalars per slice (ringed variables)
  struct WtBC { double qr, ur, wr, hr, jr; };

  BondiSolver(double rwt, int n, double u0, double du_max);

  // initialize the J slice from analytic Teukolsky data (stub path) or
  // caller-provided ringed jr values on the grid
  void init_teukolsky(double X, double rc, double tau);

  // advance to time t with internal RK4+SAT substeps; bc(u) supplies the
  // worldtube data (analytic stub now; live sphere data next stage)
  template <typename F>
  void advance(double t, const F& bc);

  // hierarchy-based psi0 t-scalar at the worldtube for the CURRENT slice
  double psi0_worldtube(const WtBC& bcnow);

  double current_u() const { return u_; }

  // analytic Teukolsky worldtube closures (ringed; expanded in 1/r)
  static WtBC teuk_bc(double u, double rwt, double X, double rc, double tau);

 private:
  int n_;
  double rwt_, u_, du_max_, taupen_;
  std::vector<double> y_, jr_;
  std::vector<double> D_;                  // n x n differentiation matrix
  std::vector<double> Aint_;               // cumulative-integral operator
  // LU factors (with pivots) of the three sweep operators
  std::vector<double> luQ_, luU_, luW_;
  std::vector<int> pivQ_, pivU_, pivW_;

  void sweep(const std::vector<double>& jr, const WtBC& bc,
             std::vector<double>* qr, std::vector<double>* ur,
             std::vector<double>* wr, std::vector<double>* hr);
  std::vector<double> rhs(double u, const std::vector<double>& jr,
                          const WtBC& bc);
};

}  // namespace z4c_ccm
#endif  // Z4C_CCM_BONDI_SOLVER_HPP_
