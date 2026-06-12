#ifndef Z4C_CCM_BONDI_SOLVER_HPP_
#define Z4C_CCM_BONDI_SOLVER_HPP_
//========================================================================================
// bondi_solver — native in-process characteristic solver (N14 stage C).
// USER DIRECTIVE (2026-06-12): implement directly in AthenaK, no embedding
// bridge. Line-by-line port of the ZccmJl-VERIFIED algorithm; the Julia
// package remains the verification harness whose admitted numbers gate
// this port (ledger N14); solver R&D stays in Julia per the migration
// constraint, production lives here.
//
// Iter 44 rewrite: UNRINGED beta-corrected scheme (ZccmJl
// bondisolver_unringed.jl). The live worldtube data arrives in the
// anchored gauge (beta != 0, non-decaying gauge tails); the ringed
// variables (jr = r J, ...) diverge at scri on such data (measured
// failure: psi0 off by O(10^3)). The unringed t-scalars are regular on
// y in [-1, 1] in BOTH gauges. Machine-derived system (exact-rational
// fits, test/derive_beta_hierarchy.jl; beta r-independent at this order):
//   2 Qt + (1-y) Qt' = -4 (1-y) Jt' + 24 bt
//   Ut' = Qt / (2 rwt)
//   2 Wt + (1-y) Wt' = (1-y) Jt/rwt + 2 Ut + (1-y) Ut'/2 + 4 bt (1-y)/rwt
//   Ht + (1-y) Ht' = ((1-y)/(4 rwt)) [(1-y)^2 Jt'' - 2 (1-y) Jt']
//                    + ((1-y)^2/(2 rwt)) Jt' + 3 bt (1-y)/rwt
//                    + Qt (1-y)/(4 rwt) + Ut ;     du Jt = Ht
//   psi0t(r) = -dr Jt/(2 r) - dr^2 Jt/4   (J-only; gauge-invariant)
//
// Cone labeling (iter 44 phase fix): the solver time u is the RETARDED
// label, u = t_Cauchy - r on the outgoing cones. The worldtube point of
// cone u is at Cauchy time t = u + rwt; the Cauchy outer boundary r_B is
// reached at t = u + r_B, so the boundary datum at Cauchy time t lives on
// the EARLIER cone u = t - r_B (causality lag, available in lockstep).
// The probe (set_probe/probe_query) records psi0t(r_B) per substep and
// serves it lag-interpolated; before the first cone it falls back to the
// analytic Teukolsky value (quiescent early times).
// Host-only, rank-0 + broadcast; never called from device code.
//========================================================================================
#include <deque>
#include <utility>
#include <vector>

namespace z4c_ccm {

class BondiSolver {
 public:
  // worldtube BC scalars per slice (RINGED values + beta t-scalar; the
  // solver converts to unringed t-scalars internally). beta = 0 for the
  // plain-gauge analytic stub; the live map supplies its Bw output.
  struct WtBC { double qr, ur, wr, hr, jr; double beta = 0.0; };

  BondiSolver(double rwt, int n, double u0, double du_max);

  // initialize the Jt slice from analytic Teukolsky data on cone u0 (also
  // stores (X, rc, tau) for the pre-cone probe fallback)
  void init_teukolsky(double X, double rc, double tau);

  // record psi0t(rstar) along the evolution (barycentric + D rows built
  // once); call before the first advance()
  void set_probe(double rstar);

  // advance to cone u = t with internal RK4+SAT substeps; bc(u) supplies
  // the worldtube data for cone u (tube sampled at Cauchy time u + rwt)
  template <typename F>
  void advance(double t, const F& bc);

  // lag-interpolated probe value on cone u (linear between substeps;
  // analytic Teukolsky fallback for u before the initial cone)
  double probe_query(double u) const;

  // psi0t at radius r >= rwt on the CURRENT cone (J-only, gauge-invariant)
  double psi0_at(double rstar) const;

  double current_u() const { return u_; }

  // analytic Teukolsky worldtube closures (ringed; expanded in 1/r)
  static WtBC teuk_bc(double u, double rwt, double X, double rc, double tau);

 private:
  int n_;
  double rwt_, u_, du_max_, taupen_;
  double X_ = 0.0, rc_ = 0.0, tau_ = 1.0;
  bool have_teuk_ = false;
  double probe_r_ = 0.0;
  std::vector<double> probe_row_, probe_rowD_, probe_rowD2_;
  std::deque<std::pair<double, double>> hist_;
  std::vector<double> y_, Jt_;
  std::vector<double> D_;                  // n x n differentiation matrix
  // LU factors (with pivots): QW = 2I + (1-y)D (shared by Q and W),
  // U = D, H = I + (1-y)D — all with the tube BC row
  std::vector<double> luQW_, luU_, luH_;
  std::vector<int> pivQW_, pivU_, pivH_;

  void sweep(const std::vector<double>& Jt, const WtBC& bc,
             std::vector<double>* Qt, std::vector<double>* Ut,
             std::vector<double>* Wt, std::vector<double>* Ht) const;
  std::vector<double> rhs(double u, const std::vector<double>& Jt,
                          const WtBC& bc) const;
  void record_probe();
  std::vector<double> bary_row(double rstar) const;
};

}  // namespace z4c_ccm
#endif  // Z4C_CCM_BONDI_SOLVER_HPP_
