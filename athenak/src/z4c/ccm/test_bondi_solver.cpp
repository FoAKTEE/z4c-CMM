//========================================================================================
// test_bondi_solver.cpp — standalone gate of the C++ unringed solver port
// against the ZccmJl-admitted numbers (iter 44):
//   plain-BC evolved psi0 probe(42.025) rel ~ 8.8e-5, tube rel ~ 0.06
//   anchored (map-BC) end-to-end: same class (8.8e-5 / 0.06)
// Build: g++ -O2 -std=c++17 -o /tmp/test_bondi test_bondi_solver.cpp bondi_solver.cpp
//========================================================================================
#include <cmath>
#include <cstdio>
#include <functional>

#include "bondi_solver.hpp"

using z4c_ccm::BondiSolver;

static double teukF(int n, double u, double X, double rc, double tau) {
  const double s = (u - rc)/tau;
  const double s2 = s*s;
  double H = 1.0;
  switch (n) {
    case 0: H = 1.0; break;
    case 1: H = 2.0*s; break;
    case 2: H = 4.0*s2 - 2.0; break;
    case 3: H = s*(8.0*s2 - 12.0); break;
    case 4: H = (16.0*s2 - 48.0)*s2 + 12.0; break;
    case 5: H = s*((32.0*s2 - 160.0)*s2 + 120.0); break;
    default: H = 0.0; break;
  }
  double pref = 1.0;
  for (int q = 0; q < n; ++q) pref *= -1.0/tau;
  return X*pref*H*std::exp(-s2);
}

int main() {
  const double X = 1.0e-5, rc = 20.0, tau = 2.0, rwt = 41.0;
  auto want = [&](double rs) {
    return -9.0/8.0*teukF(2, 20.0, X, rc, tau)/std::pow(rs, 5);
  };
  // ---- plain-BC evolution (Julia V4) ----
  {
    BondiSolver s(rwt, 65, 16.0, 0.00078125);
    s.init_teukolsky(X, rc, tau);
    s.set_probe(42.025);
    std::function<BondiSolver::WtBC(double)> bc = [&](double uu) {
      return BondiSolver::teuk_bc(uu, rwt, X, rc, tau);
    };
    s.advance(20.0, bc);
    const double pw = s.psi0_at(rwt), pp = s.psi0_at(42.025);
    std::printf("plain: tube rel = %.6e  probe rel = %.6e\n",
                std::fabs(pw - want(rwt))/std::fabs(want(rwt)),
                std::fabs(pp - want(42.025))/std::fabs(want(42.025)));
    std::printf("plain: probe_query(20) rel = %.6e\n",
                std::fabs(s.probe_query(20.0) - want(42.025))
                /std::fabs(want(42.025)));
  }
  // ---- Sec V.C pulse, quiescent BCs: gate vs the ZccmJl reference ----
  // (results/numerical/n14_pulse_ref.csv: rwt=40, rB=41, Z=1e-3, n=65,
  // du=7.8125e-4; trough psi0(t=52.50) = -1.718782e-6, peak |psi0| =
  // 2.711123e-6 at t=75.47, conditioning 7.019e-3)
  {
    const double rw = 40.0, rB = 41.0;
    BondiSolver s(rw, 65, -rw, 0.00078125);
    s.init_pulse(1.0e-3);
    s.set_probe(rB);
    std::function<BondiSolver::WtBC(double)> bc =
        [](double) { return BondiSolver::WtBC{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; };
    // production pattern: advance to cone t - rwt, then query the lagged
    // cone t - rB (the history is pruned to the lag window)
    s.advance(52.50 - rw, bc);
    const double v_trough = s.probe_query(52.50 - rB);
    s.advance(75.47 - rw, bc);
    const double v_peak = s.probe_query(75.47 - rB);
    std::printf("pulse: psi0(t=52.50) = %.6e (ref -1.718782e-6, rel %.2e)\n",
                v_trough, std::fabs(v_trough + 1.718782e-6)/1.718782e-6);
    std::printf("pulse: psi0(t=75.47) = %.6e (ref peak 2.711123e-6, rel %.2e)\n",
                v_peak, std::fabs(std::fabs(v_peak) - 2.711123e-6)/2.711123e-6);
  }
  // ---- anchored map-BC evolution (Julia V5) ----
  {
    BondiSolver s(rwt, 65, 8.0, 0.00078125);
    s.set_probe(42.025);
    std::function<BondiSolver::WtBC(double)> bc = [&](double u) {
      auto A = [&](int n) { return teukF(n, u, X, rc, tau); };
      const double r = rwt;
      const double Aser = 3*(A(2)/std::pow(r,3) + 3*A(1)/std::pow(r,4)
                             + 3*A(0)/std::pow(r,5));
      const double drA = 3*(-3*A(2)/std::pow(r,4) - 12*A(1)/std::pow(r,5)
                            - 15*A(0)/std::pow(r,6));
      const double Bser = -(A(3)/std::pow(r,2) + 3*A(2)/std::pow(r,3)
                            + 6*A(1)/std::pow(r,4) + 6*A(0)/std::pow(r,5));
      const double Cser = (A(4)/r + 2*A(3)/std::pow(r,2) + 9*A(2)/std::pow(r,3)
                           + 21*A(1)/std::pow(r,4) + 21*A(0)/std::pow(r,5))/4;
      const double dtB = -(A(4)/std::pow(r,2) + 3*A(3)/std::pow(r,3)
                           + 6*A(2)/std::pow(r,4) + 6*A(1)/std::pow(r,5));
      const double dtC = (A(5)/r + 2*A(4)/std::pow(r,2) + 9*A(3)/std::pow(r,3)
                          + 21*A(2)/std::pow(r,4) + 21*A(1)/std::pow(r,5))/4;
      const double dtA = 3*(A(3)/std::pow(r,3) + 3*A(2)/std::pow(r,4)
                            + 3*A(1)/std::pow(r,5));
      const double hTT = 1.5*(2*Cser - Aser);
      const double dt_hTT = 1.5*(2*dtC - dtA);
      const double hrr = Aser, trace = -Aser, dr_trace = -drA;
      const double dth_trace_sc = -6.0*trace;
      const double dth_dr_trace_sc = -6.0*dr_trace;
      const double hrth = -3.0*Bser*r, dt_hrth = -3.0*dtB*r;
      const double ethb_hrth_P = hrth;
      BondiSolver::WtBC b;
      b.jr = r*hTT;
      b.hr = r*dt_hTT;
      b.ur = r*r*(-dth_trace_sc/(4.0*r));
      b.qr = r*(-dth_trace_sc/4.0 - r*dth_dr_trace_sc/4.0
                - 3.0*hrth/r + dt_hrth);
      b.wr = r*r*(3.0*hrr/(4.0*r) - dr_trace/4.0);
      b.beta = hrr/4.0 - trace/8.0 - r*dr_trace/8.0 + ethb_hrth_P/(4.0*r);
      return b;
    };
    s.advance(20.0, bc);
    const double pw = s.psi0_at(rwt), pp = s.psi0_at(42.025);
    std::printf("anchored: tube rel = %.6e  probe rel = %.6e\n",
                std::fabs(pw - want(rwt))/std::fabs(want(rwt)),
                std::fabs(pp - want(42.025))/std::fabs(want(42.025)));
  }
  return 0;
}
