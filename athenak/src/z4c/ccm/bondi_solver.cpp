//========================================================================================
// bondi_solver.cpp — see bondi_solver.hpp. Port gated against the
// ZccmJl-admitted numbers (N14 ledger, iter 44 unringed scheme):
// exact-slice sweeps machine-exact on BOTH gauges; anchored end-to-end
// psi0 probe rel 8.8e-5 at du = 7.8125e-4 (7.1e-6 at du/2,
// truncation-dominated); tube value at the documented ~0.05 floor class.
//========================================================================================
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "bondi_solver.hpp"

namespace z4c_ccm {

namespace {

// physicists' Hermite F^(n) for the Gaussian pulse (matches z4c_ccm.hpp)
double Fn(int n, double u, double X, double rc, double tau) {
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
    case 6: H = ((64.0*s2 - 480.0)*s2 + 720.0)*s2 - 120.0; break;
    default: H = 0.0; break;
  }
  double pref = 1.0;
  for (int q = 0; q < n; ++q) pref *= -1.0/tau;
  return X*pref*H*std::exp(-s2);
}

// dense LU with partial pivoting (Doolittle); A is n x n row-major
void lu_factor(std::vector<double>* A, std::vector<int>* piv, int n) {
  piv->resize(n);
  for (int i = 0; i < n; ++i) (*piv)[i] = i;
  for (int k = 0; k < n; ++k) {
    int p = k;
    double mx = std::fabs((*A)[k*n + k]);
    for (int i = k+1; i < n; ++i) {
      if (std::fabs((*A)[i*n + k]) > mx) { mx = std::fabs((*A)[i*n + k]); p = i; }
    }
    if (p != k) {
      for (int j = 0; j < n; ++j) std::swap((*A)[k*n + j], (*A)[p*n + j]);
      std::swap((*piv)[k], (*piv)[p]);
    }
    for (int i = k+1; i < n; ++i) {
      (*A)[i*n + k] /= (*A)[k*n + k];
      const double lik = (*A)[i*n + k];
      for (int j = k+1; j < n; ++j) (*A)[i*n + j] -= lik*(*A)[k*n + j];
    }
  }
}

void lu_solve(const std::vector<double>& LU, const std::vector<int>& piv,
              int n, std::vector<double>* b) {
  std::vector<double> x(n);
  for (int i = 0; i < n; ++i) x[i] = (*b)[piv[i]];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) x[i] -= LU[i*n + j]*x[j];
  }
  for (int i = n-1; i >= 0; --i) {
    for (int j = i+1; j < n; ++j) x[i] -= LU[i*n + j]*x[j];
    x[i] /= LU[i*n + i];
  }
  *b = x;
}

void matvec(const std::vector<double>& A, const std::vector<double>& x,
            int n, std::vector<double>* out) {
  out->assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < n; ++j) acc += A[i*n + j]*x[j];
    (*out)[i] = acc;
  }
}

}  // namespace

BondiSolver::BondiSolver(double rwt, int n, double u0, double du_max)
    : n_(n), rwt_(rwt), u_(u0), du_max_(du_max) {
  const int N = n - 1;
  y_.resize(n);
  for (int k = 0; k <= N; ++k) y_[k] = -std::cos(M_PI*k/N);
  // Trefethen D on ascending CGL nodes
  std::vector<double> c(n, 1.0);
  c[0] = 2.0; c[n-1] = 2.0;
  for (int k = 0; k < n; ++k) c[k] *= ((k % 2 == 0) ? 1.0 : -1.0);
  D_.assign(n*n, 0.0);
  for (int i = 0; i < n; ++i) {
    double rowsum = 0.0;
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        D_[i*n + j] = (c[i]/c[j])/(y_[i] - y_[j]);
        rowsum += D_[i*n + j];
      }
    }
    D_[i*n + i] = -rowsum;
  }
  // unringed sweep operators (scri rows automatic — regular limits):
  //   QW = 2I + (1-y)D (Q and W share it), U = D, H = I + (1-y)D
  // each with the tube BC row (row 0 -> e0)
  auto build = [&](std::vector<double>* M, double diag, bool oyD) {
    M->assign(n*n, 0.0);
    for (int i = 0; i < n; ++i) {
      const double oy = 1.0 - y_[i];
      for (int j = 0; j < n; ++j) {
        (*M)[i*n + j] = (oyD ? oy : 1.0)*D_[i*n + j]
                        + ((i == j) ? diag : 0.0);
      }
    }
    for (int j = 0; j < n; ++j) (*M)[0*n + j] = 0.0;
    (*M)[0] = 1.0;
  };
  build(&luQW_, 2.0, true);  lu_factor(&luQW_, &pivQW_, n);
  build(&luU_, 0.0, false);  lu_factor(&luU_, &pivU_, n);
  build(&luH_, 1.0, true);   lu_factor(&luH_, &pivH_, n);
  taupen_ = (1.0/(2.0*rwt))*N*N;       // SAT: sigma=1, v_bdry * N^2
  Jt_.assign(n, 0.0);
}

void BondiSolver::init_teukolsky(double X, double rc, double tau) {
  X_ = X; rc_ = rc; tau_ = tau; have_teuk_ = true;
  for (int i = 0; i < n_; ++i) {
    const double ir = (1.0 - y_[i])/(2.0*rwt_);          // 1/r (0 at scri)
    Jt_[i] = 0.75*(Fn(4, u_, X, rc, tau)*ir
                   + Fn(2, u_, X, rc, tau)*ir*ir*ir);
  }
}

std::vector<double> BondiSolver::bary_row(double rstar) const {
  const int n = n_;
  const double ystar = 1.0 - 2.0*rwt_/rstar;
  std::vector<double> row(n, 0.0);
  for (int k = 0; k < n; ++k) {
    if (y_[k] == ystar) { row[k] = 1.0; return row; }
  }
  double ssum = 0.0;
  for (int k = 0; k < n; ++k) {
    double w = (k % 2 == 0) ? 1.0 : -1.0;
    if (k == 0 || k == n-1) w *= 0.5;
    row[k] = w/(ystar - y_[k]);
    ssum += row[k];
  }
  for (int k = 0; k < n; ++k) row[k] /= ssum;
  return row;
}

void BondiSolver::set_probe(double rstar) {
  probe_r_ = rstar;
  probe_row_ = bary_row(rstar);
  // rowD = row * D, rowD2 = rowD * D (vector-matrix products, built once)
  probe_rowD_.assign(n_, 0.0);
  probe_rowD2_.assign(n_, 0.0);
  for (int j = 0; j < n_; ++j) {
    double a1 = 0.0;
    for (int i = 0; i < n_; ++i) a1 += probe_row_[i]*D_[i*n_ + j];
    probe_rowD_[j] = a1;
  }
  for (int j = 0; j < n_; ++j) {
    double a2 = 0.0;
    for (int i = 0; i < n_; ++i) a2 += probe_rowD_[i]*D_[i*n_ + j];
    probe_rowD2_[j] = a2;
  }
  hist_.clear();
  record_probe();
}

double BondiSolver::psi0_at(double rstar) const {
  std::vector<double> djt, d2jt;
  matvec(D_, Jt_, n_, &djt);
  matvec(D_, djt, n_, &d2jt);
  const std::vector<double> row = bary_row(rstar);
  double j1 = 0.0, j2 = 0.0;
  for (int k = 0; k < n_; ++k) { j1 += row[k]*djt[k]; j2 += row[k]*d2jt[k]; }
  const double oy = 2.0*rwt_/rstar;                      // 1 - ystar
  const double drJ = oy*oy/(2.0*rwt_)*j1;
  const double d2rJ = oy*oy/(4.0*rwt_*rwt_)*(oy*oy*j2 - 2.0*oy*j1);
  return -drJ/(2.0*rstar) - d2rJ/4.0;
}

void BondiSolver::record_probe() {
  if (probe_r_ <= 0.0) return;
  double j1 = 0.0, j2 = 0.0;
  for (int k = 0; k < n_; ++k) {
    j1 += probe_rowD_[k]*Jt_[k];
    j2 += probe_rowD2_[k]*Jt_[k];
  }
  const double oy = 2.0*rwt_/probe_r_;
  const double drJ = oy*oy/(2.0*rwt_)*j1;
  const double d2rJ = oy*oy/(4.0*rwt_*rwt_)*(oy*oy*j2 - 2.0*oy*j1);
  hist_.emplace_back(u_, -drJ/(2.0*probe_r_) - d2rJ/4.0);
  // prune far-past entries (lag never exceeds r_B - rwt plus margin)
  const double keep = (probe_r_ - rwt_) + 16.0*du_max_ + 1.0;
  while (hist_.size() > 2 && hist_.front().first < u_ - keep) {
    hist_.pop_front();
  }
}

double BondiSolver::probe_query(double u) const {
  if (hist_.empty() || u <= hist_.front().first) {
    // before the recorded cones: analytic Teukolsky fallback (the run
    // starts quiescent; the retarded phase on cone u at any radius is u)
    if (have_teuk_ && probe_r_ > 0.0) {
      const double r5 = std::pow(probe_r_, 5);
      return -1.125*Fn(2, u, X_, rc_, tau_)/r5;
    }
    return hist_.empty() ? 0.0 : hist_.front().second;
  }
  if (u >= hist_.back().first) return hist_.back().second;
  // binary search for the bracketing pair
  size_t lo = 0, hi = hist_.size() - 1;
  while (hi - lo > 1) {
    const size_t mid = (lo + hi)/2;
    if (hist_[mid].first <= u) lo = mid; else hi = mid;
  }
  const double u0 = hist_[lo].first, u1 = hist_[hi].first;
  const double f = (u1 > u0) ? (u - u0)/(u1 - u0) : 0.0;
  return (1.0 - f)*hist_[lo].second + f*hist_[hi].second;
}

BondiSolver::WtBC BondiSolver::teuk_bc(double u, double rwt, double X,
                                       double rc, double tau) {
  const double ir = 1.0/rwt;
  WtBC b;
  b.jr = 0.75*(Fn(4, u, X, rc, tau) + Fn(2, u, X, rc, tau)*ir*ir);
  b.qr = 3.0*Fn(4, u, X, rc, tau) - 9.0*Fn(3, u, X, rc, tau)*ir
         - 9.0*Fn(2, u, X, rc, tau)*ir*ir;
  b.ur = -1.5*Fn(4, u, X, rc, tau) + 3.0*Fn(3, u, X, rc, tau)*ir
         + 2.25*Fn(2, u, X, rc, tau)*ir*ir;
  b.wr = -1.5*Fn(4, u, X, rc, tau) - 1.5*Fn(3, u, X, rc, tau)*ir
         - 0.75*Fn(2, u, X, rc, tau)*ir*ir;
  b.hr = 0.75*(Fn(5, u, X, rc, tau) + Fn(3, u, X, rc, tau)*ir*ir);
  b.beta = 0.0;                                   // plain Bondi gauge
  return b;
}

void BondiSolver::sweep(const std::vector<double>& Jt, const WtBC& bc,
                        std::vector<double>* Qt, std::vector<double>* Ut,
                        std::vector<double>* Wt,
                        std::vector<double>* Ht) const {
  const int n = n_;
  const double bt = bc.beta;
  std::vector<double> djt, d2jt;
  matvec(D_, Jt, n, &djt);
  matvec(D_, djt, n, &d2jt);
  // Q: 2 Qt + (1-y) Qt' = -4 (1-y) Jt' + 24 bt
  std::vector<double> b(n);
  for (int i = 0; i < n; ++i) {
    b[i] = -4.0*(1.0 - y_[i])*djt[i] + 24.0*bt;
  }
  b[0] = bc.qr/rwt_;
  lu_solve(luQW_, pivQW_, n, &b);
  *Qt = b;
  // U: Ut' = Qt/(2 rwt)
  for (int i = 0; i < n; ++i) b[i] = (*Qt)[i]/(2.0*rwt_);
  b[0] = bc.ur/(rwt_*rwt_);
  lu_solve(luU_, pivU_, n, &b);
  *Ut = b;
  std::vector<double> dut;
  matvec(D_, *Ut, n, &dut);
  // W: 2 Wt + (1-y) Wt' = (1-y) Jt/rwt + 2 Ut + (1-y) Ut'/2 + 4 bt (1-y)/rwt
  for (int i = 0; i < n; ++i) {
    const double oy = 1.0 - y_[i];
    b[i] = oy*Jt[i]/rwt_ + 2.0*(*Ut)[i] + oy*dut[i]/2.0 + 4.0*bt*oy/rwt_;
  }
  b[0] = bc.wr/(rwt_*rwt_);
  lu_solve(luQW_, pivQW_, n, &b);
  *Wt = b;
  // H: Ht + (1-y) Ht' = ((1-y)/(4 rwt))[(1-y)^2 Jt'' - 2 (1-y) Jt']
  //    + ((1-y)^2/(2 rwt)) Jt' + 3 bt (1-y)/rwt + Qt (1-y)/(4 rwt) + Ut
  for (int i = 0; i < n; ++i) {
    const double oy = 1.0 - y_[i];
    b[i] = oy/(4.0*rwt_)*(oy*oy*d2jt[i] - 2.0*oy*djt[i])
           + oy*oy/(2.0*rwt_)*djt[i] + 3.0*bt*oy/rwt_
           + (*Qt)[i]*oy/(4.0*rwt_) + (*Ut)[i];
  }
  b[0] = bc.hr/rwt_;
  lu_solve(luH_, pivH_, n, &b);
  *Ht = b;
}

std::vector<double> BondiSolver::rhs(double u, const std::vector<double>& Jt,
                                     const WtBC& bc) const {
  std::vector<double> Qt, Ut, Wt, Ht;
  sweep(Jt, bc, &Qt, &Ut, &Wt, &Ht);
  Ht[0] += taupen_*(bc.jr/rwt_ - Jt[0]);                  // SAT (sigma = 1)
  return Ht;
}

template <typename F>
void BondiSolver::advance(double t, const F& bc) {
  while (u_ < t - 1e-12) {
    const double du = std::fmin(du_max_, t - u_);
    const WtBC b0 = bc(u_), bh = bc(u_ + du/2.0), b1 = bc(u_ + du);
    std::vector<double> k1 = rhs(u_, Jt_, b0);
    std::vector<double> tmp(n_);
    for (int i = 0; i < n_; ++i) tmp[i] = Jt_[i] + du/2.0*k1[i];
    std::vector<double> k2 = rhs(u_ + du/2.0, tmp, bh);
    for (int i = 0; i < n_; ++i) tmp[i] = Jt_[i] + du/2.0*k2[i];
    std::vector<double> k3 = rhs(u_ + du/2.0, tmp, bh);
    for (int i = 0; i < n_; ++i) tmp[i] = Jt_[i] + du*k3[i];
    std::vector<double> k4 = rhs(u_ + du, tmp, b1);
    for (int i = 0; i < n_; ++i) {
      Jt_[i] += du/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    u_ += du;
    record_probe();
  }
}
// explicit instantiation for the std::function closure used by callers
template void BondiSolver::advance<std::function<BondiSolver::WtBC(double)>>(
    double, const std::function<BondiSolver::WtBC(double)>&);

}  // namespace z4c_ccm
