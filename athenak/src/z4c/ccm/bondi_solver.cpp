//========================================================================================
// bondi_solver.cpp — see bondi_solver.hpp. Port gated against the
// ZccmJl-admitted numbers (N14 ledger): exact-slice sweep machine-exact;
// SAT evolution ~4.3-order; Float64 psi0 conditioning floor ~0.05 of the
// r^-5-tiny linear datum (physics, not scheme — iter-35 analysis).
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
  // Aint = inv(D with row0 -> e0) * (I with row0 zeroed)
  std::vector<double> A(D_);
  for (int j = 0; j < n; ++j) A[0*n + j] = 0.0;
  A[0] = 1.0;
  std::vector<int> piv;
  lu_factor(&A, &piv, n);
  Aint_.assign(n*n, 0.0);
  for (int col = 0; col < n; ++col) {
    std::vector<double> e(n, 0.0);
    if (col != 0) e[col] = 1.0;
    lu_solve(A, piv, n, &e);
    for (int i = 0; i < n; ++i) Aint_[i*n + col] = e[i];
  }
  // sweep operators (verified linear system; BC row 1 each; W: scri row by
  // L'Hopital handled via the same structure as ZccmJl cheb_sweep)
  auto build = [&](std::vector<double>* M, int kind) {
    M->assign(n*n, 0.0);
    for (int i = 0; i < n; ++i) {
      const double oy = 1.0 - y_[i];
      for (int j = 0; j < n; ++j) {
        double v = oy*D_[i*n + j];
        if (kind == 0) v += (i == j) ? 1.0 : 0.0;        // Q: I + (1-y)D
        if (kind == 1) v += (i == j) ? -2.0 : 0.0;       // U: -2I + (1-y)D
        (*M)[i*n + j] = v;                                // W: (1-y)D
      }
    }
    for (int j = 0; j < n; ++j) (*M)[0*n + j] = 0.0;
    (*M)[0] = 1.0;                                        // BC row
    if (kind == 2) {                                      // W scri row
      for (int j = 0; j < n; ++j) (*M)[(n-1)*n + j] = -D_[(n-1)*n + j];
    }
  };
  build(&luQ_, 0); lu_factor(&luQ_, &pivQ_, n);
  build(&luU_, 1); lu_factor(&luU_, &pivU_, n);
  build(&luW_, 2); lu_factor(&luW_, &pivW_, n);
  taupen_ = (1.0/(2.0*rwt))*N*N;       // SAT: sigma=1, v_bdry * N^2
  jr_.assign(n, 0.0);
}

void BondiSolver::init_teukolsky(double X, double rc, double tau) {
  for (int i = 0; i < n_; ++i) {
    const double ir = (1.0 - y_[i])/(2.0*rwt_);          // 1/r (0 at scri)
    jr_[i] = 0.75*(Fn(4, u_, X, rc, tau)
                   + Fn(2, u_, X, rc, tau)*ir*ir);
  }
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
  return b;
}

void BondiSolver::sweep(const std::vector<double>& jr, const WtBC& bc,
                        std::vector<double>* qr, std::vector<double>* ur,
                        std::vector<double>* wr, std::vector<double>* hr) {
  const int n = n_;
  std::vector<double> djr;
  matvec(D_, jr, n, &djr);
  // Q
  std::vector<double> b(n);
  for (int i = 0; i < n; ++i) {
    b[i] = 4.0*jr[i] - 4.0*(1.0 - y_[i])*djr[i];
  }
  b[0] = bc.qr;
  lu_solve(luQ_, pivQ_, n, &b);
  *qr = b;
  // U
  b = *qr; b[0] = bc.ur;
  lu_solve(luU_, pivU_, n, &b);
  *ur = b;
  // W
  for (int i = 0; i < n; ++i) b[i] = 2.0*jr[i] + 2.0*(*ur)[i] + 0.5*(*qr)[i];
  std::vector<double> srcp;
  matvec(D_, b, n, &srcp);
  b[0] = bc.wr; b[n-1] = srcp[n-1];
  lu_solve(luW_, pivW_, n, &b);
  *wr = b;
  // H via integration by parts (first derivatives only):
  // hr = hr_wt + (G - G_1)/4 - Aint(ur)/(4 rwt) - 3 Aint(jr)/(4 rwt)
  std::vector<double> G(n), iu, ij;
  for (int i = 0; i < n; ++i) {
    const double oy = 1.0 - y_[i];
    G[i] = oy*oy/(2.0*rwt_)*djr[i] - oy*jr[i]/rwt_;
  }
  matvec(Aint_, *ur, n, &iu);
  matvec(Aint_, jr, n, &ij);
  hr->resize(n);
  for (int i = 0; i < n; ++i) {
    (*hr)[i] = bc.hr + (G[i] - G[0])/4.0 - iu[i]/(4.0*rwt_)
               - 3.0*ij[i]/(4.0*rwt_);
  }
}

std::vector<double> BondiSolver::rhs(double u, const std::vector<double>& jr,
                                     const WtBC& bc) {
  std::vector<double> qr, ur, wr, hr;
  sweep(jr, bc, &qr, &ur, &wr, &hr);
  hr[0] += taupen_*(bc.jr - jr[0]);                       // SAT (sigma = 1)
  return hr;
}

template <typename F>
void BondiSolver::advance(double t, const F& bc) {
  while (u_ < t - 1e-12) {
    const double du = std::fmin(du_max_, t - u_);
    const WtBC b0 = bc(u_), bh = bc(u_ + du/2.0), b1 = bc(u_ + du);
    std::vector<double> k1 = rhs(u_, jr_, b0);
    std::vector<double> tmp(n_);
    for (int i = 0; i < n_; ++i) tmp[i] = jr_[i] + du/2.0*k1[i];
    std::vector<double> k2 = rhs(u_ + du/2.0, tmp, bh);
    for (int i = 0; i < n_; ++i) tmp[i] = jr_[i] + du/2.0*k2[i];
    std::vector<double> k3 = rhs(u_ + du/2.0, tmp, bh);
    for (int i = 0; i < n_; ++i) tmp[i] = jr_[i] + du*k3[i];
    std::vector<double> k4 = rhs(u_ + du, tmp, b1);
    for (int i = 0; i < n_; ++i) {
      jr_[i] += du/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    u_ += du;
  }
}
// explicit instantiation for the std::function closure used by callers
template void BondiSolver::advance<std::function<BondiSolver::WtBC(double)>>(
    double, const std::function<BondiSolver::WtBC(double)>&);

double BondiSolver::psi0_worldtube(const WtBC& bcnow) {
  // hierarchy-based: dy jr|_wt from the Q equation; d2r J from the H
  // identity (ZccmJl cheb_psi0_hierarchy, verified iters 35-37)
  const int n = n_;
  std::vector<double> qr, ur, wr, hr;
  sweep(jr_, bcnow, &qr, &ur, &wr, &hr);
  std::vector<double> dq, dh;
  matvec(D_, qr, n, &dq);
  matvec(D_, hr, n, &dh);
  const double oy = 1.0 - y_[0];                          // = 2
  const double w = oy*oy/(2.0*rwt_);
  const double djr1 = (4.0*jr_[0] - qr[0] - oy*dq[0])/(4.0*oy);
  const double dJt = w*djr1/rwt_ - jr_[0]/(rwt_*rwt_);
  const double drH = w*dh[0];
  const double d2Jt = (drH + ur[0]/(rwt_*rwt_)/2.0
                       + 1.5*(jr_[0]/rwt_)/rwt_)*4.0/rwt_;
  return -dJt/(2.0*rwt_) - d2Jt/4.0;
}

}  // namespace z4c_ccm
