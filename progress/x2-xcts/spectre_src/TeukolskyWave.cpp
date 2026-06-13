// Distributed under the MIT License.
// See LICENSE.txt for details.

// Original work for the z4c-CCM mission: even-parity (l=2, m=0) Teukolsky
// gravitational-wave free data for the SpECTRE XCTS elliptic solver, providing
// the initial data for the arXiv:2308.10361 X=2 Teukolsky test (its
// eq:Teukolsky_ABC radial functions, eq:Teukolsky_angular_basis l=2 m=0
// angular tensors, and eq:Gaussian_pulse_F physicists'-Hermite Gaussian
// profile). The pointwise h_ij / hdot_ij construction is ported from the
// vacuum-verified iter-18 reference (z4c-CCM repository:
// scripts/teuk_xcts_freedata_check.py and
// athenak/src/pgen/z4c_ccm_teukolsky.cpp).

#include "PointwiseFunctions/AnalyticData/Xcts/TeukolskyWave.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <pup.h>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/CommonVariables.tpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace {
// n-th derivative of the Gaussian profile F(u) = X exp(-((u-rc)/tau)^2) via
// physicists' Hermite polynomials: F^(n)(u) = X (-1/tau)^n H_n(s) e^{-s^2}.
// (arXiv:2308.10361 eq:Gaussian_pulse_F)
double teukolsky_f(const int n, const double u, const double amplitude,
                   const double rc, const double tau) {
  const double sarg = (u - rc) / tau;
  const double s2 = sarg * sarg;
  double hermite = 1.0;
  switch (n) {
    case 0:
      hermite = 1.0;
      break;
    case 1:
      hermite = 2.0 * sarg;
      break;
    case 2:
      hermite = 4.0 * s2 - 2.0;
      break;
    case 3:
      hermite = sarg * (8.0 * s2 - 12.0);
      break;
    case 4:
      hermite = (16.0 * s2 - 48.0) * s2 + 12.0;
      break;
    case 5:
      hermite = sarg * ((32.0 * s2 - 160.0) * s2 + 120.0);
      break;
    case 6:
      hermite = ((64.0 * s2 - 480.0) * s2 + 720.0) * s2 - 120.0;
      break;
    default:
      hermite = 0.0;
      break;
  }
  double pref = 1.0;
  for (int q = 0; q < n; ++q) {
    pref *= -1.0 / tau;
  }
  return amplitude * pref * hermite * std::exp(-s2);
}

// Regular combination Q^(n)(t,r) = F^(n)(t-r) - (-1)^n F^(n)(t+r).
double teukolsky_q(const int n, const double t, const double r,
                   const double amplitude, const double rc, const double tau) {
  const double sgn = (n % 2 == 0) ? 1.0 : -1.0;
  return teukolsky_f(n, t - r, amplitude, rc, tau) -
         sgn * teukolsky_f(n, t + r, amplitude, rc, tau);
}

// Time-derivative combination R^(n)(t,r) = F^(n)(t-r) + (-1)^n F^(n)(t+r),
// with d_t Q^(n) = R^(n+1).
double teukolsky_r(const int n, const double t, const double r,
                   const double amplitude, const double rc, const double tau) {
  const double sgn = (n % 2 == 0) ? 1.0 : -1.0;
  return teukolsky_f(n, t - r, amplitude, rc, tau) +
         sgn * teukolsky_f(n, t + r, amplitude, rc, tau);
}
}  // namespace

namespace Xcts::AnalyticData {

namespace teukolsky_detail {

void teukolsky_metric_pointwise(double (*const h)[3], double (*const hdot)[3],
                                const double t, const double x, const double y,
                                const double z, const double amplitude,
                                const double rc, const double tau) {
  const double r = std::max(std::sqrt(x * x + y * y + z * z), 1.0e-12);
  const double ir = 1.0 / r;

  const double q0 = teukolsky_q(0, t, r, amplitude, rc, tau);
  const double q1 = teukolsky_q(1, t, r, amplitude, rc, tau);
  const double q2 = teukolsky_q(2, t, r, amplitude, rc, tau);
  const double q3 = teukolsky_q(3, t, r, amplitude, rc, tau);
  const double q4 = teukolsky_q(4, t, r, amplitude, rc, tau);
  const double r1 = teukolsky_r(1, t, r, amplitude, rc, tau);
  const double r2 = teukolsky_r(2, t, r, amplitude, rc, tau);
  const double r3 = teukolsky_r(3, t, r, amplitude, rc, tau);
  const double r4 = teukolsky_r(4, t, r, amplitude, rc, tau);
  const double r5 = teukolsky_r(5, t, r, amplitude, rc, tau);

  const double ir2 = ir * ir;
  const double ir3 = ir2 * ir;
  const double ir4 = ir3 * ir;
  const double ir5 = ir4 * ir;

  // Radial functions (eq:Teukolsky_ABC) for h and (with the R combinations)
  // for hdot.
  const double a_f = 3.0 * (q2 * ir3 + 3.0 * q1 * ir4 + 3.0 * q0 * ir5);
  const double b_f =
      -(q3 * ir2 + 3.0 * q2 * ir3 + 6.0 * q1 * ir4 + 6.0 * q0 * ir5);
  const double c_f = 0.25 * (q4 * ir + 2.0 * q3 * ir2 + 9.0 * q2 * ir3 +
                             21.0 * q1 * ir4 + 21.0 * q0 * ir5);
  const double a_d = 3.0 * (r3 * ir3 + 3.0 * r2 * ir4 + 3.0 * r1 * ir5);
  const double b_d =
      -(r4 * ir2 + 3.0 * r3 * ir3 + 6.0 * r2 * ir4 + 6.0 * r1 * ir5);
  const double c_d = 0.25 * (r5 * ir + 2.0 * r4 * ir2 + 9.0 * r3 * ir3 +
                             21.0 * r2 * ir4 + 21.0 * r1 * ir5);

  // Angular functions (l=2, m=0 even parity; eq:Teukolsky_angular_basis).
  const double cth = z * ir;
  const double s2 = std::max(1.0 - cth * cth, 0.0);  // sin^2 theta
  const double sth = std::sqrt(s2);
  const double frr = 0.5 * (1.0 + 3.0 * (1.0 - 2.0 * s2));  // (1+3cos2th)/2
  const double frt = -3.0 * sth * cth;
  const double f1tt = 3.0 * s2;
  const double f2tt = -1.0;
  const double f1pp = -3.0 * s2;
  const double f2pp = 3.0 * s2 - 1.0;

  // Orthonormal spherical-frame components.
  const double hrr_o = a_f * frr;
  const double hrr_d = a_d * frr;
  const double hrt_o = b_f * frt;
  const double hrt_d = b_d * frt;
  const double htt_o = c_f * f1tt + a_f * f2tt;
  const double htt_d = c_d * f1tt + a_d * f2tt;
  const double hpp_o = c_f * f1pp + a_f * f2pp;
  const double hpp_d = c_d * f1pp + a_d * f2pp;

  // Frame triad (Cartesian); axis guard: at sin(theta) -> 0 pick e_phi = y-hat.
  // e_theta = e_phi x r-hat = +theta-hat (the iter-18-corrected order).
  const std::array<double, 3> rh{{x * ir, y * ir, z * ir}};
  std::array<double, 3> eph{};
  const double rho = std::sqrt(std::max(x * x + y * y, 1.0e-300));
  if (rho > 1.0e-12 * r) {
    eph = {{-y / rho, x / rho, 0.0}};
  } else {
    eph = {{0.0, 1.0, 0.0}};
  }
  const std::array<double, 3> eth{
      {eph[1] * rh[2] - eph[2] * rh[1], eph[2] * rh[0] - eph[0] * rh[2],
       eph[0] * rh[1] - eph[1] * rh[0]}};

  for (size_t a = 0; a < 3; ++a) {
    for (size_t b = 0; b < 3; ++b) {
      h[a][b] = hrr_o * rh[a] * rh[b] +
                hrt_o * (rh[a] * eth[b] + eth[a] * rh[b]) +
                htt_o * eth[a] * eth[b] + hpp_o * eph[a] * eph[b];
      hdot[a][b] = hrr_d * rh[a] * rh[b] +
                   hrt_d * (rh[a] * eth[b] + eth[a] * rh[b]) +
                   htt_d * eth[a] * eth[b] + hpp_d * eph[a] * eph[b];
    }
  }
}

void teukolsky_h_pointwise(double (*const h)[3], const double t, const double x,
                           const double y, const double z,
                           const double amplitude, const double rc,
                           const double tau) {
  double hdot_unused[3][3];
  teukolsky_metric_pointwise(h, &hdot_unused[0], t, x, y, z, amplitude, rc, tau);
}

namespace {
// Compute h_ij at grid point `i` of the (DataType) coordinates.
template <typename DataType>
void h_at_point(double (&h)[3][3], const tnsr::I<DataType, 3>& x,
                const size_t i, const double amplitude, const double rc,
                const double tau) {
  teukolsky_h_pointwise(&h[0], 0., get_element(get<0>(x), i),
                        get_element(get<1>(x), i), get_element(get<2>(x), i),
                        amplitude, rc, tau);
}

// Compute h_ij and hdot_ij at grid point `i` of the (DataType) coordinates.
template <typename DataType>
void h_hdot_at_point(double (&h)[3][3], double (&hdot)[3][3],
                     const tnsr::I<DataType, 3>& x, const size_t i,
                     const double amplitude, const double rc,
                     const double tau) {
  teukolsky_metric_pointwise(&h[0], &hdot[0], 0., get_element(get<0>(x), i),
                             get_element(get<1>(x), i),
                             get_element(get<2>(x), i), amplitude, rc, tau);
}
}  // namespace

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ii<DataType, 3>*> conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::ConformalMetric<DataType, 3, Frame::Inertial> /*meta*/) const {
  // gammabar_ij = delta_ij + h_ij (general, non-unimodular).
  for (size_t i = 0; i < get_size(get<0>(x)); ++i) {
    double h[3][3];
    h_at_point(h, x, i, amplitude, radial_coordinate, timescale);
    for (size_t a = 0; a < 3; ++a) {
      for (size_t b = a; b < 3; ++b) {
        get_element(conformal_metric->get(a, b), i) =
            (a == b ? 1.0 : 0.0) + h[a][b];
      }
    }
  }
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::II<DataType, 3>*> inv_conformal_metric,
    const gsl::not_null<Cache*> cache,
    Tags::InverseConformalMetric<DataType, 3, Frame::Inertial> /*meta*/) const {
  // Pointwise direct cofactor inverse of gammabar_ij (no delta-h assumption).
  const auto& conformal_metric = cache->get_var(
      *this, Tags::ConformalMetric<DataType, 3, Frame::Inertial>{});
  for (size_t i = 0; i < get_size(get<0>(x)); ++i) {
    const double g00 = get_element(get<0, 0>(conformal_metric), i);
    const double g01 = get_element(get<0, 1>(conformal_metric), i);
    const double g02 = get_element(get<0, 2>(conformal_metric), i);
    const double g11 = get_element(get<1, 1>(conformal_metric), i);
    const double g12 = get_element(get<1, 2>(conformal_metric), i);
    const double g22 = get_element(get<2, 2>(conformal_metric), i);
    const double c00 = g11 * g22 - g12 * g12;
    const double c01 = g02 * g12 - g01 * g22;
    const double c02 = g01 * g12 - g02 * g11;
    const double c11 = g00 * g22 - g02 * g02;
    const double c12 = g02 * g01 - g00 * g12;
    const double c22 = g00 * g11 - g01 * g01;
    const double det = g00 * c00 + g01 * c01 + g02 * c02;
    const double inv_det = 1.0 / det;
    get_element(get<0, 0>(*inv_conformal_metric), i) = c00 * inv_det;
    get_element(get<0, 1>(*inv_conformal_metric), i) = c01 * inv_det;
    get_element(get<0, 2>(*inv_conformal_metric), i) = c02 * inv_det;
    get_element(get<1, 1>(*inv_conformal_metric), i) = c11 * inv_det;
    get_element(get<1, 2>(*inv_conformal_metric), i) = c12 * inv_det;
    get_element(get<2, 2>(*inv_conformal_metric), i) = c22 * inv_det;
  }
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ijj<DataType, 3>*> deriv_conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Tags::ConformalMetric<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  // d_k gammabar_ij = d_k h_ij, by 4th-order centered finite difference of the
  // analytic h_ij in Cartesian x, y, z (smooth field; avoids frame-vector axis
  // singularities of an analytic derivative).
  constexpr double eps = 1.0e-3;
  const double inv_12eps = 1.0 / (12.0 * eps);
  for (size_t pt = 0; pt < get_size(get<0>(x)); ++pt) {
    const double xc = get_element(get<0>(x), pt);
    const double yc = get_element(get<1>(x), pt);
    const double zc = get_element(get<2>(x), pt);
    for (size_t k = 0; k < 3; ++k) {
      const double dx = (k == 0) ? eps : 0.0;
      const double dy = (k == 1) ? eps : 0.0;
      const double dz = (k == 2) ? eps : 0.0;
      double hp1[3][3];
      double hm1[3][3];
      double hp2[3][3];
      double hm2[3][3];
      teukolsky_h_pointwise(&hp1[0], 0., xc + dx, yc + dy, zc + dz, amplitude,
                            radial_coordinate, timescale);
      teukolsky_h_pointwise(&hm1[0], 0., xc - dx, yc - dy, zc - dz, amplitude,
                            radial_coordinate, timescale);
      teukolsky_h_pointwise(&hp2[0], 0., xc + 2.0 * dx, yc + 2.0 * dy,
                            zc + 2.0 * dz, amplitude, radial_coordinate,
                            timescale);
      teukolsky_h_pointwise(&hm2[0], 0., xc - 2.0 * dx, yc - 2.0 * dy,
                            zc - 2.0 * dz, amplitude, radial_coordinate,
                            timescale);
      for (size_t a = 0; a < 3; ++a) {
        for (size_t b = a; b < 3; ++b) {
          get_element(deriv_conformal_metric->get(k, a, b), pt) =
              (-hp2[a][b] + 8.0 * hp1[a][b] - 8.0 * hm1[a][b] + hm2[a][b]) *
              inv_12eps;
        }
      }
    }
  }
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ii<DataType, 3>*> extrinsic_curvature,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::ExtrinsicCurvature<DataType, 3> /*meta*/) const {
  // K_ij = -1/2 hdot_ij (alpha = 1, beta = 0 at this order).
  for (size_t i = 0; i < get_size(get<0>(x)); ++i) {
    double h[3][3];
    double hdot[3][3];
    h_hdot_at_point(h, hdot, x, i, amplitude, radial_coordinate, timescale);
    for (size_t a = 0; a < 3; ++a) {
      for (size_t b = a; b < 3; ++b) {
        get_element(extrinsic_curvature->get(a, b), i) = -0.5 * hdot[a][b];
      }
    }
  }
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> trace_extrinsic_curvature,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::TraceExtrinsicCurvature<DataType> /*meta*/) const {
  get(*trace_extrinsic_curvature) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*> deriv_trace_extrinsic_curvature,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>, tmpl::size_t<3>,
                  Frame::Inertial> /*meta*/) const {
  std::fill(deriv_trace_extrinsic_curvature->begin(),
            deriv_trace_extrinsic_curvature->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> dt_trace_extrinsic_curvature,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataType>> /*meta*/) const {
  get(*dt_trace_extrinsic_curvature) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> conformal_factor,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::ConformalFactor<DataType> /*meta*/) const {
  get(*conformal_factor) = 1.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> conformal_factor_minus_one,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::ConformalFactorMinusOne<DataType> /*meta*/) const {
  get(*conformal_factor_minus_one) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*> deriv_conformal_factor,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Tags::ConformalFactorMinusOne<DataType>, tmpl::size_t<3>,
                  Frame::Inertial> /*meta*/) const {
  std::fill(deriv_conformal_factor->begin(), deriv_conformal_factor->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> lapse,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Lapse<DataType> /*meta*/) const {
  get(*lapse) = 1.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*> deriv_lapse,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::size_t<3>,
                  Frame::Inertial> /*meta*/) const {
  std::fill(deriv_lapse->begin(), deriv_lapse->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::LapseTimesConformalFactor<DataType> /*meta*/) const {
  get(*lapse_times_conformal_factor) = 1.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*>
        lapse_times_conformal_factor_minus_one,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::LapseTimesConformalFactorMinusOne<DataType> /*meta*/) const {
  get(*lapse_times_conformal_factor_minus_one) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*>
        deriv_lapse_times_conformal_factor,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  std::fill(deriv_lapse_times_conformal_factor->begin(),
            deriv_lapse_times_conformal_factor->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> shift_background,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::ShiftBackground<DataType, 3, Frame::Inertial> /*meta*/) const {
  std::fill(shift_background->begin(), shift_background->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::iJ<DataType, 3, Frame::Inertial>*>
        deriv_shift_background,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Xcts::Tags::ShiftBackground<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  std::fill(deriv_shift_background->begin(), deriv_shift_background->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::II<DataType, 3, Frame::Inertial>*>
        longitudinal_shift_background_minus_dt_conformal_metric,
    const gsl::not_null<Cache*> cache,
    Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
        DataType, 3, Frame::Inertial> /*meta*/) const {
  // (L_bar beta_background)^ij - ubar^ij with beta_background = 0 and
  // ubar^ij = dt(gammabar)^ij = hdot^ij raised. The wave-momentum term:
  //   result^ij = -(dt gammabar)^ij = - gammabar^ik gammabar^jl hdot_kl.
  const auto& inv_conformal_metric = cache->get_var(
      *this, Tags::InverseConformalMetric<DataType, 3, Frame::Inertial>{});
  for (size_t pt = 0; pt < get_size(get<0>(x)); ++pt) {
    double h[3][3];
    double hdot[3][3];
    h_hdot_at_point(h, hdot, x, pt, amplitude, radial_coordinate, timescale);
    double inv[3][3];
    for (size_t a = 0; a < 3; ++a) {
      for (size_t b = 0; b < 3; ++b) {
        inv[a][b] =
            get_element(inv_conformal_metric.get(std::min(a, b), std::max(a, b)),
                        pt);
      }
    }
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = i; j < 3; ++j) {
        double sum = 0.;
        for (size_t k = 0; k < 3; ++k) {
          for (size_t l = 0; l < 3; ++l) {
            sum += inv[i][k] * inv[j][l] * hdot[k][l];
          }
        }
        get_element(
            longitudinal_shift_background_minus_dt_conformal_metric->get(i, j),
            pt) = -sum;
      }
    }
  }
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> shift_excess,
    const gsl::not_null<Cache*> /*cache*/,
    Tags::ShiftExcess<DataType, 3, Frame::Inertial> /*meta*/) const {
  std::fill(shift_excess->begin(), shift_excess->end(), 0.);
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::iJ<DataType, 3>*> deriv_shift_excess,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  std::fill(deriv_shift_excess->begin(), deriv_shift_excess->end(), 0.);
}

// Matter sources: vacuum, all zero.
template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> energy_density,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>,
                        ConformalMatterScale> /*meta*/) const {
  get(*energy_density) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> stress_trace,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::StressTrace<DataType>,
                        ConformalMatterScale> /*meta*/) const {
  get(*stress_trace) = 0.;
}

template <typename DataType>
void TeukolskyWaveVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> momentum_density,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::MomentumDensity<DataType, 3>,
                        ConformalMatterScale> /*meta*/) const {
  std::fill(momentum_density->begin(), momentum_density->end(), 0.);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define INSTANTIATE(_, data) template class TeukolskyWaveVariables<DTYPE(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector))

#undef DTYPE
#undef INSTANTIATE

}  // namespace teukolsky_detail

// TeukolskyWave class implementation

TeukolskyWave::TeukolskyWave(const double amplitude,
                             const double radial_coordinate,
                             const double timescale)
    : amplitude_(amplitude),
      radial_coordinate_(radial_coordinate),
      timescale_(timescale) {}

void TeukolskyWave::pup(PUP::er& p) {
  elliptic::analytic_data::AnalyticSolution::pup(p);
  p | amplitude_;
  p | radial_coordinate_;
  p | timescale_;
}

bool operator==(const TeukolskyWave& lhs, const TeukolskyWave& rhs) {
  return lhs.amplitude_ == rhs.amplitude_ &&
         lhs.radial_coordinate_ == rhs.radial_coordinate_ &&
         lhs.timescale_ == rhs.timescale_;
}

bool operator!=(const TeukolskyWave& lhs, const TeukolskyWave& rhs) {
  return not(lhs == rhs);
}

PUP::able::PUP_ID TeukolskyWave::my_PUP_ID = 0;  // NOLINT

}  // namespace Xcts::AnalyticData

// Instantiate CommonVariables in the correct namespaces
namespace Xcts::Solutions {
template class CommonVariables<
    double, typename Xcts::AnalyticData::teukolsky_detail::
                TeukolskyWaveVariablesCache<double>>;
template class CommonVariables<
    DataVector, typename Xcts::AnalyticData::teukolsky_detail::
                    TeukolskyWaveVariablesCache<DataVector>>;
}  // namespace Xcts::Solutions

namespace Xcts::AnalyticData {
template class CommonVariables<
    double, typename teukolsky_detail::TeukolskyWaveVariablesCache<double>>;
template class CommonVariables<
    DataVector,
    typename teukolsky_detail::TeukolskyWaveVariablesCache<DataVector>>;
}  // namespace Xcts::AnalyticData
