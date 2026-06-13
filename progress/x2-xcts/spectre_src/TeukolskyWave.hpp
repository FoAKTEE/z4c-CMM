// Distributed under the MIT License.
// See LICENSE.txt for details.

// Original work for the z4c-CCM mission: even-parity (l=2, m=0) Teukolsky
// gravitational-wave free data for the SpECTRE XCTS elliptic solver. This
// provides the initial data for the arXiv:2308.10361 X=2 Teukolsky test
// (its eq:Teukolsky_ABC radial functions, eq:Teukolsky_angular_basis l=2 m=0
// angular tensors, and eq:Gaussian_pulse_F physicists'-Hermite Gaussian
// profile). The Cartesian metric perturbation h_ij(x) and its time derivative
// hdot_ij(x) at t=0 are the vacuum-verified iter-18 construction (see the
// reference implementations scripts/teuk_xcts_freedata_check.py and
// athenak/src/pgen/z4c_ccm_teukolsky.cpp in the z4c-CCM repository).

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <utility>

#include "DataStructures/CachedTempBuffer.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/CommonVariables.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "PointwiseFunctions/InitialDataUtilities/AnalyticSolution.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Xcts::AnalyticData {

namespace teukolsky_detail {

/*!
 * \brief Evaluate the even-parity (l=2, m=0) Teukolsky metric perturbation
 * \f$h_{ij}\f$ and its time derivative \f$\dot{h}_{ij}\f$ at a single Cartesian
 * point.
 *
 * Implements the construction verified in the z4c-CCM mission (see file-top
 * comment). The radial functions A, B, C are built from the regular
 * combinations \f$Q^{(n)}(t,r) = F^{(n)}(t-r) - (-1)^n F^{(n)}(t+r)\f$ and their
 * time derivatives \f$R^{(n)}(t,r) = F^{(n)}(t-r) + (-1)^n F^{(n)}(t+r)\f$ of
 * the Gaussian profile \f$F(u) = X e^{-((u-r_c)/\tau)^2}\f$ (eq:Gaussian_pulse_F),
 * assembled in the orthonormal spherical frame and rotated to Cartesian
 * components with the frame triad \f$(\hat r, \hat e_\theta, \hat e_\phi)\f$,
 * \f$\hat e_\theta = \hat e_\phi \times \hat r\f$.
 */
void teukolsky_metric_pointwise(double (*h)[3], double (*hdot)[3], double t,
                                double x, double y, double z, double amplitude,
                                double radial_coordinate, double timescale);

/// Evaluate only \f$h_{ij}\f$ at a single Cartesian point (used by the
/// finite-difference derivative).
void teukolsky_h_pointwise(double (*h)[3], double t, double x, double y,
                           double z, double amplitude, double radial_coordinate,
                           double timescale);

template <typename DataType>
using TeukolskyWaveVariablesCache =
    cached_temp_buffer_from_typelist<tmpl::push_back<
        Xcts::Solutions::common_tags<DataType>,
        ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::size_t<3>,
                      Frame::Inertial>,
        gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>, 0>,
        gr::Tags::Conformal<gr::Tags::StressTrace<DataType>, 0>,
        gr::Tags::Conformal<gr::Tags::MomentumDensity<DataType, 3>, 0>>>;

template <typename DataType>
struct TeukolskyWaveVariables
    : Xcts::Solutions::CommonVariables<DataType,
                                       TeukolskyWaveVariablesCache<DataType>> {
  static constexpr size_t Dim = 3;
  static constexpr int ConformalMatterScale = 0;
  using Cache = TeukolskyWaveVariablesCache<DataType>;
  using Base =
      Xcts::Solutions::CommonVariables<DataType,
                                       TeukolskyWaveVariablesCache<DataType>>;
  using Base::operator();

  const tnsr::I<DataType, 3>& x;
  double amplitude;
  double radial_coordinate;
  double timescale;

  TeukolskyWaveVariables(
      std::optional<std::reference_wrapper<const Mesh<Dim>>> local_mesh,
      std::optional<std::reference_wrapper<const InverseJacobian<
          DataType, Dim, Frame::ElementLogical, Frame::Inertial>>>
          local_inv_jacobian,
      const tnsr::I<DataType, 3>& local_x, double local_amplitude,
      double local_radial_coordinate, double local_timescale)
      : Base(std::move(local_mesh), std::move(local_inv_jacobian)),
        x(local_x),
        amplitude(local_amplitude),
        radial_coordinate(local_radial_coordinate),
        timescale(local_timescale) {}

  // Conformal metric gammabar_ij = delta_ij + h_ij (non-unimodular).
  void operator()(gsl::not_null<tnsr::ii<DataType, 3>*> conformal_metric,
                  gsl::not_null<Cache*> cache,
                  Tags::ConformalMetric<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::II<DataType, 3>*> inv_conformal_metric,
      gsl::not_null<Cache*> cache,
      Tags::InverseConformalMetric<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::ijj<DataType, 3>*> deriv_conformal_metric,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::ConformalMetric<DataType, 3, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;

  // Extrinsic curvature K_ij = -1/2 hdot_ij; trace K = 0 (maximal/TT).
  void operator()(
      gsl::not_null<tnsr::ii<DataType, 3>*> extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      gr::Tags::ExtrinsicCurvature<DataType, 3> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> trace_extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      gr::Tags::TraceExtrinsicCurvature<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*> deriv_trace_extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> dt_trace_extrinsic_curvature,
      gsl::not_null<Cache*> cache,
      ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataType>> /*meta*/)
      const override;

  // Conformal factor (flat initial guess).
  void operator()(gsl::not_null<Scalar<DataType>*> conformal_factor,
                  gsl::not_null<Cache*> cache,
                  Tags::ConformalFactor<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> conformal_factor_minus_one,
      gsl::not_null<Cache*> cache,
      Tags::ConformalFactorMinusOne<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*> deriv_conformal_factor,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Xcts::Tags::ConformalFactorMinusOne<DataType>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;

  // Lapse (flat initial guess).
  void operator()(gsl::not_null<Scalar<DataType>*> lapse,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Lapse<DataType> /*meta*/) const override;
  void operator()(gsl::not_null<tnsr::i<DataType, 3>*> deriv_lapse,
                  gsl::not_null<Cache*> cache,
                  ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::size_t<3>,
                                Frame::Inertial> /*meta*/) const;
  void operator()(
      gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
      gsl::not_null<Cache*> cache,
      Tags::LapseTimesConformalFactor<DataType> /*meta*/) const override;
  void operator()(
      gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor_minus_one,
      gsl::not_null<Cache*> cache,
      Tags::LapseTimesConformalFactorMinusOne<DataType> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::i<DataType, 3>*> deriv_lapse_times_conformal_factor,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::LapseTimesConformalFactorMinusOne<DataType>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;

  // Shift.
  void operator()(gsl::not_null<tnsr::I<DataType, 3>*> shift_background,
                  gsl::not_null<Cache*> cache,
                  Tags::ShiftBackground<DataType, 3, Frame::Inertial> /*meta*/)
      const override;
  void operator()(
      gsl::not_null<tnsr::iJ<DataType, 3, Frame::Inertial>*>
          deriv_shift_background,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Xcts::Tags::ShiftBackground<DataType, 3, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;
  void operator()(gsl::not_null<tnsr::II<DataType, 3, Frame::Inertial>*>
                      longitudinal_shift_background_minus_dt_conformal_metric,
                  gsl::not_null<Cache*> cache,
                  Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                      DataType, 3, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::I<DataType, 3>*> shift_excess,
      gsl::not_null<Cache*> cache,
      Tags::ShiftExcess<DataType, 3, Frame::Inertial> /*meta*/) const override;
  void operator()(
      gsl::not_null<tnsr::iJ<DataType, 3>*> deriv_shift_excess,
      gsl::not_null<Cache*> cache,
      ::Tags::deriv<Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
                    tmpl::size_t<3>, Frame::Inertial> /*meta*/) const override;

  // Matter sources (vacuum: all zero).
  void operator()(gsl::not_null<Scalar<DataType>*> energy_density,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>,
                                      ConformalMatterScale> /*meta*/) const;
  void operator()(gsl::not_null<Scalar<DataType>*> stress_trace,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<gr::Tags::StressTrace<DataType>,
                                      ConformalMatterScale> /*meta*/) const;
  void operator()(gsl::not_null<tnsr::I<DataType, 3>*> momentum_density,
                  gsl::not_null<Cache*> cache,
                  gr::Tags::Conformal<gr::Tags::MomentumDensity<DataType, 3>,
                                      ConformalMatterScale> /*meta*/) const;
};

}  // namespace teukolsky_detail

/*!
 * \brief Even-parity (l=2, m=0) Teukolsky gravitational-wave free data for the
 * XCTS equations.
 *
 * Specifies the free data for the X=2 Teukolsky test of arXiv:2308.10361
 * (Sec. V.A). The conformal metric is \f$\bar\gamma_{ij} = \delta_{ij} +
 * h_{ij}\f$ with the verified even-parity (l=2, m=0) Teukolsky perturbation
 * \f$h_{ij}\f$ (general, non-unimodular), the extrinsic curvature is
 * \f$K_{ij} = -\tfrac12 \dot{h}_{ij}\f$ with trace \f$K = 0\f$, and the
 * longitudinal-shift-minus-time-derivative background is the nonzero
 * wave-momentum term \f$-\dot{\bar\gamma}^{ij} = -\bar\gamma^{ik}\bar\gamma^{jl}
 * \dot{h}_{kl}\f$. The conformal factor, lapse, and shift are initialized flat;
 * the XCTS solver finds the constraint-satisfying initial data. This is a
 * vacuum problem: all conformal matter sources vanish.
 */
class TeukolskyWave : public elliptic::analytic_data::AnalyticSolution {
 public:
  struct Amplitude {
    using type = double;
    static constexpr Options::String help =
        "Amplitude X of the Teukolsky Gaussian profile";
    static double default_value() { return 2.0; }
  };
  struct RadialCoordinate {
    using type = double;
    static constexpr Options::String help =
        "Center r_c of the Teukolsky Gaussian pulse";
    static double default_value() { return 20.0; }
  };
  struct Timescale {
    using type = double;
    static constexpr Options::String help =
        "Width tau of the Teukolsky Gaussian pulse";
    static double default_value() { return 2.0; }
  };

  using options = tmpl::list<Amplitude, RadialCoordinate, Timescale>;
  static constexpr Options::String help =
      "Even-parity (l=2, m=0) Teukolsky gravitational-wave free data for the "
      "XCTS equations (arXiv:2308.10361 X=2 test).";

  TeukolskyWave() = default;
  TeukolskyWave(const TeukolskyWave&) = default;
  TeukolskyWave& operator=(const TeukolskyWave&) = default;
  TeukolskyWave(TeukolskyWave&&) = default;
  TeukolskyWave& operator=(TeukolskyWave&&) = default;
  ~TeukolskyWave() = default;

  TeukolskyWave(double amplitude, double radial_coordinate, double timescale);

  /// \cond
  explicit TeukolskyWave(CkMigrateMessage* m)
      : elliptic::analytic_data::AnalyticSolution(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(TeukolskyWave);
  std::unique_ptr<elliptic::analytic_data::AnalyticSolution> get_clone()
      const override {
    return std::make_unique<TeukolskyWave>(*this);
  }
  /// \endcond

  template <typename DataType>
  using tags =
      typename teukolsky_detail::TeukolskyWaveVariablesCache<DataType>::tags_list;

  template <typename DataType, typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataType, 3, Frame::Inertial>& x,
      tmpl::list<RequestedTags...> /*meta*/) const {
    return variables_impl<DataType>(x, std::nullopt, std::nullopt,
                                    tmpl::list<RequestedTags...>{});
  }

  template <typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataVector, 3, Frame::Inertial>& x, const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inv_jacobian,
      tmpl::list<RequestedTags...> /*meta*/) const {
    return variables_impl<DataVector>(x, mesh, inv_jacobian,
                                      tmpl::list<RequestedTags...>{});
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

  double amplitude() const { return amplitude_; }
  double radial_coordinate() const { return radial_coordinate_; }
  double timescale() const { return timescale_; }

 private:
  template <typename DataType, typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables_impl(
      const tnsr::I<DataType, 3, Frame::Inertial>& x,
      std::optional<std::reference_wrapper<const Mesh<3>>> mesh,
      std::optional<std::reference_wrapper<const InverseJacobian<
          DataType, 3, Frame::ElementLogical, Frame::Inertial>>>
          inv_jacobian,
      tmpl::list<RequestedTags...> /*meta*/) const {
    using VarsComputer = teukolsky_detail::TeukolskyWaveVariables<DataType>;
    typename VarsComputer::Cache cache{get_size(*x.begin())};
    const VarsComputer computer{std::move(mesh),     std::move(inv_jacobian),
                                x,                    amplitude_,
                                radial_coordinate_,   timescale_};
    return {cache.get_var(computer, RequestedTags{})...};
  }

  friend bool operator==(const TeukolskyWave& lhs, const TeukolskyWave& rhs);

  double amplitude_ = 2.0;
  double radial_coordinate_ = 20.0;
  double timescale_ = 2.0;
};

bool operator!=(const TeukolskyWave& lhs, const TeukolskyWave& rhs);

}  // namespace Xcts::AnalyticData
