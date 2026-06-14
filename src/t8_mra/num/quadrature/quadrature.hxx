#pragma once

#ifdef T8_ENABLE_MRA

#include <cstddef>
#include <vector>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/quadrature/gauss_legendre.hxx"
#include "t8_mra/num/quadrature/dunavant.hxx"

namespace t8_mra
{

// ============================================================================
// Reference-element quadrature rule
// ============================================================================
// One specialization of quadrature<TShape> per shape provides the integration
// points and weights on the reference element, behind a common interface.
// Cartesian shapes use a tensor product of 1D Gauss-Legendre; the triangle
// uses a Dunavant rule. Mirrors basis<TShape, P> (num/basis/basis.hxx).
//
// Adding a new shape: specialize quadrature<NewShape> with points/weights and
// a constructor taking the rule's accuracy parameter.

/// Common interface every quadrature specialization provides: a flat list of
/// num_points reference points (DIM coords each, point q at points[DIM*q + d])
/// and matching weights.
template <typename Q>
concept quadrature_rule = requires (const Q q) {
  { Q::DIM } -> std::convertible_to<int>;
  { q.num_points } -> std::convertible_to<std::size_t>;
  { q.points.data () } -> std::convertible_to<const double *>;
  { q.weights.data () } -> std::convertible_to<const double *>;
};

template <t8_eclass TShape, typename = void>
struct quadrature;

/// Cartesian shapes: tensor product of a 1D Gauss-Legendre rule, exact to
/// degree 2*num_points_1d - 1 per axis.
template <t8_eclass TShape>
struct quadrature<TShape, std::enable_if_t<is_cartesian<TShape>>>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;

  std::size_t num_points = 0;
  std::vector<double> points;   // flattened: point q coord d at points[DIM*q + d]
  std::vector<double> weights;

  quadrature () = default;

  explicit quadrature (int num_points_1d)
  {
    std::vector<double> p1d, w1d;
    gauss_legendre_1d (num_points_1d, p1d, w1d);

    num_points = 1;
    for (int d = 0; d < DIM; ++d)
      num_points *= num_points_1d;

    points.resize (DIM * num_points);
    weights.resize (num_points);

    // Odometer over the DIM axes (first axis fastest); order is irrelevant to
    // the integration sum, so any consistent enumeration works.
    for (std::size_t q = 0; q < num_points; ++q) {
      std::size_t rest = q;
      double w = 1.0;
      for (int d = 0; d < DIM; ++d) {
        const int id = rest % num_points_1d;
        rest /= num_points_1d;
        points[DIM * q + d] = p1d[id];
        w *= w1d[id];
      }
      weights[q] = w;
    }
  }
};

/// Triangle: a Dunavant rule on the reference triangle.
template <>
struct quadrature<T8_ECLASS_TRIANGLE>
{
  static constexpr int DIM = 2;

  std::size_t num_points = 0;
  std::vector<double> points;   // flattened: [x0, y0, x1, y1, ...]
  std::vector<double> weights;

  quadrature () = default;

  explicit quadrature (int dunavant_rule_num)
    : num_points (dunavant_order_num (dunavant_rule_num)), points (2 * num_points, 0.0), weights (num_points, 0.0)
  {
    dunavant_rule (dunavant_rule_num, num_points, points.data (), weights.data ());
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
