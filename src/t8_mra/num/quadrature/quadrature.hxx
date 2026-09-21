#pragma once

#ifdef T8_ENABLE_MRA

#include <concepts>
#include <cstddef>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"

namespace t8_mra
{

// ============================================================================
// Reference-element quadrature rule
// ============================================================================
// One specialization of quadrature<TShape> per shape provides the integration
// points and weights on the reference element, behind a common interface.
// Cartesian shapes use a tensor product of 1D Gauss-Legendre, the triangle a
// Dunavant rule and the prism their product. Mirrors basis<TShape, P>.
//
// Adding a new shape: specialize quadrature<NewShape> in num/shape/ with
// points/weights and a constructor taking the rule's accuracy parameter.

/// Common interface every quadrature specialization provides: a flat list of
/// num_points reference points (DIM coords each, point q at points[DIM*q + d])
/// and matching weights.
template <typename TQuadrature>
concept quadrature_rule = requires (const TQuadrature q) {
  { TQuadrature::DIM } -> std::convertible_to<int>;
  { q.num_points } -> std::convertible_to<std::size_t>;
  { q.points.data () } -> std::convertible_to<const double *>;
  { q.weights.data () } -> std::convertible_to<const double *>;
};

template <t8_eclass TShape, typename = void>
struct quadrature;

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/num/shape/cartesian.hxx"
#include "t8_mra/num/shape/triangle.hxx"
#endif  // T8_ENABLE_MRA
