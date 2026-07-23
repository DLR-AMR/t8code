#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <concepts>
#include <span>
#include <utility>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"

namespace t8_mra
{

// ============================================================================
// Reference-element polynomial basis
// ============================================================================
// One specialization of basis<TShape, P> per element shape is the single place
// that defines a shape's function space (in num/shape/). Everything that
// evaluates the basis (projection in dg_basis, VTK output) goes through this
// interface.

/// Minimal interface every basis specialization provides.
template <typename B>
concept reference_basis = requires (std::array<double, static_cast<std::size_t> (B::DIM)> x, double vol) {
  { B::DIM } -> std::convertible_to<int>;
  { B::DOF } -> std::convertible_to<int>;
  { B::eval (x) } -> std::same_as<std::array<double, static_cast<std::size_t> (B::DOF)>>;
  { B::normalization (vol) } -> std::convertible_to<double>;
};

template <t8_eclass TShape, int P>
struct basis;

/// Physical cell mean of a modal field (only the zeroth mode survives).
template <t8_eclass Shape, int P>
inline double
cell_mean (std::span<const double> coeffs, double vol)
{
  using basis_t = basis<Shape, P>;
  static const double phi0 = basis_t::eval ({})[0];
  return basis_t::normalization (vol) * phi0 * coeffs[0];
}

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/num/shape/cartesian.hxx"
#include "t8_mra/num/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA
