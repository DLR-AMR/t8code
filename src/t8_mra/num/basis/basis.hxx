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
template <typename TBasis>
concept reference_basis = requires (std::array<double, static_cast<std::size_t> (TBasis::DIM)> x, double vol) {
  { TBasis::DIM } -> std::convertible_to<int>;
  { TBasis::DOF } -> std::convertible_to<int>;
  { TBasis::eval (x) } -> std::same_as<std::array<double, static_cast<std::size_t> (TBasis::DOF)>>;
  { TBasis::normalization (vol) } -> std::convertible_to<double>;
};

template <t8_eclass TShape, int P>
struct basis;

/// The constant mode on the reference cell, the only mode a cell mean sees.
template <t8_eclass TShape, int P>
inline const double reference_mode0 = basis<TShape, P>::eval ({})[0];

/// Physical cell mean of a modal field (only the zeroth mode survives).
template <t8_eclass TShape, int P>
[[nodiscard]] inline double
cell_mean (std::span<const double> coeffs, double vol)
{
  return basis<TShape, P>::normalization (vol) * reference_mode0<TShape, P> * coeffs[0];
}

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/num/shape/cartesian.hxx"
#include "t8_mra/num/shape/prism.hxx"
#include "t8_mra/num/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA
