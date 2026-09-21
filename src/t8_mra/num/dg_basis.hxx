#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <type_traits>

#include "t8_eclass/t8_eclass.h"

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/quadrature/quadrature.hxx"

namespace t8_mra
{

/// Primary template left undefined: only cartesian shapes
template <t8_eclass TShape, typename = void>
struct dg_basis_base;

template <t8_eclass TShape>
struct dg_basis_base<TShape, std::enable_if_t<is_cartesian<TShape>>>
{
  static constexpr unsigned int DIM = TShape == T8_ECLASS_LINE ? 1 : (TShape == T8_ECLASS_QUAD ? 2 : 3);
  static constexpr t8_eclass Shape = TShape;

  quadrature<TShape> quad;

  dg_basis_base () = default;

  explicit dg_basis_base (int num_quad_points_1d): quad (num_quad_points_1d)
  {
  }
};

/// Simplex-based shapes take an accuracy degree rather than a point count.
template <t8_eclass TShape>
struct dg_basis_base<TShape, std::enable_if_t<TShape == T8_ECLASS_TRIANGLE || TShape == T8_ECLASS_PRISM>>
{
  static constexpr unsigned int DIM = shape_traits<TShape>::DIM;
  static constexpr t8_eclass Shape = TShape;

  quadrature<TShape> quad;

  dg_basis_base () = default;

  explicit dg_basis_base (int quadrature_rule): quad (quadrature_rule)
  {
  }
};

template <typename TElement>
class dg_basis: public dg_basis_base<TElement::Shape> {
  using Element = TElement;
  using Base = dg_basis_base<TElement::Shape>;

  static constexpr unsigned int DIM = Element::DIM;
  static constexpr auto Shape = TElement::Shape;

  static constexpr unsigned int P_DIM = Element::P_DIM;
  static constexpr unsigned int DOF = Element::DOF;
  static constexpr unsigned int W_DOF = Element::W_DOF;

  using basis_t = basis<Shape, P_DIM>;

 public:
  dg_basis () = default;

  /// Constructor for simplex-based elements (TRIANGLE, PRISM)
  explicit dg_basis (int _quadrature_rule)
    requires (Shape == T8_ECLASS_TRIANGLE || Shape == T8_ECLASS_PRISM)
    : Base (_quadrature_rule)
  {
  }

  /// Constructor for cartesian elements (LINE, QUAD, HEX)
  explicit dg_basis (int _num_quad_points_1d)
    requires is_cartesian<Shape>
    : Base (_num_quad_points_1d)
  {
  }

  /// All basis function values at a reference point.
  [[nodiscard]] std::array<double, DOF>
  basis_value (const std::array<double, DIM> &x_ref) const
  {
    return basis_t::eval (x_ref);
  }

  /// grad[dir][i] = d(phi_i)/dx_dir at a reference point.
  [[nodiscard]] std::array<std::array<double, DOF>, DIM>
  basis_gradient (const std::array<double, DIM> &x_ref) const
  {
    return basis_t::eval_gradient (x_ref);
  }
};

}  // namespace t8_mra
#endif
