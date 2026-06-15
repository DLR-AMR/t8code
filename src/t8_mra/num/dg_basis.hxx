#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <array>
#include <type_traits>

#include "t8_eclass/t8_eclass.h"
#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/quadrature/quadrature.hxx"
#include "t8_mra/num/mat.hxx"
#include "t8_mra/num/geometry.hxx"

namespace t8_mra
{

/// Primary template left undefined: only cartesian shapes and the triangle are
/// supported (specializations below). Any other shape is a compile error.
template <t8_eclass TShape, typename = void>
struct dg_basis_base;

template <t8_eclass T>
struct dg_basis_base<T, std::enable_if_t<is_cartesian<T>>>
{
  static constexpr unsigned int DIM = T == T8_ECLASS_LINE ? 1 : (T == T8_ECLASS_QUAD ? 2 : 3);
  static constexpr t8_eclass Shape = T;

  // Reference Gauss-Legendre tensor rule.
  quadrature<T> quad;

  dg_basis_base () = default;

  explicit dg_basis_base (int num_quad_points_1d): quad (num_quad_points_1d)
  {
  }

  /// Maps the reference quadrature points to the physical element (flattened).
  std::vector<double>
  deref_quad_points (const double physical_vertices[][3])
  {
    std::array<double, DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<DIM> (physical_vertices, vertices_min, vertices_max);

    return transform_quad_points<DIM> (quad.points, quad.num_points, vertices_min, vertices_max);
  }
};

template <>
struct dg_basis_base<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned int DIM = 2;
  static constexpr t8_eclass Shape = T8_ECLASS_TRIANGLE;

  quadrature<T8_ECLASS_TRIANGLE> quad;

  dg_basis_base () = default;

  explicit dg_basis_base (int dunavant_rule): quad (dunavant_rule)
  {
  }

  std::pair<t8_mra::mat, std::vector<size_t>>
  trafo_matrix_to_ref_element (const double physical_vertices[3][3])
  {
    t8_mra::mat transform_to_ref (3, 3);
    std::vector<size_t> permuation_vec (3, 0u);

    for (auto i = 0; i < 3; ++i)
      for (auto j = 0; j < 3; ++j)
        transform_to_ref (i, j) = i == 2 ? 1.0 : physical_vertices[j][i];
    t8_mra::lu_factors (transform_to_ref, permuation_vec);

    return { transform_to_ref, permuation_vec };
  }

  /// returns list of dereferenced quad points [x_0, y_0, x_1, y_1, ...]
  std::vector<double>
  deref_quad_points (const double physical_vertices[3][3])
  {
    const std::array<double, 6> corners { physical_vertices[0][0], physical_vertices[0][1], physical_vertices[1][0],
                                          physical_vertices[1][1], physical_vertices[2][0], physical_vertices[2][1] };

    return reference_to_physical_t3 (corners, quad.points);
  }

  std::vector<double>
  ref_point (const t8_mra::mat &trafo_mat, const std::vector<size_t> &permuation_vec,
             const std::vector<double> &grid_point)
  {
    std::vector<double> ret = { grid_point[0], grid_point[1], 1.0 };
    t8_mra::lu_solve (trafo_mat, permuation_vec, ret);

    return ret;
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

  // Constructor for triangular elements
  explicit dg_basis (int _dunavant_rule)
    requires (Shape == T8_ECLASS_TRIANGLE)
    : Base (_dunavant_rule)
  {
  }

  // Constructor for cartesian elements (LINE, QUAD, HEX)
  explicit dg_basis (int _num_quad_points_1d)
    requires is_cartesian<Shape>
    : Base (_num_quad_points_1d)
  {
  }

  /// All basis function values at a reference point.
  std::array<double, DOF>
  basis_value (const std::vector<double> &x_ref)
  {
    std::array<double, DIM> x;
    for (unsigned int d = 0; d < DIM; ++d)
      x[d] = x_ref[d];

    return basis_t::eval (x);
  }

  /// grad[dir][i] = d(phi_i)/dx_dir at a reference point (cartesian only).
  std::array<std::array<double, DOF>, DIM>
  basis_gradient (const std::vector<double> &x_ref)
    requires is_cartesian<Shape>
  {
    std::array<double, DIM> x;
    for (unsigned int d = 0; d < DIM; ++d)
      x[d] = x_ref[d];

    return basis_t::eval_gradient (x);
  }
};

}  // namespace t8_mra
#endif
