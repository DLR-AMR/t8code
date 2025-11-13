#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <array>

#include "t8_eclass.h"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/dunavant.hxx"
#include "t8_mra/num/mat.hpp"

namespace t8_mra
{

template <t8_eclass TShape>
struct dg_basis_base
{
  static constexpr unsigned int DIM = 0;
  static constexpr t8_eclass Shape = TShape;

  size_t num_quad_points;
  std::vector<double> ref_quad_points;
  std::vector<double> quad_weights;
  /// TODO error message
};

template <>
struct dg_basis_base<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned int DIM = 2;
  static constexpr t8_eclass Shape = T8_ECLASS_TRIANGLE;

  size_t num_quad_points;
  int dunavant_rule;

  std::vector<double> ref_quad_points;
  std::vector<double> quad_weights;

  dg_basis_base () = default;

  dg_basis_base (int _num_quad_points, int _dunavant_rule)
    : num_quad_points (t8_mra::dunavant_order_num (_dunavant_rule)), dunavant_rule (_dunavant_rule),
      ref_quad_points (2u * num_quad_points, 0.0), quad_weights (num_quad_points, 0.0)
  {
    t8_mra::dunavant_rule (dunavant_rule, num_quad_points, ref_quad_points.data (), quad_weights.data ());
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
    std::vector<double> deref_quad_points (ref_quad_points.size (), 0.0);

    std::array<double, 6> corners { physical_vertices[0][0], physical_vertices[0][1], physical_vertices[1][0],
                                    physical_vertices[1][1], physical_vertices[2][0], physical_vertices[2][1] };

    t8_mra::reference_to_physical_t3 (corners.data (), num_quad_points, ref_quad_points.data (),
                                      deref_quad_points.data ());

    return deref_quad_points;
  }

  std::vector<double>
  ref_point (const t8_mra::mat &trafo_mat, const std::vector<size_t> &permuation_vec,
             const std::vector<double> &grid_point)
  {
    std::vector<double> ret = { grid_point[0], grid_point[1], 1.0 };
    t8_mra::lu_solve (trafo_mat, permuation_vec, ret);

    // CRITICAL FIX: ret contains barycentric coordinates [λ0, λ1, λ2]
    // Reference coordinates (xi, eta) correspond to (λ1, λ2), NOT (λ0, λ1)!
    // OLD (WRONG): return ret;  // This returned [λ0, λ1, λ2]
    // NEW (CORRECT): Return [λ1, λ2, ...] as (xi, eta) coordinates
    // return { ret[1], ret[2], ret[0] };  // Reorder to [λ1, λ2, λ0]
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

 public:
  dg_basis () = default;

  dg_basis (int _num_quad_points, int _dunavant_rule)
    requires (Shape == T8_ECLASS_TRIANGLE)
    : Base (_num_quad_points, _dunavant_rule)
  {
  }

  std::array<double, DOF>
  basis_value (const std::vector<double> &x_ref)
  {
    std::array<double, DOF> res;
    ///TODO scaling_functions für allgemeine shapes machen
    for (auto i = 0u; i < DOF; ++i)
      res[i] = t8_mra::skalierungsfunktion (i, x_ref[0], x_ref[1]);

    return res;
  }
};

}  // namespace t8_mra
#endif
