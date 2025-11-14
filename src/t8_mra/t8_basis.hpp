#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <array>
#include <type_traits>

#include "t8_eclass.h"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/dunavant.hxx"
#include "t8_mra/num/mat.hpp"
#include "t8_mra/num/quadrature.hxx"
#include "t8_mra/num/legendre_basis.hxx"
#include "t8_mra/num/multiindex.hxx"
#include "t8_mra/num/geometry.hxx"

namespace t8_mra
{

template <t8_eclass TShape>
concept is_cartesian = (TShape == T8_ECLASS_LINE || TShape == T8_ECLASS_QUAD || TShape == T8_ECLASS_HEX);

template <t8_eclass TShape, typename = void>
struct dg_basis_base
{
  static constexpr unsigned int DIM = 0;
  static constexpr t8_eclass Shape = TShape;

  size_t num_quad_points;
  std::vector<double> ref_quad_points;
  std::vector<double> quad_weights;
  /// TODO error message
};

template <t8_eclass T>
struct dg_basis_base<T, std::enable_if_t<is_cartesian<T>>>
{
  static constexpr unsigned int DIM = T == T8_ECLASS_LINE ? 1 : (T == T8_ECLASS_QUAD ? 2 : 3);
  static constexpr t8_eclass Shape = T;

  int num_quad_points_1d;  // Number of 1D quadrature points
  size_t num_quad_points;  // Total number of quad points (num_quad_points_1d^DIM)
  int P;                   // Polynomial order (stored for pset generation)

  std::vector<double> ref_quad_points_1d;  // 1D quadrature points
  std::vector<double> quad_weights_1d;     // 1D quadrature weights

  std::vector<double> ref_quad_points;  // Multi-D quadrature points (flattened)
  std::vector<double> quad_weights;     // Multi-D quadrature weights

  dg_basis_base () = default;

  /**
   * @brief Constructor for cartesian elements using Gauss-Legendre quadrature
   *
   * @param _num_quad_points_1d Number of quadrature points in each dimension
   * @param _P Polynomial order for basis (used to generate pset later in dg_basis)
   */
  dg_basis_base (int _num_quad_points_1d, int _P): num_quad_points_1d (_num_quad_points_1d), P (_P)
  {
    // Generate 1D Gauss-Legendre quadrature on [0,1]
    t8_mra::gauss_legendre_1d (num_quad_points_1d, ref_quad_points_1d, quad_weights_1d);

    // Build multi-dimensional quadrature via tensor product
    build_tensor_quadrature ();
  }

 private:
  /**
   * @brief Builds multi-dimensional quadrature from 1D quadrature via tensor product
   */
  void
  build_tensor_quadrature ()
  {
    if constexpr (DIM == 1) {
      num_quad_points = num_quad_points_1d;
      ref_quad_points = ref_quad_points_1d;
      quad_weights = quad_weights_1d;
    }
    else if constexpr (DIM == 2) {
      num_quad_points = num_quad_points_1d * num_quad_points_1d;
      ref_quad_points.resize (2 * num_quad_points);
      quad_weights.resize (num_quad_points);

      size_t idx = 0;
      for (int i = 0; i < num_quad_points_1d; ++i) {
        for (int j = 0; j < num_quad_points_1d; ++j) {
          ref_quad_points[2 * idx] = ref_quad_points_1d[i];
          ref_quad_points[2 * idx + 1] = ref_quad_points_1d[j];
          quad_weights[idx] = quad_weights_1d[i] * quad_weights_1d[j];
          ++idx;
        }
      }
    }
    else if constexpr (DIM == 3) {
      num_quad_points = num_quad_points_1d * num_quad_points_1d * num_quad_points_1d;
      ref_quad_points.resize (3 * num_quad_points);
      quad_weights.resize (num_quad_points);

      size_t idx = 0;
      for (int i = 0; i < num_quad_points_1d; ++i) {
        for (int j = 0; j < num_quad_points_1d; ++j) {
          for (int k = 0; k < num_quad_points_1d; ++k) {
            ref_quad_points[3 * idx] = ref_quad_points_1d[i];
            ref_quad_points[3 * idx + 1] = ref_quad_points_1d[j];
            ref_quad_points[3 * idx + 2] = ref_quad_points_1d[k];
            quad_weights[idx] = quad_weights_1d[i] * quad_weights_1d[j] * quad_weights_1d[k];
            ++idx;
          }
        }
      }
    }
  }

 public:
  /**
   * @brief Computes Jacobian determinant for the physical element
   *
   * @param physical_vertices Vertex coordinates of the physical element
   * @return double Absolute value of Jacobian determinant
   */
  double
  jacobian_det (const double physical_vertices[][3])
  {
    std::array<double, DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<DIM> (physical_vertices, vertices_min, vertices_max);

    return jacobian_determinant<DIM> (vertices_min, vertices_max);
  }

  /**
   * @brief Maps quadrature points from reference element to physical element
   *
   * @param physical_vertices Vertex coordinates of the physical element
   * @return std::vector<double> Physical quadrature points (flattened)
   */
  std::vector<double>
  deref_quad_points (const double physical_vertices[][3])
  {
    std::array<double, DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<DIM> (physical_vertices, vertices_min, vertices_max);

    return transform_quad_points<DIM> (ref_quad_points, num_quad_points, vertices_min, vertices_max);
  }

  /**
   * @brief Maps a point from physical element to reference element [0,1]^DIM
   *
   * @param physical_vertices Vertex coordinates of the physical element
   * @param grid_point Point in physical element (must have at least DIM coordinates)
   * @return std::vector<double> Point in reference element [0,1]^DIM
   */
  std::vector<double>
  ref_point (const double physical_vertices[][3], const std::vector<double> &grid_point)
  {
    std::array<double, DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<DIM> (physical_vertices, vertices_min, vertices_max);

    std::array<double, DIM> x_phys;
    for (unsigned int d = 0; d < DIM; ++d)
      x_phys[d] = grid_point[d];

    std::array<double, DIM> x_ref = ref<DIM> (x_phys, vertices_min, vertices_max);

    std::vector<double> result (DIM);
    for (unsigned int d = 0; d < DIM; ++d)
      result[d] = x_ref[d];

    return result;
  }

  /**
   * @brief Maps a point from reference element to physical element
   *
   * @param physical_vertices Vertex coordinates of the physical element
   * @param ref_point_coords Point in reference element [0,1]^DIM
   * @return std::vector<double> Point in physical element
   */
  std::vector<double>
  deref_point (const double physical_vertices[][3], const std::vector<double> &ref_point_coords)
  {
    std::array<double, DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<DIM> (physical_vertices, vertices_min, vertices_max);

    std::array<double, DIM> x_ref;
    for (unsigned int d = 0; d < DIM; ++d) {
      x_ref[d] = ref_point_coords[d];
    }

    std::array<double, DIM> x_phys = deref<DIM> (x_ref, vertices_min, vertices_max);

    std::vector<double> result (DIM);
    for (unsigned int d = 0; d < DIM; ++d) {
      result[d] = x_phys[d];
    }
    return result;
  }
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

  // Multiindex set for tensor basis (cartesian elements only)

 public:
  std::vector<multiindex<DIM>> pset;
  dg_basis () = default;

  // Constructor for triangular elements
  dg_basis (int _num_quad_points, int _dunavant_rule)
    requires (Shape == T8_ECLASS_TRIANGLE)
    : Base (_num_quad_points, _dunavant_rule)
  {
  }

  // Constructor for cartesian elements (LINE, QUAD, HEX)
  dg_basis (int _num_quad_points_1d, int _P)
    requires is_cartesian<Shape>
    : Base (_num_quad_points_1d, _P)
  {
    // Generate multiindex set for tensor-structured basis
    pset = generate_tensor_pset<DIM> (_P);
  }

  std::array<double, DOF>
  basis_value (const std::vector<double> &x_ref)
  {
    std::array<double, DOF> res;

    if constexpr (is_cartesian<Shape>) {
      // Cartesian elements: use tensor product of Legendre polynomials
      std::array<double, DIM> x_array;
      for (unsigned int d = 0; d < DIM; ++d) {
        x_array[d] = x_ref[d];
      }

      for (auto i = 0u; i < DOF; ++i) {
        res[i] = eval_tensor_basis<DIM> (x_array, i, pset, phi_1d);
      }
    }
    else if constexpr (Shape == T8_ECLASS_TRIANGLE) {
      // Triangle elements: use existing scaling functions
      for (auto i = 0u; i < DOF; ++i)
        res[i] = t8_mra::skalierungsfunktion (i, x_ref[0], x_ref[1]);
    }

    return res;
  }

  /**
   * @brief Evaluates gradient of all basis functions at a point (cartesian elements only)
   *
   * @param x_ref Point in reference element [0,1]^DIM
   * @return std::array<std::array<double, DOF>, DIM> grad[dir][i] = d(phi_i)/dx_dir at x_ref
   */
  std::array<std::array<double, DOF>, DIM>
  basis_gradient (const std::vector<double> &x_ref)
    requires is_cartesian<Shape>
  {
    std::array<std::array<double, DOF>, DIM> grad;

    std::array<double, DIM> x_array;
    for (unsigned int d = 0; d < DIM; ++d) {
      x_array[d] = x_ref[d];
    }

    for (unsigned int dir = 0; dir < DIM; ++dir) {
      for (auto i = 0u; i < DOF; ++i) {
        grad[dir][i] = eval_tensor_basis_gradient<DIM> (x_array, i, dir, pset, phi_1d, phi_prime_1d);
      }
    }

    return grad;
  }
};

}  // namespace t8_mra
#endif
