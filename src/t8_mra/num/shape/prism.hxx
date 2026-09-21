#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <utility>
#include <vector>

#include <t8.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/basis/prism.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/mask_coefficients.hxx"
#include "t8_mra/num/quadrature/dunavant.hxx"
#include "t8_mra/num/quadrature/gauss_legendre.hxx"
#include "t8_mra/num/quadrature/quadrature.hxx"
#include "t8_mra/num/shape/triangle.hxx"

namespace t8_mra
{

/// Prism: the orthonormal Dubiner basis of the triangle factor times the
/// Legendre mode of the line factor (x = {lambda0, lambda1, z}). Orthonormal on
/// the reference prism (volume 1/2), so the physical basis scales by
/// sqrt(1/(2*vol)) exactly as the triangle does.
template <int P>
struct basis<T8_ECLASS_PRISM, P>
{
  static constexpr int DIM = 3;
  static constexpr int DOF = shape_traits<T8_ECLASS_PRISM>::dof (P);

  [[nodiscard]] static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    return [&]<std::size_t... I> (std::index_sequence<I...>) {
      return std::array<double, DOF> { prism_scaling_function<P, static_cast<int> (I)> (x[0], x[1], x[2])... };
    }(std::make_index_sequence<DOF> {});
  }

  /// grad[dir][i] = d(phi_i)/dx_dir on the reference prism (the geometric
  /// Jacobian to physical coordinates is applied by the caller).
  [[nodiscard]] static std::array<std::array<double, DOF>, DIM>
  eval_gradient (const std::array<double, DIM> &x)
  {
    std::array<std::array<double, DOF>, DIM> grad = {};
    [&]<std::size_t... I> (std::index_sequence<I...>) {
      (
        [&] {
          const auto g = prism_scaling_function_gradient<P, static_cast<int> (I)> (x[0], x[1], x[2]);
          for (int dir = 0; dir < DIM; ++dir)
            grad[dir][I] = g[dir];
        }(),
        ...);
    }(std::make_index_sequence<DOF> {});
    return grad;
  }

  [[nodiscard]] static double
  normalization (double vol)
  {
    return std::sqrt (1.0 / (2.0 * vol));
  }
};

/**
 * @brief Prism leaf geometry: affine map from six ordered vertices.
 *
 * Only extruded prisms are affine
 */
template <int P>
struct cell_geometry<T8_ECLASS_PRISM, P>
{
  static constexpr int DIM = 3;
  static constexpr int DOF = shape_traits<T8_ECLASS_PRISM>::dof (P);
  using basis_t = basis<T8_ECLASS_PRISM, P>;
  using point = std::array<double, 3>;

  point origin {};
  std::array<point, 3> edges {};    /// x_d = origin_d + sum_e edges[d][e] * ref_e
  std::array<point, 3> inv_jac {};  /// ref_e = sum_d inv_jac[e][d] * (x_d - origin_d)
  double volume = 0.0;
  double basis_scale = 0.0;
  double mass = 0.0;
  int level = 0;

  /** @brief Build from the ordered vertices (base triangle v0,v1,v2 and its extrusion v3,v4,v5). */
  [[nodiscard]] static cell_geometry
  from_prism (const point &v0, const point &v1, const point &v2, const point &v3, [[maybe_unused]] const point &v4,
              [[maybe_unused]] const point &v5, double vol)
  {
    T8_ASSERT (is_extruded (v0, v1, v2, v3, v4, v5));

    cell_geometry geom;
    geom.origin = v0;
    for (int d = 0; d < 3; ++d)
      geom.edges[d] = { v1[d] - v0[d], v2[d] - v0[d], v3[d] - v0[d] };

    const auto &J = geom.edges;
    const point cof0 { J[1][1] * J[2][2] - J[1][2] * J[2][1], J[1][2] * J[2][0] - J[1][0] * J[2][2],
                       J[1][0] * J[2][1] - J[1][1] * J[2][0] };
    const auto det = J[0][0] * cof0[0] + J[0][1] * cof0[1] + J[0][2] * cof0[2];

    geom.inv_jac = { point { cof0[0] / det, (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det,
                             (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det },
                     point { cof0[1] / det, (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det,
                             (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det },
                     point { cof0[2] / det, (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det,
                             (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det } };

    geom.volume = vol;
    geom.basis_scale = basis_t::normalization (vol);
    geom.mass = geom.basis_scale * geom.basis_scale * std::abs (det);

    return geom;
  }

  /** @brief Surface over volume, the length scale an interior-penalty face term runs on. */
  [[nodiscard]] double
  surface_to_volume () const
  {
    const point r0 { edges[0][0], edges[1][0], edges[2][0] };
    const point r1 { edges[0][1], edges[1][1], edges[2][1] };
    const point extrusion { edges[0][2], edges[1][2], edges[2][2] };
    const point diagonal { r1[0] - r0[0], r1[1] - r0[1], r1[2] - r0[2] };

    const auto triangles = norm (cross (r0, r1));
    const auto lateral
      = norm (cross (r0, extrusion)) + norm (cross (r1, extrusion)) + norm (cross (diagonal, extrusion));

    return (triangles + lateral) / volume;
  }

  /** @brief Reference (r0, r1, r2) -> basis coordinate {lambda0, lambda1, z}. */
  [[nodiscard]] static point
  basis_coord (const point &ref)
  {
    return { 1.0 - ref[0] - ref[1], ref[0], ref[2] };
  }

  /** @brief Barycentre of the reference cell. */
  [[nodiscard]] static point
  reference_centroid ()
  {
    return { 1.0 / 3.0, 1.0 / 3.0, 0.5 };
  }

  /** @brief Whether a reference point lies in the unit prism. */
  [[nodiscard]] static bool
  in_ref_cell (const point &ref)
  {
    return ref[0] >= -reference_cell_tol && ref[1] >= -reference_cell_tol && ref[0] + ref[1] <= 1.0 + reference_cell_tol
           && ref[2] >= -reference_cell_tol && ref[2] <= 1.0 + reference_cell_tol;
  }

  /** @brief Physical -> reference coordinate. */
  [[nodiscard]] point
  to_reference (const point &phys) const
  {
    const point offset { phys[0] - origin[0], phys[1] - origin[1], phys[2] - origin[2] };
    point ref {};

    for (auto e = 0; e < DIM; ++e)
      ref[e] = inv_jac[e][0] * offset[0] + inv_jac[e][1] * offset[1] + inv_jac[e][2] * offset[2];

    return ref;
  }

  /** @brief Reference -> physical coordinate. */
  [[nodiscard]] point
  to_physical (const point &ref) const
  {
    point phys {};
    for (auto d = 0; d < DIM; ++d)
      phys[d] = origin[d] + edges[d][0] * ref[0] + edges[d][1] * ref[1] + edges[d][2] * ref[2];

    return phys;
  }

  /** @brief Whether a physical point lies in the cell. */
  [[nodiscard]] bool
  contains (const point &phys) const
  {
    return in_ref_cell (to_reference (phys));
  }

  /** @brief basis_scale * sum_i coeffs_i * phi_i(basis_coord(ref)). */
  [[nodiscard]] static double
  eval_modal (std::span<const double> coeffs, const point &ref, double basis_scale)
  {
    const auto phi = basis_t::eval (basis_coord (ref));
    auto sum = 0.0;

    for (auto i = 0; i < DOF; ++i)
      sum += coeffs[i] * phi[i];

    return basis_scale * sum;
  }

  /** @brief Physical value at a reference point from the cell volume alone (no cell map). */
  [[nodiscard]] static double
  reference_value (std::span<const double> coeffs, const point &ref, double volume)
  {
    return eval_modal (coeffs, ref, basis_t::normalization (volume));
  }

  /** @brief Physical value at a reference point using the cached basis scale. */
  [[nodiscard]] double
  value (std::span<const double> coeffs, const point &ref) const
  {
    return eval_modal (coeffs, ref, basis_scale);
  }

  /** @brief Physical gradient d(u_h)/d(x_d) at a reference point. */
  [[nodiscard]] point
  gradient (std::span<const double> coeffs, const point &ref) const
  {
    const auto ref_grad = to_ref_grad (basis_t::eval_gradient (basis_coord (ref)));
    point grad {};

    for (auto d = 0; d < DIM; ++d) {
      auto sum = 0.0;

      for (auto i = 0; i < DOF; ++i)
        sum += coeffs[i]
               * (ref_grad[0][i] * inv_jac[0][d] + ref_grad[1][i] * inv_jac[1][d] + ref_grad[2][i] * inv_jac[2][d]);

      grad[d] = basis_scale * sum;
    }

    return grad;
  }

  /** @brief inv_jac * phys_dir: weights w with (phys_dir . grad_x phi) = w . grad_r phi. */
  [[nodiscard]] point
  reference_direction (const point &phys_dir) const
  {
    point ref_dir {};

    for (auto e = 0; e < DIM; ++e)
      ref_dir[e] = inv_jac[e][0] * phys_dir[0] + inv_jac[e][1] * phys_dir[1] + inv_jac[e][2] * phys_dir[2];

    return ref_dir;
  }

 private:
  /** @brief Whether the three extrusion vectors agree, the condition for an affine prism. */
  [[nodiscard]] static bool
  is_extruded (const point &v0, const point &v1, const point &v2, const point &v3, const point &v4, const point &v5)
  {
    auto scale = 0.0;
    auto defect = 0.0;

    for (auto d = 0; d < 3; ++d) {
      const auto extrusion = v3[d] - v0[d];
      scale = std::max (scale, std::abs (extrusion));
      defect = std::max ({ defect, std::abs (v4[d] - v1[d] - extrusion), std::abs (v5[d] - v2[d] - extrusion) });
    }

    return defect <= 1e-9 * std::max (scale, 1.0);
  }

  [[nodiscard]] static point
  cross (const point &a, const point &b)
  {
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
  }

  [[nodiscard]] static double
  norm (const point &a)
  {
    return std::sqrt (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  }

  /** @brief Basis gradient d/dlambda -> reference d/dr (r0=lambda1, r1=lambda2, r2=z). */
  [[nodiscard]] static std::array<std::array<double, DOF>, 3>
  to_ref_grad (const std::array<std::array<double, DOF>, 3> &basis_grad)
  {
    std::array<std::array<double, DOF>, 3> ref_grad {};

    for (auto i = 0; i < DOF; ++i) {
      ref_grad[0][i] = basis_grad[1][i] - basis_grad[0][i];
      ref_grad[1][i] = -basis_grad[0][i];
      ref_grad[2][i] = basis_grad[2][i];
    }

    return ref_grad;
  }
};

/// Prism: a Dunavant rule on the triangle factor crossed with Gauss-Legendre on
/// the line factor, both taken exact to the same degree.
template <>
struct quadrature<T8_ECLASS_PRISM>
{
  static constexpr int DIM = 3;

  std::size_t num_points = 0;
  std::vector<double> points;  /// flattened: [x0, y0, z0, x1, y1, z1, ...]
  std::vector<double> weights;

  /// Dunavant rule (accuracy degree) for the given polynomial degree, capped at the table maximum.
  [[nodiscard]] static constexpr int
  rule_for_degree (int degree)
  {
    return std::min (20, degree);
  }

  quadrature () = default;

  explicit quadrature (int rule)
  {
    const auto triangle_rule = dunavant_rule (rule);
    const int num_points_1d = rule / 2 + 1;

    std::vector<double> points_1d;
    std::vector<double> weights_1d;
    gauss_legendre_1d (num_points_1d, points_1d, weights_1d);

    const auto num_triangle_points = triangle_rule.weights.size ();
    num_points = num_triangle_points * num_points_1d;
    points.resize (DIM * num_points);
    weights.resize (num_points);

    /// Triangle index fastest
    auto q = 0u;
    for (auto iz = 0; iz < num_points_1d; ++iz) {
      for (auto it = 0u; it < num_triangle_points; ++it, ++q) {
        points[DIM * q] = triangle_rule.points[2 * it];
        points[DIM * q + 1] = triangle_rule.points[2 * it + 1];
        points[DIM * q + 2] = points_1d[iz];

        weights[q] = triangle_rule.weights[it] * weights_1d[iz];
      }
    }
  }
};

/// Prism two-scale policy, the triangle's rule with NUM_CHILDREN 8: norm = 1/(2*sqrt(8)).
template <>
struct mask_policy<T8_ECLASS_PRISM>
{
  static constexpr double norm = std::numbers::sqrt2 / 8.0;

  /// The triangle's refinement in (x,y) crossed with a bisection in z,
  /// indexed as t8_dprism_child does (triangle id fastest).
  [[nodiscard]] static auto
  child_maps ()
  {
    constexpr auto NC = shape_traits<T8_ECLASS_PRISM>::NUM_CHILDREN;
    constexpr auto TRIANGLE_NC = shape_traits<T8_ECLASS_TRIANGLE>::NUM_CHILDREN;
    const auto triangle = mask_policy<T8_ECLASS_TRIANGLE>::child_maps ();

    std::array<affine_map<3>, NC> maps {};

    for (auto k = 0; k < NC; ++k) {
      // t8_dprism_child splits the child id as tri_id + 4 * line_id, so the
      // integer quotient is the half of the z bisection this child sits in.
      const auto tri_child = k % TRIANGLE_NC;
      const auto line_child = k / TRIANGLE_NC;
      const auto &base = triangle[tri_child];

      for (auto r = 0; r < 2; ++r) {
        maps[k].b[r] = base.b[r];
        for (auto c = 0; c < 2; ++c)
          maps[k].A[r][c] = base.A[r][c];
      }

      maps[k].A[2][2] = 0.5;
      maps[k].b[2] = 0.5 * line_child;
    }

    return maps;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
