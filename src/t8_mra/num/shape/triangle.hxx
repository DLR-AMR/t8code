#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/basis/dubiner.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/mask_coefficients.hxx"
#include "t8_mra/num/quadrature/dunavant.hxx"
#include "t8_mra/num/quadrature/quadrature.hxx"

namespace t8_mra
{

/// Triangle: the orthonormal Dubiner basis in barycentric coords
/// (x = {lambda0, lambda1}). Orthonormal on the reference triangle (area 1/2),
/// so the physical basis scales by sqrt(1/(2*vol)).
template <int P>
struct basis<T8_ECLASS_TRIANGLE, P>
{
  static constexpr int DIM = 2;
  static constexpr int DOF = shape_traits<T8_ECLASS_TRIANGLE>::dof (P);

  [[nodiscard]] static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    return [&]<std::size_t... I> (std::index_sequence<I...>) {
      return std::array<double, DOF> { scaling_function<static_cast<int> (I)> (x[0], x[1])... };
    }(std::make_index_sequence<DOF> {});
  }

  /// grad[dir][i] = d(phi_i)/dx_dir on the reference triangle (the geometric
  /// Jacobian to physical coordinates is applied by the caller).
  [[nodiscard]] static std::array<std::array<double, DOF>, DIM>
  eval_gradient (const std::array<double, DIM> &x)
  {
    std::array<std::array<double, DOF>, DIM> grad = {};
    [&]<std::size_t... I> (std::index_sequence<I...>) {
      (
        [&] {
          const auto g = scaling_function_gradient<static_cast<int> (I)> (x[0], x[1]);
          grad[0][I] = g[0];
          grad[1][I] = g[1];
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

/** @brief Triangle leaf geometry: general affine map from three ordered vertices. */
template <int P>
struct cell_geometry<T8_ECLASS_TRIANGLE, P>
{
  static constexpr int DIM = 2;
  static constexpr int DOF = shape_traits<T8_ECLASS_TRIANGLE>::dof (P);
  using basis_t = basis<T8_ECLASS_TRIANGLE, P>;
  using point = std::array<double, 2>;

  point origin {};
  std::array<point, 2> edges {};    /// x_d = origin_d + sum_e edges[d][e] * ref_e
  std::array<point, 2> inv_jac {};  /// ref_e = sum_d inv_jac[e][d] * (x_d - origin_d)
  double volume = 0.0;
  double basis_scale = 0.0;
  double mass = 0.0;
  int level = 0;

  /** @brief Build from the ordered vertices (origin, r0 edge, r1 edge). */
  [[nodiscard]] static cell_geometry
  from_triangle (const point &v0, const point &v1, const point &v2, double vol)
  {
    cell_geometry geom;
    geom.origin = v0;
    geom.edges = { point { v1[0] - v0[0], v2[0] - v0[0] }, point { v1[1] - v0[1], v2[1] - v0[1] } };

    const auto J00 = geom.edges[0][0];
    const auto J01 = geom.edges[0][1];
    const auto J10 = geom.edges[1][0];
    const auto J11 = geom.edges[1][1];
    const auto det = J00 * J11 - J01 * J10;
    geom.inv_jac = { point { J11 / det, -J01 / det }, point { -J10 / det, J00 / det } };

    geom.volume = vol;
    geom.basis_scale = basis_t::normalization (vol);
    geom.mass = geom.basis_scale * geom.basis_scale * std::abs (det);

    return geom;
  }

  /** @brief Perimeter over area, the length scale an interior-penalty face term runs on. */
  [[nodiscard]] double
  surface_to_volume () const
  {
    const auto e0 = std::hypot (edges[0][0], edges[1][0]);
    const auto e1 = std::hypot (edges[0][1], edges[1][1]);
    const auto e2 = std::hypot (edges[0][1] - edges[0][0], edges[1][1] - edges[1][0]);

    return (e0 + e1 + e2) / volume;
  }

  /** @brief Reference (r0, r1) -> Dubiner coordinate {lambda0, lambda1}. */
  [[nodiscard]] static point
  basis_coord (const point &ref)
  {
    return { 1.0 - ref[0] - ref[1], ref[0] };
  }

  /** @brief Barycentre of the reference cell. */
  [[nodiscard]] static point
  reference_centroid ()
  {
    return { 1.0 / 3.0, 1.0 / 3.0 };
  }

  /** @brief Whether a reference point lies in the unit triangle. */
  [[nodiscard]] static bool
  in_ref_cell (const point &ref)
  {
    return ref[0] >= -reference_cell_tol && ref[1] >= -reference_cell_tol
           && ref[0] + ref[1] <= 1.0 + reference_cell_tol;
  }

  /** @brief Physical -> reference coordinate. */
  [[nodiscard]] point
  to_reference (const point &phys) const
  {
    const auto dx = phys[0] - origin[0];
    const auto dy = phys[1] - origin[1];

    return { inv_jac[0][0] * dx + inv_jac[0][1] * dy, inv_jac[1][0] * dx + inv_jac[1][1] * dy };
  }

  /** @brief Reference -> physical coordinate. */
  [[nodiscard]] point
  to_physical (const point &ref) const
  {
    return { origin[0] + edges[0][0] * ref[0] + edges[0][1] * ref[1],
             origin[1] + edges[1][0] * ref[0] + edges[1][1] * ref[1] };
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
    for (int d = 0; d < DIM; ++d) {
      auto sum = 0.0;
      for (auto i = 0; i < DOF; ++i)
        sum += coeffs[i] * (ref_grad[0][i] * inv_jac[0][d] + ref_grad[1][i] * inv_jac[1][d]);

      grad[d] = basis_scale * sum;
    }

    return grad;
  }

  /** @brief inv_jac * phys_dir: weights w with (phys_dir . grad_x phi) = w . grad_r phi. */
  [[nodiscard]] point
  reference_direction (const point &phys_dir) const
  {
    return { inv_jac[0][0] * phys_dir[0] + inv_jac[0][1] * phys_dir[1],
             inv_jac[1][0] * phys_dir[0] + inv_jac[1][1] * phys_dir[1] };
  }

 private:
  /** @brief Basis gradient d/dlambda -> reference d/dr (r0=lambda1, r1=lambda2). */
  [[nodiscard]] static std::array<std::array<double, DOF>, 2>
  to_ref_grad (const std::array<std::array<double, DOF>, 2> &basis_grad)
  {
    std::array<std::array<double, DOF>, 2> ref_grad {};

    for (auto i = 0; i < DOF; ++i) {
      ref_grad[0][i] = basis_grad[1][i] - basis_grad[0][i];
      ref_grad[1][i] = -basis_grad[0][i];
    }

    return ref_grad;
  }
};

/// Triangle: a Dunavant rule on the reference triangle.
template <>
struct quadrature<T8_ECLASS_TRIANGLE>
{
  static constexpr int DIM = 2;

  std::size_t num_points = 0;
  std::vector<double> points;  // flattened: [x0, y0, x1, y1, ...]
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
    auto rule_data = dunavant_rule (rule);
    num_points = rule_data.weights.size ();
    points = std::move (rule_data.points);
    weights = std::move (rule_data.weights);
  }
};

/// Triangle two-scale policy. The basis normalizes as sqrt(1/(2*vol)), so a
/// child coefficient sits at 1/sqrt(NUM_CHILDREN) of its parent's and the mask
/// factor is 1/(2*sqrt(4)).
template <>
struct mask_policy<T8_ECLASS_TRIANGLE>
{
  static constexpr double norm = 0.25;

  /// Refinement into 3 corner triangles + 1 inverted centre
  [[nodiscard]] static auto
  child_maps ()
  {
    using vertex = std::array<double, 2>;
    constexpr vertex p0 { 0.0, 0.0 };
    constexpr vertex p1 { 1.0, 0.0 };
    constexpr vertex p2 { 0.0, 1.0 };

    constexpr vertex m01 { 0.5, 0.0 };
    constexpr vertex m02 { 0.0, 0.5 };
    constexpr vertex m12 { 0.5, 0.5 };

    const std::array<std::array<vertex, 3>, 4> verts { {
      { m01, m12, m02 },  /// centre (inverted)
      { m01, p1, m12 },
      { m12, p2, m02 },
      { m02, p0, m01 },
    } };

    std::array<affine_map<2>, 4> maps {};
    for (auto k = 0; k < 4; ++k) {
      const auto &v = verts[k];
      maps[k].b = v[0];

      for (auto r = 0; r < 2; ++r) {
        maps[k].A[r][0] = v[1][r] - v[0][r];
        maps[k].A[r][1] = v[2][r] - v[0][r];
      }
    }

    return maps;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
