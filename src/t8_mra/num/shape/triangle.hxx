#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/basis/dubiner.hxx"
#include "t8_mra/num/cell_geometry.hxx"

#include <array>
#include <cmath>
#include <span>
#include <utility>

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

  static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    return [&]<std::size_t... I> (std::index_sequence<I...>) {
      return std::array<double, DOF> { scaling_function<static_cast<int> (I)> (x[0], x[1])... };
    }(std::make_index_sequence<DOF> {});
  }

  /// grad[dir][i] = d(phi_i)/dx_dir on the reference triangle (the geometric
  /// Jacobian to physical coordinates is applied by the caller).
  static std::array<std::array<double, DOF>, DIM>
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

  static double
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
  std::array<point, 2> edges {};    // x_d = origin_d + sum_e edges[d][e] * ref_e
  std::array<point, 2> inv_jac {};  // ref_e = sum_d inv_jac[e][d] * (x_d - origin_d)
  double volume = 0.0;
  double basis_scale = 0.0;
  double mass = 0.0;
  int level = 0;

  /** @brief Build from the ordered vertices (origin, r0 edge, r1 edge). */
  static cell_geometry
  from_triangle (const point &v0, const point &v1, const point &v2, double vol)
  {
    cell_geometry geom;
    geom.origin = v0;
    geom.edges = { point { v1[0] - v0[0], v2[0] - v0[0] }, point { v1[1] - v0[1], v2[1] - v0[1] } };

    const double J00 = geom.edges[0][0], J01 = geom.edges[0][1], J10 = geom.edges[1][0], J11 = geom.edges[1][1];
    const double det = J00 * J11 - J01 * J10;
    geom.inv_jac = { point { J11 / det, -J01 / det }, point { -J10 / det, J00 / det } };

    geom.volume = vol;
    geom.basis_scale = basis_t::normalization (vol);
    geom.mass = geom.basis_scale * geom.basis_scale * std::abs (det);

    return geom;
  }

  /** @brief Reference (r0, r1) -> Dubiner coordinate {lambda0, lambda1}. */
  static point
  basis_coord (const point &ref)
  {
    return { 1.0 - ref[0] - ref[1], ref[0] };
  }

  /** @brief Whether a reference point lies in the unit triangle. */
  static bool
  in_ref_cell (const point &ref)
  {
    return ref[0] >= -reference_cell_tol && ref[1] >= -reference_cell_tol
           && ref[0] + ref[1] <= 1.0 + reference_cell_tol;
  }

  /** @brief Physical -> reference coordinate. */
  point
  to_reference (const point &phys) const
  {
    const double dx = phys[0] - origin[0], dy = phys[1] - origin[1];

    return { inv_jac[0][0] * dx + inv_jac[0][1] * dy, inv_jac[1][0] * dx + inv_jac[1][1] * dy };
  }

  /** @brief Reference -> physical coordinate. */
  point
  to_physical (const point &ref) const
  {
    return { origin[0] + edges[0][0] * ref[0] + edges[0][1] * ref[1],
             origin[1] + edges[1][0] * ref[0] + edges[1][1] * ref[1] };
  }

  /** @brief Whether a physical point lies in the cell. */
  bool
  contains (const point &phys) const
  {
    return in_ref_cell (to_reference (phys));
  }

  /** @brief basis_scale * sum_i coeffs_i * phi_i(basis_coord(ref)). */
  static double
  eval_modal (std::span<const double> coeffs, const point &ref, double basis_scale)
  {
    const auto phi = basis_t::eval (basis_coord (ref));
    double sum = 0.0;
    for (int i = 0; i < DOF; ++i)
      sum += coeffs[i] * phi[i];

    return basis_scale * sum;
  }

  /** @brief Physical value at a reference point from the cell volume alone (no cell map). */
  static double
  reference_value (std::span<const double> coeffs, const point &ref, double volume)
  {
    return eval_modal (coeffs, ref, basis_t::normalization (volume));
  }

  /** @brief Physical value at a reference point using the cached basis scale. */
  double
  value (std::span<const double> coeffs, const point &ref) const
  {
    return eval_modal (coeffs, ref, basis_scale);
  }

  /** @brief Physical gradient d(u_h)/d(x_d) at a reference point. */
  point
  gradient (std::span<const double> coeffs, const point &ref) const
  {
    const auto ref_grad = to_ref_grad (basis_t::eval_gradient (basis_coord (ref)));
    point grad {};
    for (int d = 0; d < DIM; ++d) {
      double sum = 0.0;
      for (int i = 0; i < DOF; ++i)
        sum += coeffs[i] * (ref_grad[0][i] * inv_jac[0][d] + ref_grad[1][i] * inv_jac[1][d]);
      grad[d] = basis_scale * sum;
    }

    return grad;
  }

  /** @brief inv_jac * phys_dir: weights w with (phys_dir . grad_x phi) = w . grad_r phi. */
  point
  reference_direction (const point &phys_dir) const
  {
    return { inv_jac[0][0] * phys_dir[0] + inv_jac[0][1] * phys_dir[1],
             inv_jac[1][0] * phys_dir[0] + inv_jac[1][1] * phys_dir[1] };
  }

 private:
  /** @brief Basis gradient d/dlambda -> reference d/dr (r0=lambda1, r1=lambda2). */
  static std::array<std::array<double, DOF>, 2>
  to_ref_grad (const std::array<std::array<double, DOF>, 2> &basis_grad)
  {
    std::array<std::array<double, DOF>, 2> ref_grad {};
    for (int i = 0; i < DOF; ++i) {
      ref_grad[0][i] = basis_grad[1][i] - basis_grad[0][i];
      ref_grad[1][i] = -basis_grad[0][i];
    }

    return ref_grad;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
