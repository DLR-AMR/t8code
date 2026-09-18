#pragma once

#include <algorithm>
#include <numeric>
#ifdef T8_ENABLE_MRA

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/basis/legendre.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/mask_coefficients.hxx"
#include "t8_mra/num/quadrature/gauss_legendre.hxx"
#include "t8_mra/num/quadrature/quadrature.hxx"

namespace t8_mra
{

/// Cartesian shapes (LINE, QUAD, HEX): tensor product of 1D Legendre modes, the
/// basis index decomposed lexicographically (first coordinate fastest).
/// Orthonormal on the reference cell, so no volume normalization.
template <t8_eclass TShape, int P>
  requires is_cartesian<TShape>
struct basis<TShape, P>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;
  static constexpr int DOF = shape_traits<TShape>::dof (P);

  [[nodiscard]] static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    std::array<double, DOF> res = {};
    for (int p = 0; p < DOF; ++p) {
      double v = 1.0;
      int idx = p;
      for (int d = 0; d < DIM; ++d) {
        v *= phi_1d (x[d], idx % P);
        idx /= P;
      }
      res[p] = v;
    }
    return res;
  }

  [[nodiscard]] static std::array<std::array<double, DOF>, DIM>
  eval_gradient (const std::array<double, DIM> &x)
  {
    std::array<std::array<double, DOF>, DIM> grad = {};
    for (int dir = 0; dir < DIM; ++dir) {
      for (int p = 0; p < DOF; ++p) {
        double v = 1.0;
        int idx = p;
        for (int d = 0; d < DIM; ++d) {
          const int deg = idx % P;
          idx /= P;
          v *= (d == dir) ? phi_prime_1d<P> (x[d], deg) : phi_1d (x[d], deg);
        }
        grad[dir][p] = v;
      }
    }
    return grad;
  }

  [[nodiscard]] static constexpr double
  normalization (double /*vol*/) noexcept
  {
    return 1.0;
  }
};

/** @brief Cartesian leaf geometry: axis-aligned box, diagonal Jacobian. */
template <t8_eclass TShape, int P>
  requires is_cartesian<TShape>
struct cell_geometry<TShape, P>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;
  static constexpr int DOF = shape_traits<TShape>::dof (P);
  using basis_t = basis<TShape, P>;
  using point = std::array<double, DIM>;

  point origin {};
  point extent {};
  double volume = 0.0;
  double basis_scale = 1.0;
  double mass = 0.0;
  int level = 0;

  /** @brief Build from the cell's min/max corners. */
  [[nodiscard]] static cell_geometry
  from_box (const point &min_corner, const point &max_corner, double vol)
  {
    cell_geometry geom;
    geom.origin = min_corner;

    double det = 1.0;
    for (int d = 0; d < DIM; ++d) {
      geom.extent[d] = max_corner[d] - min_corner[d];
      det *= geom.extent[d];
    }

    geom.volume = vol;
    geom.basis_scale = basis_t::normalization (vol);
    geom.mass = geom.basis_scale * geom.basis_scale * std::abs (det);

    return geom;
  }

  /** @brief Surface over volume, the length scale an interior-penalty face term runs on. */
  [[nodiscard]] double
  surface_to_volume () const
  {
    const auto sum
      = std::accumulate (extent.begin (), extent.end (), 0.0, [] (double res, double x) { return res + 1.0 / x; });

    return 2.0 * sum;
  }

  /** @brief Reference coordinate -> basis coordinate (identity). */
  [[nodiscard]] static point
  basis_coord (const point &ref)
  {
    return ref;
  }

  /** @brief Barycentre of the reference cell. */
  [[nodiscard]] static point
  reference_centroid ()
  {
    point centre;
    centre.fill (0.5);

    return centre;
  }

  /** @brief Pin a shared-face point's normal reference component exactly (t8 face: axis=f>>1, side=f&1). */
  [[nodiscard]] static point
  on_face (point ref, int face)
  {
    ref[face >> 1] = ((face & 1) != 0) ? 1.0 : 0.0;
    return ref;
  }

  /** @brief Whether a reference point lies in the unit box. */
  [[nodiscard]] static bool
  in_ref_cell (const point &ref)
  {
    return !std::any_of (ref.begin (), ref.end (),
                         [&] (double x) { return x < -reference_cell_tol || x > 1.0 + reference_cell_tol; });
  }

  /** @brief Physical -> reference coordinate. */
  [[nodiscard]] point
  to_reference (const point &phys) const
  {
    point ref {};
    for (int d = 0; d < DIM; ++d)
      ref[d] = (phys[d] - origin[d]) / extent[d];

    return ref;
  }

  /** @brief Reference -> physical coordinate. */
  [[nodiscard]] point
  to_physical (const point &ref) const
  {
    point phys {};
    for (int d = 0; d < DIM; ++d)
      phys[d] = origin[d] + extent[d] * ref[d];

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
    double sum = 0.0;
    for (int i = 0; i < DOF; ++i)
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
    const auto ref_grad = basis_t::eval_gradient (ref);
    point grad {};

    for (int d = 0; d < DIM; ++d) {
      double sum = 0.0;
      for (int i = 0; i < DOF; ++i)
        sum += coeffs[i] * ref_grad[d][i];

      grad[d] = basis_scale * sum / extent[d];
    }

    return grad;
  }

  /** @brief inv_jac * phys_dir: weights w with (phys_dir . grad_x phi) = w . grad_r phi. */
  [[nodiscard]] point
  reference_direction (const point &phys_dir) const
  {
    point ref_dir {};

    for (int d = 0; d < DIM; ++d)
      ref_dir[d] = phys_dir[d] / extent[d];

    return ref_dir;
  }
};

/// Cartesian shapes: tensor product of a 1D Gauss-Legendre rule, exact to
/// degree 2*num_points_1d - 1 per axis.
template <t8_eclass TShape>
struct quadrature<TShape, std::enable_if_t<is_cartesian<TShape>>>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;

  std::size_t num_points = 0;
  std::vector<double> points;  /// flattened: point q coord d at points[DIM*q + d]
  std::vector<double> weights;

  /// 1D point count for a rule exact to the given polynomial degree (2n-1 >= degree).
  [[nodiscard]] static constexpr int
  rule_for_degree (int degree)
  {
    return degree / 2 + 1;
  }

  quadrature () = default;

  explicit quadrature (int num_points_1d)
  {
    std::vector<double> p1d;
    std::vector<double> w1d;
    gauss_legendre_1d (num_points_1d, p1d, w1d);

    num_points = 1;
    for (auto d = 0; d < DIM; ++d)
      num_points *= num_points_1d;

    points.resize (DIM * num_points);
    weights.resize (num_points);

    /// Iterating over the DIM axes (first axis fastest)
    for (auto q = 0u; q < num_points; ++q) {
      auto rest = q;
      double w = 1.0;

      for (auto d = 0; d < DIM; ++d) {
        const auto id = rest % num_points_1d;
        rest /= num_points_1d;
        points[DIM * q + d] = p1d[id];
        w *= w1d[id];
      }

      weights[q] = w;
    }
  }
};

/// Cartesian two-scale policy: the basis is orthonormal on the unit cell and
/// does not rescale with the volume
template <t8_eclass TShape>
  requires is_cartesian<TShape>
struct mask_policy<TShape>
{
  static constexpr double norm = 1.0;

  /// 2^DIM axis-aligned half-cells, Phi_k(xi) = (s_k + xi) / 2.
  [[nodiscard]] static auto
  child_maps ()
  {
    constexpr int DIM = shape_traits<TShape>::DIM;
    constexpr int NC = shape_traits<TShape>::NUM_CHILDREN;

    std::array<affine_map<DIM>, NC> maps {};
    for (auto k = 0; k < NC; ++k)
      for (auto r = 0; r < DIM; ++r) {
        maps[k].A[r][r] = 0.5;
        maps[k].b[r] = 0.5 * ((k >> r) & 1);
      }

    return maps;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
