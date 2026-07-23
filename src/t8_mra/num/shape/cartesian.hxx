#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/basis/legendre.hxx"
#include "t8_mra/num/cell_geometry.hxx"

#include <array>
#include <cmath>
#include <span>

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

  static std::array<double, DOF>
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

  static std::array<std::array<double, DOF>, DIM>
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

  static constexpr double
  normalization (double) noexcept
  {
    return 1.0;
  }
};

/** @brief Cartesian leaf geometry: axis-aligned box, diagonal Jacobian. */
template <t8_eclass Shape, int P>
  requires is_cartesian<Shape>
struct cell_geometry<Shape, P>
{
  static constexpr int DIM = shape_traits<Shape>::DIM;
  static constexpr int DOF = shape_traits<Shape>::dof (P);
  using basis_t = basis<Shape, P>;
  using point = std::array<double, DIM>;

  point origin {};
  point extent {};
  double volume = 0.0;
  double basis_scale = 1.0;
  double mass = 0.0;
  int level = 0;

  /** @brief Build from the cell's min/max corners. */
  static cell_geometry
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

  /** @brief Reference coordinate -> basis coordinate (identity). */
  static point
  basis_coord (const point &ref)
  {
    return ref;
  }

  /** @brief Pin a shared-face point's normal reference component exactly (t8 face: axis=f>>1, side=f&1). */
  static point
  on_face (point ref, int face)
  {
    ref[face >> 1] = (face & 1) ? 1.0 : 0.0;
    return ref;
  }

  /** @brief Whether a reference point lies in the unit box. */
  static bool
  in_ref_cell (const point &ref)
  {
    for (int d = 0; d < DIM; ++d)
      if (ref[d] < -reference_cell_tol || ref[d] > 1.0 + reference_cell_tol)
        return false;

    return true;
  }

  /** @brief Physical -> reference coordinate. */
  point
  to_reference (const point &phys) const
  {
    point ref {};
    for (int d = 0; d < DIM; ++d)
      ref[d] = (phys[d] - origin[d]) / extent[d];

    return ref;
  }

  /** @brief Reference -> physical coordinate. */
  point
  to_physical (const point &ref) const
  {
    point phys {};
    for (int d = 0; d < DIM; ++d)
      phys[d] = origin[d] + extent[d] * ref[d];

    return phys;
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
  point
  reference_direction (const point &phys_dir) const
  {
    point ref_dir {};

    for (int d = 0; d < DIM; ++d)
      ref_dir[d] = phys_dir[d] / extent[d];

    return ref_dir;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
