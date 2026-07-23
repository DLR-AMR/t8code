#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/dg/dg_base.hxx"
#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/dg_basis.hxx"
#include "t8_mra/num/geometry.hxx"

#include <array>
#include <span>
#include <type_traits>
#include <vector>

namespace t8_mra
{

/**
 * @brief Cartesian DG numerics (LINE, QUAD, HEX): tensor Legendre basis,
 * Gauss-Legendre projection, axis-aligned box geometry. t8code-free.
 */
template <t8_eclass Shape, int U, int P>
  requires is_cartesian<Shape>
class dg<Shape, U, P> {
 public:
  using element_t = element_data<Shape, U, P>;
  using geometry_t = cell_geometry<Shape, P>;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int DOF = element_t::DOF;

  /// Gauss-Legendre points per axis: n points integrate degree 2n-1 exactly.
  static constexpr int default_quadrature_rule = P + 1;

  dg_basis<element_t> basis;

  explicit dg (int num_quad_points_1d = default_quadrature_rule): basis (num_quad_points_1d)
  {
  }

  /** @brief Cell geometry from the leaf's t8code-order corner coords and volume. */
  geometry_t
  geometry (const double corners[T8_ECLASS_MAX_CORNERS][3], double volume, const std::array<int, 3> & = {}) const
  {
    // QUAD corners are permuted (t8code swaps 2 and 3) so index 0 is the lower
    // and the last the upper corner, as extract_cartesian_vertices expects.
    double ordered[T8_ECLASS_MAX_CORNERS][3] = {};
    for (int corner = 0; corner < shape_traits<Shape>::NUM_VERTICES; ++corner) {
      int source = corner;
      if constexpr (DIM == 2 && Shape == T8_ECLASS_QUAD) {
        constexpr int quad_corner_order[4] = { 0, 1, 3, 2 };
        source = quad_corner_order[corner];
      }
      for (int d = 0; d < 3; ++d)
        ordered[corner][d] = corners[source][d];
    }

    std::array<double, DIM> lower, upper;
    extract_cartesian_vertices<DIM> (ordered, lower, upper);
    return geometry_t::from_box (lower, upper, volume);
  }

  /** @brief Project func onto the DG basis by Gauss-Legendre quadrature. */
  template <typename Func>
  void
  project (std::span<double> coeffs, const geometry_t &geom, Func &&func)
  {
    const auto num_q = basis.quad.num_points;
    std::vector<std::array<double, DOF>> basis_at_quad (num_q);
    std::vector<std::array<double, DIM>> phys_at_quad (num_q);
    std::array<double, DIM> x_ref;
    for (auto q = 0u; q < num_q; ++q) {
      for (unsigned int d = 0; d < DIM; ++d)
        x_ref[d] = basis.quad.points[DIM * q + d];
      basis_at_quad[q] = basis.basis_value (x_ref);
      phys_at_quad[q] = geom.to_physical (x_ref);
    }

    for (auto i = 0u; i < DOF; ++i) {
      std::array<double, U_DIM> sum = {};
      for (auto q = 0u; q < num_q; ++q) {
        const auto f_val = eval_func (func, phys_at_quad[q]);
        for (auto u = 0u; u < U_DIM; ++u)
          sum[u] += basis.quad.weights[q] * f_val[u] * basis_at_quad[q][i];
      }
      for (auto u = 0u; u < U_DIM; ++u)
        coeffs[element_t::dg_idx (u, i)] = sum[u];
    }
  }

  /** @brief Solution value per component at a physical point. */
  std::array<double, U_DIM>
  evaluate (const geometry_t &geom, const element_t &data, const std::array<double, DIM> &x_phys) const
  {
    const auto x_ref = geom.to_reference (x_phys);
    std::array<double, U_DIM> res = {};
    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = geom.value (std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), x_ref);

    return res;
  }

  /** @brief Solution gradient grad[u][d] = d(u_u)/d(x_d) at a physical point. */
  std::array<std::array<double, DIM>, U_DIM>
  evaluate_gradient (const geometry_t &geom, const element_t &data, const std::array<double, DIM> &x_phys) const
  {
    const auto x_ref = geom.to_reference (x_phys);
    std::array<std::array<double, DIM>, U_DIM> grad = {};
    for (auto u = 0u; u < U_DIM; ++u)
      grad[u] = geom.gradient (std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), x_ref);

    return grad;
  }

 private:
  /// Evaluate func at a physical point; supports func(x{,y,z}) returning an
  /// array or writing into an out pointer.
  template <typename Func>
  static std::array<double, U_DIM>
  eval_func (Func &&func, const std::array<double, DIM> &x)
  {
    std::array<double, U_DIM> f_val;
    if constexpr (DIM == 1) {
      if constexpr (std::is_invocable_v<Func, double>)
        f_val = func (x[0]);
      else
        func (x[0], f_val.data ());
    }
    else if constexpr (DIM == 2) {
      if constexpr (std::is_invocable_v<Func, double, double>)
        f_val = func (x[0], x[1]);
      else
        func (x[0], x[1], f_val.data ());
    }
    else {
      if constexpr (std::is_invocable_v<Func, double, double, double>)
        f_val = func (x[0], x[1], x[2]);
      else
        func (x[0], x[1], x[2], f_val.data ());
    }
    return f_val;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
