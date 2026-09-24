#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <span>
#include <vector>

#include "t8_eclass/t8_eclass.h"

#include "t8_mra/data/element_data.hxx"
#include "t8_mra/dg/dg_base.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/dg_basis.hxx"

namespace t8_mra
{

/**
 * @brief Prism DG numerics: Dubiner-times-Legendre basis, Dunavant-times-
 * Gauss-Legendre projection, affine extruded geometry.
 */
template <int U, int P>
class dg<T8_ECLASS_PRISM, U, P> {
 public:
  static constexpr t8_eclass Shape = T8_ECLASS_PRISM;
  using element_t = element_data<Shape, U, P>;
  using geometry_t = cell_geometry<Shape, P>;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int DOF = element_t::DOF;

  /// Accuracy degree of both quadrature factors; 2P covers products of two
  /// order-P basis functions with margin for general data.
  static constexpr int default_quadrature_rule = 2 * P;

  dg_basis<element_t> basis;

  explicit dg (int quadrature_rule = default_quadrature_rule): basis (quadrature_rule)
  {
    basis_at_quad.resize (basis.quad.num_points);

    for (auto q = 0u; q < basis.quad.num_points; ++q)
      basis_at_quad[q] = basis.basis_value (geometry_t::basis_coord (quad_point (q)));
  }

  /** @brief Cell geometry from native corner coords, volume and reference vertex order. */
  [[nodiscard]] geometry_t
  geometry (const std::array<std::array<double, 3>, T8_ECLASS_MAX_CORNERS> &corners, double volume,
            const std::array<int, 3> &order) const
  {
    /// t8code numbers prism corners as (triangle vertex, extrusion end), so the
    /// triangle's order applies unchanged to both ends.
    constexpr int BASE_VERTICES = 3;
    std::array<std::array<double, 3>, 2 * BASE_VERTICES> ordered;

    for (auto i = 0; i < BASE_VERTICES; ++i) {
      ordered[order[i]] = corners[i];
      ordered[order[i] + BASE_VERTICES] = corners[i + BASE_VERTICES];
    }

    return geometry_t::from_prism (ordered[0], ordered[1], ordered[2], ordered[3], ordered[4], ordered[5], volume);
  }

  /** @brief Project func onto the DG basis by the reference prism rule. */
  template <typename Func>
  void
  project (std::span<double> coeffs, const geometry_t &geom, Func &&func)
  {
    std::ranges::fill (coeffs, 0.0);

    for (auto q = 0u; q < basis.quad.num_points; ++q) {
      const auto phys = geom.to_physical (quad_point (q));
      const auto f_val = func (phys[0], phys[1], phys[2]);
      const auto &phi = basis_at_quad[q];
      const auto weight = basis.quad.weights[q];

      for (auto i = 0u; i < DOF; ++i) {
        const auto weighted_phi = weight * phi[i];

        for (auto u = 0u; u < U_DIM; ++u)
          coeffs[element_t::dg_idx (u, i)] += weighted_phi * f_val[u];
      }
    }

    const auto scale = geom.basis_scale * geom.volume;
    for (auto &coeff : coeffs)
      coeff *= scale;
  }

  /** @brief Solution value per component at a physical point. */
  [[nodiscard]] std::array<double, U_DIM>
  evaluate (const geometry_t &geom, const element_t &data, const std::array<double, DIM> &x_phys) const
  {
    const auto x_ref = geom.to_reference (x_phys);
    std::array<double, U_DIM> res = {};

    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = geom.value (std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), x_ref);

    return res;
  }

  /** @brief Solution gradient grad[u][d] = d(u_u)/d(x_d) at a physical point. */
  [[nodiscard]] std::array<std::array<double, DIM>, U_DIM>
  evaluate_gradient (const geometry_t &geom, const element_t &data, const std::array<double, DIM> &x_phys) const
  {
    const auto x_ref = geom.to_reference (x_phys);
    std::array<std::array<double, DIM>, U_DIM> grad = {};

    for (auto u = 0u; u < U_DIM; ++u)
      grad[u] = geom.gradient (std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), x_ref);

    return grad;
  }

 private:
  /// Basis values at the fixed quadrature points, built once with the rule.
  std::vector<std::array<double, DOF>> basis_at_quad;

  /// Reference coordinates of quadrature point q.
  [[nodiscard]] std::array<double, DIM>
  quad_point (unsigned int q) const
  {
    return { basis.quad.points[DIM * q], basis.quad.points[DIM * q + 1], basis.quad.points[DIM * q + 2] };
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
