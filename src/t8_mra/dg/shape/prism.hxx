#pragma once

#ifdef T8_ENABLE_MRA

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
    const auto num_q = basis.quad.num_points;
    std::vector<std::array<double, DOF>> basis_at_quad (num_q);
    std::vector<std::array<double, U_DIM>> f_at_quad (num_q);

    for (auto j = 0u; j < num_q; ++j) {
      const std::array<double, 3> ref { basis.quad.points[3 * j], basis.quad.points[3 * j + 1],
                                        basis.quad.points[3 * j + 2] };
      const auto phys = geom.to_physical (ref);
      basis_at_quad[j] = basis.basis_value (geom.basis_coord (ref));

      f_at_quad[j] = func (phys[0], phys[1], phys[2]);
    }

    for (auto i = 0u; i < DOF; ++i) {
      std::array<double, U_DIM> sum = {};

      for (auto j = 0u; j < num_q; ++j)
        for (auto u = 0u; u < U_DIM; ++u)
          sum[u] += basis.quad.weights[j] * f_at_quad[j][u] * geom.basis_scale * basis_at_quad[j][i];

      for (auto u = 0u; u < U_DIM; ++u)
        coeffs[element_t::dg_idx (u, i)] = sum[u] * geom.volume;
    }
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
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
