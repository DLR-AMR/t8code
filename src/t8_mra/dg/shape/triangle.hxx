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
 * @brief Triangle DG numerics: orthonormal Dubiner basis, Dunavant projection,
 * affine barycentric geometry. t8code-free (vertex order supplied by the caller).
 */
template <int U, int P>
class dg<T8_ECLASS_TRIANGLE, U, P> {
 public:
  static constexpr t8_eclass Shape = T8_ECLASS_TRIANGLE;
  using element_t = element_data<Shape, U, P>;
  using geometry_t = cell_geometry<Shape, P>;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int DOF = element_t::DOF;

  /// Dunavant rule number == polynomial exactness; 2P covers products of two
  /// order-P basis functions with margin for general data.
  static constexpr int default_quadrature_rule = 2 * P;

  dg_basis<element_t> basis;

  explicit dg (int dunavant_rule = default_quadrature_rule): basis (dunavant_rule)
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
    std::array<std::array<double, 2>, 3> ordered;
    for (int i = 0; i < 3; ++i)
      ordered[order[i]] = { corners[i][0], corners[i][1] };

    return geometry_t::from_triangle (ordered[0], ordered[1], ordered[2], volume);
  }

  /** @brief Project func onto the DG basis by Dunavant quadrature. */
  template <typename Func>
  void
  project (std::span<double> coeffs, const geometry_t &geom, Func &&func)
  {
    std::ranges::fill (coeffs, 0.0);

    for (auto q = 0u; q < basis.quad.num_points; ++q) {
      const auto phys = geom.to_physical (quad_point (q));
      const auto f_val = func (phys[0], phys[1]);
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
    return { basis.quad.points[DIM * q], basis.quad.points[DIM * q + 1] };
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
