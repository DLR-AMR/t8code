#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/dg/dg_base.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/dg_basis.hxx"

#include <array>
#include <span>
#include <vector>

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
  }

  /** @brief Cell geometry from native corner coords, volume and reference vertex order. */
  geometry_t
  geometry (const double corners[T8_ECLASS_MAX_CORNERS][3], double volume, const std::array<int, 3> &order) const
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
    const auto num_q = basis.quad.num_points;
    std::vector<std::array<double, DOF>> basis_at_quad (num_q);
    std::vector<std::array<double, U_DIM>> f_at_quad (num_q);
    for (auto j = 0u; j < num_q; ++j) {
      const std::array<double, 2> ref { basis.quad.points[2 * j], basis.quad.points[2 * j + 1] };
      const auto phys = geom.to_physical (ref);
      basis_at_quad[j] = basis.basis_value (geom.basis_coord (ref));
      f_at_quad[j] = func (phys[0], phys[1]);
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
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
