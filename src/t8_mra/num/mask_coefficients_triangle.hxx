#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include <t8_eclass/t8_eclass.h>
#include <t8_mra/num/mat.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

namespace t8_mra
{

/// Two-scale low-pass masks for the triangle (4 children, dof x dof each):
///   M_k(i, j) = 1/2 * ∫_T φ_j(Φ_k ξ) φ_i(ξ) dξ
/// with φ the orthonormal Dubiner basis and Φ_k the affine map from the
/// reference triangle T onto child k of a red refinement. M_k is the parent ->
/// child-k two-scale relation consumed by the multiscale transform (mst.hxx).
template <int P>
void
initialize_triangle_mask (std::vector<t8_mra::mat> &mask_coeffs)
{
  using basis_t = basis<T8_ECLASS_TRIANGLE, P>;
  constexpr int dof = basis_t::DOF;

  mask_coeffs.assign (4, t8_mra::mat { dof, dof });

  // Child vertices in parent reference coords (red refinement); per-child vertex
  // order fixes the canonical two-scale convention.
  constexpr std::array<std::array<std::array<double, 2>, 3>, 4> child { {
    { { { 0.5, 0.0 }, { 0.5, 0.5 }, { 0.0, 0.5 } } },  // center (inverted)
    { { { 0.5, 0.0 }, { 1.0, 0.0 }, { 0.5, 0.5 } } },
    { { { 0.5, 0.5 }, { 0.0, 1.0 }, { 0.0, 0.5 } } },
    { { { 0.0, 0.5 }, { 0.0, 0.0 }, { 0.5, 0.0 } } },
  } };

  // Reference-triangle rule exact to degree 2(P-1) (product of two basis funcs).
  const quadrature<T8_ECLASS_TRIANGLE> quad (std::min (20, 2 * P));

  double wsum = 0.0;
  for (std::size_t j = 0; j < quad.num_points; ++j)
    wsum += quad.weights[j];
  // 1/2 (definition) * area(T)=1/2 / wsum, independent of the rule's weight scale.
  const double factor = 1.0 / (4.0 * wsum);

  for (int k = 0; k < 4; ++k) {
    const auto &v = child[k];
    for (std::size_t j = 0; j < quad.num_points; ++j) {
      const double x = quad.points[2 * j];
      const double y = quad.points[2 * j + 1];
      const double l0 = 1.0 - x - y;

      const std::array<double, 2> mapped { l0 * v[0][0] + x * v[1][0] + y * v[2][0],
                                           l0 * v[0][1] + x * v[1][1] + y * v[2][1] };
      const auto phi_parent = basis_t::eval (mapped);
      const auto phi_child = basis_t::eval ({ x, y });
      const double w = quad.weights[j];

      for (int i = 0; i < dof; ++i)
        for (int jj = 0; jj < dof; ++jj)
          mask_coeffs[k](i, jj) += w * phi_parent[jj] * phi_child[i];
    }

    for (int i = 0; i < dof; ++i)
      for (int jj = 0; jj < dof; ++jj)
        mask_coeffs[k](i, jj) *= factor;
  }
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
