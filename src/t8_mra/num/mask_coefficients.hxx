#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>

#include <t8_eclass/t8_eclass.h>

#include <t8_mra/core/shape_traits.hxx>

namespace t8_mra
{

/// ============================================================================
/// Two-scale low-pass mask coefficients
/// ============================================================================
/// The parent -> child-k prolongation matrix
///   M_k(i, j) = norm * ∫_ref φ_i(ξ) φ_j(Φ_k ξ) dξ
/// with φ = basis<TShape, P> (the reference function space) and Φ_k the affine
/// map of the reference element onto child k of a uniform refinement. Row i is a
/// child dof, column j a parent dof.

/// Affine map ξ -> A ξ + b on the reference element.
template <int DIM>
struct affine_map
{
  std::array<std::array<double, DIM>, DIM> A {};
  std::array<double, DIM> b {};

  [[nodiscard]] std::array<double, DIM>
  operator() (const std::array<double, DIM> &xi) const
  {
    std::array<double, DIM> out {};

    for (auto r = 0; r < DIM; ++r) {
      out[r] = b[r];
      for (auto c = 0; c < DIM; ++c)
        out[r] += A[r][c] * xi[c];
    }
    return out;
  }
};

/// Per-shape mask normalization and the reference maps of the children.
/// Specialize for a new shape.
template <t8_eclass TShape>
struct mask_policy;

}  // namespace t8_mra

// The basis and the quadrature pull in the per-shape specializations, so they
// can only be reached once the primaries above are declared.
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

namespace t8_mra
{

/**
 * @brief The NUM_CHILDREN two-scale masks of a shape at order P.
 *
 * m[k](i,j) has row i a child dof and column j a parent dof. The transpose is kept
 * alongside because the forward transform contracts over the row index
 */
template <t8_eclass TShape, int P>
struct two_scale_mask
{
  static constexpr int NC = shape_traits<TShape>::NUM_CHILDREN;
  static constexpr int DOF = shape_traits<TShape>::dof (P);

  using matrix = std::array<std::array<double, DOF>, DOF>;

  std::array<matrix, NC> m {};
  std::array<matrix, NC> transposed {};
};

/// Compute the NUM_CHILDREN two-scale masks for shape TShape at order P.
template <t8_eclass TShape, int P>
void
compute_mask (two_scale_mask<TShape, P> &mask)
{
  using basis_t = basis<TShape, P>;
  constexpr int DIM = basis_t::DIM;
  constexpr int DOF = basis_t::DOF;
  constexpr int NC = shape_traits<TShape>::NUM_CHILDREN;

  mask = {};

  const quadrature<TShape> quad (quadrature<TShape>::rule_for_degree (2 * P));
  const auto children = mask_policy<TShape>::child_maps ();

  auto wsum = 0.0;
  for (auto q = 0u; q < quad.num_points; ++q)
    wsum += quad.weights[q];

  const auto scale = mask_policy<TShape>::norm / wsum;

  for (auto k = 0; k < NC; ++k) {
    for (auto q = 0u; q < quad.num_points; ++q) {
      std::array<double, DIM> xi {};

      for (auto d = 0; d < DIM; ++d)
        xi[d] = quad.points[DIM * q + d];

      const auto phi = basis_t::eval (xi);
      const auto phi_mapped = basis_t::eval (children[k](xi));

      const auto w = scale * quad.weights[q];

      for (auto i = 0; i < DOF; ++i)
        for (auto j = 0; j < DOF; ++j)
          mask.m[k][i][j] += w * phi[i] * phi_mapped[j];
    }

    for (auto i = 0; i < DOF; ++i)
      for (auto j = 0; j < DOF; ++j)
        mask.transposed[k][j][i] = mask.m[k][i][j];
  }
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
