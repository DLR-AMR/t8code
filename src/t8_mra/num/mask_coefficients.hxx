#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include <t8_eclass/t8_eclass.h>
#include <t8_mra/num/mat.hxx>
#include <t8_mra/core/shape_traits.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

namespace t8_mra
{

// ============================================================================
// Two-scale low-pass mask coefficients
// ============================================================================
// One routine for every shape: the parent -> child-k prolongation matrix
//   M_k(i, j) = norm * ∫_ref φ_i(ξ) φ_j(Φ_k ξ) dξ
// with φ = basis<TShape, P> (the reference function space) and Φ_k the affine
// map of the reference element onto child k of a uniform refinement. Row i is a
// child dof, column j a parent dof. Consumed by the multiscale transform
// (mst.hxx). Adding a shape only needs child_maps<> + the basis/quadrature
// specializations.

/// Affine map ξ -> A ξ + b on the reference element.
template <int DIM>
struct affine_map
{
  std::array<std::array<double, DIM>, DIM> A {};
  std::array<double, DIM> b {};

  std::array<double, DIM>
  operator() (const std::array<double, DIM> &xi) const
  {
    std::array<double, DIM> out {};
    for (int r = 0; r < DIM; ++r) {
      out[r] = b[r];
      for (int c = 0; c < DIM; ++c)
        out[r] += A[r][c] * xi[c];
    }
    return out;
  }
};

/// Cartesian children: 2^DIM axis-aligned half-cells, Φ_k(ξ) = (s_k + ξ)/2.
template <t8_eclass TShape>
  requires is_cartesian<TShape>
auto
child_maps ()
{
  constexpr int DIM = shape_traits<TShape>::DIM;
  constexpr int NC = shape_traits<TShape>::NUM_CHILDREN;

  std::array<affine_map<DIM>, NC> maps {};
  for (int k = 0; k < NC; ++k) {
    for (int r = 0; r < DIM; ++r) {
      maps[k].A[r][r] = 0.5;
      maps[k].b[r] = 0.5 * ((k >> r) & 1);
    }
  }
  return maps;
}

/// Triangle children: red refinement into 3 corner triangles + 1 inverted
/// centre; per-child vertex order fixes the two-scale convention.
template <t8_eclass TShape>
  requires (TShape == T8_ECLASS_TRIANGLE)
auto
child_maps ()
{
  using vertex = std::array<double, 2>;
  constexpr vertex p0 { 0.0, 0.0 }, p1 { 1.0, 0.0 }, p2 { 0.0, 1.0 };
  constexpr vertex m01 { 0.5, 0.0 }, m02 { 0.0, 0.5 }, m12 { 0.5, 0.5 };
  const std::array<std::array<vertex, 3>, 4> verts { {
    { m01, m12, m02 },  // centre (inverted)
    { m01, p1, m12 },
    { m12, p2, m02 },
    { m02, p0, m01 },
  } };

  std::array<affine_map<2>, 4> maps {};
  for (int k = 0; k < 4; ++k) {
    const auto &v = verts[k];
    maps[k].b = v[0];
    for (int r = 0; r < 2; ++r) {
      maps[k].A[r][0] = v[1][r] - v[0][r];
      maps[k].A[r][1] = v[2][r] - v[0][r];
    }
  }
  return maps;
}

/// Mask normalization: cartesian basis is orthonormal on the unit cell (vol 1);
/// the triangle factor 1/4 = 1/2 (two-scale definition) * 1/2 (reference area).
template <t8_eclass TShape>
inline constexpr double mask_norm = is_cartesian<TShape> ? 1.0 : 0.25;

/// Reference quadrature parameter exact to degree 2(P-1): Gauss points per axis
/// (cartesian) or Dunavant rule (triangle).
template <t8_eclass TShape, int P>
inline constexpr int mask_quad_param = is_cartesian<TShape> ? (P + 1) : std::min (20, 2 * P);

/// Compute the NUM_CHILDREN two-scale masks for shape TShape at order P.
template <t8_eclass TShape, int P>
void
compute_mask (std::vector<t8_mra::mat> &mask)
{
  using basis_t = basis<TShape, P>;
  constexpr int DIM = basis_t::DIM;
  constexpr int DOF = basis_t::DOF;
  constexpr int NC = shape_traits<TShape>::NUM_CHILDREN;

  mask.assign (NC, t8_mra::mat { DOF, DOF });

  const auto children = child_maps<TShape> ();
  const quadrature<TShape> quad (mask_quad_param<TShape, P>);

  double wsum = 0.0;
  for (std::size_t q = 0; q < quad.num_points; ++q)
    wsum += quad.weights[q];
  const double scale = mask_norm<TShape> / wsum;

  for (int k = 0; k < NC; ++k) {
    for (std::size_t q = 0; q < quad.num_points; ++q) {
      std::array<double, DIM> xi {};
      for (int d = 0; d < DIM; ++d)
        xi[d] = quad.points[DIM * q + d];

      const auto child_val = basis_t::eval (xi);
      const auto parent_val = basis_t::eval (children[k] (xi));
      const double w = quad.weights[q];

      for (int i = 0; i < DOF; ++i)
        for (int j = 0; j < DOF; ++j)
          mask[k](i, j) += w * child_val[i] * parent_val[j];
    }

    for (int i = 0; i < DOF; ++i)
      for (int j = 0; j < DOF; ++j)
        mask[k](i, j) *= scale;
  }
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
