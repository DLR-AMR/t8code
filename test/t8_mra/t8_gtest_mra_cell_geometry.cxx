
#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/num/cell_geometry.hxx>
#include <t8_mra/data/element_data.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/dg_basis.hxx>

#include <array>
#include <span>
#include <vector>

namespace
{

constexpr double eps = 1e-12;

/* Cartesian: from_box gives the per-axis affine map [0,1]^DIM <-> physical box */

template <t8_eclass Shape, int P>
void
check_box (const typename t8_mra::cell_geometry<Shape, P>::point &lo,
           const typename t8_mra::cell_geometry<Shape, P>::point &hi)
{
  using geom_t = t8_mra::cell_geometry<Shape, P>;
  constexpr int DIM = geom_t::DIM;

  const auto geom = geom_t::from_box (lo, hi, 1.0);

  for (const double s : { 0.0, 0.25, 0.5, 1.0 }) {
    typename geom_t::point ref;
    for (int d = 0; d < DIM; ++d)
      ref[d] = s;

    const auto phys = geom.to_physical (ref);
    for (int d = 0; d < DIM; ++d)
      EXPECT_NEAR (phys[d], lo[d] + ref[d] * (hi[d] - lo[d]), eps);

    const auto back = geom.to_reference (phys);
    for (int d = 0; d < DIM; ++d)
      EXPECT_NEAR (back[d], ref[d], eps);

    EXPECT_TRUE (geom.contains (phys));
  }
}

TEST (mra_cell_geometry_cartesian, box_map_all_dims)
{
  check_box<T8_ECLASS_LINE, 3> ({ 2.0 }, { 5.0 });
  check_box<T8_ECLASS_QUAD, 3> ({ 1.0, 2.0 }, { 4.0, 5.0 });
  check_box<T8_ECLASS_HEX, 2> ({ 1.0, 2.0, 3.0 }, { 4.0, 6.0, 9.0 });
}

TEST (mra_cell_geometry_cartesian, box_map_negative_coordinates)
{
  check_box<T8_ECLASS_QUAD, 3> ({ -2.0, -1.0 }, { 1.0, 0.5 });
}

/* Triangle: from_triangle gives the affine map, vertices -> reference corners, and
 * physical <-> reference round trips (r0, r1) = barycentric (lambda1, lambda2). */

using tri_point = t8_mra::cell_geometry<T8_ECLASS_TRIANGLE, 2>::point;

void
check_triangle (const tri_point &v0, const tri_point &v1, const tri_point &v2)
{
  using geom_t = t8_mra::cell_geometry<T8_ECLASS_TRIANGLE, 2>;
  const auto geom = geom_t::from_triangle (v0, v1, v2, 1.0);

  // Vertices map to the reference corners (0,0), (1,0), (0,1).
  const std::array<std::pair<tri_point, tri_point>, 3> vertex_ref {
    { { v0, { 0.0, 0.0 } }, { v1, { 1.0, 0.0 } }, { v2, { 0.0, 1.0 } } }
  };
  for (const auto &[vertex, ref] : vertex_ref) {
    const auto got = geom.to_reference (vertex);
    EXPECT_NEAR (got[0], ref[0], eps);
    EXPECT_NEAR (got[1], ref[1], eps);
  }

  // A point built from known barycentric weights inverts to (r0, r1) = (w1, w2).
  const std::array<std::array<double, 3>, 3> weights {
    { { 0.5, 0.3, 0.2 }, { 0.1, 0.6, 0.3 }, { 1.0 / 3, 1.0 / 3, 1.0 / 3 } }
  };
  for (const auto &w : weights) {
    const tri_point phys { w[0] * v0[0] + w[1] * v1[0] + w[2] * v2[0], w[0] * v0[1] + w[1] * v1[1] + w[2] * v2[1] };

    const auto ref = geom.to_reference (phys);
    EXPECT_NEAR (ref[0], w[1], eps);
    EXPECT_NEAR (ref[1], w[2], eps);

    // Forward map reproduces the same physical point.
    const auto phys_back = geom.to_physical (ref);
    EXPECT_NEAR (phys_back[0], phys[0], eps);
    EXPECT_NEAR (phys_back[1], phys[1], eps);

    EXPECT_TRUE (geom.contains (phys));
  }
}

TEST (mra_cell_geometry_triangle, affine_map_over_many_triangles)
{
  check_triangle ({ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 });   // reference
  check_triangle ({ 0.5, 1.0 }, { 2.0, 0.5 }, { 1.0, 3.0 });   // generic
  check_triangle ({ 0.5, 1.0 }, { 1.0, 3.0 }, { 2.0, 0.5 });   // reversed winding
  check_triangle ({ -1.0, -1.0 }, { 1.0, -2.0 }, { 0.0, 2.0 });  // negative coordinates
}

/* Constant mode -> constant field. */
template <typename Geom>
void
check_value_constant_mode (const Geom &geom)
{
  std::array<double, Geom::DOF> coeffs {};
  coeffs[0] = 2.5;
  const std::span<const double> c (coeffs.data (), Geom::DOF);

  EXPECT_NEAR (geom.value (c, { 0.5, 0.1 }), geom.value (c, { 0.2, 0.3 }), eps);
}

/* Affine field (P1 modes `lin`): constant gradient, exact value(pA)-value(pB)=grad.(pA-pB). */
template <typename Geom>
void
check_affine_gradient (const Geom &geom, std::array<int, 2> lin)
{
  using point = typename Geom::point;
  std::array<double, Geom::DOF> coeffs {};
  coeffs[0] = 1.0;
  coeffs[lin[0]] = 0.3;
  coeffs[lin[1]] = -0.2;
  const std::span<const double> c (coeffs.data (), Geom::DOF);

  const point refA { 0.2, 0.3 }, refB { 0.5, 0.15 };
  const auto gA = geom.gradient (c, refA);
  const auto gB = geom.gradient (c, refB);
  for (int d = 0; d < Geom::DIM; ++d)
    EXPECT_NEAR (gA[d], gB[d], eps);

  const auto pA = geom.to_physical (refA), pB = geom.to_physical (refB);
  double rhs = 0.0;
  for (int d = 0; d < Geom::DIM; ++d)
    rhs += gA[d] * (pA[d] - pB[d]);
  EXPECT_NEAR (geom.value (c, refA) - geom.value (c, refB), rhs, eps);
}

/* reference_direction(v) = to_reference(origin + v). */
template <typename Geom>
void
check_reference_direction (const Geom &geom)
{
  using point = typename Geom::point;
  const point v { 0.7, -0.4 };
  const auto got = geom.reference_direction (v);
  const auto expected = geom.to_reference ({ geom.origin[0] + v[0], geom.origin[1] + v[1] });
  for (int d = 0; d < Geom::DIM; ++d)
    EXPECT_NEAR (got[d], expected[d], eps);
}

TEST (mra_cell_geometry, value_gradient_reference_direction)
{
  const auto quad = t8_mra::cell_geometry<T8_ECLASS_QUAD, 3>::from_box ({ 1.0, 2.0 }, { 4.0, 6.0 }, 12.0);
  const auto tri = t8_mra::cell_geometry<T8_ECLASS_TRIANGLE, 3>::from_triangle ({ 0.5, 1.0 }, { 2.0, 0.5 }, { 1.0, 3.0 },
                                                                                1.625);

  check_value_constant_mode (quad);
  check_value_constant_mode (tri);

  check_affine_gradient (quad, { 1, 3 });  // tensor P1 modes: x-linear = 1, y-linear = P = 3
  check_affine_gradient (tri, { 1, 2 });   // Dubiner P1 modes

  check_reference_direction (quad);
  check_reference_direction (tri);
}

/* dg_basis.basis_value / basis_gradient forward to the reference basis */

template <t8_eclass Shape, int P>
void
check_basis_forward (int quad_param, const std::vector<double> &x_ref)
{
  using elem = t8_mra::element_data<Shape, 1, P>;
  using basis_t = t8_mra::basis<Shape, P>;
  constexpr int DIM = t8_mra::shape_traits<Shape>::DIM;

  t8_mra::dg_basis<elem> dg_basis (quad_param);

  std::array<double, DIM> x {};
  for (int d = 0; d < DIM; ++d)
    x[d] = x_ref[d];

  const auto val = dg_basis.basis_value (x);
  const auto ref_val = basis_t::eval (x);
  for (std::size_t i = 0; i < val.size (); ++i)
    EXPECT_NEAR (val[i], ref_val[i], eps);

  const auto grad = dg_basis.basis_gradient (x);
  const auto ref_grad = basis_t::eval_gradient (x);
  for (int dir = 0; dir < DIM; ++dir)
    for (std::size_t i = 0; i < grad[dir].size (); ++i)
      EXPECT_NEAR (grad[dir][i], ref_grad[dir][i], eps);
}

TEST (mra_dg_basis, value_and_gradient_forward_all_shapes)
{
  check_basis_forward<T8_ECLASS_LINE, 3> (3, { 0.4 });
  check_basis_forward<T8_ECLASS_QUAD, 3> (3, { 0.3, 0.7 });
  check_basis_forward<T8_ECLASS_HEX, 2> (3, { 0.3, 0.6, 0.2 });
  check_basis_forward<T8_ECLASS_TRIANGLE, 3> (4, { 0.25, 0.4 });  // r0+r1 < 1
}

}  // namespace

#endif  // T8_ENABLE_MRA
