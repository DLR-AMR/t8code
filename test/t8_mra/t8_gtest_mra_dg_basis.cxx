
#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/data/element_data.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/dg_basis.hxx>

#include <array>
#include <vector>

namespace
{

constexpr double eps = 1e-12;

// ===========================================================================
// Cartesian: deref_quad_points = per-axis affine map [0,1]^DIM -> physical box
// ===========================================================================

TEST (mra_dg_basis_cartesian, line_deref_maps_onto_segment)
{
  using elem = t8_mra::element_data<T8_ECLASS_LINE, 1, 3>;
  t8_mra::dg_basis<elem> dg_basis (3);  // 3 Gauss points

  // index 0 = lower vertex, index 1 (= max for DIM 1) = upper vertex
  double vertices[2][3] = { { 2.0, 0.0, 0.0 }, { 5.0, 0.0, 0.0 } };
  const double lo = 2.0, hi = 5.0;

  const auto phys = dg_basis.deref_quad_points (vertices);
  ASSERT_EQ (phys.size (), dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q) {
    const double ref = dg_basis.quad.points[q];

    EXPECT_NEAR (phys[q], lo + ref * (hi - lo), eps);
    EXPECT_GE (phys[q], lo - eps);
    EXPECT_LE (phys[q], hi + eps);
  }
}

TEST (mra_dg_basis_cartesian, quad_deref_maps_onto_box)
{
  using elem = t8_mra::element_data<T8_ECLASS_QUAD, 1, 3>;
  t8_mra::dg_basis<elem> dg_basis (3);

  // index 0 = lower-left corner, index 2 = upper-right corner
  double vertices[4][3] = { { 1.0, 2.0, 0.0 }, { 4.0, 2.0, 0.0 }, { 4.0, 5.0, 0.0 }, { 1.0, 5.0, 0.0 } };
  const std::array<double, 2> lo { 1.0, 2.0 }, hi { 4.0, 5.0 };

  const auto phys = dg_basis.deref_quad_points (vertices);
  ASSERT_EQ (phys.size (), 2 * dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q)
    for (int d = 0; d < 2; ++d) {
      const double ref = dg_basis.quad.points[2 * q + d];
      EXPECT_NEAR (phys[2 * q + d], lo[d] + ref * (hi[d] - lo[d]), eps);
    }
}

TEST (mra_dg_basis_cartesian, hex_deref_maps_onto_box)
{
  using elem = t8_mra::element_data<T8_ECLASS_HEX, 1, 2>;
  t8_mra::dg_basis<elem> dg_basis (3);

  // index 0 = lower corner, index 7 = upper corner; the rest are unused
  double vertices[8][3] = {};
  const std::array<double, 3> lo { 1.0, 2.0, 3.0 }, hi { 4.0, 6.0, 9.0 };
  for (int d = 0; d < 3; ++d) {
    vertices[0][d] = lo[d];
    vertices[7][d] = hi[d];
  }

  const auto phys = dg_basis.deref_quad_points (vertices);
  ASSERT_EQ (phys.size (), 3 * dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q)
    for (int d = 0; d < 3; ++d) {
      const double ref = dg_basis.quad.points[3 * q + d];
      EXPECT_NEAR (phys[3 * q + d], lo[d] + ref * (hi[d] - lo[d]), eps);
    }
}

// ===========================================================================
// Triangle: barycentric trafo and the reference <-> physical round trip
// ===========================================================================

using tri_elem = t8_mra::element_data<T8_ECLASS_TRIANGLE, 1, 2>;

/* All trafo / ref_point / deref properties for one physical triangle:
 *   - the three vertices map to the unit barycentric coordinates
 *   - a point built from known barycentric weights inverts back to those weights
 *     (covers the centroid and partition of unity / reconstruction)
 *   - deref_quad_points -> ref_point round trips to (1-r0-r1, r0, r1), inside */
void
check_triangle (int rule, const double vertices[3][3])
{
  t8_mra::dg_basis<tri_elem> dg_basis (rule);
  auto [trafo, perm] = dg_basis.trafo_matrix_to_ref_element (vertices);

  // Each vertex maps to a unit barycentric coordinate.
  const std::array<std::array<double, 3>, 3> unit_bary { { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
  for (int v = 0; v < 3; ++v) {
    const auto bary = dg_basis.ref_point (trafo, perm, { vertices[v][0], vertices[v][1] });
    ASSERT_EQ (bary.size (), 3u);
    for (int k = 0; k < 3; ++k)
      EXPECT_NEAR (bary[k], unit_bary[v][k], eps) << "vertex " << v << " coord " << k;
  }

  // A point built from known barycentric weights inverts back to those weights
  // (covers the centroid, partition of unity and reconstruction).
  const std::array<std::array<double, 3>, 3> known_bary {
    { { 0.5, 0.3, 0.2 }, { 0.1, 0.6, 0.3 }, { 1.0 / 3, 1.0 / 3, 1.0 / 3 } }
  };
  for (const auto &expected : known_bary) {
    double px = 0.0, py = 0.0;
    for (int j = 0; j < 3; ++j) {
      px += expected[j] * vertices[j][0];
      py += expected[j] * vertices[j][1];
    }
    const auto bary = dg_basis.ref_point (trafo, perm, { px, py });
    double sum = 0.0;
    for (int k = 0; k < 3; ++k) {
      EXPECT_NEAR (bary[k], expected[k], eps);
      sum += bary[k];
    }
    EXPECT_NEAR (sum, 1.0, eps);
  }

  // Forward map (reference -> physical) and its ref_point inverse round trip.
  const auto phys = dg_basis.deref_quad_points (vertices);
  ASSERT_EQ (phys.size (), 2 * dg_basis.quad.num_points);
  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q) {
    const double r0 = dg_basis.quad.points[2 * q + 0];
    const double r1 = dg_basis.quad.points[2 * q + 1];

    const double phys_x = vertices[0][0] * (1 - r0 - r1) + vertices[1][0] * r0 + vertices[2][0] * r1;
    const double phys_y = vertices[0][1] * (1 - r0 - r1) + vertices[1][1] * r0 + vertices[2][1] * r1;
    EXPECT_NEAR (phys[2 * q + 0], phys_x, eps);
    EXPECT_NEAR (phys[2 * q + 1], phys_y, eps);

    const auto bary = dg_basis.ref_point (trafo, perm, { phys[2 * q + 0], phys[2 * q + 1] });
    EXPECT_NEAR (bary[0], 1 - r0 - r1, eps);
    EXPECT_NEAR (bary[1], r0, eps);
    EXPECT_NEAR (bary[2], r1, eps);
    for (int k = 0; k < 3; ++k)
      EXPECT_GE (bary[k], -eps) << "quad point outside the triangle";
  }
}

/* Run over several triangles: the reference triangle, a generic one, 
 * a reversed-winding one (forces a different LU pivot), and one
 * with negative coordinates. */
TEST (mra_dg_basis_triangle, trafo_and_round_trip_over_many_triangles)
{
  const double reference[3][3] = { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } };
  const double generic[3][3] = { { 0.5, 1.0, 0.0 }, { 2.0, 0.5, 0.0 }, { 1.0, 3.0, 0.0 } };
  const double reversed[3][3] = { { 0.5, 1.0, 0.0 }, { 1.0, 3.0, 0.0 }, { 2.0, 0.5, 0.0 } };
  const double negative[3][3] = { { -1.0, -1.0, 0.0 }, { 1.0, -2.0, 0.0 }, { 0.0, 2.0, 0.0 } };

  for (int rule : { 2, 4 }) {
    check_triangle (rule, reference);
    check_triangle (rule, generic);
    check_triangle (rule, reversed);
    check_triangle (rule, negative);
  }
}

/* A box with negative / shifted coordinates. */
TEST (mra_dg_basis_cartesian, quad_deref_handles_negative_coordinates)
{
  using elem = t8_mra::element_data<T8_ECLASS_QUAD, 1, 3>;
  t8_mra::dg_basis<elem> dg_basis (3);

  double vertices[4][3] = { { -2.0, -1.0, 0.0 }, { 1.0, -1.0, 0.0 }, { 1.0, 0.5, 0.0 }, { -2.0, 0.5, 0.0 } };
  const std::array<double, 2> lo { -2.0, -1.0 }, hi { 1.0, 0.5 };

  const auto phys = dg_basis.deref_quad_points (vertices);
  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q)
    for (int d = 0; d < 2; ++d) {
      const double ref = dg_basis.quad.points[2 * q + d];
      EXPECT_NEAR (phys[2 * q + d], lo[d] + ref * (hi[d] - lo[d]), eps);
    }
}

// ===========================================================================
// basis_value / basis_gradient forward to the reference basis
// ===========================================================================

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

  const auto val = dg_basis.basis_value (x_ref);
  const auto ref_val = basis_t::eval (x);
  for (std::size_t i = 0; i < val.size (); ++i)
    EXPECT_NEAR (val[i], ref_val[i], eps);

  const auto grad = dg_basis.basis_gradient (x_ref);
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
