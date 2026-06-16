
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
  double verts[2][3] = { { 2.0, 0.0, 0.0 }, { 5.0, 0.0, 0.0 } };
  const double low = 2.0, high = 5.0;

  const auto phys = dg_basis.deref_quad_points (verts);
  ASSERT_EQ (phys.size (), dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q) {
    const double ref = dg_basis.quad.points[q];

    EXPECT_NEAR (phys[q], low + ref * (high - low), eps);
    EXPECT_GE (phys[q], low - eps);
    EXPECT_LE (phys[q], high + eps);
  }
}

TEST (mra_dg_basis_cartesian, quad_deref_maps_onto_box)
{
  using elem = t8_mra::element_data<T8_ECLASS_QUAD, 1, 3>;
  t8_mra::dg_basis<elem> dg_basis (3);

  // index 0 = lower-left corner, index 2 = upper-right corner
  double verts[4][3] = { { 1.0, 2.0, 0.0 }, { 4.0, 2.0, 0.0 }, { 4.0, 5.0, 0.0 }, { 1.0, 5.0, 0.0 } };
  const std::array<double, 2> lo { 1.0, 2.0 }, hi { 4.0, 5.0 };

  const auto phys = dg_basis.deref_quad_points (verts);
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
  double verts[8][3] = {};
  const std::array<double, 3> low { 1.0, 2.0, 3.0 }, high { 4.0, 6.0, 9.0 };
  for (int d = 0; d < 3; ++d) {
    verts[0][d] = low[d];
    verts[7][d] = high[d];
  }

  const auto phys = dg_basis.deref_quad_points (verts);
  ASSERT_EQ (phys.size (), 3 * dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q)
    for (int d = 0; d < 3; ++d) {
      const double ref = dg_basis.quad.points[3 * q + d];
      EXPECT_NEAR (phys[3 * q + d], low[d] + ref * (high[d] - low[d]), eps);
    }
}

// ===========================================================================
// Triangle: barycentric trafo and the reference <-> physical round trip
// ===========================================================================

using tri_elem = t8_mra::element_data<T8_ECLASS_TRIANGLE, 1, 2>;

// A generic (non-degenerate, non-axis-aligned) physical triangle.
constexpr double TRI[3][3] = { { 0.5, 1.0, 0.0 }, { 2.0, 0.5, 0.0 }, { 1.0, 3.0, 0.0 } };

/* The three vertices map to the unit barycentric coordinates. */
TEST (mra_dg_basis_triangle, vertices_map_to_unit_barycentric)
{
  t8_mra::dg_basis<tri_elem> dg_basis (4);
  auto [trafo, perm] = dg_basis.trafo_matrix_to_ref_element (TRI);

  const std::array<std::array<double, 3>, 3> expected { { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };
  for (int v = 0; v < 3; ++v) {
    const auto lambda = dg_basis.ref_point (trafo, perm, { TRI[v][0], TRI[v][1] });
    ASSERT_EQ (lambda.size (), 3u);
    for (int k = 0; k < 3; ++k)
      EXPECT_NEAR (lambda[k], expected[v][k], eps) << "vertex " << v << " coord " << k;
  }
}

/* The centroid maps to (1/3, 1/3, 1/3). */
TEST (mra_dg_basis_triangle, centroid_maps_to_equal_barycentric)
{
  t8_mra::dg_basis<tri_elem> dg_basis (4);
  auto [trafo, perm] = dg_basis.trafo_matrix_to_ref_element (TRI);

  const double cx = (TRI[0][0] + TRI[1][0] + TRI[2][0]) / 3.0;
  const double cy = (TRI[0][1] + TRI[1][1] + TRI[2][1]) / 3.0;

  const auto lambda = dg_basis.ref_point (trafo, perm, { cx, cy });
  for (int k = 0; k < 3; ++k)
    EXPECT_NEAR (lambda[k], 1.0 / 3.0, eps);
}

/* For any point the barycentric coordinates form a partition of unity and
 * reconstruct the point. */
TEST (mra_dg_basis_triangle, barycentric_partition_of_unity_and_reconstruction)
{
  t8_mra::dg_basis<tri_elem> dg_basis (4);
  auto [trafo, perm] = dg_basis.trafo_matrix_to_ref_element (TRI);

  const std::array<std::array<double, 2>, 3> probes { { { 1.0, 1.2 }, { 0.8, 1.5 }, { 1.3, 1.0 } } };
  for (const auto &p : probes) {
    const auto lambda = dg_basis.ref_point (trafo, perm, { p[0], p[1] });

    double sum = 0.0, rx = 0.0, ry = 0.0;
    for (int j = 0; j < 3; ++j) {
      sum += lambda[j];
      rx += lambda[j] * TRI[j][0];
      ry += lambda[j] * TRI[j][1];
    }
    EXPECT_NEAR (sum, 1.0, eps);
    EXPECT_NEAR (rx, p[0], eps);
    EXPECT_NEAR (ry, p[1], eps);
  }
}

/* deref_quad_points maps reference (r0, r1) to the physical triangle, and
 * ref_point inverts it back to the barycentric coords (1-r0-r1, r0, r1). */
TEST (mra_dg_basis_triangle, deref_then_ref_point_round_trip)
{
  t8_mra::dg_basis<tri_elem> dg_basis (4);
  auto [trafo, perm] = dg_basis.trafo_matrix_to_ref_element (TRI);

  const auto phys = dg_basis.deref_quad_points (TRI);
  ASSERT_EQ (phys.size (), 2 * dg_basis.quad.num_points);

  for (std::size_t q = 0; q < dg_basis.quad.num_points; ++q) {
    const double r0 = dg_basis.quad.points[2 * q + 0];
    const double r1 = dg_basis.quad.points[2 * q + 1];

    // Forward map matches the barycentric affine combination.
    const double ex = TRI[0][0] * (1 - r0 - r1) + TRI[1][0] * r0 + TRI[2][0] * r1;
    const double ey = TRI[0][1] * (1 - r0 - r1) + TRI[1][1] * r0 + TRI[2][1] * r1;
    EXPECT_NEAR (phys[2 * q + 0], ex, eps);
    EXPECT_NEAR (phys[2 * q + 1], ey, eps);

    // Inverse recovers the barycentric coordinates and stays inside the cell.
    const auto lambda = dg_basis.ref_point (trafo, perm, { phys[2 * q + 0], phys[2 * q + 1] });
    EXPECT_NEAR (lambda[0], 1 - r0 - r1, eps);
    EXPECT_NEAR (lambda[1], r0, eps);
    EXPECT_NEAR (lambda[2], r1, eps);
    for (int k = 0; k < 3; ++k)
      EXPECT_GE (lambda[k], -eps) << "quad point outside the triangle";
  }
}

