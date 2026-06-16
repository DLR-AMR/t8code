
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

