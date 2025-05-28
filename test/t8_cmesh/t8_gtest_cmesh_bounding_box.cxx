/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/**
 * \file t8_gtest_cmesh_bounding_box.cxx
 * 
 * Test the computation of the bounding box of a t8_cmesh.
 */

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

/**
 * \brief Test fixture for testing the bounding box of a t8_cmesh. Computes a cmesh inside
 * the unit cube, computes the bounding box and checks that it is correct.
 */
class t8_cmesh_bounding_box: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_eclass eclass;
};

static void
compute_and_check_bounds (const t8_cmesh_t cmesh, const t8_eclass eclass)
{
  double bounds[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounds);

  const int dim = t8_eclass_to_dimension[eclass];

  for (int idim = 0; idim < dim; ++idim) {
    EXPECT_DOUBLE_EQ (bounds[idim * 2], 0);
    EXPECT_DOUBLE_EQ (bounds[idim * 2 + 1], 1)<< "Dimension " << idim
                                                              << " should have the same bounds, but got "
                                                              << bounds[idim * 2] << " and " << bounds[idim * 2 + 1];
  }
  for (int idim = dim; idim < 3; ++idim) {
    EXPECT_DOUBLE_EQ (bounds[idim * 2], bounds[idim * 2 + 1]);
  }
}

TEST_P (t8_cmesh_bounding_box, test_box)
{
  compute_and_check_bounds (cmesh, eclass);
}

class t8_cmesh_bounding_box_multi_trees: public testing::TestWithParam<std::tuple<int, bool, t8_eclass>> {
 protected:
  void
  SetUp () override
  {
    trees_per_dim = std::get<0> (GetParam ());
    axis_aligned_geometry = std::get<1> (GetParam ());
    eclass = std::get<2> (GetParam ());
    if (axis_aligned_geometry && (eclass != T8_ECLASS_HEX && eclass != T8_ECLASS_QUAD)) {
      GTEST_SKIP () << "Axis-aligned geometry is only supported for hex and quad elements.";
    }
  }

  int trees_per_dim;
  bool axis_aligned_geometry;
  t8_eclass eclass;
};

TEST_P (t8_cmesh_bounding_box_multi_trees, hypercube)
{
  t8_debugf ("Created cmesh with %d trees per dimension and eclass %s.\n", trees_per_dim,
                t8_eclass_to_string[eclass]);
  const double cube_bounds[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (eclass, sc_MPI_COMM_WORLD, cube_bounds, trees_per_dim, trees_per_dim,
                                                 trees_per_dim, axis_aligned_geometry);
  compute_and_check_bounds (cmesh, eclass);
  t8_cmesh_unref (&cmesh);

}

// Parameterized test cases
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_bounding_box, t8_cmesh_bounding_box, AllEclasses, print_eclass);

auto print_hypercube_test = [] (const testing::TestParamInfo<std::tuple<int, bool, t8_eclass>> &info) -> std::string
{
  int trees_per_dim = std::get<0> (info.param);
  bool axis_aligned_geometry = std::get<1> (info.param);
  t8_eclass_t eclass = std::get<2> (info.param);
  const testing::TestParamInfo<t8_eclass> dummy_eclass_info (eclass, -1);
  std::string eclass_str = print_eclass (dummy_eclass_info);
  return std::to_string (trees_per_dim) + "_trees_per_dim_" + 
         (axis_aligned_geometry ? "_linear_axis_aligned_" : "_linear_")  + eclass_str;
};

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_bounding_box, t8_cmesh_bounding_box_multi_trees,
                          testing::Combine (testing::Range (2, 5),  // trees_per_dim
                                            testing::Bool (),       // axis_aligned_geometry
                                            testing::Range (T8_ECLASS_ZERO, T8_ECLASS_PYRAMID)), print_hypercube_test);
