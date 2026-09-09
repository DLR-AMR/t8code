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

/** \file t8_gtest_geometry_negative_volume.cxx
* Provide tests to check the functionality of the negative volume check function
* of the geometries.
*/

#include <gtest/gtest.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <test/t8_gtest_memory_macros.hxx>
#include <test/t8_gtest_macros.hxx>

struct check_negative_volume: public testing::TestWithParam<t8_eclass_t>
{
 protected:
  void
  SetUp () override
  {
    tree_class = GetParam ();
    num_vertices = t8_eclass_num_vertices[tree_class];
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t tree_class;
  int num_vertices;
};

/**
 * Given an eclass fill \a vertices_ids with the corner_ids of a cube [0,1]^3, such that
 * the volume is positive
 * \param [in] eclass The eclass to use
 * \param [in, out] vertices_ids 8 ints on input, filled with the corner ids to use for i-th vertex of the element on output
 */
static void
get_vertices_ids (const t8_eclass_t eclass, int vertices_ids[T8_ECLASS_MAX_CORNERS])
{
  switch (eclass) {
  case T8_ECLASS_HEX:
    vertices_ids[4] = 4;
    vertices_ids[5] = 5;
    vertices_ids[6] = 6;
    vertices_ids[7] = 7;
    vertices_ids[3] = 3;
    vertices_ids[2] = 2;
    vertices_ids[1] = 1;
    vertices_ids[0] = 0;
    break;
  case T8_ECLASS_PRISM:
    vertices_ids[3] = 4;
    vertices_ids[4] = 5;
    vertices_ids[5] = 7;
    vertices_ids[0] = 0;
    vertices_ids[1] = 1;
    vertices_ids[2] = 3;
    break;
  case T8_ECLASS_TET:
    vertices_ids[0] = 0;
    vertices_ids[1] = 1;
    vertices_ids[2] = 5;
    vertices_ids[3] = 7;
    break;
  case T8_ECLASS_PYRAMID:
    vertices_ids[0] = 1;
    vertices_ids[1] = 3;
    vertices_ids[2] = 0;
    vertices_ids[3] = 2;
    vertices_ids[4] = 7;
    break;
  default:
    break;
  }
}

/* Test if positive volume of 3D cells is detected correctly */
TEST_P (check_negative_volume, linear_geometry_positive_volume)
{
  int vertices_ids[T8_ECLASS_MAX_CORNERS];
  get_vertices_ids (tree_class, vertices_ids);
  /* clang-format off */
  const double vertices_coords[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  /* clang-format on */
  double *elem_vertices = T8_TESTSUITE_ALLOC (double, 3 * num_vertices);
  t8_cmesh_new_translate_vertices_to_attributes (vertices_ids, vertices_coords, elem_vertices, num_vertices);
  t8_cmesh_t cmesh = t8_cmesh_new ();
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, tree_class);
  t8_cmesh_set_tree_vertices (cmesh, 0, elem_vertices, num_vertices);
  /* Disable negative volume check so that we can check the results in the test. */
  t8_cmesh_disable_negative_volume_check (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  EXPECT_TRUE (t8_cmesh_validate_geometry (cmesh, true));
  T8_TESTSUITE_FREE (elem_vertices);
  t8_cmesh_destroy (&cmesh);
}

/* Test if negative volume of 3D cells is detected correctly */
TEST_P (check_negative_volume, linear_geometry_negative_volume)
{
  int vertices_ids[T8_ECLASS_MAX_CORNERS];
  get_vertices_ids (tree_class, vertices_ids);
  /* Same nodes as above, but inverted. All 3D elements will have negative volume*/
  /* clang-format off */
  const double vertices_coords[24] = {
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0
  };
  /* clang-format on */
  double *elem_vertices = T8_TESTSUITE_ALLOC (double, 3 * num_vertices);
  t8_cmesh_new_translate_vertices_to_attributes (vertices_ids, vertices_coords, elem_vertices, num_vertices);
  t8_cmesh_t cmesh = t8_cmesh_new ();
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, tree_class);
  t8_cmesh_set_tree_vertices (cmesh, 0, elem_vertices, num_vertices);
  /* Disable negative volume check so that we can check the results in the test. */
  t8_cmesh_disable_negative_volume_check (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  if (t8_eclass_to_dimension[tree_class] <= 2) {
    EXPECT_TRUE (t8_cmesh_validate_geometry (cmesh, true));
  }
  else {
    EXPECT_FALSE (t8_cmesh_validate_geometry (cmesh, true));
  }
  T8_TESTSUITE_FREE (elem_vertices);
  t8_cmesh_destroy (&cmesh);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry_negative_volume, check_negative_volume,
                          testing::Range (T8_ECLASS_HEX, T8_ECLASS_COUNT), print_eclass);
