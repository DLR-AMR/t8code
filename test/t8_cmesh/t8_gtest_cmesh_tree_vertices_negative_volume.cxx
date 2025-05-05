/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

/** \file t8_gtest_cmesh_tree_vertices_negative_volume.cxx
* Provide tests to check the functionality of the t8_cmesh_tree_vertices_negative_volume function
* for every eclass.
*/

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <test/t8_gtest_memory_macros.hxx>

/**
 * Given an eclass fill \a vertices_ids with the corner_ids of a cube [0,1]^3, such that 
 * the volume is positive  
 * \param[in] eclass The eclass to use 
 * \param[in, out] 8 ints on input, filled with the corner ids to use for i-th vertex of the element on output 
 */
static void
get_vertices_ids (const t8_eclass_t eclass, int vertices_ids[T8_ECLASS_MAX_CORNERS])
{
  /* For Hex, Quads and Lines we set the remaining vertices by
   * not breaking at the end of the case. 
   * The same is done for Prisms and Triangles. */
  switch (eclass) {
  case T8_ECLASS_HEX:
    vertices_ids[4] = 4;
    vertices_ids[5] = 5;
    vertices_ids[6] = 6;
    vertices_ids[7] = 7;
    [[fallthrough]];
  case T8_ECLASS_QUAD:
    vertices_ids[3] = 3;
    vertices_ids[2] = 2;
    [[fallthrough]];
  case T8_ECLASS_LINE:
    vertices_ids[1] = 1;
    [[fallthrough]];
  case T8_ECLASS_VERTEX:
    vertices_ids[0] = 0;
    break;
  case T8_ECLASS_PRISM:
    vertices_ids[3] = 4;
    vertices_ids[4] = 5;
    vertices_ids[5] = 7;
    [[fallthrough]];
  case T8_ECLASS_TRIANGLE:
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

class tree_vertices_negative_volume: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    tree_class = GetParam ();
    num_vertices = t8_eclass_num_vertices[tree_class];
    get_vertices_ids (tree_class, vertices_ids);
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t tree_class;
  int num_vertices;
  int vertices_ids[T8_ECLASS_MAX_CORNERS];
};

/* Test if positive volume is detected correctly */
TEST_P (tree_vertices_negative_volume, positive_volume)
{
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

  EXPECT_FALSE (t8_cmesh_tree_vertices_negative_volume (tree_class, elem_vertices, num_vertices));
  T8_TESTSUITE_FREE (elem_vertices);
}

/* Test if negative volume is detected correctly */
TEST_P (tree_vertices_negative_volume, negative_volume)
{
  /* clang-format off */
    /* Same nodes as above, but inverted. All 3D elements will have negative volume*/
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
  if (t8_eclass_to_dimension[tree_class] <= 2) {
    EXPECT_FALSE (t8_cmesh_tree_vertices_negative_volume (tree_class, elem_vertices, num_vertices));
  }
  else {
    EXPECT_TRUE (t8_cmesh_tree_vertices_negative_volume (tree_class, elem_vertices, num_vertices));
  }
  T8_TESTSUITE_FREE (elem_vertices);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_tree_vertices_negative_volume, tree_vertices_negative_volume,
                          testing::Range (T8_ECLASS_ZERO, T8_ECLASS_PYRAMID));
