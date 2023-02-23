/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_cmesh_vtk.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

/* Test if the attribute_offsets of a partitioned cmesh with multiple  
 * attributes is calculated correctly */

/** Construct \a num_trees many cubes each of length 1 connected along the x-axis */
t8_cmesh_t
t8_cmesh_new_row_of_cubes (t8_locidx_t num_trees, const int attributes,
                           sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);
  const t8_geometry_c *linear_geom = t8_geometry_linear_new (3);
  t8_cmesh_register_geometry (cmesh, linear_geom);

  /* Vertices of first cube in row */
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
  };

  for (t8_locidx_t tree_id = 0; tree_id < num_trees; tree_id++) {
    t8_cmesh_set_tree_class (cmesh, tree_id, T8_ECLASS_HEX);
    /* Set first attribut - tree vertices */
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, 8);
    /* Update vertices_coords (x-axis) for next tree */
    for (int v_id = 0; v_id < 8; v_id++) {
      vertices[v_id * 3]++;
    }
    /* Set two more dummy attributs - tree_id & num_trees */
    if (attributes) {
      t8_cmesh_set_attribute
        (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY,
         &tree_id, sizeof (t8_locidx_t), 0);
      t8_cmesh_set_attribute
        (cmesh, tree_id, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY + 1,
         &num_trees, sizeof (t8_locidx_t), 0);
    }
  }

  for (t8_locidx_t tree_id = 0; tree_id < num_trees - 1; tree_id++) {
    t8_cmesh_set_join (cmesh, tree_id, tree_id + 1, 0, 1, 0);
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_partition_cmesh (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh_partition;
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, 0,
                                  t8_scheme_new_default_cxx ());
  t8_cmesh_commit (cmesh_partition, comm);
  return cmesh_partition;
}

/* *INDENT-OFF* */
class cmesh_multiple_attributes : public testing::TestWithParam<int>{
protected:
  void SetUp() override {
    num_trees = GetParam();
    
    cmesh_one_at = t8_cmesh_new_row_of_cubes (num_trees, 0, sc_MPI_COMM_WORLD);
    cmesh_one_at = t8_cmesh_partition_cmesh (cmesh_one_at, sc_MPI_COMM_WORLD);

    cmesh_mult_at = t8_cmesh_new_row_of_cubes (num_trees, 1, sc_MPI_COMM_WORLD);
    cmesh_mult_at = t8_cmesh_partition_cmesh (cmesh_mult_at, sc_MPI_COMM_WORLD);
  }
  void TearDown() override {
    t8_cmesh_destroy (&cmesh_one_at);
    t8_cmesh_destroy (&cmesh_mult_at);
  }

  t8_cmesh_t        cmesh_one_at;
  t8_cmesh_t        cmesh_mult_at;
  t8_locidx_t       num_trees;
};
/* *INDENT-ON* */

/** Returns 1 if two cmeshes have the same vertices.
 * \note The number of trees and the respective eclasses must be the same.
*/
int
t8_cmesh_is_equal_vertices (t8_cmesh_t cmesh_a, t8_cmesh_t cmesh_b)
{
  t8_locidx_t         num_local_trees = cmesh_a->num_local_trees;
  T8_ASSERT (num_local_trees == cmesh_b->num_local_trees);

  for (t8_locidx_t tree_id = 0; tree_id < num_local_trees; tree_id++) {
    double             *vertices_a =
      t8_cmesh_get_tree_vertices (cmesh_a, tree_id);
    double             *vertices_b =
      t8_cmesh_get_tree_vertices (cmesh_b, tree_id);

    t8_eclass_t         eclass = t8_cmesh_get_tree_class (cmesh_a, tree_id);
    T8_ASSERT (eclass == t8_cmesh_get_tree_class (cmesh_b, tree_id));
    const int           num_vertices = t8_eclass_num_vertices[eclass];
    /* Compare vertices of cmesh_a and cmesh_b */
    for (int v_id = 0; v_id < num_vertices * 3; v_id++) {
      if (vertices_a[v_id] != vertices_b[v_id]) {
        return 0;
      }
    }
  }
  return 1;
}

/** Returns 1 if the geometric volume of the trees with there coordinates 
 * is is negative. */
int
t8_cmesh_negative_volume (t8_cmesh_t cmesh)
{
  t8_locidx_t         num_local_trees = cmesh->num_local_trees;
  for (t8_locidx_t tree_id = 0; tree_id < num_local_trees; tree_id++) {
    double             *vertices =
      t8_cmesh_get_tree_vertices (cmesh, tree_id);
    t8_eclass_t         eclass = t8_cmesh_get_tree_class (cmesh, tree_id);
    const int           num_vertices = t8_eclass_num_vertices[eclass];
    /* Check if a tree does have negative volume. */
    if (t8_cmesh_tree_vertices_negative_volume
        (eclass, vertices, num_vertices)) {
      return 1;
    }
  }
  return 0;
}

/* *INDENT-OFF* */
TEST_P (cmesh_multiple_attributes, multiple_attributes) {
    EXPECT_TRUE(t8_cmesh_is_committed (cmesh_one_at));
    EXPECT_TRUE(t8_cmesh_is_committed (cmesh_mult_at));
    EXPECT_TRUE(!t8_cmesh_negative_volume (cmesh_one_at));
    EXPECT_TRUE(t8_cmesh_is_equal_vertices (cmesh_one_at, cmesh_mult_at));
}

/* Test for diffrent number of trees. */
INSTANTIATE_TEST_SUITE_P(t8_gtest_multiple_attributes, cmesh_multiple_attributes, testing::Range(1, 4));
/* *INDENT-ON* */
