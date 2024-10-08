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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* Test if multiple attributes are partitioned correctly. */

/** Return a partitioned cmesh from \a cmesh. */
static t8_cmesh_t
t8_cmesh_partition_cmesh (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh_partition;
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, 0, t8_scheme_new_default_cxx ());
  t8_cmesh_commit (cmesh_partition, comm);
  return cmesh_partition;
}

class cmesh_multiple_attributes: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    num_trees = GetParam ();

    cmesh_one_at = t8_cmesh_new_row_of_cubes (num_trees, 0, 0, sc_MPI_COMM_WORLD);
    cmesh_one_at = t8_cmesh_partition_cmesh (cmesh_one_at, sc_MPI_COMM_WORLD);

    cmesh_mult_at = t8_cmesh_new_row_of_cubes (num_trees, 1, 0, sc_MPI_COMM_WORLD);
    cmesh_mult_at = t8_cmesh_partition_cmesh (cmesh_mult_at, sc_MPI_COMM_WORLD);

    cmesh_mult_at_from_stash = t8_cmesh_new_row_of_cubes (num_trees, 1, 1, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh_one_at);
    t8_cmesh_destroy (&cmesh_mult_at);
    t8_cmesh_destroy (&cmesh_mult_at_from_stash);
  }

  t8_cmesh_t cmesh_one_at;
  t8_cmesh_t cmesh_mult_at;
  t8_cmesh_t cmesh_mult_at_from_stash;
  t8_locidx_t num_trees;
};

/** Check attribute values of cmeshes against reference values. */
TEST_P (cmesh_multiple_attributes, multiple_attributes)
{
  /* Vertices of first cube in row as reference. */
  const double vertices_ref[24] = {
    0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1,
  };

  /* Check partitioned cmesh with one attribute. */
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh_one_at));
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh_one_at);
  for (t8_locidx_t ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    t8_gloidx_t gtree_id = t8_cmesh_get_global_id (cmesh_one_at, ltree_id);
    const double *vertices_partition = t8_cmesh_get_tree_vertices (cmesh_one_at, ltree_id);

    EXPECT_EQ (T8_ECLASS_HEX, t8_cmesh_get_tree_class (cmesh_one_at, ltree_id));
    /* Compare vertices with reference vertices. */
    for (int v_id = 0; v_id < 8; v_id++) {
      EXPECT_EQ (vertices_partition[v_id * 3], vertices_ref[v_id * 3] + gtree_id);
      EXPECT_EQ (vertices_partition[v_id * 3 + 1], vertices_ref[v_id * 3 + 1]);
      EXPECT_EQ (vertices_partition[v_id * 3 + 2], vertices_ref[v_id * 3 + 2]);
    }
  }

  /* Check partitioned cmesh with three attributes. */
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh_mult_at));
  EXPECT_EQ (num_local_trees, t8_cmesh_get_num_local_trees (cmesh_mult_at));
  for (t8_locidx_t ltree_id = 0; ltree_id < num_local_trees; ltree_id++) {
    t8_gloidx_t gtree_id = t8_cmesh_get_global_id (cmesh_mult_at, ltree_id);
    const double *vertices_partition = t8_cmesh_get_tree_vertices (cmesh_one_at, ltree_id);

    EXPECT_EQ (T8_ECLASS_HEX, t8_cmesh_get_tree_class (cmesh_one_at, ltree_id));

    /* Compare vertices with reference vertices. */
    for (int v_id = 0; v_id < 8; v_id++) {
      EXPECT_EQ (vertices_partition[v_id * 3], vertices_ref[v_id * 3] + gtree_id);
      EXPECT_EQ (vertices_partition[v_id * 3 + 1], vertices_ref[v_id * 3 + 1]);
      EXPECT_EQ (vertices_partition[v_id * 3 + 2], vertices_ref[v_id * 3 + 2]);
    }
    /* Compare second attribute with global tree id. */
    t8_locidx_t att;
    att = *(t8_locidx_t *) t8_cmesh_get_attribute (cmesh_mult_at, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY,
                                                   ltree_id);
    EXPECT_EQ (gtree_id, att);
    /* Compare third attribute with global number of trees. */
    att = *(t8_locidx_t *) t8_cmesh_get_attribute (cmesh_mult_at, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY + 1,
                                                   ltree_id);
    EXPECT_EQ (att, t8_cmesh_get_num_trees (cmesh_mult_at));
  }
  /* Check partitioned cmesh from stash with three attributes. */
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh_mult_at_from_stash));
  EXPECT_EQ (num_local_trees, t8_cmesh_get_num_local_trees (cmesh_mult_at_from_stash));
  t8_locidx_t num_ghosts = t8_cmesh_get_num_ghosts (cmesh_mult_at_from_stash);
  for (t8_locidx_t ltree_id = 0; ltree_id < num_local_trees + num_ghosts; ltree_id++) {
    const t8_gloidx_t gtree_id = t8_cmesh_get_global_id (cmesh_mult_at_from_stash, ltree_id);
    const double *vertices_partition = t8_cmesh_get_tree_vertices (cmesh_mult_at_from_stash, ltree_id);
    const t8_eclass_t eclass = (ltree_id < num_local_trees)
                                 ? t8_cmesh_get_tree_class (cmesh_one_at, ltree_id)
                                 : t8_cmesh_get_ghost_class (cmesh_one_at, ltree_id - num_local_trees);
    EXPECT_EQ (T8_ECLASS_HEX, eclass);

    /* Compare vertices with reference vertices. */
    for (int v_id = 0; v_id < 8; v_id++) {
      EXPECT_EQ (vertices_partition[v_id * 3], vertices_ref[v_id * 3] + gtree_id)
        << " at tree id " << ltree_id << " and vertex " << v_id;
      EXPECT_EQ (vertices_partition[v_id * 3 + 1], vertices_ref[v_id * 3 + 1])
        << " at tree id " << ltree_id << " and rtex " << v_id;
      EXPECT_EQ (vertices_partition[v_id * 3 + 2], vertices_ref[v_id * 3 + 2])
        << " at tree id " << ltree_id << " and vertex " << v_id;
    }
    /* Compare second attribute with global tree id. */
    t8_locidx_t att;
    att = *(t8_locidx_t *) t8_cmesh_get_attribute (cmesh_mult_at_from_stash, t8_get_package_id (),
                                                   T8_CMESH_NEXT_POSSIBLE_KEY, ltree_id);
    EXPECT_EQ (gtree_id, att);
    /* Compare third attribute with global number of trees. */
    att = *(t8_locidx_t *) t8_cmesh_get_attribute (cmesh_mult_at_from_stash, t8_get_package_id (),
                                                   T8_CMESH_NEXT_POSSIBLE_KEY + 1, ltree_id);
    EXPECT_EQ (att, t8_cmesh_get_num_trees (cmesh_mult_at_from_stash));
  }
}

/* Test for different number of trees. */
INSTANTIATE_TEST_SUITE_P (t8_gtest_multiple_attributes, cmesh_multiple_attributes, testing::Range (1, 10));
