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

#include <gtest/gtest.h>
#include <t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.hxx>
#include <t8_cmesh_generator/t8_cmesh_example_sets.hxx>

/* In this file we test the t8_cmesh_vertex_conn_vertex_to_tree
 * class of the cmesh global vertex list.
 * We iterate over all cmeshes and for each case we
 * construct two global id lists.
 * 1. A single global vertex that maps to all local vertices.
 * 2. Multiple global vertices in a non geometric/semi-random pattern.
 *
 * We add the information to the list and then check whether
 * this information is maintained with the getter functions. */

static int
t8_compute_tree_num_vertices (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Get the trees class depending on whether it is a local tree or ghost. */
  const t8_eclass_t tree_class = ltreeid < num_local_trees
                                   ? t8_cmesh_get_tree_class (cmesh, ltreeid)
                                   : t8_cmesh_get_ghost_class (cmesh, ltreeid - num_local_trees);
  return t8_eclass_num_vertices[tree_class];
}

/* Given a tree id and a vertex of that tree compute a (pseudo random)
 * hash value. We use this value to assign a global vertex id.
 * This is just a number that we thought of for testing purpose. */
static t8_locidx_t
t8_compute_global_vertex_hash (t8_locidx_t itree, t8_locidx_t ivertex, t8_locidx_t num_local_trees)
{
  return (itree * ivertex) % (num_local_trees + 1);
}

class t8_test_cmesh_vertex_conn_vtt: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh = GetParam ()->cmesh_create ();
    T8_ASSERT (t8_cmesh_is_committed (cmesh));
    const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
    const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);

    t8_debugf ("Starting test with cmesh of dim %i and %" T8_GLOIDX_FORMAT " global, %i local trees.\n",
               cmesh->dimension, t8_cmesh_get_num_trees (cmesh), num_local_trees);

    /* look over all local trees */
    for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; ++itree) {
      const int num_tree_vertices = t8_compute_tree_num_vertices (cmesh, itree);

      /* loop over all vertices of this tree */
      for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
        /* Set global id of this tree and this vertex to 1 */
        t8_debugf ("Adding vertex %" T8_GLOIDX_FORMAT " to tree %i v %i\n", global_vertex_id, itree, ivertex);
        vtt_all_to_one.add_vertex_to_tree (cmesh, global_vertex_id, itree, ivertex);
        /* We assign a arbitrary but computable global id to this vertex.
         * We compute the id to be (tree_index * vertex_index) mod num_local_trees + 1 */
        const t8_gloidx_t global_id = t8_compute_global_vertex_hash (itree, ivertex, num_local_trees);
        vtt.add_vertex_to_tree (cmesh, global_id, itree, ivertex);
      }
    }
    /* We added all vertices, so we commit. */
    vtt.commit (cmesh);
    vtt_all_to_one.commit (cmesh);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;

  /* The single global vertex that everything is assigned to for
   * vtt_all_to_one. */
  const t8_gloidx_t global_vertex_id = 1;

  t8_cmesh_vertex_conn_vertex_to_tree vtt_all_to_one; /* all vertices get global id 1 */
  t8_cmesh_vertex_conn_vertex_to_tree vtt;            /* multiple global vertex ids */
};

/** Check stored global ids for the case with a single global vertex. */
TEST_P (t8_test_cmesh_vertex_conn_vtt, check_all_to_one)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);

  /* Count the number of entries (for each local tree its number of vertices). */
  size_t num_entries = 0;
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; ++itree) {
    const int num_tree_vertices = t8_compute_tree_num_vertices (cmesh, itree);

    num_entries += num_tree_vertices;
  }

  /* We can only check for the values, if we have a non-zero number of vertices.
   * Depending on the cmesh there do not have to be any vertices. */
  if (num_entries > 0) {
    /* Get the list of the single vertex. */
    auto &tree_list = vtt_all_to_one.get_tree_list_of_vertex (global_vertex_id);

    /* Asserting that the number of entries matches the computed number.
    * If it does not, we cannot continue, therefore ASSERT instead of EXPECT. */
    ASSERT_EQ (tree_list.size (), num_entries);

    /* Iterate over all entries of the tree list.
    * We expect that this single list stores all local trees, ghost trees and
    * their vertices sorted by id.
    */
    t8_locidx_t check_local_tree = 0;
    int check_vertex = 0;
    for (auto &[tree_id, tree_vertex] : tree_list) {
      T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, check_local_tree)
                 || t8_cmesh_treeid_is_ghost (cmesh, check_local_tree));

      /* Check that the current entry matches the local tree and tree vertex. */
      EXPECT_EQ (tree_id, check_local_tree);
      EXPECT_EQ (tree_vertex, check_vertex);
      /* increase check vertex */
      ++check_vertex;

      const int num_tree_vertices = t8_compute_tree_num_vertices (cmesh, check_local_tree);
      /* If we reached the end of this tree's vertices, we reset
      * the vertex counter and advance the tree counter. */
      if (check_vertex >= num_tree_vertices) {
        check_vertex = 0;
        check_local_tree++;
      }
    }
  }
}

/* Check stored global ids for the case with multiple global ids. */
TEST_P (t8_test_cmesh_vertex_conn_vtt, check_multiple_ids)
{
  /* We need to check that each local tree/ghost and each vertex
   * exists exactly once in the list.
   * We do so by setting up an indicator array storing the
   * number of vertices for each tree and count down for each occurrence.
   * At the end the values must be zero. */

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);
  const t8_locidx_t num_trees_and_ghosts = num_local_trees + num_ghost_trees;

  std::vector<int> vertex_counts (num_trees_and_ghosts);
  /* Fill each entry with the number of vertices. */
  for (t8_locidx_t itree = 0; itree < num_trees_and_ghosts; ++itree) {
    /* Compute number of vertices of this tree. */
    const int num_tree_vertices = t8_compute_tree_num_vertices (cmesh, itree);

    /* Set the number of vertices to the entry. */
    vertex_counts[itree] = num_tree_vertices;
  }

  /* Iterate over all entries in vtt.
   * Each entry corresponds to a global vertex id and
   * gives its list of tree indices and vertices. */
  for (auto &[global_vertex, tree_vertex_list] : vtt) {
    /* Iterate over the list of tree indices and vertices of this global vertex. */
    for (auto &[tree_index, tree_vertex] : tree_vertex_list) {
      const t8_locidx_t num_vertices = t8_compute_tree_num_vertices (cmesh, tree_index);
      const t8_locidx_t hash_value = t8_compute_global_vertex_hash (tree_index, tree_vertex, num_local_trees);

      /* Check that the global index matches our computed hash value. */
      EXPECT_EQ (global_vertex, hash_value);
      /* 0 <= tree_vertex < num_vertices */
      EXPECT_GE (tree_vertex, 0);
      EXPECT_LT (tree_vertex, num_vertices);

      /* remove this tree_vertex from the vertex_count */
      vertex_counts[tree_index]--;
      /* Count must be >= 0 */
      EXPECT_GE (vertex_counts[tree_index], 0);
    }
  }
  /* After we have checked all global ids, we must have matched all tree
   * vertices exactly once. If and only if this is the case, vertex_counts[i] is 0
   * for each tree. */
  for (auto &vertex_count : vertex_counts) {
    EXPECT_EQ (vertex_count, 0);
  }
}

/* TODO: Enable this test as soon as we can add attribute to
 *       derived cmeshes. */
/* Check stored global ids for the case with multiple global ids. */
TEST_P (t8_test_cmesh_vertex_conn_vtt, DISABLED_convert_to_ttv_and_back)
{
  t8_cmesh_t derived_cmesh_A, derived_cmesh_B;
  t8_cmesh_init (&derived_cmesh_A);
  t8_cmesh_init (&derived_cmesh_B);
  /* The original cmesh must survive this test to be destroyed during TearDown and
   * to be used in other tests.
   * Hence we need to ref it twice, once for each new cmesh. */
  t8_cmesh_ref (cmesh);
  t8_cmesh_ref (cmesh);
  t8_cmesh_set_derive (derived_cmesh_A, cmesh);
  t8_cmesh_set_derive (derived_cmesh_B, cmesh);
  /* Construct ttv connectivities from the two vtt connectivities. */
  t8_cmesh_vertex_conn_tree_to_vertex ttv (cmesh, derived_cmesh_A, vtt);
  t8_cmesh_vertex_conn_tree_to_vertex ttv_all_to_one (cmesh, derived_cmesh_B, vtt_all_to_one);
  /* Commit the cmeshes to actually build the ttv connectivities. */
  t8_cmesh_commit (derived_cmesh_A, sc_MPI_COMM_WORLD);
  t8_cmesh_commit (derived_cmesh_B, sc_MPI_COMM_WORLD);

  /* Now we can build vtt conns from the ttv conns.
   * They should match the original connectivities. */
  t8_cmesh_vertex_conn_vertex_to_tree vtt_new;
  vtt_new.build_from_ttv (derived_cmesh_A, ttv);
  t8_cmesh_vertex_conn_vertex_to_tree vtt_new_all_to_one;
  vtt_all_to_one.build_from_ttv (derived_cmesh_B, ttv_all_to_one);

  /* Check for equality. */
  EXPECT_EQ (vtt, vtt_new);
  EXPECT_EQ (vtt_all_to_one, vtt_new_all_to_one);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_vertex_to_tree, t8_test_cmesh_vertex_conn_vtt, AllCmeshsParam,
                          pretty_print_base_example);
