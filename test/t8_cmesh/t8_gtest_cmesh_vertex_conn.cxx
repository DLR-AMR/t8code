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
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_types.h>

/** \file In this file we test the global cmesh vertex numbers.
 *
 * We build a test cmesh consisting of two coarse triangles joined together
 * and associate global vertex numbers with the cmesh's vertices.
 * This cmesh has 4 global vertices in total.
 *
 * We then perform three tests
 *
 * 1) check_tree_to_vertex
 * Here we test the tree_to_vertex connectivity.
 * That is, given a tree id, we get a list of the global vertices of that tree
 * (in local vertex order) and check whether this list is correct.
 *
 * 2) check_vertex_to_tree
 * Here we test the vertex_to_tree connectivity.
 * Given a global vertex index, the vertex_to_tree connectivity returns a list
 * of pairs (local tree_id, local_vertex_id) of all the local trees and their local
 * vertices that are connected to the global vertex.
 * We check whether this list is correct.
 *
 * 3) check_global_vertex_number
 * We verify that the number of global vertices is 4.
 * We additionally verify that the process local number of global vertices is 4 as well.
 * This is true, since the cmesh is not partitioned.
 *
 * Additionally, t8_test_cmesh_vertex_conn_partitioned is the start of a test
 * suite with partitioned cmesh that is currently disabled and could be enabled and extended
 * when cmesh vertex connectivity supports partitioned cmeshes.
 * Note that the test itself then has to be set to parallel in the CMake file.
 */

class t8_test_cmesh_vertex_conn: public testing::Test {
 protected:
  void
  SetUp () override
  {
    t8_cmesh_init (&cmesh);

    /*

    We build a replicated cmesh of 2 triangle trees that are connected via one face.
    The triangles are connected via faces 0 and 1 in positive orientation.
    The local vertex numbers look like this:

                 x ----- x
                /2\2   1/
               /   \   /
              /0  1 \0/
             x ----- x

    We manually set global vertex numbers:

                   2       3
                   x ---- x
                  /\     /
                 /  \   /
                /    \ /
             0 x ---- x 1

    */

    const t8_eclass_t tree_class = T8_ECLASS_TRIANGLE;
    /* Set two triangle trees and join them */
    t8_cmesh_set_tree_class (cmesh, 0, tree_class);
    t8_cmesh_set_tree_class (cmesh, 1, tree_class);
    t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
    /* Define and set the global vertices of the trees */
    constexpr t8_gloidx_t global_vertices_of_tree_0[testcase_num_vertices_per_tree] = { 0, 1, 2 };
    constexpr t8_gloidx_t global_vertices_of_tree_1[testcase_num_vertices_per_tree] = { 1, 3, 2 };
    t8_cmesh_set_global_vertices_of_tree (cmesh, 0, global_vertices_of_tree_0, testcase_num_vertices_per_tree);
    t8_cmesh_set_global_vertices_of_tree (cmesh, 1, global_vertices_of_tree_1, testcase_num_vertices_per_tree);
    /* TODO: When cmesh becomes a class implement a cmesh member function to do this instead:
    cmesh.set_global_vertices_of_tree (0, global_vertices_of_tree_0, testcase_num_vertices_per_tree);
    cmesh.set_global_vertices_of_tree (1, global_vertices_of_tree_1, testcase_num_vertices_per_tree);
    */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  /* This test cmesh has 2 trees and 4 global vertices. */
  static constexpr t8_gloidx_t testcase_num_global_trees = 2;
  static constexpr int testcase_num_global_vertices = 4;
  static constexpr int testcase_num_vertices_per_tree = 3;

  t8_cmesh_t cmesh;
};

/** In this test we retrieve the global vertex indices for each of the trees
 * and verify that they are correct.
*/
TEST_F (t8_test_cmesh_vertex_conn, check_tree_to_vertex)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  ASSERT_FALSE (t8_cmesh_is_partitioned (cmesh));

  /* Get the vertices of the trees and check their values. */
  int returned_num_vertices_per_tree;
  const t8_gloidx_t *check_global_vertices_tree_0
    = t8_cmesh_get_global_vertices_of_tree (cmesh, 0, &returned_num_vertices_per_tree);
  EXPECT_EQ (testcase_num_vertices_per_tree, returned_num_vertices_per_tree);
  const t8_gloidx_t *check_global_vertices_tree_1
    = t8_cmesh_get_global_vertices_of_tree (cmesh, 1, &returned_num_vertices_per_tree);
  EXPECT_EQ (testcase_num_vertices_per_tree, returned_num_vertices_per_tree);

  EXPECT_EQ (check_global_vertices_tree_0[0], 0);
  EXPECT_EQ (check_global_vertices_tree_0[1], 1);
  EXPECT_EQ (check_global_vertices_tree_0[2], 2);
  EXPECT_EQ (check_global_vertices_tree_1[0], 1);
  EXPECT_EQ (check_global_vertices_tree_1[1], 3);
  EXPECT_EQ (check_global_vertices_tree_1[2], 2);
}

/** In this test we iterate over all global vertices and for
 * each one check that the list of connected local trees and local tree vertices
 * is correct.
*/
TEST_F (t8_test_cmesh_vertex_conn, check_vertex_to_tree)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  ASSERT_FALSE (t8_cmesh_is_partitioned (cmesh));

  /*

  The vertex to tree lists should have been built up with
  the following entries ( x:y corresponding to tree x, local tree vertex y)

  0: 0:0
  1: 0:1, 1:0
  2: 0:2, 1:2
  3: 1:1

  */
  /* Build a test array to check against.
  * The entry [vertex][tree] corresponds to the local tree vertex
  * of 'tree' that maps to the global vertex 'vertex'.
  * -1 is an invalid value */
  const int check_local_vertices[testcase_num_global_vertices][testcase_num_global_trees] = { {
                                                                                                0, -1 /* Vertex 0 */
                                                                                              },
                                                                                              {
                                                                                                1, 0 /* Vertex 1 */
                                                                                              },
                                                                                              {
                                                                                                2, 2 /* Vertex 2 */
                                                                                              },
                                                                                              {
                                                                                                -1, 1 /* Vertex 3 */
                                                                                              } };
  const int check_num_trees_at_vertex[testcase_num_global_vertices] = { 1, 2, 2, 1 };

  /* Get the vertex to tree list for each vertex. */
  const int num_local_vertices = t8_cmesh_get_num_local_vertices (cmesh);

  for (int ivertex = 0; ivertex < num_local_vertices; ++ivertex) {
    auto &vertex_to_tree_list = t8_cmesh_get_vertex_to_tree_list (cmesh, ivertex);
    /* Check the values via the iterator */
    for (auto &[local_tree, local_vertex] : vertex_to_tree_list) {
      EXPECT_EQ (local_vertex, check_local_vertices[ivertex][local_tree]);
    }

    /* Check that the number of trees at the vertex is correct. */
    const t8_locidx_t num_trees_at_vertex = t8_cmesh_get_num_trees_at_vertex (cmesh, ivertex);
    EXPECT_EQ (num_trees_at_vertex, check_num_trees_at_vertex[ivertex]);
  }
}

/** Check that the number of global/local unique vertices is correct.
 * Since the cmesh is not partitioned, both numbers should be equal to 4.
 */
TEST_F (t8_test_cmesh_vertex_conn, check_global_vertex_number)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  ASSERT_FALSE (t8_cmesh_is_partitioned (cmesh));

  const int num_global_vertices = t8_cmesh_get_num_global_vertices (cmesh);
  const int num_local_vertices = t8_cmesh_get_num_local_vertices (cmesh);
  EXPECT_EQ (num_global_vertices, testcase_num_global_vertices);
  EXPECT_EQ (num_local_vertices, testcase_num_global_vertices);
}

/* Parallel test suite (to be extended) is currently disabled since
* the cmesh vertex connecticity does not support partitioned cmeshes. */
class t8_test_cmesh_vertex_conn_partitioned: public testing::Test {
 protected:
  void
  SetUp () override
  {

    /*

    We build a partitioned cmesh of mpirank many hex trees that are all connected in a long chain.
    Each Hex is connected to the next via face 1 to face 0.
    The last Hex is connected to the first.
    All vertices are mapped to the single global vertex id 42.

    */

    /* We want to have as many trees as mpiranks, so that
     * each rank has at least one tree when partitioned.
     * Compute the mpisize and assign it to the number of trees:
     * */
    const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
    int mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);

    const t8_locidx_t num_trees = mpisize;
    /* Set the trees and join consecutive trees together,
     * also joining the last tree to the first tree, all via faces 1 and 0.
     *
     * Then set all local vertices to point to the same global vertex. */
    ASSERT_EQ (testcase_num_vertices_per_tree, t8_eclass_num_vertices[tree_class]);
    t8_gloidx_t global_vertices_of_tree[testcase_num_vertices_per_tree];
    // Fill array of global vertex ids per tree with unique vertex id.
    std::fill_n (global_vertices_of_tree, testcase_num_vertices_per_tree, global_vertex_id);

    t8_cmesh_init (&cmesh);
    for (t8_gloidx_t itree = 0; itree < num_trees; ++itree) {
      t8_cmesh_set_tree_class (cmesh, itree, tree_class);
      const t8_gloidx_t join_with_tree = (itree + 1) % num_trees;
      const int face_of_this_tree = 1;
      const int face_of_join_tree = 0;
      const int orientation = 0;
      // Join this tree with the next tree
      t8_debugf ("Adding join %li %li [%i %i]\n", static_cast<long>(itree), static_cast<long>(join_with_tree), face_of_this_tree, face_of_join_tree);
      t8_cmesh_set_join (cmesh, itree, join_with_tree, face_of_this_tree, face_of_join_tree, orientation);
      // Set all vertices of this tree to the same single global index.
      t8_cmesh_set_global_vertices_of_tree (cmesh, itree, global_vertices_of_tree, testcase_num_vertices_per_tree);
    }
    // Partition the cmesh such that
    // each process should have one tree and two ghost tree.
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    t8_cmesh_set_partition_range (cmesh, 3, mpirank, mpirank);
    t8_cmesh_commit (cmesh, comm);
    ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), 1) << "Cmesh was not partitioned correctly.";
    if (mpisize > 1) {
      ASSERT_EQ (t8_cmesh_get_num_ghosts (cmesh), 2) << "CMesh was not partitioned correctly.";
    }
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  /* This test cmesh has 1 global vertex. */
  static constexpr t8_eclass_t tree_class = T8_ECLASS_HEX;
  static constexpr int testcase_num_global_vertices = 1;
  static constexpr int testcase_num_vertices_per_tree = 8;
  static constexpr t8_gloidx_t global_vertex_id = 42;
  int mpisize;
  int mpirank;
  t8_cmesh_t cmesh;
};

/** Check that the number of global/local unique vertices is correct.
 * Both numbers should be equal to 1.
 */
TEST_F (t8_test_cmesh_vertex_conn_partitioned, DISABLED_check_global_vertex_number)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  ASSERT_TRUE (t8_cmesh_is_partitioned (cmesh));

  const int num_global_vertices = t8_cmesh_get_num_global_vertices (cmesh);
  const int num_local_vertices = t8_cmesh_get_num_local_vertices (cmesh);
  EXPECT_EQ (num_global_vertices, testcase_num_global_vertices);
  EXPECT_EQ (num_local_vertices, testcase_num_global_vertices);
}
