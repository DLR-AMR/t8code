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
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* In this file we test TODO: document */

class t8_test_cmesh_vertex_conn: public testing::Test {
 protected:
  void
  SetUp () override
  {
    t8_cmesh_init (&cmesh);

    /*

    We build a replicated cmesh of 2 triangle trees that are coinnected via one face.
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

/** */
TEST_F (t8_test_cmesh_vertex_conn, check_tree_to_vertex)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  ASSERT_FALSE (t8_cmesh_is_partitioned (cmesh));

  /* Get the vertices of the trees and check their values. */
  const t8_gloidx_t *check_global_vertices_tree_0
    = t8_cmesh_get_global_vertices_of_tree (cmesh, 0, testcase_num_vertices_per_tree);
  const t8_gloidx_t *check_global_vertices_tree_1
    = t8_cmesh_get_global_vertices_of_tree (cmesh, 0, testcase_num_vertices_per_tree);
  EXPECT_EQ (check_global_vertices_tree_0[0], 0);
  EXPECT_EQ (check_global_vertices_tree_0[1], 1);
  EXPECT_EQ (check_global_vertices_tree_0[2], 2);
  EXPECT_EQ (check_global_vertices_tree_1[0], 1);
  EXPECT_EQ (check_global_vertices_tree_1[1], 3);
  EXPECT_EQ (check_global_vertices_tree_1[2], 2);
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

/** */
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