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
#include <t8_cmesh/t8_cmesh_vertex_conn_vertex_to_tree.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

class cmesh_vertex_conn_vtt: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    const int cmesh_id = GetParam ();
    const t8_cmesh_t cmesh = t8_test_create_cmesh (cmesh_id);
    T8_ASSERT (t8_cmesh_is_committed (cmesh));
    const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

    t8_debugf ("Starting test with cmesh of dim %i and %li global, %i local trees.\n", cmesh->dimension,
               t8_cmesh_get_num_trees (cmesh), num_local_trees);

    /* look over all local trees */
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

      const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

      /* loop over all vertices of this tree */
      for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
        /* Set global id of this tree and this vertex to 1 */
        vtt_all_to_one.add_vertex_to_tree (cmesh, 1, itree, ivertex);
        /* We assign a arbitrary but computable global id to this vertex.
         * We comput the id to be (tree_index * vertex_index) mod num_local_trees + 1 */
        const t8_gloidx_t global_id = (itree * ivertex) % (num_local_trees + 1);
        vtt.add_vertex_to_tree (cmesh, global_id, itree, ivertex);
      }
    }
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;

  t8_cmesh_vertex_conn_vertex_to_tree_c vtt_all_to_one; /* all vertices get global id 1 */
  t8_cmesh_vertex_conn_vertex_to_tree_c vtt;            /* multiple global vertex ids */
};

/** Check stored global ids. */
TEST_P (cmesh_vertex_conn_vtt, check_all_to_one)
{
  for (auto& [global_id, tree_vertex_list] : vtt_all_to_one) {
    EXPECT_EQ (global_id, 1);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_vertex_to_tree, cmesh_vertex_conn_vtt, AllCmeshs);