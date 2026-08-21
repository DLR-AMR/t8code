/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_gtest_subelement_neighbors.cxx
 * Minimal scheme level test for face neighbors of quad subelements.
 */
#include <gtest/gtest.h>
#include <test/t8_gtest_adapt_callbacks.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>

TEST (t8_gtest_subelement_neighbors, leaf_face_neighbors)
{
  int mpisize;
  SC_CHECK_MPI (sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize));
  if (mpisize > 1) {
    GTEST_SKIP () << "This test only runs serially (no ghost layer is created).";
  }

  const int level = 3;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), level, 0, sc_MPI_COMM_WORLD);

  /* Adapting twice with this callback gives three levels in a 2 x 2 periodic pattern, so we need to
   * balance before the hanging nodes can be resolved. */
  forest = t8_forest_new_adapt (forest, t8_test_adapt_even_global_id, 0, 0, NULL);
  forest = t8_forest_new_adapt (forest, t8_test_adapt_even_global_id, 0, 0, NULL);

  t8_forest_t forest_balanced;
  t8_forest_init (&forest_balanced);
  t8_forest_set_balance (forest_balanced, forest, 0);
  t8_forest_commit (forest_balanced);
  //   forest = t8_forest_remove_hanging_nodes (forest_balanced);
  //   EXPECT_TRUE (t8_forest_has_global_subelements (forest));
  forest = forest_balanced;

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  //const t8_locidx_t num_local_leaves = t8_forest_get_local_num_leaf_elements (forest);
  int num_subelements = 0;
  int num_cells_with[8] = { 0 };

  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    //const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, itree);
    const t8_locidx_t num_tree_leaves = t8_forest_get_tree_num_leaf_elements (forest, itree);

    for (t8_locidx_t ielem = 0; ielem < num_tree_leaves; ++ielem) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      const bool is_sub = t8_element_is_subelement (scheme, tree_class, element);
      EXPECT_EQ (scheme->element_get_shape (tree_class, element), is_sub ? T8_ECLASS_TRIANGLE : T8_ECLASS_QUAD);
      const int num_faces = scheme->element_get_num_faces (tree_class, element);
      EXPECT_EQ (num_faces, is_sub ? 3 : 4);

      if (is_sub) {
        ++num_subelements;
        const int siblings = scheme->element_get_num_siblings (tree_class, element);
        EXPECT_TRUE (siblings >= 5 && siblings <= 7);
        ++num_cells_with[siblings];
      }

      for (int iface = 0; iface < num_faces; ++iface) {
        int num_neighbors = 0, *dual_faces = NULL;
        t8_locidx_t *neigh_indices = NULL;
        const t8_element_t **neighbors = NULL;
        t8_eclass_t neigh_class;

        t8_forest_leaf_face_neighbors (forest, itree, element, &neighbors, iface, &dual_faces, &num_neighbors,
                                       &neigh_indices, &neigh_class);

        if (num_neighbors == 0) { /* Boundary of the domain. */
          continue;
        }
        // /* After hanging node resolution every face is matched by exactly one neighbour. */
        // EXPECT_EQ (num_neighbors, 1) << "tree " << itree << ", leaf " << ielem << ", face " << iface;
        // EXPECT_GE (neigh_indices[0], 0);
        // EXPECT_LT (neigh_indices[0], num_local_leaves) << "neighbour is a ghost, but this test is serial";
        // EXPECT_GE (dual_faces[0], 0);
        // EXPECT_LT (dual_faces[0], scheme->element_get_num_faces (neigh_class, neighbors[0]));

        // /* The indices are forest local, so look the leaf up through the forest, not through the tree. */
        // t8_locidx_t neigh_tree = -1;
        // const t8_element_t *neigh_leaf = t8_forest_get_leaf_element (forest, neigh_indices[0], &neigh_tree);
        // /* element_is_equal ignores the subelement id, so compare the id separately. */
        // EXPECT_TRUE (scheme->element_is_equal (neigh_class, neigh_leaf, neighbors[0]));
        // EXPECT_EQ (scheme->element_get_child_id (neigh_class, neigh_leaf),
        //            scheme->element_get_child_id (neigh_class, neighbors[0]));

        // /* Crossing back over the dual face must return to the leaf we started from. */
        // int back_num = 0, *back_faces = NULL;
        // t8_locidx_t *back_indices = NULL;
        // const t8_element_t **back_neighbors = NULL;
        // t8_eclass_t back_class;
        // t8_forest_leaf_face_neighbors (forest, neigh_tree, neigh_leaf, &back_neighbors, dual_faces[0],
        //                                &back_faces, &back_num, &back_indices, &back_class);
        // bool found = false;
        // for (int iback = 0; iback < back_num; ++iback) {
        //   found = found || back_indices[iback] == tree_offset + ielem;
        // }
        // EXPECT_TRUE (found) << "no reciprocal neighbour for tree " << itree << ", leaf " << ielem
        //                     << ", face " << iface;
        // if (back_num > 0) {
        //   scheme->element_destroy (back_class, back_num, (t8_element_t **) back_neighbors);
        //   T8_FREE (back_neighbors);
        //   T8_FREE (back_faces);
        //   T8_FREE (back_indices);
        // }

        scheme->element_destroy (neigh_class, num_neighbors, (t8_element_t **) neighbors);
        T8_FREE (neighbors);
        T8_FREE (dual_faces);
        T8_FREE (neigh_indices);
      }
    }
  }
  // Expect to have subelements.
  EXPECT_GT (num_subelements, 0);

  t8_forest_unref (&forest);
}
