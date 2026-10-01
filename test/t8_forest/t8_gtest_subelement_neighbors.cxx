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
 * Minimal test for face neighbors of subelements.
 */
#include <gtest/gtest.h>
#include <test/t8_gtest_adapt_callbacks.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>

TEST (t8_gtest_subelement_neighbors, face_neighbors_single_tree)
{
  /* A single quad tree, so no face neighbor crosses a tree boundary. */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), 3, 0, sc_MPI_COMM_WORLD);

  /* A single adapt pass from a uniform forest is balanced, so we can resolve the hanging nodes directly. */
  forest = t8_forest_new_adapt (forest, refine_every_nth_element_callback<3>, 0, 0, NULL);
  forest = t8_forest_remove_hanging_nodes (forest);
  ASSERT_TRUE (t8_forest_has_subelements (forest));

  const t8_scheme *scheme = t8_forest_get_scheme (forest);

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
        /* A transition cell has 4 subelements plus one per hanging face. */
        const int siblings = scheme->element_get_num_siblings (tree_class, element);
        EXPECT_TRUE (siblings >= 5 && siblings <= 8);
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
        const std::string where
          = "tree " + std::to_string (itree) + ", leaf " + std::to_string (ielem) + ", face " + std::to_string (iface);

        /* After hanging node resolution every face is matched by exactly one neighbor. */
        EXPECT_EQ (num_neighbors, 1) << where;

        /* The reported element and the leaf stored at the reported index must be the same. */
        t8_locidx_t neigh_tree = -1;
        const t8_element_t *neigh_leaf = t8_forest_get_leaf_element (forest, neigh_indices[0], &neigh_tree);
        EXPECT_TRUE (scheme->element_is_equal (neigh_class, neigh_leaf, neighbors[0])) << where;

        /* The face and the dual face of the neighbor are the same segment, so their centroids coincide. */
        double centroid[3], dual_centroid[3];
        t8_forest_element_face_centroid (forest, itree, element, iface, centroid);
        t8_forest_element_face_centroid (forest, neigh_tree, neigh_leaf, dual_faces[0], dual_centroid);
        EXPECT_LT (std::fabs (centroid[0] - dual_centroid[0]) + std::fabs (centroid[1] - dual_centroid[1]), 1e-10)
          << where;

        /* Crossing back over the dual face must return to the leaf we started from. */
        // int back_num = 0, *back_faces = NULL;
        // t8_locidx_t *back_indices = NULL;
        // const t8_element_t **back_neighbors = NULL;
        // t8_eclass_t back_class;
        // t8_forest_leaf_face_neighbors (forest, neigh_tree, neigh_leaf, &back_neighbors, dual_faces[0], &back_faces,
        //                                &back_num, &back_indices, &back_class);
        // bool found = false;
        // for (int iback = 0; iback < back_num; ++iback) {
        //   found = found || back_indices[iback] == tree_offset + ielem;
        // }
        // EXPECT_TRUE (found) << where << ": no reciprocal neighbor";

        // /* The neighbors are the forest's own leaves, so only the arrays are freed. */
        // T8_FREE (back_neighbors);
        // T8_FREE (back_faces);
        // T8_FREE (back_indices);
        T8_FREE (neighbors);
        T8_FREE (dual_faces);
        T8_FREE (neigh_indices);
      }
    }
  }
  t8_forest_unref (&forest);
}

// TEST (t8_gtest_subelement_geometry, known_face_centroid)
// {
//   t8_subelem_scheme_hanging_nodes_quad scheme;
//   t8_element_t *root;
//   scheme.element_new (1, &root);
//   scheme.set_to_root (root);

//   const int type = 8; /* only f0 (left) is hanging: 5 subelements */
//   const int num_sub = scheme.element_get_num_children (root, type);
//   ASSERT_EQ (num_sub, 5);
//   t8_element_t **cell = T8_ALLOC (t8_element_t *, num_sub);
//   scheme.element_new (num_sub, cell);
//   scheme.element_get_children (root, num_sub, cell, type);

//   /* Subelement 0 is the lower half of the left face, so its three vertices are the centre of the
//    * cell, the lower left corner and the midpoint of the left face. */
//   double v[3][2];
//   for (int ivertex = 0; ivertex < 3; ++ivertex) {
//     scheme.element_get_vertex_reference_coords (cell[0], ivertex, v[ivertex]);
//   }
//   EXPECT_DOUBLE_EQ (v[0][0], 0.5);
//   EXPECT_DOUBLE_EQ (v[0][1], 0.5);

//   scheme.element_destroy (num_sub, cell);
//   T8_FREE (cell);
//   scheme.element_destroy (1, &root);
// }
