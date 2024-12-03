/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <test/t8_gtest_macros.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

class forest_face_neighbors: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    t8_cmesh_t cmesh = GetParam ()->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* we skip empty cmeshes case */
      t8_cmesh_unref (&cmesh);
      GTEST_SKIP ();
    }
    t8_scheme_cxx_t *default_scheme = t8_scheme_new_default_cxx ();
    const int level = 1;
    const bool do_ghost = true;
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, do_ghost, sc_MPI_COMM_WORLD);
    cmesh = t8_forest_get_cmesh (forest);
  }

  void
  TearDown () override
  {
    if (forest != nullptr) {
      t8_forest_unref (&forest);
    }
  }

  t8_forest_t forest { nullptr };
};

TEST_P (forest_face_neighbors, test_face_neighbors)
{
  /* iterate over all elements */
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {
    const bool is_ghost = itree >= num_local_trees;
    const t8_locidx_t ghost_tree_id = itree - num_local_trees;
    /* Get the leaf element array */
    const t8_element_array_t *leaf_elements = !is_ghost ? t8_forest_get_tree_element_array (forest, itree)
                                                        : t8_forest_ghost_get_tree_elements (forest, ghost_tree_id);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
    const t8_locidx_t num_leafs = t8_element_array_get_count (leaf_elements);
    for (t8_locidx_t ileaf = 0; ileaf < num_leafs; ++ileaf) {
      // Iterate over each leaf element
      const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);
      const int num_faces = scheme->t8_element_num_faces (element);
      for (int iface = 0; iface < num_faces; ++iface) {
        // Iterate over all faces and compute the face neighbors

        // preparation
        t8_element_t **neighbor_leaves;
        int *dual_faces;
        int num_neighbors;
        t8_locidx_t *element_indices;
        t8_eclass_scheme_c *neigh_scheme;
        t8_gloidx_t gneigh_tree;
        int orientation;
        // Actual computation
        t8_forest_leaf_face_neighbors_ext (forest, itree, element, &neighbor_leaves, iface, &dual_faces, &num_neighbors,
                                           &element_indices, &neigh_scheme, &gneigh_tree, &orientation);

        // clean-up
        if (num_neighbors > 0) {
          scheme->t8_element_destroy (num_neighbors, neighbor_leaves);
          T8_FREE (neighbor_leaves);
          T8_FREE (element_indices);
          T8_FREE (dual_faces);
        }
      }
    }
  }
}

#if 0
//// OLD CODE

    for (t8_locidx_t ielement = 0; ielement < t8_forest_get_tree_num_elements (forest, itree); ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      /* iterate over the faces */
      for (int face = 0; face < ts->t8_element_num_faces (element); face++) {
        /* Get the eclass of the face neighbor and get the scheme */
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, element, face);
        t8_eclass_scheme_c *neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
        const int num_face_neighs = ts->t8_element_num_face_children (element, face);
        t8_element_t **half_neighbors = T8_ALLOC (t8_element_t *, num_face_neighs);
        ts->t8_element_new (num_face_neighs, half_neighbors);
        t8_forest_element_half_face_neighbors (forest, itree, element, half_neighbors, neigh_scheme, face,
                                               num_face_neighs, NULL);
        /* allocate memory for element's neighbor and construct it */
        neigh_scheme->t8_element_new (1, &neighbor);
        const t8_locidx_t neigh_tree
          = t8_forest_element_face_neighbor (forest, itree, element, neighbor, neigh_scheme, face, &dual_face);
        if (neigh_tree > 0) {
          /* We now check whether the face children of neighbor are the half neighbors. */
          T8_ASSERT (num_face_neighs == neigh_scheme->t8_element_num_face_children (neighbor, dual_face));
          t8_element_t **neighbor_face_children = T8_ALLOC (t8_element_t *, num_face_neighs);
          neigh_scheme->t8_element_new (num_face_neighs, neighbor_face_children);
          int *child_ids = T8_ALLOC (int, num_face_neighs);
          neigh_scheme->t8_element_children_at_face (neighbor, dual_face, neighbor_face_children, num_face_neighs,
                                                     child_ids);
          /* Check that the children at face of the neighbor are the half neighbors of the element */
          for (int ineigh = 0; ineigh < num_face_neighs; ineigh++) {
            EXPECT_ELEM_EQ (neigh_scheme, neighbor_face_children[ineigh], half_neighbors[ineigh])
              << "ineigh = " << ineigh << " face = " << face;
          }
          neigh_scheme->t8_element_destroy (num_face_neighs, neighbor_face_children);
          T8_FREE (child_ids);
          T8_FREE (neighbor_face_children);
        }
        neigh_scheme->t8_element_destroy (1, &neighbor);
        neigh_scheme->t8_element_destroy (num_face_neighs, half_neighbors);
        T8_FREE (half_neighbors);
      }
    }
  }
  t8_forest_unref (&forest);
  sc_array_reset (&owners);
}
#endif

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neighbors, forest_face_neighbors, AllCmeshsParam);
