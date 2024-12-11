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
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <test/t8_gtest_macros.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

bool
test_face_neighbors_skip_cmesh (const t8_cmesh_t cmesh)
{
  // We only allow cmeshes with pure quad or hex elements.
  // So we check all eclass and if the cmesh contains any of those
  // we skip the cmesh.
  for (int eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {
    if (eclass != T8_ECLASS_QUAD && eclass != T8_ECLASS_HEX) {
      if (cmesh->num_trees_per_eclass[eclass] > 0) {
        return true;
      }
    }
  }
  // Additionally, we skip empty cmeshes.
  return t8_cmesh_is_empty (cmesh);
}

class forest_face_neighbors: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    t8_cmesh_t cmesh = GetParam ()->cmesh_create ();
    if (test_face_neighbors_skip_cmesh (cmesh)) {
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
  const t8_forest_t forest_uniform = forest;
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest_uniform);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest_uniform);
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest_uniform);
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements (forest_uniform);
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {
    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_uniform, itree);
    const bool is_ghost = itree >= num_local_trees;
    const t8_locidx_t ghost_tree_id = itree - num_local_trees;
    /* Get the leaf element array */
    const t8_element_array_t *leaf_elements = !is_ghost
                                                ? t8_forest_get_tree_element_array (forest_uniform, itree)
                                                : t8_forest_ghost_get_tree_elements (forest_uniform, ghost_tree_id);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest_uniform, itree);
    const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest_uniform, tree_class);
    const t8_locidx_t num_leafs = t8_element_array_get_count (leaf_elements);
    const t8_locidx_t cmesh_tree = t8_forest_ltreeid_to_cmesh_ltreeid (forest_uniform, itree);
    for (t8_locidx_t ileaf = 0; ileaf < num_leafs; ++ileaf) {
      // Iterate over each leaf element
      const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);
      const int num_faces = scheme->t8_element_num_faces (element);
      for (int iface = 0; iface < num_faces; ++iface) {
        // Iterate over all faces and compute the face neighbors

        // preparation
        t8_element_t **neighbor_leaves;
        int *dual_faces;
        int num_neighbors = 0;
        t8_locidx_t *element_indices;
        t8_eclass_scheme_c *neigh_scheme;
        t8_gloidx_t gneigh_tree;
        int orientation;
        // Actual computation of the face neighbors
        t8_forest_leaf_face_neighbors_ext (forest_uniform, itree, element, &neighbor_leaves, iface, &dual_faces,
                                           &num_neighbors, &element_indices, &neigh_scheme, &gneigh_tree, &orientation);

        t8_debugf ("Tree %i element %i at face %i has %i face neighbors.\n", itree, ileaf, iface, num_neighbors);

        if (gneigh_tree < 0) {
          // No neighbors are found, check for correctly set return values
          ASSERT_EQ (num_neighbors, 0);
          ASSERT_EQ (neighbor_leaves, nullptr);
          ASSERT_EQ (element_indices, nullptr);
          ASSERT_EQ (dual_faces, nullptr);
        }
        else {
          ASSERT_GE (num_neighbors, 0);
          ASSERT_NE (neighbor_leaves, nullptr);
          ASSERT_NE (element_indices, nullptr);
          ASSERT_NE (dual_faces, nullptr);
        }
        // Checking for:
        //      uniform forest:
        //        - inner element has 1 face neighbors
        //        - boundary element has 0 face neighbors // TODO: Add this test for all forests
        //        - If E face f has neighbor E' face f', then
        //             E' face f' must have neighbor E face f.

        // Now checking for inner and boundary elements.

        // Compute whether this element is a boundary element or not.
        // An element is a boundary element if it lies on the tree boundary
        // and if the corresponding tree face is at the domain boundary.
        const bool is_root_boundary = scheme->t8_element_is_root_boundary (element, iface);
        const int tree_face = scheme->t8_element_tree_face (element, iface);
        const bool is_boundary_element
          = is_root_boundary && t8_cmesh_tree_face_is_boundary (cmesh, cmesh_tree, tree_face);

        if (!is_boundary_element) {
          EXPECT_EQ (num_neighbors, 1) << "Inner element should have exactly 1 neighbor, has " << num_neighbors << ".";
        }
        else {
          EXPECT_EQ (num_neighbors, 0) << "Boundary element should have exactly 0 neighbors, has " << num_neighbors
                                       << ".";
          ;
        }

        // Check that the neighbor of the neighbor is the original element.
        for (int ineigh = 0; ineigh < num_neighbors; ++ineigh) {
          const t8_element_t *neighbor = neighbor_leaves[ineigh];
          const int dual_face = dual_faces[ineigh];
          const t8_locidx_t neigh_index = element_indices[ineigh];

          t8_debugf ("Checking neighbor element %p.\n", neighbor);
          t8_debugf ("dual face is %i, index is %i\n", dual_face, neigh_index);
          ASSERT_TRUE (neigh_scheme->t8_element_is_valid (neighbor))
            << "Neighbor element " << ineigh << " is not valid";
          // Compute the local tree id of the neighbors tree depending on whether
          // it is a local tree or a ghost tree.
          const t8_locidx_t neigh_ltreeid
            = neigh_index < num_local_elements
                ? gneigh_tree - t8_forest_get_first_local_tree_id (forest_uniform)
                : t8_forest_ghost_get_ghost_treeid (forest_uniform, gneigh_tree) + num_local_trees;
          // preparation
          t8_element_t **neigh_neighbor_leaves;
          int *neigh_dual_faces;
          int neigh_num_neighbors = 0;
          t8_locidx_t *neigh_element_indices;
          t8_eclass_scheme_c *neigh_neigh_scheme;
          t8_gloidx_t neigh_gneigh_tree;
          int neigh_orientation;
          // Actual computation of the neighbor's face neighbors
          t8_forest_leaf_face_neighbors_ext (forest_uniform, neigh_ltreeid, neighbor, &neigh_neighbor_leaves, dual_face,
                                             &neigh_dual_faces, &neigh_num_neighbors, &neigh_element_indices,
                                             &neigh_neigh_scheme, &neigh_gneigh_tree, &neigh_orientation);

          // We must have found at least one face neighbor, namely the original element.
          EXPECT_GE (neigh_num_neighbors, 1);
          // The neighbor's neighbor tree must be the current tree
          EXPECT_EQ (gtree_id, neigh_gneigh_tree);
          // The neighbor's scheme must be the current scheme
          EXPECT_EQ (scheme, neigh_neigh_scheme);
          // The neighbor's orientation must be the orientation
          EXPECT_EQ (orientation, neigh_orientation);

          // TODO Check for the specific face neighbor that would match the original element

          // clean-up neighbor's neighbors
          if (neigh_num_neighbors > 0) {
            T8_FREE (neigh_neighbor_leaves);
            T8_FREE (neigh_element_indices);
            T8_FREE (neigh_dual_faces);
          }
        }

        // clean-up original element neighbors
        if (num_neighbors > 0) {
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
