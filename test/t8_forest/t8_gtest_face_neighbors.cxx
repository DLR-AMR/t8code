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
#include <test/t8_gtest_adapt_callbacks.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_schemes.hxx>

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

class forest_face_neighbors: public testing::TestWithParam<std::tuple<int, cmesh_example_base *> > {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    t8_cmesh_t cmesh = std::get<1> (GetParam ())->cmesh_create ();
    if (test_face_neighbors_skip_cmesh (cmesh)) {
      /* we skip empty cmeshes case */
      t8_cmesh_unref (&cmesh);
      GTEST_SKIP ();
    }
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);
    const int level = 1;
    const int adapt_levels = 2;
    const int max_adapt_level = level + adapt_levels;
    const bool do_ghost = true;
    const bool do_recursive_adapt = true;
    forests[0] = t8_forest_new_uniform (cmesh, scheme, level, do_ghost, sc_MPI_COMM_WORLD);
    cmesh = t8_forest_get_cmesh (forests[0]);
    t8_forest_ref (forests[0]);
    forests[1] = t8_forest_new_adapt (forests[0], t8_test_adapt_first_child, do_recursive_adapt, do_ghost,
                                      (void *) &max_adapt_level);
  }

  void
  TearDown () override
  {
    for (auto &forest : forests) {
      if (forest != nullptr) {
        t8_forest_unref (&forest);
      }
    }
  }

  t8_forest_t forests[2] { nullptr, nullptr };
};

TEST_P (forest_face_neighbors, test_face_neighbors)
{
  /* iterate over all elements */
  bool forest_is_uniform = true;  // The first forest is uniform. We set this to false at the end of the for loop.
  for (auto &forest : forests) {
    const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
#if T8_ENABLE_DEBUG
    if (t8_cmesh_get_tree_geometry (cmesh, 0) != NULL) {
      // Debug vtk output, only if cmesh has a registered geometry
      t8_forest_write_vtk (forest, "debug_face_neigh");
    }
#endif
    if (!forest_is_uniform) {
      // Currently, adaptive forest is not working properly due to a bug in the
      // face neighbor computation.
      // We hence currently skip the test for adaptive forests.
      // break;
    }
    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
    const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements (forest);
    t8_locidx_t ielement_index = 0;
    for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {
      const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, itree);
      const bool is_ghost = itree >= num_local_trees;
      const t8_locidx_t ghost_tree_id = itree - num_local_trees;
      /* Get the leaf element array */
      const t8_element_array_t *leaf_elements = !is_ghost ? t8_forest_get_tree_element_array (forest, itree)
                                                          : t8_forest_ghost_get_tree_elements (forest, ghost_tree_id);
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
      const t8_scheme *scheme = t8_forest_get_scheme (forest);
      const t8_locidx_t num_leaves = t8_element_array_get_count (leaf_elements);
      const t8_locidx_t cmesh_tree = t8_forest_ltreeid_to_cmesh_ltreeid (forest, itree);
      for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf, ++ielement_index) {
        // Iterate over each leaf element
        const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);
        const int num_faces = scheme->element_get_num_faces (tree_class, element);
        for (int iface = 0; iface < num_faces; ++iface) {
          // Iterate over all faces and compute the face neighbors

          // preparation
          t8_element_t **neighbor_leaves;
          int *dual_faces;
          int num_neighbors = 0;
          t8_locidx_t *element_indices;
          t8_eclass_t neigh_class;
          t8_gloidx_t gneigh_tree;
          int orientation;

          t8_debugf ("Compute face neighbor for tree %i (%s) element %i (index %i), at face %i.\n", itree,
                     is_ghost ? "ghost" : "local", ileaf, ielement_index, iface);

          // Actual computation of the face neighbors
          t8_forest_leaf_face_neighbors_ext (forest, itree, element, &neighbor_leaves, iface, &dual_faces,
                                             &num_neighbors, &element_indices, &neigh_class, &gneigh_tree,
                                             &orientation);

          t8_debugf ("Tree %i element %i at face %i has %i face neighbors.\n", itree, ileaf, iface, num_neighbors);

          if (gneigh_tree < 0) {
            // If there is no neighbor tree then there cannot be any face neighbors.
            // Note that there can also be no face neighbors computed if a neighbor tree exists, but
            // the element is a ghost and the neighbor would is neither a local element nor ghost.
            ASSERT_EQ (num_neighbors, 0);
          }
          if (num_neighbors == 0) {
            // No neighbors are found, check for correctly set return values
            ASSERT_TRUE (element_indices == NULL);
            ASSERT_TRUE (neighbor_leaves == NULL);
            ASSERT_TRUE (dual_faces == NULL);
          }
          else {
            ASSERT_GE (num_neighbors, 0);
            ASSERT_TRUE (neighbor_leaves != NULL);
            ASSERT_TRUE (element_indices != NULL);
            ASSERT_TRUE (dual_faces != NULL);
          }

          // Checking for:
          //      uniform and adapted forest:
          //        - inner local element has >= 1 face neighbors (= 1 for uniform)
          //        - inner ghost element has 0 or 1 face neighbors (= 1 for uniform)
          //        - boundary element has 0 face neighbors
          //        - If E face f has neighbor E' face f', then
          //             E' face f' must have neighbor E face f.

          // Now checking for inner and boundary elements.

          // Compute whether this element is a boundary element or not.
          // An element is a boundary element if it lies on the tree boundary
          // and if the corresponding tree face is at the domain boundary.
          const bool is_root_boundary = scheme->element_is_root_boundary (tree_class, element, iface);
          const int tree_face = scheme->element_get_tree_face (tree_class, element, iface);
          const bool is_boundary_element
            = is_root_boundary && t8_cmesh_tree_face_is_boundary (cmesh, cmesh_tree, tree_face);

          if (!is_boundary_element) {
            if (!is_ghost) {  // Local element
              if (forest_is_uniform) {
                // In a uniform forest we must have exactly 1 neighbor.
                EXPECT_EQ (num_neighbors, 1)
                  << "Inner local element should have exactly 1 neighbor, has " << num_neighbors << ".";
              }
              else {
                // In an adaptive forest we have 1 or more neighbors.
                EXPECT_GE (num_neighbors, 1)
                  << "Inner local element should have at least 1 neighbor, has " << num_neighbors << ".";
              }
            }
            else {  // Ghost element
              if (forest_is_uniform) {
                // In a uniform forest a ghost element has none or one neighbor.
                EXPECT_TRUE (num_neighbors == 0 || num_neighbors == 1)
                  << "Inner ghost element should have exactly 1 or 0 neighbors, has " << num_neighbors << ".";
              }
              else {
                // In an adaptive forest a ghost element has 0 or more neighbors.
                EXPECT_GE (num_neighbors, 0)
                  << "Inner ghost element should have 0 or more neighbors, has " << num_neighbors << ".";
              }
            }
          }
          else {
            EXPECT_EQ (num_neighbors, 0) << "Boundary element should have exactly 0 neighbors, has " << num_neighbors
                                         << ".";
          }

          if (forest_is_uniform) {
            ASSERT_TRUE (num_neighbors == 0 || num_neighbors == 1);
            // Check the index computation function and that it computes the correct neighbor index.
            int check_dual_face;
            const t8_locidx_t check_same_level_index = t8_forest_same_level_leaf_face_neighbor_index (
              forest, ielement_index, iface, gtree_id, &check_dual_face);

            if (check_dual_face < 0) {
              EXPECT_EQ (num_neighbors, 0);
            }
            if (check_dual_face >= 0) {
              EXPECT_EQ (dual_faces[0], check_dual_face);
              EXPECT_EQ (element_indices[0], check_same_level_index);
            }
          }

          // Check that the neighbor of the neighbor is the original element.
          for (int ineigh = 0; ineigh < num_neighbors; ++ineigh) {
            const t8_element_t *neighbor = neighbor_leaves[ineigh];
            const int dual_face = dual_faces[ineigh];
            const t8_locidx_t neigh_index = element_indices[ineigh];

            t8_debugf ("Checking neighbor element %p in (global) tree %li.\n", (void *) neighbor, gneigh_tree);
            t8_debugf ("dual face is %i, index is %i\n", dual_face, neigh_index);

#if T8_ENABLE_DEBUG
            ASSERT_TRUE (scheme->element_is_valid (neigh_class, neighbor))
              << "Neighbor element " << ineigh << " is not valid";
#endif
            t8_locidx_t neigh_ltreeid_from_index;
            // Check that neighbor index correctly yields neighbor element.
            if (neigh_index < num_local_elements) {
              const t8_element_t *neighbor_from_index
                = t8_forest_get_element (forest, neigh_index, &neigh_ltreeid_from_index);
              EXPECT_TRUE (scheme->element_is_equal (neigh_class, neighbor_from_index, neighbor));
            }
            // TODO: Check neighbor index if the element is a ghost element

            // Compute the local tree id of the neighbors tree depending on whether
            // it is a local tree or a ghost tree.
            const t8_locidx_t neigh_ltreeid
              = neigh_index < num_local_elements
                  ? gneigh_tree - t8_forest_get_first_local_tree_id (forest)
                  : t8_forest_ghost_get_ghost_treeid (forest, gneigh_tree) + num_local_trees;
            if (neigh_index < num_local_elements) {
              EXPECT_EQ (neigh_ltreeid, neigh_ltreeid_from_index);
            }  // TODO: Check neighbor ltreeid if ghost tree
            // preparation
            t8_element_t **neigh_neighbor_leaves;
            int *neigh_dual_faces;
            int neigh_num_neighbors = 0;
            t8_locidx_t *neigh_element_indices;
            t8_eclass_t neigh_neigh_class;
            t8_gloidx_t neigh_gneigh_tree;
            int neigh_orientation;
            // Actual computation of the neighbor's face neighbors
            t8_forest_leaf_face_neighbors_ext (forest, neigh_ltreeid, neighbor, &neigh_neighbor_leaves, dual_face,
                                               &neigh_dual_faces, &neigh_num_neighbors, &neigh_element_indices,
                                               &neigh_neigh_class, &neigh_gneigh_tree, &neigh_orientation);

            t8_debugf ("original element\n");
            scheme->element_debug_print (neigh_class, element);
            t8_debugf ("neighbor element\n");
            scheme->element_debug_print (neigh_class, neighbor);
            fflush (stdout);
            // We must have found at least one face neighbor, namely the original element.
            EXPECT_GE (neigh_num_neighbors, 1);
            // The neighbor's neighbor tree must be the current tree
            EXPECT_EQ (gtree_id, neigh_gneigh_tree);
            // The neighbor's scheme must be the current scheme
            EXPECT_EQ (tree_class, neigh_neigh_class);
            // The neighbor's orientation must be the orientation
            EXPECT_EQ (orientation, neigh_orientation);

            // We now (try to) find the original element among the neighbors.
            // If it does not exist there was an error.
            // If it exists we check that dual face and index were computed correctly.

            int position_of_original_element = -1;
            bool found_original = false;
            for (int ineighneigh = 0; ineighneigh < neigh_num_neighbors && !found_original; ++ineighneigh) {
              // Check that the neighbor of the neighbor element is the original element
              const t8_element_t *neigh_of_neigh = neigh_neighbor_leaves[ineighneigh];
              if (scheme->element_is_equal (tree_class, element, neigh_of_neigh)) {
                position_of_original_element = ineighneigh;
                found_original = true;  // Stop the for loop
              }
            }
            // We must have found the original element among the neighbors.
            ASSERT_TRUE (found_original) << "The original element was not a neighbor of its neighbor.";

            // Check that the dual face of the dual face is the original face
            const int neigh_dual_face = neigh_dual_faces[position_of_original_element];
            EXPECT_EQ (neigh_dual_face, iface);

            // Check that the index is correct, i.e. when getting the neighbor neighbor element from the index
            // we retrieve the original element.
            const t8_locidx_t element_index = neigh_element_indices[position_of_original_element];
            EXPECT_GE (element_index, 0);

            if (element_index < num_local_elements) {
              const t8_element_t *element_from_index = t8_forest_get_element (forest, element_index, NULL);
              EXPECT_EQ (element_from_index, element)
                << "Neighbor neighbor element at index " << element_index << " is not original element.";
            }
            // TODO: Check element index if original element is a ghost element

            // clean-up neighbor's neighbors
            if (neigh_num_neighbors > 0) {
              scheme->element_destroy (neigh_class, neigh_num_neighbors, neigh_neighbor_leaves);
              T8_FREE (neigh_neighbor_leaves);
              T8_FREE (neigh_element_indices);
              T8_FREE (neigh_dual_faces);
            }
          }

          // clean-up original element neighbors
          if (num_neighbors > 0) {
            scheme->element_destroy (neigh_class, num_neighbors, neighbor_leaves);
            T8_FREE (neighbor_leaves);
            T8_FREE (element_indices);
            T8_FREE (dual_faces);
          }
        }
      }
    }
    forest_is_uniform = false;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neighbors, forest_face_neighbors,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
