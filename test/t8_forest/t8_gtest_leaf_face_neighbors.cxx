/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>

#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_adapt_callbacks.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include "t8_test_data_dir.h"
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_adapt_callbacks.hxx>

bool
test_face_neighbors_skip_cmesh (const t8_cmesh_t cmesh)
{
  // We skip empty cmeshes.
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
    const int uniform_level = 0;
    const int adapt_levels = 2;
    const int max_adapt_level = uniform_level + adapt_levels;
    const bool do_ghost = true;
    const bool do_recursive_adapt = true;
    forests[0] = t8_forest_new_uniform (cmesh, scheme, uniform_level, do_ghost, sc_MPI_COMM_WORLD);
    cmesh = t8_forest_get_cmesh (forests[0]);
    t8_forest_ref (forests[0]);
    forests[1] = t8_forest_new_adapt (forests[0], t8_test_adapt_first_child, do_recursive_adapt, do_ghost,
                                      (void *) &max_adapt_level);
    t8_forest_ref (forests[0]);
    // Add another adapted forest that does not create a ghost layer to test
    // the face neighbor algorithm on forests without ghost.
    const bool dont_do_ghost = false;
    forests[2] = t8_forest_new_adapt (forests[0], t8_test_adapt_first_child, do_recursive_adapt, dont_do_ghost,
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

  t8_forest_t forests[3] { nullptr, nullptr, nullptr };
};

/* Check that a leaf/ghost element index matches a given leaf/ghost element.
 * Returns true if the element with index "element_index" (0<= element_index < num_leaves + num_ghosts)
 * lies in the tree "gtreeid" and matches the element "element". */
static void
verify_leaf_element_index (const t8_forest_t forest, const t8_gloidx_t gtreeid, const t8_locidx_t element_index,
                           const t8_element_t *element)
{
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);

  if (element_index < num_local_elements) {
    // The element index belongs to a local leaf element.
    // Check that the  tree is a local tree.
    const t8_locidx_t ltreeid = t8_forest_get_local_id (forest, gtreeid);
    EXPECT_GE (ltreeid, 0) << "Tree of local element is not a local tree.";
    t8_locidx_t check_local_tree_id;
    const t8_element_t *element_from_index = t8_forest_get_leaf_element (forest, element_index, &check_local_tree_id);
    EXPECT_EQ (check_local_tree_id, ltreeid) << "Element index does not match local tree index.";
    EXPECT_EQ (element_from_index, element)
      << "Element at index " << element_index << " does not match the given element.";
  }
  else {
    // The element index belongs to a ghost leaf element.
    // Check that the neighbor tree is a ghost tree.
    const t8_locidx_t ghost_tree_id = t8_forest_ghost_get_ghost_treeid (forest, gtreeid);
    EXPECT_GE (ghost_tree_id, 0) << "Tree of ghost element is not a ghost tree.";
    // There is no get_ghost_leaf_element function that takes only the element index.
    // So we have to convert the element index into a tree local index first.
    // Subtract the number of local elements and the number of ghosts in previous trees.
    const t8_locidx_t ghost_in_tree_index
      = element_index - num_local_elements - t8_forest_ghost_get_tree_element_offset (forest, ghost_tree_id);
    const t8_element_t *ghost_element_from_index
      = t8_forest_ghost_get_leaf_element (forest, ghost_tree_id, ghost_in_tree_index);
    EXPECT_EQ (ghost_element_from_index, element)
      << "Neighbor neighbor ghost element at index " << element_index << " is not original element.";
  }
}

TEST_P (forest_face_neighbors, test_face_neighbors)
{
  /* iterate over all elements */
  bool forest_is_uniform = true;  // The first forest is uniform. We set this to false at the end of the for loop.
  for (auto &forest : forests) {
    const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
    const bool has_ghost = forest->ghosts != NULL;
#if T8_ENABLE_DEBUG
    if (t8_cmesh_get_tree_geometry (cmesh, 0) != NULL) {
      // Debug vtk output, only if cmesh has a registered geometry
      t8_forest_write_vtk (forest, "debug_face_neigh");
      t8_debugf ("writing forest to \'debug_face_neigh\'");
    }
#endif
    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
    const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);
    t8_locidx_t ielement_index = 0;
    for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {
      const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, itree);
      const bool is_ghost = itree >= num_local_trees;
      const t8_locidx_t ghost_tree_id = itree - num_local_trees;
      /* Get the leaf element array */
      const t8_element_array_t *leaf_elements = is_ghost
                                                  ? t8_forest_ghost_get_tree_leaf_elements (forest, ghost_tree_id)
                                                  : t8_forest_get_tree_leaf_element_array (forest, itree);
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
          const t8_element_t **neighbor_leaves;
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
            // the element is a ghost and the neighbor element is neither a local element nor ghost.
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
          // TODO: Use t8_forest_leaf_is_boundary when https://github.com/DLR-AMR/t8code/pull/1081  is integrated
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
              else if (has_ghost) {
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
              else if (has_ghost) {
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
            scheme->element_debug_print (neigh_class, neighbor);
            ASSERT_TRUE (scheme->element_is_valid (neigh_class, neighbor))
              << "Neighbor element " << ineigh << " is not valid";
#endif
            // Check that neighbor index correctly yields neighbor element.
            verify_leaf_element_index (forest, gneigh_tree, neigh_index, neighbor);

            // Compute the local treeid of the neighbor tree.
            const t8_locidx_t neigh_ltreeid
              = neigh_index < num_local_elements
                  ? gneigh_tree - t8_forest_get_first_local_tree_id (forest)
                  : t8_forest_ghost_get_ghost_treeid (forest, gneigh_tree) + num_local_trees;

            // preparation
            const t8_element_t **neigh_neighbor_leaves;
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

#if T8_ENABLE_DEBUG
            t8_debugf ("original element\n");
            scheme->element_debug_print (tree_class, element);
            t8_debugf ("neighbor element\n");
            scheme->element_debug_print (neigh_class, neighbor);
#endif
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
            t8_debugf ("Checking all %i neighbors of neighbor for original element:\n", neigh_num_neighbors);
            for (int ineighneigh = 0; ineighneigh < neigh_num_neighbors && !found_original; ++ineighneigh) {
              // Check that the neighbor of the neighbor element is the original element
              const t8_element_t *neigh_of_neigh = neigh_neighbor_leaves[ineighneigh];
              t8_debugf ("Checking neighbor of neighbor #%i:\n", ineighneigh);
#if T8_ENABLE_DEBUG
              scheme->element_debug_print (tree_class, neigh_of_neigh);
#endif
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

            verify_leaf_element_index (forest, neigh_gneigh_tree, element_index, element);

            // Devnote: We are not using T8_TESTSUITE_FREE here, since the memory is allocated inside of t8_forest_leaf_face_neighbors_ext
            //          by t8code itself using T8_ALLOC. Thus, we need to call T8_FREE to free it.
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
    forest_is_uniform = false;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neighbors, forest_face_neighbors,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);

/** Adapt callback that refines all elements in the tree with global id 0.
 * This callback belong to the test case forest_face_neighbors_two_quad_mesh.
 */
int
t8_test_adapt_first_tree (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                          [[maybe_unused]] t8_locidx_t which_tree, const t8_eclass_t eclass,
                          [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                          [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                          t8_element_t *elements[])
{
  T8_ASSERT (!is_family || (is_family && num_elements == scheme->element_get_num_children (eclass, elements[0])));

  const int level = scheme->element_get_level (eclass, elements[0]);

  /* we set a maximum refinement level as forest user data */
  const int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  const t8_gloidx_t global_tree = t8_forest_global_tree_id (forest_from, which_tree);
  if (global_tree == 0) {
    return 1;
  }
  return 0;
}

/** This test case loads a specific 2D .msh file and iterates over it to compute
 * the leaf face neighbors.
 * It originated from debugging session by Benedict Geihe using Trixi.jl where
 * this setting resulted in errors.
 * This, we implemented it as a test case in order to catch possible errors again in the future.
 * The mesh consists of 2 quad trees that are connected via nontrivial orientation.
 *
 * - Load msh file with 2 quad trees
 * - Adapt mesh such that the first tree is refined once uniformly.
 * - Iterate over all mesh elements and compute the leaf face neighbors.
 *
 *   __ __ _____
 *  |__|__|     |
 *  |__|__|_____|
 *  Tree_0 Tree_1
 *
 *
*/
class forest_face_neighbors_two_quad_mesh: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {

    /* Read our specific mesh file into a cmesh and build a forest. */
    auto meshfile_prefix = t8_test_data_dir / "test_twosquares_twisted";
    const int partition_mesh = 0;
    const int mesh_dim = 2;
    const int main_proc = 0;
    const int use_cad = 0;
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_from_msh_file (&cmesh, meshfile_prefix.c_str (), partition_mesh, sc_MPI_COMM_WORLD, mesh_dim, main_proc,
                            use_cad);
    ASSERT_NE (cmesh, nullptr) << "Could not open mesh file.";

    /* Build a uniform forest on it. */
    const int scheme_id = GetParam ();
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);
    const int level = 0;
    const int do_ghost = 1;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, do_ghost, sc_MPI_COMM_WORLD);

    /* Build an adaptive forest */
    const int do_recursive_adapt = 0;
    const int max_adapt_level = 2;
    adaptive_forest
      = t8_forest_new_adapt (forest, t8_test_adapt_first_tree, do_recursive_adapt, do_ghost, (void *) &max_adapt_level);
  }

  void
  TearDown () override
  {
    if (adaptive_forest != nullptr) {
      t8_forest_unref (&adaptive_forest);
    }
    else if (cmesh != nullptr) {
      t8_cmesh_unref (&cmesh);
    }
  }

  t8_cmesh_t cmesh = nullptr;
  t8_forest_t adaptive_forest = nullptr;
};

/* Perform the actual test for the forest_face_neighbors_two_quad_mesh.
 * Iterate over all leaves and ghosts and call the leaf face neighbor function. */
TEST_P (forest_face_neighbors_two_quad_mesh, check_neighbors)
{
  // For debug purpoese, we write the forest to vtk
  t8_forest_write_vtk (adaptive_forest, "lfn_test_twoquads");

  // Iterate over all elements and compute their neighbors.
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (adaptive_forest);
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (adaptive_forest);
  t8_locidx_t ielement_index = 0;
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {
    const t8_gloidx_t gtreeid = t8_forest_global_tree_id (adaptive_forest, itree);
    const bool is_ghost = itree >= num_local_trees;
    const t8_locidx_t ghost_tree_id = itree - num_local_trees;
    /* Get the leaf element array */
    const t8_element_array_t *leaf_elements
      = is_ghost ? t8_forest_ghost_get_tree_leaf_elements (adaptive_forest, ghost_tree_id)
                 : t8_forest_get_tree_leaf_element_array (adaptive_forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (adaptive_forest, itree);
    const t8_scheme *scheme = t8_forest_get_scheme (adaptive_forest);
    const t8_locidx_t num_leaves = t8_element_array_get_count (leaf_elements);
    int num_faces_with_2_neighbors = 0;
    int num_faces_with_1_neighbor = 0;
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf, ++ielement_index) {
      // Iterate over each leaf element
      const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);
      const int num_faces = scheme->element_get_num_faces (tree_class, element);
      for (int iface = 0; iface < num_faces; ++iface) {
        // Iterate over all faces and compute the face neighbors

        // preparation
        const t8_element_t **neighbor_leaves;
        int *dual_faces;
        int num_neighbors = 0;
        t8_locidx_t *element_indices;
        t8_eclass_t neigh_class;
        t8_gloidx_t gneigh_tree;
        int orientation;

        t8_debugf ("Compute face neighbor for tree %i (%s) element %i (index %i), at face %i.\n", itree,
                   is_ghost ? "ghost" : "local", ileaf, ielement_index, iface);

        // Actual computation of the face neighbors
        t8_forest_leaf_face_neighbors_ext (adaptive_forest, itree, element, &neighbor_leaves, iface, &dual_faces,
                                           &num_neighbors, &element_indices, &neigh_class, &gneigh_tree, &orientation);

        std::string buffer = "";
        for (int ineigh = 0; ineigh < num_neighbors; ++ineigh) {
          buffer += std::to_string (element_indices[ineigh]) + " ";
        }
        t8_debugf ("Tree %i, Element %i, Face %i has %i neighbors:\t%s\n", itree, ileaf, iface, num_neighbors,
                   buffer.c_str ());

        // We are now checking explicit face neighbor values.
        // The number of neighbors must be >= 0
        EXPECT_GE (num_neighbors, 0);
        // Count how many faces have 1 and 2 neighbors.
        // We expect tree 1 to have 1 face with 2 neighbors
        // We cannot state how many faces with 1 neighbor we count for tree 0, since we
        // do not know how the elements are distributed across processes.
        // But we know, that each element has at least one face with 1 neighbor, so the total count must be >0.
        if (num_neighbors == 1) {
          num_faces_with_1_neighbor++;
        }
        else if (num_neighbors == 2) {
          num_faces_with_2_neighbors++;
        }
        if (gtreeid == 0) {
          // Every face in global tree 0 has 0 or 1 face neighbors.
          EXPECT_LE (num_neighbors, 1) << "Tree 0 faces must have 0 or 1 neighbor.";
        }
        else {
          T8_ASSERT (gtreeid == 1);
          // Tree 1 has only one element, it has either 0 neighbors, or 2
          if (num_neighbors > 0) {
            EXPECT_EQ (num_neighbors, 2) << "Tree 2 faces must have 0 or 2 neighbors.";
          }
        }

        // clean-up original element neighbors
        if (num_neighbors > 0) {
          // Devnote: We are not using T8_TESTSUITE_FREE here, since the memory is allocated inside of t8_forest_leaf_face_neighbors_ext
          //          by t8code itself using T8_ALLOC. Thus, we need to call T8_FREE to free it.
          T8_FREE (neighbor_leaves);
          T8_FREE (element_indices);
          T8_FREE (dual_faces);
        }
      }  // End face loop
    }    // End leaf in tree loop
    if (gtreeid == 0) {
      // Tree 0 must have >0 faces with 1 neighbor.
      EXPECT_GE (num_faces_with_1_neighbor, 1) << "Tree 1 must have 1 face with 2 neighbors.";
    }
    if (gtreeid == 1) {
      // Tree 1 must have exactly 1 face with 2 neighbors.
      EXPECT_EQ (num_faces_with_2_neighbors, 1) << "Tree 1 must have 1 face with 2 neighbors.";
    }
  }  // End tree loop
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neighbors, forest_face_neighbors_two_quad_mesh, AllSchemeCollections);

/**
 * \brief Class to test the functionality of \b t8_forest_leaf_neighbor_subface.
 */
struct forest_face_neighbors_subface: public testing::TestWithParam<std::tuple<int, cmesh_example_base *> >
{
 protected:
  /**
   * \brief Set up the test suite.
   */
  void
  SetUp () override
  {
    // Read test parameters.
    const int scheme_id = std::get<0> (GetParam ());
    t8_cmesh_t cmesh = std::get<1> (GetParam ())->cmesh_create ();

    // Skip empty cmeshes.
    if (test_face_neighbors_skip_cmesh (cmesh)) {
      t8_cmesh_unref (&cmesh);
      GTEST_SKIP ();
    }

    // Create uniform base forest from cmesh.
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);
    const bool do_ghost = true;
    const bool do_recursive_adapt = false;
    t8_forest_t base_forest = t8_forest_new_uniform (cmesh, scheme, base_level, do_ghost, sc_MPI_COMM_WORLD);

    // Refine some elements of forest once (and store the resulting forest as member variable).
    forest = t8_forest_new_adapt (base_forest, refine_every_nth_element_callback<3, 1>, do_recursive_adapt, do_ghost,
                                  nullptr);
  }

  /**
   * \brief Tear down the test suite.
   */
  void
  TearDown () override
  {
    // Unref the forest.
    if (forest != nullptr)
      t8_forest_unref (&forest);
  }

  // Member variables.
  t8_forest_t forest = nullptr;  // The forest used within the tests.
  const int base_level = 1;      // The coarse base level in the forest.
};

// Test the functionality of \ref t8_forest_leaf_neighbor_subface.
// The testing is outlined as follows:
// - Loop over all trees, all its elements and all its faces
// - For each face, determine whether it has a one-level coarser neighbor
// - If it does, determine the subface ID via the function t8_forest_leaf_neighbor_subface
// - To validate its output, we ensure that the associated (virtual) child of the neighbor
//   is in fact the corresponding neighbor of the original element.
TEST_P (forest_face_neighbors_subface, test_face_neighbor_subface)
{

  // ----------------------------------------------------------------------------------------------
  // -------------------------------- (0.) Preparation --------------------------------------------
  // ----------------------------------------------------------------------------------------------

  // Explicitly store some forest information for readability.
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements (forest);

  // ----------------------------------------------------------
  // -------- LOOP 1: Iterate over all local trees.
  // ----------------------------------------------------------
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    // Get eclass and number of leaf elements in this tree.
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    t8_locidx_t num_elems_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);

    // ----------------------------------------------------------
    // -------- LOOP 2: Iterate over all local elements in the tree.
    // ----------------------------------------------------------
    for (t8_locidx_t ielem_tree = 0; ielem_tree < num_elems_tree; ielem_tree++) {
      // Get current element and its level.
      const t8_element *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielem_tree);
      const int level = scheme->element_get_level (tree_class, element);

      // Skip if it is on the base level, because it cannot not have a coarser neighbor.
      if (level == base_level)
        continue;

      // Get number of faces.
      const int num_faces = scheme->element_get_num_faces (tree_class, element);

      // ----------------------------------------------------------
      // -------- LOOP 3: Iterate over the element faces
      // ----------------------------------------------------------
      for (int iface = 0; iface < num_faces; iface++) {
        // Compute face neighbor along that face
        // -------------------------------------
        // (i) Define variables
        const t8_element_t **neighbor_leaves;
        int *dual_faces;
        int num_neighbors = 0;
        t8_locidx_t *element_indices;
        t8_eclass_t neigh_class;
        t8_gloidx_t gneigh_tree;
        int orientation;

        // (ii) Actual computation of the face neighbors.
        t8_forest_leaf_face_neighbors_ext (forest, itree, element, &neighbor_leaves, iface, &dual_faces, &num_neighbors,
                                           &element_indices, &neigh_class, &gneigh_tree, &orientation);

        // Only investigate deeper if the face has exactly one neighbor, because otherwise it cannot be a coarser neighbor.
        if (num_neighbors == 1) {

          // Compute level of neighbor.
          const int neighbor_level = scheme->element_get_level (neigh_class, neighbor_leaves[0]);

          // Is the neighbor one level coarser?
          if (neighbor_level == level - 1) {
            // For code readability:
            const t8_element_t *neigh = neighbor_leaves[0];
            const int neigh_face = dual_faces[0];

            // (iii.) Compute subface ID (via the function t8_forest_leaf_neighbor_subface this test is all about).
            int neighbor_subface_id
              = t8_forest_leaf_neighbor_subface (forest, itree, element, iface, neigh_class, neigh, neigh_face);

            // VALIDATION: Iterate over all the neighbor's children and check whether the one associated with neighbor_subface_id
            //             has the original element as face neighbor.
            // --------------------------------------------------------------------
            //
            // (0.) Determine local tree ID of neighbor.
            const t8_locidx_t neigh_ltreeid
              = element_indices[0] < num_local_elements
                  ? gneigh_tree - t8_forest_get_first_local_tree_id (forest)
                  : t8_forest_ghost_get_ghost_treeid (forest, gneigh_tree) + num_local_trees;

            // (1.) Compute the neighbor's children.
            const int num_children = scheme->element_get_num_children (neigh_class, neigh);
            t8_element_t **neigh_children = T8_TESTSUITE_ALLOC (t8_element_t *, num_children);
            scheme->element_new (neigh_class, num_children, neigh_children);
            scheme->element_get_children (neigh_class, neigh, num_children, neigh_children);

            // (2.) Iterate over children and check whether they touch the considered face.
            //      NOTE: This could be done in a simpler way via \ref t8_element_get_children_at_face.
            //            However, the tested function \ref t8_forest_leaf_neighbor_subface itself
            //            already relies on \ref t8_element_get_children_at_face, so it seemed more
            //            reasonable to at least use a slightly different way to get there.
            int check_face = -1;
            int neigh_child_to_check = -1;
            int i_child_at_face = -1;

            for (int ichild = 0; ichild < num_children; ichild++) {

              // (3.) Iterate over all faces of the child.
              const int num_child_faces = scheme->element_get_num_faces (neigh_class, neigh_children[ichild]);
              for (int ichildface = 0; ichildface < num_child_faces; ichildface++) {
                // Get corresponding face of parent (if applicable).
                int parent_face
                  = scheme->element_face_get_parent_face (neigh_class, neigh_children[ichild], ichildface);

                // (4.) If the parent face matches the face we are interested in, we know this child touches the face at its
                // local face iface.
                if (parent_face == neigh_face) {
                  i_child_at_face++;
                  check_face = ichildface;
                  break;
                }
              }
              // (5.) If we have reached the (neighbor_subface_id)-th child touching neigh_face, we can leave the iteration
              // and go on to the checks.
              if (i_child_at_face == neighbor_subface_id) {
                neigh_child_to_check = ichild;
                break;
              }
            }

            // Assertions making sure we actually found the correct child and face to examine.
            ASSERT_GE (check_face, 0);
            ASSERT_GE (neigh_child_to_check, 0);

            // (6.) Compute face neighbor of the neighbor's child.
            t8_element_t *check_elem;
            scheme->element_new (tree_class, 1, &check_elem);
            int dual_dual_face;
            t8_forest_element_face_neighbor (forest, neigh_ltreeid, neigh_children[neigh_child_to_check], check_elem,
                                             tree_class, check_face, &dual_dual_face);

            // Make sure it matches the original element.
            EXPECT_ELEM_EQ (scheme, tree_class, element, check_elem);

            // Make sure the face relation is also correct.
            EXPECT_EQ (iface, dual_dual_face);

            // (7.) Free memory of neigh_children and check_elem
            scheme->element_destroy (neigh_class, num_children, neigh_children);
            scheme->element_destroy (tree_class, 1, &check_elem);
            T8_TESTSUITE_FREE (neigh_children);

          }  // end if(neighbor_level == level - 1)
        }    // end if(num_neighbors == 1)

        // Manually free memory allocated by t8code inside t8_forest_leaf_face_neighbors_ext.
        if (num_neighbors > 0) {
          T8_FREE (neighbor_leaves);
          T8_FREE (element_indices);
          T8_FREE (dual_faces);
        }
      }  // end face loop
    }    // end element loop
  }      // end tree loop
}

// We check for all cmesh examples and scheme collections.
INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neighbors, forest_face_neighbors_subface,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
