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

#include <test/t8_gtest_schemes.hxx>

#include <t8_forest/t8_forest_ghost.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_private.h>

struct element_neighbor_eclass: public testing::TestWithParam<std::tuple<int, int>>
{
 protected:
  void
  SetUp () override
  {
    const bool do_ghost = true;
    const int scheme_id = std::get<0> (GetParam ());
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);
    const int level = std::get<1> (GetParam ());

    // Construct a hybrid coarse mesh
    const t8_cmesh_t cmesh = t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);

    // Build a uniform forest
    forest = t8_forest_new_uniform (cmesh, scheme, level, do_ghost, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }

  t8_forest_t forest;
};

TEST_P (element_neighbor_eclass, test_half_neighbors)
{
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  // Iterate over all trees (local and ghosts)
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {

    const t8_eclass_t tree_eclass = t8_forest_get_tree_class (forest, itree);
    const bool is_ghost = itree >= num_local_trees;
    const t8_locidx_t ghost_tree_id = itree - num_local_trees;
    const t8_element_array_t *leaf_elements = is_ghost ? t8_forest_ghost_get_tree_leaf_elements (forest, ghost_tree_id)
                                                       : t8_forest_get_tree_leaf_element_array (forest, itree);
    const t8_locidx_t num_leaves = t8_element_array_get_count (leaf_elements);
    const t8_locidx_t ltreeid_in_cmesh = t8_forest_ltreeid_to_cmesh_ltreeid (forest, itree);

    // Iterate over all leaf elements
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf) {

      const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);

      // Iterate over all faces
      for (int iface = 0; iface < scheme->element_get_num_faces (tree_eclass, element); iface++) {

        // Get the eclass of the face neighbor (function to test)
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, element, iface);

        // Inside the current tree
        if (!scheme->element_is_root_boundary (tree_eclass, element, iface)) {
          // We expect to get the same eclass
          EXPECT_EQ (neigh_class, tree_eclass);
        }
        // Crossing a tree boundary
        else {
          // Get neighbor tree id
          const int tree_face = scheme->element_get_tree_face (tree_eclass, element, iface);
          const t8_locidx_t neighbor_id
            = t8_cmesh_get_face_neighbor (cmesh, ltreeid_in_cmesh, tree_face, nullptr, nullptr);

          if (neighbor_id < 0) {
            // No neighbor
            EXPECT_EQ (neigh_class, T8_ECLASS_INVALID);
          }
          else {
            // Neighbor
            const bool neighbor_is_ghost = t8_cmesh_treeid_is_ghost (cmesh, neighbor_id);
            if (!neighbor_is_ghost) {
              EXPECT_EQ (neigh_class, t8_cmesh_get_tree_class (cmesh, neighbor_id));
            }
            else {
              const t8_locidx_t lghost_neighbor_id = neighbor_id - t8_cmesh_get_num_local_trees (cmesh);
              EXPECT_EQ (neigh_class, t8_cmesh_get_ghost_class (cmesh, lghost_neighbor_id));
            }
          }
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_neighbor_eclass, element_neighbor_eclass,
                          testing::Combine (AllSchemeCollections, testing::Range (0, 4)));
