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

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>

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
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_full_hybrid (cmesh, sc_MPI_COMM_WORLD);

    // Partition the cmesh
    t8_cmesh_t cmesh_partition;
    t8_cmesh_init (&cmesh_partition);
    // The cmesh takes ownership of the scheme, but we still pass it to the forest below.
    scheme->ref ();
    t8_cmesh_set_partition_uniform (cmesh_partition, level, scheme);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);

    // Build a uniform forest
    forest = t8_forest_new_uniform (cmesh_partition, scheme, level, do_ghost, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }

  t8_forest_t forest;
};

TEST_P (element_neighbor_eclass, neighbor_eclass_matches_neighbor_leaves)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  // Iterate over all trees (local and ghosts)
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees + num_ghost_trees; itree++) {

    const t8_eclass_t tree_eclass = t8_forest_get_tree_class (forest, itree);
    const bool is_ghost = itree >= num_local_trees;
    const t8_locidx_t ghost_tree_id = itree - num_local_trees;
    const t8_element_array_t *leaf_elements
      = is_ghost ? t8_forest_ghost_get_tree_leaf_elements (forest, ghost_tree_id)
                 : t8_forest_tree_get_leaf_elements (forest, itree);
    const t8_locidx_t num_leaves = t8_element_array_get_count (leaf_elements);

    // Iterate over all leaf elements
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf) {

      const t8_element_t *element = t8_element_array_index_locidx (leaf_elements, ileaf);
      const int num_faces = scheme->element_get_num_faces (tree_eclass, element);

      // Iterate over all faces
      for (int iface = 0; iface < num_faces; iface++) {

        // Get the eclass of the face neighbor (function to test)
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, element, iface);

        // Compute the actual face neighbor leaves
        const t8_element_t **neighbors = nullptr;
        int *dual_faces = nullptr;
        int num_neighbors = 0;
        t8_locidx_t *element_indices = nullptr;
        t8_eclass_t dummy_eclass;
        t8_gloidx_t gneigh_tree = -1;
        t8_forest_leaf_face_neighbors_ext (forest, itree, element, &neighbors, iface, &dual_faces, &num_neighbors,
                                           &element_indices, &dummy_eclass, &gneigh_tree, nullptr);

        if (num_neighbors == 0) {
          if (!is_ghost) {
            // With a ghost layer, every local leaf finds its neighbors if they exist.
            // No neighbors thus means that the face is on the domain boundary.
            EXPECT_EQ (neigh_class, T8_ECLASS_INVALID)
              << "tree " << itree << " leaf " << ileaf << " face " << iface << ": no neighbor leaves found";
          }
          // For ghost leaves the neighbor may simply not be present on this process.
        }
        else  {
          // The neighbor class must be the class of the forest tree containing the neighbor leaves
          const t8_locidx_t neigh_ltree = t8_forest_get_local_or_ghost_id (forest, gneigh_tree);
          ASSERT_GE (neigh_ltree, 0) << "Neighbor tree " << gneigh_tree << " is neither local nor ghost.";
          EXPECT_EQ (neigh_class, t8_forest_get_tree_class (forest, neigh_ltree))
            << "tree " << itree << " leaf " << ileaf << " face " << iface;

          // Looking back from each neighbor leaf must yield our own eclass
          for (int ineigh = 0; ineigh < num_neighbors; ++ineigh) {
            EXPECT_EQ (t8_forest_element_neighbor_eclass (forest, neigh_ltree, neighbors[ineigh], dual_faces[ineigh]),
                      tree_eclass)
              << "tree " << itree << " leaf " << ileaf << " face " << iface << " neighbor " << ineigh;
          }
        }
        T8_FREE (neighbors);
        T8_FREE (dual_faces);
        T8_FREE (element_indices);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_neighbor_eclass, element_neighbor_eclass,
                          testing::Combine (AllSchemeCollections, testing::Range (0, 4)));
