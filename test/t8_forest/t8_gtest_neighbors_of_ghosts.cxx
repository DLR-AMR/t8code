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

/** \file t8_gtest_neighbors_of_ghosts.cxx
 * Check that neighbors of ghosts work as expected. This test is very simple but caught a bug.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_adapt_callbacks.hxx>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_schemes/t8_default/t8_default.hxx>

TEST (t8_gtest_neighbors_of_ghosts, check_neighbors)
{
  // Test does not make sense without ghosts.
  int mpisize;
  int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  if (!(mpisize > 1)) {
    GTEST_SKIP () << "Skipping test as no ghost elements can be created.";
  }

  // Construct adapted and partitioned forest with ghost layer.
  const int level = 2;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, 0);
  t8_forest_t forest_uniform = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, true, sc_MPI_COMM_WORLD);
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_adapt (forest, forest_uniform, refine_every_nth_element_callback<2, 0>, 0);
  t8_forest_set_partition (forest, NULL, 0);
  t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  t8_forest_commit (forest);

  // Check neighbors of ghost.
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_locidx_t num_ghost_trees = t8_forest_ghost_num_trees (forest);

  for (t8_locidx_t ighost_tree = 0; ighost_tree < num_ghost_trees; ++ighost_tree) {
    const t8_eclass_t eclass = t8_forest_ghost_get_tree_class (forest, ighost_tree);
    const t8_locidx_t num_elems = t8_forest_ghost_tree_num_leaf_elements (forest, ighost_tree);

    const t8_locidx_t ltreeid = num_local_trees + ighost_tree;

    for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
      const t8_element_t *ghost = t8_forest_ghost_get_leaf_element (forest, ighost_tree, ielem);
      const int num_faces = scheme->element_get_num_faces (eclass, ghost);

      for (int iface = 0; iface < num_faces; ++iface) {
        const t8_element_t **neighbors = nullptr;
        int *dual_faces = nullptr;
        int num_neighbors = 0;
        t8_locidx_t *element_indices = nullptr;
        t8_eclass_t neigh_eclass;

        t8_forest_leaf_face_neighbors (forest, ltreeid, ghost, &neighbors, iface, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_eclass);

        // A face can have at most as many neighbors as it has face children.
        const int max_neighbors = scheme->element_get_num_face_children (eclass, ghost, iface);
        EXPECT_LE (num_neighbors, max_neighbors)
          << "ghost tree " << ighost_tree << " element " << ielem << " face " << iface << ": " << num_neighbors
          << " neighbors, but the face has only " << max_neighbors << " face children.";

        // No neighbour may be reported twice.
        const std::set<t8_locidx_t> unique_indices (element_indices, element_indices + num_neighbors);
        EXPECT_EQ (static_cast<int> (unique_indices.size ()), num_neighbors) << "Duplicated neighbors found.";

        if (num_neighbors > 0) {
          T8_FREE (neighbors);
          T8_FREE (dual_faces);
          T8_FREE (element_indices);
        }
      }
    }
  }

  t8_forest_unref (&forest);
}
