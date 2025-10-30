/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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

/**
 * \file t8_gtest_ghost.cxx
 * Tests if the ghost elements and the get_face_neighbor implementation work as intended.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include <t8.h>

#include <t8_mesh_handle/mesh.hxx>
#include <t8_mesh_handle/competences.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <vector>

/** Parametrized test fixture for the mesh handle tests. */
class t8_mesh_ghost_test: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
 protected:
  void
  SetUp () override
  {
    scheme = t8_scheme_new_default ();
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 1, 0);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 1, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_forest_t forest;
  const t8_scheme* scheme;
  t8_eclass_t eclass;
  int level;
};

/** Check the implementation of ghosts and all functions accessible by ghosts. */
TEST_P (t8_mesh_ghost_test, check_get_face_neighbors)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));
  t8_forest_ghost_print (forest);

  t8_mesh_handle::mesh<> mesh = t8_mesh_handle::mesh<> (forest);
  EXPECT_EQ (mesh.get_local_num_ghosts (), t8_forest_get_num_ghosts (forest));

  // Check that functionality of mesh elements is still valid.
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_FALSE (it->is_ghost_element ());
    EXPECT_EQ (level, it->get_level ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    auto vertex_coordinates = (*it).get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }

  // Check functions for ghost elements.
  const t8_locidx_t num_local_elements = mesh.get_local_num_elements ();
  const t8_locidx_t num_ghost_elements = mesh.get_local_num_ghosts ();
  for (t8_locidx_t ighost = num_local_elements; ighost < num_local_elements + num_ghost_elements; ++ighost) {
    EXPECT_TRUE (mesh[ighost].is_ghost_element ());
    EXPECT_EQ (level, mesh[ighost].get_level ());
    auto centroid = mesh[ighost].get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    auto vertex_coordinates = mesh[ighost].get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }
}

/** Check that the function \ref get_face_neighbors of the handle works as intended (equal results to forest).*/
TEST_P (t8_mesh_ghost_test, compare_neighbors_to_forest)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));

  t8_mesh_handle::mesh<> mesh = t8_mesh_handle::mesh<> (forest);
  EXPECT_EQ (mesh.get_local_num_ghosts (), t8_forest_get_num_ghosts (forest));

  const t8_scheme* scheme = t8_forest_get_scheme (forest);
  auto mesh_iterator = mesh.cbegin ();
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielem) {
      // --- Compare elements. ---
      EXPECT_EQ (mesh_iterator->get_tree_id (), itree);
      EXPECT_EQ (mesh_iterator->get_element_id (), ielem);
      // --- Compare neighbors. ---
      const t8_element_t* elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      const int num_faces = scheme->element_get_num_faces (tree_class, elem);
      EXPECT_EQ (mesh_iterator->get_num_faces (), num_faces);
      for (int iface = 0; iface < num_faces; iface++) {
        // --- Get neighbors from forest. ---
        t8_element_t** neighbors;
        int num_neighbors;
        const int forest_is_balanced = t8_forest_is_balanced (forest);
        t8_eclass_t neigh_eclass;
        int* dual_faces;
        t8_locidx_t* neigh_ids;
        t8_forest_leaf_face_neighbors (forest, itree, elem, &neighbors, iface, &dual_faces, &num_neighbors, &neigh_ids,
                                       &neigh_eclass, forest_is_balanced);
        // --- Get neighbors from mesh element. ---
        int num_neighbors_handle;
        int* dual_faces_handle;
        auto neighbor_ids_handle = mesh_iterator->get_face_neighbors (iface, &num_neighbors_handle, &dual_faces_handle);
        // --- Compare results. ---
        EXPECT_EQ ((int) neighbor_ids_handle.size (), num_neighbors_handle);
        EXPECT_EQ (num_neighbors, num_neighbors_handle);
        EXPECT_TRUE (std::equal (dual_faces, dual_faces + num_neighbors, dual_faces_handle));
        EXPECT_EQ (neighbor_ids_handle, std::vector<t8_locidx_t> (neigh_ids, neigh_ids + num_neighbors));
        // Free memory.
        if (num_neighbors > 0) {
          scheme->element_destroy (neigh_eclass, num_neighbors, neighbors);
          T8_FREE (neigh_ids);
          T8_FREE (neighbors);
          T8_FREE (dual_faces);
          T8_FREE (dual_faces_handle);
        }
      }
      // Evolve mesh iterator.
      mesh_iterator++;
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost, t8_mesh_ghost_test, testing::Combine (AllEclasses, testing::Range (1, 2)));
