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
 * \file t8_gtest_compare_handle_to_forest.cxx
 * Tests that the functionality of the handle gives the same results as if worked with the forest directly.
 */

#include <gtest/gtest.h>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** Tests that the functionality of the handle gives the same results as if worked with the forest directly.
 * Therefore, we compare the results of the mesh handle with the one accessed with the forest directly.
 */
TEST (t8_gtest_compare_handle_to_forest, compare_handle_to_forest)
{
  // Define forest and mesh handle.
  const int level = 2;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *init_scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, init_scheme, level, 0, sc_MPI_COMM_WORLD);

  t8_mesh_handle::mesh<> mesh = t8_mesh_handle::mesh<> (forest);

  // Iterate over the elements of the forest and of the mesh handle simultaneously and compare results of different functions.
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  // Use iterator for the mesh handle and two nested for-loops for the forest.
  auto mesh_iterator = mesh.begin ();
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielem) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      // --- Compare elements. ---
      EXPECT_EQ (mesh_iterator->get_local_tree_id (), itree);
      EXPECT_EQ (mesh_iterator->get_local_element_id (), ielem);
      // --- Compare basics. ---
      EXPECT_EQ (mesh_iterator->get_level (), scheme->element_get_level (tree_class, elem));
      EXPECT_EQ (mesh_iterator->get_shape (), scheme->element_get_shape (tree_class, elem));
      EXPECT_EQ (mesh_iterator->get_volume (), t8_forest_element_volume (forest, itree, elem));
      // --- Compare centroid. ---
      t8_3D_point centroid;
      t8_forest_element_centroid (forest, itree, elem, centroid.data ());
      EXPECT_EQ (mesh_iterator->get_centroid (), centroid);
      // --- Compare vertex coordinates. ---
      auto vertex_coordinates = mesh_iterator->get_vertex_coordinates ();
      for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
        t8_3D_point vertex_forest;
        t8_forest_element_coordinate (forest, itree, elem, ivertex, vertex_forest.data ());
        EXPECT_EQ (vertex_forest, vertex_coordinates[ivertex]);
      }
      // --- Evolve mesh iterator. ---
      mesh_iterator++;
    }
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}
