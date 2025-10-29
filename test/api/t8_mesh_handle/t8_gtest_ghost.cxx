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
 * \file t8_gtest_mesh_handle.cxx
 * Tests if the mesh class of the handle works as intended for different types of predefined template parameter classes. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include <t8.h>

#include <t8_mesh_handle/mesh.hxx>
#include <t8_mesh_handle/element.hxx>
#include <t8_mesh_handle/competences.hxx>
#include <t8_mesh_handle/ghost.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** Parametrized test fixture for the mesh handle tests. */
class t8_mesh_handle_test: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
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

/** Test some default functionality and the iterator of \ref t8_mesh_handle::mesh class. */
TEST_P (t8_mesh_handle_test, test_iterator)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));
  t8_forest_ghost_print (forest);

  // --- Check default functionality. ---
  t8_mesh_handle::mesh<> mesh = t8_mesh_handle::mesh<> (forest);
  EXPECT_EQ (mesh.get_local_num_ghosts (), t8_forest_get_num_ghosts (forest));

  // Iterate with the iterator over all mesh elements and check some functionality.
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_FALSE (it->is_ghost_element ());
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

  for (int ighost = mesh.get_local_num_elements ();
       ighost < mesh.get_local_num_elements () + mesh.get_local_num_ghosts (); ++ighost) {

    if (auto* ghost = static_cast<t8_mesh_handle::ghost_element<>*> (&mesh[ighost])) {
      EXPECT_TRUE (ghost->is_ghost_element ());
    }

    // Idea: virtual class element with children mesh_element and ghost_element.
    // Alternative is to let the user choose the right function ghost or not ghost with [].

    // auto centroid = it->get_centroid ();
    // for (const auto& coordinate : centroid) {
    //   EXPECT_GE (1, coordinate);
    //   EXPECT_LE (0, coordinate);
    // }
    // auto vertex_coordinates = (*it).get_vertex_coordinates ();
    // for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
    //   for (const auto& coordinate : vertex_coordinates[ivertex]) {
    //     EXPECT_GE (1, coordinate);
    //     EXPECT_LE (0, coordinate);
    //   }
    // }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_mesh, t8_mesh_handle_test, testing::Combine (AllEclasses, testing::Range (2, 3)));
