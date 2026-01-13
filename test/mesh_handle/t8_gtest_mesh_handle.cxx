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
 * Tests if \ref t8_mesh_handle::mesh works as intended for different types of predefined template parameters. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competences.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** Parametrized test fixture for the mesh handle tests. */
struct t8_mesh_handle_test: public testing::TestWithParam<std::tuple<t8_eclass_t, int>>
{
 protected:
  void
  SetUp () override
  {
    const t8_scheme* scheme = t8_scheme_new_default ();
    t8_eclass_t eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 1, 0);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }

  t8_forest_t forest;
  int level;
};

/** Test some default functionality and the iterator of \ref t8_mesh_handle::mesh. */
TEST_P (t8_mesh_handle_test, test_iterator)
{
  // --- Check default functionality. ---
  using mesh_type = t8_mesh_handle::mesh<>;
  using element_type = mesh_type::element_type;
  mesh_type mesh = mesh_type (forest);
  EXPECT_FALSE (element_type::has_vertex_cache ());
  EXPECT_FALSE (element_type::has_centroid_cache ());

  // Iterate with the iterator over all mesh elements and check some functionality.
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_FALSE (it->has_vertex_cache ());
    EXPECT_FALSE (it->has_centroid_cache ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    // Test dereference operator.
    auto vertex_coordinates = (*it).get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }

  // Check loop with indices.
  for (int ielement = 0; ielement < mesh.get_num_local_elements (); ielement++) {
    EXPECT_EQ (level, mesh[ielement].get_level ());
  }
}

/** Test competences. */
TEST_P (t8_mesh_handle_test, test_competences)
{
  // --- Version with cached vertex coordinates. ---
  using mesh_type_vertex = t8_mesh_handle::mesh<t8_mesh_handle::cache_vertex_coordinates>;
  using element_type_vertex = mesh_type_vertex::element_type;
  mesh_type_vertex mesh_vertex = mesh_type_vertex (forest);
  EXPECT_TRUE (element_type_vertex::has_vertex_cache ());
  EXPECT_FALSE (element_type_vertex::has_centroid_cache ());

  // Iterate with the iterator over all mesh elements and check functionality.
  for (auto it = mesh_vertex.begin (); it != mesh_vertex.end (); ++it) {
    EXPECT_FALSE (it->vertex_cache_filled ());
    EXPECT_EQ (level, it->get_level ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }
  // Check cached value.
  for (auto it = mesh_vertex.begin (); it != mesh_vertex.end (); ++it) {
    EXPECT_TRUE (it->vertex_cache_filled ());
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }

  // --- Version with cached centroid variable. ---
  using mesh_type_centroid = t8_mesh_handle::mesh<t8_mesh_handle::cache_centroid>;
  using element_type_centroid = mesh_type_centroid::element_type;
  mesh_type_centroid mesh_centroid = mesh_type_centroid (forest);
  EXPECT_FALSE (element_type_centroid::has_vertex_cache ());
  EXPECT_TRUE (element_type_centroid::has_centroid_cache ());

  // Iterate with the iterator over all mesh elements.
  for (auto it = mesh_centroid.begin (); it != mesh_centroid.end (); ++it) {
    EXPECT_FALSE (it->centroid_cache_filled ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    EXPECT_TRUE (it->centroid_cache_filled ());
    auto centroid_cached = it->get_centroid ();
    for (const auto& coordinate : centroid_cached) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
  }
}

/** Test \ref t8_mesh_handle::mesh with more than one competence. */
TEST_P (t8_mesh_handle_test, test_2_competences)
{
  // --- Use competences to cache level and centroid. ---
  using mesh_type = t8_mesh_handle::mesh<t8_mesh_handle::cache_vertex_coordinates, t8_mesh_handle::cache_centroid>;
  using element_type = mesh_type::element_type;
  mesh_type mesh = mesh_type (forest);
  EXPECT_TRUE (element_type::has_vertex_cache ());
  EXPECT_TRUE (element_type::has_centroid_cache ());

  // Iterate with the iterator over all mesh elements.
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_FALSE (it->centroid_cache_filled ());
    EXPECT_FALSE (it->vertex_cache_filled ());
    EXPECT_EQ (level, it->get_level ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }
  // Test dereference operator. (Here the cached values should be used.)
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_TRUE (it->centroid_cache_filled ());
    EXPECT_TRUE (it->vertex_cache_filled ());
    EXPECT_EQ (level, (*it).get_level ());
    auto centroid = (*it).get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_mesh, t8_mesh_handle_test, testing::Combine (AllEclasses, testing::Range (2, 3)));
