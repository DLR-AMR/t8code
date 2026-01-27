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
 * Tests if the \ref t8_mesh_handle::mesh class of the handle works as intended for different types of template parameter classes. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competences.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>

/** Parametrized test fixture for the mesh handle tests. */
struct t8_mesh_handle_test: public testing::TestWithParam<std::tuple<t8_eclass_t, int>>
{
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());
  }

  t8_eclass_t eclass; /**< Element class used for testing.*/
  int level;          /**< Refinement level used for testing.*/
};

/** Test default \ref t8_mesh_handle::mesh handle, the iterator and some exemplary functionality. */
TEST_P (t8_mesh_handle_test, test_default_mesh_handle)
{
  using mesh_class = t8_mesh_handle::mesh<>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (eclass, level, sc_MPI_COMM_WORLD, true,
                                                                            true, false);
  EXPECT_FALSE (element_class::has_vertex_cache ());
  EXPECT_FALSE (element_class::has_centroid_cache ());

  // Iterate with the constant iterator over all mesh elements and check some functionality.
  for (auto it = mesh->cbegin (); it != mesh->cend (); ++it) {
    EXPECT_FALSE (it->has_vertex_cache ());
    EXPECT_FALSE (it->has_centroid_cache ());
    auto centroid = it->get_centroid ();
    for (const auto& coordinate : centroid) {
      EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
    }
    auto vertex_coordinates = (*it).get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto& coordinate : vertex_coordinates[ivertex]) {
        EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
      }
    }
  }

  // Check loop with indices.
  for (int ielement = 0; ielement < mesh->get_num_local_elements (); ielement++) {
    EXPECT_EQ (level, (*mesh)[ielement].get_level ());
    EXPECT_EQ (ielement, (*mesh)[ielement].get_element_handle_id ());
  }
}

/** Test mesh class with all competence available using some exemplary functionality. */
TEST_P (t8_mesh_handle_test, test_all_cache_competence)
{
  // --- Use predefined competences to use all available caching competences. ---
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::all_cache_competences>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (eclass, level, sc_MPI_COMM_WORLD, true,
                                                                            true, false);
  EXPECT_TRUE (element_class::has_volume_cache ());
  EXPECT_TRUE (element_class::has_diameter_cache ());
  EXPECT_TRUE (element_class::has_vertex_cache ());
  EXPECT_TRUE (element_class::has_centroid_cache ());
  EXPECT_TRUE (element_class::has_face_area_cache ());
  EXPECT_TRUE (element_class::has_face_centroid_cache ());
  EXPECT_TRUE (element_class::has_face_normal_cache ());
  EXPECT_TRUE (element_class::has_face_neighbor_cache ());

  // Iterate over all mesh elements and access some exemplary functionality which sets the caches.
  for (auto it = mesh->cbegin (); it != mesh->cend (); ++it) {
    EXPECT_FALSE (it->volume_cache_filled ());
    EXPECT_FALSE (it->centroid_cache_filled ());
    EXPECT_FALSE (it->vertex_cache_filled ());
    EXPECT_EQ (level, it->get_level ());
    EXPECT_GE (it->get_volume (), 0);
    for (const auto& coordinate : it->get_centroid ()) {
      EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
    }
    for (int ivertex = 0; ivertex < it->get_num_vertices (); ++ivertex) {
      for (const auto& coordinate : it->get_vertex_coordinates (ivertex)) {
        EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
      }
    }
  }
  // Check if caches are set. If caches are accessed correctly is checked in another test.
  for (auto it = mesh->cbegin (); it != mesh->cend (); ++it) {
    EXPECT_TRUE (it->volume_cache_filled ());
    EXPECT_TRUE (it->centroid_cache_filled ());
    EXPECT_TRUE (it->vertex_cache_filled ());
    EXPECT_FALSE (it->diameter_cache_filled ());
    if (it->get_num_faces () > 0) {
      EXPECT_FALSE (it->face_area_cache_filled (0));
      EXPECT_FALSE (it->face_centroid_cache_filled (0));
      EXPECT_FALSE (it->face_normal_cache_filled (0));
    }
    EXPECT_FALSE (it->neighbor_cache_filled_any ());
  }
}

/** Test mesh class with all predefined face competences using some exemplary functionality. */
TEST_P (t8_mesh_handle_test, test_cache_face_competences)
{
  // --- Use all predefined competences. ---
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::cache_face_competences>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (eclass, level, sc_MPI_COMM_WORLD, true,
                                                                            true, false);
  EXPECT_FALSE (element_class::has_volume_cache ());
  EXPECT_FALSE (element_class::has_diameter_cache ());
  EXPECT_FALSE (element_class::has_vertex_cache ());
  EXPECT_FALSE (element_class::has_centroid_cache ());
  EXPECT_TRUE (element_class::has_face_area_cache ());
  EXPECT_TRUE (element_class::has_face_centroid_cache ());
  EXPECT_TRUE (element_class::has_face_normal_cache ());
  EXPECT_TRUE (element_class::has_face_neighbor_cache ());

  // Iterate over all mesh elements and access some exemplary functionality which sets the caches.
  for (auto it = mesh->cbegin (); it != mesh->cend (); ++it) {
    if (it->get_num_faces () == 0) {
      GTEST_SKIP () << "Skipping test as there are no faces to be tested.";
    }
    for (int iface = 0; iface < it->get_num_faces (); ++iface) {
      EXPECT_FALSE (it->face_area_cache_filled (iface));
      EXPECT_FALSE (it->face_normal_cache_filled (iface));
      EXPECT_GE (it->get_face_area (iface), 0);
      for (const auto& coordinate : it->get_face_normal (iface)) {
        EXPECT_TRUE (coordinate >= -1 && coordinate <= 1);
      }
    }
  }
  // Check if caches are set. If caches are accessed correctly is checked in another test.
  for (auto it = mesh->cbegin (); it != mesh->cend (); ++it) {
    for (int iface = 0; iface < it->get_num_faces (); ++iface) {
      EXPECT_TRUE (it->face_area_cache_filled (iface));
      EXPECT_FALSE (it->face_centroid_cache_filled (iface));
      EXPECT_TRUE (it->face_normal_cache_filled (iface));
    }
    EXPECT_FALSE (it->neighbor_cache_filled_any ());
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_mesh, t8_mesh_handle_test, testing::Combine (AllEclasses, testing::Range (2, 3)));
