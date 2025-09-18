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

/** Tests for the unstructured mesh class. */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8.h>

#include <t8_unstructured_mesh/t8_unstructured_mesh.hxx>
#include <t8_unstructured_mesh/t8_element_competences.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** Construct a forest of a hypercube with volume 1. If the elements are refined uniformly,
 * all elements have volume 1/global_num_elements. 
 */
class t8_unstructured_mesh_test: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_forest_t forest;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
  int level;
};

/** Test some default functionality and the iterator class of the unstructured mesh. */
TEST_P (t8_unstructured_mesh_test, test_iterator)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Check default functionality. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<>> unstructured_mesh_calculate
    = t8_unstructured_mesh<t8_unstructured_mesh_element<>> (forest);

  // Iterate with the iterator over all unstructured mesh elements and check some functionality.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
        EXPECT_GE (1, vertex_coordinates[ivertex][coord]);
        EXPECT_LE (0, vertex_coordinates[ivertex][coord]);
      }
    }
  }
  // Test dereference operator.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    EXPECT_EQ (level, (*it).get_level ());
  }
}

/** Test competences. */
TEST_P (t8_unstructured_mesh_test, test_competences)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Version with cached vertex coordinates. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_vertex_coordinates>> unstructured_mesh_vertex_coordinates
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_vertex_coordinates>> (forest);

  // Iterate with the iterator over all unstructured mesh elements and check functionality.
  for (auto it = unstructured_mesh_vertex_coordinates.begin (); it != unstructured_mesh_vertex_coordinates.end ();
       ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
        EXPECT_GE (1, vertex_coordinates[ivertex][coord]);
        EXPECT_LE (0, vertex_coordinates[ivertex][coord]);
      }
    }
  }
  // Test dereference operator. (Here the cached value should be used.)
  for (auto it = unstructured_mesh_vertex_coordinates.begin (); it != unstructured_mesh_vertex_coordinates.end ();
       ++it) {
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
        EXPECT_GE (1, vertex_coordinates[ivertex][coord]);
        EXPECT_LE (0, vertex_coordinates[ivertex][coord]);
      }
    }
  }

  // --- Version with cached centroid variable. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_centroid>> unstructured_mesh_centroid
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_centroid>> (forest);

  // Iterate with the iterator over all unstructured mesh elements.
  for (auto it = unstructured_mesh_centroid.begin (); it != unstructured_mesh_centroid.end (); ++it) {
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      // Second call (here cached value should be used).
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
  }
}

/** Test unstructured mesh (element) class with more than one competence. */
TEST_P (t8_unstructured_mesh_test, test_2_competences)
{
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Use competences to cache level and centroid. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_vertex_coordinates, t8_cache_centroid>> unstructured_mesh
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_vertex_coordinates, t8_cache_centroid>> (forest);

  // Iterate with the iterator over all unstructured mesh elements.
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
        EXPECT_GE (1, vertex_coordinates[ivertex][coord]);
        EXPECT_LE (0, vertex_coordinates[ivertex][coord]);
      }
    }
  }
  // Test dereference operator. (Here the cached values should be used.)
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    EXPECT_EQ (level, (*it).get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, (*it).get_centroid ()[coord]);
      EXPECT_LE (0, (*it).get_centroid ()[coord]);
    }
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
        EXPECT_GE (1, vertex_coordinates[ivertex][coord]);
        EXPECT_LE (0, vertex_coordinates[ivertex][coord]);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_unstructured_mesh, t8_unstructured_mesh_test,
                          testing::Combine (AllSchemes, testing::Range (0, 4)));
