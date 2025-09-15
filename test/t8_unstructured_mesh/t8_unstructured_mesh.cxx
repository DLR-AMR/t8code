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
#include <t8.h>

#include <t8_unstructured_mesh/t8_unstructured_mesh.hxx>
#include <t8_unstructured_mesh/t8_element_competences.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** Test some default functionality and the iterator class of the unstructured mesh. */
TEST (t8_unstructured_mesh, test_iterator)
{
  // Simple example forest.
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Check default functionality. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<>> unstructured_mesh_calculate
    = t8_unstructured_mesh<t8_unstructured_mesh_element<>> (forest);

  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
  }
  // Test dereference operator.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    EXPECT_EQ (level, (*it).get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, (*it).get_centroid ()[coord]);
      EXPECT_LE (0, (*it).get_centroid ()[coord]);
    }
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}

/** Test competences. */
TEST (t8_unstructured_mesh, test_competences)
{
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Version with cached level variable. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_level>> unstructured_mesh_level
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_level>> (forest);

  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh_level.begin (); it != unstructured_mesh_level.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
  }
  // Test dereference operator. (Here the cached value should be used.)
  for (auto it = unstructured_mesh_level.begin (); it != unstructured_mesh_level.end (); ++it) {
    EXPECT_EQ (level, (*it).get_level ());
  }

  // --- Version with cached centroid variable. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_centroid>> unstructured_mesh_centroid
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_centroid>> (forest);

  // Iterate with the iterator over all unstructured mesh elements.
  for (auto it = unstructured_mesh_centroid.begin (); it != unstructured_mesh_centroid.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      // Second call (here cached value should be used).
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}

/** Test unstructured mesh (element) class with more than one competence. */
TEST (t8_unstructured_mesh, test_2_competences)
{
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_TRUE (t8_forest_is_committed (forest));

  // --- Use competences to cache level and centroid. ---
  t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_level, t8_cache_centroid>> unstructured_mesh
    = t8_unstructured_mesh<t8_unstructured_mesh_element<t8_cache_level, t8_cache_centroid>> (forest);

  // Iterate with the iterator over all unstructured mesh elements.
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    EXPECT_EQ (level, it->get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, it->get_centroid ()[coord]);
      EXPECT_LE (0, it->get_centroid ()[coord]);
    }
  }
  // Test dereference operator. (Here the cached values should be used.)
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    EXPECT_EQ (level, (*it).get_level ());
    for (int coord = 0; coord < T8_ECLASS_MAX_DIM; ++coord) {
      EXPECT_GE (1, (*it).get_centroid ()[coord]);
      EXPECT_LE (0, (*it).get_centroid ()[coord]);
    }
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}
