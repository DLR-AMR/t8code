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

#include "t8_unstructured_mesh/t8_unstructured_mesh.hxx"
#include "t8_unstructured_mesh/t8_element_competences.hxx"
#include <gtest/gtest.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

TEST (t8_unstructured_mesh, test_iterator)
{
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();

  // Start with a uniform forest.
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_EQ (true, t8_forest_is_committed (forest));

  // Define an unstructured mesh for the forest.
  t8_unstructured_mesh<t8_unstructured_mesh_element<CacheLevel>> unstructured_mesh
    = t8_unstructured_mesh<t8_unstructured_mesh_element<CacheLevel>> (forest);

  // Version without cache.
  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    ASSERT_EQ (level, it->get_level ());
  }

  // Version with cached level variable.
  //unstructured_mesh.cache_level ();
  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh.begin (); it != unstructured_mesh.end (); ++it) {
    ASSERT_EQ (level, (*it).get_level ());
  }

  // Define an unstructured mesh for the forest.
  t8_unstructured_mesh<t8_unstructured_mesh_element<ComputeLevel>> unstructured_mesh_calculate
    = t8_unstructured_mesh<t8_unstructured_mesh_element<ComputeLevel>> (forest);

  // Version without cache.
  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    ASSERT_EQ (level, it->get_level ());
  }

  // Version with cached level variable.
  //unstructured_mesh.cache_level ();
  // Iterate with the iterator over all unstructured mesh elements and check the level.
  for (auto it = unstructured_mesh_calculate.begin (); it != unstructured_mesh_calculate.end (); ++it) {
    ASSERT_EQ (level, (*it).get_level ());
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}
