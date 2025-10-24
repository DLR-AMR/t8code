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
 * \file t8_gtest_cache_competence.cxx
 * Checks that the competences to cache functionality actually use the cache. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8.h>

#include <t8_mesh_interface/t8_interface_mesh.hxx>
#include <t8_mesh_interface/t8_interface_element.hxx>
#include <t8_mesh_interface/t8_interface_competences.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

/** TODO. */
template <typename TUnderlying>
struct t8_cache_vertex_coordinates_overwrite: public t8_cache_vertex_coordinates<TUnderlying>
{
 public:
  void
  overwrite_cache (std::vector<t8_3D_vec> new_vec)
  {
    this->m_vertex_coordinates = new_vec;
  }
};

/** TODO */
TEST (t8_gtest_cache_competence, cache_competence)
{
  // Define forest to construct mesh.
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_EQ (true, t8_forest_is_committed (forest));

  // Check mesh with custom defined competence.
  using mesh_element = t8_interface_element<t8_cache_vertex_coordinates_overwrite>;
  t8_interface_mesh<mesh_element> mesh = t8_interface_mesh<mesh_element> (forest);
  EXPECT_TRUE (mesh_element::has_vertex_cache ());
  EXPECT_FALSE (mesh_element::has_centroid_cache ());

  std::vector<t8_3D_vec> unrealistic_vertex = { t8_3D_vec ({ 41, 42, 43 }), t8_3D_vec ({ 99, 100, 101 }) };
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_FALSE (it->vertex_cache_filled ());
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto &coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
    EXPECT_TRUE (it->vertex_cache_filled ());
    // Overwrite the cache with unrealistic values.
    it->overwrite_cache (unrealistic_vertex);
    EXPECT_TRUE (it->vertex_cache_filled ());
    EXPECT_EQ (it->get_vertex_coordinates (), unrealistic_vertex);
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}
