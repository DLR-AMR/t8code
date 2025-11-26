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
 * Checks that the cache is actually used if the element gets a cache competence as template parameter. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competences.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

/** Child class of \ref t8_mesh_handle::cache_vertex_coordinates that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_vertex_coordinates_overwrite: public t8_mesh_handle::cache_vertex_coordinates<TUnderlying>
{
 public:
  /** Overwrites the cache variable for the vertex coordinates.
   * \param [in] new_vec New cache vector. 
   */
  void
  overwrite_cache (std::vector<t8_3D_point> new_vec) const
  {
    this->m_vertex_coordinates = new_vec;
  }
};

/** Child class of \ref t8_mesh_handle::cache_centroid that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_centroid_overwrite: public t8_mesh_handle::cache_centroid<TUnderlying>
{
 public:
  /** Overwrites the cache variable for the centroid.
   * \param [in] new_vec New point for the cache. 
   */
  void
  overwrite_cache (t8_3D_point new_vec) const
  {
    this->m_centroid = new_vec;
  }
};

/** Test fixture for cache competence tests. */
class t8_gtest_cache_competence: public testing::Test {
 protected:
  void
  SetUp () override
  {
    level = 1;
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
    const t8_scheme *scheme = t8_scheme_new_default ();
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

/** Use child class of \ref t8_mesh_handle::cache_vertex_coordinates class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST_F (t8_gtest_cache_competence, cache_vertex_coordinates)
{
  using mesh_class = t8_mesh_handle::mesh<cache_vertex_coordinates_overwrite>;
  using element_class = mesh_class::abstract_element_class;
  mesh_class mesh = mesh_class (forest);
  EXPECT_TRUE (element_class::has_vertex_cache ());
  EXPECT_FALSE (element_class::has_centroid_cache ());

  std::vector<t8_3D_point> unrealistic_vertex = { t8_3D_point ({ 41, 42, 43 }), t8_3D_point ({ 99, 100, 101 }) };
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    // Check that cache is empty at the beginning.
    EXPECT_FALSE (it->vertex_cache_filled ());
    // Check that values are valid.
    auto vertex_coordinates = it->get_vertex_coordinates ();
    for (int ivertex = 0; ivertex < (int) vertex_coordinates.size (); ++ivertex) {
      for (const auto &coordinate : vertex_coordinates[ivertex]) {
        EXPECT_GE (1, coordinate);
        EXPECT_LE (0, coordinate);
      }
    }
    // Check that cache is set.
    EXPECT_TRUE (it->vertex_cache_filled ());
    // Overwrite the cache with unrealistic values.
    it->overwrite_cache (unrealistic_vertex);
    // Check that the cache is actually used.
    EXPECT_EQ (it->get_vertex_coordinates (), unrealistic_vertex);
  }
}

/** Use child class of \ref t8_mesh_handle::cache_centroid class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST_F (t8_gtest_cache_competence, cache_centroid)
{
  using mesh_class = t8_mesh_handle::mesh<cache_centroid_overwrite>;
  using element_class = mesh_class::abstract_element_class;
  mesh_class mesh = mesh_class (forest);
  EXPECT_FALSE (element_class::has_vertex_cache ());
  EXPECT_TRUE (element_class::has_centroid_cache ());

  t8_3D_point unrealistic_centroid ({ 999, 1000, 998 });
  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    // Check that cache is empty at the beginning.
    EXPECT_FALSE (it->centroid_cache_filled ());
    // Check that values are valid.
    auto centroid = it->get_centroid ();
    for (const auto &coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    // Check that cache is set.
    EXPECT_TRUE (it->centroid_cache_filled ());
    // Overwrite the cache with unrealistic values.
    it->overwrite_cache (unrealistic_centroid);
    // Check that the cache is actually used.
    EXPECT_EQ (it->get_centroid (), unrealistic_centroid);
  }
}
