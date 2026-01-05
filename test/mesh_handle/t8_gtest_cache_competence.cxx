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
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competences.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrapper.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

/** Child class of \ref t8_mesh_handle::cache_volume that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_volume_overwrite: public t8_mesh_handle::cache_volume<TUnderlying>
{
 public:
  /** Overwrites the cache variable for the volume.
   * \param [in] new_volume New volume. 
   */
  void
  overwrite_cache (double new_volume) const
  {
    this->m_volume = new_volume;
  }
};

/** Child class of \ref t8_mesh_handle::cache_diameter that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_diameter_overwrite: public t8_mesh_handle::cache_diameter<TUnderlying>
{
 public:
  /** Overwrites the cache variable for the diameter.
   * \param [in] new_diameter New diameter. 
   */
  void
  overwrite_cache (double new_diameter) const
  {
    this->m_diameter = new_diameter;
  }
};

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

/** Use child class of \ref t8_mesh_handle::cache_volume class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_volume)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<cache_volume_overwrite>>;
  using element_class = typename mesh_class::abstract_element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_volume_cache ());

  double unrealistic_volume = -3000;
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    // Check that cache is empty at the beginning.
    EXPECT_FALSE (it->volume_cache_filled ());
    // Fill cache and check that volume is valid.
    EXPECT_GE (it->get_volume (), 0);
    // Check that cache is set.
    EXPECT_TRUE (it->volume_cache_filled ());
    // Overwrite the cache with unrealistic values.
    it->overwrite_cache (unrealistic_volume);
    // Check that the cache is actually used.
    EXPECT_EQ (it->get_volume (), unrealistic_volume);
  }
}

/** Use child class of \ref t8_mesh_handle::cache_diameter class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_diameter)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<cache_diameter_overwrite>>;
  using element_class = typename mesh_class::abstract_element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_diameter_cache ());

  double unrealistic_diameter = -3000;
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    EXPECT_FALSE (it->diameter_cache_filled ());
    EXPECT_GE (it->get_diameter (), 0);
    EXPECT_TRUE (it->diameter_cache_filled ());
    it->overwrite_cache (unrealistic_diameter);
    EXPECT_EQ (it->get_diameter (), unrealistic_diameter);
  }
}

/** Use child class of \ref t8_mesh_handle::cache_vertex_coordinates class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_vertex_coordinates)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<cache_vertex_coordinates_overwrite>>;
  using element_class = typename mesh_class::abstract_element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_vertex_cache ());

  std::vector<t8_3D_point> unrealistic_vertex = { t8_3D_point ({ 41, 42, 43 }), t8_3D_point ({ 99, 100, 101 }) };
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
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
TEST (t8_gtest_cache_competence, cache_centroid)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<cache_centroid_overwrite>>;
  using element_class = mesh_class::abstract_element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_centroid_cache ());

  t8_3D_point unrealistic_centroid ({ 999, 1000, 998 });
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
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

// --- Face related caches. ---

/** Child class of \ref t8_mesh_handle::cache_face_area that allows to modify the cache variables for test purposes. */
template <typename TUnderlying>
struct cache_face_area_overwrite: public t8_mesh_handle::cache_face_area<TUnderlying>
{
 public:
  /** Overwrites the cache variable for the face area.
   * \param [in] face The face for which the cache should be set.
   * \param [in] new_face_area New face area. 
   */
  void
  overwrite_cache (int face, double new_face_area) const
  {
    this->m_face_area[face] = new_face_area;
  }
};

/** Use child class of \ref t8_mesh_handle::cache_face_area class to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_face_area)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<cache_face_area_overwrite>>;
  using element_class = typename mesh_class::abstract_element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_face_area_cache ());

  double unrealistic_face_area = 41.1;
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    for (int iface = 0; iface < it->get_num_faces (); ++iface) {
      EXPECT_FALSE (it->face_area_cache_filled (iface));
      auto face_area = it->get_face_area (iface);
      EXPECT_LE (0, face_area);
      EXPECT_TRUE (it->face_area_cache_filled (iface));
      it->overwrite_cache (iface, unrealistic_face_area + iface);
      EXPECT_EQ (it->get_face_area (iface), unrealistic_face_area + iface);
      std::cout << "here";
    }
  }
}

//TODO: Add tests for face centroid and face normal caches. Also add these to compare handle to forest.
