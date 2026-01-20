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
#include <mesh_handle/constructor_wrappers.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

/** Child of \ref t8_mesh_handle::cache_volume that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_volume_overwrite: public t8_mesh_handle::cache_volume<TUnderlying>
{
  /** Overwrites the cache variable for the volume.
   * \param [in] new_volume New volume. 
   */
  void
  overwrite_cache (double new_volume) const
  {
    this->m_volume = new_volume;
  }
};

/** Child of \ref t8_mesh_handle::cache_diameter that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_diameter_overwrite: public t8_mesh_handle::cache_diameter<TUnderlying>
{
  /** Overwrites the cache variable for the diameter.
   * \param [in] new_diameter New diameter. 
   */
  void
  overwrite_cache (double new_diameter) const
  {
    this->m_diameter = new_diameter;
  }
};

/** Child of \ref t8_mesh_handle::cache_vertex_coordinates that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_vertex_coordinates_overwrite: public t8_mesh_handle::cache_vertex_coordinates<TUnderlying>
{
  /** Overwrites the cache variable for the vertex coordinates.
   * \param [in] new_vec New cache vector. 
   */
  void
  overwrite_cache (std::vector<t8_3D_point> new_vec) const
  {
    this->m_vertex_coordinates = new_vec;
  }
};

/** Child of \ref t8_mesh_handle::cache_centroid that allows to modify the cache variable for test purposes. */
template <typename TUnderlying>
struct cache_centroid_overwrite: public t8_mesh_handle::cache_centroid<TUnderlying>
{
  /** Overwrites the cache variable for the centroid.
   * \param [in] new_vec New point for the cache. 
   */
  void
  overwrite_cache (t8_3D_point new_vec) const
  {
    this->m_centroid = new_vec;
  }
};

// --- Face related caches. ---
/** Child of \ref t8_mesh_handle::cache_face_area that allows to modify the cache variables for test purposes. */
template <typename TUnderlying>
struct cache_face_area_overwrite: public t8_mesh_handle::cache_face_area<TUnderlying>
{
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

/** Child of \ref t8_mesh_handle::cache_face_centroid that allows to modify the cache variables for test purposes. */
template <typename TUnderlying>
struct cache_face_centroid_overwrite: public t8_mesh_handle::cache_face_centroid<TUnderlying>
{
  /** Overwrites the cache variable for the face centroid.
   * \param [in] face The face for which the cache should be set.
   * \param [in] new_face_centroid New face centroid. 
   */
  void
  overwrite_cache (int face, t8_3D_point new_face_centroid) const
  {
    this->m_face_centroid[face] = new_face_centroid;
  }
};

/** Child of \ref t8_mesh_handle::cache_face_normal that allows to modify the cache variables for test purposes. */
template <typename TUnderlying>
struct cache_face_normal_overwrite: public t8_mesh_handle::cache_face_normal<TUnderlying>
{
  /** Overwrites the cache variable for the face normal.
   * \param [in] face The face for which the cache should be set.
   * \param [in] new_face_normal New face normal. 
   */
  void
  overwrite_cache (int face, t8_3D_vec new_face_normal) const
  {
    this->m_face_normal[face] = new_face_normal;
  }
};

// --- Element related cache tests. ---
/** Use child of \ref t8_mesh_handle::cache_volume to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_volume)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_volume_overwrite>;
  using element_class = typename mesh_class::element_class;
  const auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
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

/** Use child of \ref t8_mesh_handle::cache_diameter to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_diameter)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_diameter_overwrite>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
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

/** Use child of \ref t8_mesh_handle::cache_vertex_coordinates to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_vertex_coordinates)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_vertex_coordinates_overwrite>;
  using element_class = typename mesh_class::element_class;
  const auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_vertex_cache ());

  std::vector<t8_3D_point> unrealistic_vertex = { t8_3D_point ({ 41, 42, 43 }), t8_3D_point ({ 99, 100, 101 }) };
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    EXPECT_FALSE (it->vertex_cache_filled ());
    for (int ivertex = 0; ivertex < it->get_num_vertices (); ++ivertex) {
      for (const auto &coordinate : it->get_vertex_coordinates (ivertex)) {
        EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
      }
    }
    EXPECT_TRUE (it->vertex_cache_filled ());
    it->overwrite_cache (unrealistic_vertex);
    EXPECT_EQ (it->get_vertex_coordinates (), unrealistic_vertex);
  }
}

/** Use child of \ref t8_mesh_handle::cache_centroid to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_centroid)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_centroid_overwrite>;
  using element_class = mesh_class::element_class;
  const auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_centroid_cache ());

  t8_3D_point unrealistic_centroid ({ 999, 1000, 998 });
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    EXPECT_FALSE (it->centroid_cache_filled ());
    for (const auto &coordinate : it->get_centroid ()) {
      EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
    }
    EXPECT_TRUE (it->centroid_cache_filled ());
    it->overwrite_cache (unrealistic_centroid);
    EXPECT_EQ (it->get_centroid (), unrealistic_centroid);
  }
}

// --- Face related cache tests. ---
/** Use child of \ref t8_mesh_handle::cache_face_area to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_face_area)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_face_area_overwrite>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
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
    }
  }
}

/** Use child of \ref t8_mesh_handle::cache_face_centroid to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_face_centroid)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_face_centroid_overwrite>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_face_centroid_cache ());

  t8_3D_point unrealistic_face_centroid ({ 999, 1000, 998 });
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    for (int iface = 0; iface < it->get_num_faces (); ++iface) {
      EXPECT_FALSE (it->face_centroid_cache_filled (iface));
      for (const auto &coordinate : it->get_face_centroid (iface)) {
        EXPECT_TRUE (coordinate >= 0 && coordinate <= 1);
      }
      EXPECT_TRUE (it->face_centroid_cache_filled (iface));
      it->overwrite_cache (iface, unrealistic_face_centroid);
      EXPECT_EQ (it->get_face_centroid (iface), unrealistic_face_centroid);
    }
  }
}

/** Use child of \ref t8_mesh_handle::cache_face_normal to check that the cache is actually set 
 * and accessed correctly. This is done by modifying the cache to an unrealistic value and 
 * checking that the functionality actually outputs this unrealistic value.
 */
TEST (t8_gtest_cache_competence, cache_face_normal)
{
  const int level = 1;
  using mesh_class = t8_mesh_handle::mesh<cache_face_normal_overwrite>;
  using element_class = typename mesh_class::element_class;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);
  EXPECT_TRUE (element_class::has_face_normal_cache ());

  t8_3D_vec unrealistic_face_normal ({ 41, 42, 43 });
  for (auto it = mesh->begin (); it != mesh->end (); ++it) {
    for (int iface = 0; iface < it->get_num_faces (); ++iface) {
      EXPECT_FALSE (it->face_normal_cache_filled (iface));
      for (const auto &coordinate : it->get_face_normal (iface)) {
        EXPECT_TRUE (coordinate >= -1 && coordinate <= 1);
      }
      EXPECT_TRUE (it->face_normal_cache_filled (iface));
      it->overwrite_cache (iface, unrealistic_face_normal);
      EXPECT_EQ (it->get_face_normal (iface), unrealistic_face_normal);
    }
  }
}
