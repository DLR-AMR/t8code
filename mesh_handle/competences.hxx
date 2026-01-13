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

/** \file competences.hxx
 * Definition of the additional competences/functionalities that can be used for the mesh class.
 * Especially, competences to cache functionalities of elements instead of calculating them each time a function
 * is called are provided.
 *
 * All competences have the same inheritance pattern: 
 * We use the CRTP pattern as we may need to access members of the derived classes like 
 * \ref t8_mesh_handle::element. 
 * The t8_crtp_operator is used for convenience/clear code (avoid to type a static cast explicitly each time 
 * we need functionality of TUnderlying).
 * Especially for the competences to cache functionality, the access of members is not necessary, 
 * such that it is not obvious why we use the crtp. For competences that extend the functionality of the element, 
 * this is required. 
 * We use it for all competences for consistency and compatibility with the \ref t8_mesh_handle::element class.
 */

#ifndef T8_COMPETENCES_HXX
#define T8_COMPETENCES_HXX

#include <t8.h>
#include <t8_types/t8_operators.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>
#include <optional>

namespace t8_mesh_handle
{

/**
 * Competence to cache the volume of an element at the first function call.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct cache_volume: public t8_crtp_operator<TUnderlying, cache_volume>
{
 public:
  /**
   * Function that checks if the cache for the volume has been filled.
   * \return true if the cache has been filled, false otherwise.
   */
  bool
  volume_cache_filled () const
  {
    return m_volume.has_value ();
  }

 protected:
  mutable std::optional<double>
    m_volume; /**< Cache for the volume. Use optional to allow no value if cache is not filled. */
};

/**
 * Competence to cache the vertex coordinates of an element at the first function call.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct cache_vertex_coordinates: public t8_crtp_operator<TUnderlying, cache_vertex_coordinates>
{
 public:
  /**
   * Function that checks if the cache for the vertex coordinates has been filled.
   * \return true if the cache for the vertex coordinates has been filled, false otherwise.
   */
  bool
  vertex_cache_filled () const
  {
    return !m_vertex_coordinates.empty ();
  }

 protected:
  mutable std::vector<t8_3D_point>
    m_vertex_coordinates; /**< Cache for the vector of vertex coordinate arrays. Empty vector if not filled. */
};

/**
 * Competence to cache the centroid of an element at the first function call.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct cache_centroid: public t8_crtp_operator<TUnderlying, cache_centroid>
{
 public:
  /**
   * Function that checks if the cache for the centroid has been filled.
   * \return true if the cache for the centroid has been filled, false otherwise.
   */
  bool
  centroid_cache_filled () const
  {
    return m_centroid.has_value ();
  }

 protected:
  mutable std::optional<t8_3D_point>
    m_centroid; /**< Cache for the coordinates of the centroid. Use optional to allow no value if cache is not filled. */
};

/**
 * Competence to cache the neighbors of an element at a specific face at the first function call.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct cache_neighbors: t8_crtp_operator<TUnderlying, cache_neighbors>
{
 public:
  /**
   * Function that checks if the neighbor cache for a face has been filled.
   * \param [in] face The face for which the cache should be checked.
   * \return true if the cache has been filled, false otherwise.
   */
  bool
  neighbor_cache_filled (int face) const
  {
    return m_num_neighbors[face].has_value ();
  }

  /**
   * Function that checks if the neighbor cache for any face has been filled.
   * \return true if the cache has been filled, false otherwise.
   */
  bool
  neighbor_cache_filled_any () const
  {
    for (int iface = 0; iface < this->underlying ().get_num_faces (); ++iface) {
      if (neighbor_cache_filled (iface)) {
        return true;
      }
    }
    return false;
  }

 protected:
  mutable std::vector<std::vector<const TUnderlying *>>
    m_neighbors; /**< Neighboring elements at each face. The length of the vectors is stored in \ref m_num_neighbors. */
  mutable std::vector<std::optional<int>>
    m_num_neighbors; /**< Vector with the numbers of neighbor elements at each face. 
                        num_neighbors is stored to indicate that the cache is filled if a face does not have any neighbor. */
  mutable std::vector<std::vector<int>>
    m_dual_faces; /**< Face id's of the neighboring elements' faces for each face. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_COMPETENCES_HXX */
