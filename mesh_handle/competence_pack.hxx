/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file competence_pack.hxx
 * Define classes to pack different competences into one template parameter for the \ref t8_mesh_handle::mesh class.
 * We have one class for element competences and one for mesh competences.
 * There are several predefined competence packs for different use cases, e.g. for all caching competences or
 * all data related competences.
 * Use the union operator to combine different competence packs.
 */

#pragma once

#include "competences.hxx"
#include "data_handler.hxx"
#include "internal/competence_pack_union.hxx"
namespace t8_mesh_handle
{
// --- Element competence pack. ---
/** Class to pack different element competences into one template parameter for the \ref mesh class.
 * \tparam TElementCompetence The competences to be packed.
 */
template <template <typename> class... TElementCompetence>
struct element_competence_pack
{
  /** Apply the competence pack to a template class, e.g. the \ref element class.
   * \tparam TMeshClass The mesh class given to the element class.
   * \tparam Target The target template class to apply the \a TElementCompetence pack to.
   */
  template <typename TMeshClass, template <typename, template <typename> class...> class Target>
  using apply = Target<TMeshClass, TElementCompetence...>;

  using is_element_competence_pack = void; /**< Tag to identify this class. */
};

/** Empty competence pack. */
using empty_element_competences = element_competence_pack<>;

/** Predefined element competence pack combining all caching competences. */
using all_cache_element_competences
  = element_competence_pack<cache_volume, cache_diameter, cache_vertex_coordinates, cache_centroid, cache_face_areas,
                            cache_face_centroids, cache_face_normals, cache_neighbors>;

/** Predefined element competence pack combining all competences related to faces. */
using cache_face_element_competences
  = element_competence_pack<cache_face_areas, cache_face_centroids, cache_face_normals, cache_neighbors>;

/** Predefined element data competence pack. 
 *  Please note that you must combine this with \ref t8_mesh_handle::data_mesh_competences_basic. */
using data_element_competences_basic = element_competence_pack<element_data_element_competence>;

/** Predefined element competence pack combining the element data competence and competence to set new element data. 
 *  Please note that you must combine this with \ref t8_mesh_handle::new_data_mesh_competences. */
using new_data_element_competences
  = element_competence_pack<element_data_element_competence, new_element_data_element_competence>;

// --- Mesh competence pack. ---
/** Class to pack different mesh competences into one template parameter for the \ref mesh class.
 * \tparam TMeshCompetence The mesh competences to be packed.
 */
template <template <typename> class... TMeshCompetence>
struct mesh_competence_pack
{
  /** Apply the mesh competence pack to a mesh type.
   *  By inheriting from all mesh competences, the functionality of the competences gets added to the mesh type.
   *  \tparam TMesh The mesh type to which the competences are applied.
   */
  template <typename TMesh>
  struct apply: public TMeshCompetence<TMesh>...
  {
  };

  using is_mesh_competence_pack = void; /**< Tag to identify this class. */
};

/** Empty competence pack. */
using empty_mesh_competences = mesh_competence_pack<>;

/** Predefined mesh competence pack to handle element data. 
 * If you want to access the data also via the elements, combine this with \ref t8_mesh_handle::data_element_competences_basic.
 */
template <T8MPISafeType TElementDataType>
using data_mesh_competences_basic = mesh_competence_pack<element_data_mesh_competence<TElementDataType>::template type>;

/** Predefined mesh competence pack combining all competences related to data and competence to set new element data. 
 * If you want to access the data also via the elements, combine this with \ref t8_mesh_handle::new_data_element_competences.
 */
template <T8MPISafeType TElementDataType>
using new_data_mesh_competences
  = mesh_competence_pack<element_data_mesh_competence<TElementDataType>::template type,
                         new_element_data_mesh_competence<TElementDataType>::template type>;

// --- Compute union of competence packs. ---
/** Compute the unique union of the competences of several competence_pack. This could be
 * \ref t8_mesh_handle::element_competence_pack or \ref t8_mesh_handle::mesh_competence_pack.
 *  A new competence pack is produced containing all competences of the given competence packs with duplicates removed.
 * \tparam TPacks The competence packs for which we should compute the unique union of the competences.
 *         Each competence pack is expected to be of the same type and of type 
 *         \ref t8_mesh_handle::element_competence_pack or \ref t8_mesh_handle::mesh_competence_pack.
 */
template <typename... TPacks>
using union_competence_packs_type = typename detail::union_competence_packs<TPacks...>::type;

}  // namespace t8_mesh_handle
