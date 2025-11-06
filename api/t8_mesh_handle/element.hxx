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

/** \file element.hxx
 * Definition of the elements used in the mesh class.
 */

#ifndef T8_ELEMENT_HXX
#define T8_ELEMENT_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_eclass.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_types/t8_vec.hxx>
#include <array>
#include <vector>

namespace t8_mesh_handle
{
/* Forward declaration of the mesh class of the handle.
 */
template <class TMeshElement>
class mesh;

/** 
 * Class for the elements of the mesh handle. 
 * The element without specified template parameters provides default implementations for basic functionality 
 * as accessing the refinement level or the centroid. With this implementation, the functionality is calculated each time
 * the function is called. 
 * Use the competences defined in \ref competences.hxx as template parameter to cache the functionality instead of 
 * recalculating in every function call.
 * To add functionality to the element, you can also simply write your own competence class and give it as a template parameter.
 * You can access the functions implemented in your competence via the element. 
 *
 * The inheritance pattern is inspired by the \ref T8Type class (which also uses the CRTP).
 * We decided to use this structure 1.) to be able to add new functionality easily and 
 *    2.) for the cached options to keep the number of class member variables of the default to a minimum to save memory.
 * The choice between calculate and cache is a tradeoff between runtime and memory usage. 
 *
 * \tparam TCompetence The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class element: public TCompetence<element<TCompetence...>>... {
  using SelfType = element<TCompetence...>;

 private:
  // --- Variables to check which functionality is defined in TCompetence. ---
  /** Helper function to check if class T implements the function vertex_cache_filled.
   * \tparam T The competence to be checked.
   * \return true if T implements the function, false if not.
   */
  template <template <typename> class T>
  static constexpr bool
  vertex_cache_defined ()
  {
    return requires (T<SelfType>& competence) { competence.vertex_cache_filled (); };
  }
  /* This variable is true if any of the given competences \ref TCompetence implements 
  a function vertex_cache_filled */
  static constexpr bool vertex_cache_exists = (false || ... || vertex_cache_defined<TCompetence> ());

  /** Helper function to check if class T implements the function centroid_cache_filled.
   * \tparam T The competence to be checked.
   * \return true if T implements the function, false if not.
   */
  template <template <typename> class T>
  static constexpr bool
  centroid_cache_defined ()
  {
    return requires (T<SelfType>& competence) { competence.centroid_cache_filled (); };
  }
  /* This variable is true if any of the given competences \ref TCompetence implements 
  a function centroid_cache_filled. */
  static constexpr bool centroid_cache_exists = (false || ... || centroid_cache_defined<TCompetence> ());

 public:
  /**
   * Constructor for an element of a mesh.
   * \param [in] mesh           Pointer to the mesh the element should belong to.
   * \param [in] tree_id        The tree id of the element in the forest defining the mesh.
   * \param [in] element_id     The element id of the element in the forest defining the mesh.
   */
  element (mesh<SelfType>* mesh, t8_locidx_t tree_id, t8_locidx_t element_id)
    : m_mesh (mesh), m_tree_id (tree_id), m_element_id (element_id)
  {
  }

  // --- Functions to check if caches exist. ---
  /**
   * Function that checks if a cache for the vertex coordinates exists.
   * \return true if a cache for the vertex coordinates exists, false otherwise.
   */
  static constexpr bool
  has_vertex_cache ()
  {
    return vertex_cache_exists;
  }

  /**
   * Function that checks if a cache for the centroid exists.
   * \return true if a cache for the centroid exists, false otherwise.
   */
  static constexpr bool
  has_centroid_cache ()
  {
    return centroid_cache_exists;
  }

  // --- Functionality of the element. In each function, it is checked if a cached version exists (and is used then). ---
  /**
   * Getter for the refinement level of the mesh element.
   * For this easily accessible variable, it makes no sense to provide a cached version.
   * \return Refinement level of the mesh element.
   */
  t8_element_level
  get_level () const
  {
    const t8_eclass_t tree_class = get_tree_class ();
    const t8_element_t* element = get_element ();
    return t8_forest_get_scheme (m_mesh->m_forest)->element_get_level (tree_class, element);
  }

  /**
   * Getter for the vertex coordinates of the mesh element.
   * This function uses or sets the cached version defined in TCompetence if available and calculates if not.
   * \return Vector with one coordinate array for each vertex of the element.
   */
  std::vector<t8_3D_vec>
  get_vertex_coordinates () const
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (vertex_cache_exists) {
      if (this->vertex_cache_filled ()) {
        return this->m_vertex_coordinates;
      }
    }
    // Calculate the vertex coordinates.
    const t8_element_t* element = get_element ();
    const int num_corners
      = t8_forest_get_scheme (m_mesh->m_forest)->element_get_num_corners (get_tree_class (), element);
    std::vector<t8_3D_vec> vertex_coordinates;
    vertex_coordinates.reserve (num_corners);
    for (int icorner = 0; icorner < num_corners; ++icorner) {
      t8_3D_vec vertex;
      t8_forest_element_coordinate (m_mesh->m_forest, m_tree_id, element, icorner, vertex.data ());
      vertex_coordinates.push_back (vertex);
    }
    // Fill the cache in the cached version.
    if constexpr (vertex_cache_exists) {
      this->m_vertex_coordinates = std::move (vertex_coordinates);
      return this->m_vertex_coordinates;
    }
    return vertex_coordinates;
  }

  /**
   * Getter for the center of mass of the mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \return Coordinates of the center.
   */
  t8_3D_vec
  get_centroid () const
  {
    // Check if we have a cached version and if the cache has already been filled.
    if constexpr (centroid_cache_exists) {
      if (this->centroid_cache_filled ()) {
        return this->m_centroid.value ();
      }
    }
    t8_3D_vec coordinates;
    t8_forest_element_centroid (m_mesh->m_forest, m_tree_id, get_element (), coordinates.data ());
    // Fill the cache in the cached version.
    if constexpr (centroid_cache_exists) {
      this->m_centroid = coordinates;
    }
    return coordinates;
  }

  //--- Getter for the member variables. ---
  /**
   * Getter for the tree id of the mesh element.
   * \return The element's tree id.
   */
  t8_locidx_t
  get_tree_id () const
  {
    return m_tree_id;
  }

  /**
   * Getter for the element id of the mesh element.
   * \return The element id of the mesh element.
   */
  t8_locidx_t
  get_element_id () const
  {
    return m_element_id;
  }

  /**
   * Getter for the mesh to which the mesh element is belonging.
   * \return Reference to the mesh.
   */
  const mesh<SelfType>*
  get_mesh () const
  {
    return m_mesh;
  }

 private:
  //--- Private getter for internal use. ---
  /**
   * Getter for the leaf element of the mesh element.
   * \return The leaf element.
   */
  const t8_element_t*
  get_element () const
  {
    return t8_forest_get_leaf_element_in_tree (m_mesh->m_forest, m_tree_id, m_element_id);
  }

  /**
   * Getter for the eclass of the mesh element.
   * \return The element's eclass.
   */
  t8_eclass_t
  get_tree_class () const
  {
    return t8_forest_get_tree_class (m_mesh->m_forest, m_tree_id);
  }

  mesh<SelfType>* m_mesh;   /**< Pointer to the mesh the element is defined for. */
  t8_locidx_t m_tree_id;    /**< The tree id of the element in the forest defined in the mesh. */
  t8_locidx_t m_element_id; /**< The element id of the element in the forest defined in the mesh. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_ELEMENT_HXX */
