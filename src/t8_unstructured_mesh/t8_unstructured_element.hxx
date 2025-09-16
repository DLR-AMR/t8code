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

/** \file t8_unstructured_element.hxx
 * Definition of the elements used in the unstructured mesh element.
 */

#ifndef T8_UNSTRUCTURED_ELEMENT_HXX
#define T8_UNSTRUCTURED_ELEMENT_HXX
#include <t8.h>

#include <t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <array>
#include <t8_schemes/t8_scheme.hxx>

/* Forward declaration of the unstructured mesh class.
 */
template <class TUnstructuredMeshElement>
class t8_unstructured_mesh;

/** 
 * Unstructured mesh element class. 
 * The unstructured element without specified template parameters provides default implementations for basic functionality 
 * as accessing the refinement level or the centroid. With this implementation, the functionality is calculated each time
 * the function is called. 
 * Use the competences defined in t8_element_competences.hxx as template parameter to cache the functionality instead of 
 * calculating them each time. 
 * To add functionality to the element, you can simply write you own competence class and give it as a template parameter.
 * You can access the functions implemented in your competence via the element. 
 *
 * The inheritance pattern is inspired by the \ref T8Type class.
 * We decided to use this structure 1.) to be able to add new functionality easily and 
 *    2.) for the cached options to keep the number of class member variables of the default to a minimum to safe memory.
 * The choice between calculate and cache is a tradeoff between runtime and memory usage. 
 *
 * \tparam The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class t8_unstructured_mesh_element: public TCompetence<t8_unstructured_mesh_element<TCompetence...>>... {
  using SelfType = t8_unstructured_mesh_element<TCompetence...>;

  // --- Variables to check which functionality is defined in TCompetence. ---
  // Checks if one of the competences (like t8_cache_level) defines the function get_level_cached().
  // Helper function.
  template <template <typename> class T>
  static constexpr bool
  has_get_level_cached ()
  {
    return requires (T<SelfType>& competence) { competence.get_level_cached (); };
  }
  static constexpr bool get_level_defined = (false || ... || has_get_level_cached<TCompetence> ());

  template <template <typename> class T>
  static constexpr bool
  has_get_centroid_cached ()
  {
    return requires (T<SelfType>& competence) { competence.get_centroid_cached (); };
  }
  static constexpr bool get_centroid_defined = (false || ... || has_get_centroid_cached<TCompetence> ());

 public:
  /**
   * Constructor of the unstructured mesh element.
   * \param unstructured_mesh Reference to the unstructured mesh the element should belong to.
   * \param tree_id The tree id of the element in the forest defining the unstructured mesh.
   * \param element_id The element id of the element in the forest defining the unstructured mesh.
   */
  t8_unstructured_mesh_element (t8_unstructured_mesh<SelfType>* unstructured_mesh, t8_locidx_t tree_id,
                                t8_locidx_t element_id)
    : m_tree_id (tree_id), m_element_id (element_id), m_unstructured_mesh (unstructured_mesh)
  {
  }

  // --- Functionality of the element. In each function, it is checked if a cached version exists (and is used then). ---

  /**
   * Getter for the refinement level of the unstructured mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates the refinement level if not.
   * \return Refinement level of the unstructured mesh element.
   */
  t8_element_level
  get_level ()
  {
    if constexpr (get_level_defined) {
      return this->get_level_cached ();
    }
    else {
      const t8_eclass_t tree_class = t8_forest_get_tree_class (m_unstructured_mesh->m_forest, m_tree_id);
      const t8_element_t* element
        = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
      return t8_forest_get_scheme (m_unstructured_mesh->m_forest)->element_get_level (tree_class, element);
    }
  }

  /**
   * Getter for the center of mass of the unstructured mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \return Coordinates of the center.
   */
  std::array<double, T8_ECLASS_MAX_DIM>
  get_centroid ()
  {
    if constexpr (get_centroid_defined) {
      return this->get_centroid_cached ();
    }
    else {
      std::array<double, T8_ECLASS_MAX_DIM> coordinates;
      const t8_element_t* element
        = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
      t8_forest_element_centroid (m_unstructured_mesh->m_forest, m_tree_id, element, coordinates.data ());
      return coordinates;
    }
  }

  //--- Getter for the member variables. ---
  /**
   * Getter for the tree id of the unstructured mesh element.
   */
  t8_locidx_t
  get_tree_id ()
  {
    return m_tree_id;
  }

  /**
   * Getter for the element id of the unstructured mesh element.
   */
  t8_locidx_t
  get_element_id ()
  {
    return m_element_id;
  }

  /**
   * Getter for the unstructured mesh to which the unstructured mesh element is belonging.
   * \return Reference to the unstructured mesh.
   */
  t8_unstructured_mesh<SelfType>*
  get_unstructured_mesh ()
  {
    return m_unstructured_mesh;
  }

 private:
  t8_locidx_t m_tree_id,
    m_element_id; /**< The tree id and the element id of the element in the forest defined in the unstructured mesh. */
  t8_unstructured_mesh<SelfType>*
    m_unstructured_mesh; /**< Pointer to the unstructured mesh the element is defined for. */
};

#endif /* !T8_UNSTRUCTURED_ELEMENT_HXX */
