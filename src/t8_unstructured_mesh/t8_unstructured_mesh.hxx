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

/** \file t8_unstructured_mesh.hxx
 * Definition of the unstructured mesh class.
 */

#ifndef T8_UNSTRUCTURED_MESH_HXX
#define T8_UNSTRUCTURED_MESH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_unstructured_mesh/t8_unstructured_element.hxx>
#include <iterator>
#include <memory>
#include <vector>

/**
 * Wrapper for a forest that enables it to be handled like an unstructured mesh object.
 * \tparam TUnstructuredMeshElement The element class that should be used for the unstructured mesh elements. 
 *                                  This template parameter defines which element functionality is available and if it is cached or calculated.
 */
template <class TUnstructuredMeshElement = t8_unstructured_mesh_element<>>
class t8_unstructured_mesh {
 public:
  // Declare unstructured mesh element as friend such that private members (e.g. the forest) can be accessed.
  friend TUnstructuredMeshElement;

  using t8_unstructured_iterator = typename std::vector<TUnstructuredMeshElement>::iterator;

  /** 
   * Constructor for an unstructured mesh. 
   * \param [in] input_forest The forest from which the unstructured mesh should be created. 
   */
  t8_unstructured_mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
    update_elements ();
  }

  /**
  * Getter for the number of local elements in the unstructured mesh.
  * \return Number of local elements in the unstructured mesh.
  */
  t8_locidx_t
  get_local_num_elements () const
  {
    return t8_forest_get_local_num_leaf_elements (m_forest);
  }

  /**
   * Returns an iterator to the first (local) unstructured mesh element.
   * \return Iterator to the first (local) unstructured mesh element.
   */
  t8_unstructured_iterator
  begin ()
  {
    return m_elements.begin ();
  }

  /**
   * Returns an iterator to an unstructured mesh element following the last (local) element of the unstructured mesh.
   * \return Iterator to the mesh element following the last (local) element of the unstructured mesh.
   */
  t8_unstructured_iterator
  end ()
  {
    return m_elements.end ();
  }

  /**
   * Getter for an unstructured mesh element given its local index.
   *\param[in] local_index The local index of the element to access.
   * \return Reference to the unstructured mesh element.
   */
  TUnstructuredMeshElement&
  operator[] (t8_locidx_t local_index)
  {
    return m_elements[local_index];
  }

  /**
   * Getter for the forest the unstructured mesh is defined for.
   * \return The forest  the unstructured mesh is defined for.
   */
  t8_forest_t
  get_forest () const
  {
    return m_forest;
  }

  /** 
  * Setter for the forest. 
  * \param [in] input_forest The forest from which the unstructured mesh should be a wrapper. 
  */
  void
  set_forest (t8_forest_t input_forest)
  {
    m_forest = input_forest;
    update_elements ();
  }

 private:
  /** 
   * Update the storage of the unstructured mesh elements according to the current forest. 
   */
  void
  update_elements ()
  {
    // Clear the element vector if already created.
    if (!m_elements.empty ()) {
      m_elements.clear ();
    }
    m_elements.reserve (get_local_num_elements ());
    // Iterate through forest elements and fill the element vector with newly created unstructured mesh elements.
    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (m_forest); ++itree) {
      const t8_locidx_t num_elems = t8_forest_get_tree_num_leaf_elements (m_forest, itree);
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        m_elements.emplace_back (this, itree, ielem);
      }
    }
  }

  t8_forest_t m_forest;                             /**< The forest the unstructured mesh should be defined for. */
  std::vector<TUnstructuredMeshElement> m_elements; /**< Vector storing the unstructured mesh elements. */
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */
