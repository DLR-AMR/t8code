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

/** \file mesh.hxx
 * Definition of the mesh class of the handle.
 */

#ifndef T8_MESH_HXX
#define T8_MESH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include "element.hxx"
#include <t8_forest/t8_forest_ghost.h>
#include <vector>

namespace t8_mesh_handle
{

/**
 * Wrapper for a forest that enables it to be handled as a simple mesh object.
 * \tparam TCompetence The competences you want to add to the default functionality of the mesh.
 *         \see element for more details on the choice of the template parameter.   
 */
template <template <typename> class... TCompetence>
class mesh {
 public:
  using element_class = element<TCompetence...>; /**< The element class of the mesh with given competences. */
  friend element_class; /**< Element class as friend such that private members (e.g. the forest) can be accessed. */
  using mesh_const_iterator =
    typename std::vector<element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */

  /** 
   * Constructor for a mesh of the handle. 
   * \param [in] input_forest The forest from which the mesh should be created. 
   */
  mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
    T8_ASSERT (t8_forest_is_committed (m_forest));
    update_elements ();
  }

  /** 
   * Destructor for a mesh of the handle. 
   * The forest in use will be unreferenced. 
   * Call \ref t8_forest_ref before if you want to keep it alive.
   */
  ~mesh ()
  {
    t8_forest_unref (&m_forest);
  }

  /**
   * Getter for the number of local elements in the mesh.
   * \return Number of local elements in the mesh.
   */
  t8_locidx_t
  get_num_local_elements () const
  {
    return t8_forest_get_local_num_leaf_elements (m_forest);
  }

  /**
   * Getter for the number of ghost elements.
   * \return Number of ghost elements in the mesh.
   */
  t8_locidx_t
  get_num_ghosts () const
  {
    return t8_forest_get_num_ghosts (m_forest);
  }

  /** 
   * Getter for the dimension of the mesh.
   * \return The dimension.
   */
  int
  get_dimension () const
  {
    return t8_forest_get_dimension (m_forest);
  }

  /**
   * Returns a constant iterator to the first (local) mesh element.
   * \return Constant iterator to the first (local) mesh element.
   */
  mesh_const_iterator
  begin () const
  {
    return m_elements.cbegin ();
  }

  /**
   * Returns a constant iterator to a mesh element following the last (local) element of the mesh.
   * \return Constant iterator to the mesh element following the last (local) element of the mesh.
   */
  mesh_const_iterator
  end () const
  {
    return m_elements.cend ();
  }

  /**
   * Getter for an element given its local index. This could be a (local) mesh element or 
   *  a ghost element. 
   * The indices 0, 1, ... num_local_el - 1 refer to local mesh elements and 
   *    num_local_el , ... , num_local_el + num_ghosts - 1 refer to ghost elements.
   * \param [in] local_index The local index of the element to access.
   * \return Reference to the element.
   */
  const element_class&
  operator[] (t8_locidx_t local_index) const
  {
    T8_ASSERT (0 <= local_index && local_index < get_num_local_elements () + get_num_ghosts ());
    if (local_index < get_num_local_elements ()) {
      return m_elements[local_index];
    }
    else {
      return m_ghosts[local_index - get_num_local_elements ()];
    }
  }

  /**
   * Getter for the forest the mesh is defined for.
   * \return The forest the mesh is defined for.
   */
  t8_forest_t
  get_forest () const
  {
    return m_forest;
  }

  /** 
  * Setter for the forest. 
  * \param [in] input_forest The forest from which the mesh should be a wrapper. 
  */
  void
  set_forest (t8_forest_t input_forest)
  {
    T8_ASSERT (t8_forest_is_committed (input_forest));
    m_forest = input_forest;
    update_elements ();
  }

 private:
  /** 
   * Update the storage of the mesh elements according to the current forest. 
   */
  void
  update_elements ()
  {
    // Clear the element vector if already created.
    m_elements.clear ();
    m_elements.reserve (get_num_local_elements ());
    // Iterate through forest elements and fill the element vector with newly created mesh elements.
    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (m_forest); ++itree) {
      const t8_locidx_t num_elems = t8_forest_get_tree_num_leaf_elements (m_forest, itree);
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        m_elements.push_back (element_class (this, itree, ielem));
      }
    }
    update_ghost_elements ();
  }

  /** 
   * Update the storage of the ghost elements according to the current forest. 
   */
  void
  update_ghost_elements ()
  {
    // Clear the ghost vector if already created.
    m_ghosts.clear ();
    m_ghosts.reserve (get_num_ghosts ());
    t8_locidx_t num_loc_trees = t8_forest_get_num_local_trees (m_forest);

    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_ghost_trees (m_forest); ++itree) {
      const t8_locidx_t num_elems = t8_forest_ghost_tree_num_leaf_elements (m_forest, itree);
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        m_ghosts.push_back (element_class (this, num_loc_trees + itree, ielem, true));
      }
    }
  }

  t8_forest_t m_forest;                  /**< The forest the mesh should be defined for. */
  std::vector<element_class> m_elements; /**< Vector storing the (local) mesh elements. */
  std::vector<element_class> m_ghosts;   /**< Vector storing the (local) ghost elements. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_MESH_HXX */
