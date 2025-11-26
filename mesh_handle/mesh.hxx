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
#include "abstract_element.hxx"
#include "mesh_element.hxx"
#include "ghost_element.hxx"
#include "competence_pack.hxx"
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <type_traits>

namespace t8_mesh_handle
{

/**
 * Wrapper for a forest that enables it to be handled as a simple mesh object.
 * \tparam TCompetencePack The competences you want to add to the default functionality of the mesh.
 *         \see abstract_element for more details on the choice of the template parameter.   
 *         \note Please pack your competences using the \ref competence_pack class.
 */
template <typename TCompetencePack = competence_pack<>, typename TUserData = void>
class mesh {
 public:
  using SelfType = mesh<TCompetencePack, TUserData>;
  /** Type definitions of the element classes with given competences. */
  using abstract_element_class = TCompetencePack::template apply<
    SelfType, abstract_element>; /**< The abstract element class of the mesh (could be a mesh element of ghost). */
  using mesh_element_class
    = TCompetencePack::template apply<SelfType, mesh_element>; /**< The mesh element class of the mesh. */
  using ghost_element_class
    = TCompetencePack::template apply<SelfType, ghost_element>; /**< The ghost element class of the mesh. */

  // Declare all element classes as friend such that private members (e.g. the forest) can be accessed.
  friend abstract_element_class;
  friend mesh_element_class;
  friend ghost_element_class;

  using mesh_const_iterator =
    typename std::vector<mesh_element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */

  /** 
   * Constructor for a mesh of the handle. 
   * \param [in] forest The forest from which the mesh should be created. 
   */
  mesh (t8_forest_t forest): m_forest (forest)
  {
    T8_ASSERT ((std::is_same<typename TCompetencePack::is_competence_pack, void>::value));
    T8_ASSERT (t8_forest_is_committed (m_forest));
    update_elements ();
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
   * Getter for the number of local ghost elements.
   * \return Number of local ghost elements in the mesh.
   */
  t8_locidx_t
  get_num_local_ghosts () const
  {
    return t8_forest_get_num_ghosts (m_forest);
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
  const abstract_element_class&
  operator[] (t8_locidx_t local_index) const
  {
    T8_ASSERT (0 <= local_index && local_index < get_num_local_elements () + get_num_local_ghosts ());
    if (local_index < get_num_local_elements ()) {
      return m_elements[local_index];
    }
    else {
      return m_ghosts[local_index - get_num_local_elements ()];
    }
  }

  /**
   * Getter for a mesh element given its local index. 
   * This function could be used instead of the operator[] if a mesh element is 
   * specifically required to access their functionality.
   * \param [in] local_index The local index of the element to access.
   * \return Reference to the mesh element.
   */
  const mesh_element_class&
  get_mesh_element (t8_locidx_t local_index) const
  {
    T8_ASSERT (0 <= local_index && local_index < get_num_local_elements ());
    return m_elements[local_index];
  }

  /**
   * Getter for a ghost element given its local flat index. 
   * This function could be used instead of the operator[] if a ghost element is 
   * specifically required to access ghost functionality.
   * \param [in] flat_index The local flat index of the element to access.
   *             Therefore, flat_index should lie in between num_local_el and
   *              num_local_el + num_ghosts - 1.
   * \return Reference to the ghost element.
   */
  const ghost_element_class&
  get_ghost_element (t8_locidx_t flat_index) const
  {
    T8_ASSERT (get_num_local_elements () <= flat_index
               && flat_index < get_num_local_elements () + get_num_local_ghosts ());
    return m_ghosts[flat_index - get_num_local_elements ()];
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
    m_forest = input_forest;
    T8_ASSERT (t8_forest_is_committed (m_forest));
    update_elements ();
  }

  /** 
    * Set the user data of the mesh. This can i.e. be used to pass user defined arguments to the adapt routine.
    * \param [in] data The user data. Data will never be touched by mesh handling routines.
    */
  template <typename U = TUserData, typename = std::enable_if_t<!std::is_void<U>::value>>
  void
  set_user_data (const U& data)
  {
    t8_forest_set_user_data (m_forest, data);
  }

  /** 
    * TODO
    */
  template <typename U = TUserData, typename = std::enable_if_t<!std::is_void<U>::value>>
  const U&
  get_user_data ()
  {
    return *static_cast<const U*> (t8_forest_get_user_data (m_forest));
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
        m_elements.push_back (mesh_element_class (this, itree, ielem));
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
    m_ghosts.reserve (get_num_local_ghosts ());
    t8_locidx_t num_loc_trees = t8_forest_get_num_local_trees (m_forest);

    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_ghost_trees (m_forest); ++itree) {
      const t8_locidx_t num_elems = t8_forest_ghost_tree_num_leaf_elements (m_forest, itree);
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        m_ghosts.push_back (ghost_element_class (this, num_loc_trees + itree, itree, ielem));
      }
    }
  }

  t8_forest_t m_forest;                       /**< The forest the mesh should be defined for. */
  std::vector<mesh_element_class> m_elements; /**< Vector storing the (local) mesh elements. */
  std::vector<ghost_element_class> m_ghosts;  /**< Vector storing the (local) ghost elements. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_MESH_HXX */
