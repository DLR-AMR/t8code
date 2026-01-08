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
#include "competence_pack.hxx"
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <type_traits>

namespace t8_mesh_handle
{

/**
 * Wrapper for a forest that enables it to be handled as a simple mesh object.
 * \tparam TCompetence The competences you want to add to the default functionality of the mesh.
 *         \see element for more details on the choice of the template parameter.   
 *         \note Please pack your competences using the \ref competence_pack class.
 * \tparam TUserData The user data type you want to associate with the mesh. Use void (this is also the default) if you do not want to set user data.
 * \tparam TElementData The element data type you want to use for each element of the mesh. 
 *         Use void (this is also the default) if you do not want to set element data.
 */
template <typename TCompetencePack = competence_pack<>, typename TUserData = void, typename TElementData = void>
class mesh {
 public:
  using SelfType = mesh<TCompetencePack, TUserData,
                        TElementData>;  /**< Type of the current class with all template parameters specified. */
  using UserDataType = TUserData;       /**< Make Type of the user data accessible. */
  using ElementDataType = TElementData; /**< Make Type of the element data accessible. */
  using element_class
    = TCompetencePack::template apply<SelfType, element>; /**< The element class of the mesh with given competences. */
  friend element_class; /**< Element class as friend such that private members (e.g. the forest) can be accessed. */
  using mesh_const_iterator =
    typename std::vector<element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */
  using mesh_iterator =
    typename std::vector<element_class>::iterator; /**< Non-const iterator type for the mesh elements. */

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
  cbegin () const
  {
    return m_elements.cbegin ();
  }

  /**
   * Returns a constant iterator to a mesh element following the last (local) element of the mesh.
   * \return Constant iterator to the mesh element following the last (local) element of the mesh.
   */
  mesh_const_iterator
  cend () const
  {
    return m_elements.cend ();
  }

  /**
   * Non-const version of \ref cbegin.
   * \return Iterator to the first (local) mesh element.
   */
  mesh_iterator
  begin ()
  {
    return m_elements.begin ();
  }

  /**
   * Non-const version of \ref end.
   * \return Iterator to the mesh element following the last (local) element of the mesh.
   */
  mesh_iterator
  end ()
  {
    return m_elements.end ();
  }

  /**
   * Getter for an element given its local index. This could be a (local) mesh element or 
   *  a ghost element. 
   * The indices 0, 1, ... num_local_el - 1 refer to local mesh elements and 
   *    num_local_el , ... , num_local_el + num_ghosts - 1 refer to ghost elements.
   * \param [in] local_index The local index of the element to access.
   * \return Constant reference to the element.
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
   * Non const version of operator above.
   * \param [in] local_index The local index of the element to access.
   * \return Reference to the element.
   */
  element_class&
  operator[] (t8_locidx_t local_index)
  {
    return const_cast<element_class&> (static_cast<const mesh*> (this)->operator[] (local_index));
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
    if constexpr (!std::is_void<TElementData>::value) {
      t8_global_infof (
        "The elements of the mesh handle have been updated. Please note that the element data is not interpolated "
        "automatically. Use the function set_element_data() to provide new adapted element data.\n");
    }
  }

  /** 
   * Set the user data of the mesh. This can i.e. be used to pass user defined arguments to the adapt routine.
   * \param [in] data The user data of class TUserData. Data will never be touched by mesh handling routines.
   */
  template <typename U = TUserData, typename = std::enable_if_t<!std::is_void<U>::value>>
  void
  set_user_data (const U& data)
  {
    t8_forest_set_user_data (m_forest, data);
  }

  /** 
   * Get the user data of the mesh. 
   * \return The user data previously set using \ref set_user_data.   
   */
  template <typename U = TUserData, typename = std::enable_if_t<!std::is_void<U>::value>>
  const U&
  get_user_data () const
  {
    return *static_cast<const U*> (t8_forest_get_user_data (m_forest));
  }

  /** 
   * Set the element data vector. The vector should have the length of num_local_elements.
   * \param [in] element_data The element data vector to set with one entry of class TElementData 
   *            for each local mesh element (excluding ghosts).
   */
  template <typename E = TElementData, typename = std::enable_if_t<!std::is_void<E>::value>>
  void
  set_element_data (const std::vector<E>& element_data)
  {
    T8_ASSERT (element_data.size () == get_num_local_elements ());
    m_element_data.clear ();
    m_element_data.resize (get_num_local_elements () + get_num_ghosts ());
    std::copy (element_data.begin (), element_data.end (), m_element_data.begin ());
  }

  /** 
   * Get the element data vector.
   * The element data of the local mesh elements can be set using \ref set_element_data.
   * If ghost entries should be filled, one should call \ref exchange_ghost_data on each process first.
   * \return Element data vector with data of Type TElementData.
   */
  template <typename E = TElementData, typename = std::enable_if_t<!std::is_void<E>::value>>
  const std::vector<E>&
  get_element_data () const
  {
    return m_element_data;
  }

  /** 
  * Exchange the element data for ghost elements between processes.
  * This routine has to be called on each process after setting the element data for all local elements.
  * \return The element data vector of size num_local_elements + num_local_ghosts with data of Type TElementData.
  */
  template <typename E = TElementData, typename = std::enable_if_t<!std::is_void<E>::value>>
  const std::vector<E>&
  exchange_ghost_data ()
  {
    // t8_forest_ghost_exchange_data expects an sc_array, so we need to wrap our data array to one.
    sc_array* sc_array_wrapper;
    T8_ASSERT (m_element_data.size () == get_num_local_elements () + get_num_ghosts ());
    sc_array_wrapper = sc_array_new_data (m_element_data.data (), sizeof (TElementData),
                                          get_num_local_elements () + get_num_ghosts ());

    // Data exchange: entries with indices > num_local_elements will get overwritten.
    t8_forest_ghost_exchange_data (m_forest, sc_array_wrapper);

    sc_array_destroy (sc_array_wrapper);
    return m_element_data;
  }

 private:
  /** 
   * Update the storage of the mesh elements according to the current forest. 
   */
  void
  update_elements ()
  {
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
  std::conditional_t<!std::is_void_v<TElementData>, std::vector<TElementData>, std::nullptr_t>
    m_element_data; /**< Vector storing the (local) element data. */
};

}  // namespace t8_mesh_handle
#endif /* !T8_MESH_HXX */
