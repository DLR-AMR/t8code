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

#pragma once

#include <t8.h>
#include "element.hxx"
#include "competence_pack.hxx"
#include "adapt.hxx"
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <type_traits>
#include <functional>
#include <memory>

namespace t8_mesh_handle
{

/** Concept to ensure that a type is MPI safe.
 */
template <typename TType>
concept T8MPISafeType
  = std::is_void_v<TType> || (std::is_trivially_copyable_v<TType> && std::is_standard_layout_v<TType>);

/**
 * Wrapper for a forest that enables it to be handled as a simple mesh object.
 * \tparam TCompetencePack The competences you want to add to the default functionality of the mesh.
 *         \see element for more details on the choice of the template parameter.   
 *         \note Please pack your competences using the \ref competence_pack class.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 *         Use void (this is also the default) if you do not want to set element data.
 */
template <typename TCompetencePack = competence_pack<>, T8MPISafeType TElementDataType = void>
class mesh {
 public:
  using SelfType
    = mesh<TCompetencePack, TElementDataType>; /**< Type of the current class with all template parameters specified. */
  using ElementDataType = TElementDataType;    /**< Make Type of the element data accessible. */
  using element_class =
    typename TCompetencePack::template apply<SelfType,
                                             element>; /**< The element class of the mesh with given competences. */
  friend element_class; /**< Element class as friend such that private members (e.g. the forest) can be accessed. */
  using mesh_const_iterator =
    typename std::vector<element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */
  using mesh_iterator =
    typename std::vector<element_class>::iterator; /**< Non-const iterator type for the mesh elements. */

  /** Callback function prototype to decide for refining and coarsening of a family of elements
   * or one element in a mesh handle.
   * If \a elements contains more than one element, they must form a family and we decide whether this family should be coarsened
   * or only the first element should be refined.
   * Family means multiple elements that can be coarsened into one parent element.
   * \see set_adapt for the usage of this callback.
   * \param [in] mesh     The mesh that should be adapted.
   * \param [in] elements One element or a family of elements to consider for adaptation.
   * \return 1 if the first entry in \a elements should be refined,
   *        -1 if the family \a elements shall be coarsened,
   *         0 else.
   */
  using adapt_callback_type = std::function<int (const SelfType& mesh, const std::vector<element_class>& elements)>;

  /** Templated callback function prototype to decide for refining and coarsening of a family of elements
   * or one element in a mesh handle including user data.
   * See the version without user_data \ref adapt_callback_type for more details.
   * Use \ref mesh_adapt_callback_wrapper to convert this type into \ref adapt_callback_type 
   * to be able to pass the callback to \ref set_adapt.
   * \tparam TUserDataType The type of the user data to be passed to the callback.
   * \param [in] mesh       The mesh that should be adapted.
   * \param [in] elements   One element or a family of elements to consider for adaptation.
    * \param [in] user_data The user data to be used during the adaptation process.
   * \return 1 if the first entry in \a elements should be refined,
   *        -1 if the family \a elements shall be coarsened,
   *         0 else.
   */
  template <typename TUserDataType>
  using adapt_callback_type_with_userdata
    = std::function<int (const SelfType& mesh, const std::vector<element_class>& elements, TUserDataType user_data)>;

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
   * Destructor for a mesh of the handle. 
   * The forest in use will be unreferenced. 
   * Call \ref t8_forest_ref before if you want to keep it alive.
   */
  ~mesh ()
  {
    t8_forest_unref (&m_forest);
  }

  // --- Getter for mesh related information. ---
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
   * Getter for the forest the mesh is defined for.
   * \return The forest the mesh is defined for.
   */
  t8_forest_t
  get_forest () const
  {
    return m_forest;
  }

  // --- Methods to access elements. ---
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
   * Returns an iterator to the first (local) mesh element.
   * \return Iterator to the first (local) mesh element.
   */
  mesh_iterator
  begin ()
  {
    return m_elements.begin ();
  }

  /**
   * Returns an iterator to a mesh element following the last (local) element of the mesh.
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

  // --- Methods to change the mesh, e.g. adapt, partition, balance, ... ---
  /** Wrapper to convert an adapt callback with user data of type \ref adapt_callback_type_with_userdata
   * into a callback without user data of type \ref adapt_callback_type using the defined user data \a user_data.
   * This is required to pass an adapt callback with user data to \ref set_adapt.
   * \tparam TUserDataType The type of the user data to be passed to the callback.
   * \param [in] adapt_callback_with_userdata The adapt callback including user data.
   * \param [in] user_data The user data to be used during the adaptation process.
   * \return An adapt callback without user data parameter that can be passed to \ref set_adapt.
  */
  template <typename TUserDataType>
  static adapt_callback_type
  mesh_adapt_callback_wrapper (adapt_callback_type_with_userdata<TUserDataType> adapt_callback_with_userdata,
                               const TUserDataType& user_data)
  {
    return [=] (const SelfType& mesh, const std::vector<element_class>& elements) {
      return adapt_callback_with_userdata (mesh, elements, user_data);
    };
  }

  /** Set an adapt function to be used to adapt the mesh on committing.
   * \param [in] adapt_callback    The adapt callback used on committing.
   * \param [in] recursive         Specifying whether adaptation is to be done recursively or not. 
   * \note The adaptation is carried out only when \ref commit is called.
   * \note This setting can be combined with set_partition and set_balance. The order in which
   * these operations are executed is always 1) Adapt 2) Partition 3) Balance.
   */
  void
  set_adapt (adapt_callback_type adapt_callback, bool recursive)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    // Create and register adaptation context holding the mesh handle and the user defined callback.
    detail::AdaptRegistry::register_context (
      m_forest, std::make_unique<detail::MeshAdaptContext<SelfType>> (*this, std::move (adapt_callback)));

    // Set up the forest for adaptation using the wrapper callback.
    t8_forest_set_adapt (m_uncommitted_forest.value (), m_forest, detail::mesh_adapt_callback_wrapper, recursive);
  }

  /** Enable or disable the creation of a layer of ghost elements.
   * \param [in]      do_ghost  If true a ghost layer will be created.
   * \param [in]      ghost_type Controls which neighbors count as ghost elements,
   *                             currently only T8_GHOST_FACES is supported. This value
   *                             is ignored if \a do_ghost = false.
   */
  void
  set_ghost (bool do_ghost = true, t8_ghost_type_t ghost_type = T8_GHOST_FACES)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    t8_forest_set_ghost (m_uncommitted_forest.value (), do_ghost, ghost_type);
  }

  /** After allocating and adding properties to the mesh, commit the changes.
   * This call updates the internal state of the mesh.
   * The forest used to define the mesh handle is replaced in this function.
   * The previous forest is unreferenced. Call \ref t8_forest_ref before if you want to keep it alive.
   * Specialize the update with calls like \ref set_adapt first.
   */
  void
  commit ()
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    /* It can happen that the user only calls set_ghost before commit. 
    This does not set the set_from member of the forest and we copy the current forest in this case. */
    if (m_uncommitted_forest.value ()->set_from == NULL) {
      t8_forest_set_copy (m_uncommitted_forest.value (), m_forest);
    }
    t8_forest_ref (m_forest);
    t8_forest_commit (m_uncommitted_forest.value ());
    // Check if we adapted and unregister the adapt context if so.
    if (detail::AdaptRegistry::get (m_uncommitted_forest.value ()) != nullptr) {
      detail::AdaptRegistry::unregister_context (m_forest);
      if (!std::is_void<TElementDataType>::value) {
        t8_global_infof (
          "Please note that the element data is not interpolated automatically during adaptation. Use the "
          "function set_element_data() to provide new adapted element data.\n");
      }
    }
    t8_forest_unref (&m_forest);
    // Update underlying forest of the mesh.
    m_forest = m_uncommitted_forest.value ();
    m_uncommitted_forest = std::nullopt;
    update_elements ();
  }

  // --- Methods to set and get user and element data and exchange data between processes. ---
  /** 
   * Set the element data vector. The vector should have the length of num_local_elements.
   * \param [in] element_data The element data vector to set with one entry of class TElementDataType 
   *            for each local mesh element (excluding ghosts).
   */
  template <typename ElementDataType = TElementDataType,
            typename = std::enable_if_t<!std::is_void<ElementDataType>::value>>
  void
  set_element_data (std::vector<ElementDataType> element_data)
  {
    T8_ASSERT (element_data.size () == static_cast<size_t> (get_num_local_elements ()));
    m_element_data = std::move (element_data);
    m_element_data.reserve (get_num_local_elements () + get_num_ghosts ());
    m_element_data.resize (get_num_local_elements ());
  }

  /** 
   * Get the element data vector.
   * The element data of the local mesh elements can be set using \ref set_element_data.
   * If ghost entries should be filled, one should call \ref exchange_ghost_data on each process first.
   * \return Element data vector with data of Type TElementDataType.
   */
  template <typename ElementDataType = TElementDataType,
            typename = std::enable_if_t<!std::is_void<ElementDataType>::value>>
  const std::vector<ElementDataType>&
  get_element_data () const
  {
    return m_element_data;
  }

  /** 
  * Exchange the element data for ghost elements between processes.
  * This routine has to be called on each process after setting the element data for all local elements.
  */
  template <typename ElementDataType = TElementDataType,
            typename = std::enable_if_t<!std::is_void<ElementDataType>::value>>
  void
  exchange_ghost_data ()
  {
    // t8_forest_ghost_exchange_data expects an sc_array, so we need to wrap our data array to one.
    sc_array* sc_array_wrapper;
    m_element_data.resize (get_num_local_elements () + get_num_ghosts ());
    sc_array_wrapper = sc_array_new_data (m_element_data.data (), sizeof (ElementDataType),
                                          get_num_local_elements () + get_num_ghosts ());

    // Data exchange: entries with indices > num_local_elements will get overwritten.
    t8_forest_ghost_exchange_data (m_forest, sc_array_wrapper);

    sc_array_destroy (sc_array_wrapper);
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
    if (get_num_ghosts () == 0) {
      m_ghosts.clear ();
      return;
    }
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
  std::conditional_t<!std::is_void_v<TElementDataType>, std::vector<TElementDataType>, std::nullptr_t>
    m_element_data; /**< Vector storing the (local) element data. */
  std::optional<t8_forest_t>
    m_uncommitted_forest; /**< Forest in which the set flags are set for a new forest before committing. */
};

}  // namespace t8_mesh_handle
