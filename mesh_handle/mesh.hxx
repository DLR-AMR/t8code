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
#include "competences.hxx"
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <functional>
#include <memory>

namespace t8_mesh_handle
{

/** Concept to ensure that a type is an element competence pack.
 */
template <typename TType>
concept ElementCompetencePack = requires { typename TType::is_element_competence_pack; };
/** Concept to ensure that a type is a mesh competence pack.
 */
template <typename TType>
concept MeshCompetencePack = requires { typename TType::is_mesh_competence_pack; };

/**
 * Wrapper for a forest that enables it to be handled as a simple mesh object.
 * \tparam TElementCompetencePack The competences you want to add to the default functionality of the elements.
 *         \see element for more details on the choice of the template parameter.   
 *         \note Please pack your competences using the \ref element_competence_pack class.
 * \tparam TMeshCompetences The competences you want to add to the default functionality of the mesh.  
 *         \note Please pack your competences using the \ref mesh_competence_pack class.
 *         One of the most important competences to add is \ref handle_element_data.
 */
template <ElementCompetencePack TElementCompetencePack = element_competence_pack<>,
          MeshCompetencePack TMeshCompetencePack = mesh_competence_pack<>>
class mesh: public TMeshCompetencePack::template apply<mesh<TElementCompetencePack, TMeshCompetencePack>> {
 public:
  using SelfType = mesh<TElementCompetencePack, TMeshCompetencePack>; /**< Type of the current class. */
  using element_class
    = TElementCompetencePack::template apply<SelfType,
                                             element>; /**< The element class of the mesh with given competences. */
  friend element_class; /**< Element class as friend such that private members (e.g. the forest) can be accessed. */
  using mesh_const_iterator =
    typename std::vector<element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */
  using mesh_iterator =
    typename std::vector<element_class>::iterator;  /**< Non-const iterator type for the mesh elements. */
  friend struct access_element_data<element_class>; /**< Friend struct to access its element data vector. */

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

  /** Check if the local elements of the mesh are balanced. 
  * The mesh is said to be balanced if each element has face neighbors of level
  * at most +1 or -1 of the element's level.
  * \return true if the local elements are balanced, false otherwise.
  */
  bool
  is_balanced ()
  {
    return t8_forest_is_balanced (m_forest);
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

  /** If this function is called, the mesh will be partitioned on committing.
   * The partitioning is done according to the SFC and each rank is assigned
   * the same (maybe +1) number of elements.
   * \note The partition is carried out only when \ref commit is called.
   * \note This setting can be combined with \ref set_adapt and \ref set_balance. The order in which
   * these operations are executed is always 1) Adapt 2) Partition 3) Balance.
   * \param [in] set_for_coarsening If true, the partitions are choose such that coarsening 
   *        an element once is a process local operation. Default is false.
   */
  void
  set_partition (bool set_for_coarsening = false)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    t8_forest_set_partition (m_uncommitted_forest.value (), m_forest, set_for_coarsening);
  }

  /** If this function is called, the mesh will be balanced on committing.
   * The mesh is said to be balanced if each element has face neighbors of level
   * at most +1 or -1 of the element's level.
   * \note The balance is carried out only when \ref commit is called.
   * \param [in] no_repartition Balance constructs several intermediate steps that
   *       are refined from each other. In order to maintain a balanced load, a repartitioning is performed in each 
   *       round and the resulting mesh is load-balanced per default. 
   *       Set \a no_repartition to true if this behaviour is not desired.
   *       If \a no_repartition is false (default), an additional call of \ref set_partition is not necessary.
   * \note This setting can be combined with \ref set_adapt and \ref set_partition. The order in which
   * these operations are executed is always 1) Adapt 2) Partition 3) Balance.
   */
  void
  set_balance (bool no_repartition = false)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    // Disable repartitioning and let the user call set_partition if desired.
    t8_forest_set_balance (m_uncommitted_forest.value (), m_forest, no_repartition);
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
      if (has_element_data_handler_competence ()) {
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

  // --- Methods to check for mesh competences. ---
  /** Function that checks if a competence for element data handling is given.
   * \return true if mesh has a data handler, false otherwise.
   */
  static constexpr bool
  has_element_data_handler_competence ()
  {
    return requires (SelfType& mesh) { mesh.get_element_data (); };
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
  std::optional<t8_forest_t>
    m_uncommitted_forest; /**< Forest in which the set flags are set for a new forest before committing. */
};

}  // namespace t8_mesh_handle
