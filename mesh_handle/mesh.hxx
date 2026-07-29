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

/** \file mesh.hxx
 * Definition of the mesh class of the handle.
 */

#pragma once

#include <t8.h>
#include "element.hxx"
#include "competence_pack.hxx"
#include "internal/adapt.hxx"
#include "internal/interpolate.hxx"
#include "competences/element_data_competences.hxx"
#include "concepts.hxx"
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_iterate.h>
#include <vector>
#include <functional>
#include <memory>
#include <span>

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
 *         \note Please pack your competences using the \ref t8_mesh_handle::mesh_competence_pack class.
 *         One of the most important competences to add is \ref element_data_mesh_competence.
 */
template <ElementCompetencePack TElementCompetencePack = element_competence_pack<>,
          MeshCompetencePack TMeshCompetencePack = mesh_competence_pack<>>
class mesh: public TMeshCompetencePack::template apply<mesh<TElementCompetencePack, TMeshCompetencePack>> {
 public:
  using SelfType = mesh<TElementCompetencePack, TMeshCompetencePack>; /**< Type of the current class. */
  using element_class = typename TElementCompetencePack::template apply<
    SelfType, element>;  /**< The element class of the mesh with given competences. */
  friend element_class;  /**< Element class as friend such that private members (e.g. the forest) can be accessed. */
  using mesh_tag = void; /**< Mesh tag for identification in concept. */
  using mesh_const_iterator =
    typename std::vector<element_class>::const_iterator; /**< Constant iterator type for the mesh elements. */
  using mesh_iterator =
    typename std::vector<element_class>::iterator;              /**< Non-const iterator type for the mesh elements. */
  friend struct element_data_element_competence<element_class>; /**< Friend struct to access its element data vector. */

  // --- Definition of callback types. ---
  /** Callback function prototype to decide for refining and coarsening of a family of elements
   * or one element in a mesh handle.
   * If \a elements contains more than one element, they must form a family and we decide whether this family should be
   * coarsened or only the first element should be refined.
   * Family means multiple elements that can be coarsened into one parent element.
   * \see set_adapt for the usage of this callback.
   * \param [in] mesh     The mesh that should be adapted.
   * \param [in] elements One element or a family of elements (if more than one element) to consider for adaptation.
   * \return 1 if the first entry in \a elements should be refined,
   *        -1 if the family \a elements shall be coarsened,
   *         0 else.
   * \note We currently do not provide functionality to delete elements.
   */
  using adapt_callback_type = std::function<int (const SelfType& mesh, std::span<const element_class> elements)>;

  /** Templated callback function prototype to decide for refining and coarsening of a family of elements
   * or one element in a mesh handle including user data.
   * See the version without user_data \ref adapt_callback_type for more details.
   * Use \ref mesh_adapt_callback_wrapper to convert this type into \ref adapt_callback_type 
   * to be able to pass the callback to \ref set_adapt.
   * \tparam TUserDataType The type of the user data to be passed to the callback.
   * \param [in] mesh       The mesh that should be adapted.
   * \param [in] elements One element or a family of elements (if more than one element) to consider for adaptation.
   * \param [in] user_data The user data to be used during the adaptation process.
   * \return 1 if the first entry in \a elements should be refined,
   *        -1 if the family \a elements shall be coarsened,
   *         0 else.
   * \note We currently do not provide functionality to delete elements.
   */
  template <typename TUserDataType>
  using adapt_callback_type_with_userdata
    = std::function<int (const SelfType& mesh, std::span<const element_class> elements, TUserDataType user_data)>;

  /** Callback function prototype to interpolate the element data after refining or coarsening.
   * \note You need to include \ref interpolate_element_data_mesh_competence to you competences to be able to
   * interpolate. The best way to do this is via the predefined pack \ref interpolate_data_mesh_competence_pack 
   * defined in \ref competence_pack.hxx.
   *
   * For each group of elements that changed during adaption, the outgoing elements of the old mesh are passed in 
   * \a old_elements and the incoming elements of the new mesh in \a new_elements; the callback reads the old data 
   * and writes the interpolated data onto the new elements. \a refine is the value \ref adapt_callback_type returned
   * for this group.
   * \see set_interpolate_callback for the usage of this callback.
   * \param [in]     mesh_old     The old mesh that is adapted from.
   * \param [in,out] mesh_new     The new mesh constructed from \a mesh_old.
   * \param [in]     refine       -1 if the family \a old_elements got coarsened, 0 if the element was not touched,
   *                            1 if the element got refined. Same convention as the return of \ref adapt_callback_type.
   * \param [in]     old_elements Span over the outgoing elements: the whole family on coarsening, 
   *                              a single element if refined or untouched.
   * \param [in,out] new_elements Span over the incoming elements to write the interpolated data to: the children on
   *                              refinement, a single element if coarsened or untouched.
   */
  using interpolate_callback_type
    = std::function<void (const SelfType& mesh_old, SelfType& mesh_new, const int refine,
                          std::span<const element_class> old_elements, std::span<element_class> new_elements)>;

  /** Templated callback function prototype to interpolate the element data after refining or coarsening, 
   * including user data.
   * See the version without user_data \ref interpolate_callback_type for more details!
   * Use \ref mesh_interpolate_callback_wrapper to convert this type into \ref interpolate_callback_type
   * to be able to pass the callback to \ref set_interpolate_callback (see \ref element_data_competence.hxx).
   * \tparam TUserDataType The type of the user data to be passed to the callback.
   * \param [in]     mesh_old     The old mesh that is adapted from.
   * \param [in,out] mesh_new     The new mesh constructed from \a mesh_old.
   * \param [in]     refine       -1 if the family got coarsened, 0 if the element was not touched, 1 if it got refined.
   * \param [in]     old_elements Span over the outgoing elements from \a mesh_old.
   * \param [in,out] new_elements Span over the incoming elements to write the interpolated data to from \a mesh_new.
   * \param [in]     user_data    The user data to be used during the interpolation.
   */
  template <typename TUserDataType>
  using interpolate_callback_type_with_userdata = std::function<void (
    const SelfType& mesh_old, SelfType& mesh_new, const int refine, std::span<const element_class> old_elements,
    std::span<element_class> new_elements, TUserDataType user_data)>;

  // --- Constructor and destructor. ---
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
  * The mesh is said to be balanced if the level difference between face neighbors is at most 1.
  * at most +1 or -1 of the element's level.
  * \return true if the local elements are balanced, false otherwise.
  */
  bool
  is_balanced () const
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
   * Returns a constant iterator to the first (local) mesh element.
   * \return Constant iterator to the first (local) mesh element.
   */
  mesh_const_iterator
  begin () const
  {
    return this->cbegin ();
  }

  /**
   * Returns a constant iterator to a mesh element following the last (local) element of the mesh.
   * \return Constant iterator to the mesh element following the last (local) element of the mesh.
   */
  mesh_const_iterator
  end () const
  {
    return this->cend ();
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
  /** Wrapper to convert an interpolate callback with user data of type \ref interpolate_callback_type_with_userdata
   * into a callback without user data of type \ref interpolate_callback_type using the defined user data \a user_data.
   * The returned callback can then be passed to \ref set_interpolate_callback.
   * See also \ref element_data_competences.hxx for the interpolation competence.
   * \tparam TUserDataType The type of the user data to be passed to the callback.
   * \param [in] interpolate_callback_with_userdata The interpolate callback including user data.
   * \param [in] user_data The user data to be used during the interpolation process.
   * \return An interpolate callback without user data parameter that can be passed to \ref set_interpolate_callback.
   */
  template <typename TUserDataType>
  static interpolate_callback_type
  mesh_interpolate_callback_wrapper (
    interpolate_callback_type_with_userdata<TUserDataType> interpolate_callback_with_userdata,
    const TUserDataType& user_data)
  {
    return [=] (const SelfType& mesh_old, SelfType& mesh_new, const int refine,
                std::span<const element_class> old_elements, std::span<element_class> new_elements) {
      return interpolate_callback_with_userdata (mesh_old, mesh_new, refine, old_elements, new_elements, user_data);
    };
  }

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
    return [=] (const SelfType& mesh, std::span<const element_class> elements) {
      return adapt_callback_with_userdata (mesh, elements, user_data);
    };
  }

  /** Set an adapt function to be used to adapt the mesh on committing.
   * \param [in] adapt_callback The adapt callback used on committing.
   * \note The adaptation is carried out only when \ref commit is called.
   * \note We currently do not provide the functionality to delete elements.
   * \note This setting can be combined with set_partition and set_balance. The order in which
   * these operations are executed is always 1) Adapt 2) Balance 3) Partition.
   */
  void
  set_adapt (adapt_callback_type adapt_callback)
  {
    SC_CHECK_ABORT (m_forest->incomplete_trees == 0, "The mesh handle can't adapt forests with incomplete trees.\n");
    if (!m_uncommitted_forest.has_value ()) {
      m_uncommitted_forest.emplace ();
      t8_forest_init (&*m_uncommitted_forest);
    }
    // Create and register adaptation context holding the mesh handle and the user defined callback.
    detail::adapt_registry::register_context (
      m_forest, std::make_unique<detail::mesh_adapt_context<SelfType>> (*this, std::move (adapt_callback)));

    // Set up the forest for adaptation using the wrapper callback.
    // Recursive adaptation is currently not supported.
    t8_forest_set_adapt (m_uncommitted_forest.value (), m_forest, detail::mesh_adapt_callback_wrapper, false);
  }

  /** If this function is called, the mesh will be partitioned on committing.
   * The partitioning is done according to the SFC and each rank is assigned
   * the same (maybe +1) number of elements.
   * \note The partition is carried out only when \ref commit is called.
   * \note This setting can be combined with \ref set_adapt and \ref set_balance. The order in which
   * these operations are executed is always 1) Adapt 2) Balance 3) Partition.
   * \param [in] set_for_coarsening If true, the partitions are chosen such that coarsening 
   *        an element once is a process local operation. Default is false.
   */
  void
  set_partition (bool set_for_coarsening = false)
  {
    // If the mesh has an interpolate callback, we partition the mesh after the interpolation step
    // (and the first committing of the forest), such that we store the partition information for later.
    if constexpr (has_interpolate_data_competence ()) {
      this->m_partition_for_coarsening = set_for_coarsening;
      return;
    }
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
    t8_forest_set_partition (m_uncommitted_forest.value (), m_forest, set_for_coarsening);
  }

  /** If this function is called, the mesh will be balanced on committing.
   * The mesh is said to be balanced if the element level between face neighbors differs by at most 1.
   * \note The balance is carried out only when \ref commit is called.
   * \param [in] no_repartition Balance constructs several intermediate steps that
   *       are refined from each other. In order to maintain a balanced load, a repartitioning is performed in each 
   *       round and the resulting mesh is load-balanced per default. 
   *       Set \a no_repartition to true if this behaviour is not desired.
   *       If \a no_repartition is false (default), an additional call of \ref set_partition is not necessary.
   * \note This setting can be combined with \ref set_adapt and \ref set_partition. The order in which
   * these operations are executed is always 1) Adapt 2) Balance 3) Partition.
   */
  void
  set_balance (bool no_repartition = false)
  {
    if (!m_uncommitted_forest.has_value ()) {
      t8_forest_t new_forest;
      t8_forest_init (&new_forest);
      m_uncommitted_forest = new_forest;
    }
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
      m_uncommitted_forest.emplace ();
      t8_forest_init (&*m_uncommitted_forest);
    }
    t8_forest_set_ghost (m_uncommitted_forest.value (), do_ghost, ghost_type);
  }

  /** After allocating and adding properties to the mesh, commit the changes.
   * This call updates the internal state of the mesh.
   * The forest used to define the mesh handle is replaced in this function.
   * The previous forest is unreferenced. Call \ref t8_forest_ref before if you want to keep it alive.
   * Specialize the update with calls like \ref set_adapt first.
   * The order of the calls is always 1) Adapt 2) Balance 3) Data Interpolation 4) Partition 5) Ghost, 
   * where calls, that are not set beforehand, are skipped.
   * The order of the calls does not matter, the operations are always executed in this order.
   */
  void
  commit ()
  {
    if (!m_uncommitted_forest.has_value ()) {
      m_uncommitted_forest.emplace ();
      t8_forest_init (&*m_uncommitted_forest);
    }
    /* It can happen that the user only calls set_ghost before commit. 
    This does not set the set_from member of the forest and we copy the current forest in this case. */
    if (m_uncommitted_forest.value ()->set_from == NULL) {
      t8_forest_set_copy (m_uncommitted_forest.value (), m_forest);
    }
    t8_forest_ref (m_forest);
    t8_forest_commit (m_uncommitted_forest.value ());
    // Check if we adapted and unregister the adapt context if so.
    if (detail::adapt_registry::get (m_forest) != nullptr) {
      detail::adapt_registry::unregister_context (m_forest);

      // If data are set for the mesh, we now try to interpolate them after adaptation.
      if constexpr (has_element_data_handler_competence ()) {
        if constexpr (has_interpolate_data_competence ()) {
          if (this->m_interpolate_callback) {
            // Create new intermediate mesh to interpolate the data from the current mesh to the new mesh.
            SelfType new_mesh (m_uncommitted_forest.value ());
            t8_forest_ref (m_uncommitted_forest.value ());
            // Register the interpolate context with the callback for the new mesh. With this, the standard
            // iterate replace can be called.
            detail::interpolate_registry::register_context (
              m_forest, std::make_unique<detail::mesh_interpolate_context<SelfType>> (
                          *this, new_mesh, std::move (this->m_interpolate_callback)));
            t8_forest_iterate_replace (m_uncommitted_forest.value (), m_forest, detail::mesh_replace_callback_wrapper);
            detail::interpolate_registry::unregister_context (m_forest);
            // Override the element data of the current mesh with the interpolated data from the "new mesh".
            this->m_element_data = new_mesh.take_element_data ();
            // Now we update the forest of the current mesh with the new forest and partition it if required.
            t8_forest_unref (&m_forest);
            if (this->set_partition_called ()) {
              t8_forest_init (&m_forest);
              t8_forest_set_partition (m_forest, m_uncommitted_forest.value (),
                                       this->m_partition_for_coarsening.value ());
              if (t8_forest_get_num_ghosts (m_uncommitted_forest.value ()) > 0) {
                t8_forest_set_ghost (m_forest, true, T8_GHOST_FACES);
              }
              t8_forest_commit (m_forest);

              /* Now we repartition also the data: The interpolated data follows m_uncommitted_forest. 
               * We align it now with the partitioned m_forest. */
              this->repartition_element_data (m_uncommitted_forest.value (), m_forest);
            }
            else {
              // Update underlying forest of the mesh for the case where we do not repartition.
              m_forest = m_uncommitted_forest.value ();
            }
            // Cleanup and update the elements of the mesh.
            m_uncommitted_forest.reset ();
            update_elements ();
            return;
          }
          else {
            t8_global_infof ("No interpolation context set.\n");
          }
        }
        else {
          t8_global_infof ("The element data was not interpolated during adaptation. Use set_element_data() to provide "
                           "new data or use the mesh competence interpolate_element_data_mesh_competence.\n");
        }
      }
    }
    t8_forest_unref (&m_forest);
    // Update underlying forest of the mesh.
    m_forest = m_uncommitted_forest.value ();
    m_uncommitted_forest.reset ();
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

  /** Function that checks if a competence for the interpolation of element data is given.
   * \return true if mesh has the competence, false otherwise.
   */
  static constexpr bool
  has_interpolate_data_competence ()
  {
    return requires (SelfType& mesh) { mesh.set_partition_called (); };
  }

  /** Function that checks if a competence to determine the ranks of the elements is given.
   * \return true if mesh has the competence, false otherwise.
   */
  static constexpr bool
  has_remote_ranks_mesh_competence ()
  {
    return requires (SelfType& mesh) { mesh.fill_rank_vector (); };
  }

  /** Function that checks if a competence to determine a unique vector of the faces is given.
   * \return true if mesh has the competence, false otherwise.
   */
  static constexpr bool
  has_face_vector_mesh_competence ()
  {
    return requires (SelfType& mesh) { mesh.fill_unique_face_vector (); };
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
