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

/** \file element_data_competences.hxx
 * Handler for the element data of a \ref t8_mesh_handle::mesh.
 * The file defines mesh and element competences for element data handling.
 * The mesh competences make it possible to manage element data and exchange it for ghost elements between processes. 
 * The element competences makes it possible to access these element data directly for each element of the mesh.
 * A competence to interpolate data after adaptation using a user defined callback is provided.
 */
#pragma once

#include <t8.h>
#include <mesh_handle/concepts.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_types/t8_crtp.hxx>
#include <t8_types/t8_operators.hxx>
#include <type_traits>
#include <vector>
#include <functional>
#include <optional>
#include <span>

namespace t8_mesh_handle
{

/** Namespace detail to hide implementation details from the user. */
namespace detail
{
/** Helper function to wrap a span based interpolation callback (see \ref mesh::interpolate_callback_type) into
 * the element-index based \ref interpolate_element_data_mesh_competence::internal_interpolate_callback_type.
 * The returned wrapper receives the index/count pairs, builds the element spans, and forwards them to \a callback. 
 * The spans are built here and not in \ref interpolate_element_data_mesh_competence because \a TMesh is complete, 
 * so element_class is nameable — which it is not inside the competence (see note on 
 * \ref interpolate_element_data_mesh_competence::internal_interpolate_callback_type).
 * This is used in \ref interpolate_element_data_mesh_competence::set_interpolate_callback
 * \tparam TMesh The (complete) mesh handle type.
 * \param [in] callback The span based user callback of type \ref mesh::interpolate_callback_type. Taken by value and
 *                      moved into the returned wrapper, which owns it.
 * \return Callback of type \ref interpolate_element_data_mesh_competence::internal_interpolate_callback_type.
 */
template <T8MeshType TMesh>
auto
to_replace_callback (typename TMesh::interpolate_callback_type callback)
{
  return
    [callback = std::move (callback)] (const TMesh& mesh_old, TMesh& mesh_new, const int refine, const int num_old,
                                       const t8_locidx_t first_old, const int num_new, const t8_locidx_t first_new) {
      callback (mesh_old, mesh_new, refine, std::span (&mesh_old[first_old], num_old),
                std::span (&mesh_new[first_new], num_new));
    };
}
}  // namespace detail

// --- Mesh competence for element data management. ---
/** Handler for the element data of a \ref mesh.
 * Use this competence if you want to manage element data for the elements of the mesh.
 * Use the helper \ref element_data_mesh_competence to get this competence with the correct template parameters form. 
 * If you want to access the data not only in vector form but also directly for each element, 
 * you can combine this competence with \ref element_data_element_competence.
 * In summary you can use the competences like this: 
 *    mesh<element_competence_pack<element_data_element_competence>,
 *         mesh_competence_pack<element_data_mesh_competence<YourElementDataType>::template type>>;
 * Some predefined competences are also defined in \ref competence_pack.hxx.
 *
 * \tparam TUnderlying Use the \ref mesh class here.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class element_data_mesh_competence_impl: public t8_crtp_basic<TUnderlying> {
 public:
  using ElementDataType = TElementDataType; /**< Make Type of the element data publicly accessible. */

  /** Set the element data vector. The vector should have the length of num_local_elements.
   * \param [in] element_data The element data vector to set with one entry of class TElementDataType 
   *            for each local mesh element (excluding ghosts). The vector is moved, not copied.
   */
  void
  set_element_data (std::vector<TElementDataType> element_data)
  {
    T8_ASSERT (element_data.size () == static_cast<size_t> (this->underlying ().get_num_local_elements ()));
    // element_data is moved, not copied.
    m_element_data = std::move (element_data);
  }

  /** Get the element data vector.
   * The element data of the local mesh elements can be set using \ref set_element_data.
   * If ghost entries should be filled, one should call \ref exchange_ghost_data on each process first.
   * \return Element data vector with data of Type TElementDataType.
   */
  const std::vector<TElementDataType>&
  get_element_data () const
  {
    return m_element_data;
  }

  /** Get the element data vector by moving it out of the competence.
   * In contrast to \ref get_element_data, this transfers ownership of the data instead of
   * returning a reference. After this call the internal element data vector is left empty (moved-from),
   * so \ref set_element_data should be used before accessing the data again.
   * \return Element data vector with data of Type TElementDataType, moved out of the competence.
   */
  std::vector<TElementDataType>
  take_element_data ()
  {
    return std::move (m_element_data);
  }

  /** Exchange the element data for ghost elements between processes.
  * This routine has to be called on each process after setting the element data for all local elements.
  */
  void
  exchange_ghost_data ()
  {
    // Extend element data array to hold also the ghost elements.
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    m_element_data.resize (num_local_elements + num_ghosts);
    // t8_forest_ghost_exchange_data expects an sc_array, so we need to wrap our data array to one.
    sc_array* sc_array_wrapper
      = sc_array_new_data (m_element_data.data (), sizeof (ElementDataType), num_local_elements + num_ghosts);

    // Data exchange: entries with indices > num_local_elements will get overwritten.
    t8_forest_ghost_exchange_data (this->underlying ().get_forest (), sc_array_wrapper);

    sc_array_destroy (sc_array_wrapper);
  }

 protected:
  std::vector<TElementDataType> m_element_data; /**< Vector storing the (local) element data. */
};

/** Wrapper for \ref element_data_mesh_competence_impl to hide TUnderlying and provide the form needed to 
 * pass it as a mesh competence.
 * Use mesh_competence_pack<element_data_mesh_competence<YourElementDataType>::template type> 
 * to get this competence with the correct template parameter form for the mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <T8MPISafeType TElementDataType>
struct element_data_mesh_competence
{
  /** Type to provide the form needed for the mesh competence pack. 
  * \tparam TUnderlying Use the \ref mesh class here.
  */
  template <typename TUnderlying>
  using type = element_data_mesh_competence_impl<TUnderlying, TElementDataType>;
};

// --- Element competence for element data management. ---
/** Element competence to enable that element data can be accessed directly for each element of the mesh.
 * \note This competence requires that the mesh has the \ref element_data_mesh_competence_impl 
 *     (or \ref element_data_mesh_competence) competence that defines the element data vector and the element data type.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct element_data_element_competence: public t8_crtp_operator<TUnderlying, element_data_element_competence>
{
 public:
  /** Set the element data for the element. 
   * \note You can only set element data for non-ghost elements.
   * \param [in] element_data The element data to be set of Type TMeshClass::ElementDataType. 
   *              element_data is moved, not copied.
   */
  void
  set_element_data (auto element_data)
  {
    T8_ASSERT (this->underlying ().m_mesh->has_element_data_handler_competence ());
    SC_CHECK_ABORT (!this->underlying ().is_ghost_element (), "Element data cannot be set for ghost elements.\n");
    // Resize to num_local_elements on first use so that operator[] is valid.S
    if (static_cast<t8_locidx_t> (this->underlying ().m_mesh->m_element_data.size ())
        <= this->underlying ().get_element_handle_id ()) {
      this->underlying ().m_mesh->m_element_data.resize (this->underlying ().m_mesh->get_num_local_elements ());
    }
    // element_data is moved, not copied.
    this->underlying ().m_mesh->m_element_data[this->underlying ().get_element_handle_id ()] = std::move (element_data);
  }

  /** Getter for the element data.
   * For ghost elements ensure that \ref element_data_mesh_competence_impl::exchange_ghost_data is called on 
   * each process first.
   * Element data for non-ghost elements can be accessed (if set) directly.
   * \return Element data with data of Type TMeshClass::ElementDataType.
   */
  const auto&
  get_element_data () const
  {
    T8_ASSERT (this->underlying ().m_mesh->has_element_data_handler_competence ());
    const t8_locidx_t handle_id = this->underlying ().get_element_handle_id ();
    T8_ASSERTF (static_cast<size_t> (handle_id) < this->underlying ().m_mesh->m_element_data.size (),
                "Element data not set.\n");
    return this->underlying ().m_mesh->m_element_data[handle_id];
  }
};

// --- Mesh competence to interpolate data. ---
/** Mesh competence to interpolate the element data after an adaptation step.
 * The \ref element_data_mesh_competence stores a vector of element data, but that data has to be updated if
 * the mesh is adapted, since the elements it refers to are refined, coarsened or reordered. This competence adds the
 * ability to interpolate the data after the adaptation via a user defined callback set using \ref set_interpolate_callback.
 * The next \ref mesh::commit applies it. 
 * \note It therefore only makes sense in combination with the element data competence 
 *       (see \ref interpolate_data_mesh_competence_pack, which bundles the two).
 * \tparam TUnderlying Use the \ref mesh class here.
 */
template <typename TUnderlying>
class interpolate_element_data_mesh_competence:
  public t8_crtp_operator<TUnderlying, interpolate_element_data_mesh_competence> {
 public:
  /** Mesh internal, element-index based storage for the interpolation callback.
   * Users should use the easier span based callback type \ref mesh::interpolate_callback_type.
   * \see set_interpolate_callback for registering a callback.
   * The span based callback is automatically wrapped in \ref set_interpolate_callback to match this type and stored as
   * \ref m_interpolate_callback to be used in the next \ref mesh::commit.
   * \note We can not store or the span based \ref mesh::interpolate_callback_type directly. This competence uses the 
   *       CRTP pattern, so while the competence is instantiated, the mesh (\a TUnderlying) is still an incomplete type.
   *       A data member of type \ref mesh::interpolate_callback_type, would require \c TUnderlying::element_class, 
   *       which is not available for the incomplete type. Therefore we use this index based callback type for storage 
   *       without the need for element_class. Using \ref set_interpolate_callback, we move the element_class lookup 
   *       to the call site.
   * \param [in]     mesh_old  The old mesh that is adapted from.
   * \param [in,out] mesh_new  The new mesh constructed from \a mesh_old.
   * \param [in]     refine    -1 if a family in the old mesh got coarsened, 0 if the element was not touched,
   *                           1 if the element got refined.
   * \param [in]     num_old   The number of outgoing elements.
   * \param [in]     first_old The local mesh handle index of the first outgoing element in the old mesh.
   * \param [in]     num_new   The number of incoming elements.
   * \param [in]     first_new The local mesh handle index of the first incoming element in the new mesh.
  
   */
  using internal_interpolate_callback_type
    = std::function<void (const TUnderlying& mesh_old, TUnderlying& mesh_new, const int refine, const int num_old,
                          const t8_locidx_t first_old, const int num_new, const t8_locidx_t first_new)>;

  /** Register a user callback to interpolate the element data after adaptation.
   * Note that data can only be interpolated for a level difference of at most one.
   * Please use the type \ref mesh::interpolate_callback_type for the callback. 
   * \see mesh::interpolate_callback_type for the expected callback shape and the meaning of its arguments.
   * \note This function is templated on purpose: it is only instantiated at the call site, where the mesh is a
   *       complete type and \ref mesh::element_class is nameable. The element class is needed in the definition 
   *       of the interpolate_callback_type. You do not have to provide this template, it is normally auto deduced.
   * \note This function is templated on purpose: it is only instantiated at the call site, where the mesh is a
   *       complete type and \ref mesh::element_class (needed to name \ref mesh::interpolate_callback_type) is
   *       nameable, which it is not inside this competence. 
   *       The template parameter is deduced from the passed callback, so you do not have to provide it explicitly!
   * \tparam TInterpolateCallback The user callback type \ref mesh::interpolate_callback_type.
   * \param [in] interpolate_callback The span based interpolation callback.
   * 
   */
  template <typename TInterpolateCallback>
  void
  set_interpolate_callback (TInterpolateCallback&& interpolate_callback)
  {
    /* We wrap the user defined, span based callback using \ref detail::to_replace_callback to the index based callback
     * type \ref internal_interpolate_callback_type to be able to store the callback without the need of knowing 
     * \ref mesh::element_class.*/
    m_interpolate_callback
      = detail::to_replace_callback<TUnderlying> (std::forward<TInterpolateCallback> (interpolate_callback));
  }

 protected:
  /** Decide whether \ref mesh::set_partition has been requested for the upcoming \ref mesh::commit.
   * With the interpolation competence the partition step is postponed so that it runs after the element data has been
   * interpolated onto the new mesh; \ref mesh::set_partition therefore records its choice in
   * \ref m_partition_for_coarsening instead of building the partitioned forest directly.
   * \return true if \ref mesh::set_partition has been called (i.e. \ref m_partition_for_coarsening holds a value),
   *         false otherwise.
   */
  bool
  set_partition_called ()
  {
    return m_partition_for_coarsening.has_value ();
  }

  /** Repartition the element data so it follows a newly partitioned forest.
   * The element data currently belongs to \a forest_from; this moves it to the layout of \a forest_to, which must
   * have been created by partitioning \a forest_from. 
   * Analogous to \ref element_data_mesh_competence_impl::exchange_ghost_data, but using \ref t8_forest_partition_data. 
   * This function is called from \ref mesh::commit after the interpolated data has been produced and the partitioned 
   * forest has been committed.
   * \param [in] forest_from The (committed) forest the current element data belongs to.
   * \param [in] forest_to   The committed forest that was partitioned from \a forest_from.
   * \note Both forests could also be accessed directly (by this->underlying()) but this requires that the function is
   * called on the exact right states of m_forest and m_uncommitted_forest. 
   * Providing the variables is the saver implementation.
   */
  void
  repartition_element_data (t8_forest_t forest_from, t8_forest_t forest_to)
  {
    using element_data_type = typename TUnderlying::ElementDataType;
    // Take ownership of old data and wrap into sc_array. This is because the forest functions expect an sc_array.
    std::vector<element_data_type> old_data = this->underlying ().take_element_data ();
    sc_array* data_in = sc_array_new_data (old_data.data (), sizeof (element_data_type),
                                           t8_forest_get_local_num_leaf_elements (forest_from));
    const t8_locidx_t num_new_local = t8_forest_get_local_num_leaf_elements (forest_to);
    // Define vector for the new data and wrap it.
    std::vector<element_data_type> partitioned_data (num_new_local);
    sc_array* data_out = sc_array_new_data (partitioned_data.data (), sizeof (element_data_type), num_new_local);
    // Partition magic using the forest function.
    t8_forest_partition_data (forest_from, forest_to, data_in, data_out);
    // Clean up.
    sc_array_destroy (data_in);
    sc_array_destroy (data_out);
    // Set the partitioned element data to the mesh.
    this->underlying ().set_element_data (std::move (partitioned_data));
  }

  internal_interpolate_callback_type
    m_interpolate_callback; /**< The wrapped element-index based interpolation callback, 
                              * applied on the next \ref mesh::commit. */
  std::optional<bool>
    m_partition_for_coarsening; /**< Postponed \ref mesh::set_partition request: a value means partition on the next 
                          * commit (with value passed to \ref mesh::set_partition); no value means do not partition. */
};

}  // namespace t8_mesh_handle
