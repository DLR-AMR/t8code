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

/** \file data_handler.hxx
 * Handler for the element data of a \ref t8_mesh_handle::mesh.
 * The file defines a mesh and an element competence for element data handling.
 * Use both competences together if you want to manage element data for the elements of the mesh and access it directly for each element.
 */
#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_crtp.hxx>
#include <type_traits>
#include <vector>

namespace t8_mesh_handle
{
// --- Mesh competence for element data management. ---
/** Concept to ensure that a type is MPI safe. */
template <typename TType>
concept T8MPISafeType
  = std::is_void_v<TType> || (std::is_trivially_copyable_v<TType> && std::is_standard_layout_v<TType>);

/** Handler for the element data of a \ref mesh.
 * Use this competence if you want to manage element data for the elements of the mesh.
 * Use the helper \ref element_data_competence to get this competence with the correct template parameters form for the mesh. 
 * If you want to access the data not only in vector form but also directly for each element, 
 * you can combine this competence with the \ref access_element_data competence.
 * In summary you can use the competences like this: 
 *    mesh<element_competence_pack<access_element_data>,
 *         mesh_competence_pack<element_data_competence<YourElementDataType>::template type>>;
 *
 * \tparam TUnderlying Use the \ref mesh class here.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class handle_element_data: public t8_crtp_basic<TUnderlying> {
 public:
  using ElementDataType = TElementDataType; /**< Make Type of the element data publicly accessible. */

  /** Set the element data vector. The vector should have the length of num_local_elements.
   * \param [in] element_data The element data vector to set with one entry of class TElementDataType 
   *            for each local mesh element (excluding ghosts).
   */
  void
  set_element_data (std::vector<TElementDataType> element_data)
  {
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    T8_ASSERT (element_data.size () == static_cast<size_t> (num_local_elements));
    m_element_data = std::move (element_data);
    m_element_data.reserve (num_local_elements + num_ghosts);
    m_element_data.resize (num_local_elements);
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

  /** Exchange the element data for ghost elements between processes.
  * This routine has to be called on each process after setting the element data for all local elements.
  */
  void
  exchange_ghost_data ()
  {
    // t8_forest_ghost_exchange_data expects an sc_array, so we need to wrap our data array to one.
    sc_array* sc_array_wrapper;
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    m_element_data.resize (num_local_elements + num_ghosts);
    sc_array_wrapper
      = sc_array_new_data (m_element_data.data (), sizeof (ElementDataType), num_local_elements + num_ghosts);

    // Data exchange: entries with indices > num_local_elements will get overwritten.
    t8_forest_ghost_exchange_data (this->underlying ().get_forest (), sc_array_wrapper);

    sc_array_destroy (sc_array_wrapper);
  }

 protected:
  std::vector<TElementDataType> m_element_data; /**< Vector storing the (local) element data. */
};

/** Wrapper for \ref handle_element_data to hide TUnderlying and provide the form needed to pass it as a mesh competence.
 * Use mesh_competence_pack<element_data_competence<YourElementDataType>::template type> 
 * to get this competence with the correct template parameter form for the mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <T8MPISafeType TElementDataType>
struct element_data_competence
{
  /** Type to provide the form needed for the mesh competence pack. 
  * \tparam TUnderlying Use the \ref mesh class here.
  */
  template <typename TUnderlying>
  using type = handle_element_data<TUnderlying, TElementDataType>;
};

// --- Element competence for element data management. ---
/** Element competence to enable that element data can be accessed directly for each element of the mesh.
 * \note This competence requires that the mesh has the \ref handle_element_data competence to manage the
 *   element data vector and exchange ghost data.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct access_element_data: public t8_crtp_basic<TUnderlying>
{
 public:
  // --- Getter and setter for element data. ---
  /** Getter for the element data.
   * For ghost elements ensure that \ref handle_element_data::exchange_ghost_data is called on each process first.
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

  /** Set the element data for the element. 
   * \note You can only set element data for non-ghost elements.
   * \param [in] element_data The element data to be set of Type TMeshClass::ElementDataType.
   */
  void
  set_element_data (auto element_data)
  {
    T8_ASSERT (this->underlying ().m_mesh->has_element_data_handler_competence ());
    SC_CHECK_ABORT (!this->underlying ().is_ghost_element (), "Element data cannot be set for ghost elements.\n");
    // Resize for the case that no data vector has been set previously.
    this->underlying ().m_mesh->m_element_data.reserve (this->underlying ().m_mesh->get_num_local_elements ()
                                                        + this->underlying ().m_mesh->get_num_ghosts ());
    this->underlying ().m_mesh->m_element_data.resize (this->underlying ().m_mesh->get_num_local_elements ());
    this->underlying ().m_mesh->m_element_data[this->underlying ().get_element_handle_id ()] = std::move (element_data);
  }
};

}  // namespace t8_mesh_handle
