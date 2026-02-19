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
 */
#pragma once

#include <t8_types/t8_crtp.hxx>
#include <type_traits>
#include <vector>

namespace t8_mesh_handle
{

/** Concept to ensure that a type is MPI safe.
 */
template <typename TType>
concept T8MPISafeType
  = std::is_void_v<TType> || (std::is_trivially_copyable_v<TType> && std::is_standard_layout_v<TType>);

/**
 * Handler for the element data of a \ref mesh.
 * \tparam TElementDataType The element data type you want to use for each element of the mesh. 
 *         The data type has to be MPI safe as the data for ghost elements will be exchanged via MPI.
 */
template <typename TUnderlying, T8MPISafeType TElementDataType>
class handle_element_data: public t8_crtp_basic<TUnderlying> {
 public:
  using ElementDataType = TElementDataType; /**< Make Type of the element data accessible. */
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
    const auto num_local_elements = this->underlying ().get_num_local_elements ();
    const auto num_ghosts = this->underlying ().get_num_ghosts ();
    T8_ASSERT (element_data.size () == static_cast<size_t> (num_local_elements));
    m_element_data = std::move (element_data);
    m_element_data.reserve (num_local_elements + num_ghosts);
    m_element_data.resize (num_local_elements);
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

/**
 * Helper alias to create a mesh competence for element data
 * without exposing TUnderlying to the user.
 *
 * Usage:
 *   using mesh_type =
 *     mesh<competence_pack<...>, element_data_competence<MyData>>;
 */
// template <T8MPISafeType TElementDataType>
// struct element_data_wrapper
// {
//   template <typename TUnderlying>
//   using type = handle_element_data<TUnderlying, TElementDataType>;
// };

// template <T8MPISafeType TElementDataType>
// using element_data_mesh_competence = typename element_data_wrapper<TElementDataType>::template type;

}  // namespace t8_mesh_handle
