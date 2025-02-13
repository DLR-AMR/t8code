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

#ifndef T8_DATA_HANDLER_HXX
#define T8_DATA_HANDLER_HXX

/**
 * \file t8_data_handler.hxx
 * 
 * This file provides a templated class for handling single data items.
 * 
 */
#include <t8.h>
#include <t8_types/t8_data_handler_type.hxx>

constexpr bool
is_internal_data (const t8_data_handler_type type)
{
  switch (type.underlying ().get ()) {
  // Placeholder for future internal data, which is handled here.
  default:
    return false;
  }
};

/**
 * Template class for handling single data items.
 *
 * This class implements methods for packing, unpacking, and determining the size of single data items.
 *
 * \tparam TType The type of data to be handled.
 */
template <typename TType>
class t8_data_handler {
 public:
  /**
   * Construct a new t8 single data handler.
   * 
   */
  t8_data_handler () {};

  /** 
   * Returns the size of the data.
   * 
   * \param[in] data The data to compute the size of. 
   * \param[in] comm The MPI communicator used for communication.
   * \return An integer representing the size of the data.
   */
  int
  size (const TType &data, sc_MPI_Comm comm);

  /**
   * Packs the given data into a buffer for communication.
   *
   * \param[in] data The data to be packed.
   * \param[in] pos The current position in the buffer where the data should be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] num_bytes The number of bytes available in the buffer.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  pack (const TType &data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm);

  /**
   * Unpacks data from a buffer.
   *
   * This function unpacks data from a given buffer into the provided data structure.
   *
   * \param[in] buffer A pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos A reference to an integer representing the current position in the buffer.
   * \param[in, out] data A pointer to the data structure where the unpacked data will be stored.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  unpack (const void *buffer, const int num_bytes, int &pos, TType &data, sc_MPI_Comm comm);

  /**
   * Returns the type of the data handler.
   *
   * This function returns the type of the data handler.
   *
   * \return An integer representing the type.
   */
  t8_data_handler_type
  type ();

  ~t8_data_handler () {};
};

#endif /* T8_DATA_HANDLER_HXX */
