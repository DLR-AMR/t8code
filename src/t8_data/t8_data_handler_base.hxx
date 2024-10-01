/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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

#ifndef T8_DATA_HANDLER_BASE
#define T8_DATA_HANDLER_BASE

#include <t8.h>

template <typename T>
class t8_single_data_handler {
 public:
  t8_single_data_handler () {};

  int
  size (const T &data, sc_MPI_Comm comm);

  /**
   * Packs the given data into a buffer for communication.
   *
   * \tparam T The type of the data to be packed.
   * \param[in] data The data to be packed.
   * \param[in] pos The current position in the buffer where the data should be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] num_bytes The number of bytes available in the buffer.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  pack (const T &data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm);

  /**
   * Unpacks data from a buffer.
   *
   * This function unpacks data from a given buffer into the provided data structure.
   *
   * \tparam T The type of the data to be unpacked.
   * \param[in] buffer A pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos A reference to an integer representing the current position in the buffer.
   * \param[in, out] data A pointer to the data structure where the unpacked data will be stored.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  unpack (const void *buffer, const int num_bytes, int &pos, T &data, sc_MPI_Comm comm);

  /**
   * Returns the type of the data handler.
   *
   * This function returns the type of the data handler.
   *
   * \return An integer representing the type.
   */
  int
  type ();

  ~t8_single_data_handler () {};
};

/**
 * \class t8_single_data_handler_c
 * @\rief A class for handling single data operations in an MPI environment.
 *
 * This class provides methods to determine the size of data, pack data into a buffer,
 * unpack data from a buffer, and get the type of data.
 */
struct t8_single_data_handler_c
{

  int
  size (const void *data, sc_MPI_Comm comm);

  void
  pack (const void *data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm);

  void
  unpack (const void *buffer, const int num_bytes, int &pos, void *data, sc_MPI_Comm comm);

  int
  type ();
};

#endif /* T8_DATA_HANDLER_BASE */
