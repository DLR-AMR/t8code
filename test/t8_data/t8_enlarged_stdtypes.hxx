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

/**
 * \file t8_enlarged_stdtypes.hxx 
 * 
 * This file provides a templated class for testing single data items.
 * 
 */

#ifndef T8_ENLARGED_STDTYPES_HXX
#define T8_ENLARGED_STDTYPES_HXX

#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler.hxx>

/**
 * pseudo_types for testing
 * 
 */
enum pseudo_types { T8_ENLARGED_INT = 0, T8_ENLARGED_DOUBLE = 1 };

/**
 * A template specialisation for handling single enlarged ints (int plus an additional int, this 
 * data-type is meant for testing).
 * 
 * This class implements methods for packing, unpacking, and determining the size of single data items.
 * 
 */
template <>
class t8_data_handler<enlarged_data<int>> {
 public:
  /**
   * Returns the size of an enlarged_int.
   * 
   * \param[in] data The data to compute the size of. 
   * \param[in] comm The MPI communicator used for communication.
   * \return An integer representing the size of the data.
   */
  inline int
  size ([[maybe_unused]]const enlarged_data<int> &item, sc_MPI_Comm comm)
  {
    int size;
    const int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &size);
    SC_CHECK_MPI (mpiret);
    return size;
  }

  /**
   * Packs an enlarged_int into a buffer for communication.
   * 
   * \param[in] data The data to be packed.
   * \param[in] pos The current position in the buffer where the data should be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] num_bytes The number of bytes available in the buffer.
   * \param[in] comm The MPI communicator used for communication.
   */
  inline void
  pack (const enlarged_data<int> &data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&(data.data), 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&(data.check), 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Unpacks an enlarged_int from a buffer.
   * 
   * \param[in] buffer A pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos A reference to an integer representing the current position in the buffer.
   * \param[in, out] data A pointer to the data structure where the unpacked data will be stored.
   * \param[in] comm The MPI communicator used for communication.
   */
  inline void
  unpack (const void *buffer, const int num_bytes, int &pos, enlarged_data<int> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &(data.data), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &(data.check), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Returns the type of the data handler.
   * 
   * This function returns the type of the data handler.
   * 
   * \return An integer representing the type.
   */
  constexpr t8_data_handler_type
  type ()
  {
    return t8_data_handler_type (T8_ENLARGED_INT);
  }
};

/**
 * t8_data_handler
 * A template specialization for handling single enlarged doubles (double plus an additional int, this 
 * data-type is meant for testing).
 * 
 * This class implements methods for packing, unpacking, and determining the size of single data items.
 * 
 */
template <>
class t8_data_handler<enlarged_data<double>> {
 public:
  /**
   * Returns the size of an enlarged_double.
   * 
   * \param[in] data The data to compute the size of. 
   * \param[in] comm The MPI communicator used for communication.
   * \return An integer representing the size of the data.
   */
  inline int
  size ([[maybe_unused]] const enlarged_data<double> &item, sc_MPI_Comm comm)
  {
    int int_size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &int_size);
    SC_CHECK_MPI (mpiret);
    int double_size;
    mpiret = sc_MPI_Pack_size (1, sc_MPI_DOUBLE, comm, &double_size);
    SC_CHECK_MPI (mpiret);
    return int_size + double_size;
  }

  /**
   * Packs an enlarged_double into a buffer for communication.
   * 
   * \param[in] data The data to be packed.
   * \param[in] pos The current position in the buffer where the data should be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] num_bytes The number of bytes available in the buffer.
   * \param[in] comm The MPI communicator used for communication.
   */
  inline void
  pack (const enlarged_data<double> &data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&(data.data), 1, sc_MPI_DOUBLE, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&(data.check), 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Unpacks an enlarged_double from a buffer.
   * 
   * \param[in] buffer A pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos A reference to an integer representing the current position in the buffer.
   * \param[in, out] data A pointer to the data structure where the unpacked data will be stored.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  unpack (const void *buffer, const int num_bytes, int &pos, enlarged_data<double> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &(data.data), 1, sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &(data.check), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Returns the type of the data handler.
   * 
   * This function returns the type of the data handler.
   * 
   * \return An integer representing the type.
   */
  constexpr t8_data_handler_type
  type ()
  {
    return t8_data_handler_type (T8_ENLARGED_DOUBLE);
  }
};

#endif /* T8_ENLARGED_STDTYPES_HXX */
