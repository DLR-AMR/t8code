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

/**
 * \file t8_data_handler_specs.hxx
 * This file provides data types for testing data handlers: a templated class that enlarges a data type by a
 * check-int, and a type that has to be serialized by a passed data handler.
 * Should only be used for testing purposes.
 */

#ifndef T8_DATA_HANDLER_SPECS_HXX
#define T8_DATA_HANDLER_SPECS_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_vector_handler.hxx>

/**
 * A template class to handle data with an additional check value.
 *
 * This class is designed to store a piece of data along with an integer
 * value that can be used for additional checks or validations.
 *
 * \tparam TType The type of the data to be stored.
 *
 */
template <typename TType>
struct enlarged_data
{
 public:
  /**
   * Default constructor.
   * Creates an empty enlarged_data object.
   */
  enlarged_data ()
  {
  }

  /**
   * Constructor that initializes the data and check value.
   *
   * \param[in] data The original data to be stored.
   * \param[in] check An integer value for additional checks.
   */
  enlarged_data (TType data, int check): data (data), check (check)
  {
  }

  TType data; /**< original data */
  int check;  /**< additional data to check against */
};

/**
 * A template class to create data of type T.
 *
 * \tparam TType The type of the data to be stored.
 */
template <typename TType>
struct data_creator
{
 public:
  data_creator ()
  {
    large_data = std::vector<TType> (0);
  };

  /**
   * Create several data items.
   *
   * \param[in] num_data  The number of data items to create.
   */
  void
  create (const int num_data);

  std::vector<TType> large_data; /**< A vector to hold the created data items. */
};

/**
 * An item that has to be serialized, since it holds a std::vector. It has no t8_data_handler specialization.
 */
struct t8_test_cell
{
  std::vector<double> values; /**< The values of the cell, e.g. coefficients of a polynomial. */

  /** Compare two cells value by value. */
  bool
  operator== (const t8_test_cell &other) const
    = default;
};

/**
 * A data handler for \ref t8_test_cell that is passed to the vector handler instead of specializing
 * t8_data_handler. Every cell has the same number of values, which is only known at run time, e.g. from a
 * polynomial degree. The handler therefore carries it as state and does not pack it per cell.
 */
struct t8_test_cell_handler
{
  int num_values; /**< The number of values of every cell. */

  /**
   * Returns the number of bytes needed to pack a cell.
   * \param[in] cell The cell.
   * \param[in] comm The MPI communicator.
   * \return The number of bytes.
   */
  int
  size ([[maybe_unused]] const t8_test_cell &cell, sc_MPI_Comm comm) const
  {
    int pack_size = 0;
    const int mpiret = sc_MPI_Pack_size (num_values, sc_MPI_DOUBLE, comm, &pack_size);
    SC_CHECK_MPI (mpiret);
    return pack_size;
  }

  /**
   * Packs a cell into a buffer.
   * \param[in] cell The cell to pack.
   * \param[in, out] buffer The buffer.
   * \param[in] buffer_size The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator.
   */
  void
  pack (const t8_test_cell &cell, void *buffer, const int buffer_size, int &position, sc_MPI_Comm comm) const
  {
    T8_ASSERT (cell.values.size () == static_cast<size_t> (num_values));
    const int mpiret
      = sc_MPI_Pack (cell.values.data (), num_values, sc_MPI_DOUBLE, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Unpacks a cell from a buffer.
   * \param[in] buffer The buffer.
   * \param[in] buffer_size The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the unpacked bytes.
   * \param[out] cell The unpacked cell.
   * \param[in] comm The MPI communicator.
   */
  void
  unpack (const void *buffer, const int buffer_size, int &position, t8_test_cell &cell, sc_MPI_Comm comm) const
  {
    cell.values.resize (num_values);
    const int mpiret
      = sc_MPI_Unpack (buffer, buffer_size, &position, cell.values.data (), num_values, sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);
  }
};

#endif /* T8_DATA_HANDLER_SPECS_HXX */
