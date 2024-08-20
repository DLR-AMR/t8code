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

#ifndef T8_DATA_HANDLER_HXX
#define T8_DATA_HANDLER_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_data_handler_base.hxx>
#include <t8_data/t8_data_packs/t8_packed_types.hxx>

template <typename T>
class t8_data_handler: public t8_single_data_handler<T> {
 public:
  /**
   * Compute the size of a buffer for \a num_data items of type T
   * 
   * \param[in] num_data Number of items that will be packed into the buffer
   * \param[in] comm The communicator that will be used. 
   * \return The size of the buffer in bytes. 
   */
  int
  buffer_size (const int num_data, sc_MPI_Comm comm)
  {
    const int single_size = this->data_size (comm);
    int num_data_size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &num_data_size);
    SC_CHECK_MPI (mpiret);
    return num_data_size + num_data * single_size;
  }

  /**
   * Pack a vector of items into a buffer. 
   * 
   * \param[in] data A vector of items to pack
   * \param[in, out] buffer A vector that will be filled with the packed data. 
   * \param[in] comm The used communicator
   */
  void
  data_pack_vector (std::vector<T> &data, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int pos = 0;
    T8_ASSERT (buffer.size () == buffer_size (data.size (), comm));
    const int num_data = data.size ();
    sc_MPI_Pack (&num_data, 1, MPI_INT, buffer.data (), buffer.size (), &pos, comm);

    for (int idata = 0; idata < num_data; idata++) {
      this->data_pack (data[idata], pos, buffer, comm);
    }
  }

  /**
   * Unpack a buffer into a vector of items. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in, out] data Vector of type T, that will be filled with the unpacked data.
   * \param[in, out] outcount Number of items that were packed
   * \param[in] comm The communicator to use. 
   */
  void
  data_unpack_vector (std::vector<char> &buffer, std::vector<T> &data, int &outcount, sc_MPI_Comm comm)
  {
    int pos = 0;

    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    data.resize (outcount);

    for (int ipack = 0; ipack < outcount; ++ipack) {
      this->data_unpack (buffer, pos, data[ipack], comm);
    }
  }
};

#endif /* T8_DATA_HANDLER_HXX */
