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

template <typename T>
class t8_data_handler {
 public:
  t8_data_handler ();

  int
  t8_data_size (sc_MPI_Comm comm);

  int
  t8_buffer_size (const int num_data, sc_MPI_Comm comm)
  {
    const int single_size = t8_data_size ();
    const int num_data_size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &num_data_size);
    SC_CHECK_MPI (mpiret);
    return num_data_size + num_data * single_size;
  }

  /**
         * Overwrite this routine to describe how data of type T should be packed
         * 
         * \param[in] data Data to be packed via MPI_Pack
         * \return the size of the packed data in number of bytes. 
         */
  void
  t8_data_pack (const T &data, const int pos, std::vector<char> &buffer, sc_MPI_Comm comm);

  void
  t8_data_pack_vector (const std::vector<T> &data, const int num_data, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int pos = 0;
    int size;
    T8_ASSERT (buffer.size == t8_buffer_size (num_data, comm));

    sc_MPI_Pack_size (1, sc_MPI_Int, comm, &size);
    sc_MPI_Pack (&num_data, 1, MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    pos += size;

    const int single_size = t8_data_size (comm);
    for (int idata = 0; idata < num_data; idata++) {
      t8_data_pack (data[idata], pos, buffer, comm);
      pos += single_size;
    }
  }

  /**
         * Overwrite this routine to describe how data of type T should be unpacked 
         * 
         * \param packed_data A void-pointer to the packed Data
         * \return T* the unpacked data. 
         */
  void
  t8_data_unpack (const void *packed_data, const int insize, const int outcount, std::vector<T> &unpacked_data,
                  sc_MPI_Comm comm);

  void
  t8_data_unpack_vector (const std::vector<char> &packed_data, std::vector<T> &data, sc_MPI_Comm comm)
  {
  }
};

#endif /* T8_DATA_HANDLER_HXX */
