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
#include <test/t8_data/t8_data_handler_specs.hxx>

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
  t8_data_pack (T &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm);

  void
  t8_data_pack_vector (const std::vector<T> &data, const int num_data, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int pos = 0;
    T8_ASSERT (buffer.size == t8_buffer_size (num_data, comm));

    sc_MPI_Pack (&num_data, 1, MPI_INT, buffer.data (), buffer.size (), &pos, comm);

    for (int idata = 0; idata < num_data; idata++) {
      t8_data_pack (data[idata], pos, buffer, comm);
    }
  }

  /**
         * Overwrite this routine to describe how data of type T should be unpacked 
         * 
         * \param packed_data A void-pointer to the packed Data
         * \return T* the unpacked data. 
         */
  void
  t8_data_unpack (std::vector<char> &buffer, int &pos, T &data, sc_MPI_Comm comm);

  void
  t8_data_unpack_vector (std::vector<char> &buffer, std::vector<T> &data, int &outcount, sc_MPI_Comm comm)
  {
    int pos = 0;

    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    data.resize (outcount);

    for (int ipack = 0; ipack < outcount; ++ipack) {
      t8_data_unpack (buffer, pos, data[ipack], comm);
    }
  }
};

template <>
class t8_data_handler<enlarged_data<int>>;

#endif /* T8_DATA_HANDLER_HXX */
