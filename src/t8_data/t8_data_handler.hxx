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

    for (T item : data) {
      this->data_pack (item, pos, buffer, comm);
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

    /* Get the number of items we received. */
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    data.resize (outcount);

    for (T item : data) {
      this->data_unpack (buffer, pos, item, comm);
    }
  }

  /**
   * Wrapper around a \a data_pack_vector and an sc_MPI_Send. 
   * Packs the \a data and sends it to rank \a dest using \a tag via \a comm
   * 
   * \param[in] data The data to pack and send
   * \param[in] dest The rank we send to. 
   * \param[in] tag The tag to use during communication
   * \param[in] comm The communicator to use. 
   * \return The result of the mpi-communication
   */
  int
  data_send (std::vector<T> &data, int dest, int tag, sc_MPI_Comm comm)
  {
    std::vector<char> buffer (buffer_size (data.size (), comm));
    data_pack_vector (data, buffer, comm);

    const int mpiret = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, dest, tag, comm);

    return mpiret;
  }

  /**
   * Wrapper around an \a sc_MPI_Recv and \a data_unpack. 
   * Receives and unpackes data coming from \a source. 
   * 
   * \param[in, out] data The output buffer. Will be filled with the unpacked data.
   * \param[in] source The rank we receive data from
   * \param[in] tag The tag used during communication
   * \param[in] comm The communicator to use. 
   * \param[in] status Status of the MPI-communication
   * \param[in, out] outcount After execution it is the number of items of type \a T received. 
   * \return The result of the mpi communication.  
   */
  int
  data_recv (std::vector<T> &data, int source, int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
  {
    int mpiret = sc_MPI_Probe (source, tag, comm, status);
    SC_CHECK_MPI (mpiret);

    int size;
    mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &size);
    SC_CHECK_MPI (mpiret);

    std::vector<char> buffer (size);
    int pos = 0;
    mpiret = sc_MPI_Recv (buffer.data (), buffer.size (), sc_MPI_PACKED, source, pos, comm, status);
    SC_CHECK_MPI (mpiret);

    data_unpack_vector (buffer, data, outcount, comm);

    return mpiret;
  }
};

#endif /* T8_DATA_HANDLER_HXX */
