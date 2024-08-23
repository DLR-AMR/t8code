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
    const int single_size = this->size (comm);
    int num_data_size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &num_data_size);
    SC_CHECK_MPI (mpiret);
    return num_data_size + num_data * single_size;
  }

  /**
   * Pack a vector of items into a buffer. The first integer of the packed data tells how many
   * items were packed. 
   * 
   * \param[in] data        A vector of items to pack
   * \param[in, out] buffer A vector that will be filled with the packed data. 
   *                        Adds a prefix-int for the size.
   * \param[in] comm        The used communicator
   */
  void
  pack_vector_prefix (const std::vector<T> &data, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int pos = 0;
    T8_ASSERT (buffer.size () == (long unsigned int) buffer_size (data.size (), comm));
    const int num_data = data.size ();
    sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);

    for (const T item : data) {
      this->pack (item, pos, buffer, comm);
    }
  }

  /**
   * Pack a vector of items into a buffer. In contrast to \a pack_vector_sizeprefix no prefix
   * to tell how many items have been packed is used. 
   * 
   * \param[in] data A vector of items to pack
   * \param[in, out] buffer A vector that will be filled with the packed data.
   * \param[in] comm The used communicator
   */
  void
  pack_vector_no_prefix (const std::vector<T> &data, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int pos = 0;
    for (const T item : data) {
      this->pack (item, pos, buffer, comm);
    }
  }

  /**
   * Unpack a buffer into a vector of items. Expects a prefix telling how many items of type T have been packed. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in, out] data Vector of type T, that will be filled with the unpacked data.
   * \param[in, out] outcount Number of items that were packed
   * \param[in] comm The communicator to use. 
   */
  void
  unpack_vector_prefix (std::vector<char> &buffer, std::vector<T> &data, int &outcount, sc_MPI_Comm comm)
  {
    int pos = 0;

    /* Get the number of items we received. */
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    data.resize (outcount);

    for (T item : data) {
      this->unpack (buffer, pos, item, comm);
    }
  }

  /**
   * Unpack a buffer into a vector of items. Does not expect a prefix giving any metadata. 
   * 
   * \param[in] buffer  The input buffer
   * \param[in, out] data Vector of type T, that will be filled with the unpacked data. The vector is already
   *                      allocated to the proper size. 
   * \param[in] comm The communicator to use. 
   */
  void
  unpack_vector_no_prefix (const std::vector<char> &buffer, std::vector<T> &data, sc_MPI_Comm comm)
  {
    int pos = 0;

    for (T item : data) {
      this->unpack (buffer, pos, item, comm);
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
  send (std::vector<T> &data, int dest, int tag, sc_MPI_Comm comm)
  {
#if T8_ENABLE_MPI
    std::vector<char> buffer (buffer_size (data.size (), comm));
    pack_vector_prefix (data, buffer, comm);

    const int mpiret = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, dest, tag, comm);

    return mpiret;
#else
    t8_infof ("send only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
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
  recv (std::vector<T> &data, int source, int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
  {
#if T8_ENABLE_MPI
    int mpiret = sc_MPI_Probe (source, tag, comm, status);
    SC_CHECK_MPI (mpiret);

    int size;
    mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &size);
    SC_CHECK_MPI (mpiret);

    std::vector<char> buffer (size);
    int pos = 0;
    mpiret = sc_MPI_Recv (buffer.data (), buffer.size (), sc_MPI_PACKED, source, pos, comm, status);
    SC_CHECK_MPI (mpiret);

    unpack_vector_prefix (buffer, data, outcount, comm);

    return mpiret;
#else
    t8_infof ("recv only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
  }

  /**
   * Wrapper around an sc_MPI_Allgather, packing the data before sending it and 
   * unpacking it afterwards. 
   * 
   * \param[in] send        The data that should be send to all other ranks in the communicator.
   * \param[in, out] recv   A vector of type \a T to gather the data. Will be set to proper size. 
   * \param[in] comm        The communicator to use. 
   * \return int 
   */
  int
  allgather (const std::vector<T> &send, std::vector<T> &recv, sc_MPI_Comm comm)
  {
#if T8_ENABLE_MPI
    const int bsize = buffer_size (send.size (), comm);
    std::vector<char> buffer (bsize);
    pack_vector_no_prefix (send, buffer, comm);

    int mpisize;
    int mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);

    const int recv_size = mpisize * bsize;
    std::vector<char> recv_buffer (recv_size);
    mpiret = sc_MPI_Allgather (buffer.data (), buffer.size (), sc_MPI_PACKED, recv_buffer.data (), buffer.size (),
                               sc_MPI_PACKED, comm);
    SC_CHECK_MPI (mpiret);

    recv.resize (send.size () * mpisize);

    unpack_vector_no_prefix (recv_buffer, recv, comm);

    return mpiret;
#else
    t8_infof ("Allgather only available when configured with --enable-mpi\n");
    return sc_MPI_ERR_OTHER;
#endif
  }
};

#endif /* T8_DATA_HANDLER_HXX */
