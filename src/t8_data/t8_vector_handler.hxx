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
 * \file t8_vector_handler.hxx
 * 
 * This file provides functions to pack and unpack data for communication and to send and receive it using MPI.
 * 
 */

#ifndef T8_VECTOR_HANDLER_HXX
#define T8_VECTOR_HANDLER_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_data_handler.hxx>
#include <algorithm>
#include <memory>
#include <numeric>

class t8_abstract_vector_handler {
 public:
  /**
   * Pure virtual function to determine the buffer size.
   *
   * This function must be implemented by derived classes to calculate
   * the size of the buffer required for communication.
   *
   * \param[in] comm The MPI communicator.
   * \return The size of the buffer.
   */
  virtual int
  buffer_size (sc_MPI_Comm comm)
    = 0;

  /**
   * Packs a vector into a buffer. The vector data will be prefixed with the number of elements in the vector.
   *
   * This pure virtual function is responsible for packing a vector prefix into the provided buffer.
   *
   * \param[in, out] buffer A pointer to the buffer where the vector prefix will be packed.
   * \param[in] num_bytes The number of bytes to be packed.
   * \param[in, out] pos A reference to an integer representing the current position in the buffer. This will be updated as bytes are packed.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm)
    = 0;

  /**
   * Unpacks a vector from a buffer. Expected to be prefixed with the number of elements in the vector.
   *
   * This pure virtual function is responsible for unpacking a vector prefix from the provided buffer.
   *
   * \param[in] buffer Pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos Reference to an integer representing the current position in the buffer. This will be updated as data is unpacked.
   * \param[in] outcount Reference to an integer where the count of unpacked elements will be stored.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm)
    = 0;

  /**
   * Pure virtual function to send data to a specified destination.
   *
   * This function is responsible for packing and sending data to a given destination
   * with a specific tag using the provided MPI communicator.
   *
   * \param[in] dest The destination rank to which the data will be sent.
   * \param[in] tag The tag associated with the message to be sent.
   * \param[in] comm The MPI communicator used for the communication.
   * \return An integer indicating the status of the send operation.
   */
  virtual int
  send (const int dest, const int tag, sc_MPI_Comm comm)
    = 0;

  /**
   * Receives a message from a specified source.
   *
   * This pure virtual function is responsible for receiving and unpacking a message from a given source
   * with a specific tag within the provided MPI communicator. The function will also
   * update the status and output count of the received message.
   *
   * \param[in] source The rank of the source process from which the message is received.
   * \param[in] tag The tag of the message to be received.
   * \param[in] comm The MPI communicator within which the message is received.
   * \param[in] status A pointer to an MPI status object that will be updated with the status of the received message.
   * \param[in] outcount A reference to an integer that will be updated with the count of received elements.
   * \return An integer indicating the success or failure of the receive operation.
   */
  virtual int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
    = 0;

  /**
   * Pure virtual function to get the type.
   * 
   * This function must be overridden in derived classes to return the type.
   * 
   * \return An integer representing the type.
   * TODO: implement a proper type/enum for this. 
   */
  virtual int
  type ()
    = 0;

  virtual ~t8_abstract_vector_handler () {};
};

/**
 * A template class for handling data in a distributed environment.
 *
 * This class inherits from t8_abstract_data_handler and provides methods for
 * packing, unpacking, sending, and receiving data using MPI.
 *
 * \tparam TType The type of data to be handled.
 */
template <typename TType>
class t8_vector_handler: public t8_abstract_vector_handler {
 public:
  /**
   * Construct a new t8 data handler.
   * m_data is initialized to nullptr.
   */
  t8_vector_handler (): single_handler ()
  {
    m_data = nullptr;
  }

  /**
   * Construct a new t8 data handler with the given data.
   *
   * \param[in] data The data to be handled.
   */
  t8_vector_handler (const std::vector<TType> &data)
    : m_data (std::make_shared<std::vector<TType>> (data)), single_handler ()
  {
  }

  /**
   * Get the data.
   * 
   * \return std::shared_ptr<std::vector<TType>> 
   */
  std::shared_ptr<std::vector<TType>>
  get_data () const
  {
    return m_data;
  }

  int
  buffer_size (sc_MPI_Comm comm) override
  {
    int total_size = 0;
    const int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &total_size);
    SC_CHECK_MPI (mpiret);
    if (m_data) {
      total_size += std::accumulate (m_data->begin (), m_data->end (), 0, [&] (int sum, const TType &item) {
        return sum + single_handler.size (item, comm);
      });
    }
    return total_size;
  }

  void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm) override
  {
    const int num_data = m_data->size ();
    const int mpiret = sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    std::for_each (m_data->begin (), m_data->end (),
                   [&] (const TType &item) { single_handler.pack (item, pos, buffer, num_bytes, comm); });
  }

  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    const int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    if (!m_data) {
      m_data = std::make_shared<std::vector<TType>> (outcount);
    }
    else {
      m_data->clear ();
      m_data->resize (outcount);
    }
    std::for_each (m_data->begin (), m_data->end (),
                   [&] (TType &item) { single_handler.unpack (buffer, num_bytes, pos, item, comm); });
  }

  int
  send (const int dest, const int tag, sc_MPI_Comm comm) override
  {
    int pos = 0;
    const int num_bytes = buffer_size (comm);
    std::vector<char> buffer (num_bytes);
    pack_vector_prefix (buffer.data (), num_bytes, pos, comm);

    const int mpiret = sc_MPI_Send (buffer.data (), num_bytes, sc_MPI_PACKED, dest, tag, comm);
    SC_CHECK_MPI (mpiret);
    return mpiret;
  }

  int
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount) override
  {
    int pos = 0;
    int mpiret = sc_MPI_Probe (source, tag, comm, status);
    SC_CHECK_MPI (mpiret);

    int num_bytes;
    mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &num_bytes);
    SC_CHECK_MPI (mpiret);
    std::vector<char> buffer (num_bytes);

    mpiret = sc_MPI_Recv (buffer.data (), num_bytes, sc_MPI_PACKED, source, tag, comm, status);
    SC_CHECK_MPI (mpiret);
    unpack_vector_prefix (buffer.data (), num_bytes, pos, outcount, comm);
    return mpiret;
  }

  int
  type () override
  {
    return single_handler.type ();
  }

 private:
  /**
  * A shared pointer to a vector of data. 
  * This data will be packed, unpacked, and communicated via MPI.
  */
  std::shared_ptr<std::vector<TType>> m_data;
  /**
   * A single data handler for the data type T.
   * This handler will be used to pack and unpack individual data items.
   */
  t8_data_handler<TType> single_handler;
};

inline t8_abstract_vector_handler *
create_internal_handler (const int type)
{
  switch (type) {
    /* Place holder to create a handler for internal data structures */
  default:
    return nullptr;
  }
}

#endif /* T8_VECTOR_HANDLER_HXX */
