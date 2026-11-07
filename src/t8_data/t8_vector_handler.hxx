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
#include <memory>
#include <utility>
#include <limits>

/**
 * Convert a count of bytes or items to the int that MPI expects.
 * Aborts if the count does not fit into an int.
 *
 * \param[in] count The count to convert.
 * \return \a count as int.
 */
inline int
t8_mpi_int (const size_t count)
{
  SC_CHECK_ABORTF (count <= static_cast<size_t> (std::numeric_limits<int>::max ()),
                   "Count %zu exceeds the int range of MPI.", count);
  return static_cast<int> (count);
}

/**
 * A base class for vector handlers.
 */
struct t8_abstract_vector_handler
{
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
   * Sends and receives messages to and from specified destinations.
   *
   * This pure virtual function is responsible for packing, sending, receiving and unpacking a message from and
   * to a given source and destination with a specific tag within the provided MPI communicator. The function will also
   * update the status and output count of the received message.
   *
   * \param[in] dest The destination rank to which the data will be sent.
   * \param[in] source The rank of the source process from which the message is received.
   * \param[in] tag The tag associated with the messages to be sent.
   * \param[in] comm The MPI communicator used for the communication.
   * \param[in] status A pointer to an MPI status object that will be updated with the status of the received message.
   * \param[in] outcount A reference to an integer that will be updated with the count of received elements.
   * \return An integer indicating the status of the send operation.
   */
  virtual int
  sendrecv (const int dest, const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount)
    = 0;

  /**
   * Pure virtual function to get the type.
   *
   * This function must be overridden in derived classes to return the type.
   *
   * \return An integer representing the type.
   */
  virtual t8_data_handler_type
  type ()
    = 0;

  virtual ~t8_abstract_vector_handler () {};
};

/**
 * A vector handler for a std::vector with one item of type TType per entry.
 *
 * Single items are packed and unpacked by a data handler of type TDataHandler. By default this is the
 * \ref t8_data_handler of TType, which handles trivially copyable types without user code. For data that has
 * to be serialized, e.g. items holding a std::vector, pass a handler that describes how one item is packed.
 *
 * \tparam TType The type of the data items.
 * \tparam TDataHandler The type of the handler for single items, see \ref T8DataHandler.
 */
template <typename TType, typename TDataHandler = t8_data_handler<TType>>
  requires T8DataHandler<TDataHandler, TType>
struct t8_vector_handler: public t8_abstract_vector_handler
{
 public:
  /**
   * Construct a t8_vector_handler that owns the given data.
   *
   * \param[in] data The data to be handled. Pass an rvalue (e.g. via std::move) to avoid a copy.
   * \param[in] data_handler The handler for single items, e.g. one that depends on run time parameters.
   */
  explicit t8_vector_handler (std::vector<TType> data = {}, TDataHandler data_handler = {})
    : m_data (std::move (data)), m_data_handler (std::move (data_handler))
  {
  }

  /**
   * Get the data.
   *
   * \return A reference to the data.
   */
  const std::vector<TType> &
  get_data () const
  {
    return m_data;
  }

  /**
   * Replace the data.
   *
   * \param[in] data The new data. Pass an rvalue (e.g. via std::move) to avoid a copy.
   */
  void
  set_data (std::vector<TType> data)
  {
    m_data = std::move (data);
  }

  /**
   * Move the data out of the handler. The handler is empty afterwards.
   *
   * \return The data.
   */
  std::vector<TType>
  take_data ()
  {
    return std::exchange (m_data, {});
  }

  /**
   * Compute the size of the buffer that is needed to pack the data.
   *
   * \param[in] comm The MPI communicator used for the operation.
   * \return The size of the buffer.
   */
  int
  buffer_size (sc_MPI_Comm comm) override
  {
    int prefix_size = 0;
    const int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &prefix_size);
    SC_CHECK_MPI (mpiret);
    return t8_mpi_int (prefix_size + range_buffer_size (0, m_data.size (), comm));
  }

  /**
   * Packs a vector into a buffer. The vector data will be prefixed with the number of elements in the vector.
   *
   * \param[in, out] buffer A pointer to the buffer where the vector prefix will be packed.
   * \param[in] num_bytes The number of bytes to be packed.
   * \param[in, out] pos A reference to an integer representing the current position in the buffer. This will be updated as bytes are packed.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  pack_vector_prefix (void *buffer, const int num_bytes, int &pos, sc_MPI_Comm comm) override
  {
    const int num_data = t8_mpi_int (m_data.size ());
    const int mpiret = sc_MPI_Pack (&num_data, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
    pack_range (0, m_data.size (), buffer, num_bytes, pos, comm);
  }

  /**
   * Unpacks a vector from a buffer. Expected to be prefixed with the number of elements in the vector.
   *
   * \param[in] buffer Pointer to the buffer containing the packed data.
   * \param[in] num_bytes The number of bytes in the buffer.
   * \param[in] pos Reference to an integer representing the current position in the buffer. This will be updated as data is unpacked.
   * \param[in] outcount Reference to an integer where the count of unpacked elements will be stored.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &pos, int &outcount, sc_MPI_Comm comm) override
  {
    const int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    m_data.resize (outcount);
    unpack_range (buffer, num_bytes, pos, m_data.data (), m_data.size (), comm);
  }

  /**
   * Sends a message to a specified destination.
   *
   * This function is responsible for packing and sending data to a given destination
   * with a specific tag using the provided MPI communicator.
   *
   * \param[in] dest The destination rank to which the data will be sent.
   * \param[in] tag The tag associated with the message to be sent.
   * \param[in] comm The MPI communicator used for the communication.
   * \return An integer indicating the status of the send operation.
   */
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

  /**
   * Receives a message from a specified source.
   *
   * This function is responsible for receiving and unpacking a message from a given source
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

  /**
   * Sends and receives messages to and from specified destinations.
   *
   * This function is responsible for packing, sending and receiving data to and from given destinations
   * with a specific tag using the provided MPI communicator.
   *
   * \param[in] dest The destination rank to which the data will be sent.
   * \param[in] source The rank of the source process from which the message is received.
   * \param[in] tag The tag associated with the messages to be sent.
   * \param[in] comm The MPI communicator used for the communication.
   * \param[in] status A pointer to an MPI status object that will be updated with the status of the received message.
   * \param[in] outcount A reference to an integer that will be updated with the count of received elements.
   * \return An integer indicating the status of the send operation.
   */
  int
  sendrecv ([[maybe_unused]] const int dest, [[maybe_unused]] const int source, [[maybe_unused]] const int tag,
            [[maybe_unused]] sc_MPI_Comm comm, [[maybe_unused]] sc_MPI_Status *status,
            [[maybe_unused]] int &outcount) override
  {
    {
#if T8_ENABLE_MPI
      /* Compute sizes of outgoing and incoming messages */
      const int send_num_bytes = buffer_size (comm);
      int recv_num_bytes = 0;
      int mpiret = sc_MPI_Sendrecv (&send_num_bytes, 1, sc_MPI_INT, dest, tag, &recv_num_bytes, 1, sc_MPI_INT, source,
                                    tag, comm, status);
      SC_CHECK_MPI (mpiret);

      /* Prepare send buffer */
      std::vector<char> send_buffer (send_num_bytes);
      int send_pos = 0;
      pack_vector_prefix (send_buffer.data (), send_num_bytes, send_pos, comm);

      /* Prepare receive buffer using size received */
      std::vector<char> recv_buffer (recv_num_bytes);

      /* Actual sendrecv of data */
      mpiret = sc_MPI_Sendrecv (
        /* send payload */
        send_buffer.data (), send_num_bytes, sc_MPI_PACKED, dest, tag,
        /* receive payload */
        recv_buffer.data (), recv_num_bytes, sc_MPI_PACKED, source, tag,
        /* communicator + status */
        comm, status);
      SC_CHECK_MPI (mpiret);

      /* Unpack received data */
      int recv_pos = 0;
      unpack_vector_prefix (recv_buffer.data (), recv_num_bytes, recv_pos, outcount, comm);

      return mpiret;
#else
      SC_ABORT ("MPI_Sendrecv not available without MPI");
      return 0;
#endif
    }
  }

  /**
   * Get the type of the data handler.
   *
   * \return The type of the data handler, \ref t8_data_handler_type_unregistered if the data handler has no
   *         member type ().
   */
  t8_data_handler_type
  type () override
  {
    if constexpr (requires { m_data_handler.type (); }) {
      return m_data_handler.type ();
    }
    else {
      return t8_data_handler_type_unregistered;
    }
  }

 private:
  /**
   * Compute the number of bytes needed to pack a range of items (without count prefix).
   *
   * \param[in] first The index of the first item of the range.
   * \param[in] count The number of items in the range.
   * \param[in] comm The MPI communicator used for the operation.
   * \return The number of bytes needed to pack the items [first, first + count).
   */
  size_t
  range_buffer_size (const size_t first, const size_t count, sc_MPI_Comm comm) const
  {
    T8_ASSERT (first + count <= m_data.size ());
    if constexpr (T8BytewiseDataHandler<TDataHandler, TType>) {
      /* Bytewise items are packed in bulk. */
      int pack_size = 0;
      const int mpiret = sc_MPI_Pack_size (t8_mpi_int (count * sizeof (TType)), sc_MPI_BYTE, comm, &pack_size);
      SC_CHECK_MPI (mpiret);
      return pack_size;
    }
    else {
      size_t pack_size = 0;
      for (size_t iitem = first; iitem < first + count; ++iitem) {
        pack_size += m_data_handler.size (m_data[iitem], comm);
      }
      return pack_size;
    }
  }

  /**
   * Pack a range of items into a buffer (without count prefix).
   *
   * \param[in] first The index of the first item of the range.
   * \param[in] count The number of items in the range.
   * \param[in, out] buffer The buffer to pack the items into.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  pack_range (const size_t first, const size_t count, void *buffer, const int num_bytes, int &position,
              sc_MPI_Comm comm) const
  {
    T8_ASSERT (first + count <= m_data.size ());
    if constexpr (T8BytewiseDataHandler<TDataHandler, TType>) {
      /* Bytewise items are packed in bulk. */
      const int mpiret = sc_MPI_Pack (m_data.data () + first, t8_mpi_int (count * sizeof (TType)), sc_MPI_BYTE, buffer,
                                      num_bytes, &position, comm);
      SC_CHECK_MPI (mpiret);
    }
    else {
      for (size_t iitem = first; iitem < first + count; ++iitem) {
        m_data_handler.pack (m_data[iitem], buffer, num_bytes, position, comm);
      }
    }
  }

  /**
   * Unpack a range of items from a buffer (without count prefix).
   *
   * \param[in] buffer The buffer containing the packed items.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the unpacked bytes.
   * \param[out] items Array of at least \a count items to unpack into.
   * \param[in] count The number of items to unpack.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  unpack_range (const void *buffer, const int num_bytes, int &position, TType *items, const size_t count,
                sc_MPI_Comm comm) const
  {
    if constexpr (T8BytewiseDataHandler<TDataHandler, TType>) {
      /* Bytewise items are unpacked in bulk. */
      const int mpiret
        = sc_MPI_Unpack (buffer, num_bytes, &position, items, t8_mpi_int (count * sizeof (TType)), sc_MPI_BYTE, comm);
      SC_CHECK_MPI (mpiret);
    }
    else {
      for (size_t iitem = 0; iitem < count; ++iitem) {
        m_data_handler.unpack (buffer, num_bytes, position, items[iitem], comm);
      }
    }
  }

  std::vector<TType> m_data;   /**< The data that is packed, unpacked and communicated via MPI. */
  TDataHandler m_data_handler; /**< Handler to pack and unpack single items. */
};

/**
 * Create an internal handler for the given data type. Can be used for complex data types that use
 * data structures implemented by t8code.
 *
 * \param[in] type The type of the data handler to create.
 * \return A handler for the given data type.
 */
inline t8_abstract_vector_handler *
create_internal_handler (const t8_data_handler_type type)
{
  switch (type.get ()) {
    /* Place holder to create a handler for internal data structures */
  default:
    return nullptr;
  }
}

#endif /* T8_VECTOR_HANDLER_HXX */
