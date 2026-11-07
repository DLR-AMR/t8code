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
 * Vector handlers own a vector of data items (e.g. one item per element) and pack, unpack and
 * communicate it via MPI. Single items are handled by a data handler, see \ref t8_data_handler.
 */

#ifndef T8_VECTOR_HANDLER_HXX
#define T8_VECTOR_HANDLER_HXX

#include <t8.h>
#include <t8_data/t8_data_exchange.hxx>
#include <t8_data/t8_data_handler.hxx>
#include <memory>
#include <utility>
#include <vector>

/**
 * A base class for vector handlers.
 *
 * A vector handler owns a vector of data items and can pack, unpack and communicate it via MPI.
 * Derived classes implement the item-type dependent packing; sending and receiving is implemented
 * once in this base class on top of it.
 */
struct t8_abstract_vector_handler
{
  /** Default destructor. */
  virtual ~t8_abstract_vector_handler () = default;

  /**
   * Compute an upper bound of the number of bytes needed by \ref pack_vector_prefix.
   *
   * \param[in] comm The MPI communicator.
   * \return The size of the buffer.
   */
  virtual int
  buffer_size (sc_MPI_Comm comm) const
    = 0;

  /**
   * Packs the vector into a buffer. The data is prefixed with the number of items in the vector.
   *
   * \param[in, out] buffer The buffer to pack the vector into.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  pack_vector_prefix (void *buffer, const int num_bytes, int &position, sc_MPI_Comm comm) const
    = 0;

  /**
   * Unpacks a vector from a buffer and replaces the data of this handler.
   * The data is expected to be prefixed with the number of items in the vector.
   *
   * \param[in] buffer The buffer containing the packed data.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the unpacked bytes.
   * \param[out] outcount The number of unpacked items.
   * \param[in] comm The MPI communicator used for the operation.
   */
  virtual void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &position, int &outcount, sc_MPI_Comm comm)
    = 0;

  /**
   * Get the type of the handled data.
   *
   * \return The type of the data handler.
   */
  virtual t8_data_handler_type
  type () const
    = 0;

  /**
   * Pack the vector and send it to a process (blocking). Requires MPI.
   *
   * \param[in] dest The rank to send to.
   * \param[in] tag The tag of the message.
   * \param[in] comm The MPI communicator used for the communication.
   */
  void
  send (const int dest, const int tag, sc_MPI_Comm comm) const;

  /**
   * Receive a vector from a process (blocking) and replace the data of this handler with it. Requires MPI.
   *
   * \param[in] source The rank to receive from.
   * \param[in] tag The tag of the message.
   * \param[in] comm The MPI communicator used for the communication.
   * \param[out] status The MPI status of the received message.
   * \param[out] outcount The number of received items.
   */
  void
  recv (const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount);

  /**
   * Send the vector to one process and replace the data of this handler with the vector received from
   * another process. Requires MPI.
   *
   * \param[in] dest The rank to send to.
   * \param[in] source The rank to receive from.
   * \param[in] tag The tag of the messages.
   * \param[in] comm The MPI communicator used for the communication.
   * \param[out] status The MPI status of the received message.
   * \param[out] outcount The number of received items.
   */
  void
  sendrecv (const int dest, const int source, const int tag, sc_MPI_Comm comm, sc_MPI_Status *status, int &outcount);
};

/**
 * A vector handler for a std::vector with one item of type TType per entry.
 *
 * Single items are packed and unpacked by a data handler of type TDataHandler. By default this is the
 * \ref t8_data_handler of TType, which handles trivially copyable types without user code. For data that has
 * to be serialized, e.g. items holding a std::vector, pass a handler that describes how one item is packed.
 * Ranges of items with a bytewise data handler (\ref T8BytewiseDataHandler) are packed in bulk, all other items
 * are serialized one by one.
 *
 * \tparam TType The type of the data items.
 * \tparam TDataHandler The type of the handler for single items, see \ref T8DataHandler.
 */
template <typename TType, typename TDataHandler = t8_data_handler<TType>>
  requires T8DataHandler<TDataHandler, TType>
struct t8_vector_handler: public t8_abstract_vector_handler
{
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
   * Compute an upper bound of the number of bytes needed by \ref pack_vector_prefix.
   *
   * \param[in] comm The MPI communicator used for the operation.
   * \return The size of the buffer.
   */
  int
  buffer_size (sc_MPI_Comm comm) const override
  {
    int prefix_size = 0;
    const int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &prefix_size);
    SC_CHECK_MPI (mpiret);
    return t8_mpi_int (prefix_size + range_buffer_size (0, m_data.size (), comm));
  }

  /**
   * Packs the vector into a buffer. The data is prefixed with the number of items in the vector.
   *
   * \param[in, out] buffer The buffer to pack the vector into.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  pack_vector_prefix (void *buffer, const int num_bytes, int &position, sc_MPI_Comm comm) const override
  {
    const int num_items = t8_mpi_int (m_data.size ());
    const int mpiret = sc_MPI_Pack (&num_items, 1, sc_MPI_INT, buffer, num_bytes, &position, comm);
    SC_CHECK_MPI (mpiret);
    pack_range (0, m_data.size (), buffer, num_bytes, position, comm);
  }

  /**
   * Unpacks a vector from a buffer and replaces the data of this handler.
   * The data is expected to be prefixed with the number of items in the vector.
   *
   * \param[in] buffer The buffer containing the packed data.
   * \param[in] num_bytes The number of bytes of \a buffer.
   * \param[in, out] position The current position in \a buffer. Advanced by the unpacked bytes.
   * \param[out] outcount The number of unpacked items.
   * \param[in] comm The MPI communicator used for the operation.
   */
  void
  unpack_vector_prefix (const void *buffer, const int num_bytes, int &position, int &outcount,
                        sc_MPI_Comm comm) override
  {
    const int mpiret = sc_MPI_Unpack (buffer, num_bytes, &position, &outcount, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (outcount >= 0);

    m_data.resize (outcount);
    unpack_range (buffer, num_bytes, position, m_data.data (), m_data.size (), comm);
  }

  /**
   * Get the type of the data handler.
   *
   * \return The type of the data handler, \ref t8_data_handler_type_unregistered if the data handler has no
   *         member type ().
   */
  t8_data_handler_type
  type () const override
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
 * \return A handler for the given data type, nullptr if \a type is not an internal type.
 */
inline std::unique_ptr<t8_abstract_vector_handler>
t8_create_internal_handler (const t8_data_handler_type type)
{
  switch (type.get ()) {
    /* Place holder to create a handler for internal data structures */
  default:
    return nullptr;
  }
}

#endif /* T8_VECTOR_HANDLER_HXX */
