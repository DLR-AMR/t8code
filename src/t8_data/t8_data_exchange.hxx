/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
 * \file t8_data_exchange.hxx
 *
 * Building blocks to exchange distributed arrays of items between processes: the partition pattern of a
 * repartitioning, computed from global offsets, nonblocking byte messages, the repartitioning of plain bytes
 * without intermediate buffers, and the ghost pattern with the ghost exchange of plain bytes.
 * They depend neither on the type of the items nor on the forest, so containers like \ref t8_vector_handler and
 * forest-internal arrays can share them.
 */

#ifndef T8_DATA_EXCHANGE_HXX
#define T8_DATA_EXCHANGE_HXX

#include <t8.h>
#include <limits>
#include <memory>
#include <vector>

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
 * A contiguous range of items that this process exchanges with one other process, e.g. during a repartitioning or
 * a ghost exchange.
 */
struct t8_data_partition_chunk
{
  int rank;     /**< The process the items are sent to or received from. */
  size_t first; /**< Local index of the first item: in the old data for sends, in the new data for receives. */
  size_t count; /**< The number of items. Always positive. */
};

/**
 * Describes how this process takes part in repartitioning a distributed array whose items keep their global
 * order, e.g. one item per leaf element when a forest is repartitioned.
 *
 * Both sides of every message decide from the same global offsets whether the message exists and how many items
 * it contains. Hence all messages can be posted at once and no message size has to be probed for plain bytes.
 */
struct t8_data_partition_pattern
{
  std::vector<t8_data_partition_chunk> sends; /**< Chunks of the old data to send, ascending ranks, without self. */
  std::vector<t8_data_partition_chunk> recvs; /**< Chunks of the new data to receive, ascending ranks, without self. */
  size_t self_first_old = 0; /**< Local index in the old data of the first item that stays on this process. */
  size_t self_first_new = 0; /**< Local index in the new data of the first item that stays on this process. */
  size_t self_count = 0;     /**< The number of items that stay on this process. */
  size_t num_old = 0;        /**< The number of local items before the repartitioning. */
  size_t num_new = 0;        /**< The number of local items after the repartitioning. */
  /** The number of items that new storage for the new data has room for at least, e.g. to append the ghost items
   * of the new forest without reallocation. Containers reserve max (num_new, min_capacity_new) items. */
  size_t min_capacity_new = 0;
  sc_MPI_Comm comm = sc_MPI_COMM_NULL; /**< The communicator of the repartitioning. */
  int tag = 0;                         /**< The MPI tag of all messages. */
};

/**
 * Compute the partition pattern of this process from the old and new global offsets of a distributed array.
 * The items of process p are [offsets[p], offsets[p + 1]), offsets[0] = 0 and offsets[num_ranks] is the global
 * number of items, which has to be equal in both distributions. Empty processes are allowed.
 * The cost is O(log (num_ranks) + number of partners), independent of the number of items.
 *
 * \param[in] offsets_old The old offsets, \a num_ranks + 1 entries.
 * \param[in] offsets_new The new offsets, \a num_ranks + 1 entries.
 * \param[in] rank The rank of this process.
 * \param[in] num_ranks The number of processes.
 * \param[in] comm The communicator of the repartitioning.
 * \param[in] tag The MPI tag of all messages.
 * \return The partition pattern of \a rank.
 */
t8_data_partition_pattern
t8_data_partition_pattern_from_offsets (const t8_gloidx_t *offsets_old, const t8_gloidx_t *offsets_new, const int rank,
                                        const int num_ranks, sc_MPI_Comm comm, const int tag);

/**
 * Check that the chunks of a partition pattern cover the old and the new data exactly once.
 * Meant for assertions at the start of a repartitioning.
 *
 * \param[in] pattern A partition pattern.
 * \return True if the chunks are consistent with num_old and num_new, false otherwise.
 */
bool
t8_data_partition_pattern_is_valid (const t8_data_partition_pattern &pattern);

/**
 * The nonblocking messages of one communication step with plain bytes, e.g. a repartitioning.
 * A chunk with more than the maximum number of bytes per message is split into consecutive messages; both sides
 * split it identically, since both know its size. All messages have to be posted before \ref wait_all, so that no
 * process blocks before its messages exist.
 */
struct t8_data_messages
{
  /**
   * Prepare the messages of a communication step.
   *
   * \param[in] comm The communicator of the messages.
   * \param[in] tag The MPI tag of all messages.
   * \param[in] max_message_bytes Chunks with more bytes are split into several messages. Only meant to be
   *                              changed by tests; the default is the largest message size MPI can describe.
   */
  t8_data_messages (sc_MPI_Comm comm, const int tag,
                    const size_t max_message_bytes = static_cast<size_t> (std::numeric_limits<int>::max ()));

  /** The messages are not copyable, since MPI requests cannot be shared. */
  t8_data_messages (const t8_data_messages &) = delete;

  /** The messages are not copyable, since MPI requests cannot be shared. \return This object. */
  t8_data_messages &
  operator= (const t8_data_messages &)
    = delete;

  /** All posted messages have to be completed with \ref wait_all before destruction. */
  ~t8_data_messages ();

  /**
   * Post the receive of a chunk of bytes.
   *
   * \param[in] rank The process to receive from.
   * \param[out] buffer Space for \a num_bytes bytes, filled when \ref wait_all returns.
   * \param[in] num_bytes The number of bytes of the chunk.
   */
  void
  receive (const int rank, void *buffer, const size_t num_bytes);

  /**
   * Post the send of a chunk of bytes.
   *
   * \param[in] rank The process to send to.
   * \param[in] buffer The \a num_bytes bytes to send. Must not change until \ref wait_all returns.
   * \param[in] num_bytes The number of bytes of the chunk.
   */
  void
  send (const int rank, const void *buffer, const size_t num_bytes);

  /**
   * Wait until all posted messages are complete. In debug mode, check that every receive got the expected bytes.
   */
  void
  wait_all ();

 private:
  sc_MPI_Comm m_comm;                             /**< The communicator of the messages. */
  int m_tag;                                      /**< The MPI tag of all messages. */
  size_t m_max_message_bytes;                     /**< The maximum number of bytes of one message. */
  std::vector<sc_MPI_Request> m_receive_requests; /**< The requests of the posted receives. */
  std::vector<int> m_receive_bytes;               /**< The expected number of bytes of each posted receive. */
  std::vector<sc_MPI_Request> m_send_requests;    /**< The requests of the posted sends. */
};

/**
 * Repartition an array of items that are plain bytes, without packing and without intermediate buffers.
 * All receives are posted directly into \a data_out first, then all sends directly from \a data_in. The items
 * that stay on this process are copied while the messages are in flight.
 *
 * \param[in] pattern The partition pattern of this process.
 * \param[in] data_in The old data, \a pattern.num_old items of \a item_size bytes. Must not overlap \a data_out.
 * \param[out] data_out The new data, space for \a pattern.num_new items of \a item_size bytes.
 * \param[in] item_size The number of bytes of one item. Must be equal on all processes.
 * \param[in] max_message_bytes Chunks with more bytes are split into several messages. Only meant to be
 *                              changed by tests; the default is the largest message size MPI can describe.
 */
void
t8_data_partition_bytes (const t8_data_partition_pattern &pattern, const void *data_in, void *data_out,
                         const size_t item_size,
                         const size_t max_message_bytes = static_cast<size_t> (std::numeric_limits<int>::max ()));

/**
 * Describes how this process takes part in exchanging the data of ghost items, e.g. of the ghost elements of a forest.
 * The data holds the local items, followed by the ghost items. Every process sends some of its local items to the
 * processes for which they are ghosts, and receives its ghost items from their owners into contiguous ranges.
 */
struct t8_data_ghost_pattern
{
  std::vector<int> send_ranks;      /**< The processes to send to, ascending, without empty sends. */
  std::vector<size_t> send_offsets; /**< The items sent to send_ranks[i] are send_indices[send_offsets[i]] to
                                         send_indices[send_offsets[i + 1] - 1]. One entry more than send_ranks,
                                         first entry 0. */
  std::vector<size_t> send_indices; /**< The local indices of the items to send, for each process in the order of
                                         its ghost items. */
  std::vector<t8_data_partition_chunk> recvs; /**< The ghost items to receive: the sending process, the index of the
                                                   first one in the data and their number. Ascending ranks, without
                                                   empty chunks. */
  size_t num_local = 0;                       /**< The number of local items. */
  size_t num_ghosts = 0;                      /**< The number of ghost items. */
  sc_MPI_Comm comm = sc_MPI_COMM_NULL;        /**< The communicator of the exchange. */
  int tag = 0;                                /**< The MPI tag of all messages. */
};

/**
 * Check that a ghost pattern is consistent: the send offsets match the send ranks and indices, only local items are
 * sent, and the received chunks lie within the ghost items.
 * Meant for assertions at the start of a ghost exchange.
 *
 * \param[in] pattern A ghost pattern.
 * \return True if the pattern is consistent, false otherwise.
 */
bool
t8_data_ghost_pattern_is_valid (const t8_data_ghost_pattern &pattern);

/**
 * A byte buffer for messages that is reused, e.g. across the ghost exchanges of one forest, so that the messages do
 * not pay for fresh memory each time. It only grows and keeps its largest size until destruction.
 */
struct t8_data_buffer
{
  /**
   * Get the buffer with room for at least \a num_bytes bytes. Growing the buffer discards its content and
   * invalidates pointers returned before.
   *
   * \param[in] num_bytes The number of bytes needed.
   * \return The buffer. Its content is uninitialized after growing.
   */
  char *
  get (const size_t num_bytes)
  {
    if (num_bytes > m_num_bytes) {
      /* Free the old buffer first, so that the old and the new one do not occupy memory at the same time. */
      m_bytes.reset ();
      m_bytes = std::make_unique_for_overwrite<char[]> (num_bytes);
      m_num_bytes = num_bytes;
    }
    return m_bytes.get ();
  }

 private:
  std::unique_ptr<char[]> m_bytes; /**< The buffer, nullptr before the first call of \ref get. */
  size_t m_num_bytes = 0;          /**< The number of bytes of the buffer. */
};

/**
 * A nonblocking ghost exchange of items that are plain bytes.
 * The constructor posts all messages: the receives directly into the ghost slots of the data, the sends from a
 * gather buffer that collects the items for all processes. Each process gets its items as soon as they are gathered.
 * \ref wait completes the exchange; in between, the local items may be read but no item may be modified.
 */
struct t8_data_ghost_bytes_exchange
{
  /**
   * Post the ghost exchange.
   *
   * \param[in] pattern The ghost pattern of this process.
   * \param[in,out] data The \a pattern.num_local local items followed by space for the \a pattern.num_ghosts ghost
   *                     items, each of \a item_size bytes. The ghost items are filled when \ref wait returns.
   * \param[in] item_size The number of bytes of one item. Must be equal on all processes.
   * \param[in,out] send_buffer The gather buffer of the items to send. It must not be used otherwise until
   *                            \ref wait returns.
   * \param[in] max_message_bytes Larger messages are split into several messages. Only meant to be changed by
   *                              tests; the default is the largest message size MPI can describe.
   */
  t8_data_ghost_bytes_exchange (const t8_data_ghost_pattern &pattern, void *data, const size_t item_size,
                                t8_data_buffer &send_buffer,
                                const size_t max_message_bytes
                                = static_cast<size_t> (std::numeric_limits<int>::max ()));

  /**
   * Wait until the exchange is complete.
   */
  void
  wait ();

 private:
  t8_data_messages m_messages; /**< The posted messages. */
};

#endif /* T8_DATA_EXCHANGE_HXX */
