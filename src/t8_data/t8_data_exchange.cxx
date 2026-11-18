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
 * \file t8_data_exchange.cxx
 *
 * Implementation of the type-independent exchange of distributed arrays, see \ref t8_data_exchange.hxx.
 */

#include <t8_data/t8_data_exchange.hxx>
#include <algorithm>
#include <cstring>

t8_data_partition_pattern
t8_data_partition_pattern_from_offsets (const t8_gloidx_t *offsets_old, const t8_gloidx_t *offsets_new, const int rank,
                                        const int num_ranks, sc_MPI_Comm comm, const int tag)
{
  T8_ASSERT (0 <= rank && rank < num_ranks);
  T8_ASSERT (offsets_old[0] == 0 && offsets_new[0] == 0);
  T8_ASSERT (offsets_old[num_ranks] == offsets_new[num_ranks]);

  const t8_gloidx_t old_first = offsets_old[rank];
  const t8_gloidx_t old_end = offsets_old[rank + 1];
  const t8_gloidx_t new_first = offsets_new[rank];
  const t8_gloidx_t new_end = offsets_new[rank + 1];

  t8_data_partition_pattern pattern;
  pattern.num_old = static_cast<size_t> (old_end - old_first);
  pattern.num_new = static_cast<size_t> (new_end - new_first);
  pattern.comm = comm;
  pattern.tag = tag;

  /* Send the overlap of the old local range with the new range of each process. The first process to consider
   * is the new owner of our first old item; the offsets are sorted, so it is found by binary search. */
  if (old_first < old_end) {
    int iproc
      = static_cast<int> (std::upper_bound (offsets_new, offsets_new + num_ranks + 1, old_first) - offsets_new) - 1;
    for (; iproc < num_ranks && offsets_new[iproc] < old_end; ++iproc) {
      const t8_gloidx_t begin = std::max (old_first, offsets_new[iproc]);
      const t8_gloidx_t end = std::min (old_end, offsets_new[iproc + 1]);
      if (begin >= end) {
        /* iproc is empty in the new partition. */
        continue;
      }
      if (iproc == rank) {
        pattern.self_first_old = static_cast<size_t> (begin - old_first);
        pattern.self_first_new = static_cast<size_t> (begin - new_first);
        pattern.self_count = static_cast<size_t> (end - begin);
      }
      else {
        pattern.sends.push_back ({ iproc, static_cast<size_t> (begin - old_first), static_cast<size_t> (end - begin) });
      }
    }
  }

  /* Receive the overlap of the new local range with the old range of each process. */
  if (new_first < new_end) {
    int iproc
      = static_cast<int> (std::upper_bound (offsets_old, offsets_old + num_ranks + 1, new_first) - offsets_old) - 1;
    for (; iproc < num_ranks && offsets_old[iproc] < new_end; ++iproc) {
      const t8_gloidx_t begin = std::max (new_first, offsets_old[iproc]);
      const t8_gloidx_t end = std::min (new_end, offsets_old[iproc + 1]);
      if (begin >= end || iproc == rank) {
        /* iproc is empty in the old partition, or this is the self chunk computed above. */
        continue;
      }
      pattern.recvs.push_back ({ iproc, static_cast<size_t> (begin - new_first), static_cast<size_t> (end - begin) });
    }
  }
  return pattern;
}

bool
t8_data_partition_pattern_is_valid (const t8_data_partition_pattern &pattern)
{
  size_t num_sent = pattern.self_count;
  for (const auto &chunk : pattern.sends) {
    if (chunk.count == 0 || chunk.first + chunk.count > pattern.num_old) {
      return false;
    }
    num_sent += chunk.count;
  }
  size_t num_received = pattern.self_count;
  for (const auto &chunk : pattern.recvs) {
    if (chunk.count == 0 || chunk.first + chunk.count > pattern.num_new) {
      return false;
    }
    num_received += chunk.count;
  }
  return num_sent == pattern.num_old && num_received == pattern.num_new;
}

t8_data_messages::t8_data_messages (sc_MPI_Comm comm, const int tag, const size_t max_message_bytes)
  : m_comm (comm), m_tag (tag), m_max_message_bytes (max_message_bytes)
{
  T8_ASSERT (0 < max_message_bytes && max_message_bytes <= static_cast<size_t> (std::numeric_limits<int>::max ()));
}

t8_data_messages::~t8_data_messages ()
{
  T8_ASSERT (m_receive_requests.empty () && m_send_requests.empty ());
}

void
t8_data_messages::receive (const int rank, void *buffer, const size_t num_bytes)
{
  char *const bytes = static_cast<char *> (buffer);
  for (size_t offset = 0; offset < num_bytes; offset += m_max_message_bytes) {
    const int message_bytes = t8_mpi_int (std::min (m_max_message_bytes, num_bytes - offset));
    m_receive_requests.emplace_back ();
    const int mpiret
      = sc_MPI_Irecv (bytes + offset, message_bytes, sc_MPI_BYTE, rank, m_tag, m_comm, &m_receive_requests.back ());
    SC_CHECK_MPI (mpiret);
    m_receive_bytes.push_back (message_bytes);
  }
}

void
t8_data_messages::send (const int rank, const void *buffer, const size_t num_bytes)
{
  /* The cast is needed because the non-MPI prototype takes a non-const buffer. */
  char *const bytes = const_cast<char *> (static_cast<const char *> (buffer));
  for (size_t offset = 0; offset < num_bytes; offset += m_max_message_bytes) {
    const int message_bytes = t8_mpi_int (std::min (m_max_message_bytes, num_bytes - offset));
    m_send_requests.emplace_back ();
    const int mpiret
      = sc_MPI_Isend (bytes + offset, message_bytes, sc_MPI_BYTE, rank, m_tag, m_comm, &m_send_requests.back ());
    SC_CHECK_MPI (mpiret);
  }
}

void
t8_data_messages::wait_all ()
{
  if (!m_receive_requests.empty ()) {
#if T8_ENABLE_DEBUG
    /* Check that every receive got exactly the expected number of bytes. */
    std::vector<sc_MPI_Status> statuses (m_receive_requests.size ());
    int mpiret = sc_MPI_Waitall (t8_mpi_int (m_receive_requests.size ()), m_receive_requests.data (), statuses.data ());
    SC_CHECK_MPI (mpiret);
    for (size_t irequest = 0; irequest < statuses.size (); ++irequest) {
      int received_bytes = 0;
      mpiret = sc_MPI_Get_count (&statuses[irequest], sc_MPI_BYTE, &received_bytes);
      SC_CHECK_MPI (mpiret);
      T8_ASSERT (received_bytes == m_receive_bytes[irequest]);
    }
#else
    const int mpiret
      = sc_MPI_Waitall (t8_mpi_int (m_receive_requests.size ()), m_receive_requests.data (), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
#endif
  }
  if (!m_send_requests.empty ()) {
    const int mpiret
      = sc_MPI_Waitall (t8_mpi_int (m_send_requests.size ()), m_send_requests.data (), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  m_receive_requests.clear ();
  m_receive_bytes.clear ();
  m_send_requests.clear ();
}

void
t8_data_partition_bytes (const t8_data_partition_pattern &pattern, const void *data_in, void *data_out,
                         const size_t item_size, const size_t max_message_bytes)
{
  T8_ASSERT (t8_data_partition_pattern_is_valid (pattern));
  T8_ASSERT (item_size > 0);
  const char *const bytes_in = static_cast<const char *> (data_in);
  char *const bytes_out = static_cast<char *> (data_out);
  T8_ASSERT (pattern.num_old == 0 || pattern.num_new == 0 || bytes_in + pattern.num_old * item_size <= bytes_out
             || bytes_out + pattern.num_new * item_size <= bytes_in);
  t8_data_messages messages (pattern.comm, pattern.tag, max_message_bytes);

  /* Post the receives first, so that incoming messages land directly at their final position. */
  for (const auto &chunk : pattern.recvs) {
    messages.receive (chunk.rank, bytes_out + chunk.first * item_size, chunk.count * item_size);
  }
  /* Send directly from the old data. */
  for (const auto &chunk : pattern.sends) {
    messages.send (chunk.rank, bytes_in + chunk.first * item_size, chunk.count * item_size);
  }

  /* Copy the items that stay on this process while the messages are in flight. */
  if (pattern.self_count > 0) {
    std::memcpy (bytes_out + pattern.self_first_new * item_size, bytes_in + pattern.self_first_old * item_size,
                 pattern.self_count * item_size);
  }
  messages.wait_all ();
}

bool
t8_data_ghost_pattern_is_valid (const t8_data_ghost_pattern &pattern)
{
  if (pattern.send_offsets.size () != pattern.send_ranks.size () + 1 || pattern.send_offsets.front () != 0
      || pattern.send_offsets.back () != pattern.send_indices.size ()) {
    return false;
  }
  for (size_t isend = 0; isend < pattern.send_ranks.size (); ++isend) {
    if (pattern.send_offsets[isend] >= pattern.send_offsets[isend + 1]) {
      return false;
    }
  }
  for (const size_t index : pattern.send_indices) {
    if (index >= pattern.num_local) {
      return false;
    }
  }
  for (const auto &chunk : pattern.recvs) {
    if (chunk.count == 0 || chunk.first < pattern.num_local
        || chunk.first + chunk.count > pattern.num_local + pattern.num_ghosts) {
      return false;
    }
  }
  return true;
}

t8_data_ghost_bytes_exchange::t8_data_ghost_bytes_exchange (const t8_data_ghost_pattern &pattern, void *data,
                                                            const size_t item_size, t8_data_buffer &send_buffer,
                                                            const size_t max_message_bytes)
  : m_messages (pattern.comm, pattern.tag, max_message_bytes)
{
  T8_ASSERT (t8_data_ghost_pattern_is_valid (pattern));
  T8_ASSERT (item_size > 0);
  char *const bytes = static_cast<char *> (data);

  /* Post the receives first, so that incoming messages land directly in the ghost slots. */
  for (const auto &chunk : pattern.recvs) {
    m_messages.receive (chunk.rank, bytes + chunk.first * item_size, chunk.count * item_size);
  }

  /* Gather the items for one process after another and send them right away. */
  char *const gathered = send_buffer.get (pattern.send_indices.size () * item_size);
  for (size_t isend = 0; isend < pattern.send_ranks.size (); ++isend) {
    const size_t first = pattern.send_offsets[isend];
    const size_t end = pattern.send_offsets[isend + 1];
    for (size_t iitem = first; iitem < end; ++iitem) {
      std::memcpy (gathered + iitem * item_size, bytes + pattern.send_indices[iitem] * item_size, item_size);
    }
    m_messages.send (pattern.send_ranks[isend], gathered + first * item_size, (end - first) * item_size);
  }
}

void
t8_data_ghost_bytes_exchange::wait ()
{
  m_messages.wait_all ();
}
