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
 * \file t8_vector_handler.cxx
 *
 * Implementation of the item-type independent communication of \ref t8_abstract_vector_handler.
 */

#include <t8_data/t8_vector_handler.hxx>
#include <memory>

void
t8_abstract_vector_handler::send ([[maybe_unused]] const int dest, [[maybe_unused]] const int tag,
                                  [[maybe_unused]] sc_MPI_Comm comm) const
{
#if T8_ENABLE_MPI
  const int num_bytes = buffer_size (comm);
  /* The buffer is completely overwritten by packing, so it is not initialized. */
  const std::unique_ptr<char[]> buffer (new char[num_bytes]);
  int position = 0;
  pack_vector_prefix (buffer.get (), num_bytes, position, comm);

  /* buffer_size is an upper bound, only the packed bytes are sent. */
  const int mpiret = sc_MPI_Send (buffer.get (), position, sc_MPI_PACKED, dest, tag, comm);
  SC_CHECK_MPI (mpiret);
#else
  SC_ABORT ("t8_abstract_vector_handler::send requires MPI.");
#endif
}

void
t8_abstract_vector_handler::recv ([[maybe_unused]] const int source, [[maybe_unused]] const int tag,
                                  [[maybe_unused]] sc_MPI_Comm comm, [[maybe_unused]] sc_MPI_Status *status,
                                  [[maybe_unused]] int &outcount)
{
#if T8_ENABLE_MPI
  int mpiret = sc_MPI_Probe (source, tag, comm, status);
  SC_CHECK_MPI (mpiret);
  int num_bytes = 0;
  mpiret = sc_MPI_Get_count (status, sc_MPI_PACKED, &num_bytes);
  SC_CHECK_MPI (mpiret);

  const std::unique_ptr<char[]> buffer (new char[num_bytes]);
  mpiret = sc_MPI_Recv (buffer.get (), num_bytes, sc_MPI_PACKED, source, tag, comm, status);
  SC_CHECK_MPI (mpiret);

  int position = 0;
  unpack_vector_prefix (buffer.get (), num_bytes, position, outcount, comm);
#else
  SC_ABORT ("t8_abstract_vector_handler::recv requires MPI.");
#endif
}

void
t8_abstract_vector_handler::sendrecv ([[maybe_unused]] const int dest, [[maybe_unused]] const int source,
                                      [[maybe_unused]] const int tag, [[maybe_unused]] sc_MPI_Comm comm,
                                      [[maybe_unused]] sc_MPI_Status *status, [[maybe_unused]] int &outcount)
{
#if T8_ENABLE_MPI
  /* Pack the outgoing message first, so that only the packed bytes are announced and sent. */
  const int send_buffer_size = buffer_size (comm);
  const std::unique_ptr<char[]> send_buffer (new char[send_buffer_size]);
  int send_num_bytes = 0;
  pack_vector_prefix (send_buffer.get (), send_buffer_size, send_num_bytes, comm);

  /* Exchange the sizes of the outgoing and incoming messages. */
  int recv_num_bytes = 0;
  int mpiret = sc_MPI_Sendrecv (&send_num_bytes, 1, sc_MPI_INT, dest, tag, &recv_num_bytes, 1, sc_MPI_INT, source, tag,
                                comm, status);
  SC_CHECK_MPI (mpiret);

  /* Exchange the payload. */
  const std::unique_ptr<char[]> recv_buffer (new char[recv_num_bytes]);
  mpiret = sc_MPI_Sendrecv (send_buffer.get (), send_num_bytes, sc_MPI_PACKED, dest, tag, recv_buffer.get (),
                            recv_num_bytes, sc_MPI_PACKED, source, tag, comm, status);
  SC_CHECK_MPI (mpiret);

  int position = 0;
  unpack_vector_prefix (recv_buffer.get (), recv_num_bytes, position, outcount, comm);
#else
  SC_ABORT ("t8_abstract_vector_handler::sendrecv requires MPI.");
#endif
}
