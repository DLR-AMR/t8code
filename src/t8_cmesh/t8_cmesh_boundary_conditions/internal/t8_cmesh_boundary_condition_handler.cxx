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
 * \file t8_cmesh_boundary_condition_handler.cxx
 * Implementation context of \ref t8_cmesh_boundary_condition_handler.hxx
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

#include <cstring>

using namespace detail;

#if T8_ENABLE_DEBUG
int
t8_cmesh_boundary_condition_handler::verify () const
{
  T8_ASSERT (!t8_cmesh_is_committed (m_cmesh, 0));

  /* Retrieve eclasses and boundary conditions from stash */
  std::vector<std::pair<t8_gloidx_t, std::span<boundary_condition_hash>>> boundary_conditions
    = t8_stash_extract_attribute_list<boundary_condition_hash> (m_cmesh->stash, t8_get_package_id (),
                                                                T8_CMESH_BOUNDARY_CONDITION_ATTRIBUTE_KEY);
  std::vector<std::pair<t8_gloidx_t, t8_eclass>> eclasses = t8_stash_extract_eclasses (m_cmesh->stash);

  /* Check that every tree has a bc */
  if (eclasses.size () != boundary_conditions.size ()) {
    t8_errorf ("ERROR: Number of cmesh eclasses does not match the number of cmesh boundary conditions.\n"
               "If boundary conditions are applied, all trees have to get them assigned. They can be left empty.\n");
    return false;
  }

  /* Sort the bcs and eclasses by tree id. */
  std::sort (boundary_conditions.begin (), boundary_conditions.end (),
             [] (const auto &a, const auto &b) { return a.first < b.first; });
  std::sort (eclasses.begin (), eclasses.end (), [] (const auto &a, const auto &b) { return a.first < b.first; });

  auto eclasses_it = eclasses.cbegin ();
  auto bc_it = boundary_conditions.cbegin ();
  for (; eclasses_it != eclasses.end (); ++eclasses_it, ++bc_it) {
    /* Make sure, that the eclasses and bcs tree id are the same.  */
    if (eclasses_it->first != bc_it->first) {
      t8_errorf ("ERROR: eclass tree id and boundary condition tree id do not match.\nProbably the boundary condition "
                 "got a wrong tree id. Boundary condition tree id: %li \n",
                 bc_it->first);
      return false;
    }
    /* Make sure that the amount of registered bcs matches th number of tree faces. */
    if (static_cast<size_t> (t8_eclass_num_faces[eclasses_it->second]) != bc_it->second.size ()) {
      t8_errorf ("ERROR: Tree %li has a mismatch in its face count and boundary condition count.\n"
                 "Face count: %i, boundary condition count: %li\n",
                 eclasses_it->first, t8_eclass_num_faces[eclasses_it->second], bc_it->second.size ());
      return false;
    }
    /* Check, that the hash of the bc attribute is present in the bc map. This failure cannot be caused by the user. */
    for (const auto bc : bc_it->second) {
      if (!get_boundary_condition_name_safe (bc).has_value ()) {
        t8_errorf ("ERROR: Could not find name of applied boundary condition. You probably found a bug.\n");
        return false;
      }
    }
  }
  return true;
}

#endif /* T8_ENABLE_DEBUG */

int
t8_cmesh_boundary_condition_handler::get_boundary_condition_attribute_key () const
{
  return T8_CMESH_BOUNDARY_CONDITION_ATTRIBUTE_KEY;
}

std::vector<char>
t8_cmesh_boundary_condition_handler::serialize_map () const
{
  std::vector<char> serial_data;
  /* Fill serial_data with strings only. The hashes can be re-generated locally. */
  for (const auto &[key, string] : t8_cmesh_boundary_condition_handler::m_boundary_conditions) {
    serial_data.insert (serial_data.end (), string.begin (), string.end ());
    /* Null terminate strings to be able to split them again later. */
    serial_data.push_back ('\0');
  }
  serial_data.shrink_to_fit ();
  return serial_data;
}

void
t8_cmesh_boundary_condition_handler::unpack_map (std::vector<char> &serial_data, bool overwrite)
{
  if (overwrite) {
    m_boundary_conditions.clear ();
  }

  /* Iterate over stings. */
  for (const char *string = serial_data.data (); string < serial_data.data () + serial_data.size ();
       string += std::strlen (string) + 1) {
    /* Interpret c string as std::string, rehash it and insert it. */
    std::string value (string);
    t8_cmesh_boundary_condition_handler::m_boundary_conditions.try_emplace (
      t8_cmesh_boundary_condition_handler::hash_boundary_condition_name (value), std::move (value));
  }
}

void
t8_cmesh_boundary_condition_handler::synchronize (sc_MPI_Comm comm)
{
  /* Use a bottom-up binomial tree merge approach instead of an MPI_Allgatherv to secure O(log(p)) scaling.
   * In every level of the merge tree each rank with rank = rank & ~mask merges all information of rank = rank | mask.
   * The rank = rank && ~mask then drops out of the communication pattern.
   * The receiving rank also checks if there is a sender in the first place (src >= mpisize) for non-power of 2 mpisizes.
   *
   * After all data is collected, rank 0 broadcasts the collected data.
   */

  T8_ASSERT (comm != sc_MPI_COMM_NULL);

  int rank, mpisize;
  sc_MPI_Comm_rank (comm, &rank);
  sc_MPI_Comm_size (comm, &mpisize);

  if (mpisize == 1) {
    /* Nothing to do. */
    return;
  }

  /* Iterate over all levels of the merge tree. */
  for (int mask = 1; mask < mpisize; mask <<= 1) {
    /* This rank receives a message. */
    if ((rank & mask) == 0) {
      const int src = rank | mask;
      if (src >= mpisize) {
        /* There is no sender, so we have nothing to do on this level. */
        continue;
      }
      /* Probe the size of the incoming message. */
      sc_MPI_Status status;
      sc_MPI_Probe (src, T8_MPI_BOUNDARY_CONDITION_SYNC_TAG, comm, &status);
      int incoming_bytes;
      sc_MPI_Get_count (&status, sc_MPI_BYTE, &incoming_bytes);

      /* Prepare buffer and receive data. */
      std::vector<char> incoming (incoming_bytes);
      sc_MPI_Recv (incoming.data (), incoming_bytes, sc_MPI_BYTE, src, T8_MPI_BOUNDARY_CONDITION_SYNC_TAG, comm,
                   sc_MPI_STATUS_IGNORE);
      t8_cmesh_boundary_condition_handler::unpack_map (incoming, false);
    }
    /* This process sends a message and then drops out. */
    else {
      const int dst = rank & ~mask;
      const std::vector<char> outgoing = serialize_map ();
      sc_MPI_Send (const_cast<char *> (outgoing.data ()), static_cast<int> (outgoing.size ()), sc_MPI_BYTE, dst,
                   T8_MPI_BOUNDARY_CONDITION_SYNC_TAG, comm);
      /* Drop out. */
      break;
    }
  }

  t8_cmesh_boundary_condition_handler::bcast (0, comm);
}

void
t8_cmesh_boundary_condition_handler::bcast (int main_rank, sc_MPI_Comm comm)
{
  int rank;
  sc_MPI_Comm_rank (comm, &rank);

  /* Buffer for sending and receiving */
  std::vector<char> buffer;

  /* Serialize data. */
  if (rank == main_rank) {
    buffer = t8_cmesh_boundary_condition_handler::serialize_map ();
  }

  /* Prepare buffer length. */
  int size = static_cast<int> (buffer.size ());
  sc_MPI_Bcast (&size, 1, sc_MPI_INT, main_rank, comm);
  if (rank != main_rank) {
    buffer.resize (size);
  }

  /* Broadcast data */
  sc_MPI_Bcast (buffer.data (), size, sc_MPI_BYTE, main_rank, comm);

  /* Add data to local map. */
  if (rank != main_rank) {
    t8_cmesh_boundary_condition_handler::unpack_map (buffer, true);
  }
}
