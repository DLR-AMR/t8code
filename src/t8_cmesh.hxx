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

/** \file t8_cmesh.hxx
 * We define the coarse mesh of trees in this file.
 */

#ifndef T8_CMESH_HXX
#define T8_CMESH_HXX

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_handler.hxx>

/**
 * Compute the first element of a process in a partitioned mesh, via floor(process * mpisize / global_num_elements).
 * Prevents overflow by splitting global_num_elements into two 32 bit integers.
 * global_num_elements = a_0 * tau + a_1, tau = 2^32.
 * 
 * floor(a_0 * tau + a_1) * process / mpisize
 * = floor (a_0 * tau * process / mpisize + a_1 * process / mpisize)
 * (add the fractional part of the summands and take the floor) to split the computation into separate parts.
 * = floor (a_0 * tau * process / mpisize) + floor (a_1 * process / mpisize) + floor ((a_o*tau*process) % mpisize + a_1*process % mpisize) / mpisize)
 * the first term can computed as:
 * floor (a_0 * tau * process / mpisize)
 * = floor (floor (a_0 * tau / mpisize * process ) + (a_0 * tau % mpisize) * process / mpisize)
 * We divide a_0 * tau by mpisize first to prevent an overflow. The second term corrects for the rounding error.
 *
 * Adding all summands gives the result.
 *
 * This enables us to partition 2^64-1 elements over 2^32-1 processes.
 * Update this function if we have supercomputers with more than 2^32-1 processes, or need larger meshes.
 *
 * \param[in] process   The number of processes
 * \param[in] mpisize   The size of the MPI communicator
 * \param[in] global_num_elements   The number of elements in the global mesh
 * \return The first element of the process in the partitioned mesh
 * 
 * \warning This function assumes that process <= mpisize. mpisize has to be greater than 0.
 * Otherwise the result of process * global_num_elements / mpisize can exceed the range of uint64_t.
 */
constexpr uint64_t
t8_cmesh_get_first_element_of_process (const uint32_t process, const uint32_t mpisize,
                                       const uint64_t global_num_elements)
{
  T8_ASSERT (mpisize > 0);
  T8_ASSERT (process <= mpisize);
  /* Split the uint64_t */
  const uint64_t a_0 = (global_num_elements >> 32);
  const uint64_t a_1 = (global_num_elements << 32) >> 32;

  /* Cast everything into uint64_t */
  const uint64_t process_64 = static_cast<uint64_t> (process);
  const uint64_t mpisize_64 = static_cast<uint64_t> (mpisize);

  /* sum_0 = floor(a_0 * tau / mpisize * process) */
  const uint64_t sum_0 = (a_0 << 32) / mpisize_64 * process_64;

  /* sum_1 = floor(a_0 * tau % mpisize * process / mpisize) */
  const uint64_t sum_1 = (((a_0 << 32) % mpisize_64) * process_64) / mpisize_64;

  /* sum_2 = floor(a_1 * process / mpisize) */
  const uint64_t sum_2 = (a_1 * process_64) / mpisize_64;

  /* sum_3 is the correction term of splitting the floor operation of the sum into two parts.  */
  const uint64_t sum_3_1 = ((a_0 << 32) % mpisize_64 * (process_64 % mpisize_64)) % mpisize_64;
  const uint64_t sum_3_2 = (a_1 * process_64) % mpisize_64;
  const uint64_t sum_3 = (sum_3_1 + sum_3_2) / mpisize_64;

  return (sum_0 + sum_1 + sum_2 + sum_3);
}

/**
 * Create and register a geometry with the coarse mesh. The coarse mesh takes the ownership of the geometry.
 * @tparam geometry_type 
 * \param [in,out] cmesh The cmesh.
 * \param [in,out] args The constructor arguments of the geometry.
 * \return         A pointer to the geometry.
 */

template <typename geometry_type, typename... _args>
inline geometry_type *
t8_cmesh_register_geometry (t8_cmesh_t cmesh, _args &&...args)
{
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not constructed, do it now. */
    cmesh->geometry_handler = new t8_geometry_handler ();
  }
  return cmesh->geometry_handler->register_geometry<geometry_type> (std::forward<_args> (args)...);
}

#endif /* T8_CMESH_HXX */
