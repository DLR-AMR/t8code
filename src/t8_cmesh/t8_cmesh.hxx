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

#include <vector>
#include <unordered_map>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity_types.hxx>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_connectivity_types.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

/**
 * Compute the first element of a process in a partitioned mesh, via floor(process * global_num_elements / mpisize).
 * Prevents overflowing by using division with remainder to store partial results in uint64_t.
 * We also take into account that process <= mpisize.
 * Both process and global_num_elements can be written as:
 * 
 * global_num_elements = elem_over_size * mpisize + remainder_0
 * process = proc_over_size * mpisize + remainder_1
 * with:
 * elem_over_size = global_num_elements / mpisize
 * proc_over_size = process / mpisize
 * and remainders:
 * remainder_0 = global_num_elements % mpisize
 * remainder_1 = process % mpisize
 * 
 * Putting this together in the computation of first_element = global_num_elements * process / mpisize gives:
 * 
 * = elem_over_size * proc_over_size * mpisize + elem_over_size * process % mpisize
 *   + proc_over_size * global_num_elements % mpisize + (remainder_0 * remainder_1) / mpisize
 *
 * Each variable is less than 2^32-1 taking into account that process <= mpisize 
 * we can assure that the first summand does not overflow
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

  /* Cast everything into uint64_t */
  const uint64_t process_64 = static_cast<uint64_t> (process);
  const uint64_t mpisize_64 = static_cast<uint64_t> (mpisize);

  /* Split the uint64_t */
  const uint64_t elem_over_size = global_num_elements / mpisize_64;
  const uint64_t remainder_0 = global_num_elements % mpisize_64;

  const uint64_t proc_over_size = process_64 / mpisize_64;
  const uint64_t remainder_1 = process_64 % mpisize_64;

  const uint64_t sum_0 = (elem_over_size * proc_over_size) * mpisize_64;
  const uint64_t sum_1 = elem_over_size * (process_64 % mpisize_64);
  const uint64_t sum_2 = proc_over_size * (global_num_elements % mpisize_64);
  const uint64_t sum_3 = (remainder_0 * remainder_1) / mpisize_64;

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
geometry_type *
t8_cmesh_register_geometry (t8_cmesh_t cmesh, _args &&...args)
{
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not constructed, do it now. */
    cmesh->geometry_handler = new t8_geometry_handler ();
  }
  return cmesh->geometry_handler->register_geometry<geometry_type> (std::forward<_args> (args)...);
}


/** Get the list of global trees and local vertex ids a global vertex is connected to.
 * Cmesh Interface function.
 *  
 * \param [in] cmesh A committed cmesh.
 * \param [in] global_vertex The global id of a vertex in \a cmesh.
 * \return The list of global tree ids and local vertex ids of \a global_vertex_id.
 */
const tree_vertex_list &
t8_cmesh_get_vertex_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_vertex);

const tree_edge_list &
t8_cmesh_get_edge_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_edge);

typedef class t8_neigh_info {
 public:
  t8_gloidx_t neighid;
  int neigh_bdy_id;
  int orientation;
  int sign;
} t8_neigh_info;

std::vector<t8_neigh_info>
t8_cmesh_get_neighs (t8_cmesh_t cmesh, t8_locidx_t treeid, int bdy_dim, int bdy_id);

#endif /* T8_CMESH_HXX */
