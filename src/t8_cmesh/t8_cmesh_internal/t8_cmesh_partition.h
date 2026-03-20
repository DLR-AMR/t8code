/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_partition.h
 * Functionality related to the partitioning of a cmesh.
 */

#ifndef T8_CMESH_PARTITION_H
#define T8_CMESH_PARTITION_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>

T8_EXTERN_C_BEGIN ();

/** Given a cmesh which is to be partitioned, execute the partition task.
 *  This includes partitioning by uniform level and partitioning from a second cmesh
 *  \param [in,out]  cmesh   The cmesh to be partitioned
 *  \param [in]      comm    The MPI communicator
 *  This function is usually called within \ref t8_cmesh_commit */
void
t8_cmesh_partition (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** From num_local_trees_per_eclass compute num_trees_per_eclass.
 *  This function is collective.
 * \param [in,out]    cmesh  The cmesh whose num_trees_per_eclass values should be created.
 *                           Must be partitioned and committed.
 * \param [in]        comm   Mpi communicator used to create the offset array.
 * \warning This function does not perform a check whether \a cmesh is committed.
 * Use with caution.
 */
void
t8_cmesh_gather_trees_per_eclass (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Create the offset array for a partitioned cmesh.
 * This function is collective.
 * \param [in,out]    cmesh  The cmesh whose array should be created.
 *                           Must be partitioned and committed.
 * \param [in]        comm   Mpi communicator used to create the offset array.
 * \note if the offset array (cmesh->tree_offsets) already exists, it is not changed.
 */
void
t8_cmesh_gather_treecount (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Perform the same task as \ref t8_cmesh_gather_treecount, but do
 * not perform the debugging check whether cmesh is committed.
 * \warning Use with caution and only if you know what you are doing.
 * Prefer \ref t8_cmesh_gather_treecount.
 * This function is collective.
 * \param [in,out]    cmesh  The cmesh whose array should be created.
 *                           Must be partitioned and first and last local tree
 *                           as well as the total number of tree must be set.
 * \param [in]        comm   Mpi communicator used to create the offset array.
 * \note if the offset array (cmesh->tree_offsets) already exists, it is not changed.
 */
void
t8_cmesh_gather_treecount_nocommit (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/**
 * Print the offset array of a cmesh to stdout.
 * This function is used for debugging purposes only.
 * \warning The function definition will be empty if not in debug mode!
 * \param [in]      cmesh   A cmesh that is committed and partitioned.
 * \param [in]      comm    A valid MPI communicator for cmesh.
 */
void
t8_cmesh_offset_print (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Create a valid partition table that concentrates all trees at a given
 *  process.
 * \param[in]        proc    The processor that should get all trees.
 * \param[in]        comm    The communicator to use.
 * \param[in]        num_trees The number of global trees in the partition.
 * \return                   A valid partition table for a cmesh with \a num_trees trees
 *                           and communicator \a comm, where each tree is on process \a proc.
 * \note This function is MPI collective.
 */
t8_shmem_array_t
t8_cmesh_offset_concentrate (int proc, sc_MPI_Comm comm, t8_gloidx_t num_trees);

/** Create a random partition table.
 * The use of this function is only reasonable for debugging.
 * \param[in]        comm    The communicator to use.
 * \param[in]        num_trees The number of global trees in the partition.
 * \param[in]        shared  If true than there will be shared trees in the generated partition table.
 * \param[in]        seed    A seed to be used for the random number generator. If zero, a random
 *                           seed is chosen. \a seed has to be the same on each process.
 * \return                   A valid partition table for a cmesh with \a num_trees trees
 *                           and communicator \a comm, where each processor gets a random number
 *                           of trees. The number of trees per processor is roughly uniformly distributed.
 * \note This function is MPI collective. 
 */
t8_shmem_array_t
t8_cmesh_offset_random (sc_MPI_Comm comm, t8_gloidx_t num_trees, int shared, unsigned seed);

/** Create a repartition array, where each process sends half of its
 * trees to the next process. The last process does not send any trees.
 * \param [in]      cmesh   A cmesh that is committed and partitioned.
 * \param [in]      comm    A valid MPI communicator for cmesh.
 * \return                  A shared memory offset array storing the new offsets.
 *                          Its associated communicator is \a comm.
 * \note This function is MPI collective.
 */
t8_shmem_array_t
t8_cmesh_offset_half (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Create a repartition array, where each process sends a percentage of its
 * trees to the next process. The last process does not send any trees.
 * \param [in]      cmesh   A cmesh that is committed and partitioned.
 * \param [in]      comm    A valid MPI communicator for cmesh.
 * \param [in]      percent The percentage of trees that this process should send
 *                          to the next one. Must satisfy 0 <= \a percent <= 100 and
 *                          be the same on each process.
 * \return                  A shared memory offset array storing the new offsets.
 *                          Its associated communicator is \a comm.
 * \note This function is MPI collective.
 */
t8_shmem_array_t
t8_cmesh_offset_percent (t8_cmesh_t cmesh, sc_MPI_Comm comm, int percent);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PARTITION_H */
