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
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_PARTITION_H
#define T8_CMESH_PARTITION_H

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"

T8_EXTERN_C_BEGIN ();

/** Given a cmesh which is to be partitioned, execute the partition task.
 *  This includes partitioning by uiniform level and partitioning from a second cmesh
 *  \param [in,out]  cmesh   The cmesh to be partitioned
 *  This function is usually called within \ref t8_cmesh_commit */
void                t8_cmesh_partition (t8_cmesh_t cmesh);

/** Create a valid partition table that concentrates all trees at a given
 *  process.
 * \param[in]        proc    The processor that should get all trees.
 * \param[in]        comm    The communicator to use.
 * \param[in]        num_trees The number of global trees in the partition.
 * \return                   A valid partition table for a mesh with \a num_trees trees
 *                           and communicator \a comm, where each tree is on process \a proc.
 */
t8_gloidx_t        *t8_cmesh_offset_concentrate (int proc, sc_MPI_Comm comm,
                                                 t8_gloidx_t num_trees);

/** Create a random partition table.
 * The use of this function is only reasonable for debugging.
 * \param[in]        comm    The communicator to use.
 * \param[in]        num_trees The number of global trees in the partition.
 * \param[in]        shared  If true than there will be shared trees in the generated partition table.
 * \return                   A valid partition table for a mesh with \a num_trees trees
 *                           and communicator \a comm, where each processor gets a random number
 *                           of trees. The number of trees per processor is roughly uniformly distributed.
 */
t8_gloidx_t        *t8_cmesh_offset_random (sc_MPI_Comm comm,
                                            t8_gloidx_t num_trees,
                                            int shared);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PARTITION_H */
