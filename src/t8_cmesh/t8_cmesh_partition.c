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

/** \file t8_cmesh_partition.c
 *
 * TODO: document this file
 */

#include <t8_shmem.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_partition.h"

/* Given the trees_per_proc array compute the first local tree */
static void
t8_cmesh_scan_trees (t8_cmesh_t cmesh)
{
  int                 iproc;
  int                 last_tree_shared;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->tree_per_proc != NULL);

  cmesh->first_tree = 0;
  for (iproc = 0; iproc < cmesh->mpirank; iproc++) {
    last_tree_shared = cmesh->tree_per_proc[iproc] < 0;
    cmesh->first_tree += last_tree_shared ? -cmesh->tree_per_proc[iproc] - 1
      : cmesh->tree_per_proc[iproc];
  }
  /* TODO :An alternative is to use MPI_Scan with a self-defined operation.
   *       This could be faster since MPI_Scan is O(log(P)) instead of O(P)
   *       but it introduces a new synchronization point. */

}

static void
t8_cmesh_partition_given (t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from,
                          t8_locidx_t * trees_per_proc)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (cmesh_from != NULL);
  T8_ASSERT (cmesh_from->committed);
}

void
t8_cmesh_partition (t8_cmesh_t cmesh)
{
  t8_cmesh_t          cmesh_from;
  t8_gloidx_t         last_tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->set_from->committed);

  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  /**********************************************/
  /*      Compute local number of trees         */
  /*         and trees per proc array           */
  /**********************************************/
  if (cmesh->set_level >= 0) {
    /* Compute first and last tree index */
    T8_ASSERT (cmesh->tree_per_proc == NULL);
    t8_cmesh_uniform_bounds (cmesh_from, cmesh->set_level, &cmesh->first_tree,
                             NULL, &last_tree, NULL,
                             &cmesh->last_tree_shared);
    cmesh->num_local_trees = last_tree - cmesh->first_tree;
    /* To compute the tree_per_proc correctly we have to invert the sign on the
     * local tree count if last tree is shared */
    if (cmesh->last_tree_shared) {
      cmesh->num_local_trees = -cmesh->num_local_trees;
    }
    /* allocate and fill tree_per_proc array with number of trees per process */
    /* We have to be careful with shared memory in the case where cmesh comm is
     * duplicated, since we have to use the same communicator for allocating and freeing memory.
     * Thus this function must only be called after cmesh communicator was duplicated. */
    T8_ASSERT (cmesh->mpisize > 0);
    cmesh->tree_per_proc = SC_SHMEM_ALLOC (t8_locidx_t, cmesh->mpisize,
                                           cmesh->mpicomm);
    sc_shmem_allgather (&cmesh->num_local_trees, 1, T8_MPI_LOCIDX,
                        cmesh->tree_per_proc, 1, T8_MPI_LOCIDX,
                        cmesh->mpicomm);
    /* tree_to_proc was computed, reinvert the sign */
    if (cmesh->last_tree_shared) {
      T8_ASSERT (cmesh->num_local_trees < 0);
      cmesh->num_local_trees = -cmesh->num_local_trees;
    }
  }
  else {
    T8_ASSERT (cmesh->tree_per_proc != NULL);
    cmesh->last_tree_shared = cmesh->tree_per_proc[cmesh->mpirank] < 0;
    cmesh->num_local_trees = cmesh->last_tree_shared ?
      -cmesh->tree_per_proc[cmesh->mpirank] :
      cmesh->tree_per_proc[cmesh->mpirank];
    t8_cmesh_scan_trees (cmesh);
  }
  /***************************************************/
  /*        Done with local num and trees_per_proc   */
  /***************************************************/
  t8_cmesh_partition_given (cmesh, cmesh->set_from, cmesh->tree_per_proc);
}
