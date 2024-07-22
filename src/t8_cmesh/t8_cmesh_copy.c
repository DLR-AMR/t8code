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

/** \file t8_cmesh_copy.c
 *
 * TODO: document this file
 */

#include <t8_data/t8_shmem.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_partition.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_copy.h"

void
t8_cmesh_copy (t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, sc_MPI_Comm comm)
{
  size_t num_parts, iz;
  t8_locidx_t first_tree, num_trees, first_ghost, num_ghosts;

  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (t8_cmesh_is_committed (cmesh_from));

  /* Copy all variables */
  cmesh->dimension = cmesh_from->dimension;
  cmesh->face_knowledge = cmesh_from->face_knowledge;
  cmesh->first_tree = cmesh_from->first_tree;
  cmesh->first_tree_shared = cmesh_from->first_tree_shared;
  cmesh->mpirank = cmesh_from->mpirank;
  cmesh->mpisize = cmesh_from->mpisize;
  cmesh->num_ghosts = cmesh_from->num_ghosts;
  cmesh->num_local_trees = cmesh_from->num_local_trees;
  cmesh->num_trees = cmesh_from->num_trees;
  cmesh->set_partition = cmesh_from->set_partition;
  cmesh->set_partition_level = cmesh_from->set_partition_level;
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  /* Copy the tree_offsets */
  if (cmesh_from->tree_offsets != NULL) {
    T8_ASSERT (cmesh->tree_offsets == NULL);
    cmesh->tree_offsets = t8_cmesh_alloc_offsets (cmesh->mpisize, comm);
    t8_shmem_array_copy (cmesh->tree_offsets, cmesh_from->tree_offsets);
  }
  /* Copy the numbers of trees */
  memcpy (cmesh->num_trees_per_eclass, cmesh_from->num_trees_per_eclass, T8_ECLASS_COUNT * sizeof (t8_gloidx_t));
  /* Copy the numbers of local trees */
  memcpy (cmesh->num_local_trees_per_eclass, cmesh_from->num_local_trees_per_eclass,
          T8_ECLASS_COUNT * sizeof (t8_locidx_t));

  /* Copy the tree info */
  cmesh->trees = NULL;
  if (cmesh_from->trees) {
    num_parts = t8_cmesh_trees_get_numproc (cmesh_from->trees);
    t8_cmesh_trees_init (&cmesh->trees, num_parts, cmesh_from->num_local_trees, cmesh_from->num_ghosts);
    t8_cmesh_trees_copy_toproc (cmesh->trees, cmesh_from->trees, cmesh_from->num_local_trees, cmesh_from->num_ghosts);
    for (iz = 0; iz < num_parts; iz++) {
      t8_cmesh_trees_get_part_data (cmesh_from->trees, iz, &first_tree, &num_trees, &first_ghost, &num_ghosts);
      t8_cmesh_trees_start_part (cmesh->trees, iz, first_tree, num_trees, first_ghost, num_ghosts, 0);
      t8_cmesh_trees_copy_part (cmesh->trees, iz, cmesh_from->trees, iz);
    }
  }
}
