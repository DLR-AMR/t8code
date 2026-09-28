/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_ghost_definition_base.cxx
 * Implementation details for t8_forest_ghost_definition_base.hxx
 */

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_base.hxx>

int
t8_forest_ghost_definition::do_ghost (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  if (ghost_get_type () == T8_GHOST_NONE) {
    t8_errorf ("WARNING: Trying to construct ghosts with ghost_type NONE. "
               "Ghost layer is not constructed.\n");
    return 0;
  }

  const int memory_flag = communicate_ownerships (forest);

  /* Processes without local elements also get an (empty) ghost structure, so that
   * all processes agree on whether a ghost layer exists. */
  t8_forest_ghost_init (&forest->ghosts, ghost_type);

  if (t8_forest_get_local_num_leaf_elements (forest) > 0) {
    fill_remote_ghosts (forest);

    communicate_ghost_elements (forest);
  }
  clean_up (forest, memory_flag);

  return 1;
}

int
t8_forest_ghost_definition::communicate_ownerships (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  int memory_flag = 0;

  if (forest->element_offsets == nullptr) {
    /* create element offset array if not done already */
    memory_flag = memory_flag | CREATE_ELEMENT_ARRAY;
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == nullptr) {
    /* Create tree offset array if not done already */
    memory_flag = memory_flag | CREATE_TREE_ARRAY;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == nullptr) {
    /* Create global first desc array if not done already */
    memory_flag = memory_flag | CREATE_GFIRST_DESC_ARRAY;
    t8_forest_partition_create_first_desc (forest);
  }
  return memory_flag;
}

void
t8_forest_ghost_definition::communicate_ghost_elements (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  t8_forest_ghost_t const ghost = forest->ghosts;
  sc_MPI_Request *requests;

  /* Start sending the remote elements */
  t8_ghost_mpi_send_info_t *const send_info = t8_forest_ghost_send_start (forest, ghost, &requests);

  /* Receive the ghost elements from the remote processes */
  t8_forest_ghost_receive (forest, ghost);

  /* End sending the remote elements */
  t8_forest_ghost_send_end (forest, ghost, send_info, requests);
}

void
t8_forest_ghost_definition::clean_up (t8_forest_t forest, const int memory_flag)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  if (memory_flag & CREATE_ELEMENT_ARRAY) {
    /* Free the offset memory, if allocated */
    t8_shmem_array_destroy (&forest->element_offsets);
  }
  if (memory_flag & CREATE_TREE_ARRAY) {
    /* Free the offset memory, if allocated */
    t8_shmem_array_destroy (&forest->tree_offsets);
  }
  if (memory_flag & CREATE_GFIRST_DESC_ARRAY) {
    /* Free the offset memory, if allocated */
    t8_shmem_array_destroy (&forest->global_first_desc);
  }
}
