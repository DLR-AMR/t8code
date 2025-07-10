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

void
t8_forest_ghost_definition::communicate_ownerships (t8_forest_t forest)
{
  if (forest->element_offsets == NULL) {
    /* create element offset array if not done already */
    memory_flag = memory_flag | CREATE_ELEMENT_ARRAY;
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == NULL) {
    /* Create tree offset array if not done already */
    memory_flag = memory_flag | CREATE_TREE_ARRAY;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    /* Create global first desc array if not done already */
    memory_flag = memory_flag | CREATE_GFIRST_DESC_ARRAY;
    t8_forest_partition_create_first_desc (forest);
  }
}

void
t8_forest_ghost_definition::communicate_ghost_elements (t8_forest_t forest)
{
  t8_forest_ghost_t ghost = forest->ghosts;
  t8_ghost_mpi_send_info_t *send_info;
  sc_MPI_Request *requests;

  /* Start sending the remote elements */
  send_info = t8_forest_ghost_send_start (forest, ghost, &requests);

  /* Receive the ghost elements from the remote processes */
  t8_forest_ghost_receive (forest, ghost);

  /* End sending the remote elements */
  t8_forest_ghost_send_end (forest, ghost, send_info, requests);
}

void
t8_forest_ghost_definition::clean_up (t8_forest_t forest)
{
  if (memory_flag & CREATE_GFIRST_DESC_ARRAY) {
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
