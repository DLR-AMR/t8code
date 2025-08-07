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

/** \file t8_forest_ghost_definition_helpers.hxx
 * Implements helper functions for different ghost definitions.
 */

#ifndef T8_FOREST_GHOST_DEFINITION_HELPERS_HXX
#define T8_FOREST_GHOST_DEFINITION_HELPERS_HXX

#include <t8_forest/t8_forest_general.h>

/* The information stored for the ghost trees */
typedef struct
{
  t8_gloidx_t global_id;       /* global id of the tree */
  t8_locidx_t element_offset;  /* The count of all ghost elements in all smaller ghost trees */
  t8_element_array_t elements; /* The ghost elements of that tree */
  t8_eclass_t eclass;          /* The trees element class */
} t8_ghost_tree_t;

/* The data structure stored in the global_tree_to_ghost_tree hash table. */
typedef struct
{
  t8_gloidx_t global_id; /* global tree id */
  size_t index;          /* the index of that global tree in the ghost_trees array. */
} t8_ghost_gtree_hash_t;

/* The data structure stored in the process_offsets array. */
typedef struct
{
  int mpirank;              /* rank of the process */
  t8_locidx_t ghost_offset; /* The number of ghost elements for all previous ranks */
  size_t tree_index;        /* index of first ghost tree of this process in ghost_trees */
  size_t first_element;     /* the index of the first element in the elements array of the ghost tree. */
} t8_ghost_process_hash_t;

/* The information for a remote process, what data we have to send to them.
 */
typedef struct
{
  int recv_rank;           /* The rank to which we send. */
  size_t num_bytes;        /* The number of bytes that we send. */
  sc_MPI_Request *request; /* Communication request, not owned by this struct. */
  char *buffer;            /* The send buffer. */
} t8_ghost_mpi_send_info_t;

/**
 * Initializes the forest ghost structure and allocates memory.
 * \param [in, out] pghost  The ghost structure
 * \param [in] ghost_type   The type of ghost to use.
 */
void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type);

/**
  * Add a new element to the remote hash table (if not already in it).
  * Must be called for elements in linear order
  * element_index is the tree local index of this element
  * \param [in, out] forest         The forest
  * \param [in, out] ghost          The ghosts of the forest
  * \param [in] remote_rank         The other rank
  * \param [in] ltreeid             The local tree id of \a elem
  * \param [in] elem                The element to add
  * \param [in] element_index       The index of the \a elem
  */
void
t8_ghost_add_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int remote_rank, t8_locidx_t ltreeid,
                     const t8_element_t *elem, t8_locidx_t element_index);

/* Begin sending the ghost elements from the remote ranks
 * using non-blocking communication.
 * Afterwards,
 *  t8_forest_ghost_send_end
 * must be called to end the communication.
 * Returns an array of mpi_send_info_t, one for each remote rank.
 */
t8_ghost_mpi_send_info_t *
t8_forest_ghost_send_start (t8_forest_t forest, t8_forest_ghost_t ghost, sc_MPI_Request **requests);

void
t8_forest_ghost_send_end (t8_forest_t forest, t8_forest_ghost_t ghost, t8_ghost_mpi_send_info_t *send_info,
                          sc_MPI_Request *requests);

/* Probe for all incoming messages from the remote ranks and receive them.
 * We receive the message in the order in which they arrive. To achieve this,
 * we have to use polling. */
void
t8_forest_ghost_receive (t8_forest_t forest, t8_forest_ghost_t ghost);

#endif /* T8_FOREST_GHOST_DEFINITION_HELPERS_HXX */
