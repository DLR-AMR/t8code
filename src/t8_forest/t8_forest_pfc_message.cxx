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
 * \file t8_forest_pfc_message.cxx
 * This file implements the functions for \ref t8_forest_pfc_message.hxx
*/

#include <t8_forest/t8_forest_pfc_message.hxx>
#include <t8_forest/t8_forest_types.h>

void
t8_forest_pfc_message::pack (void *buf, int buf_size, int *position)
{
  /* pack: itree */
  sc_MPI_Pack (&itree, 1, T8_MPI_GLOIDX, buf, buf_size, position, comm);

  /* pack: eclass */
  int eclass_int = (int) eclass;
  sc_MPI_Pack (&eclass_int, 1, sc_MPI_INT, buf, buf_size, position, comm);

  /* pack: parent */
  t8_element_MPI_Pack (scheme, eclass, &parent, 1, buf, buf_size, position, comm);

  /* pack: num_siblings */
  sc_MPI_Pack (&num_siblings, 1, sc_MPI_INT, buf, buf_size, position, comm);
}

void
t8_forest_pfc_message::unpack (void *buf, int buf_size, int *position)
{
  /* unpack: itree */
  sc_MPI_Unpack (buf, buf_size, position, &itree, 1, T8_MPI_GLOIDX, comm);

  /* unpack: eclass */
  int eclass_int;
  sc_MPI_Unpack (buf, buf_size, position, &eclass_int, 1, sc_MPI_INT, comm);
  eclass = (t8_eclass_t) eclass_int;

  /* unpack: parent */
  t8_element_new (scheme, eclass, 1, &parent);
  t8_element_MPI_Unpack (scheme, eclass, buf, buf_size, position, &parent, 1, comm);

  /* unpack: num_siblings */
  sc_MPI_Unpack (buf, buf_size, position, &num_siblings, 1, sc_MPI_INT, comm);
}

int
t8_forest_pfc_message::pack_size ()
{
  /* initialize sum*/
  int message_size = 0;
  int datasize;

  /* add size: itree */
  sc_MPI_Pack_size (1, T8_MPI_GLOIDX, comm, &datasize);
  message_size += datasize;

  /* add size: eclass */
  sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  message_size += datasize;

  /* add size: parent */
  if (parent != nullptr) {
    t8_element_MPI_Pack_size (scheme, eclass, 1, comm, &datasize);
    message_size += datasize;
  }

  /* add size: num_siblings */
  sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  message_size += datasize;
  return message_size;
}

void
t8_forest_pfc_message::mpi_Isend (const t8_forest_t forest, sc_MPI_Request &request)
{
  /* Allocate buffer */
  int buffer_size = pack_size ();
  char *send_buffer = T8_ALLOC (char, buffer_size);

  /* Pack message to buffer */
  int position = 0;
  pack (send_buffer, buffer_size, &position);

  /* Send buffer */
  const int mpiret = sc_MPI_Isend (send_buffer, position, sc_MPI_PACKED, iproc, T8_MPI_PFC_TAG,
                                   t8_forest_get_mpicomm (forest), &request);
  SC_CHECK_MPI (mpiret);
  T8_FREE (send_buffer);
}

void
t8_forest_pfc_message::mpi_Recv (char *&recv_buf, int &buf_size)
{
  /** Get needed size of message via MPI probe and allocate buffer. */
  sc_MPI_Status status;
  sc_MPI_Probe (iproc, T8_MPI_PFC_TAG, comm, &status);
  sc_MPI_Get_count (&status, sc_MPI_PACKED, &buf_size);
  recv_buf = T8_ALLOC (char, buf_size);

  /* Actually receive buffer. */
  const int mpiret = sc_MPI_Recv (recv_buf, buf_size, sc_MPI_PACKED, iproc, T8_MPI_PFC_TAG, comm, sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);
}

/**
   * Fill this instance of class t8_forest_pfc_message, i.e., setting the member viarables based on the given forest.
   *
   * \param[in] forest the forest
  */
void
t8_forest_pfc_message::fill (t8_forest_t forest)
{
  // Set mpi rank and partition element offsets.
  const t8_procidx_t rank = t8_forest_get_mpirank (forest);
  const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

  // Determine the global element ID of the current process ("rank") that is closest to the elements of process iproc:
  // If iproc < rank, the first element is closest to iproc; otherwise, the last one is.
  const t8_gloidx_t closest_to_rank_gid = (iproc <= rank) ? partition[rank] : partition[rank + 1] - 1;

  // Use helper function to get various element and tree indices.
  t8_tree_t tree;
  t8_locidx_t index_in_tree;
  t8_element_t *element_closest_to_receiver;
  t8_forest_pfc_helper_index_in_tree_from_globalid (forest, closest_to_rank_gid, itree, tree, index_in_tree,
                                                    element_closest_to_receiver);
  // Set scheme and eclass.
  scheme = t8_forest_get_scheme (forest);
  eclass = t8_forest_get_eclass (forest, t8_forest_get_local_id (forest, itree));

  // Allocate memory for the parent element.
  t8_element_new (scheme, eclass, 1, &parent);

  // If we are already the root element, we cannot be part of a split family, so we send any(the root) element and no num_siblings.
  if (scheme->element_get_level (eclass, element_closest_to_receiver) == 0) {
    // Copy the root element to the allocated parent element.
    scheme->element_copy (eclass, element_closest_to_receiver, parent);
    // Set the number of siblings to zero, since there are none on the root level.
    num_siblings = 0;
  }
  else {
    // Compute parent.
    scheme->element_get_parent (eclass, element_closest_to_receiver, parent);

    // Distinguish send "direction"
    if (iproc > rank) {
      // Sending will be towards higher-rank process, so count siblings in decreasing index direction.
      const t8_locidx_t min_id = t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, -1);
      num_siblings = index_in_tree - min_id + 1;
    }
    else {
      T8_ASSERT (iproc != rank);
      // Sending will be towards lower-rank process, so count siblings in increasing index direction.
      const t8_locidx_t max_id = t8_forest_pfc_extreme_local_sibling (scheme, tree, index_in_tree, +1);
      num_siblings = max_id - index_in_tree + 1;
    }
  }
}
