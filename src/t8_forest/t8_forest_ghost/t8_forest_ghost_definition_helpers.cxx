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
 * Implementations for t8_forest_ghost_definition_helpers.hxx
 */

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>

/* The hash function for the global tree hash.
 * As hash value we just return the global tree id. */
static unsigned
t8_ghost_gtree_hash_function (const void *ghost_gtree_hash, [[maybe_unused]] const void *data)
{
  const t8_ghost_gtree_hash_t *object = (const t8_ghost_gtree_hash_t *) ghost_gtree_hash;

  return (unsigned) object->global_id;
}

/* The equal function for two global tree hash objects.
 * Two t8_ghost_gtree_hash_t are considered equal if the global tree ids are the same.
 */
static int
t8_ghost_gtree_equal_function (const void *ghost_gtreea, const void *ghost_gtreeb, [[maybe_unused]] const void *user)
{
  const t8_ghost_gtree_hash_t *objecta = (const t8_ghost_gtree_hash_t *) ghost_gtreea;
  const t8_ghost_gtree_hash_t *objectb = (const t8_ghost_gtree_hash_t *) ghost_gtreeb;

  /* return true if and only if the global_ids are the same */
  return objecta->global_id == objectb->global_id;
}

/* The hash value for an entry of the process_offsets hash is the processes mpirank. */
static unsigned
t8_ghost_process_hash_function (const void *process_data, [[maybe_unused]] const void *user_data)
{
  const t8_ghost_process_hash_t *process = (const t8_ghost_process_hash_t *) process_data;

  return process->mpirank;
}

/* The equal function for the process_offsets array.
 * Two entries are the same if their mpiranks are equal. */
static int
t8_ghost_process_equal_function (const void *process_dataa, const void *process_datab,
                                 [[maybe_unused]] const void *user)
{
  const t8_ghost_process_hash_t *processa = (const t8_ghost_process_hash_t *) process_dataa;
  const t8_ghost_process_hash_t *processb = (const t8_ghost_process_hash_t *) process_datab;

  return processa->mpirank == processb->mpirank;
}

/* The hash function for the remote_ghosts hash table.
 * The hash value for an mpirank is just the rank */
static unsigned
t8_ghost_remote_hash_function (const void *remote_data, [[maybe_unused]] const void *user_data)
{
  const t8_ghost_remote_t *remote = (const t8_ghost_remote_t *) remote_data;

  return remote->remote_rank;
}

/* The equal function for the remote hash table.
 * Two entries are the same if they have the same rank. */
static int
t8_ghost_remote_equal_function (const void *remote_dataa, const void *remote_datab, [[maybe_unused]] const void *user)
{
  const t8_ghost_remote_t *remotea = (const t8_ghost_remote_t *) remote_dataa;
  const t8_ghost_remote_t *remoteb = (const t8_ghost_remote_t *) remote_datab;

  return remotea->remote_rank == remoteb->remote_rank;
}

/* Initialize a t8_ghost_remote_tree_t */
static void
t8_ghost_init_remote_tree (t8_forest_t forest, t8_gloidx_t gtreeid, int remote_rank, t8_eclass_t tree_class,
                           t8_ghost_remote_tree_t *remote_tree)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_locidx_t local_treeid;

  T8_ASSERT (remote_tree != NULL);

  local_treeid = gtreeid - t8_forest_get_first_local_tree_id (forest);
  /* Set the entries of the new remote tree */
  remote_tree->global_id = gtreeid;
  remote_tree->mpirank = remote_rank;
  remote_tree->eclass = t8_forest_get_eclass (forest, local_treeid);
  /* Initialize the array to store the element */
  t8_element_array_init (&remote_tree->elements, scheme, tree_class);
  /* Initialize the array to store the element indices. */
  sc_array_init (&remote_tree->element_indices, sizeof (t8_locidx_t));
}

void
t8_ghost_add_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int remote_rank, t8_locidx_t ltreeid,
                     const t8_element_t *elem, t8_locidx_t element_index)
{
  t8_ghost_remote_t remote_entry_lookup, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_element_t *elem_copy;
  sc_array_t *remote_array;
  size_t index, element_count;
  int *remote_process_entry;
  int level, copy_level = 0;

  /* Get the tree's element class and the scheme */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_gloidx_t gtreeid = t8_forest_get_first_local_tree_id (forest) + ltreeid;

  /* Check whether the remote_rank is already present in the remote ghosts
   * array. */
  remote_entry_lookup.remote_rank = remote_rank;
  /* clang-format off */
  remote_entry = (t8_ghost_remote_t *) sc_hash_array_insert_unique (ghost->remote_ghosts, (void *) &remote_entry_lookup,
                                                                    &index);
  /* clang-format on */
  if (remote_entry != NULL) {
    /* The remote rank was not in the array and was inserted now */
    remote_entry->remote_rank = remote_rank;
    remote_entry->num_elements = 0;
    /* Initialize the tree array of the new entry */
    sc_array_init_size (&remote_entry->remote_trees, sizeof (t8_ghost_remote_tree_t), 1);
    /* Get a pointer to the new entry */
    remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees, 0);
    /* initialize the remote_tree */
    t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, tree_class, remote_tree);
    /* Since the rank is a new remote rank, we also add it to the remote ranks array */
    remote_process_entry = (int *) sc_array_push (ghost->remote_processes);
    *remote_process_entry = remote_rank;
  }
  else {
    /* The remote rank already is contained in the remotes array at position index. */
    remote_array = &ghost->remote_ghosts->a;
    remote_entry = (t8_ghost_remote_t *) sc_array_index (remote_array, index);
    T8_ASSERT (remote_entry->remote_rank == remote_rank);
    /* Check whether the tree has already an entry for this process.
     * Since we only add in local tree order the current tree is either
     * the last entry or does not have an entry yet. */
    remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees,
                                                             remote_entry->remote_trees.elem_count - 1);
    if (remote_tree->global_id != gtreeid) {
      /* The tree does not exist in the array. We thus need to add it and
       * initialize it. */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_push (&remote_entry->remote_trees);
      t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, tree_class, remote_tree);
    }
  }
  /* remote_tree now points to a valid entry for the tree.
   * We can add a copy of the element to the elements array
   * if it does not exist already. If it exists it is the last entry in the array. */
#if T8_ENABLE_DEBUG
  {
    /* debugging assertion that the element is really not contained already */
    int ielem;
    int elem_count = t8_element_array_get_count (&remote_tree->elements);
    for (ielem = 0; ielem < elem_count - 1; ielem++) {
      const t8_element_t *test_el = t8_element_array_index_int (&remote_tree->elements, ielem);
      SC_CHECK_ABORTF (!scheme->element_is_equal (tree_class, test_el, elem),
                       "Local element %i already in remote ghosts at pos %i\n", element_index, ielem);
    }
  }
#endif
  elem_copy = NULL;
  level = scheme->element_get_level (tree_class, elem);
  element_count = t8_element_array_get_count (&remote_tree->elements);
  if (element_count > 0) {
    elem_copy = t8_element_array_index_locidx_mutable (&remote_tree->elements, element_count - 1);
    copy_level = scheme->element_get_level (tree_class, elem_copy);
  }
  /* Check if the element was not contained in the array.
   * If so, we add a copy of elem to the array.
   * Otherwise, we do nothing. */
  if (elem_copy == NULL || level != copy_level
      || scheme->element_get_linear_id (tree_class, elem_copy, copy_level)
           != scheme->element_get_linear_id (tree_class, elem, level)) {
    /* Add the element */
    elem_copy = t8_element_array_push (&remote_tree->elements);
    scheme->element_copy (tree_class, elem, elem_copy);
    /* Add the index of the element */
    *(t8_locidx_t *) sc_array_push (&remote_tree->element_indices) = element_index;
    remote_entry->num_elements++;
  }
}

void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type)
{
  t8_forest_ghost_t ghost;

  T8_ASSERT (ghost_type != T8_GHOST_NONE);

  /* Allocate memory for ghost */
  ghost = *pghost = T8_ALLOC_ZERO (t8_forest_ghost_struct_t, 1);
  /* initialize the reference counter */
  t8_refcount_init (&ghost->rc);
  /* Set the ghost type */
  ghost->ghost_type = ghost_type;

  /* Allocate the trees array */
  ghost->ghost_trees = sc_array_new (sizeof (t8_ghost_tree_t));

  /* initialize the global_tree_to_ghost_tree hash table */
  ghost->glo_tree_mempool = sc_mempool_new (sizeof (t8_ghost_gtree_hash_t));
  ghost->global_tree_to_ghost_tree
    = sc_hash_new (t8_ghost_gtree_hash_function, t8_ghost_gtree_equal_function, NULL, NULL);

  /* initialize the process_offset hash table */
  ghost->proc_offset_mempool = sc_mempool_new (sizeof (t8_ghost_process_hash_t));
  ghost->process_offsets = sc_hash_new (t8_ghost_process_hash_function, t8_ghost_process_equal_function, NULL, NULL);
  /* initialize the remote ghosts hash table */
  ghost->remote_ghosts = sc_hash_array_new (sizeof (t8_ghost_remote_t), t8_ghost_remote_hash_function,
                                            t8_ghost_remote_equal_function, NULL);
  /* initialize the remote processes array */
  ghost->remote_processes = sc_array_new (sizeof (int));
}

/* Begin sending the ghost elements from the remote ranks
 * using non-blocking communication.
 * Afterwards,
 *  t8_forest_ghost_send_end
 * must be called to end the communication.
 * Returns an array of mpi_send_info_t, one for each remote rank.
 */
t8_ghost_mpi_send_info_t *
t8_forest_ghost_send_start (t8_forest_t forest, t8_forest_ghost_t ghost, sc_MPI_Request **requests)
{
  int proc_index, remote_rank;
  int num_remotes;
  size_t remote_index;
  t8_ghost_remote_t *remote_entry;
  sc_array_t *remote_trees;
  t8_ghost_remote_tree_t *remote_tree = NULL;
  t8_ghost_mpi_send_info_t *send_info, *current_send_info;
  char *current_buffer;
  size_t bytes_written, element_bytes, element_count, element_size;
#if T8_ENABLE_DEBUG
  size_t acc_el_count = 0;
#endif
  int mpiret;

  /* Allocate a send_buffer for each remote rank */
  num_remotes = ghost->remote_processes->elem_count;
  send_info = T8_ALLOC (t8_ghost_mpi_send_info_t, num_remotes);
  *requests = T8_ALLOC (sc_MPI_Request, num_remotes);

  /* Loop over all remote processes */
  for (proc_index = 0; proc_index < (int) ghost->remote_processes->elem_count; proc_index++) {
    current_send_info = send_info + proc_index;
    /* Get the rank of the current remote process. */
    remote_rank = *(int *) sc_array_index_int (ghost->remote_processes, proc_index);
    t8_debugf ("Filling send buffer for process %i\n", remote_rank);
    /* initialize the send_info for the current rank */
    current_send_info->recv_rank = remote_rank;
    current_send_info->num_bytes = 0;
    current_send_info->request = *requests + proc_index;
    /* Lookup the ghost elements for the first tree of this remote */
    remote_entry = t8_forest_ghost_get_remote (forest, remote_rank);
    T8_ASSERT (remote_entry->remote_rank == remote_rank);
    /* Loop over all trees of the remote rank and count the bytes */
    /* At first we store the number of remote trees in the buffer */
    current_send_info->num_bytes += sizeof (size_t);
    /* add padding before the eclass */
    current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
    /* TODO: put this in a function */
    /* TODO: Use remote_entry to count the number of bytes while inserting
     *        the remote ghosts. */
    remote_trees = &remote_entry->remote_trees;
    for (remote_index = 0; remote_index < remote_trees->elem_count; remote_index++) {
      /* Get the next remote tree. */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (remote_trees, remote_index);
      /* We will store the global tree id, the element class and the list
       * of elements in the send_buffer. */
      current_send_info->num_bytes += sizeof (t8_gloidx_t);
      /* add padding before the eclass */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += sizeof (t8_eclass_t);
      /* add padding before the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_count = t8_element_array_get_count (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* We will store the number of elements */
      current_send_info->num_bytes += sizeof (size_t);
      /* add padding before the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += element_bytes;
      /* add padding after the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
    }

    /* We now now the number of bytes for our send_buffer and thus allocate it. */
    current_send_info->buffer = T8_ALLOC_ZERO (char, current_send_info->num_bytes);

    /* We iterate through the tree again and store the tree info and the elements into the send_buffer. */
    current_buffer = current_send_info->buffer;
    bytes_written = 0;
    /* Start with the number of remote trees in the buffer */
    memcpy (current_buffer + bytes_written, &remote_trees->elem_count, sizeof (size_t));
    bytes_written += sizeof (size_t);
    bytes_written += T8_ADD_PADDING (bytes_written);
#if T8_ENABLE_DEBUG
    acc_el_count = 0;
#endif
    for (remote_index = 0; remote_index < remote_trees->elem_count; remote_index++) {
      /* Get a pointer to the tree */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (remote_trees, remote_index);
      T8_ASSERT (remote_tree->mpirank == remote_rank);

      /* Copy the global tree id */
      memcpy (current_buffer + bytes_written, &remote_tree->global_id, sizeof (t8_gloidx_t));
      bytes_written += sizeof (t8_gloidx_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Copy the trees element class */
      memcpy (current_buffer + bytes_written, &remote_tree->eclass, sizeof (t8_eclass_t));
      bytes_written += sizeof (t8_eclass_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Store the number of elements in the buffer */
      element_count = t8_element_array_get_count (&remote_tree->elements);
      memcpy (current_buffer + bytes_written, &element_count, sizeof (size_t));
      bytes_written += sizeof (size_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* Copy the elements into the send buffer */
      memcpy (current_buffer + bytes_written, t8_element_array_get_data (&remote_tree->elements), element_bytes);
      bytes_written += element_bytes;
      /* add padding after the elements */
      bytes_written += T8_ADD_PADDING (bytes_written);

      /* Add to the counter of remote elements. */
      ghost->num_remote_elements += element_count;
#if T8_ENABLE_DEBUG
      acc_el_count += element_count;
#endif
    } /* End tree loop */

    T8_ASSERT (bytes_written == current_send_info->num_bytes);
    /* We can now post the MPI_Isend for the remote process */
    mpiret = sc_MPI_Isend (current_buffer, bytes_written, sc_MPI_BYTE, remote_rank, T8_MPI_GHOST_FOREST,
                           forest->mpicomm, *requests + proc_index);
    SC_CHECK_MPI (mpiret);
  } /* end process loop */
  return send_info;
}

void
t8_forest_ghost_send_end ([[maybe_unused]] t8_forest_t forest, t8_forest_ghost_t ghost,
                          t8_ghost_mpi_send_info_t *send_info, sc_MPI_Request *requests)
{
  int num_remotes;
  int proc_pos, mpiret;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (ghost != NULL);

  /* Get the number of remote processes */
  num_remotes = ghost->remote_processes->elem_count;

  /* We wait for all communication to end. */
  mpiret = sc_MPI_Waitall (num_remotes, requests, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* Clean-up */
  for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
    T8_FREE (send_info[proc_pos].buffer);
  }
  T8_FREE (send_info);
  T8_FREE (requests);
}

/* Receive a single message from a remote process, after the message was successfully probed.
 * Returns the allocated receive buffer and the number of bytes received */
static char *
t8_forest_ghost_receive_message (int recv_rank, sc_MPI_Comm comm, sc_MPI_Status status, int *recv_bytes)
{
  char *recv_buffer;
  int mpiret;

  T8_ASSERT (recv_rank == status.MPI_SOURCE);
  T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);

  /* Get the number of bytes in the message */
  mpiret = sc_MPI_Get_count (&status, sc_MPI_BYTE, recv_bytes);

  /* Allocate receive buffer */
  recv_buffer = T8_ALLOC_ZERO (char, *recv_bytes);
  /* receive the message */
  mpiret
    = sc_MPI_Recv (recv_buffer, *recv_bytes, sc_MPI_BYTE, recv_rank, T8_MPI_GHOST_FOREST, comm, sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);

  return recv_buffer;
}

/* Parse a message from a remote process and correctly include the received
 * elements in the ghost structure.
 * The message looks like:
 * num_trees | pad | treeid 0 | pad | eclass 0 | pad | num_elems 0 | pad | elements | pad | treeid 1 | ...
 *  size_t   |     |t8_gloidx |     |t8_eclass |     | size_t      |     | t8_element_t |
 *
 * pad is paddind, see T8_ADD_PADDING
 *
 * current_element_offset is updated in each step to store the element offset
 * of the next ghost tree to be inserted.
 * When called with the first message, current_element_offset must be set to 0.
 */
/* Currently we expect that the messages arrive in order of the sender's rank. */
static void
t8_forest_ghost_parse_received_message (t8_forest_t forest, t8_forest_ghost_t ghost,
                                        t8_locidx_t *current_element_offset, int recv_rank, char *recv_buffer,
                                        int recv_bytes)
{
  size_t bytes_read, first_tree_index = 0, first_element_index = 0;
  t8_locidx_t num_trees, itree;
  t8_gloidx_t global_id;
  t8_eclass_t eclass;
  size_t num_elements, old_elem_count, ghosts_offset;
  t8_ghost_gtree_hash_t *tree_hash, **pfound_tree, *found_tree;
  t8_ghost_tree_t *ghost_tree;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_element_t *element_insert;
  t8_ghost_process_hash_t *process_hash;
#if T8_ENABLE_DEBUG
  int added_process;
#endif

  bytes_read = 0;
  /* read the number of trees */
  num_trees = *(size_t *) recv_buffer;
  bytes_read += sizeof (size_t);
  bytes_read += T8_ADD_PADDING (bytes_read);

  t8_debugf ("Received %li trees from %i (%i bytes)\n", (long) num_trees, recv_rank, recv_bytes);

  /* Count the total number of ghosts that we receive from this rank */
  ghosts_offset = ghost->num_ghosts_elements;
  for (itree = 0; itree < num_trees; itree++) {
    /* Get tree id */
    /* search if tree was inserted */
    /* if not: new entry, add elements */
    /* if yes: add the elements to the end of the tree's element array. */

    /* read the global id of this tree. */
    global_id = *(t8_gloidx_t *) (recv_buffer + bytes_read);
    bytes_read += sizeof (t8_gloidx_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* read the element class of the tree */
    eclass = *(t8_eclass_t *) (recv_buffer + bytes_read);
    bytes_read += sizeof (t8_eclass_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* read the number of elements sent */
    num_elements = *(size_t *) (recv_buffer + bytes_read);

    /* Add to the counter of ghost elements. */
    ghost->num_ghosts_elements += num_elements;

    bytes_read += sizeof (size_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* Search for the tree in the ghost_trees array */
    tree_hash = (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    tree_hash->global_id = global_id;

    /* Get the scheme for this tree */
    if (sc_hash_insert_unique (ghost->global_tree_to_ghost_tree, tree_hash, (void ***) &pfound_tree)) {
      /* The tree was not stored already, tree_hash is now an entry in the hash table. */
      /* If the tree was not contained, it is the newest tree in the array and
       * thus has as index the number of currently inserted trees. */
      tree_hash->index = ghost->ghost_trees->elem_count;
      found_tree = tree_hash;
      /* We grow the array by one and initialize the entry */
      ghost_tree = (t8_ghost_tree_t *) sc_array_push (ghost->ghost_trees);
      ghost_tree->global_id = global_id;
      ghost_tree->eclass = eclass;
      /* Initialize the element array */
      t8_element_array_init_size (&ghost_tree->elements, scheme, eclass, num_elements);
      /* pointer to where the elements are to be inserted */
      element_insert = t8_element_array_get_data_mutable (&ghost_tree->elements);
      /* Compute the element offset of this new tree by adding the offset
       * of the previous tree to the element count of the previous tree. */
      ghost_tree->element_offset = *current_element_offset;
      /* Allocate a new tree_hash for the next search */
      old_elem_count = 0;
      tree_hash = (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    }
    else {
      /* The entry was found in the trees array */
      found_tree = *pfound_tree;
      T8_ASSERT (found_tree->global_id == global_id);
      /* Get a pointer to the tree */
      ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees, found_tree->index);
      T8_ASSERT (ghost_tree->eclass == eclass);
      T8_ASSERT (ghost_tree->global_id == global_id);
      T8_ASSERT (ghost_tree->elements.tree_class == eclass);

      old_elem_count = t8_element_array_get_count (&ghost_tree->elements);

      /* Grow the elements array of the tree to fit the new elements */
      t8_element_array_resize (&ghost_tree->elements, old_elem_count + num_elements);
      /* Get a pointer to where the new elements are to be inserted */
      element_insert = t8_element_array_index_locidx_mutable (&ghost_tree->elements, old_elem_count);
    }

    if (itree == 0) {
      /* We store the index of the first tree and the first element of this
       * rank */
      first_tree_index = found_tree->index;
      first_element_index = old_elem_count;
    }
    /* Insert the new elements */
    memcpy (element_insert, recv_buffer + bytes_read, num_elements * scheme->get_element_size (eclass));

    bytes_read += num_elements * scheme->get_element_size (eclass);
    bytes_read += T8_ADD_PADDING (bytes_read);
    *current_element_offset += num_elements;
  }
  T8_ASSERT (bytes_read == (size_t) recv_bytes);
  T8_FREE (recv_buffer);

  /* At last we add the receiving rank to the ghosts process_offset hash table */
  process_hash = (t8_ghost_process_hash_t *) sc_mempool_alloc (ghost->proc_offset_mempool);
  process_hash->mpirank = recv_rank;
  process_hash->tree_index = first_tree_index;
  process_hash->first_element = first_element_index;
  process_hash->ghost_offset = ghosts_offset;
  /* Insert this rank into the hash table. We assert if the rank was not already contained. */
#if T8_ENABLE_DEBUG
  added_process =
#else
  (void)
#endif
    sc_hash_insert_unique (ghost->process_offsets, process_hash, NULL);
  T8_ASSERT (added_process);
}

/* In forest_ghost_receive we need a lookup table to give us the position
 * of a process in the ghost->remote_processes array, given the rank of a process.
 * We implement this via a hash table with the following struct as entry. */
typedef struct t8_recv_list_entry_struct
{
  int rank;                    /* The rank of this process */
  int pos_in_remote_processes; /* The position of this process in the remote_processes array */
} t8_recv_list_entry_t;

/* We hash these entries by their rank */
unsigned
t8_recv_list_entry_hash (const void *v1, [[maybe_unused]] const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;

  return e1->rank;
}

/* two entries are considered equal if they have the same rank. */
int
t8_recv_list_entry_equal (const void *v1, const void *v2, [[maybe_unused]] const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;
  const t8_recv_list_entry_t *e2 = (const t8_recv_list_entry_t *) v2;

  return e1->rank == e2->rank;
}

/* Probe for all incoming messages from the remote ranks and receive them.
 * We receive the message in the order in which they arrive. To achieve this,
 * we have to use polling. */
void
t8_forest_ghost_receive (t8_forest_t forest, t8_forest_ghost_t ghost)
{
  int num_remotes;
  int proc_pos;
  int recv_rank;
  int mpiret;
  sc_MPI_Comm comm;
  sc_MPI_Status status;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (ghost != NULL);

  comm = forest->mpicomm;
  /* Get the number of remote processes */
  num_remotes = ghost->remote_processes->elem_count;

  if (num_remotes == 0) {
    /* There is nothing to do */
    return;
  }

  {
    /* This code receives the message in order of their arrival.
     * This is effective in terms of runtime, but makes it more difficult
     * to store the received data, since the data has to be stored in order of
     * ascending ranks.
     * We include the received data into the ghost structure in order of the
     * ranks of the receivers and we do this as soon as the message from
     * the next rank that we can include was received.
     */
    char **buffer;
    int *recv_bytes;
    int received_messages = 0;
    int *received_flag;
    int last_rank_parsed = -1, parse_it;
#undef T8_POLLING /* activates polling for Mpi messages, currently used for testing */
#ifdef T8_POLLING
    sc_link_t *proc_it, *prev;
    int iprobe_flag;
    sc_list_t *receivers;
#else
    t8_recv_list_entry_t **pfound, *found;
#ifdef T8_ENABLE_DEBUG
    int ret;
#endif
    sc_hash_t *recv_list_entries_hash;
#endif
    t8_recv_list_entry_t recv_list_entry, *recv_list_entries;
    t8_locidx_t current_element_offset = 0;

    buffer = T8_ALLOC (char *, num_remotes);
    recv_bytes = T8_ALLOC (int, num_remotes);
    received_flag = T8_ALLOC_ZERO (int, num_remotes);
    recv_list_entries = T8_ALLOC (t8_recv_list_entry_t, num_remotes);

    /* Sort the array of remote processes, such that the ranks are in
     * ascending order. */
    sc_array_sort (ghost->remote_processes, sc_int_compare);

    /* We build a hash table of all ranks from which we receive and their position
     * in the remote_processes array. */
#ifdef T8_POLLING /* polling */
    receivers = sc_list_new (NULL);
#else
    recv_list_entries_hash = sc_hash_new (t8_recv_list_entry_hash, t8_recv_list_entry_equal, NULL, NULL);
#endif
    for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
      recv_list_entries[proc_pos].rank = *(int *) sc_array_index_int (ghost->remote_processes, proc_pos);
      recv_list_entries[proc_pos].pos_in_remote_processes = proc_pos;
#ifndef T8_POLLING
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_insert_unique (recv_list_entries_hash, recv_list_entries + proc_pos, NULL);
      T8_ASSERT (ret == 1);
#else /* polling */
      sc_list_append (receivers, recv_list_entries + proc_pos);
#endif
    }

    /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for a message for each
     * sender and if there is one we receive it and remove the sender from the list.
     * The last message can be received via probe */
#ifdef T8_POLLING
    while (received_messages < num_remotes - 1) {
      /* TODO: This part of the code using polling and IProbe to receive the
       *       messages. We replaced with a non-polling version that uses the
       *       blocking Probe. */
      iprobe_flag = 0;
      prev = NULL; /* ensure that if the first receive entry is matched first,
                                it is removed properly. */
      for (proc_it = receivers->first; proc_it != NULL && iprobe_flag == 0;) {
        /* pointer to the rank of a receiver */
        recv_rank = ((t8_recv_list_entry_t *) proc_it->data)->rank;
        proc_pos = ((t8_recv_list_entry_t *) proc_it->data)->pos_in_remote_processes;

        mpiret = sc_MPI_Iprobe (recv_rank, T8_MPI_GHOST_FOREST, comm, &iprobe_flag, &status);
        SC_CHECK_MPI (mpiret);
#else
    while (received_messages < num_remotes) {
      /* blocking probe for a message. */
      mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_GHOST_FOREST, comm, &status);
      SC_CHECK_MPI (mpiret);
#endif
#ifdef T8_POLLING
        if (iprobe_flag == 0) {
          /* There is no message to receive, we continue */
          prev = proc_it;
          proc_it = proc_it->next;
        }
        else {
#else
      /* There is a message to receive, we receive it. */
      recv_rank = status.MPI_SOURCE;
      /* Get the position of this rank in the remote processes array */
      recv_list_entry.rank = recv_rank;
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (recv_list_entries_hash, &recv_list_entry, (void ***) &pfound);
      T8_ASSERT (ret != 0);
      found = *pfound;
      proc_pos = found->pos_in_remote_processes;
#endif
          T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);
          buffer[proc_pos] = t8_forest_ghost_receive_message (recv_rank, comm, status, recv_bytes + proc_pos);
          /* mark this entry as received. */
          T8_ASSERT (received_flag[proc_pos] == 0);
          received_flag[proc_pos] = 1;
          received_messages++;
          /* Parse all messages that we can parse now.
           * We have to parse the messages in order of their rank. */
          T8_ASSERT (last_rank_parsed < proc_pos);
          /* For all ranks that we haven't parsed yet, but can be parsed in order */
          for (parse_it = last_rank_parsed + 1; parse_it < num_remotes && received_flag[parse_it] == 1; parse_it++) {
            recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
            t8_forest_ghost_parse_received_message (forest, ghost, &current_element_offset, recv_rank, buffer[parse_it],
                                                    recv_bytes[parse_it]);
            last_rank_parsed++;
          }

#ifdef T8_POLLING /* polling */
          /* Remove the process from the list of receivers. */
          proc_it = proc_it->next;
          sc_list_remove (receivers, prev);
        }
      } /* end for */
#endif
    } /* end while */

#ifdef T8_POLLING
    /* polling */
    T8_ASSERT (receivers->elem_count == 1);
    /* Get the last rank from which we didnt receive yet */
    recv_list_entry = *(t8_recv_list_entry_t *) sc_list_pop (receivers);
    recv_rank = recv_list_entry.rank;
    proc_pos = recv_list_entry.pos_in_remote_processes;
    /* destroy the list */
    sc_list_destroy (receivers);
    /* Blocking probe for the last message */
    mpiret = sc_MPI_Probe (recv_rank, T8_MPI_GHOST_FOREST, comm, &status);
    SC_CHECK_MPI (mpiret);
    /* Receive the message */
    T8_ASSERT (received_flag[proc_pos] == 0);
    buffer[proc_pos] = t8_forest_ghost_receive_message (recv_rank, comm, status, recv_bytes + proc_pos);
    received_flag[proc_pos] = 1;
    received_messages++;
    T8_ASSERT (received_messages == num_remotes);
    /* parse all messages that are left */
    /* For all ranks that we haven't parsed yet, but can be parsed in order */
    for (parse_it = last_rank_parsed + 1; parse_it < num_remotes && received_flag[parse_it] == 1; parse_it++) {
      recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
      t8_forest_ghost_parse_received_message (forest, ghost, &current_element_offset, recv_rank, buffer[parse_it],
                                              recv_bytes[parse_it]);
      last_rank_parsed++;
    }
#endif
#ifdef T8_ENABLE_DEBUG
    for (parse_it = 0; parse_it < num_remotes; parse_it++) {
      T8_ASSERT (received_flag[parse_it] == 1);
    }
#endif
    T8_ASSERT (last_rank_parsed == num_remotes - 1);

    /* clean-up */
#ifndef T8_POLLING
    sc_hash_destroy (recv_list_entries_hash);
#endif
    T8_FREE (buffer);
    T8_FREE (received_flag);
    T8_FREE (recv_list_entries);
    T8_FREE (recv_bytes);
  }
}
