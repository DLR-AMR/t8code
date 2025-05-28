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

#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_data/t8_containers.h>
#include <sc_statistics.h>

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_base.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_w_search.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** This struct is used during a ghost data exchange.
 * Since we use asynchronous communication, we store the
 * send buffers and mpi requests until we end the communication.
 */
typedef struct
{
  int num_remotes;
  /** The number of processes, we send to */
  char **send_buffers;
  /** For each remote the send buffer */
  sc_MPI_Request *send_requests;
  /** For each process we send to, the MPI request used */
  sc_MPI_Request *recv_requests;
  /** For each process we receive from, the MPI request used */
} t8_ghost_data_exchange_t;

/* return the number of trees in a ghost */
t8_locidx_t
t8_forest_ghost_num_trees (const t8_forest_t forest)
{
  if (forest->ghosts == NULL) {
    return 0;
  }
  T8_ASSERT (forest->ghosts != NULL);
  if (forest->ghosts->num_ghosts_elements <= 0) {
    return 0;
  }
  T8_ASSERT (forest->ghosts->ghost_trees != NULL);

  return forest->ghosts->ghost_trees->elem_count;
}

/* Given an index into the ghost_trees array return the ghost tree */
static t8_ghost_tree_t *
t8_forest_ghost_get_tree (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  t8_forest_ghost_t ghost;

  T8_ASSERT (t8_forest_is_committed (forest));
  ghost = forest->ghosts;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->ghost_trees != NULL);
  T8_ASSERT (0 <= lghost_tree && lghost_tree < t8_forest_ghost_num_trees (forest));

  ghost_tree = (t8_ghost_tree_t *) t8_sc_array_index_locidx (ghost->ghost_trees, lghost_tree);
  return ghost_tree;
}

t8_locidx_t
t8_forest_ghost_get_tree_element_offset (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return t8_forest_ghost_get_tree (forest, lghost_tree)->element_offset;
}

/* Given an index in the ghost_tree array, return this tree's number of elements */
t8_locidx_t
t8_forest_ghost_tree_num_elements (t8_forest_t forest, t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return t8_element_array_get_count (&ghost_tree->elements);
}

t8_element_array_t *
t8_forest_ghost_get_tree_elements (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  return &t8_forest_ghost_get_tree (forest, lghost_tree)->elements;
}

t8_locidx_t
t8_forest_ghost_get_ghost_treeid (t8_forest_t forest, t8_gloidx_t gtreeid)
{
  t8_ghost_gtree_hash_t query, *found, **pfound;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  query.global_id = gtreeid;
  if (sc_hash_lookup (forest->ghosts->global_tree_to_ghost_tree, &query, (void ***) &pfound)) {
    /* The tree was found */
    found = *pfound;
    return found->index;
  }
  else {
    /* The tree was not found */
    return -1;
  }
}

/* Given an index in the ghost_tree array, return this tree's element class */
t8_eclass_t
t8_forest_ghost_get_tree_class (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->eclass;
}

/* Given an index in the ghost_tree array, return this tree's global id */
t8_gloidx_t
t8_forest_ghost_get_global_treeid (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->global_id;
}

/* Given an index into the ghost_trees array and for that tree an element index,
 * return the corresponding element. */
t8_element_t *
t8_forest_ghost_get_element (t8_forest_t forest, t8_locidx_t lghost_tree, t8_locidx_t lelement)
{
  t8_ghost_tree_t *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  T8_ASSERT (0 <= lelement && lelement < t8_forest_ghost_tree_num_elements (forest, lghost_tree));
  /* TODO: In future, make return type const (and offer additional mutable version) and call t8_element_array_index_locidx (the const version). */
  return t8_element_array_index_locidx_mutable (&ghost_tree->elements, lelement);
}

void
t8_forest_ghost_create_ext (t8_forest_t forest)
{
  t8_forest_ghost_t ghost;
  t8_forest_ghost_definition_c *ghost_definition;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghost_definition != NULL);

  ghost_definition = forest->ghost_definition;

  t8_productionf ("Into t8_forest_ghost with %i local elements.\n", t8_forest_get_local_num_elements (forest));

  /* In parallel, check forest for deleted elements. The ghost algorithm currently
  * does not work on forests with deleted elements.
  * See also: https://github.com/DLR-AMR/t8code/issues/825
  * See also the test case: TODO Add a test case that currently fails. */
  SC_CHECK_ABORT (
    !forest->incomplete_trees || forest->mpisize == 1,
    "ERROR: Cannot compute ghost layer for forest with deleted elements (incomplete trees/holes in the mesh).\n");

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start ghost at %f  %f\n", sc_MPI_Wtime (), forest->profile->ghost_runtime);
  }
  /* Call the dot_ghost function on the ghost_definition class of the forest to compute the ghost layer */
  ghost_definition->do_ghost (forest);

  ghost = forest->ghosts;

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime += sc_MPI_Wtime ();
    /* We also store the number of ghosts and remotes */
    if (ghost != NULL) {
      forest->profile->ghosts_received = ghost->num_ghosts_elements;
      forest->profile->ghosts_shipped = ghost->num_remote_elements;
    }
    else {
      forest->profile->ghosts_received = 0;
      forest->profile->ghosts_shipped = 0;
    }
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End ghost at %f  %f\n", sc_MPI_Wtime (), forest->profile->ghost_runtime);
  }

  t8_productionf ("Done t8_forest_ghost with %i local elements and %i ghost elements.\n",
                  t8_forest_get_local_num_elements (forest), t8_forest_get_num_ghosts (forest));
}

void
t8_forest_ghost_create_topdown (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghost_definition != NULL);
  T8_ASSERT (t8_forest_ghost_definition_face_get_version (forest->ghost_definition) == 3);
  t8_forest_ghost_create_ext (forest);
}

/** Return the array of remote ranks.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in,out] num_remotes On output the number of remote ranks is stored here.
 * \return              The array of remote ranks in ascending order.
 */
int *
t8_forest_ghost_get_remotes (t8_forest_t forest, int *num_remotes)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->ghosts == NULL) {
    *num_remotes = 0;
    return NULL;
  }
  T8_ASSERT (forest->ghosts != NULL);

  *num_remotes = forest->ghosts->remote_processes->elem_count;
  return (int *) forest->ghosts->remote_processes->array;
}

/* Return the remote struct of a given remote rank */
t8_ghost_remote_t *
t8_forest_ghost_get_remote (t8_forest_t forest, int remote)
{
  t8_ghost_remote_t remote_search;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  size_t index;

  T8_ASSERT (t8_forest_is_committed (forest));

  remote_search.remote_rank = remote;
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (forest->ghosts->remote_ghosts, &remote_search, &index);
  T8_ASSERT (ret);
  return (t8_ghost_remote_t *) sc_array_index (&forest->ghosts->remote_ghosts->a, index);
}

/* Return a remote processes info about the stored ghost elements */
static t8_ghost_process_hash_t *
t8_forest_ghost_get_proc_info (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t proc_hash_search, **pproc_hash_found, *proc_hash_found;
#if T8_ENABLE_DEBUG
  int ret;
#endif

  T8_ASSERT (t8_forest_is_committed (forest));

  proc_hash_search.mpirank = remote;
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_lookup (forest->ghosts->process_offsets, &proc_hash_search, (void ***) &pproc_hash_found);
  T8_ASSERT (ret);
  proc_hash_found = *pproc_hash_found;
  T8_ASSERT (proc_hash_found->mpirank == remote);
  return proc_hash_found;
}

/** Return the first local ghost tree of a remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The ghost tree id of the first ghost tree that stores ghost
 *                      elements of \a remote.
 */
t8_locidx_t
t8_forest_ghost_remote_first_tree (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t *proc_entry;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  proc_entry = t8_forest_ghost_get_proc_info (forest, remote);
  T8_ASSERT (proc_entry->mpirank == remote);
  return proc_entry->tree_index;
}

/** Return the local index of the first ghost element that belongs to a given remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The index i in the ghost elements of the first element of rank \a remote
 */
t8_locidx_t
t8_forest_ghost_remote_first_elem (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t *proc_entry;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  proc_entry = t8_forest_ghost_get_proc_info (forest, remote);
  T8_ASSERT (proc_entry->mpirank == remote);

  return proc_entry->ghost_offset;
}

/* Fill the send buffer for a ghost data exchange for on remote rank.
 * returns the number of bytes in the buffer. */
static size_t
t8_forest_ghost_exchange_fill_send_buffer (t8_forest_t forest, int remote, char **pbuffer, sc_array_t *element_data)
{
  char *buffer;
  t8_ghost_remote_t lookup_rank, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_forest_ghost_t ghost;
  size_t index, element_index, data_size;
  size_t elements_inserted, byte_count;
  t8_tree_t local_tree;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  t8_locidx_t itree, ielement, element_pos;
  t8_locidx_t ltreeid;
  size_t elem_count;

  ghost = forest->ghosts;
  data_size = element_data->elem_size;
  elements_inserted = 0;
  lookup_rank.remote_rank = remote;

  /* Lookup the remote entry of this remote process */
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (ghost->remote_ghosts, &lookup_rank, &index);
  T8_ASSERT (ret != 0);
  remote_entry = (t8_ghost_remote_t *) sc_array_index (&ghost->remote_ghosts->a, index);
  T8_ASSERT (remote_entry->remote_rank == remote);

  /* allocate memory for the send buffer */
  byte_count = data_size * remote_entry->num_elements;
  buffer = *pbuffer = T8_ALLOC (char, byte_count);

  /* We now iterate over the remote trees and their elements to find the
   * local element indices of the remote elements */
  for (itree = 0; itree < (t8_locidx_t) remote_entry->remote_trees.elem_count; itree++) {
    /* tree loop */
    remote_tree = (t8_ghost_remote_tree_t *) t8_sc_array_index_locidx (&remote_entry->remote_trees, itree);
    /* Get the local id of this tree */
    /* TODO: Why does remote_tree store the global id? could be local instead */
    ltreeid = t8_forest_get_local_id (forest, remote_tree->global_id);
    /* Get a pointer to the forest tree */
    local_tree = t8_forest_get_tree (forest, ltreeid);
    elem_count = t8_element_array_get_count (&remote_tree->elements);
    for (ielement = 0; ielement < (t8_locidx_t) elem_count; ielement++) {
      /* element loop */
      /* Get the index of this remote element in its local tree */
      element_pos = *(t8_locidx_t *) t8_sc_array_index_locidx (&remote_tree->element_indices, ielement);
      T8_ASSERT (0 <= element_pos);
      /* Compute the index of this element in the element_data array */
      element_index = local_tree->elements_offset + element_pos;
      /* Copy the data of this element from the element_data array to the send buffer */
      memcpy (buffer + elements_inserted * data_size, sc_array_index (element_data, element_index), data_size);
      elements_inserted++;
    }
  }
  return byte_count;
}

static t8_ghost_data_exchange_t *
t8_forest_ghost_exchange_begin (t8_forest_t forest, sc_array_t *element_data)
{
  t8_ghost_data_exchange_t *data_exchange;
  t8_forest_ghost_t ghost;
  size_t bytes_to_send, ghost_start;
  int iremote, remote_rank;
  int mpiret, recv_rank, bytes_recv;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  char **send_buffers;
  t8_ghost_process_hash_t lookup_proc, *process_entry, **pfound;
  t8_locidx_t remote_offset, next_offset;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (element_data != NULL);
  T8_ASSERT (forest->ghosts != NULL);

  ghost = forest->ghosts;

  /* Allocate the new exchange context */
  data_exchange = T8_ALLOC (t8_ghost_data_exchange_t, 1);
  /* The number of processes we need to send to */
  data_exchange->num_remotes = ghost->remote_processes->elem_count;
  /* Allocate MPI requests */
  data_exchange->send_requests = T8_ALLOC (sc_MPI_Request, data_exchange->num_remotes);
  data_exchange->recv_requests = T8_ALLOC (sc_MPI_Request, data_exchange->num_remotes);
  /* Allocate pointers to send buffers */
  send_buffers = data_exchange->send_buffers = T8_ALLOC (char *, data_exchange->num_remotes);

  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* Iterate over all remote processes and fill their send buffers */
    remote_rank = *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    /* Fill the send buffers and compute the number of bytes to send */
    bytes_to_send
      = t8_forest_ghost_exchange_fill_send_buffer (forest, remote_rank, send_buffers + iremote, element_data);

    /* Post the asynchronuos send */
    mpiret = sc_MPI_Isend (send_buffers[iremote], bytes_to_send, sc_MPI_BYTE, remote_rank, T8_MPI_GHOST_EXC_FOREST,
                           forest->mpicomm, data_exchange->send_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }

  /* The index in element_data at which the ghost elements start */
  ghost_start = t8_forest_get_local_num_elements (forest);
  /* Receive the incoming messages */
  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* We need to compute the offset in element_data to which we can receive the message */
    /* Search for this processes' entry in the ghost struct */
    recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    lookup_proc.mpirank = recv_rank;
#if T8_ENABLE_DEBUG
    ret =
#else
    (void)
#endif
      sc_hash_lookup (ghost->process_offsets, &lookup_proc, (void ***) &pfound);
    T8_ASSERT (ret);
    process_entry = *pfound;
    /* In process_entry we stored the offset of this ranks ghosts under all
     * ghosts. Thus in element_data we look at the position
     *  ghost_start + offset
     */
    remote_offset = process_entry->ghost_offset;
    /* Compute the offset of the next remote rank */
    if (iremote + 1 < data_exchange->num_remotes) {
      lookup_proc.mpirank = *(int *) sc_array_index_int (ghost->remote_processes, iremote + 1);
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &lookup_proc, (void ***) &pfound);
      T8_ASSERT (ret);
      process_entry = *pfound;
      next_offset = process_entry->ghost_offset;
    }
    else {
      /* We are the last rank, the next offset is the total number of ghosts */
      next_offset = ghost->num_ghosts_elements;
    }
    /* Calculate the number of bytes to receive */
    bytes_recv = (next_offset - remote_offset) * element_data->elem_size;
    /* receive the message */
    mpiret = sc_MPI_Irecv (sc_array_index (element_data, ghost_start + remote_offset), bytes_recv, sc_MPI_BYTE,
                           recv_rank, T8_MPI_GHOST_EXC_FOREST, forest->mpicomm, data_exchange->recv_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }
  return data_exchange;
}

static void
t8_forest_ghost_exchange_end (t8_ghost_data_exchange_t *data_exchange)
{
  int iproc;

  T8_ASSERT (data_exchange != NULL);
  /* Wait for all communications to end */
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->recv_requests, sc_MPI_STATUSES_IGNORE);
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->send_requests, sc_MPI_STATUSES_IGNORE);

  /* Free the send buffers */
  for (iproc = 0; iproc < data_exchange->num_remotes; iproc++) {
    T8_FREE (data_exchange->send_buffers[iproc]);
  }
  T8_FREE (data_exchange->send_buffers);
  /* free requests */
  T8_FREE (data_exchange->send_requests);
  T8_FREE (data_exchange->recv_requests);
  T8_FREE (data_exchange);
}

void
t8_forest_ghost_exchange_data (t8_forest_t forest, sc_array_t *element_data)
{
  t8_ghost_data_exchange_t *data_exchange;

  t8_debugf ("Entering ghost_exchange_data\n");
  T8_ASSERT (t8_forest_is_committed (forest));

  if (forest->ghosts == NULL) {
    /* This process has no ghosts */
    return;
  }

  T8_ASSERT (forest->ghosts != NULL);
  T8_ASSERT (element_data != NULL);
  T8_ASSERT ((t8_locidx_t) element_data->elem_count
             == t8_forest_get_local_num_elements (forest) + t8_forest_get_num_ghosts (forest));

  data_exchange = t8_forest_ghost_exchange_begin (forest, element_data);
  if (forest->profile != NULL) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime = -sc_MPI_Wtime ();
  }
  t8_forest_ghost_exchange_end (data_exchange);
  if (forest->profile != NULL) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime += sc_MPI_Wtime ();
  }
  t8_debugf ("Finished ghost_exchange_data\n");
}

/* Print a forest ghost structure */
void
t8_forest_ghost_print (t8_forest_t forest)
{
  t8_forest_ghost_t ghost;
  t8_ghost_remote_t *remote_found;
  t8_ghost_remote_tree_t *remote_tree;
  t8_ghost_process_hash_t proc_hash, **pfound, *found;
  size_t iremote, itree;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  int remote_rank;
  char remote_buffer[BUFSIZ] = "";
  char buffer[BUFSIZ] = "";

  if (forest->ghosts == NULL) {
    return;
  }
  T8_ASSERT (forest->ghosts != NULL);
  ghost = forest->ghosts;
  snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer), "\tRemotes:\n");
  snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "\tReceived:\n");

  if (ghost->num_ghosts_elements > 0) {
    for (iremote = 0; iremote < ghost->remote_processes->elem_count; iremote++) {
      /* Get the rank of the remote process */
      remote_rank = *(int *) sc_array_index (ghost->remote_processes, iremote);
      /* Get this remote's entry */
      remote_found = t8_forest_ghost_get_remote (forest, remote_rank);
      /* investigate the entry of this remote process */
      snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer), "\t[Rank %i] (%li trees):\n",
                remote_found->remote_rank, remote_found->remote_trees.elem_count);
      for (itree = 0; itree < remote_found->remote_trees.elem_count; itree++) {
        remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_found->remote_trees, itree);
        snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer),
                  "\t\t[id: %lli, class: %s, #elem: %li]\n", (long long) remote_tree->global_id,
                  t8_eclass_to_string[remote_tree->eclass], (long) t8_element_array_get_count (&remote_tree->elements));
      }

      /* Investigate the elements that we received from this process */
      proc_hash.mpirank = remote_rank;
      /* look up this rank in the hash table */
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &proc_hash, (void ***) &pfound);

      T8_ASSERT (ret);
      found = *pfound;
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                "\t[Rank %i] First tree: %li\n\t\t First element: %li\n", remote_rank, (long) found->tree_index,
                (long) found->first_element);
    }
  }
  t8_debugf ("Ghost structure:\n%s\n%s\n", remote_buffer, buffer);
}

/* Completely destroy a ghost structure */
static void
t8_forest_ghost_reset (t8_forest_ghost_t *pghost)
{
  t8_forest_ghost_t ghost;
  size_t it, it_trees;
  t8_ghost_tree_t *ghost_tree;
  t8_ghost_remote_t *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->rc.refcount == 0);

  /* Clean-up the arrays */
  for (it_trees = 0; it_trees < ghost->ghost_trees->elem_count; it_trees++) {
    ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees, it_trees);
    t8_element_array_reset (&ghost_tree->elements);
  }

  sc_array_destroy (ghost->ghost_trees);
  sc_array_destroy (ghost->remote_processes);
  /* Clean-up the hashtables */
  sc_hash_destroy (ghost->global_tree_to_ghost_tree);
  sc_hash_destroy (ghost->process_offsets);
  /* Clean-up the remote ghost entries */
  for (it = 0; it < ghost->remote_ghosts->a.elem_count; it++) {
    remote_entry = (t8_ghost_remote_t *) sc_array_index (&ghost->remote_ghosts->a, it);
    for (it_trees = 0; it_trees < remote_entry->remote_trees.elem_count; it_trees++) {
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees, it_trees);
      t8_element_array_reset (&remote_tree->elements);
      sc_array_reset (&remote_tree->element_indices);
    }
    sc_array_reset (&remote_entry->remote_trees);
  }
  sc_hash_array_destroy (ghost->remote_ghosts);

  /* Clean-up the memory pools for the data inside
   * the hash tables */
  sc_mempool_destroy (ghost->glo_tree_mempool);
  sc_mempool_destroy (ghost->proc_offset_mempool);

  /* Free the ghost */
  T8_FREE (ghost);
  pghost = NULL;
}

void
t8_forest_ghost_ref (t8_forest_ghost_t ghost)
{
  T8_ASSERT (ghost != NULL);

  t8_refcount_ref (&ghost->rc);
}

void
t8_forest_ghost_unref (t8_forest_ghost_t *pghost)
{
  t8_forest_ghost_t ghost;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);

  if (t8_refcount_unref (&ghost->rc)) {
    t8_forest_ghost_reset (pghost);
  }
}

void
t8_forest_ghost_destroy (t8_forest_ghost_t *pghost)
{
  T8_ASSERT (pghost != NULL && *pghost != NULL && t8_refcount_is_last (&(*pghost)->rc));
  t8_forest_ghost_unref (pghost);
  T8_ASSERT (*pghost == NULL);
}

T8_EXTERN_C_END ();
