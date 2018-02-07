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

#include <t8_data/t8_shmem.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"
#include "t8_cmesh_partition.h"
#include "t8_cmesh_offset.h"

#if 0
/* Return the minimum of two t8_gloidx_t's */
static              t8_gloidx_t
t8_glo_min (t8_gloidx_t A, t8_gloidx_t B)
{
  return A < B ? A : B;
}
#endif

/* Change the neighbor entry of a tree to match the new partition.
 * Input: A face_neighbor entry in cmesh_from and a process to which the corresponding tree will be send
 * Output: The face neighbor entry is changed to match its new id in cmesh.
 */
static void
t8_cmesh_partition_send_change_neighbor (t8_cmesh_t cmesh,
                                         t8_cmesh_t cmesh_from,
                                         t8_locidx_t * neighbor, int to_proc)
{
  t8_gloidx_t         temp;
  t8_gloidx_t        *tree_offset;

  tree_offset = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (0 <= *neighbor && *neighbor < cmesh_from->num_local_trees) {
    /* Neighbor is a local tree in cmesh */
    temp = cmesh_from->first_tree - t8_offset_first (to_proc, tree_offset);
    /* Assert for possible overflow du to gloidx computation */
    T8_ASSERT ((t8_locidx_t) (temp + *neighbor) == temp + *neighbor);
    *neighbor = temp + *neighbor;
  }
  else {
    t8_cghost_t         ghost;
    /* neighbor is a ghost in cmesh_from */
    T8_ASSERT (*neighbor >= cmesh_from->num_local_trees &&
               *neighbor <
               cmesh_from->num_local_trees + cmesh_from->num_ghosts);
    ghost =
      t8_cmesh_trees_get_ghost (cmesh_from->trees,
                                *neighbor - cmesh_from->num_local_trees);
    /* We only do something if the neighbor will be a local tree
     * of to_proc in cmesh */
    if (t8_offset_in_range (ghost->treeid, to_proc, tree_offset)) {
      /* The new neighbor id is The global index of ghost - first tree of cmesh
       */
      temp = ghost->treeid - t8_offset_first (to_proc, tree_offset);
      /* assert for gloidx overflow */
      T8_ASSERT ((t8_locidx_t) temp == temp);
      *neighbor = temp;
    }
  }
}

/* After we received all parts from the other processes we have to compute
 * the new local ids of the ghosts and set those in the neighbor fields
 * of all local trees.
 * We also insert their global ids into the hash table of global_id -> local_id
 */
static void
t8_partition_new_ghost_ids (t8_cmesh_t cmesh,
                            t8_part_tree_t recv_part,
                            t8_locidx_t first_ghost, int proc)
{
  t8_locidx_t         ghost_it;
  t8_cghost_t         ghost;
  t8_gloidx_t        *face_neighbors, tree_id_glo;
  t8_locidx_t        *tree_neighbors;
  int8_t             *ttf;
  int                 iface, face_tree;
  t8_trees_glo_lo_hash_t *new_hash;
  int                 ret;

  for (ghost_it = 0; ghost_it < recv_part->num_ghosts; ghost_it++) {
    /* loop over all ghosts of recv_part */
    cmesh->trees->ghost_to_proc[ghost_it + first_ghost] = proc;
    ghost =
      t8_cmesh_trees_get_ghost_ext (cmesh->trees, first_ghost + ghost_it,
                                    &face_neighbors, &ttf);
    for (iface = 0; iface < t8_eclass_num_faces[ghost->eclass]; iface++) {
      /* loop over all faces of ghost */
      tree_id_glo = face_neighbors[iface];
      if (cmesh->first_tree <= tree_id_glo &&
          tree_id_glo < cmesh->first_tree + cmesh->num_local_trees) {
        /* the face neighbor is a local tree */
        (void) t8_cmesh_trees_get_tree_ext (cmesh->trees,
                                            tree_id_glo - cmesh->first_tree,
                                            &tree_neighbors, NULL);
        /* Get the number of the face of tree that is connected with ghost
         * and set the new local ghost id */
        face_tree = ttf[iface] % t8_eclass_max_num_faces[cmesh->dimension];
        tree_neighbors[face_tree] =
          ghost_it + first_ghost + cmesh->num_local_trees;
      }
    }
    /* Insert this ghost's global and local id into the hash table */
    new_hash = (t8_trees_glo_lo_hash_t *)
      sc_mempool_alloc (cmesh->trees->global_local_mempool);
    new_hash->global_id = ghost->treeid;
    /* The new local ghost id is the concurrent id of this ghost plus the
     * number of local trees */
    new_hash->local_id = ghost_it + first_ghost + cmesh->num_local_trees;
    ret = sc_hash_insert_unique (cmesh->trees->ghost_globalid_to_local_id,
                                 new_hash, NULL);
    t8_debugf ("[H] Added global id %li local id %i %p\n", ghost->treeid,
               new_hash->local_id, cmesh->trees->ghost_globalid_to_local_id);
    /* The entry must not have existed before */
    T8_ASSERT (ret);
  }
}

/* From num_local_trees_per_eclass compute num_trees_per_eclass.
 * collective function */
void
t8_cmesh_gather_trees_per_eclass (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_gloidx_t         temp_trees_per_eclass[T8_ECLASS_COUNT];
  int                 ieclass;

  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  if (cmesh->set_partition) {
    /* Copy the local values */
    /* We need to do it in a loop since we convert from locidx to gloidx.
     * memcpy is thus not possible */
    for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
      temp_trees_per_eclass[ieclass] =
        cmesh->num_local_trees_per_eclass[ieclass];
    }

    if (cmesh->first_tree_shared) {
      t8_eclass_t         eclass;
      T8_ASSERT (cmesh->num_local_trees > 0);
      /* If our first tree is shared, we must not count it in the
       * global trees_per_eclass field */
      /* We need to use t8_cmesh_trees_get_tree instead of t8_cmesh_get_tree
       * since the latter performs a is_committed check on cmesh, but we may
       * call t8_cmesh_gather_trees_per_eclass when the cmesh is not
       * yet fully committed.
       * (Since a cmesh only counts as committed if its trees_per_eclass
       * values are set.)
       */
      eclass = t8_cmesh_trees_get_tree (cmesh->trees, 0)->eclass;
      temp_trees_per_eclass[eclass]--;
    }
    sc_MPI_Allreduce (temp_trees_per_eclass, cmesh->num_trees_per_eclass,
                      T8_ECLASS_COUNT, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  }
  else {
    /* The cmesh is not partitioned, we can just copy local_num_trees_per_eclass */
    for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
      cmesh->num_trees_per_eclass[ieclass] =
        cmesh->num_local_trees_per_eclass[ieclass];
    }
  }
#ifdef T8_ENABLE_DEBUG
  /* Count the number of trees and check if it matches cmesh->num_trees */
  {
    t8_gloidx_t         num_trees = 0;
    for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
      num_trees += cmesh->num_trees_per_eclass[ieclass];
    }
    T8_ASSERT (num_trees == cmesh->num_trees);
  }
#endif
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process,
 * additional flag whether we compute the trees per eclass or not
 * additional flag whether to check if cmesh is committed.
 * Warning: use with caution with check_commit = 0 */
static void
t8_cmesh_gather_treecount_ext (t8_cmesh_t cmesh, sc_MPI_Comm comm,
                               int check_commit)
{
  t8_gloidx_t         tree_offset;
  t8_gloidx_t        *tree_offset_array;
  int                 is_empty, has_empty;

  if (check_commit) {
    T8_ASSERT (t8_cmesh_is_committed (cmesh));
  }
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  tree_offset = cmesh->first_tree_shared ? -cmesh->first_tree - 1 :
    cmesh->first_tree;
  if (cmesh->tree_offsets == NULL) {
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    /* Only allocate the shmem array, if it is not already allocated */
    cmesh->tree_offsets = t8_cmesh_alloc_offsets (cmesh->mpisize, comm);
    t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX,
                              cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
    t8_shmem_array_set_gloidx (cmesh->tree_offsets, cmesh->mpisize,
                               cmesh->num_trees);

    if (cmesh->num_local_trees <= 0) {
      /* This process is empty */
      is_empty = 1;
    }
    else {
      is_empty = 0;
    }
    /* Communicate whether we have empty processes */
    sc_MPI_Allreduce (&is_empty, &has_empty, 1, sc_MPI_INT, sc_MPI_LOR, comm);
    if (has_empty) {
      int                 next_nonempty;

      tree_offset_array =
        t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
      /* there exist empty ranks, we have to recalculate the offset.
       * Each empty rank stores the offset of the next nonempty rank */
      if (is_empty) {
        next_nonempty =
          t8_offset_next_nonempty_rank (cmesh->mpirank, cmesh->mpisize,
                                        tree_offset_array);
        /* Set the tree offset to the first nonshared tree of the next rank */
        tree_offset = t8_offset_first (next_nonempty, tree_offset_array);
        if (tree_offset_array[next_nonempty] < 0) {
          tree_offset++;
        }
      }
      /* Communicate the new tree offsets */
      t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX,
                                cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
    }
  }
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process */
void
t8_cmesh_gather_treecount (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_gather_treecount_ext (cmesh, comm, 1);
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process */
void
t8_cmesh_gather_treecount_nocommit (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_gather_treecount_ext (cmesh, comm, 0);
}

#if 0
/* If for a cmesh set_partition_range was called, create
 * from this information the complete partition table */
/* TODO: For this function each process needs to communicate with
 *       both of its nearest nonempty neighbors.
 *       This is too much communication.
 *       Usually if we have shared trees then these come from an
 *       undelying forest, and it is thus the responsibility of the
 *       forest to set the correct offsets.
 *       Therefore, if set_partition_range was called, we will now allways assume
 *       that there are no trees.
 *       This function is thus obsolete.
 */
static void
t8_cmesh_partition_create_offsets (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_gloidx_t         first_tree, *offset;
  t8_locidx_t         num_trees_prev, recv_data[2], send_data[2];
  int                 mpiret, first_tree_shared, proc_prev;
  sc_MPI_Request      mpi_request;

  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  T8_ASSERT (cmesh->face_knowledge == -1 || cmesh->face_knowledge == 3);
  T8_ASSERT (cmesh->first_tree >= 0);
  T8_ASSERT (cmesh->num_local_trees >= 0);
  T8_ASSERT (cmesh->tree_offsets == NULL);

  /* Initialize the tree offset array */
  t8_shmem_array_init (&cmesh->tree_offsets, sizeof (t8_gloidx_t),
                       cmesh->mpisize + 1, comm);
  /* Store the first local tree of each process in the trees array */
  first_tree = cmesh->first_tree;
  t8_shmem_array_allgather (&first_tree, 1, T8_MPI_GLOIDX,
                            cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
  t8_shmem_array_set_gloidx (cmesh->tree_offsets, cmesh->mpisize,
                             cmesh->num_trees);
  offset = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);

  /* Now we need to find out if our first tree is shared with other processes */
  /* To achieve this, we need to now the number of trees from the previous,
   * nonempty process */
  num_trees_prev = 0;

  if (cmesh->mpirank != 0) {
    /* We receive the number of trees of the next nonempty process */
    mpiret = sc_MPI_Irecv (recv_data, 2, T8_MPI_LOCIDX,
                           cmesh->mpirank - 1, 0, comm, &mpi_request);
    SC_CHECK_MPI (mpiret);
  }
  if (cmesh->num_local_trees == 0 && cmesh->mpirank > 0) {
    /* If we do not have any tree, we wait to receive the number of trees
     * from the previous process and send them to the next one */
    /* Wait until the message is here */
    mpiret = sc_MPI_Wait (&mpi_request, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    mpi_request = sc_MPI_REQUEST_NULL;
    send_data[0] = recv_data[0];
    send_data[1] = recv_data[1];
  }
  else {
    /* If we do have local trees, we send this number to the next process */
    send_data[0] = cmesh->mpirank;
    send_data[1] = cmesh->num_local_trees;
  }
  if (cmesh->mpirank != cmesh->mpisize - 1) {
    mpiret = sc_MPI_Send (&send_data, 2, T8_MPI_LOCIDX, cmesh->mpirank + 1, 0,
                          comm);
    SC_CHECK_MPI (mpiret);
  }
  if (cmesh->mpirank > 0) {
    /* Wait until the message is here */
    mpiret = sc_MPI_Wait (&mpi_request, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    proc_prev = recv_data[0];
    num_trees_prev = recv_data[1];
  }
  else {
    proc_prev = -1;
    num_trees_prev = 0;
  }
  /* Calculate whether our first tree is shared */
  first_tree_shared = 0;
  t8_debugf ("proc prev %i, num_trees_prev %i\n", proc_prev, num_trees_prev);
  if (num_trees_prev != 0) {
    if (t8_offset_first (proc_prev, offset) + num_trees_prev > first_tree) {
      T8_ASSERT (t8_offset_first (proc_prev, offset) + num_trees_prev
                 == first_tree + 1);
      first_tree_shared = 1;
    }
  }
  cmesh->first_tree_shared = first_tree_shared;
  /* Given the first_tree_shared info calculate the correct entry for
   * this process in the offset array */
  first_tree = t8_offset_first_tree_to_entry (first_tree, first_tree_shared);
  /* Allgather the new first tree entries */
  t8_shmem_array_allgather (&first_tree, 1, T8_MPI_GLOIDX,
                            cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
}
#endif

/* TODO: currently this function is unused.
 *        Also it better fits to cmesh_offset.c/h */
#if 0
/* Return a process that a given process definitely sends to/receives from */
static int
t8_cmesh_partition_send_any (int proc, t8_cmesh_t cmesh_from,
                             t8_gloidx_t * offset_from,
                             t8_gloidx_t * offset_to, int receive)
{
  int                 lookhere, range[2], *sender, *receiver;
  int                 search_dir;
  int                 last;
  int                 done = 0;

  T8_ASSERT (proc >= 0 && proc < cmesh_from->mpisize);

  if ((!receive && t8_offset_nosend (proc, cmesh_from->mpisize, offset_from,
                                     offset_to)) ||
      (receive
       && t8_offset_empty (proc, receive ? offset_to : offset_from))) {
    /* We do not send/recv anything */
    return -1;
  }

  range[0] = 0;
  range[1] = cmesh_from->mpisize;
  lookhere = (range[0] + range[1]) / 2;

  if (!receive) {
    sender = &proc;
    receiver = &lookhere;
  }
  else {
    sender = &lookhere;
    receiver = &proc;
  }
  search_dir = +1;
  /* If we send our reveive from ourselves, we use it as our guess */
  if (t8_offset_sendsto (cmesh_from->mpirank, cmesh_from->mpirank,
                         offset_from, offset_to)) {
    return cmesh_from->mpirank;
  }
  while (!done
         && !t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
    last = lookhere;
    while ((receive && t8_offset_nosend (lookhere, cmesh_from->mpisize,
                                         offset_from, offset_to))
           || (!receive && t8_offset_empty (lookhere, offset_to))) {
      /* Skip processes that do not send/recv anything */
      lookhere += search_dir;
      if (lookhere < 0) {
        /* There is no candidate to the left of last */
        range[0] = last + 1;
        lookhere = (range[0] + range[1]) / 2;
        last = lookhere;
        search_dir = 1;
      }
      else if (lookhere >= cmesh_from->mpisize) {
        /* There is no candidate to the right of last */
        range[1] = last - 1;
        lookhere = (range[0] + range[1]) / 2;
        last = lookhere;
        search_dir = -1;
      }
    }
    if (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      if (t8_offset_last (*sender, offset_from)
          < t8_offset_first (*receiver, offset_to)
          + (offset_from[*receiver] < 0 && *sender != *receiver
             && t8_offset_in_range (t8_offset_first (*receiver, offset_to),
                                    *receiver, offset_from))) {

        /* If the last tree we could send is smaller than the first we could
         * receive then we have to make receiver smaller or sender bigger */
        if (!receive) {
          range[1] = SC_MIN (lookhere - 1, range[1] - 1);
          search_dir = -1;
        }
        else {
          range[0] = SC_MAX (lookhere + 1, range[0] + 1);
          search_dir = +1;
        }
      }
      else {
        if (!receive) {
          range[0] = SC_MAX (lookhere, range[0] + (search_dir > 0));
          search_dir = +1;
        }
        else {
          range[1] = SC_MIN (lookhere, range[1] - (search_dir < 0));
          search_dir = -1;
        }
      }
    }
    else {                      /* t8_offset_sendsto */
      range[0] = range[1] = lookhere;
    }
    if (range[0] <= range[1]) {
      lookhere = (range[0] + range[1]) / 2;
    }
    if (range[0] >= range[1]) {
      done = 1;
    }
    t8_debugf ("%i\n ", lookhere);
    T8_ASSERT (0 <= lookhere && lookhere < cmesh_from->mpisize);
#if 0
    if (last == lookhere) {
      /* We did not change in this round which means, that we found the rank */
      done = 1;
    }
#endif
  }
  if (done && !t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
    /* We were not able to find a process, so we do a linear search instead */
    /* TODO: Currently the function above has some flaws such that in rare cases,
     *       no process is found. We should eventually repair it. */
    int                 searcher[2];    /* We search from mpirank in both directions */
    int                 pos = 0;
    t8_debugf
      ("[H] Could not find any process in bin search. Start linear search.\n");
    lookhere = cmesh_from->mpirank;
    search_dir = 1;
    searcher[0] = cmesh_from->mpirank - 1;
    searcher[1] = cmesh_from->mpirank + 1;
    while (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      T8_ASSERT (0 <= searcher[0] || searcher[1] < cmesh_from->mpisize);
      if (searcher[0] < 0) {
        /* Always search below, if the top search reached its end */
        pos = 1;
      }
      else if (searcher[1] >= cmesh_from->mpirank) {
        /* Always search on top, if the below search reached its end */
        pos = 0;
      }
      else {
        /* In each step we change the search direction */
        pos = 1 - pos;
      }
      lookhere = searcher[pos];
      searcher[pos] += 2 * pos - 1;     /* if pos = 0, substract one, if pos = 1, add 1 */
    }
  }
  return lookhere;
}
#endif

#if 0
/* Compute first and last process to which we will send data */
/* Returns the local tree_id of last local tree on send_first. */
static              t8_locidx_t
t8_cmesh_partition_sendrange (t8_cmesh_t cmesh_to, t8_cmesh_t cmesh_from,
                              int *send_first, int *send_last, int receive)
{
  t8_gloidx_t         ret = -1;
  int                 range[2], lookhere, first_guess, search_dir = -1, last;
  int                *sender, *receiver;
  t8_gloidx_t        *offset_from, *offset_to;

  if (!cmesh_from->set_partition) {
    /* If cmesh_from is not partitioned we can send/recv all information
     * to/from ourself */
    *send_first = *send_last = cmesh_from->mpirank;
    offset_to = t8_shmem_array_get_gloidx_array (cmesh_to->tree_offsets);
    return t8_offset_last (cmesh_from->mpirank, offset_to);
  }

  *send_first = -1;
  *send_last = -1;
  range[0] = 0;
  range[1] = cmesh_from->mpisize - 1;

  /* Binary search the smallest process to which we send something */

  offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  offset_to = t8_shmem_array_get_gloidx_array (cmesh_to->tree_offsets);

  /* We start with any process that we send to/recv from */
  first_guess = t8_cmesh_partition_send_any (cmesh_from->mpirank, cmesh_from,
                                             offset_from, offset_to, receive);
  t8_debugf ("first guess %i\n", first_guess);
  T8_ASSERT (first_guess == -1
             || (!receive
                 && t8_offset_sendsto (cmesh_from->mpirank, first_guess,
                                       offset_from, offset_to))
             || (receive
                 && t8_offset_sendsto (first_guess, cmesh_from->mpirank,
                                       offset_from, offset_to)));

  if (first_guess == -1 || (!receive && cmesh_from->num_local_trees == 0)
      || (receive && cmesh_to->num_local_trees == 0)) {
    /* This partition is empty and we can not send/recv anything */
    *send_last = -2;
    return -1;
  }
  T8_ASSERT (first_guess >= 0 && first_guess < cmesh_from->mpisize);
  if (!receive) {
    /* in sendmode we are the sender and look for receivers */
    sender = &cmesh_from->mpirank;
    receiver = &lookhere;
  }
  else {
    /* in receive mode we are receiver and look for senders */
    sender = &lookhere;
    receiver = &cmesh_from->mpirank;
  }

  range[1] = first_guess;
  /* We start by probing our first_guess */
  lookhere = first_guess;
  while (*send_first == -1) {
    /* However it could be empty in this case we search linearly to
     * find the next nonempty receiver in search direction */
    /* TODO: replace nonempty with nonsending? */
    while (1 <= lookhere && lookhere < cmesh_from->mpisize - 1 &&
           (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)
            /* TODO: The two next lines can be removed? */
            || (receive && t8_offset_num_trees (lookhere, offset_from) == 1
                && offset_from[lookhere] < 0 && lookhere != *receiver))) {
      lookhere += search_dir;
    }
    if (t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      /* We send to this process */
      /* We have to look further left */
      range[1] = SC_MIN (lookhere, range[1]);
      search_dir = -1;
    }
    else {
      /* We do not send to this process */
      /* We have to look further right */
      range[0] = SC_MAX (lookhere + 1, range[0] + 1);
      search_dir = +1;
    }
    if (range[0] == range[1]) {
      /* The only possibility left is send_first */
      *send_first = range[0];
    }
    /* Our next guess is the midpoint between range[0] and [1] */
    lookhere = (range[0] + range[1]) / 2;
  }
  range[0] = first_guess;
  range[1] = cmesh_from->mpisize - 1;
  search_dir = +1;
  lookhere = first_guess;

  last = -2;                    /* If this does not change we found our process */
  while (*send_last == -1) {
    /* Our guess could be empty in this case we search linearly to
     * find the next nonempty receiver in search direction */
    while (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      t8_debugf ("lookhere %i\n", lookhere);
      /* We may be already at the end or beginning in which
       * case the search direction is the opposite one */
      if (lookhere == cmesh_from->mpisize - 1) {
        search_dir = -1;
      }
      if (lookhere == 0) {
        search_dir = +1;
      }
      lookhere += search_dir;
    }
    if (lookhere == last) {
      /* We did not change, so we found the right process */
      *send_last = lookhere;
    }
    else if (t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      /* We have to look further right */
      range[0] = SC_MAX (lookhere, range[0]);
      search_dir = +1;
    }
    else {
      /* We have to look further left */
      range[1] = SC_MIN (lookhere - 1, range[1] - 1);
      search_dir = -1;
    }
    if (range[0] == range[1]) {
      /* The only possibility left is send_last */
      T8_ASSERT (0 <= range[0] && range[0] < cmesh_from->mpisize);
      *send_last = range[0];
    }
    last = lookhere;
    /* Our next guess is the midpoint between range[0] and [1] */
    /* We use ceil, since we always want to look more to the right than
     * the left */
    lookhere = ceil ((range[0] + range[1]) / 2.);
  }

  if (!receive) {
    /* Calculate the last local tree that we need to send to send_first */
    /* Set it to the last tree on send_first */
    ret = t8_offset_last (*send_first, offset_to) - cmesh_from->first_tree;
    /* If there are actually more trees on send_first than we have, we need to send
     * all our local trees to send_first */
    ret = SC_MIN (ret, cmesh_from->num_local_trees - 1);
    if (cmesh_from->mpirank != *send_first &&
        t8_offset_in_range (t8_offset_last (cmesh_from->mpirank, offset_from),
                            *send_first, offset_from)
        && ret == cmesh_from->num_local_trees - 1) {
      /* Substract one if our last tree already belonged to send_first,
       * and we counted this last tree. */
      ret--;
    }
  }
}

t8_debugf ("%s_first = %i, %s_last = %i, last_tree = %li\n",
           receive ? "recv" : "send", *send_first,
           receive ? "recv" : "send", *send_last, ret);

T8_ASSERT (*send_first >= 0);
//TODO:reactivate  T8_ASSERT (*send_last >= 0);
T8_ASSERT (receive || (ret >= 0 && ret < cmesh_from->num_local_trees));
T8_ASSERT (receive || ret == (t8_locidx_t) ret);
return (t8_locidx_t) ret;
}
#endif

/* TODO: deprecated, can be removed */
#if 0
/* Compute first and last process to which we will send data */
/* Returns the local tree_id of last local tree on send_first. */
static              t8_locidx_t
t8_cmesh_partition_sendrange_old (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                                  t8_gloidx_t * tree_offsets,
                                  int *send_first, int *send_last,
                                  int receive)
{
  t8_gloidx_t         first_tree = 0, last_tree, last_local_tree,
    first_local_tree, ret;
  int                 range[2], lookhere;
  int                 was_empty;

  /* p_i is first process we send to if and only if
   *
   * tree_first[p_i]    <=    tree_first    <=    tree_last[p_i]
   *      |                       ||                    |
   *  new_partition     cmesh_from->first_tree     new_partition
   */
  /* p_i is last process we send to if and only if
   *
   * tree_first[p_i]    <=    tree_last    <=    tree_last[p_i]
   *        |                     ||                    |
   *  new_partition     cmesh_from->first_tree     new_partition
   *                       + num_local_trees-1
   */
  /* We have to be careful with empty processes and first tree shared */
  /* If the first tree of a process is shared than it will not send this tree,
     but the smallest process that owns this tree will do. */

  t8_debugf ("Determining %srange\n", receive ? "receive" : "send");

  T8_ASSERT (tree_offsets != NULL);

  if ((!receive && !cmesh_from->set_partition)
      || (receive && !cmesh->set_partition)) {
    /* If cmesh_from is not partitioned we can send/recv all information
     * to ourself */
    *send_first = *send_last = cmesh_from->mpirank;
    return t8_offset_last (cmesh_from->mpirank, tree_offsets);
  }

  *send_first = -1;
  *send_last = -1;
  range[0] = 0;
  range[1] = cmesh->mpisize - 1;
  last_local_tree = cmesh_from->first_tree + cmesh_from->num_local_trees - 1;
  /* We exclude the first tree rom sending if it is shared */
  first_local_tree = !cmesh_from->first_tree_shared ?
    cmesh_from->first_tree : cmesh_from->first_tree + !receive;

  if (cmesh_from->num_local_trees == 0) {
    /* This partition is empty and we can not send anything */
    *send_last = -2;
    return -1;
  }

  /* This variable is true if we reveive and the cmesh we reveive from
   * had no local trees */
  was_empty = receive && cmesh->num_local_trees <= 0;

#if 0
  if (!receive
      && t8_offset_in_range (cmesh_from->first_tree, cmesh->mpirank,
                             tree_offsets)) {
    /* We keep the first tree, thus we are the first process we send to */
    *send_first = cmesh->mpirank;
  }
  else
#endif
  if (receive && cmesh->first_tree == cmesh_from->first_tree) {
    /* We receive the first tree from ourselves.
     * If this is our smallest local tree there can't be any smaller rank
     * sending to us. */
    if (!was_empty) {
      *send_first = cmesh->mpirank;
    }
  }

  /* Determine send_first */
  while (*send_first == -1) {
    lookhere = (range[0] + range[1]) / 2;
    /* first tree stores new first_tree of process lookhere */
    /* last tree the new last tree of process lookhere */
    first_tree = t8_offset_first (lookhere, tree_offsets);
    if (lookhere == cmesh->mpirank) {
      first_local_tree = cmesh_from->first_tree;
    }
    else {
      first_local_tree =
        cmesh_from->first_tree + cmesh_from->first_tree_shared;
    }
    /* If this first tree is shared we exclude it when we receive */
    if (receive && tree_offsets[lookhere] < 0) {
      first_tree++;
    }
    last_tree = t8_offset_last (lookhere, tree_offsets);

    if (first_tree == first_local_tree) {
      *send_first = lookhere;
      while (last_tree >= first_local_tree && lookhere >= 0) {
        lookhere--;
        last_tree = t8_offset_last (lookhere, tree_offsets);
      }
      lookhere++;
      /* Now we found the correct process, however the proc found may be
       * the first one and empty */
      while (t8_offset_empty (lookhere, tree_offsets)) {
        lookhere++;
      }
      *send_first = lookhere;
      /* We are done now */
    }
    else if (first_tree < first_local_tree && first_local_tree <= last_tree) {
      /* This has to be the correct process */
      *send_first = lookhere;
    }
    else if (last_tree < first_local_tree) {
      /* We have to look further right */
      range[0] = lookhere + 1;
    }
    else {
      T8_ASSERT (first_local_tree < first_tree);
      /* We have to look further left */
      range[1] = lookhere - 1;
    }
  }

  t8_debugf ("Done send first %i\n", *send_first);

  /* Determine send_last */
  range[0] = *send_first;
  range[1] = cmesh->mpisize;
  while (*send_last == -1) {
    /* TODO: This is our initial guess. Maybe mpirank is a better choice? */
    lookhere = (range[0] + range[1]) / 2;
    /* first tree stores new first_tree of process lookhere */
    /* last tree the new last tree of process lookhere */
    first_tree = t8_offset_first (lookhere, tree_offsets);
    /* If the first tree of lookhere is shared and we determine receive range
     * then we ignore this first tree */
    if (receive && tree_offsets[lookhere] < 0) {
      first_tree++;
    }

    last_tree = t8_offset_last (lookhere, tree_offsets);

    if (last_tree == last_local_tree) {
      while (first_tree <= last_local_tree && lookhere < cmesh->mpisize) {
        lookhere++;
        first_tree = t8_offset_first (lookhere, tree_offsets);
        /* If the first tree of lookhere is shared and we determine receive range
         * then we ignore this first tree */
        if (receive && tree_offsets[lookhere] < 0) {
          first_tree++;
        }
      }
      lookhere--;
      /* We have found our proc, but it may be the last and empty */
      /* It could also be that it has only one tree and this is shared,  then
       * we also do not consider it when we receive */
      while (t8_offset_empty (lookhere, tree_offsets)
             || (receive
                 && t8_offset_num_trees (lookhere, tree_offsets) == 1
                 && tree_offsets[lookhere] < 0)) {
        lookhere--;
      }
      *send_last = lookhere;
    }
    else if (first_tree <= last_local_tree && last_local_tree < last_tree) {
      *send_last = lookhere;
    }
    else if (last_tree < last_local_tree) {
      /* We have to look further right */
      range[0] = lookhere + 1;
    }
    else {
      T8_ASSERT (last_local_tree < first_tree);
      /* We have to look further left */
      range[1] = lookhere - 1;
    }
  }
  t8_debugf ("prelim sendlast = %i\n", *send_last);
  while (!receive && *send_last != cmesh->mpirank && *send_last >= 0
         && t8_offset_in_range (last_local_tree, *send_last, tree_offsets)
         && (last_local_tree == t8_offset_first (*send_last, tree_offsets)
             || last_local_tree == first_local_tree)) {
    /* The last process we want to send to only needs last_local_tree
     * but in fact it already has it, so we exclude it from the list.
     * We repeat this until we find a nonempty process that does not own the
     * last tree yet. */
    /* TODO: Will this cause trouble with the example
     *   Old    |0|1|-2| 4|4|5|
     *
     *   New    |0|1|-2|-2|4|5|
     *      since now both p1 and p2 want to send to p3?
     */
    (*send_last)--;
    while (*send_last >= 0 && t8_offset_empty (*send_last, tree_offsets)) {
      (*send_last)--;
    }
  }
  if (receive && t8_offset_in_range (cmesh->first_tree, cmesh->mpirank,
                                     tree_offsets)) {
    T8_ASSERT (*send_first <= cmesh->mpirank);
    if (cmesh->mpirank > *send_last) {
      *send_last = cmesh->mpirank;
    }
  }

  ret = t8_glo_min (t8_offset_last (*send_first, tree_offsets) -
                    cmesh_from->first_tree, cmesh_from->num_local_trees - 1);
  t8_debugf ("%s_first = %i, %s_last = %i, first_tree = %li\n",
             receive ? "recv" : "send", *send_first,
             receive ? "recv" : "send", *send_last, ret);

  T8_ASSERT (*send_first >= 0);
  T8_ASSERT (*send_last >= 0);
  T8_ASSERT (receive || (ret >= 0 && ret <= cmesh_from->num_local_trees));
  T8_ASSERT (ret == (t8_locidx_t) ret);
  return (t8_locidx_t) ret;
}
#endif

/* TODO: deprecated. Can be removed */
#if 0
/* Compute the first and last process from which we will receive local trees */
static void
t8_cmesh_partition_recvrange (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                              int *recv_first, int *recv_last)
{
  /*  */

  if (!cmesh_from->set_partition) {
    /* If cmesh_from is not partitioned we can receive all information
     * from ourself */
    *recv_first = cmesh_from->mpirank;
    *recv_last = cmesh_from->mpirank;
    return;
  }

  (void) t8_cmesh_partition_sendrange (cmesh, cmesh_from,
                                       recv_first, recv_last, 1);
}
#endif

/* A much faster version to compute the sendrange */
static              t8_locidx_t
t8_cmesh_partition_alternative_sendrange (t8_cmesh_t cmesh,
                                          t8_cmesh_t cmesh_from,
                                          int *send_first, int *send_last)
{
  t8_gloidx_t         first_tree = t8_cmesh_get_first_treeid (cmesh_from),
    last_tree;
  t8_gloidx_t         ret;
  t8_gloidx_t        *offset_to;
  t8_gloidx_t        *offset_from;
  int                 alternative_sendfirst, alternative_sendlast = -2;
  int                 flag, count;
  int                 some_owner = -1;  /* Passes as argument to first/last owner functions */

  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (!cmesh_from->set_partition) {
    /* If cmesh_from is not partitioned we can send/recv all information
     * to/from ourself */
    *send_first = *send_last = cmesh_from->mpirank;
    offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
    return t8_offset_last (cmesh_from->mpirank, offset_to);
  }

  offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);

  /* TODO: try to work around calling this function, since it has log(P) runtime */
  if (t8_offset_nosend (cmesh->mpirank, cmesh->mpisize,
                        offset_from, offset_to)) {
    /* We do not send at all */
    *send_last = -2;
    *send_first = -1;
    return -1;
  }
  flag = 0;
  if (cmesh_from->first_tree_shared == 1) {
    /* If the first tree was shared then we only send it, if
     * we own it in the new partition. In this case we send to ourself first. */
    if (t8_offset_in_range (first_tree, cmesh->mpirank, offset_to)) {
      alternative_sendfirst = cmesh->mpirank;
      flag = 1;                 /* We found the process */
      t8_debugf ("[H] I send to myself first.\n");
    }
    else {
      /* We do not own the first tree in the new partition */
      if (cmesh_from->num_local_trees == 1) {
        /* We do not send at all */
        *send_last = -2;
        *send_first = -1;
        return -1;
      }
      /* The second local tree of cmesh_from is the first that is sent */
      first_tree += 1;
    }
  }
  /* If we have not found the process yet, our first tree is not shared
   * or we consider our the second local tree */
  if (flag != 1) {
    /* search for smallest new owner that did not own this tree before */
    /* Get the smallest owner of first_tree in the new partition */
    alternative_sendfirst =
      t8_offset_first_owner_of_tree (cmesh->mpisize, first_tree, offset_to,
                                     &some_owner);
    while (alternative_sendfirst >= 0 && flag == 0) {
      if (alternative_sendfirst == cmesh->mpirank
          || t8_offset_empty (alternative_sendfirst, offset_from)
          || t8_offset_first (alternative_sendfirst,
                              offset_from) != first_tree) {
        /* We found the process if it is either ourself or it did not own first_tree
         * before */
        flag = 1;
      }
      else {
        /* Compute the next bigger process that has first_tree as local tree in
         * the new partition. */
        alternative_sendfirst =
          t8_offset_next_owner_of_tree (cmesh->mpisize, first_tree,
                                        offset_to, alternative_sendfirst);
      }
#if 0
      /* If the first tree is shared, we consider the second tree here
       * and if we did not find an owner in this round,
       * then we do not send any trees to any processes. */
      if (cmesh_from->first_tree_shared && flag == 0) {
        t8_debugf ("[H] New no send\n");
        *send_last = -2;
        *send_first = -1;
        return -1;
      }
      /* TODO: deprecated */
      /* If the first tree is shared, we consider the  second tree here and
       * must have found the process in the first round */
      T8_ASSERT (!cmesh_from->first_tree_shared || flag == 1);
#endif
    }
    T8_ASSERT (flag == 1);      /* We must have found the process by now */
    t8_debugf ("[H] Alternative send first = %i\n", alternative_sendfirst);
  }
  /* Get the last local tree on cmesh_from */
  last_tree = t8_cmesh_get_first_treeid (cmesh_from) +
    t8_cmesh_get_num_local_trees (cmesh_from) - 1;
  flag = 0;
  count = 0;
  if (last_tree == t8_cmesh_get_first_treeid (cmesh_from)
      && cmesh_from->first_tree_shared) {
    /* If our last tree is the first tree and it is shared then we only send
     * to ourselves */
    alternative_sendlast = alternative_sendfirst;
    flag = 1;
    T8_ASSERT (alternative_sendfirst == cmesh->mpirank);
  }
  else {
    /* Check the last tree and maybe the second last tree.
     * If we do not send our last tree, we have to check the second last one.
     * This while loop is executed once for the last tree and, if unsuccessfull
     * once for the second last tree. */
    while (last_tree >= first_tree && count < 2 && flag == 0) {
      count++;
      flag = 0;
      some_owner = -1;          /* Reset some_owner, since we do not know an owner of last_tree */
      /* Parse the new owners from the top and stop at the first process
       * that did not own last_tree */ ;
      alternative_sendlast =
        t8_offset_last_owner_of_tree (cmesh->mpisize, last_tree, offset_to,
                                      &some_owner);
      t8_debugf ("[H] Last owner of %lli is %i\n", (long long) last_tree,
                 alternative_sendlast);
      while (alternative_sendlast >= 0 && flag == 0) {
        if (alternative_sendlast == cmesh->mpirank ||
            t8_offset_empty (alternative_sendlast, offset_from) ||
            t8_offset_first (alternative_sendlast,
                             offset_from) != last_tree) {
          /* alternative_sendlast is either this process or it did not have
           * last_tree before */
          flag = 1;
        }
        else {
          /* Compute next smaller process that has last_tree as local tree */
          alternative_sendlast =
            t8_offset_prev_owner_of_tree (cmesh->mpisize, last_tree,
                                          offset_to, alternative_sendlast);
          t8_debugf ("[H] Prev owner of %lli is %i\n", (long long) last_tree,
                     alternative_sendlast);
        }
      }
      /* If we did not found the alternative send here, then all procs
       * owned the last tree before and we have to check the second last tree.
       * Here a success is guaranteed. */
      last_tree--;
      /* If this second last tree is the first tree and the first tree is shared,
       * then we do not send to any other processes than ourselves. */
      if (flag == 0 && last_tree == t8_cmesh_get_first_treeid (cmesh_from)
          && cmesh_from->first_tree_shared) {
        alternative_sendlast = alternative_sendfirst;
        T8_ASSERT (alternative_sendfirst == cmesh->mpirank);
        flag = 1;
      }
    }
  }
  if (flag == 0) {
    /* This should never happen */
    t8_debugf ("[H] Could not find alternative sendlast.\n");
    T8_ASSERT (0 == 1);
  }
  else {
    t8_debugf ("[H] Found alternative last send %i\n", alternative_sendlast);
  }
  *send_first = alternative_sendfirst;
  *send_last = alternative_sendlast;

  /* Calculate the last local tree that we need to send to send_first */
  /* Set it to the last tree on send_first */
  ret = t8_offset_last (*send_first, offset_to) - cmesh_from->first_tree;
  /* If there are actually more trees on send_first than we have, we need to send
   * all our local trees to send_first */
  ret = SC_MIN (ret, cmesh_from->num_local_trees - 1);
  if (cmesh_from->mpirank != *send_first &&
      t8_offset_in_range (t8_offset_last (cmesh_from->mpirank, offset_from),
                          *send_first, offset_from)
      && ret == cmesh_from->num_local_trees - 1) {
    /* Substract one if our last tree already belonged to send_first,
     * and we counted this last tree. */
    ret--;
  }

  t8_debugf ("%s_first = %i, %s_last = %i, last_tree = %li\n",
             "send", *send_first, "send", *send_last, ret);

  T8_ASSERT (*send_first >= 0);
//TODO:reactivate  T8_ASSERT (*send_last >= 0);
  T8_ASSERT ((ret >= 0 && ret < cmesh_from->num_local_trees));
  T8_ASSERT (ret == (t8_locidx_t) ret);
  return (t8_locidx_t) ret;
}

/* A much faster version to compute the receive range */
static void
t8_cmesh_partition_alternative_recvrange (t8_cmesh_t cmesh,
                                          t8_cmesh_t cmesh_from,
                                          int *recv_first, int *recv_last)
{
  t8_gloidx_t         first_tree, last_tree;
  t8_gloidx_t        *offset_to;
  t8_gloidx_t        *offset_from;
  int                 alternative_recvfirst, alternative_recvlast;
  int                 some_owner = -1;  /* Passes as argument to first/last owner functions */

  if (cmesh->num_local_trees == 0) {
    /* This process is empty and therefore does not receive anything */
    *recv_last = -2;
    *recv_first = -1;
    return;
  }
  if (!cmesh_from->set_partition) {
    /* If cmesh_from is not partitioned we can receive all information
     * from ourself */
    *recv_first = cmesh_from->mpirank;
    *recv_last = cmesh_from->mpirank;
    return;
  }
  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  /* Get the new first local tree */
  first_tree = t8_cmesh_get_first_treeid (cmesh);
  if (t8_offset_in_range (first_tree, cmesh->mpirank, offset_from)) {
    /* It it already was a local tree then we received it from ourselves
     * and are thus the first process we receive from */
    alternative_recvfirst = cmesh->mpirank;
  }
  else {
    /* Otherwise the first process we receive from is the smallest process that
     * had our new first tree as a local tree. */
    some_owner = -1;
    alternative_recvfirst =
      t8_offset_first_owner_of_tree (cmesh->mpisize, first_tree, offset_from,
                                     &some_owner);
  }
  /* Get the new last local tree */
  last_tree = t8_offset_last (cmesh->mpirank, offset_to);
  if (t8_offset_in_range (last_tree, cmesh->mpirank, offset_from)) {
    /* We had our last local tree as a local tree before and thus
     * we are the last process that we receive from */
    alternative_recvlast = cmesh->mpirank;
  }
  else {
    /* We receive from the smallest process that had our new last local tree
     * as a local tree. */
    if (first_tree != last_tree) {
      /* We can reuse the computed owner if first_tree == last_tree,
       * otherwise we have to reset it */
      some_owner = -1;
    }
    alternative_recvlast = t8_offset_first_owner_of_tree (cmesh->mpisize,
                                                          last_tree,
                                                          offset_from,
                                                          &some_owner);
  }
  *recv_first = alternative_recvfirst;
  *recv_last = alternative_recvlast;
}

/* Compute the number of bytes that need to be allocated in the send buffer
 * for the neighbor entries of ghost */
static              size_t
t8_partition_compute_gnb (t8_cmesh_t cmesh_from, sc_array_t * send_as_ghost)
{
  size_t              ghost_neighbor_bytes = 0, ighost;
  t8_locidx_t         ghost_id;
  t8_eclass_t         eclass;

  for (ighost = 0; ighost < send_as_ghost->elem_count; ighost++) {
    ghost_id = *(t8_locidx_t *) sc_array_index (send_as_ghost, ighost);
    T8_ASSERT (ghost_id >= 0);
    if (ghost_id < cmesh_from->num_local_trees) {
      /* The ghost to send is currently a tree */
      eclass = t8_cmesh_get_tree_class (cmesh_from, ghost_id);
    }
    else {
      /* The ghost to send is currently a ghost */
      T8_ASSERT (ghost_id < cmesh_from->num_local_trees
                 + cmesh_from->num_ghosts);
      eclass = t8_cmesh_get_ghost_class (cmesh_from, ghost_id
                                         - cmesh_from->num_local_trees);
    }
    ghost_neighbor_bytes += t8_eclass_num_faces[eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))      /* offset */
      +T8_ADD_PADDING (t8_eclass_num_faces[eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))); /* padding */
  }
  return ghost_neighbor_bytes;
}

/* Compute the number of bytes that need to be allocated in the send buffer
 * for the attribute entries of all ghosts. */
static size_t
t8_partition_compute_gab (t8_cmesh_t cmesh_from, sc_array_t * send_as_ghost,
                          size_t * attr_info_bytes)
{
  size_t              ghost_attribute_bytes = 0, ighost;
  t8_locidx_t         ghost_id, ghost_id_min_offset;
  t8_cghost_t         ghost;
  t8_ctree_t          tree;

  T8_ASSERT (attr_info_bytes != NULL);
  *attr_info_bytes = 0;
  for (ighost = 0; ighost < send_as_ghost->elem_count; ighost++) {
    ghost_id = *(t8_locidx_t *) sc_array_index (send_as_ghost, ighost);
    t8_debugf ("[H] reading ghost_id %i\n", (int) ghost_id);
    T8_ASSERT (ghost_id >= 0);
    if (ghost_id < cmesh_from->num_local_trees) {
      /* This ghost is currently a local tree */
      tree = t8_cmesh_get_tree (cmesh_from, ghost_id);
      ghost_attribute_bytes += t8_cmesh_trees_attribute_size (tree);
      *attr_info_bytes += tree->num_attributes
        * sizeof (t8_attribute_info_struct_t);
    }
    else {
      /* This ghost is currently a ghost */
      ghost_id_min_offset =
        ghost_id - t8_cmesh_get_num_local_trees (cmesh_from);
      t8_debugf ("[H] ghost_id %i is ghost number %i\n", ghost_id,
                 ghost_id_min_offset);
      T8_ASSERT (0 <= ghost_id_min_offset
                 && ghost_id_min_offset <
                 t8_cmesh_get_num_ghosts (cmesh_from));
      ghost =
        t8_cmesh_trees_get_ghost (cmesh_from->trees, ghost_id_min_offset);
      ghost_attribute_bytes += t8_cmesh_trees_ghost_attribute_size (ghost);
      *attr_info_bytes += ghost->num_attributes
        * sizeof (t8_attribute_info_struct_t);
    }
  }
  return ghost_attribute_bytes;
}

/* Determine whether a local tree or ghost should be send to process p as a ghost.
 * This is the case if and only if:
 *  - tree will not be a local tree on p
 * and
 *  - If p sends to itself then tree is currently neither
 *    a local tree nor a ghost on p
 * and
 *  - we are the smallest rank under all procs sending a tree to p that
 *    has this tree as ghost or local tree . */
static int
t8_cmesh_send_ghost (t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from,
                     int p, t8_locidx_t tree)
{
  t8_gloidx_t         tree_id, *ghost_neighbors, neighbor;
  t8_gloidx_t        *from_offsets;
  t8_locidx_t        *tree_neighbors;
  t8_cghost_t         ghost = NULL;
  t8_ctree_t          ctree = NULL;
  t8_eclass_t         eclass;
  int                 proc, iface;
  size_t              left = 0, right = cmesh->mpirank;
  int                 found = 0;
  t8_gloidx_t        *offset_to;

  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);

  if (cmesh_from->set_partition) {
    T8_ASSERT (cmesh_from->tree_offsets != NULL);
    from_offsets = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  }
  else {
    from_offsets = NULL;
  }
  if (tree < cmesh_from->num_local_trees) {
    /* Given local id belongs to a tree. We compute its global id */
    ctree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, tree,
                                         &tree_neighbors, NULL);
    tree_id = tree + cmesh_from->first_tree;
    eclass = ctree->eclass;
  }
  else {
    /* Given local id belongs to a ghost. We store the ghost and its
     * global id */
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees,
                                          tree - cmesh_from->num_local_trees,
                                          &ghost_neighbors, NULL);
    tree_id = ghost->treeid;
    eclass = ghost->eclass;
  }

  if (t8_offset_in_range (tree_id, p, offset_to)) {
    /* The tree/ghost will be a local tree on p */
    return 0;
  }

  /* Parse all processes that own neighbors of tree/ghost and check for each if
   * it sends this neighbor to p. If the current rank is the smallest such rank,
   * then we send the tree/ghost. */
  for (iface = 0; iface < t8_eclass_num_faces[eclass]; iface++) {
    /* Get the global id of the considered neighbor */
    neighbor = ctree != NULL ?
      t8_cmesh_get_global_id ((t8_cmesh_t) cmesh_from,
                              tree_neighbors[iface]) : ghost_neighbors[iface];
    if (neighbor == tree_id) {
      /* There is no neighbor at this face */
      continue;
    }
    if (!t8_offset_in_range (neighbor, p, offset_to)) {
      /* This neighbor will not be send to p */
      continue;
    }
    /* If the receiving rank is the sending rank, we definetely send the ghost */
    if (p == cmesh_from->mpirank) {
      return 1;
    }
    /* If the receiving rank did own neighbor in cmesh_from, then it will
     * send it to itself and thus also this ghost. */
    if (t8_offset_in_range (neighbor, p, from_offsets)) {
      return 0;
    }
    /* We perform a binary search in the offset array to find the
     * process proc that owns the neighbor iface.
     * We start with checking mpirank, since it is most probable.
     */
    found = 0;
    proc = cmesh->mpirank;
    if (!cmesh_from->set_partition) {
      /* If the original cmesh is replicated then we already own
       * all trees */
      found = 1;
    }
    left = 0, right = cmesh->mpisize;
    while (!found) {
      if (t8_offset_first (proc, from_offsets) + (from_offsets[proc] < 0)
          > neighbor) {
        /* look left */
        right = proc;
        proc = left + (proc - left) / 2;
      }
      else if (t8_offset_last (proc, from_offsets)
               < neighbor) {
        /* look right */
        left = proc;
        proc++;
        proc += (right - proc) / 2;
      }
      else {
        found = 1;
      }
    }
    if (cmesh_from->set_partition) {
      /* We have to consider the special case where neighbor is owned by more
       * than one process. We find the smallest one of them. */
      while (neighbor == t8_offset_first (proc, from_offsets)
             && from_offsets[proc] < 0) {
        proc--;
        while (t8_offset_empty (proc, from_offsets)) {
          /* Skip empty processes */
          proc--;
        }
      }
    }
    T8_ASSERT (0 <= proc && proc < cmesh->mpisize);
    T8_ASSERT (t8_offset_in_range (neighbor, proc, from_offsets));
    if (proc < cmesh->mpirank && t8_offset_sendstree (proc, p, neighbor,
                                                      from_offsets,
                                                      offset_to)) {
      /* We do not send the ghost */
      return 0;
    }
  }
  /* Under all processes that have our tree/ghost as ghost and send a neighbor
   * of it to p there was none that
   * has smaller rank than myrank and sends a neighbor to p */
  return 1;
}

/* copy all tree/ghost/attribute data to the send buffer */
static void
t8_cmesh_partition_copy_data (char *send_buffer, t8_cmesh_t cmesh,
                              const struct t8_cmesh *cmesh_from,
                              t8_locidx_t num_trees, size_t attr_info_bytes,
                              size_t ghost_attr_info_bytes,
                              size_t ghost_neighbor_bytes,
                              size_t tree_neighbor_bytes,
                              size_t tree_attribute_bytes,
                              sc_array_t * send_as_ghost,
                              t8_locidx_t send_first, t8_locidx_t send_last,
                              size_t total_alloc, int to_proc)
{
  t8_ctree_t          tree, tree_cpy;
  int                 num_attributes;
  size_t              temp_offset_tree, temp_offset_att, iz, temp_offset,
    temp_offset_data, last_offset, last_num_att, last_size,
    temp_offset_ghost_att, temp_offset_ghost_data, temp_offset_ghost,
    ghost_attr_info_bytes_sofar;
  size_t              ghost_att_size;
  //ssize_t             last_attribute_diff;
  t8_attribute_info_struct_t *attr_info;
  void               *first_attribute;
  t8_locidx_t         num_ghost_send = send_as_ghost->elem_count;
  t8_locidx_t         ghosts_left;
  t8_locidx_t        *face_neighbor, ghost_id, itree;
  t8_gloidx_t        *face_neighbor_g, *face_neighbor_gnew, new_neighbor;
  t8_cghost_t         ghost, ghost_cpy;
  int                 iface, iatt;
  int8_t             *ttf_ghost, *ttf;

  /* Copy all trees to the send buffer */
  /* TODO: This is currently inefficient since we copy each tree for itself.
   *       Best practive is to copy chunks of trees out of the different part
   *       arrays of cmesh_from */
  if (total_alloc == 0 || send_buffer == NULL) {
    t8_debugf ("No data to store in buffer.\n");
    return;
  }
  /* offset from the last ghost_faces entry to the current tree_faces entry */
  temp_offset = 0;
  /* offset from the last tree_faces entry to the current tree_att_info entry */
  temp_offset_att = 0;
  /* offset from the last tree_att_info entry to the current attribute's data */
  temp_offset_data = 0;
  /* offset from the beginning to the current tree */
  temp_offset_tree = 0;
  for (itree = send_first; itree <= send_last; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree,
                                        &face_neighbor, NULL);

    (void) memcpy (send_buffer + temp_offset_tree, tree,
                   sizeof (t8_ctree_struct_t));
    temp_offset_tree += sizeof (t8_ctree_struct_t);
    /* Copy all face neighbor information to send_buffer */
    (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                   num_ghost_send * sizeof (t8_cghost_struct_t) +
                   ghost_neighbor_bytes + temp_offset, face_neighbor,
                   t8_eclass_num_faces[tree->eclass] *
                   (sizeof (t8_locidx_t) + sizeof (int8_t)));
    temp_offset += t8_eclass_num_faces[tree->eclass] *
      (sizeof (t8_locidx_t) + sizeof (int8_t))
      + T8_ADD_PADDING (t8_eclass_num_faces[tree->eclass] *
                        (sizeof (t8_locidx_t) + sizeof (int8_t)));
    /* TODO: ??? temp_offset += T8_ADD_PADDING (temp_offset) instead of the last 2 lines? */
    if (tree->num_attributes > 0) {
      /* Copy all attribute infos to send_buffer */
      (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                     num_ghost_send * sizeof (t8_cghost_struct_t) +
                     ghost_neighbor_bytes + tree_neighbor_bytes +
                     temp_offset_att, T8_TREE_ATTR_INFO (tree, 0),
                     tree->num_attributes *
                     sizeof (t8_attribute_info_struct_t));
      temp_offset_att +=
        tree->num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute data to send_buffer */
      (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                     num_ghost_send * sizeof (t8_cghost_struct_t) +
                     ghost_neighbor_bytes + tree_neighbor_bytes +
                     attr_info_bytes + temp_offset_data,
                     T8_TREE_ATTR (tree, T8_TREE_ATTR_INFO (tree, 0)),
                     t8_cmesh_trees_attribute_size (tree));
      temp_offset_data += t8_cmesh_trees_attribute_size (tree);
    }
  }
  T8_ASSERT (tree_attribute_bytes == temp_offset_data);
  /* Set new face_neighbor offsets */
  /* TODO: indent bug? */
  /* Computes the offset of the face neighbors of the new trees */
  temp_offset = num_trees * sizeof (t8_ctree_struct_t) +
    num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes;
  /* Compute the offset of the new attribute infos */
  temp_offset_att = temp_offset + tree_neighbor_bytes;
  temp_offset_tree = 0;

#if 0
  /* This has to be add to each attribute info offset to calculated the new offset */
  last_attribute_diff = attr_info_bytes - attr_info->attribute_offset;
#endif
  /* Set attribute offsets of trees and attribute data offsets of info objects */
  tree_cpy = NULL;
  last_num_att = 0;
  last_size = 0;
  last_offset = attr_info_bytes;

  for (itree = send_first; itree <= send_last; itree++) {
    /* Get the current tree */
    tree_cpy = (t8_ctree_t) (send_buffer + temp_offset_tree);

    /* new neighbor offset of tree */
    tree_cpy->neigh_offset = temp_offset - temp_offset_tree;

    /* Set new face neighbor entries, since we store local ids we have to adapt
     * to the local ids of the new process */
    face_neighbor = (t8_locidx_t *) T8_TREE_FACE (tree_cpy);
    for (iface = 0; iface < t8_eclass_num_faces[tree_cpy->eclass]; iface++) {
      t8_cmesh_partition_send_change_neighbor (cmesh, (t8_cmesh_t) cmesh_from,
                                               face_neighbor + iface,
                                               to_proc);
    }

    /* compute neighbor offset for next tree */
    temp_offset += t8_eclass_num_faces[tree_cpy->eclass] *
      (sizeof (t8_locidx_t) + sizeof (int8_t))
      + T8_ADD_PADDING (t8_eclass_num_faces[tree_cpy->eclass] *
                        (sizeof (t8_locidx_t) + sizeof (int8_t)));

    /* new attribute offset for tree */
    tree_cpy->att_offset = temp_offset_att - temp_offset_tree;
    if (tree_cpy->num_attributes > 0) {
      attr_info = T8_TREE_ATTR_INFO (tree_cpy, 0);
      attr_info->attribute_offset = last_offset + last_size -
        last_num_att * sizeof (t8_attribute_info_struct_t);
      last_offset = attr_info->attribute_offset;
      last_size = attr_info->attribute_size;

      /* set new attribtue data offsets */
      for (iz = 1; iz < (size_t) tree_cpy->num_attributes; iz++) {
        attr_info->attribute_offset = last_offset + last_size;
        last_size = attr_info->attribute_size;
        attr_info++;
      }
      temp_offset_att += tree_cpy->num_attributes *
        sizeof (t8_attribute_info_struct_t);
    }
    temp_offset_tree += sizeof (t8_ctree_struct_t);
    last_num_att = tree_cpy->num_attributes;
  }

  /* Copy all ghosts and set their face entries and offsets */
  /* Offset of ghost face_neighbor from first ghost */
  temp_offset = num_ghost_send * sizeof (t8_cghost_struct_t);
  /* Offset of current ghost from first ghost */
  temp_offset_ghost = 0;
  /* offset from the last tree_attribute entry to the current tree_att_info entry
   * for the ghost */
  temp_offset_ghost_att = 0;
  /* offset from the last tree_att_info entry to the current attribute */
  temp_offset_ghost_data = 0;
  /* number of bytes of attribute_infos that we already counted */
  ghost_attr_info_bytes_sofar = 0;
  for (iz = 0; iz < send_as_ghost->elem_count; iz++) {
    ghost_id = *((t8_locidx_t *) sc_array_index (send_as_ghost, iz));
    ghost_cpy = (t8_cghost_t) (send_buffer +
                               num_trees * sizeof (t8_ctree_struct_t)
                               + iz * sizeof (t8_cghost_struct_t));
    /* Get the correct element class from the ghost_id.
     * For this we have to check whether ghost_id points to a local tree or ghost */
    ghost_cpy->eclass = ghost_id < cmesh_from->num_local_trees ?
      t8_cmesh_get_tree_class ((t8_cmesh_t) cmesh_from, ghost_id) :
      t8_cmesh_get_ghost_class ((t8_cmesh_t) cmesh_from, ghost_id
                                - cmesh_from->num_local_trees);
    ghost_cpy->neigh_offset = temp_offset - temp_offset_ghost;  /* New face neighbor offset */
    face_neighbor_gnew = (t8_gloidx_t *) T8_GHOST_FACE (ghost_cpy);
    ttf_ghost = T8_GHOST_TTF (ghost_cpy);
    if (ghost_id >= cmesh_from->num_local_trees) {
      /* The ghost that we send was a local ghost */
      ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees,
                                            ghost_id -
                                            cmesh_from->num_local_trees,
                                            &face_neighbor_g, &ttf);
      tree = NULL;
      /* Set entries of the new ghost */
      ghost_cpy->eclass = ghost->eclass;
      ghost_cpy->treeid = ghost->treeid;
      ghost_cpy->num_attributes = ghost->num_attributes;
      num_attributes = ghost_cpy->num_attributes;
      if (num_attributes > 0) {
        /* get a pointer to the first att_info and the first attribute */
        attr_info = T8_GHOST_ATTR_INFO (ghost, 0);
        first_attribute = T8_GHOST_ATTR (ghost, attr_info);
        ghost_att_size = t8_cmesh_trees_ghost_attribute_size (ghost);
      }
      else {
        attr_info = NULL;
        first_attribute = NULL;
        ghost_att_size = 0;
      }
      /* Copy face_neighbor entries and ttf entries */
      memcpy (face_neighbor_gnew, face_neighbor_g,
              t8_eclass_num_faces[ghost_cpy->eclass] * (sizeof (t8_gloidx_t)
                                                        + sizeof (int8_t)));
    }
    else {
      /* The ghost we send was a local tree */
      T8_ASSERT (0 <= ghost_id && ghost_id < cmesh_from->num_local_trees);
      tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, ghost_id,
                                          &face_neighbor, &ttf);
      ghost = NULL;
      T8_ASSERT (ghost_id == tree->treeid);
      /* Set entries of the new ghost */
      ghost_cpy->eclass = tree->eclass;
      ghost_cpy->treeid = ghost_id + cmesh_from->first_tree;
      ghost_cpy->num_attributes = tree->num_attributes;
      num_attributes = ghost_cpy->num_attributes;
      if (num_attributes > 0) {
        /* get a pointer to the first att_info and the first attribute */
        attr_info = T8_TREE_ATTR_INFO (tree, 0);
        first_attribute = T8_TREE_ATTR (tree, attr_info);
        ghost_att_size = t8_cmesh_trees_attribute_size (tree);
      }
      else {
        attr_info = NULL;
        first_attribute = NULL;
        ghost_att_size = 0;
      }
      /* copy face_neighbor entries, since the ones on the tree are local and
       * we need global, we have to compute each one */
      for (iface = 0; iface < t8_eclass_num_faces[ghost_cpy->eclass]; iface++) {
        if (face_neighbor[iface] < 0) {
          /* TODO: think about this */
          new_neighbor = -1;    /* boundary indicator */
        }
        else {
          /* Compute global index from local index */
          new_neighbor = t8_cmesh_get_global_id ((t8_cmesh_t) cmesh_from,
                                                 face_neighbor[iface]);
        }
        face_neighbor_gnew[iface] = new_neighbor;
      }
      /* Copy tree_to_face entries */
      memcpy (ttf_ghost, ttf, t8_eclass_num_faces[ghost_cpy->eclass]
              * sizeof (int8_t));
    }                           /* Done distinction between from tree and from ghost */

    /* Compute and store new attribute offset of this ghost */
    ghosts_left = send_as_ghost->elem_count - iz;
    ghost_cpy->att_offset =
      ghosts_left * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes +
      tree_neighbor_bytes + attr_info_bytes + tree_attribute_bytes +
      temp_offset_ghost_att;
    if (num_attributes > 0) {
      size_t              this_ghosts_att_info_size;
      t8_attribute_info_struct_t *first_attr_info;

      t8_debugf ("[H] Copy %i atts of total size %lu\n",
                 num_attributes, ghost_att_size);

      /* The byte count of this ghosts attribute info structs */
      this_ghosts_att_info_size =
        num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute info data of this ghost */
      first_attr_info =
        (t8_attribute_info_struct_t *) T8_GHOST_FIRST_ATT (ghost_cpy);
      memcpy (first_attr_info, attr_info, this_ghosts_att_info_size);
      temp_offset_ghost_att += this_ghosts_att_info_size;

      /* Compute all new attribute data offsets */
      for (iatt = 0; iatt < num_attributes; iatt++) {
        /* Get the current attribute info */
        attr_info = T8_GHOST_ATTR_INFO (ghost_cpy, iatt);
        /* The new attribute offset is the offset from the first att_info to the data.
         * Thus, the count of the bytes occupied by the att_info (ghosts_attr_info_bytes)
         * plus the count of all attributes before this attribute (this_data_temp_offset).*/
        /* all att info from this ghost and after. + all attributes before this attribute  */
        attr_info->attribute_offset = ghost_attr_info_bytes -
          ghost_attr_info_bytes_sofar + temp_offset_ghost_data;
        temp_offset_ghost_data += attr_info->attribute_size;
      }
      ghost_attr_info_bytes_sofar +=
        num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute data of this ghost */
      memcpy (T8_GHOST_ATTR (ghost_cpy, first_attr_info),
              first_attribute, ghost_att_size);
    }                           /* end num_attributes > 0 */
    /* compute new offsets */
    temp_offset += t8_eclass_num_faces[ghost_cpy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))    /* offset */
      +T8_ADD_PADDING (t8_eclass_num_faces[ghost_cpy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t)));      /* padding */
    temp_offset_ghost += sizeof (t8_cghost_struct_t);
  }
  /* Store number of trees and ghosts at the end of send buffer */
  *(t8_locidx_t *) (send_buffer + total_alloc - 2 * sizeof (t8_locidx_t))
    = num_trees;
  *(t8_locidx_t *) (send_buffer + total_alloc - sizeof (t8_locidx_t))
    = num_ghost_send;
}

/* For each process that we send trees to, we loop over all these trees and
 * determine
 *  - The number of bytes to send for the neighbor entries,
 *  - the attribute entries,
 *  - the attribute info entries.
 *  - The local trees and ghosts that have to be send as a ghost.
 *  - The local trees that we will be a ghost in the new partition and
 *    therefore need to be kept local.
 */
static void
t8_cmesh_partition_sendtreeloop (t8_cmesh_t cmesh,
                                 const struct t8_cmesh *cmesh_from,
                                 t8_locidx_t range_start,
                                 t8_locidx_t range_end,
                                 size_t * tree_neighbor_bytes,
                                 size_t * attr_bytes,
                                 size_t * attr_info_bytes,
                                 int8_t * ghost_flag_send,
                                 t8_gloidx_t * tree_offset, int iproc,
                                 sc_array_t * send_as_ghost)
{
  t8_ctree_t          tree;
  t8_locidx_t         neighbor, *face_neighbor, itree;
  int8_t             *ttf;
  int                 iface;
  t8_gloidx_t        *offset_from, *offset_to;

  if (cmesh_from->set_partition) {
    offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  }
  else {
    offset_from = NULL;
  }
  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  /* loop over all trees that will be send */
  for (itree = range_start; itree <= range_end; itree++) {
    /* test if we really send the tree itree to the process iproc */
    T8_ASSERT ((!cmesh_from->set_partition && iproc == cmesh_from->mpirank) ||
               t8_offset_sendstree (cmesh_from->mpirank, iproc,
                                    itree + cmesh_from->first_tree,
                                    offset_from, offset_to));
    tree =
      t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree, &face_neighbor,
                                   &ttf);
    /* Count the additional memory needed per tree from neighbors */
    *tree_neighbor_bytes += t8_eclass_num_faces[tree->eclass] *
      (sizeof (*face_neighbor) + sizeof (*ttf));
    /* TODO: change padding to sizeof (void*) */
    *tree_neighbor_bytes += T8_ADD_PADDING (*tree_neighbor_bytes);      /* padding to make number of bytes per tree
                                                                           a multiple of 4 */
    /*  Compute number of attribute bytes in this tree range.
     *       Not every tree has an attribute */
    *attr_info_bytes += tree->num_attributes *
      sizeof (t8_attribute_info_struct_t);
    *attr_bytes += t8_cmesh_trees_attribute_size (tree);

    /* loop over all faces of each tree to determine ghost to send */
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      neighbor = face_neighbor[iface];
      if (neighbor >= 0 && neighbor != itree) { /* Consider only non-boundary neighbors */
        if ((neighbor < cmesh_from->num_local_trees &&
             (neighbor < range_start || neighbor > range_end))
            || neighbor >= cmesh_from->num_local_trees) {
          /* neighbor is a local tree or local ghost and will be ghost on iproc */
          if (ghost_flag_send[neighbor] == 0
              && t8_cmesh_send_ghost (cmesh, cmesh_from, iproc, neighbor)) {
            /* we did not add this neighbor yet and it should be send to iproc */
            ghost_flag_send[neighbor] = 1;
            *((t8_locidx_t *) sc_array_push (send_as_ghost)) = neighbor;
          }
        }
      }
    }
  }
}

/*********************************************/
/*        Send loop                          */
/*********************************************/
/* For each process we have to send to, we fill the send buffers and start
 * the communication */
/* Returns the number of processes we send data to */
static int
t8_cmesh_partition_sendloop (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                             int *num_request_alloc, int *send_first,
                             int *send_last,
                             char ***send_buffer, char **my_buffer,
                             size_t * my_buffer_bytes,
                             sc_MPI_Request ** requests, sc_MPI_Comm comm)
{
  size_t              attr_bytes = 0, tree_neighbor_bytes,
    ghost_neighbor_bytes, attr_info_bytes, ghost_attribute_bytes,
    ghost_attr_info_bytes;
  size_t              total_alloc;
  int                 iproc, flag;
  int                 mpiret, num_send_mpi = 0;
  char               *buffer;
  t8_locidx_t         num_trees, num_ghost_send, range_start, range_end;
  sc_array_t          send_as_ghost;    /* Stores local id's of trees and ghosts that will be send as ghosts */
  int8_t             *ghost_flag_send;  /* For each local tree and ghost set to 1 if it is in send_as_ghost */
  t8_gloidx_t        *offset_from, *offset_to;

  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  /* We use two different flag arrays here, since ghost_flag_send needs to
   * be reset for each process we send to, while ghost_flag_keep keeps is entries.
   * Otherwise we could have set bitflags and only use one array */
  ghost_flag_send = T8_ALLOC (int8_t, cmesh_from->num_local_trees +
                              cmesh_from->num_ghosts);

  sc_array_init (&send_as_ghost, sizeof (t8_locidx_t));

  range_end =
    t8_cmesh_partition_alternative_sendrange (cmesh,
                                              (t8_cmesh_t) cmesh_from,
                                              send_first, send_last);

  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (cmesh_from->set_partition) {
    range_start = cmesh_from->first_tree_shared;        /* Stores the first local tree that was not send yet */
    /* We do not send the first tree if it is shared, so we start with the second tree then */
    if (t8_offset_in_range
        (cmesh_from->first_tree, cmesh->mpirank, offset_to)) {
      /* In the special case that we send our first tree to ourselves, we
       * set it as start even if it is shared */
      range_start = 0;
    }
    offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  }
  else {
    offset_from = NULL;
    range_start = t8_offset_first (cmesh_from->mpirank, offset_to);
  }

  /* Determine the number of request for MPI communication.
   * If the current process is included in the sendrange we do not allocate
   * a request and buffer for it. */
  flag = cmesh->mpirank < *send_first || cmesh->mpirank > *send_last ? 1 : 0;
  *num_request_alloc = *send_last - *send_first + flag;
  /* range_end stores (my rank) local tree_id of last tree on *send_first */

  *send_buffer = T8_ALLOC (char *, *num_request_alloc);
  *requests = T8_ALLOC (sc_MPI_Request, *num_request_alloc);

  flag = 0;

  /* iproc is incremented below such that we skip those processes we do not
   * send to */
  for (iproc = *send_first; iproc <= *send_last;) {
    /* Since we do not allocate a request and buffer for my rank if it is included in
     * the sendrange, we have to offset the index into these array by -1 when the
     * current irpoc is bigger than my rank */
    if (*send_first <= cmesh->mpirank && iproc > cmesh->mpirank) {
      flag = 1;
    }
    while (cmesh_from->set_partition &&
           !t8_offset_sendstree (cmesh_from->mpirank, iproc,
                                 range_start + cmesh_from->first_tree,
                                 offset_from, offset_to)) {
      /* We skip trees that we do not send */
      t8_debugf ("[H] Skipping local tree %i\n", range_start);
      range_start++;
    }

    t8_debugf ("send to %i [%i,%i]\n", iproc, range_start, range_end);
    t8_debugf ("Start send loop\n");
    attr_bytes = 0;
    tree_neighbor_bytes = 0;
    attr_info_bytes = 0;
    ghost_neighbor_bytes = 0;
    ghost_attribute_bytes = 0;
    ghost_attr_info_bytes = 0;

    memset (ghost_flag_send, 0,
            (cmesh_from->num_local_trees + cmesh_from->num_ghosts)
            * sizeof (int8_t)); /* Yes, i know that sizeof(int8_t) is always 1, so what?! */
    sc_array_truncate (&send_as_ghost);
    /* loop over trees to calculate buffersize and which trees and ghosts to send */
    t8_cmesh_partition_sendtreeloop (cmesh, cmesh_from, range_start,
                                     range_end, &tree_neighbor_bytes,
                                     &attr_bytes,
                                     &attr_info_bytes, ghost_flag_send,
                                     offset_to, iproc, &send_as_ghost);
    /* loop over trees ends here */
    /* Calculate and allocate memory for send buffer */
      /**********************************************************/
    /*
     *      The data that we send has the layout
     *
     * | tree 0,1...| ghost 0,1...| face_neighbors0, tree_to_face0,...| ghost_neighbors 0,1,...| Attinfo00,01,...,Attinfoend | attdata00,01,... |
     *
     */
#if 0
    /* TODO: deprecated , remove? */
    attr_info_bytes += sizeof (t8_attribute_info_struct_t);     /* We always put one attribute info structat the end */
#endif
    num_trees = SC_MAX (range_end - range_start + 1, 0);

#if 0
    /* Output of all ghosts that are send */
    {
      size_t              count;
      t8_debugf ("Send as ghosts (global):\n");
      for (count = 0; count < send_as_ghost.elem_count; count++) {
        t8_debugf ("%li\n", t8_cmesh_get_global_id (cmesh_from,
                                                    *(t8_locidx_t
                                                      *) (sc_array_index_int
                                                          (&send_as_ghost,
                                                           count))));
      }
    }
#endif

    num_ghost_send = send_as_ghost.elem_count;
    /* parse through send_as_ghost to compute ghost_neighbor_bytes */
    ghost_neighbor_bytes = t8_partition_compute_gnb (cmesh_from,
                                                     &send_as_ghost);
    /* parse through send_as_ghost to compute ghost_attribute_bytes and attr_info_bytes */
    ghost_attribute_bytes = t8_partition_compute_gab (cmesh_from,
                                                      &send_as_ghost,
                                                      &ghost_attr_info_bytes);
    /* Total number of bytes that we send to the other process */
    total_alloc = num_trees * sizeof (t8_ctree_struct_t) +
      num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes +
      tree_neighbor_bytes + attr_info_bytes + attr_bytes +
      ghost_attribute_bytes + ghost_attr_info_bytes;
    /* Extra space to store the number of trees and ghosts in the send buffer */
    total_alloc += 2 * sizeof (t8_locidx_t);

    /* Debugging output about shipped trees and ghosts. */
    t8_debugf ("NT %i -- each %zdB\n", num_trees, sizeof (t8_ctree_struct_t));
    t8_debugf ("NGS %i - each %zdB\n", num_ghost_send,
               sizeof (t8_cghost_struct_t));
    t8_debugf ("GNB %zd\n", ghost_neighbor_bytes);
    t8_debugf ("TNB %zd\n", tree_neighbor_bytes);
    t8_debugf ("AIB %zd\n", attr_info_bytes);
    t8_debugf ("AB %zd\n", attr_bytes);
    t8_debugf ("GAIB %zd\n", ghost_attr_info_bytes);
    t8_debugf ("GAB %zd\n", ghost_attribute_bytes);
    t8_debugf ("Ta %zd\n", total_alloc);
    /* If profiling is enabled, we count the number of shipped trees/ghosts
     * and processes we ship to. */
    if (cmesh->profile) {
      if (iproc != cmesh->mpirank) {
        /* We do not count if we send to ourself. */
        cmesh->profile->partition_ghosts_shipped += num_ghost_send;
        cmesh->profile->partition_trees_shipped += num_trees;
        cmesh->profile->partition_procs_sent++;
        cmesh->profile->partition_bytes_sent += total_alloc;
      }
    }
    if (iproc != cmesh->mpirank) {
      buffer = (*send_buffer)[iproc - *send_first - flag] = T8_ALLOC (char,
                                                                      total_alloc);
    }
    else if (num_trees > 0 || num_ghost_send > 0) {
      *my_buffer = buffer = T8_ALLOC (char, total_alloc);
      *my_buffer_bytes = total_alloc;
    }
    else {
      my_buffer = NULL;
      buffer = NULL;
    }
    T8_ASSERT (num_trees + num_ghost_send == 0 || (!cmesh_from->set_partition
                                                   && iproc ==
                                                   cmesh_from->mpirank)
               || t8_offset_sendsto (cmesh->mpirank, iproc, offset_from,
                                     offset_to));

    /* Copy all data to the send buffer */
    t8_cmesh_partition_copy_data (buffer, cmesh,
                                  cmesh_from, num_trees, attr_info_bytes,
                                  ghost_attr_info_bytes,
                                  ghost_neighbor_bytes,
                                  tree_neighbor_bytes, attr_bytes,
                                  &send_as_ghost,
                                  range_start, range_end, total_alloc, iproc);

    /* If we send to a remote process we post the MPI_Isend here */
    if (iproc != cmesh->mpirank) {
      if (num_trees + num_ghost_send > 0) {
        /* send buffer to remote process */
        t8_debugf ("Post send of %i trees/%zd bytes to %i\n",
                   *(t8_locidx_t *) (buffer +
                                     total_alloc - 2 * sizeof (t8_locidx_t)),
                   total_alloc, iproc - flag);
        mpiret =
          sc_MPI_Isend (buffer, total_alloc,
                        sc_MPI_BYTE, iproc, T8_MPI_PARTITION_CMESH,
                        comm, *requests + iproc - flag - *send_first);
        SC_CHECK_MPI (mpiret);
        num_send_mpi++;
      }
      else {
        /* If num_trees + num_ghost_send = 0 we do not post an MPI_Send
         * and we set the associated request to NULL in order for it to
         * be ignored in the MPI_Waitall call later */
        t8_debugf ("Set request %i to NULL\n", iproc - flag - *send_first);
        *(*requests + iproc - flag - *send_first) = sc_MPI_REQUEST_NULL;
      }
    }
    /* Calculate the next process we send to */
    iproc++;
    while (iproc < *send_last && !t8_offset_sendsto (cmesh_from->mpirank,
                                                     iproc, offset_from,
                                                     offset_to)) {
      /* Skip processes we do not send to and set their requests to NULL */

      if (*send_first <= cmesh->mpirank && iproc > cmesh->mpirank) {
        flag = 1;
      }
      t8_debugf ("Set request %i to NULL\n", iproc - flag - *send_first);
      *(*requests + iproc - flag - *send_first) = sc_MPI_REQUEST_NULL;
      (*send_buffer)[iproc - flag - *send_first] = NULL;
      iproc++;
    }
    if (iproc <= *send_last) {
      /* compute new ranges here */
      /* If tree_offset of iproc is < 0 the last tree is shared and has to
       * be send to the next process as well, except when this process already
       * has the last tree */
      /* Set range_start to the first tree that iproc needs */
      t8_gloidx_t         first_send;

      first_send = t8_offset_first (iproc, offset_to);
      while (first_send < t8_cmesh_get_num_local_trees (cmesh_from)
             && !t8_offset_sendstree (cmesh_from->mpirank, iproc,
                                      first_send, offset_from, offset_to)) {

        first_send++;
      }
      range_start = first_send - cmesh_from->first_tree;
      t8_debugf ("RS: %i\n", range_start);

      /* add number of trees in iproc to send to range_end */
      /* We have to be careful with locidx overflow when we go out of bounds
       * of our process */

      /* Set range end to the last tree of the receiver or to our last tree,
       * whichever is smaller */
      range_end =
        SC_MIN (t8_offset_last (iproc, offset_to) - cmesh_from->first_tree,
                cmesh_from->num_local_trees - 1);
      t8_debugf ("RE: %i\n", range_end);
      /* If the last tree was already present on the receiving process
       * we exclude it from sending */
      if (cmesh_from->mpirank != iproc &&
          range_end == cmesh_from->num_local_trees - 1
          && !t8_offset_empty (iproc, offset_from)
          && t8_offset_first (iproc, offset_from)
          == cmesh_from->first_tree + cmesh_from->num_local_trees - 1) {
        range_end--;
        t8_debugf ("RE: %i\n", range_end);
      }
      if (range_end < range_start
          || range_start >= t8_cmesh_get_num_local_trees (cmesh_from)) {
        /* We do not send to this process and are finished */
        iproc = *send_last + 1;
        *send_last = -2;
      }
    }
  }                             /* sending loop ends here */
  T8_FREE (ghost_flag_send);
  sc_array_reset (&send_as_ghost);
  t8_debugf ("End send loop\n");
  return num_send_mpi;
}

void
t8_cmesh_partition_receive_message (t8_cmesh_t cmesh, sc_MPI_Comm comm,
                                    int proc_recv, sc_MPI_Status * status,
                                    int *local_procid, int recv_first,
                                    t8_locidx_t * num_ghosts)
{
  int                 mpiret;
  int                 recv_bytes;
  t8_part_tree_t      recv_part;

  T8_ASSERT (proc_recv == status->MPI_SOURCE);
  T8_ASSERT (status->MPI_TAG == T8_MPI_PARTITION_CMESH);

  mpiret = sc_MPI_Get_count (status, sc_MPI_BYTE, &recv_bytes);
  SC_CHECK_MPI (mpiret);
  /* Allocate receive buffer */
  recv_part =
    t8_cmesh_trees_get_part (cmesh->trees,
                             local_procid[proc_recv - recv_first]);
  /* take first tree of part and allocate recv_bytes */
  recv_part->first_tree = T8_ALLOC (char, recv_bytes);
  /* Receive message */
  mpiret = sc_MPI_Recv (recv_part->first_tree, recv_bytes, sc_MPI_BYTE,
                        proc_recv, T8_MPI_PARTITION_CMESH, comm,
                        sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);
  /* Read num trees and num ghosts */
  recv_part->num_trees =
    *((t8_locidx_t *) (recv_part->first_tree + recv_bytes -
                       2 * sizeof (t8_locidx_t)));
  recv_part->num_ghosts =
    *((t8_locidx_t *) (recv_part->first_tree + recv_bytes -
                       sizeof (t8_locidx_t)));
  *num_ghosts += recv_part->num_ghosts;

  t8_debugf ("Received %i trees/%i ghosts/%i bytes from %i to %i\n",
             recv_part->num_trees, recv_part->num_ghosts, recv_bytes,
             proc_recv, local_procid[proc_recv - recv_first]);
  /* If we are profiling, we count the number of trees and ghosts that
   * we received. */
  if (cmesh->profile != NULL && proc_recv != cmesh->mpirank) {
    cmesh->profile->partition_ghosts_recv += recv_part->num_ghosts;
    cmesh->profile->partition_trees_recv += recv_part->num_trees;
  }
#if 0
  /* Free memory that was only used to store num_trees/ghosts */
  /* TODO: If we want to do this properly we have to realloc the memory,
   * if realloc is not enabled this means copying all the bytes.
   * Right now we do not free this memory here, but rely on it to be freed when
   * the trees_part is destroyed */
  T8_FREE (recv_part->first_tree + recv_bytes - 2 * sizeof (t8_locidx_t));
#endif
}

/* Loop over all the processes that we receive trees from and get the
 * MPI messages. */
/* stores the number of received ghosts in num_ghosts */
/* fr and lr are only for debugging, see t8_cmesh_partition_debug_listprocs */
/* TODO: Remove the const qualifier at the cmesh_from parameter */
static void
t8_cmesh_partition_recvloop (t8_cmesh_t cmesh,
                             const struct t8_cmesh *cmesh_from,
                             t8_gloidx_t * tree_offset,
                             char *my_buffer, size_t my_buffer_bytes,
                             sc_MPI_Comm comm, int fr, int lr)
{
  int                 num_receive, *local_procid;       /* ranks of the processor from which we will receive */
  int                 mpiret, proc_recv, iproc;
  int                 recv_first, recv_last;    /* ranks of the processor from which we will receive */
  int                 num_parts, myrank_part;
  t8_locidx_t         num_trees, num_ghosts;
  t8_gloidx_t        *from_offsets;
  t8_part_tree_t      recv_part;
  sc_MPI_Status       status;
#if 0
  sc_list_t          *possible_receivers;       /* Linked list storing the ranks from which
                                                   we still expect a message */
  /* TODO: These variables belong to the polling code */
  int                *poss_recv_data;   /* Data storage for these ranks */
  sc_link_t          *iterate;  /* Iterator through the linked list */
  sc_link_t          *prev;     /* We need to store iterates predecessor in order to
                                   be able to remove iterate when needed */
  int                 iprobe_flag;
#endif

  num_trees = t8_offset_num_trees (cmesh->mpirank, tree_offset);
  /* Receive from other processes */
  if (cmesh_from->set_partition) {
    T8_ASSERT (cmesh_from->tree_offsets != NULL);
    from_offsets = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
    t8_cmesh_partition_alternative_recvrange (cmesh, (t8_cmesh_t) cmesh_from,
                                              &recv_first, &recv_last);
#if T8_ENABLE_DEBUG
    if (recv_first <= recv_last) {
      T8_ASSERT (fr == recv_first);
      T8_ASSERT (lr == recv_last);
    }
#endif

    num_parts = t8_offset_range_send (recv_first, recv_last, cmesh->mpirank,
                                      from_offsets, tree_offset);
  }
  else {
    recv_first = cmesh->mpirank;
    recv_last = cmesh->mpirank;
    num_parts = my_buffer != NULL ? 1 : 0;
    from_offsets = NULL;
  }
  /* Initialize trees structure with yet unknown number of ghosts */
  t8_cmesh_trees_init (&cmesh->trees, num_parts, num_trees, 0);
  num_ghosts = 0;

  /* Total number of processors from which we receive a MPI message is one
   * less if we counted our own rank. */
  if (cmesh_from->set_partition) {
    if (t8_offset_sendsto (cmesh->mpirank, cmesh->mpirank, from_offsets,
                           tree_offset)) {
      num_receive = num_parts - 1;
    }
    else {
      num_receive = num_parts;
    }
  }
  else {
    num_receive = 0;
  }

  num_receive = SC_MAX (num_receive, 0);        /* num_receive could get negative, so we set it
                                                   to 0 in that case, since we need to pass it to
                                                   malloc */
  if (recv_first < 0 || recv_last < recv_first) {
    recv_last = -2;
    recv_first = -1;
  }
  /* stores at i-th position the new local proc id of the respective process.
   * This is important to account for empty processes */
  local_procid = T8_ALLOC_ZERO (int, recv_last - recv_first + 1);
  /* Allocate memory */
  for (iproc = 1; iproc < recv_last - recv_first + 1; iproc++) {
    local_procid[iproc] =
      !t8_offset_sendsto (iproc + recv_first, cmesh->mpirank, from_offsets,
                          tree_offset) ?
      local_procid[iproc - 1] : local_procid[iproc - 1] + 1;
  }
  if (cmesh_from->set_partition) {
    if (recv_first <= cmesh->mpirank && cmesh->mpirank <= recv_last
        && t8_offset_sendsto (cmesh->mpirank, cmesh->mpirank, from_offsets,
                              tree_offset)) {
      myrank_part = local_procid[cmesh->mpirank - recv_first];
    }
    else {
      myrank_part = -1;
    }
  }
  else {
    /* cmesh_from is replicated */
    myrank_part = my_buffer != NULL ? 0 : -1;
  }
  T8_ASSERT (my_buffer == NULL || myrank_part >= 0);

  T8_ASSERT (recv_first == -1 || !cmesh_from->set_partition ||
             !t8_offset_empty (recv_first, from_offsets));

  /**************************************************/
  /*            Receive MPI messages                */
  /**************************************************/

  /****     Setup     ****/

  if (num_receive > 0) {
#if 0
    /* TODO: This belongs to the polling MPI communication, see below */
    /* Find first process from which we will receive */
    proc_recv = recv_first;
    /* Check whether we expect an MPI message from this process */
    while (proc_recv == cmesh->mpirank ||
           !t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                               tree_offset)) {
      proc_recv++;
    }

    poss_recv_data = T8_ALLOC (int, num_receive);
    possible_receivers = sc_list_new (NULL);

    /* Loop over all processes that we receive from and create an entry in a linked
     * list for each one */
    for (iproc = 0; iproc < num_receive && proc_recv <= recv_last; iproc++) {
      /* For each process we expect a message from we create an entry in a linked list */
      poss_recv_data[iproc] = proc_recv;
      sc_list_append (possible_receivers, (void *) (poss_recv_data + iproc));
      t8_debugf ("Added %i to list\n", proc_recv);
      /* Calculate next rank from which we receive a message */
      proc_recv++;
      /* Check whether we expect an MPI message from this process, otherwise increment */
      while (proc_recv <= recv_last &&
             (proc_recv == cmesh->mpirank
              || !t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                                     tree_offset))) {
        proc_recv++;
      }
    }
#endif

    /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for an message for each
     * sender and if there is one we receive it and remove the sender from
     * the list.
     * The last message can be received via probe */
    while (num_receive > 0) {
      t8_debugf ("Probing for %i messages.\n", num_receive);
      mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_PARTITION_CMESH, comm,
                             &status);
      SC_CHECK_MPI (mpiret);
      num_receive--;
#if 0
      prev = NULL;
      iprobe_flag = 0;
      /* TODO: This part of the code uses polling to receive the
       *       messages. Remove if the unpolling version has prooved to run
       *       well.
       */
      for (iterate = possible_receivers->first;
           iterate != NULL && iprobe_flag == 0;) {
        /* Iterate through all entries and probe for message */
        proc_recv = *(int *) iterate->data;
        mpiret = sc_MPI_Iprobe (proc_recv, T8_MPI_PARTITION_CMESH, comm,
                                &iprobe_flag, &status);
        SC_CHECK_MPI (mpiret);
        if (iprobe_flag == 0) {
          /* We do only continue if there is no message to receive */
          prev = iterate;
          iterate = iterate->next;
        }
      }
#endif
#if 0
      /* polling code, see above */
      if (iprobe_flag != 0) {
#endif
        /* There is a message to receive */
        proc_recv = status.MPI_SOURCE;
        /* TODO: assert that proc_recv is still contained in the list of receivers. */
        T8_ASSERT (status.MPI_TAG == T8_MPI_PARTITION_CMESH);
        T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last &&
                   t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                                      tree_offset));
        t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status,
                                            local_procid, recv_first,
                                            &num_ghosts);
#if 0
        sc_list_remove (possible_receivers, prev);
        /* polling code, see above */
      }
#endif
    }

#if 0
    /* TODO: polling code */
    /* We have one message left and does we do not need to ause ANY_SOURCE
     * but call the blocking Probe with the remaining process instead. */
    T8_ASSERT (possible_receivers->elem_count == 1);
    proc_recv = *(int *) sc_list_pop (possible_receivers);
    sc_list_destroy (possible_receivers);
    T8_FREE (poss_recv_data);
    mpiret = sc_MPI_Probe (proc_recv, T8_MPI_PARTITION_CMESH, comm, &status);
    SC_CHECK_MPI (mpiret);
    t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status,
                                        local_procid, recv_first,
                                        &num_ghosts);
#endif
  }

#if 0
  for (iproc = 0; iproc < num_receive; iproc++) {
    t8_debugf ("Start receive\n");

    /* Blocking test for message. */
    /* TODO: Right now we are waiting for each message in order to be received.
     *       We do this, since when we use Iprobe and ANY_SOURCE, we could receive
     *       messages that other processors send in later calls to partition, which
     *       then mess up the current partition.
     *       It would certainly be more efficient not to receive the messages in order
     *       of the receiving processes but in temporal order of receivement.
     */
    /* TODO: Proposed solution is to do an MPI_IProbe for each process that we
     *       receive from and then do MPI_Waitany/some until all messages are
     *       received.
     *       Problem if one process sends in the next round to the same one?
     */
    t8_debugf ("Probing for message from %i\n", proc_recv);
    mpiret = sc_MPI_Probe (proc_recv, T8_MPI_PARTITION_CMESH, comm, &status);
    SC_CHECK_MPI (mpiret);

    T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last &&
               t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                                  tree_offset));

    t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status,
                                        local_procid, recv_first,
                                        &num_ghosts);

    /* Calculate next rank from which we receive a message */
    proc_recv++;
    /* Check whether we expect an MPI message from this process */
    while (proc_recv <= recv_last &&
           (proc_recv == cmesh->mpirank
            || !t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                                   tree_offset))) {
      proc_recv++;
    }
  }
#endif
  t8_debugf ("End receive\n");

  /**************************************************/
  /*            End receiving MPI messages          */
  /**************************************************/

  /* Got trees and ghosts from myself */
  if (my_buffer != NULL) {
    T8_ASSERT (myrank_part >= 0);
    /* TODO: Get ghosts from myself! */
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, myrank_part);
    recv_part->first_tree = my_buffer;
    /* Read num trees and num ghosts */
    recv_part->num_trees = *(t8_locidx_t *) (recv_part->first_tree +
                                             my_buffer_bytes
                                             - 2 * sizeof (t8_locidx_t));
    recv_part->num_ghosts = *(t8_locidx_t *) (recv_part->first_tree +
                                              my_buffer_bytes
                                              - sizeof (t8_locidx_t));
    num_ghosts += recv_part->num_ghosts;

    t8_debugf ("Received %i trees/%i ghosts from myself to %i\n",
               recv_part->num_trees, recv_part->num_ghosts, myrank_part);
  }
  else if (myrank_part >= 0) {
    /* We did not receive trees from ourself, however we allocated the part,
     * so we set everything to 0 there */
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, myrank_part);
    recv_part->num_ghosts = recv_part->num_trees = 0;
    recv_part->first_tree = NULL;
    /* TODO: This case should not happen, eventually remove */
    SC_ABORT_NOT_REACHED ();
  }
  t8_debugf ("Total number of ghosts in new partition: %i\n", num_ghosts);
  cmesh->trees->ghost_to_proc = T8_ALLOC (int, num_ghosts);
  cmesh->num_ghosts = num_ghosts;
  T8_FREE (local_procid);
}

static void
t8_cmesh_partition_debug_listprocs (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                                    sc_MPI_Comm comm, int *fs, int *ls,
                                    int *fr, int *lr)
{
  int                 mpiret, mpisize, mpirank, p;
  char                out[BUFSIZ] = { };
  t8_gloidx_t        *from, *to;

  if (cmesh_from->set_partition) {
    from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  }
  else {
    from = NULL;
  }
  to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  *fs = *fr = mpisize;
  *ls = *lr = 0;
  for (p = 0; p < mpisize; p++) {
    if (t8_offset_sendsto (mpirank, p, from, to)) {
      snprintf (out + strlen (out), BUFSIZ - strlen (out), "%i%c ", p,
                p == mpisize - 1 ? '!' : ',');
      *fs = SC_MIN (*fs, p);
      *ls = SC_MAX (*ls, p);
    }
  }
  t8_debugf ("I send to: %s\n", out);
  sprintf (out, " ");
  if (cmesh_from->set_partition) {
    for (p = 0; p < mpisize; p++) {
      if (t8_offset_sendsto (p, mpirank, from, to)) {
        snprintf (out + strlen (out), BUFSIZ - strlen (out), "%i%c ", p,
                  p == mpisize - 1 ? '!' : ',');
        *fr = SC_MIN (*fr, p);
        *lr = SC_MAX (*lr, p);
      }
    }
  }
  else {
    *fr = *lr = cmesh_from->mpirank;
    snprintf (out, BUFSIZ, "%i", cmesh_from->mpirank);
  }
  t8_debugf ("I receive from: %s\n", out);
}

/* Given an initial cmesh (cmesh_from) and a new partition table (tree_offset)
 * create the new partition on the destination cmesh (cmesh) */
/* TODO: remove offset argument and use cmesh_from.tree_offsets */
/* TODO: remove const */
static void
t8_cmesh_partition_given (t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from,
                          t8_gloidx_t * tree_offset, sc_MPI_Comm comm)
{
  int                 send_first, send_last, num_request_alloc; /* ranks of the processor to which we will send */
  int                 iproc, num_send_mpi, mpiret;
  size_t              my_buffer_bytes;
  char              **send_buffer = NULL, *my_buffer = NULL;

  int                 fs, ls, fr, lr;

  sc_MPI_Request     *requests = NULL;
  t8_locidx_t         num_ghosts, itree, num_trees;
  t8_part_tree_t      recv_part;
  t8_ctree_t          tree;

  /* TODO: computing recv information needs the shared array of the old partition in cmesh_from,
   *       thus, two of those huge arrays need to exist at the same time.
   *       Using MPI-2.1 One-sided communication we could resolve this, since
   *       pi \in pj_recv if and only if pj \in pi_send. The latter can be computed
   *       using only the new partition array w/o communication and then via one-sided
   *       communication we can compute the number of receiving processes on each process, which
   *       should be enough to receive all messages in a while loop. */

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partition);
  T8_ASSERT (cmesh_from != NULL);
  T8_ASSERT (cmesh_from->committed);

  /* determine send and receive range. temp_tree is last local tree of send_first in new partition */
  cmesh->first_tree = t8_offset_first (cmesh->mpirank, tree_offset);
  cmesh->num_local_trees = t8_offset_num_trees (cmesh->mpirank, tree_offset);

  if (cmesh_from->set_partition) {
    t8_cmesh_partition_debug_listprocs (cmesh, (t8_cmesh_t) cmesh_from, comm,
                                        &fs, &ls, &fr, &lr);
  }

  /*********************************************/
  /*        Done with setup                    */
  /*********************************************/

  /* Send all trees and ghosts */
  num_send_mpi =
    t8_cmesh_partition_sendloop (cmesh, (t8_cmesh_t) cmesh_from,
                                 &num_request_alloc, &send_first, &send_last,
                                 &send_buffer, &my_buffer,
                                 &my_buffer_bytes, &requests, comm);
  T8_ASSERT (!cmesh_from->set_partition || send_first == -1
             || send_first == fs);
  T8_ASSERT (!cmesh_from->set_partition || send_last == -2
             || send_last == ls);

  /* receive all trees and ghosts */
  t8_cmesh_partition_recvloop (cmesh, cmesh_from, tree_offset, my_buffer,
                               my_buffer_bytes, comm, fr, lr);
  if (num_send_mpi > 0) {
    mpiret = sc_MPI_Waitall (num_request_alloc, requests,
                             sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean-up */
  for (iproc = 0; iproc < send_last - send_first +
       !(cmesh->mpirank >= send_first && cmesh->mpirank <= send_last);
       iproc++) {
    T8_FREE (send_buffer[iproc]);
  }
  T8_FREE (send_buffer);
  T8_FREE (requests);
  /* Done with Clean-up */

  /* set recv_part->first_tree_id/first_ghost_id */
  /* Use as temporary variables to store the first tree_id/ghost_id of the new parts */
  num_trees = 0;
  num_ghosts = 0;
  for (iproc = 0; iproc < (int) cmesh->trees->from_proc->elem_count; iproc++) {
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, iproc);
    /* Set first ids */
    recv_part->first_tree_id = num_trees;
    recv_part->first_ghost_id = num_ghosts;
    /* count trees_per_eclass */
    for (itree = recv_part->first_tree_id;
         itree < recv_part->first_tree_id + recv_part->num_trees; itree++) {
      cmesh->trees->tree_to_proc[itree] = iproc;
      tree = t8_cmesh_trees_get_tree (cmesh->trees, itree);
      tree->treeid = itree;
      cmesh->num_local_trees_per_eclass[tree->eclass]++;
    }
    /* Calculate num_trees and num_ghosts for next part */
    num_trees += recv_part->num_trees;
    num_ghosts += recv_part->num_ghosts;
  }
  num_ghosts = 0;
  for (iproc = 0; iproc < (int) cmesh->trees->from_proc->elem_count; iproc++) {
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, iproc);
    /* Assing new local ids to the ghosts of this part, also set ghost_to_proc */
    t8_partition_new_ghost_ids (cmesh, recv_part, num_ghosts, iproc);
    /* We need to do this in a second loop, since all the tree_to_proc entries have
     * to be set before. This is because we may be accessing any tree in the new cmesh.  */
    num_ghosts += recv_part->num_ghosts;
  }
  T8_ASSERT (cmesh->num_local_trees == num_trees);
  cmesh->num_ghosts = num_ghosts;
}

/* Given a cmesh which is to be partitioned, execute the partition task.
 * This includes partitioning by uniform level and partitioning from a second cmesh */
/* TODO: Check whether the input data is consistent.
 *       If tree_offset is set on one process it has to be set on each process.
 *       If first_tree is set   "       "         "
 *       what else? */
void
t8_cmesh_partition (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh_from;
  t8_gloidx_t         last_tree, *tree_offsets;

  T8_ASSERT (t8_cmesh_is_committed (cmesh->set_from));
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partition);

  t8_global_productionf ("Enter cmesh partition\n");
  /* If profiling is enabled, we measure the runtime of this routine. */
  if (cmesh->profile != NULL) {
    cmesh->profile->partition_runtime = sc_MPI_Wtime ();
  }
  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  cmesh->num_trees = cmesh_from->num_trees;

  /**********************************************/
  /*      Compute local number of trees         */
  /*         and trees per proc array           */
  /**********************************************/
  if (cmesh->set_partition_level >= 0) {
    /* Compute first and last tree index */
    T8_ASSERT (cmesh->tree_offsets == NULL);
    t8_cmesh_uniform_bounds (cmesh_from, cmesh->set_partition_level,
                             &cmesh->first_tree, NULL, &last_tree, NULL,
                             &cmesh->first_tree_shared);
    cmesh->num_local_trees = last_tree - cmesh->first_tree + 1;
    /* Compute the tree offset */
    t8_cmesh_gather_treecount_nocommit (cmesh, comm);
    /* Set the tree offsets to the cmesh's offset array */
    tree_offsets = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  }
  else {
    /* We compute the partition after a given partition table in cmesh->tree_offsets */
    T8_ASSERT (cmesh->tree_offsets != NULL);
    tree_offsets = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
    /* Check whether the offset at mpirank is smaller 0, if so then the
     * first local tree is also the last local tree of the next smaller nonempty process. */
    cmesh->first_tree_shared = tree_offsets[cmesh->mpirank] < 0;
    /* compute local first tree */
    cmesh->first_tree = t8_offset_first (cmesh->mpirank, tree_offsets);
    /* compute local num trees */
    cmesh->num_local_trees =
      t8_offset_num_trees (cmesh->mpirank, tree_offsets);
  }
  if (cmesh->set_from->set_partition && cmesh->set_from->tree_offsets == NULL) {
    /* Create the partition table for cmesh_from */
    t8_cmesh_gather_treecount (cmesh->set_from, comm);
  }
  t8_debugf ("Partition from:\n");
  t8_cmesh_offset_print (cmesh->set_from, comm);
  t8_debugf ("To:\n");
  t8_cmesh_offset_print (cmesh, comm);
  /***************************************************/
  /*        Done with local num and tree_offset      */
  /***************************************************/
  t8_cmesh_partition_given (cmesh, cmesh->set_from, tree_offsets, comm);
  /* If profiling is enabled, we measure the runtime of this routine. */
  if (cmesh->profile) {
    /* Runtime = current_time - start_time */
    cmesh->profile->partition_runtime = sc_MPI_Wtime ()
      - cmesh->profile->partition_runtime;
  }
  t8_global_productionf ("Done cmesh partition\n");
}

void
t8_cmesh_offset_print (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
#if T8_ENABLE_DEBUG
  int                 offset_isnew = 0;

  if (cmesh->set_partition) {
    if (cmesh->tree_offsets == NULL) {
      t8_cmesh_gather_treecount (cmesh, comm);
      offset_isnew = 1;
    }
    t8_offset_print (cmesh->tree_offsets, comm);
    T8_ASSERT (t8_offset_consistent (cmesh->mpisize,
                                     t8_shmem_array_get_gloidx_array
                                     (cmesh->tree_offsets),
                                     cmesh->num_trees));
    if (offset_isnew == 1) {
      t8_shmem_array_destroy (&cmesh->tree_offsets);
      T8_ASSERT (cmesh->tree_offsets == NULL);
    }
    else {
      t8_debugf ("Replicated cmesh with %lli trees.\n",
                 (long long) cmesh->num_trees);
    }
  }
#endif
}

/* Create a partition that concentrates everything at a given proc */
t8_shmem_array_t
t8_cmesh_offset_concentrate (int proc, sc_MPI_Comm comm,
                             t8_gloidx_t num_trees)
{
  int                 mpirank, mpiret, mpisize, iproc;
  t8_shmem_array_t    shmem_array;
  t8_gloidx_t        *offsets;
#ifdef T8_ENABLE_DEBUG
  char                out[BUFSIZ] = "";
#endif

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  shmem_array = t8_cmesh_alloc_offsets (mpisize, comm);
  offsets = t8_shmem_array_get_gloidx_array (shmem_array);
  offsets[0] = 0;
  for (iproc = 1; iproc <= mpisize; iproc++) {
    if (iproc == proc + 1) {
      offsets[iproc] = num_trees;
    }
    else {
      offsets[iproc] = offsets[iproc - 1];
    }
#ifdef T8_ENABLE_DEBUG
    snprintf (out + strlen (out), BUFSIZ - strlen (out), "%li,",
              offsets[iproc]);
#endif
  }
#ifdef T8_ENABLE_DEBUG
  t8_debugf ("Partition with offsets:0,%s\n", out);
#endif

  T8_ASSERT (t8_offset_consistent (mpisize, offsets, num_trees));
  return shmem_array;
}

/* Create a random partition */
/* if shared is nonzero than first trees can be shared */
t8_shmem_array_t
t8_cmesh_offset_random (sc_MPI_Comm comm, t8_gloidx_t num_trees, int shared,
                        unsigned seed)
{
  int                 iproc, mpisize, mpiret, random_number, mpirank;
  int                 first_shared;
  int                 i;
  unsigned            u_seed;
  t8_gloidx_t        *offsets;
  t8_shmem_array_t    shmem_array;
  t8_gloidx_t         new_first;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  shmem_array = t8_cmesh_alloc_offsets (mpisize, comm);

  offsets = t8_shmem_array_get_gloidx_array (shmem_array);

  if (seed == 0) {
    u_seed = sc_MPI_Wtime () * 10000;
  }
  else {
    u_seed = seed;
  }

  if (mpirank == 0) {
    t8_debugf ("Random number seed = %u\n", u_seed);
  }
  mpiret = sc_MPI_Bcast (&u_seed, 1, sc_MPI_INT, 0, comm);
  SC_CHECK_MPI (mpiret);
  srand (u_seed);

  offsets[0] = 0;
  first_shared = 0;
  for (iproc = 1; iproc < mpisize; iproc++) {
    offsets[iproc] = 0;
    /* Create a random number between 0 and 200% of an ideal partition */
    /* This is the number of trees on process iproc-1. */
    if ((int) (num_trees * 2. / mpisize) == 0) {
      /* This case prevents division by 0 */
      random_number = 1;
    }
    else {
      random_number = rand () % (int) (num_trees * 2. / mpisize);
    }

    if (random_number == 0 && first_shared) {
      /* The previous proc is empty but set its first tree to be shared. */
      /* We have to manually reset the shared flag. */
      offsets[iproc - 1] = -offsets[iproc - 1] - 1;
      first_shared = 0;
    }
    random_number += first_shared;
    /* If we would exceed the number of trees we cut the random number */
    new_first = t8_offset_first (iproc - 1, offsets) + random_number;
    if (new_first > num_trees) {
      random_number = num_trees - t8_offset_first (iproc - 1, offsets);
      new_first = num_trees;
    }
    if (shared && new_first < num_trees) {      /* new first is num_trees, this process must be empty */
      first_shared = rand () % 2;
    }
    else {
      first_shared = 0;
    }

    offsets[iproc] = random_number + t8_offset_first (iproc - 1, offsets);
    if (first_shared && offsets[iproc] != num_trees) {
      offsets[iproc] = -offsets[iproc] - 1;
    }
  }
  offsets[mpisize] = num_trees;

  for (i = 0; i < mpisize; i++) {
    sc_MPI_Barrier (comm);
    if (i == mpirank) {
      for (iproc = 0; iproc < mpisize + 1; iproc++) {
        t8_debugf ("[H] %li\n", offsets[iproc]);
      }
    }
    seed = -1;
  }
  T8_ASSERT (t8_offset_consistent (mpisize, offsets, num_trees));
  return shmem_array;
}

/* TODO: Check that percent is the same on each process */
t8_shmem_array_t
t8_cmesh_offset_percent (t8_cmesh_t cmesh, sc_MPI_Comm comm, int percent)
{
  t8_gloidx_t         new_first_tree, old_first_tree;
  t8_locidx_t         old_num_trees_pm1;
  t8_shmem_array_t    partition_array;
  t8_gloidx_t        *old_partition;
  int                 mpirank, mpisize, mpiret;
  int                 created = 0;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Get old partition table and check for existence. */
  /* TODO: If it does not exists (NULL is returned), we could compute
   * the number of trees and first tree of the next smaller rank by hand. */
  if (cmesh->tree_offsets == NULL) {
    /* We need to create the old partition array. */
    t8_cmesh_gather_treecount (cmesh, comm);
    created = 1;
  }
  old_partition = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  T8_ASSERT (old_partition != NULL);

  /* Allocate new partition array. */
  partition_array = t8_cmesh_alloc_offsets (mpisize, comm);
  /* get the first local tree from the current cmesh */
  old_first_tree = t8_cmesh_get_first_treeid (cmesh);
  /* Get the number of local trees of next smaller process, or
   * 0 if we are rank 0. */
  old_num_trees_pm1 = t8_offset_num_trees (mpirank > 0 ? mpirank - 1 : 0,
                                           old_partition);
  /* Compute the new first local tree */
  new_first_tree = old_first_tree - old_num_trees_pm1 * percent / 100;
  /* Compute the new entry in the offset array.
   * If the old first tree was shared, then the new one will be as well
   * and if not it will not be shared. */
  if (mpirank != 0) {
    new_first_tree = t8_offset_first_tree_to_entry (new_first_tree,
                                                    cmesh->first_tree_shared);
  }
  else {
    new_first_tree = 0;
  }
  t8_shmem_array_allgather (&new_first_tree, 1, T8_MPI_GLOIDX,
                            partition_array, 1, T8_MPI_GLOIDX);
  t8_shmem_array_set_gloidx (partition_array, mpisize,
                             t8_cmesh_get_num_trees (cmesh));
  if (created) {
    /* We needed to create the old partition array and thus we clean it up
     * again. */
    t8_shmem_array_destroy (&cmesh->tree_offsets);
  }
  T8_ASSERT (t8_offset_consistent (mpisize,
                                   t8_shmem_array_get_gloidx_array
                                   (partition_array),
                                   t8_cmesh_get_num_trees (cmesh)));
  return partition_array;
}

/* Create a repartition array, where each process sends half of its
 * trees to the next process. The last process does not send any trees. */
/* TODO: This function was not tested with shared trees yet. */
t8_shmem_array_t
t8_cmesh_offset_half (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  return t8_cmesh_offset_percent (cmesh, comm, 50);
}
