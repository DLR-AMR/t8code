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

/* Return the minimum of two t8_gloidx_t's */
static              t8_gloidx_t
t8_glo_min (t8_gloidx_t A, t8_gloidx_t B)
{
  return A < B ? A : B;
}

/* Return -1 if A is smaller 0
 * Return 0 if not. */
static int
t8_glo_kl0 (t8_gloidx_t A)
{
  return A < 0 ? -1 : 0;
}

/* The first tree of a given process in a partition */
static              t8_gloidx_t
t8_offset_first (int proc, t8_gloidx_t * offset)
{
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);
  return T8_GLOIDX_ABS (offset[proc]) + t8_glo_kl0 (offset[proc]);
}

/* The number of trees of a given process in a partition */
/* Can get negative for empty processes */
static              t8_gloidx_t
t8_offset_num_trees (int proc, t8_gloidx_t * offset)
{
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);

  return T8_GLOIDX_ABS (offset[proc + 1]) - t8_offset_first (proc, offset);
}

/* The last local tree of a given process in a partition */
static              t8_gloidx_t
t8_offset_last (int proc, t8_gloidx_t * offset)
{
  T8_ASSERT (proc >= -1);
  T8_ASSERT (offset != NULL);

  return T8_GLOIDX_ABS (offset[proc + 1]) - 1;
}

/* Return 1 if the process has no trees in the partition.
 * Return 0 if the process has at least one tree */
static int
t8_offset_empty (int proc, t8_gloidx_t * offset)
{
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);

  /* proc is empty if the first tree of the next process is smaller than
   * the "first tree" of proc.
   * In this case the first tree of proc+1 is shared. */
  if (t8_offset_first (proc + 1, offset) < t8_offset_first (proc, offset))
    return 1;
  /* Or the "first tree" of proc equals the first tree of proc+1 but the latter
   * is not shared */
  if (t8_offset_first (proc + 1, offset) == t8_offset_first (proc, offset)
      && offset[proc + 1] >= 0) {
    return 1;
  }
  return 0;
}

/* Determine whether a given global tree id is in the range of a given process */
static int
t8_offset_in_range (t8_gloidx_t tree_id, int proc, t8_gloidx_t * offset)
{
  return t8_offset_first (proc, offset) <= tree_id
    && tree_id <= t8_offset_last (proc, offset);
}

/* Compute a list of all processes that own a specific tree */
static void
t8_offset_all_owners_of_tree (int mpisize, t8_gloidx_t gtree, t8_gloidx_t * offset,
                              sc_array_t * owners)
{
  int proc, range[2], found, proc_temp;
  int *entry;
  T8_ASSERT (owners != NULL);
  T8_ASSERT (owners->elem_count == 0);
  T8_ASSERT (owners->elem_size == sizeof (int));
  T8_ASSERT (0 <= gtree && gtree <= t8_offset_last (mpisize - 1, offset));

  range[0] = 0;
  range[1] = mpisize - 1;
  /* At first we find any process that owns the tree with a binary search */
  found = 0;
  while (!found) {
    proc = (range[0] + range[1]) /2;
    if (t8_offset_in_range (gtree, proc, offset)) {
      found = 1;
    }
    else if (t8_offset_last (proc, offset) < gtree){
      /* look further right */
      range[0] = proc + 1;
    }
    else {
      range[1] = proc - 1;
    }
  }
  /* Now we find the smallest process that has the tree */
  proc_temp = proc;
  while (proc_temp >= 0 && t8_offset_in_range (gtree, proc_temp, offset)) {
    /* decrement temp as long as it holds the tree */
    proc_temp--;
    while (0 <= proc_temp && t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp--;
    }
  }
  /* We are now one nonempty proccess below the smallest process having the tree */
  if (proc_temp >= 0) {
    proc_temp++;
    while (t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp++;
    }
    T8_ASSERT (t8_offset_in_range (gtree, proc_temp, offset));
  }
  else {
    proc_temp = proc;
  }
  /* Add the first process to the array */
  proc = proc_temp;
  entry = (int*) sc_array_push (owners);
  *entry = proc;
  found = 0;
  /* Now we parse through all processes bigger than the firt one until
   * they do not have the tree anymore. */
  while (!found) {
    proc++;
    while (proc < mpisize && t8_offset_empty (proc, offset)) {
      /* Skip empty processes */
      proc++;
    }
    if (proc < mpisize && t8_offset_in_range (gtree, proc, offset)) {
      entry = (int*) sc_array_push (owners);
      *entry = proc;
    }
    else {
      found = 1;
    }
  }
}

/* Return 1 if the process will not send any trees, that is if it is
 * empty or has only one shared tree. */
static int
t8_offset_nosend (int proc,int mpisize, t8_gloidx_t * offset_from, t8_gloidx_t * offset_to)
{
  t8_gloidx_t num_trees;

  num_trees = t8_offset_num_trees (proc, offset_from);
  if (t8_offset_empty (proc, offset_from)) {
    return 1;
  }
  else if (num_trees <= 2) {
    int first_not_send, last_not_send = 0;
    /* We only have one or two trees */
    /* The first tree is not send if it is shared and will not be on
     * this process in the new partition */
    first_not_send = offset_from[proc] < 0
        && !t8_offset_in_range (t8_offset_first (proc, offset_from),proc, offset_to);

    if ((first_not_send && num_trees == 2)
        || (!first_not_send && num_trees == 1)) {
      int temp_proc, i;
      t8_gloidx_t last_tree = t8_offset_last (proc, offset_from);
      sc_array_t receivers;
      /* It could be that our last tree is not send either, this is the case
       * if all processes in the new partition that have it, already have it shared
       * in the current partition. */
      /* At first we find out if last tree is shared by anybody, since if not
       * we do not need to find all receivers. */
      if (t8_offset_in_range (last_tree, proc, offset_to)) {
        /* We keep our last tree, thus it is send by us */
          return 0;
      }

      temp_proc = proc + 1;
      /* Find next nonempty process */
      while (temp_proc < mpisize && t8_offset_empty (temp_proc, offset_from)) {
        temp_proc++;
      }
      if (temp_proc >= mpisize ||
          t8_offset_first (temp_proc, offset_from) != last_tree) {
        /* Our last tree is unique and thus we need to send it */
        return 0;
      }
      sc_array_init (&receivers, sizeof (int));
      /* Now we need to find out all process in the new partition that have last_tree
       * and check if any of them did not have it before. */
      t8_offset_all_owners_of_tree (mpisize, last_tree, offset_to, &receivers);
      for (i = 0; i < receivers.elem_count;i++) {
        temp_proc = *(int *) sc_array_index_int (&receivers, i);
        if (!t8_offset_in_range (last_tree, temp_proc, offset_from)) {
          /* We found at least one process that needs our last tree */
          sc_array_reset (&receivers);
          return 0;
        }
      }
      sc_array_reset (&receivers);
      last_not_send = 1;
    }
    if (num_trees - first_not_send - last_not_send <= 0) {
      /* We only have one tree, it is shared and we do not keep it.
       * Or we have only two trees and both are also on other procs and only
       * persist on these procs */
      return 1;
    }
  }
  return 0;
}


/* Return one if proca sends trees to procb when partitioning from
 * offset_from to offset_to */
static int
t8_offset_sendsto (int proca, int procb, t8_gloidx_t * t8_offset_from,
                   t8_gloidx_t * t8_offset_to)
{
  t8_gloidx_t         proca_first, proca_last;
  t8_gloidx_t         procb_first, procb_last;
  int                 keeps_first;

  T8_ASSERT (t8_offset_from != NULL && t8_offset_to != NULL);
  /* proca sends to procb if proca's first tree (plus 1 if it is shared)
   * is smaller than procb's last tree and
   * proca's last tree is bigger than procb's first tree
   * and proca has trees to send */
  /*  true if procb will keep its first tree and it is shared */
  keeps_first = t8_offset_from[procb] < 0 &&
      t8_offset_in_range (t8_offset_first (procb, t8_offset_from),
                          procb, t8_offset_to)
      && !t8_offset_empty (procb, t8_offset_from);
  if (proca == procb && keeps_first) {
    /* If proca = procb and the first tree is kept, we definetely send */
    return 1;
  }
  proca_first = t8_offset_first (proca, t8_offset_from) +
    (t8_offset_from[proca] < 0);
  proca_last = t8_offset_last (proca, t8_offset_from);
  procb_first = t8_offset_first (procb, t8_offset_to);
  procb_last = t8_offset_last (procb, t8_offset_to);
  if (keeps_first && proca_last == t8_offset_first (procb, t8_offset_from)) {
    proca_last--;
  }
  if (proca_first <= proca_last &&      /* There are trees to send  and... */
      proca_first <= procb_last       /* The first tree on a before is smaller than
                                       * the last on b after partitioning and... */
      && proca_last >= procb_first    /* The last tree on before is bigger than */
                    + (keeps_first    /* the first on b after partitioning */
                      && procb_first == t8_offset_first(procb, t8_offset_from))

    ) {

    return 1;
  }
  /* We also have to cover the case where a process sends the first tree to
   * itself */
  if (proca == procb &&
      t8_offset_in_range (t8_offset_first (proca, t8_offset_from), procb,
                          t8_offset_to)
      && !t8_offset_empty (proca, t8_offset_from)) {
    return 1;
  }
  return 0;
}

/* Determine whether a process A will send a tree T to a process B.
 */
static int
t8_offset_sendstree (int proc_send, int proc_to, t8_gloidx_t gtree,
                     t8_gloidx_t * offset_from, t8_gloidx_t * offset_to)
{
  /* If the tree is not the last tree on proc_send it is send if
   * first_tree <= tree <= last_tree on the new partition for process proc_to.
   * If the tree is the last tree is is send if the same condition holds and
   * it is not already the first tree of proc_to in the old partition */
  if (!t8_offset_in_range (gtree, proc_send, offset_from)) {
    /* The tree is not in proc_send's part of the partition */
      return 0;
  }

  if (!t8_offset_in_range (gtree, proc_to, offset_to)) {
    /* The tree will not be in proc_to's part of the new partition */
    return 0;
  }
  if (!t8_offset_empty (proc_to, offset_from) &&
      gtree == t8_offset_first (proc_to, offset_from) &&
      proc_send != proc_to) {
    /* The tree is already a tree of proc_to. This catches the case where
     * tree is shared between proc_send and proc_to. */
    return 0;
  }
  if (proc_send != proc_to &&
      gtree == t8_offset_first (proc_send, offset_from) &&
      offset_from[proc_send] < 0) {
    /* proc_send is not proc_to and tree is the first of proc send and shared
     * then a process smaller then proc_send will send tree */
    return 0;
  }
  /* The tree is on proc_send and will be on proc_to.
   * If proc_send is not proc_to then
   * tree is not shared by a smaller process than proc_send. */
  return 1;
}

/* Count the number of sending procs from start to end sending to mpirank
 * A process counts as sending if it has at least one non-shared local tree */
static int
t8_offset_range_send (int start, int end, int mpirank,
                      t8_gloidx_t * offset_from, t8_gloidx_t * offset_to)
{
  int                 count = 0, i;

  for (i = start; i <= end; i++) {
    if (t8_offset_sendsto (i, mpirank, offset_from, offset_to)) {
      count++;
    }
  }
  return count;
}

#if 0
/* TODO: Do we need this function? */
/* Determine whether a given global tree id is in the range of a given
 * process without considering the first tree if it is shared */
static int
t8_offset_in_range_wofirstshared (t8_gloidx_t tree_id, int proc,
                                  t8_gloidx_t * offset)
{
  return t8_offset_first (proc, offset) + (offset[proc] < 0) <= tree_id
    && tree_id <= t8_offset_last (proc, offset);
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
  }
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process */
void
t8_cmesh_gather_treecount (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_gloidx_t         tree_offset;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  T8_ASSERT (cmesh->tree_offsets == NULL);

  tree_offset = cmesh->first_tree_shared ? -cmesh->first_tree - 1 :
    cmesh->first_tree;
  cmesh->tree_offsets = t8_cmesh_alloc_offsets (cmesh->mpisize, comm);
  t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX,
                            cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
  t8_shmem_array_set_gloidx (cmesh->tree_offsets, cmesh->mpisize,
                             cmesh->num_trees);
}

void
t8_offset_print (t8_cmesh_t cmesh, sc_MPI_Comm comm){
#if T8_ENABLE_DEBUG
  char                buf[BUFSIZ] = "| ";
  int                 i, offset_isnew = 0;

  if (cmesh->tree_offsets == NULL) {
    t8_cmesh_gather_treecount (cmesh, comm);
    offset_isnew = 1;
  }
  for (i = 0; i <= cmesh->mpisize; i++) {
    snprintf (buf + strlen (buf), BUFSIZ - strlen (buf), " % lli |",
              (long long) t8_shmem_array_get_gloidx (cmesh->tree_offsets,
                                                     i));
  }
  if (offset_isnew == 1) {
    t8_shmem_array_destroy (&cmesh->tree_offsets);
    T8_ASSERT (cmesh->tree_offsets == NULL);
  }
  t8_debugf ("Offsets = %s\n", buf);
#endif
}

/* Return a process that a given process definitely sends to/receives from */
static int
t8_cmesh_partition_send_any (int proc, t8_cmesh_t cmesh_from,
                             t8_gloidx_t * offset_from,
                             t8_gloidx_t * offset_to, int receive)
{
  int                 lookhere, range[2], *sender, *receiver;
  int                 search_dir;
  int                 last;


  T8_ASSERT (proc >= 0 && proc < cmesh_from->mpisize);


  if ((!receive && t8_offset_nosend (proc, cmesh_from->mpisize, offset_from,
                                     offset_to)) ||
      (receive && t8_offset_empty (proc, receive ? offset_to : offset_from))) {
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
  while (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
    last = lookhere;
    while ((receive && t8_offset_nosend (lookhere, cmesh_from->mpisize,
                                          offset_from, offset_to))
            || (!receive && t8_offset_empty (lookhere, offset_to))) {
      /* Skip processes that do not send/recv anything */
      lookhere += search_dir;
      if (lookhere < 0) {
        /* There is no candidate to the left of last */
        range[0] = last + 1;
        lookhere = (range[0] + range[1])/2;
        last = lookhere;
        search_dir = 1;
      }
      else if (lookhere >= cmesh_from->mpisize) {
        /* There is no candidate to the right of last */
        range[1] = last - 1;
        lookhere = (range[0] + range[1])/2;
        last = lookhere;
        search_dir = -1;
      }
    }
    if (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)) {
      if (t8_offset_last (*sender, offset_from)
          < t8_offset_first (*receiver, offset_to)
          + (offset_to[*receiver] < 0 && *sender != *receiver
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
    else { /* t8_offset_sendsto */
        range[0] = range[1] = lookhere;
      }
    if (range[0] <= range[1]) {
      lookhere = (range[0] + range[1]) / 2;
    }
    t8_debugf ("%i\n ", lookhere);
    T8_ASSERT (0 <= lookhere && lookhere < cmesh_from->mpisize);
#if 0
    if (last == lookhere) {
      /* We did not change in this round which means, that we found the rank */
      return lookhere;
    }
#endif
  }

  return lookhere;
}

/* Compute first and last process to which we will send data */
/* Returns the local tree_id of last local tree on send_first. */
static              t8_locidx_t
t8_cmesh_partition_sendrange (t8_cmesh_t cmesh_to, t8_cmesh_t cmesh_from,
                              int *send_first, int *send_last, int receive)
{
  t8_gloidx_t         ret;
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
  /* We start be probing our first_guess */
  lookhere = first_guess;
  while (*send_first == -1) {
    /* However it could be empty in this case we search linearly to
     * find the next nonempty receiver in search direction */
    /* TODO: replace nonempty with nonsending? */
    while (1 <= lookhere && lookhere < cmesh_from->mpisize - 1 &&
           (!t8_offset_sendsto (*sender, *receiver, offset_from, offset_to)
            /* TODO: The two next lines can be removed?*/
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

  last = -2; /* If this does not change we found our process */
  while (*send_last == -1) {
    /* Our guess could be empty in this case we search linearly to
     * find the next nonempty receiver in search direction */
    /* TODO: replace nonempty with nonsending? */
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

  /* Calculate the last local tree that we need to send to send_first */
  /* Set it to the last tree on send_first */
  ret =  t8_offset_last (*send_first, offset_to) -
      cmesh_from->first_tree;
  /* If there are actually more trees on send_first than we have, we need to send
   * all our local trees to send_first */
  ret = SC_MIN (ret, cmesh_from->num_local_trees - 1);
  if (cmesh_from->mpirank != *send_first &&
      t8_offset_in_range (t8_offset_last (cmesh_from->mpirank, offset_from),
                          *send_first, offset_from)) {
    /* Substract one if our last tree already belonged to send_first */
    ret--;
  }

  t8_debugf ("%s_first = %i, %s_last = %i, last_tree = %li\n",
             receive ? "recv" : "send", *send_first,
             receive ? "recv" : "send", *send_last, ret);

  T8_ASSERT (*send_first >= 0);
//TODO:reactivate  T8_ASSERT (*send_last >= 0);
  T8_ASSERT (receive || (ret >= 0 && ret < cmesh_from->num_local_trees));
  T8_ASSERT (ret == (t8_locidx_t) ret);
  return (t8_locidx_t) ret;
}

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
      +((4 - t8_eclass_num_faces[eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t)) % 4) % 4);  /* padding */
  }
  return ghost_neighbor_bytes;
}

/* Determine whether a local tree or ghost should be send to process p as a ghost.
 * This is the case if and only if:
 *  - tree will not be a local tree on p
 * and
 *  - we are the smallest rank under all procs sending a tree to p that
 *    has this tree as ghost or local tree. */
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
                                                     from_offsets, offset_to)) {
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
                              size_t ghost_neighbor_bytes,
                              size_t tree_neighbor_bytes,
                              sc_array_t * send_as_ghost,
                              t8_locidx_t send_first, t8_locidx_t send_last,
                              size_t total_alloc, int to_proc)
{
  t8_ctree_t          tree, tree_cpy;
  size_t              temp_offset_tree, temp_offset_att, iz, temp_offset,
    temp_offset_data, last_offset, last_num_att, last_size;
  //ssize_t             last_attribute_diff;
  t8_attribute_info_struct_t *attr_info;
  t8_locidx_t         num_ghost_send = send_as_ghost->elem_count;
  t8_locidx_t        *face_neighbor, ghost_id, itree;
  t8_gloidx_t        *face_neighbor_g, *face_neighbor_gnew, new_neighbor;
  t8_cghost_t         ghost, ghost_cpy;
  int                 iface;
  int8_t             *ttf_ghost, *ttf;

  /* Copy all trees to the send buffer */
  /* TODO: This is currently inefficient since we copy each tree for itself.
   *       Best practive is to copy chunks of trees out of the different part
   *       arrays of cmesh_from */
  if (total_alloc == 0 || send_buffer == NULL) {
    t8_debugf ("No data to store in buffer.\n");
    return;
  }
  temp_offset = 0;
  temp_offset_att = 0;
  temp_offset_data = 0;
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
      + ((4 - t8_eclass_num_faces[tree->eclass] *
          (sizeof (t8_locidx_t) + sizeof (int8_t)) % 4) % 4);
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
  /* Set new face_neighbor offsets */
  /* TODO: indent bug? */
  temp_offset = num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes;   /* Computes the offset of the face neighbors of the new trees */
  temp_offset_att = temp_offset + tree_neighbor_bytes;  /* Compute the offset of the new attribute infos */
  temp_offset_tree = 0;
  attr_info = (t8_attribute_info_struct_t *) (send_buffer +
                                              num_trees *
                                              sizeof (t8_ctree_struct_t) +
                                              num_ghost_send *
                                              sizeof (t8_cghost_struct_t) +
                                              ghost_neighbor_bytes +
                                              tree_neighbor_bytes);
#if 0
  /* This has to be add to each attribute info offset to calculated the new offset */
  last_attribute_diff = attr_info_bytes - attr_info->attribute_offset;
#endif
  /* Set attribute offsets of trees and attribute data offsets of info objects */
  tree_cpy = NULL;
  last_num_att = 0;
  last_size = 0;
  last_offset = attr_info_bytes;

  /* TODO: The last changes here did not make it better */
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
      + ((4 - t8_eclass_num_faces[tree_cpy->eclass] *
          (sizeof (t8_locidx_t) + sizeof (int8_t)) % 4) % 4);

    /* new attribute offset for tree */
    tree_cpy->att_offset = temp_offset_att - temp_offset_tree;
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
    temp_offset_tree += sizeof (t8_ctree_struct_t);
    last_num_att = tree_cpy->num_attributes;
  }

  /* Copy all ghosts and set their face entries and offsets */
  /* Offset of ghost face_neighbor from first ghost */
  temp_offset = num_ghost_send * sizeof (t8_cghost_struct_t);
  /* Offset of current ghost from first ghost */
  temp_offset_tree = 0;
  for (iz = 0; iz < send_as_ghost->elem_count; iz++) {
    ghost_id = *((t8_locidx_t *) sc_array_index (send_as_ghost, iz));
    ghost_cpy = (t8_cghost_t) (send_buffer +
                               num_trees * sizeof (t8_ctree_struct_t)
                               + iz * sizeof (t8_cghost_struct_t));
    /* Get the correc element class from the ghost_id.
     * For this we have to check whether ghost_id points to a local tree or ghost */
    ghost_cpy->eclass = ghost_id < cmesh_from->num_local_trees ?
      t8_cmesh_get_tree_class ((t8_cmesh_t) cmesh_from, ghost_id) :
      t8_cmesh_get_ghost_class ((t8_cmesh_t) cmesh_from, ghost_id
                                - cmesh_from->num_local_trees);
    ghost_cpy->neigh_offset = temp_offset - temp_offset_tree;   /* New face neighbor offset */
    face_neighbor_gnew = (t8_gloidx_t *) T8_GHOST_FACE (ghost_cpy);
    ttf_ghost = T8_GHOST_TTF (ghost_cpy);
    if (ghost_id >= cmesh_from->num_local_trees) {
      /* The ghost that we send was a local ghost */
      ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees,
                                            ghost_id -
                                            cmesh_from->num_local_trees,
                                            &face_neighbor_g, &ttf);
      ghost_cpy->eclass = ghost->eclass;
      ghost_cpy->treeid = ghost->treeid;
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
      ghost_cpy->eclass = tree->eclass;
      ghost_cpy->treeid = ghost_id + cmesh_from->first_tree;
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
    }
    /* compute new offsets */
    temp_offset += t8_eclass_num_faces[ghost_cpy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))    /* offset */
      +((4 - t8_eclass_num_faces[ghost_cpy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t)) % 4) % 4);       /* padding */
    temp_offset_tree += sizeof (t8_cghost_struct_t);
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

  /* loop over all trees that will be send */
  for (itree = range_start; itree <= range_end; itree++) {
    T8_ASSERT ((!cmesh_from->set_partition && iproc == cmesh_from->mpirank) ||
    t8_offset_sendstree (cmesh_from->mpirank, iproc, itree + cmesh_from->first_tree,
                                  t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets),
                                  t8_shmem_array_get_gloidx_array (cmesh->tree_offsets)));
    tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree,
                                        &face_neighbor, &ttf);
    /* Count the additional memory needed per tree from neighbors */
    *tree_neighbor_bytes += t8_eclass_num_faces[tree->eclass] *
      (sizeof (*face_neighbor) + sizeof (*ttf));
    /* TODO: change padding to sizeof (void*) */
    *tree_neighbor_bytes += (4 - *tree_neighbor_bytes % 4) % 4; /* padding to make number of bytes per tree
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
                             int *my_buffer_bytes, sc_MPI_Request ** requests,
                             sc_MPI_Comm comm)
{
  size_t              attr_bytes = 0, tree_neighbor_bytes,
    ghost_neighbor_bytes, total_alloc, attr_info_bytes;
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

  range_end = t8_cmesh_partition_sendrange (cmesh, (t8_cmesh_t) cmesh_from,
                                            send_first, send_last, 0);

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

  for (iproc = *send_first; iproc <= *send_last; iproc++) {
    /* Since we do not allocate a request and buffer for my rank if it is included in
     * the sendrange, we have to offset the index into these array by -1 when the
     * current irpoc is bigger than my rank */
    if (*send_first <= cmesh->mpirank && iproc > cmesh->mpirank) {
      flag = 1;
    }

    t8_debugf ("send to %i [%i,%i]\n", iproc, range_start, range_end);
    t8_debugf ("Start send loop\n");
    attr_bytes = 0;
    tree_neighbor_bytes = 0;
    attr_info_bytes = 0;
    ghost_neighbor_bytes = 0;

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

#ifdef T8_ENABLE_DEBUG
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
    /* Total number of bytes that we send to the other process */
    total_alloc = num_trees * sizeof (t8_ctree_struct_t) +
      num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes +
      tree_neighbor_bytes + attr_info_bytes + attr_bytes;
    /* Extra space to store the number of trees and ghosts in the send buffer */
    total_alloc += 2 * sizeof (t8_locidx_t);

    t8_debugf ("NT %i -- each %zdB\n", num_trees, sizeof (t8_ctree_struct_t));
    t8_debugf ("NGS %i - each %zdB\n", num_ghost_send,
               sizeof (t8_cghost_struct_t));
    t8_debugf ("GNB %zd\n", ghost_neighbor_bytes);
    t8_debugf ("TNB %zd\n", tree_neighbor_bytes);
    t8_debugf ("AIB %zd\n", attr_info_bytes);
    t8_debugf ("AB %zd\n", attr_bytes);
    t8_debugf ("Ta %zd\n", total_alloc);
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
                                                   && iproc == cmesh_from->mpirank) ||
               t8_offset_sendsto (cmesh->mpirank, iproc, offset_from,
                                  offset_to));

    /* Copy all data to the send buffer */
    t8_cmesh_partition_copy_data (buffer, cmesh,
                                  cmesh_from, num_trees, attr_info_bytes,
                                  ghost_neighbor_bytes,
                                  tree_neighbor_bytes, &send_as_ghost,
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
        *(*requests + iproc - flag - *send_first) = sc_MPI_REQUEST_NULL;
      }
    }
    if (iproc < *send_last) {
      /* compute new ranges here */
      /* If tree_offset of iproc + 1 is < 0 the last tree is shared and has to
       * be send to the next process as well, except when this process already
       * has the last tree */
      range_start = range_end + 1 - (offset_to[iproc + 1] < 0);
      /* add number of trees in iproc + 1 to send to range_end */
      /* We have to be careful with locidx overflow when we go out of bounds
       * of our process */

      range_end =
        t8_glo_min (range_end + t8_offset_num_trees (iproc + 1, offset_to)
                    - (offset_to[iproc + 1] < 0),
                    cmesh_from->num_local_trees - 1);
      t8_debugf ("RE: %i\n", range_end);
      /* If the last tree was already present on the receiving process
       * we exclude it from sending */
      if (cmesh_from->mpirank != iproc + 1 &&
          range_end == cmesh_from->num_local_trees - 1
          && !t8_offset_empty (iproc + 1, offset_from)
          && t8_offset_first (iproc + 1, offset_from)
          == cmesh_from->first_tree + cmesh_from->num_local_trees - 1) {
        range_end--;
        t8_debugf ("RE: %i\n", range_end);
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
static void
t8_cmesh_partition_recvloop (t8_cmesh_t cmesh,
                             const struct t8_cmesh *cmesh_from,
                             t8_gloidx_t * tree_offset,
                             char *my_buffer, int my_buffer_bytes,
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
  sc_list_t          *possible_receivers;       /* Linked list storing the ranks from which
                                                   we still expect a message */
  int                *poss_recv_data;   /* Data storage for these ranks */
  sc_link_t          *iterate;  /* Iterator through the linked list */
  sc_link_t          *prev;     /* We need to store iterates predecessor in order to
                                   be able to remove iterate when needed */
  int                 iprobe_flag;

  num_trees = t8_offset_num_trees (cmesh->mpirank, tree_offset);
  /* Receive from other processes */
  if (cmesh_from->set_partition) {
    T8_ASSERT (cmesh_from->tree_offsets != NULL);
    from_offsets = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
    t8_cmesh_partition_recvrange (cmesh, (t8_cmesh_t) cmesh_from,
                                  &recv_first, &recv_last);
    if (recv_first <= recv_last) {
      T8_ASSERT (fr == recv_first);
      T8_ASSERT (lr == recv_last);
    }

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

    /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for an message for each
     * sender and if there is one we receive it and remove the sender from
     * the list.
     * The last message can be reveived via probe */
    while (possible_receivers->elem_count > 1) {
      iprobe_flag = 0;
      prev = NULL;
      for (iterate = possible_receivers->first;
           iterate != NULL && iprobe_flag == 0;) {
        /* Iterate through all entries and probe for message */
        proc_recv = *(int *) iterate->data;
        t8_debugf ("Probing for message from %i\n", proc_recv);
        mpiret = sc_MPI_Iprobe (proc_recv, T8_MPI_PARTITION_CMESH, comm,
                                &iprobe_flag, &status);
        SC_CHECK_MPI (mpiret);
        if (iprobe_flag == 0) {
          /* We do only continue if there is no message to receive */
          prev = iterate;
          iterate = iterate->next;
        }
      }
      if (iprobe_flag != 0) {
        /* There is a message to receive */
        T8_ASSERT (proc_recv == status.MPI_SOURCE);
        T8_ASSERT (status.MPI_TAG == T8_MPI_PARTITION_CMESH);
        T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last &&
                   t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets,
                                      tree_offset));
        t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status,
                                            local_procid, recv_first,
                                            &num_ghosts);
        sc_list_remove (possible_receivers, prev);
      }
    }

    /* We have one message left and does we do not need to IProbe
     * but call the blocking Probe instead. */
    T8_ASSERT (possible_receivers->elem_count == 1);
    proc_recv = *(int *) sc_list_pop (possible_receivers);
    sc_list_destroy (possible_receivers);
    T8_FREE (poss_recv_data);
    mpiret = sc_MPI_Probe (proc_recv, T8_MPI_PARTITION_CMESH, comm, &status);
    SC_CHECK_MPI (mpiret);
    t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status,
                                        local_procid, recv_first,
                                        &num_ghosts);
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
  int                 iproc, my_buffer_bytes, num_send_mpi, mpiret;
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
  T8_ASSERT (!cmesh_from->set_partition || send_first == -1 || send_first == fs);
  T8_ASSERT (!cmesh_from->set_partition || send_last == -2 || send_last == ls);

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
      cmesh->num_trees_per_eclass[tree->eclass]++;
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
 * This includes partitioning by uiniform level and partitioning from a second cmesh */
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
    /* To compute the tree_offsets correctly we have to invert the sign on the
     * first tree if last tree is shared, since we use it for MPI to allgather the tree_offsets */
    if (cmesh->first_tree_shared) {
      cmesh->first_tree = -cmesh->first_tree - 1;
    }
    /* allocate and fill tree_offset array with number of trees per process */
    /* We have to be careful with shared memory in the case where cmesh comm is
     * duplicated, since we have to use the same communicator for allocating and freeing memory.
     * Thus this function must only be called after cmesh communicator was duplicated,
     * so we check whether mpisize has been set */
    T8_ASSERT (cmesh->mpisize > 0);
    t8_shmem_array_init (&cmesh->tree_offsets, sizeof (t8_gloidx_t),
                         cmesh->mpisize + 1, comm);
    t8_shmem_array_allgather (&cmesh->first_tree, 1, T8_MPI_GLOIDX,
                              cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
    t8_shmem_array_set_gloidx (cmesh->tree_offsets, cmesh->mpisize,
                               cmesh_from->num_trees);
    /* tree_offsets was computed, reinvert the sign */
    if (cmesh->first_tree_shared) {
      T8_ASSERT (cmesh->first_tree <= 0);
      cmesh->first_tree = -cmesh->first_tree - 1;
    }
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
  t8_offset_print (cmesh->set_from, comm);
  t8_debugf ("To:\n");
  t8_offset_print (cmesh, comm);
  /***************************************************/
  /*        Done with local num and tree_offset      */
  /***************************************************/
  t8_cmesh_partition_given (cmesh, cmesh->set_from, tree_offsets, comm);
  t8_global_productionf ("Done cmesh partition\n");
}

static int
t8_offset_consistent (int mpisize, t8_gloidx_t * offset,
                      t8_gloidx_t num_trees)
{
  int                 i, ret = 1;
  t8_gloidx_t         temp;

  ret = offset[0] == 0;
  temp = 0;
  for (i = 1; i < mpisize && ret; i++) {
    if (offset[i] < 0) {
      ret &= (fabs (offset[i] + 1) >= temp);
      temp = fabs (offset[i] + 1);
    }
    else {
      ret &= (offset[i] >= temp);
      temp = offset[i];
    }
    ret &= (temp <= num_trees);
  }
  ret &= (offset[mpisize] == num_trees);
  t8_debugf ("Offset is %s\n", ret ? "consistend." : "not consistend!");
  return ret;
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
  return shmem_array;
}

/* Create a random partition */
/* if shared is nonzero than first trees can be shared */
t8_shmem_array_t
t8_cmesh_offset_random (sc_MPI_Comm comm, t8_gloidx_t num_trees, int shared,
                        int seed)
{
  int                 iproc, mpisize, mpiret, random_number, mpirank;
  int                 first_shared;
  int                 i;
  unsigned            u_seed;
  t8_gloidx_t        *offsets;
  t8_shmem_array_t    shmem_array;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  shmem_array = t8_cmesh_alloc_offsets (mpisize, comm);

  offsets = t8_shmem_array_get_gloidx_array (shmem_array);

  if (seed < 0) {
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
    if ((int) (num_trees * 2. / mpisize) == 0) {
      /* This case prevents division by 0 */
      random_number = 1;
    }
    else {
      random_number = rand_r (&u_seed) % (int) (num_trees * 2. / mpisize);
    }

    random_number += first_shared;
    /* If we would exceed the number of trees we cut the random number */
    if (t8_offset_first (iproc - 1, offsets) + random_number > num_trees) {
      random_number = num_trees - t8_offset_first (iproc - 1, offsets);
    }
    if (shared) {
      first_shared = rand_r (&u_seed) % 2;
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
