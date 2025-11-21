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

/** \file t8_cmesh_offset.c
 *
 * TODO: document this file
 */

#include "t8_cmesh_offset.h"

/* Return -1 if A is smaller 0
 * Return 0 if not. */
static int
t8_glo_kl0 (const t8_gloidx_t A)
{
  return A < 0 ? -1 : 0;
}

/* The first tree of a given process in a partition */
t8_gloidx_t
t8_offset_first (const int proc, const t8_gloidx_t *offset)
{
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);
  return T8_GLOIDX_ABS (offset[proc]) + t8_glo_kl0 (offset[proc]);
}

/* Given the global tree id of the first local tree of a process and
 * the flag whether it is shared or not, compute the entry in the offset array. */
/* This entry is the first_tree if it is not shared and
 * -first_tree - 1  if it is shared.
 */
t8_gloidx_t
t8_offset_first_tree_to_entry (const t8_gloidx_t first_tree, const int shared)
{
  return shared ? -first_tree - 1 : first_tree;
}

/* The number of trees of a given process in a partition */
t8_gloidx_t
t8_offset_num_trees (const int proc, const t8_gloidx_t *offset)
{
  t8_gloidx_t num_global_trees;
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);

  num_global_trees = T8_GLOIDX_ABS (offset[proc + 1]) - t8_offset_first (proc, offset);

  /* Fails if both offset[proc] and offset[proc+1] are -1 */
  T8_ASSERT (num_global_trees >= 0);

  return num_global_trees;
}

/* The last local tree of a given process in a partition */
t8_gloidx_t
t8_offset_last (const int proc, const t8_gloidx_t *offset)
{
  T8_ASSERT (proc >= -1);
  T8_ASSERT (offset != NULL);

  return T8_GLOIDX_ABS (offset[proc + 1]) - 1;
}

#if T8_ENABLE_DEBUG
/* Query whether a given global tree is in a valid range of a partition */
static int
t8_offset_valid_tree (const t8_gloidx_t gtree, const int mpisize, const t8_gloidx_t *offset)
{
  T8_ASSERT (offset != NULL);

  return 0 <= gtree && gtree <= t8_offset_last (mpisize - 1, offset);
}
#endif

/* Return 1 if the process has no trees in the partition.
 * Return 0 if the process has at least one tree */
int
t8_offset_empty (const int proc, const t8_gloidx_t *offset)
{
  T8_ASSERT (proc >= 0);
  T8_ASSERT (offset != NULL);

  if (t8_offset_num_trees (proc, offset) <= 0) {
    return 1;
  }
  return 0;
}

/* Find the next higher rank that is not empty.
 * returns mpisize if this rank does not exist. */
int
t8_offset_next_nonempty_rank (const int rank, const int mpisize, const t8_gloidx_t *offset)
{
  int next_nonempty = rank + 1;
  while (next_nonempty < mpisize && t8_offset_empty (next_nonempty, offset)) {
    next_nonempty++;
  }
  return next_nonempty;
}

#if T8_ENABLE_DEBUG
/* Check whether a given offset array represents a valid
 * partition. */
int
t8_offset_consistent (const int mpisize, const t8_shmem_array_t offset_shmem, const t8_gloidx_t num_trees)
{
  int i, ret = 1;
  t8_gloidx_t last_tree;
  T8_ASSERT (t8_shmem_array_is_initialized (offset_shmem));

  const t8_gloidx_t *offset = t8_shmem_array_get_gloidx_array (offset_shmem);

  ret = offset[0] == 0;
  last_tree = t8_offset_last (0, offset); /* stores the last tree of process i-1 */
  for (i = 1; i < mpisize && ret; i++) {
    if (t8_offset_empty (i, offset)) {
      /* If the process is empty, then its first tree must not be shared */
      ret &= offset[i] >= 0;
    }
    else {
      /* If the process is not empty its first local tree must be bigger or
       * equal to the last tree of the last nonempty process.
       * Equality must only hold, when the first tree is shared, thus offset[i] < 0 */
      if (offset[i] < 0) {
        ret &= t8_offset_first (i, offset) == last_tree;
      }
      else {
        ret &= t8_offset_first (i, offset) > last_tree;
      }
      last_tree = t8_offset_last (i, offset);
      ret &= (last_tree <= num_trees);
    }
  }
  ret &= (offset[mpisize] == num_trees);
  t8_debugf ("Offset is %s %i\n", ret ? "consistent." : "not consistent!", i - 1);
  return ret;
}
#endif

/* Determine whether a given global tree id is in the range of a given process */
int
t8_offset_in_range (const t8_gloidx_t tree_id, const int proc, const t8_gloidx_t *offset)
{
  return t8_offset_first (proc, offset) <= tree_id && tree_id <= t8_offset_last (proc, offset);
}

int
t8_offset_any_owner_of_tree_ext (const int mpisize, const int start_proc, const t8_gloidx_t gtree,
                                 const t8_gloidx_t *offset)
{
  int proc = start_proc;
  int range[2] = { 0, mpisize - 1 };

  range[0] = 0;
  range[1] = mpisize - 1;
  /* find any process that owns the tree with a binary search */
  int found = 0;
  while (!found) {
    if (t8_offset_in_range (gtree, proc, offset)) {
      found = 1;
      break;
    }
    else if (t8_offset_last (proc, offset) < gtree) {
      /* look further right */
      range[0] = proc + 1;
    }
    else {
      range[1] = proc - 1;
    }
    proc = (range[0] + range[1]) / 2;
  }
  return proc;
}
/* Find any owner of a given tree.
 */
int
t8_offset_any_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset)
{
  return t8_offset_any_owner_of_tree_ext (mpisize, (mpisize - 1) / 2, gtree, offset);
}

/* Find the smallest process that owns a given tree.
 * To increase the runtime, some_owner can be a process that
 * already owns the tree. Otherwise (some_owner < 0), the function will compute one. */
int
t8_offset_first_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, int *some_owner)
{
  int proc, proc_temp;

  T8_ASSERT (t8_offset_valid_tree (gtree, mpisize, offset));
  if (*some_owner < 0) {
    *some_owner = t8_offset_any_owner_of_tree (mpisize, gtree, offset);
  }
  T8_ASSERT (*some_owner < mpisize);
  proc = *some_owner;
  /* Here, proc stores an arbitrary owner of gtree */
  /* Now we find the smallest process that owns the tree */
  proc_temp = proc;
  while (proc_temp >= 0 && t8_offset_in_range (gtree, proc_temp, offset)) {
    /* decrement temp as long as it holds the tree */
    proc_temp--;
    while (0 <= proc_temp && t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp--;
    }
  }
  /* We are now one nonempty process below the smallest process having the tree */
  if (proc_temp >= -1) {
    proc_temp++;
    while (t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp++;
    }
    T8_ASSERT (t8_offset_in_range (gtree, proc_temp, offset));
  }
  else {
    /* This should never happen */
    SC_ABORT ("ERROR: proc_temp ran out of bounds");
  }
  proc = proc_temp;
  return proc;
}

static int
t8_offset_next_prev_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset,
                                   const int current_owner, const int search_dir)
{
  int proc;

  int search_left_or_right = search_dir > 0 ? 1 : -1; /* Adjust search dir. Positive means next process,
                                           negative previous. */
  proc = current_owner + search_left_or_right;
  while (proc >= 0 && proc < mpisize && t8_offset_empty (proc, offset)) {
    /* Skip empty processes */
    proc += search_left_or_right;
  }
  if (proc >= 0 && proc < mpisize && t8_offset_in_range (gtree, proc, offset)) {
    /* proc is still in the range and it owns gtree */
    return proc;
  }
  /* Either proc is greater than mpisize or smaller 0 or it does not
   * own gtree any more. */
  return -1;
}

/* Given a process current_owner that has the local tree gtree,
 * find the next bigger rank that also has this tree.
 * If none is found, we return -1.
 */
int
t8_offset_next_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset,
                              const int current_owner)
{
  return t8_offset_next_prev_owner_of_tree (mpisize, gtree, offset, current_owner, +1);
}

/* Given a process current_owner that has the local tree gtree,
 * find the next smaller rank that also has this tree.
 * If none is found, we return -1.
 */
int
t8_offset_prev_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset,
                              const int current_owner)
{
  return t8_offset_next_prev_owner_of_tree (mpisize, gtree, offset, current_owner, -1);
}

/* Find the biggest process that owns a given tree.
 * To increase the runtime, some_owner can be a process that
 * already owns the tree. Otherwise (some_owner < 0), the function will compute one. */
int
t8_offset_last_owner_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, int *some_owner)
{
  int proc, proc_temp;

  T8_ASSERT (t8_offset_valid_tree (gtree, mpisize, offset));
  if (*some_owner < 0) {
    *some_owner = t8_offset_any_owner_of_tree (mpisize, gtree, offset);
  }
  T8_ASSERT (*some_owner < mpisize);
  proc = *some_owner;
  /* Here, proc stores an arbitrary owner of gtree */
  /* Now we find the smallest process that owns the tree */
  proc_temp = proc;
  while (proc_temp < mpisize && t8_offset_in_range (gtree, proc_temp, offset)) {
    /* increment temp as long as it holds the tree */
    proc_temp++;
    while (proc_temp < mpisize && t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp++;
    }
  }
  /* We are now one nonempty process above the biggest process having the tree */
  if (proc_temp <= mpisize) {
    proc_temp--;
    while (t8_offset_empty (proc_temp, offset)) {
      /* Skip empty processes */
      proc_temp--;
    }
    T8_ASSERT (t8_offset_in_range (gtree, proc_temp, offset));
  }
  else {
    /* This should never happen */
    SC_ABORT ("ERROR: proc_temp ran out of bounds");
  }
  proc = proc_temp;
  return proc;
}

/* Compute a list of all processes that own a specific tree */
/* Owners must be an initialized sc_array with int elements and
 * element count 0 */
void
t8_offset_all_owners_of_tree (const int mpisize, const t8_gloidx_t gtree, const t8_gloidx_t *offset, sc_array_t *owners)
{
  int proc;
  int *entry;
  int some_owner = -1;
  T8_ASSERT (owners != NULL);
  T8_ASSERT (owners->elem_count == 0);
  T8_ASSERT (owners->elem_size == sizeof (int));
  T8_ASSERT (t8_offset_valid_tree (gtree, mpisize, offset));

  /* Find the smallest process that has the tree */
  proc = t8_offset_first_owner_of_tree (mpisize, gtree, offset, &some_owner);
  /* Add the first process to the array */
  entry = (int *) sc_array_push (owners);
  *entry = proc;
  /* Now we parse through all processes bigger than the first one until
   * they do not have the tree anymore. */
  while (proc >= 0) {
    /* Find the next owner of gtree */
    proc = t8_offset_next_owner_of_tree (mpisize, gtree, offset, proc);
    if (proc >= 0) {
      /* If we found one, we add it to the array */
      /* Otherwise we found all owners */
      entry = (int *) sc_array_push (owners);
      *entry = proc;
    }
  }
}

/* Return 1 if the process will not send any trees, that is if it is
 * empty or has only one shared tree. */
int
t8_offset_nosend (const int proc, const int mpisize, const t8_gloidx_t *offset_from, const t8_gloidx_t *offset_to)
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
    first_not_send
      = offset_from[proc] < 0 && !t8_offset_in_range (t8_offset_first (proc, offset_from), proc, offset_to);

    if ((first_not_send && num_trees == 2) || (!first_not_send && num_trees == 1)) {
      int temp_proc;
      size_t iz;
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
      if (temp_proc >= mpisize || t8_offset_first (temp_proc, offset_from) != last_tree) {
        /* Our last tree is unique and thus we need to send it */
        return 0;
      }
      sc_array_init (&receivers, sizeof (int));
      /* Now we need to find out all process in the new partition that have last_tree
       * and check if any of them did not have it before. */
      t8_offset_all_owners_of_tree (mpisize, last_tree, offset_to, &receivers);
      for (iz = 0; iz < receivers.elem_count; iz++) {
        temp_proc = *(int *) sc_array_index (&receivers, iz);
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
int
t8_offset_sendsto (const int proca, const int procb, const t8_gloidx_t *t8_offset_from, const t8_gloidx_t *t8_offset_to)
{
  t8_gloidx_t proca_first, proca_last;
  t8_gloidx_t procb_first, procb_last;
  int keeps_first;

  T8_ASSERT (t8_offset_from != NULL && t8_offset_to != NULL);
  /* proca sends to procb if proca's first tree (plus 1 if it is shared)
   * is smaller than procb's last tree and
   * proca's last tree is bigger than procb's first tree
   * and proca has trees to send */
  /*  true if procb will keep its first tree and it is shared */
  if (t8_offset_empty (proca, t8_offset_from) || t8_offset_empty (procb, t8_offset_to)) {
    /* If the sender or the receiver is empty we do not send anything. */
    return 0;
  }
  keeps_first = t8_offset_from[procb] < 0
                && t8_offset_in_range (t8_offset_first (procb, t8_offset_from), procb, t8_offset_to)
                && !t8_offset_empty (procb, t8_offset_from);
  if (proca == procb && keeps_first) {
    /* If proca = procb and the first tree is kept, we definitely send */
    return 1;
  }
  proca_first = t8_offset_first (proca, t8_offset_from) + (t8_offset_from[proca] < 0);
  proca_last = t8_offset_last (proca, t8_offset_from);
  procb_first = t8_offset_first (procb, t8_offset_to);
  procb_last = t8_offset_last (procb, t8_offset_to);
  if (keeps_first && proca_last == t8_offset_first (procb, t8_offset_from)) {
    proca_last--;
  }
  if (proca_first <= proca_last &&      /* There are trees to send  and... */
      proca_first <= procb_last         /* The first tree on a before is smaller than
                                 * the last on b after partitioning and... */
      && proca_last >= procb_first      /* The last tree on a before is bigger than */
                         + (keeps_first /* the first on b after partitioning */
                            && procb_first == t8_offset_first (procb, t8_offset_from))) {
    return 1;
  }
  return 0;
}

/* Determine whether a process A will send a tree T to a process B.
 */
int
t8_offset_sendstree (const int proc_send, const int proc_to, const t8_gloidx_t gtree, const t8_gloidx_t *offset_from,
                     const t8_gloidx_t *offset_to)
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
  if (!t8_offset_empty (proc_to, offset_from) && gtree == t8_offset_first (proc_to, offset_from)
      && proc_send != proc_to) {
    /* The tree is already a tree of proc_to. This catches the case where
     * tree is shared between proc_send and proc_to. */
    return 0;
  }
  if (proc_send != proc_to && gtree == t8_offset_first (proc_send, offset_from) && offset_from[proc_send] < 0) {
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
int
t8_offset_range_send (const int start, const int end, const int mpirank, const t8_gloidx_t *offset_from,
                      const t8_gloidx_t *offset_to)
{
  int count = 0, i;

  for (i = start; i <= end; i++) {
    if (t8_offset_sendsto (i, mpirank, offset_from, offset_to)) {
      count++;
    }
  }
  return count;
}

/**
 * Print the offsets of a partition.
 * 
 * This function prints the offsets of a partition in a debug message.
 * 
 * \param [in] offset  The offsets to print.
 * \param [in] comm    The MPI communicator to use for printing.
 */
void
t8_offset_print (__attribute__ ((unused)) const t8_shmem_array_t offset, __attribute__ ((unused)) sc_MPI_Comm comm)
{
#if T8_ENABLE_DEBUG
  char buf[BUFSIZ] = "| ";
  int i, mpiret, mpisize;

  if (offset == NULL) {
    t8_debugf ("Offsets = NULL\n");
    return;
  }

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  for (i = 0; i <= mpisize; i++) {
    snprintf (buf + strlen (buf), BUFSIZ - strlen (buf), " % lli |", (long long) t8_shmem_array_get_gloidx (offset, i));
  }
  t8_debugf ("Offsets = %s\n", buf);
#endif
}
