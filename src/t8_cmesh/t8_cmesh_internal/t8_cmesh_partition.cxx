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

/** \file t8_cmesh_partition.cxx
 * Implementation of functionality related to the partitioning of a cmesh.
 */

#include <cstring>

#include <t8_data/t8_shmem.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_element/t8_element.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_geometry/t8_geometry_handler.hxx>

/* Change the neighbor entry of a tree to match the new partition.
 * Input: A face_neighbor entry in cmesh_from and a process to which the corresponding tree will be send
 * Output: The face neighbor entry is changed to match its new id in cmesh.
 */
static void
t8_cmesh_partition_send_change_neighbor (const t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, t8_locidx_t *neighbor,
                                         const int to_proc)
{
  t8_gloidx_t temp;
  const t8_gloidx_t *tree_offset = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (0 <= *neighbor && *neighbor < cmesh_from->num_local_trees) {
    /* Neighbor is a local tree in cmesh */
    temp = cmesh_from->first_tree - t8_offset_first (to_proc, tree_offset);
    /* Assert for possible overflow du to gloidx computation */
    T8_ASSERT ((t8_locidx_t) (temp + *neighbor) == temp + *neighbor);
    *neighbor = temp + *neighbor;
  }
  else {
    t8_cghost_t ghost;
    /* neighbor is a ghost in cmesh_from */
    T8_ASSERT (*neighbor >= cmesh_from->num_local_trees
               && *neighbor < cmesh_from->num_local_trees + cmesh_from->num_ghosts);
    ghost = t8_cmesh_trees_get_ghost (cmesh_from->trees, *neighbor - cmesh_from->num_local_trees);
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
t8_partition_new_ghost_ids (const t8_cmesh_t cmesh, const t8_part_tree_t recv_part, const t8_locidx_t first_ghost,
                            const int proc)
{
  t8_locidx_t ghost_it;
  t8_cghost_t ghost;
  t8_gloidx_t *face_neighbors, tree_id_glo;
  t8_locidx_t *tree_neighbors;
  int8_t *ttf;
  int iface, face_tree;
  t8_trees_glo_lo_hash_t *new_hash;
#if T8_ENABLE_DEBUG
  int ret;
#endif

  for (ghost_it = 0; ghost_it < recv_part->num_ghosts; ghost_it++) {
    /* loop over all ghosts of recv_part */
    cmesh->trees->ghost_to_proc[ghost_it + first_ghost] = proc;
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees, first_ghost + ghost_it, &face_neighbors, &ttf);
    for (iface = 0; iface < t8_eclass_num_faces[ghost->eclass]; iface++) {
      /* loop over all faces of ghost */
      tree_id_glo = face_neighbors[iface];
      if (cmesh->first_tree <= tree_id_glo && tree_id_glo < cmesh->first_tree + cmesh->num_local_trees) {
        /* the face neighbor is a local tree */
        (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, tree_id_glo - cmesh->first_tree, &tree_neighbors, NULL);
        /* Get the number of the face of tree that is connected with ghost
         * and set the new local ghost id */
        face_tree = ttf[iface] % t8_eclass_max_num_faces[cmesh->dimension];
        tree_neighbors[face_tree] = ghost_it + first_ghost + cmesh->num_local_trees;
      }
    }
    /* Insert this ghost's global and local id into the hash table */
    new_hash = (t8_trees_glo_lo_hash_t *) sc_mempool_alloc (cmesh->trees->global_local_mempool);
    new_hash->global_id = ghost->treeid;
    /* The new local ghost id is the concurrent id of this ghost plus the
     * number of local trees */
    new_hash->local_id = ghost_it + first_ghost + cmesh->num_local_trees;
#if T8_ENABLE_DEBUG
    ret =
#endif
      sc_hash_insert_unique (cmesh->trees->ghost_globalid_to_local_id, new_hash, NULL);
    /* The entry must not have existed before */
    T8_ASSERT (ret);
  }
}

/* From num_local_trees_per_eclass compute num_trees_per_eclass.
 * collective function */
void
t8_cmesh_gather_trees_per_eclass (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_gloidx_t temp_trees_per_eclass[T8_ECLASS_COUNT];
  int ieclass;

  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  if (cmesh->set_partition) {
    /* Copy the local values */
    /* We need to do it in a loop since we convert from locidx to gloidx.
     * memcpy is thus not possible */
    for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
      temp_trees_per_eclass[ieclass] = cmesh->num_local_trees_per_eclass[ieclass];
    }

    if (cmesh->first_tree_shared) {
      t8_eclass_t eclass;
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
    sc_MPI_Allreduce (temp_trees_per_eclass, cmesh->num_trees_per_eclass, T8_ECLASS_COUNT, T8_MPI_GLOIDX, sc_MPI_SUM,
                      comm);
  }
  else {
    /* The cmesh is not partitioned, we can just copy local_num_trees_per_eclass */
    for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
      cmesh->num_trees_per_eclass[ieclass] = cmesh->num_local_trees_per_eclass[ieclass];
    }
  }
#if T8_ENABLE_DEBUG
  /* Count the number of trees and check if it matches cmesh->num_trees */
  {
    t8_gloidx_t num_trees = 0;
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
t8_cmesh_gather_treecount_ext (const t8_cmesh_t cmesh, sc_MPI_Comm comm, const int check_commit)
{
  t8_gloidx_t tree_offset;
  int is_empty, has_empty;

  if (check_commit) {
    T8_ASSERT (t8_cmesh_is_committed (cmesh));
  }
  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  tree_offset = cmesh->first_tree_shared ? -cmesh->first_tree - 1 : cmesh->first_tree;
  if (cmesh->tree_offsets == NULL) {
    SC_CHECK_ABORT (t8_shmem_init (comm) > 0, "Error in shared memory setup.");
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    /* Only allocate the shmem array, if it is not already allocated */
    cmesh->tree_offsets = t8_cmesh_alloc_offsets (cmesh->mpisize, comm);
    t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX, cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
    /* Set the last entry to the total number of trees. */
    if (t8_shmem_array_start_writing (cmesh->tree_offsets)) {
      t8_shmem_array_set_gloidx (cmesh->tree_offsets, cmesh->mpisize, cmesh->num_trees);
    }
    t8_shmem_array_end_writing (cmesh->tree_offsets);

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
      int next_nonempty;

      const t8_gloidx_t *tree_offset_array = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
      /* there exist empty ranks, we have to recalculate the offset.
       * Each empty rank stores the offset of the next nonempty rank */
      if (is_empty) {
        next_nonempty = t8_offset_next_nonempty_rank (cmesh->mpirank, cmesh->mpisize, tree_offset_array);
        /* Set the tree offset to the first nonshared tree of the next rank */
        tree_offset = t8_offset_first (next_nonempty, tree_offset_array);
        if (tree_offset_array[next_nonempty] < 0) {
          tree_offset++;
        }
      }
      /* Communicate the new tree offsets */
      t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX, cmesh->tree_offsets, 1, T8_MPI_GLOIDX);
    }
  }
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process */
void
t8_cmesh_gather_treecount (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_gather_treecount_ext (cmesh, comm, 1);
}

/* Given a cmesh create its tree_offsets from the local number of
 * trees on each process */
void
t8_cmesh_gather_treecount_nocommit (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_gather_treecount_ext (cmesh, comm, 0);
}

/* A fast way to compute the sendrange */
static t8_locidx_t
t8_cmesh_partition_sendrange (const t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, int *send_first, int *send_last)
{
  t8_gloidx_t first_tree = t8_cmesh_get_first_treeid (cmesh_from);
  int sendfirst;
  int sendlast = -2;
  int some_owner = -1; /* Passes as argument to first/last owner functions */

  const t8_gloidx_t *offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (!cmesh_from->set_partition) {
    /* If cmesh_from is not partitioned we can send/recv all information
     * to/from ourself */
    *send_first = *send_last = cmesh_from->mpirank;
    offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
    return t8_offset_last (cmesh_from->mpirank, offset_to);
  }

  const t8_gloidx_t *offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);

  /* TODO: try to work around calling this function, since it has log(P) runtime */
  if (t8_offset_nosend (cmesh->mpirank, cmesh->mpisize, offset_from, offset_to)) {
    /* We do not send at all */
    *send_last = -2;
    *send_first = -1;
    return -1;
  }
  int flag = 0;
  if (cmesh_from->first_tree_shared == 1) {
    /* If the first tree was shared then we only send it, if
     * we own it in the new partition. In this case we send to ourself first. */
    if (t8_offset_in_range (first_tree, cmesh->mpirank, offset_to)) {
      sendfirst = cmesh->mpirank;
      flag = 1; /* We found the process */
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
    sendfirst = t8_offset_first_owner_of_tree (cmesh->mpisize, first_tree, offset_to, &some_owner);
    while (sendfirst >= 0 && flag == 0) {
      if (sendfirst == cmesh->mpirank || t8_offset_empty (sendfirst, offset_from)
          || t8_offset_first (sendfirst, offset_from) != first_tree) {
        /* We found the process if it is either ourself or it did not own first_tree
         * before */
        flag = 1;
      }
      else {
        /* Compute the next bigger process that has first_tree as local tree in
         * the new partition. */
        sendfirst = t8_offset_next_owner_of_tree (cmesh->mpisize, first_tree, offset_to, sendfirst);
      }
    }
    T8_ASSERT (flag == 1); /* We must have found the process by now */
  }
  /* Get the last local tree on cmesh_from */
  t8_gloidx_t last_tree = t8_cmesh_get_first_treeid (cmesh_from) + t8_cmesh_get_num_local_trees (cmesh_from) - 1;
  flag = 0;
  int count = 0;
  if (last_tree == t8_cmesh_get_first_treeid (cmesh_from) && cmesh_from->first_tree_shared) {
    /* If our last tree is the first tree and it is shared then we only send
     * to ourselves */
    sendlast = sendfirst;
    flag = 1;
    T8_ASSERT (sendfirst == cmesh->mpirank);
  }
  else {
    /* Check the last tree and maybe the second last tree.
     * If we do not send our last tree, we have to check the second last one.
     * This while loop is executed once for the last tree and, if unsuccessful
     * once for the second last tree. */
    while (last_tree >= first_tree && count < 2 && flag == 0) {
      count++;
      flag = 0;
      some_owner = -1; /* Reset some_owner, since we do not know an owner of last_tree */
                       /* Parse the new owners from the top and stop at the first process
       * that did not own last_tree */

      sendlast = t8_offset_last_owner_of_tree (cmesh->mpisize, last_tree, offset_to, &some_owner);
      while (sendlast >= 0 && flag == 0) {
        if (sendlast == cmesh->mpirank || t8_offset_empty (sendlast, offset_from)
            || t8_offset_first (sendlast, offset_from) != last_tree) {
          /* sendlast is either this process or it did not have
           * last_tree before */
          flag = 1;
        }
        else {
          /* Compute next smaller process that has last_tree as local tree */
          sendlast = t8_offset_prev_owner_of_tree (cmesh->mpisize, last_tree, offset_to, sendlast);
        }
      }
      /* If we did not found the alternative send here, then all procs
       * owned the last tree before and we have to check the second last tree.
       * Here a success is guaranteed. */
      last_tree--;
      /* If this second last tree is the first tree and the first tree is shared,
       * then we do not send to any other processes than ourselves. */
      if (flag == 0 && last_tree == t8_cmesh_get_first_treeid (cmesh_from) && cmesh_from->first_tree_shared) {
        sendlast = sendfirst;
        T8_ASSERT (sendfirst == cmesh->mpirank);
        flag = 1;
      }
    }
  }
  if (flag == 0) {
    /* This should never happen */
    T8_ASSERT (0 == 1);
  }
  else {
  }
  *send_first = sendfirst;
  *send_last = sendlast;

  /* Calculate the last local tree that we need to send to send_first */
  /* Set it to the last tree on send_first */
  t8_gloidx_t ret = t8_offset_last (*send_first, offset_to) - cmesh_from->first_tree;
  /* If there are actually more trees on send_first than we have, we need to send
   * all our local trees to send_first */
  ret = SC_MIN (ret, cmesh_from->num_local_trees - 1);
  if (cmesh_from->mpirank != *send_first
      && t8_offset_in_range (t8_offset_last (cmesh_from->mpirank, offset_from), *send_first, offset_from)
      && ret == cmesh_from->num_local_trees - 1) {
    /* Subtract one if our last tree already belonged to send_first,
     * and we counted this last tree. */
    ret--;
  }

  t8_debugf ("%s_first = %i, %s_last = %i, last_tree = %" T8_GLOIDX_FORMAT "\n", "send", *send_first, "send",
             *send_last, ret);

  T8_ASSERT (*send_first >= 0);
  //TODO:reactivate  T8_ASSERT (*send_last >= 0);
  T8_ASSERT ((ret >= 0 && ret < cmesh_from->num_local_trees));
  T8_ASSERT (ret == (t8_locidx_t) ret);
  return (t8_locidx_t) ret;
}

/* A fast way to compute the receive range */
static void
t8_cmesh_partition_recvrange (const t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, int *recv_first, int *recv_last)
{
  int recvfirst;
  int recvlast;
  int some_owner = -1; /* Passes as argument to first/last owner functions */

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
  const t8_gloidx_t *offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  const t8_gloidx_t *offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  /* Get the new first local tree */
  const t8_gloidx_t first_tree = t8_cmesh_get_first_treeid (cmesh);
  if (t8_offset_in_range (first_tree, cmesh->mpirank, offset_from)) {
    /* It it already was a local tree then we received it from ourselves
     * and are thus the first process we receive from */
    recvfirst = cmesh->mpirank;
  }
  else {
    /* Otherwise the first process we receive from is the smallest process that
     * had our new first tree as a local tree. */
    some_owner = -1;
    recvfirst = t8_offset_first_owner_of_tree (cmesh->mpisize, first_tree, offset_from, &some_owner);
  }
  /* Get the new last local tree */
  const t8_gloidx_t last_tree = t8_offset_last (cmesh->mpirank, offset_to);
  if (t8_offset_in_range (last_tree, cmesh->mpirank, offset_from)) {
    /* We had our last local tree as a local tree before and thus
     * we are the last process that we receive from */
    recvlast = cmesh->mpirank;
  }
  else {
    /* We receive from the smallest process that had our new last local tree
     * as a local tree. */
    if (first_tree != last_tree) {
      /* We can reuse the computed owner if first_tree == last_tree,
       * otherwise we have to reset it */
      some_owner = -1;
    }
    recvlast = t8_offset_first_owner_of_tree (cmesh->mpisize, last_tree, offset_from, &some_owner);
  }
  *recv_first = recvfirst;
  *recv_last = recvlast;
}

/* Compute the number of bytes that need to be allocated in the send buffer
 * for the neighbor entries of ghost */
static size_t
t8_partition_compute_gnb (const t8_cmesh_t cmesh_from, sc_array_t *send_as_ghost)
{
  size_t ghost_neighbor_bytes = 0, ighost;
  t8_locidx_t ghost_id;
  t8_eclass_t eclass;

  for (ighost = 0; ighost < send_as_ghost->elem_count; ighost++) {
    ghost_id = *(t8_locidx_t *) sc_array_index (send_as_ghost, ighost);
    T8_ASSERT (ghost_id >= 0);
    if (ghost_id < cmesh_from->num_local_trees) {
      /* The ghost to send is currently a tree */
      eclass = t8_cmesh_get_tree_class (cmesh_from, ghost_id);
    }
    else {
      /* The ghost to send is currently a ghost */
      T8_ASSERT (ghost_id < cmesh_from->num_local_trees + cmesh_from->num_ghosts);
      eclass = t8_cmesh_get_ghost_class (cmesh_from, ghost_id - cmesh_from->num_local_trees);
    }
    ghost_neighbor_bytes
      += t8_eclass_num_faces[eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))                     /* offset */
         + T8_ADD_PADDING (t8_eclass_num_faces[eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t))); /* padding */
  }
  return ghost_neighbor_bytes;
}

/* Compute the number of bytes that need to be allocated in the send buffer
 * for the attribute entries of all ghosts. */
static size_t
t8_partition_compute_gab (const t8_cmesh_t cmesh_from, sc_array_t *send_as_ghost, size_t *attr_info_bytes)
{
  size_t ghost_attribute_bytes = 0, ighost;
  t8_locidx_t ghost_id, ghost_id_min_offset;
  t8_cghost_t ghost;
  t8_ctree_t tree;

  T8_ASSERT (attr_info_bytes != NULL);
  *attr_info_bytes = 0;
  for (ighost = 0; ighost < send_as_ghost->elem_count; ighost++) {
    ghost_id = *(t8_locidx_t *) sc_array_index (send_as_ghost, ighost);
    T8_ASSERT (ghost_id >= 0);
    if (ghost_id < cmesh_from->num_local_trees) {
      /* This ghost is currently a local tree */
      tree = t8_cmesh_get_tree (cmesh_from, ghost_id);
      ghost_attribute_bytes += t8_cmesh_trees_attribute_size (tree);
      *attr_info_bytes += tree->num_attributes * sizeof (t8_attribute_info_struct_t);
    }
    else {
      /* This ghost is currently a ghost */
      ghost_id_min_offset = ghost_id - t8_cmesh_get_num_local_trees (cmesh_from);
      T8_ASSERT (0 <= ghost_id_min_offset && ghost_id_min_offset < t8_cmesh_get_num_ghosts (cmesh_from));
      ghost = t8_cmesh_trees_get_ghost (cmesh_from->trees, ghost_id_min_offset);
      ghost_attribute_bytes += t8_cmesh_trees_ghost_attribute_size (ghost);
      *attr_info_bytes += ghost->num_attributes * sizeof (t8_attribute_info_struct_t);
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
t8_cmesh_send_ghost (const t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from, const int p, const t8_locidx_t tree)
{
  t8_gloidx_t tree_id, *ghost_neighbors, neighbor;
  const t8_gloidx_t *from_offsets;
  t8_locidx_t *tree_neighbors;
  t8_cghost_t ghost = NULL;
  t8_ctree_t ctree = NULL;
  t8_eclass_t eclass;
  int proc, iface;
  size_t left = 0, right = cmesh->mpirank;
  int found = 0;
  const t8_gloidx_t *offset_to;

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
    ctree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, tree, &tree_neighbors, NULL);
    tree_id = tree + cmesh_from->first_tree;
    eclass = ctree->eclass;
  }
  else {
    /* Given local id belongs to a ghost. We store the ghost and its
     * global id */
    ghost
      = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees, tree - cmesh_from->num_local_trees, &ghost_neighbors, NULL);
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
    neighbor = ctree != NULL ? t8_cmesh_get_global_id ((t8_cmesh_t) cmesh_from, tree_neighbors[iface])
                             : ghost_neighbors[iface];
    if (neighbor == tree_id) {
      /* There is no neighbor at this face */
      continue;
    }
    if (!t8_offset_in_range (neighbor, p, offset_to)) {
      /* This neighbor will not be send to p */
      continue;
    }
    /* If the receiving rank is the sending rank, we definitely send the ghost */
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
      if (t8_offset_first (proc, from_offsets) + (from_offsets[proc] < 0) > neighbor) {
        /* look left */
        right = proc;
        proc = left + (proc - left) / 2;
      }
      else if (t8_offset_last (proc, from_offsets) < neighbor) {
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
      while (neighbor == t8_offset_first (proc, from_offsets) && from_offsets[proc] < 0) {
        proc--;
        while (t8_offset_empty (proc, from_offsets)) {
          /* Skip empty processes */
          proc--;
        }
      }
    }
    T8_ASSERT (0 <= proc && proc < cmesh->mpisize);
    T8_ASSERT (t8_offset_in_range (neighbor, proc, from_offsets));
    if (proc < cmesh->mpirank && t8_offset_sendstree (proc, p, neighbor, from_offsets, offset_to)) {
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
t8_cmesh_partition_copy_data (char *send_buffer, t8_cmesh_t cmesh, const t8_cmesh *cmesh_from,
                              const t8_locidx_t num_trees, const size_t attr_info_bytes,
                              const size_t ghost_attr_info_bytes, const size_t ghost_neighbor_bytes,
                              const size_t tree_neighbor_bytes, const size_t tree_attribute_bytes,
                              sc_array_t *send_as_ghost, const t8_locidx_t send_first, const t8_locidx_t send_last,
                              const size_t total_alloc, const int to_proc)
{
  t8_ctree_t tree, tree_copy;
  int num_attributes;
  size_t temp_offset_tree, temp_offset_att, temp_offset, temp_offset_data, last_offset, last_num_att, last_size,
    temp_offset_ghost_att, temp_offset_ghost_data, temp_offset_ghost, ghost_attr_info_bytes_sofar;
  size_t ghost_att_size;
  //ssize_t             last_attribute_diff;
  t8_attribute_info_struct_t *attr_info;
  void *first_attribute;
  t8_locidx_t num_ghost_send = send_as_ghost->elem_count;
  t8_locidx_t ghosts_left;
  t8_locidx_t *face_neighbor, ghost_id, itree;
  t8_gloidx_t *face_neighbor_g, *face_neighbor_gnew, new_neighbor;
  t8_cghost_t ghost, ghost_copy;
  int iface, iatt;
  int8_t *ttf_ghost, *ttf;

  /* Copy all trees to the send buffer */
  /* TODO: This is currently inefficient since we copy each tree for itself.
   *       Best practice is to copy chunks of trees out of the different part
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
    tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree, &face_neighbor, NULL);

    (void) memcpy (send_buffer + temp_offset_tree, tree, sizeof (t8_ctree_struct_t));
    temp_offset_tree += sizeof (t8_ctree_struct_t);
    /* Copy all face neighbor information to send_buffer */
    (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t)
                     + ghost_neighbor_bytes + temp_offset,
                   face_neighbor, t8_eclass_num_faces[tree->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t)));
    temp_offset += t8_eclass_num_faces[tree->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t))
                   + T8_ADD_PADDING (t8_eclass_num_faces[tree->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t)));
    /* TODO: ??? temp_offset += T8_ADD_PADDING (temp_offset) instead of the last 2 lines? */
    if (tree->num_attributes > 0) {
      /* Copy all attribute infos to send_buffer */
      (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t)
                       + ghost_neighbor_bytes + tree_neighbor_bytes + temp_offset_att,
                     T8_TREE_ATTR_INFO (tree, 0), tree->num_attributes * sizeof (t8_attribute_info_struct_t));
      temp_offset_att += tree->num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute data to send_buffer */
      (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t)
                       + ghost_neighbor_bytes + tree_neighbor_bytes + attr_info_bytes + temp_offset_data,
                     T8_TREE_ATTR (tree, T8_TREE_ATTR_INFO (tree, 0)), t8_cmesh_trees_attribute_size (tree));
      temp_offset_data += t8_cmesh_trees_attribute_size (tree);
    }
  }
  T8_ASSERT (tree_attribute_bytes == temp_offset_data);
  /* Set new face_neighbor offsets */
  /* TODO: indent bug? */
  /* Computes the offset of the face neighbors of the new trees */
  temp_offset
    = num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes;
  /* Compute the offset of the new attribute infos */
  temp_offset_att = temp_offset + tree_neighbor_bytes;
  temp_offset_tree = 0;

  /* Set attribute offsets of trees and attribute data offsets of info objects */
  tree_copy = NULL;
  last_num_att = 0;
  last_size = 0;
  last_offset = attr_info_bytes;

  for (itree = send_first; itree <= send_last; itree++) {
    /* Get the current tree */
    tree_copy = (t8_ctree_t) (send_buffer + temp_offset_tree);

    /* new neighbor offset of tree */
    tree_copy->neigh_offset = temp_offset - temp_offset_tree;

    /* Set new face neighbor entries, since we store local ids we have to adapt
     * to the local ids of the new process */
    face_neighbor = (t8_locidx_t *) T8_TREE_FACE (tree_copy);
    for (iface = 0; iface < t8_eclass_num_faces[tree_copy->eclass]; iface++) {
      t8_cmesh_partition_send_change_neighbor (cmesh, (t8_cmesh_t) cmesh_from, face_neighbor + iface, to_proc);
    }

    /* compute neighbor offset for next tree */
    temp_offset += t8_eclass_num_faces[tree_copy->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t))
                   + T8_ADD_PADDING (t8_eclass_num_faces[tree_copy->eclass] * (sizeof (t8_locidx_t) + sizeof (int8_t)));

    /* new attribute offset for tree */
    tree_copy->att_offset = temp_offset_att - temp_offset_tree;
    if (tree_copy->num_attributes > 0) {
      attr_info = T8_TREE_ATTR_INFO (tree_copy, 0);
      attr_info->attribute_offset = last_offset + last_size - last_num_att * sizeof (t8_attribute_info_struct_t);
      last_offset = attr_info->attribute_offset;
      last_size = attr_info->attribute_size;

      /* set new attribute data offsets */
      for (size_t iattribute = 1; iattribute < (size_t) tree_copy->num_attributes; iattribute++) {
        attr_info++;
        attr_info->attribute_offset = last_offset + last_size;
        last_offset = attr_info->attribute_offset;
        last_size = attr_info->attribute_size;
      }
      temp_offset_att += tree_copy->num_attributes * sizeof (t8_attribute_info_struct_t);
    }
    temp_offset_tree += sizeof (t8_ctree_struct_t);
    last_num_att = tree_copy->num_attributes;
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
  for (size_t isendghost = 0; isendghost < send_as_ghost->elem_count; isendghost++) {
    ghost_id = *((t8_locidx_t *) sc_array_index (send_as_ghost, isendghost));
    ghost_copy
      = (t8_cghost_t) (send_buffer + num_trees * sizeof (t8_ctree_struct_t) + isendghost * sizeof (t8_cghost_struct_t));
    /* Get the correct element class from the ghost_id.
     * For this we have to check whether ghost_id points to a local tree or ghost */
    ghost_copy->eclass = ghost_id < cmesh_from->num_local_trees
                           ? t8_cmesh_get_tree_class ((t8_cmesh_t) cmesh_from, ghost_id)
                           : t8_cmesh_get_ghost_class ((t8_cmesh_t) cmesh_from, ghost_id - cmesh_from->num_local_trees);
    ghost_copy->neigh_offset = temp_offset - temp_offset_ghost; /* New face neighbor offset */
    face_neighbor_gnew = (t8_gloidx_t *) T8_GHOST_FACE (ghost_copy);
    ttf_ghost = T8_GHOST_TTF (ghost_copy);
    if (ghost_id >= cmesh_from->num_local_trees) {
      /* The ghost that we send was a local ghost */
      ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees, ghost_id - cmesh_from->num_local_trees, &face_neighbor_g,
                                            &ttf);
      tree = NULL;
      /* Set entries of the new ghost */
      ghost_copy->eclass = ghost->eclass;
      ghost_copy->treeid = ghost->treeid;
      ghost_copy->num_attributes = ghost->num_attributes;
      num_attributes = ghost_copy->num_attributes;
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
              t8_eclass_num_faces[ghost_copy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t)));
    }
    else {
      /* The ghost we send was a local tree */
      T8_ASSERT (0 <= ghost_id && ghost_id < cmesh_from->num_local_trees);
      tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, ghost_id, &face_neighbor, &ttf);
      ghost = NULL;
      T8_ASSERT (ghost_id == tree->treeid);
      /* Set entries of the new ghost */
      ghost_copy->eclass = tree->eclass;
      ghost_copy->treeid = ghost_id + cmesh_from->first_tree;
      ghost_copy->num_attributes = tree->num_attributes;
      num_attributes = ghost_copy->num_attributes;
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
      for (iface = 0; iface < t8_eclass_num_faces[ghost_copy->eclass]; iface++) {
        if (face_neighbor[iface] < 0) {
          /* TODO: think about this */
          new_neighbor = -1; /* boundary indicator */
        }
        else {
          /* Compute global index from local index */
          new_neighbor = t8_cmesh_get_global_id ((t8_cmesh_t) cmesh_from, face_neighbor[iface]);
        }
        face_neighbor_gnew[iface] = new_neighbor;
      }
      /* Copy tree_to_face entries */
      memcpy (ttf_ghost, ttf, t8_eclass_num_faces[ghost_copy->eclass] * sizeof (int8_t));
    } /* Done distinction between from tree and from ghost */

    /* Compute and store new attribute offset of this ghost */
    ghosts_left = send_as_ghost->elem_count - isendghost;
    ghost_copy->att_offset = ghosts_left * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes + tree_neighbor_bytes
                             + attr_info_bytes + tree_attribute_bytes + temp_offset_ghost_att;
    if (num_attributes > 0) {
      size_t this_ghosts_att_info_size;
      t8_attribute_info_struct_t *first_attr_info;

      /* The byte count of this ghosts attribute info structs */
      this_ghosts_att_info_size = num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute info data of this ghost */
      first_attr_info = (t8_attribute_info_struct_t *) T8_GHOST_FIRST_ATT_INFO (ghost_copy);
      memcpy (first_attr_info, attr_info, this_ghosts_att_info_size);
      temp_offset_ghost_att += this_ghosts_att_info_size;

      /* Compute all new attribute data offsets */
      for (iatt = 0; iatt < num_attributes; iatt++) {
        /* Get the current attribute info */
        attr_info = T8_GHOST_ATTR_INFO (ghost_copy, iatt);
        /* The new attribute offset is the offset from the first att_info to the data.
         * Thus, the count of the bytes occupied by the att_info (ghosts_attr_info_bytes)
         * plus the count of all attributes before this attribute (this_data_temp_offset).*/
        /* all att info from this ghost and after. + all attributes before this attribute  */
        attr_info->attribute_offset = ghost_attr_info_bytes - ghost_attr_info_bytes_sofar + temp_offset_ghost_data;
        temp_offset_ghost_data += attr_info->attribute_size;
      }
      ghost_attr_info_bytes_sofar += num_attributes * sizeof (t8_attribute_info_struct_t);
      /* Copy all attribute data of this ghost */
      memcpy (T8_GHOST_ATTR (ghost_copy, first_attr_info), first_attribute, ghost_att_size);
    } /* end num_attributes > 0 */
    /* compute new offsets */
    temp_offset += t8_eclass_num_faces[ghost_copy->eclass] * (sizeof (t8_gloidx_t) + sizeof (int8_t)) /* offset */
                   + T8_ADD_PADDING (t8_eclass_num_faces[ghost_copy->eclass]
                                     * (sizeof (t8_gloidx_t) + sizeof (int8_t))); /* padding */
    temp_offset_ghost += sizeof (t8_cghost_struct_t);
  }
  /* Store number of trees and ghosts at the end of send buffer */
  *(t8_locidx_t *) (send_buffer + total_alloc - 2 * sizeof (t8_locidx_t)) = num_trees;
  *(t8_locidx_t *) (send_buffer + total_alloc - sizeof (t8_locidx_t)) = num_ghost_send;
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
t8_cmesh_partition_sendtreeloop (t8_cmesh_t cmesh, const t8_cmesh *cmesh_from, const t8_locidx_t range_start,
                                 const t8_locidx_t range_end, size_t *tree_neighbor_bytes, size_t *attr_bytes,
                                 size_t *attr_info_bytes, int8_t *ghost_flag_send, const int iproc,
                                 sc_array_t *send_as_ghost)
{
  t8_ctree_t tree;
  t8_locidx_t neighbor, *face_neighbor, itree;
  int8_t *ttf;
  int iface;
#if T8_ENABLE_DEBUG
  const t8_gloidx_t *offset_from, *offset_to;

  if (cmesh_from->set_partition) {
    offset_from = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
  }
  else {
    offset_from = NULL;
  }
  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
#endif
  /* loop over all trees that will be send */
  for (itree = range_start; itree <= range_end; itree++) {
    /* test if we really send the tree itree to the process iproc */
    T8_ASSERT (
      (!cmesh_from->set_partition && iproc == cmesh_from->mpirank)
      || t8_offset_sendstree (cmesh_from->mpirank, iproc, itree + cmesh_from->first_tree, offset_from, offset_to));
    tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree, &face_neighbor, &ttf);
    /* Count the additional memory needed per tree from neighbors */
    *tree_neighbor_bytes += t8_eclass_num_faces[tree->eclass] * (sizeof (*face_neighbor) + sizeof (*ttf));
    /* TODO: change padding to sizeof (void*) */
    *tree_neighbor_bytes += T8_ADD_PADDING (*tree_neighbor_bytes); /* padding to make number of bytes per tree
                                                                           a multiple of 4 */
    /*  Compute number of attribute bytes in this tree range.
     *       Not every tree has an attribute */
    *attr_info_bytes += tree->num_attributes * sizeof (t8_attribute_info_struct_t);
    *attr_bytes += t8_cmesh_trees_attribute_size (tree);

    /* loop over all faces of each tree to determine ghost to send */
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      neighbor = face_neighbor[iface];
      if (neighbor >= 0 && neighbor != itree) { /* Consider only non-boundary neighbors */
        if ((neighbor < cmesh_from->num_local_trees && (neighbor < range_start || neighbor > range_end))
            || neighbor >= cmesh_from->num_local_trees) {
          /* neighbor is a local tree or local ghost and will be ghost on iproc */
          if (ghost_flag_send[neighbor] == 0 && t8_cmesh_send_ghost (cmesh, cmesh_from, iproc, neighbor)) {
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
t8_cmesh_partition_sendloop (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from, int *num_request_alloc, int *send_first,
                             int *send_last, char ***send_buffer, char **my_buffer, size_t *my_buffer_bytes,
                             sc_MPI_Request **requests, sc_MPI_Comm comm)
{
  size_t attr_bytes = 0, tree_neighbor_bytes, ghost_neighbor_bytes, attr_info_bytes, ghost_attribute_bytes,
         ghost_attr_info_bytes;
  size_t total_alloc;
  int iproc, flag;
  int mpiret, num_send_mpi = 0;
  char *buffer;
  t8_locidx_t num_trees, num_ghost_send, range_start, range_end;
  sc_array_t send_as_ghost; /* Stores local id's of trees and ghosts that will be send as ghosts */
  int8_t *ghost_flag_send;  /* For each local tree and ghost set to 1 if it is in send_as_ghost */
  const t8_gloidx_t *offset_from, *offset_to;

  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  /* We use two different flag arrays here, since ghost_flag_send needs to
   * be reset for each process we send to, while ghost_flag_keep keeps is entries.
   * Otherwise we could have set bitflags and only use one array */
  ghost_flag_send = T8_ALLOC (int8_t, cmesh_from->num_local_trees + cmesh_from->num_ghosts);

  sc_array_init (&send_as_ghost, sizeof (t8_locidx_t));

  range_end = t8_cmesh_partition_sendrange (cmesh, (t8_cmesh_t) cmesh_from, send_first, send_last);

  offset_to = t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
  if (cmesh_from->set_partition) {
    range_start = cmesh_from->first_tree_shared; /* Stores the first local tree that was not send yet */
    /* We do not send the first tree if it is shared, so we start with the second tree then */
    if (t8_offset_in_range (cmesh_from->first_tree, cmesh->mpirank, offset_to)) {
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
    while (cmesh_from->set_partition
           && !t8_offset_sendstree (cmesh_from->mpirank, iproc, range_start + cmesh_from->first_tree, offset_from,
                                    offset_to)) {
      /* We skip trees that we do not send */
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
    t8_cmesh_partition_sendtreeloop (cmesh, cmesh_from, range_start, range_end, &tree_neighbor_bytes, &attr_bytes,
                                     &attr_info_bytes, ghost_flag_send, iproc, &send_as_ghost);
    /* loop over trees ends here */
    /* Calculate and allocate memory for send buffer */
    /**********************************************************/
    /*
     *      The data that we send has the layout
     *
     * | tree 0,1...| ghost 0,1...| face_neighbors0, tree_to_face0,...| ghost_neighbors 0,1,...| Attinfo00,01,...,Attinfoend | attdata00,01,... |
     *
     */
    num_trees = SC_MAX (range_end - range_start + 1, 0);

    num_ghost_send = send_as_ghost.elem_count;
    /* parse through send_as_ghost to compute ghost_neighbor_bytes */
    ghost_neighbor_bytes = t8_partition_compute_gnb (cmesh_from, &send_as_ghost);
    /* parse through send_as_ghost to compute ghost_attribute_bytes and attr_info_bytes */
    ghost_attribute_bytes = t8_partition_compute_gab (cmesh_from, &send_as_ghost, &ghost_attr_info_bytes);
    /* Total number of bytes that we send to the other process */
    total_alloc = num_trees * sizeof (t8_ctree_struct_t) + num_ghost_send * sizeof (t8_cghost_struct_t)
                  + ghost_neighbor_bytes + tree_neighbor_bytes + attr_info_bytes + attr_bytes + ghost_attribute_bytes
                  + ghost_attr_info_bytes;
    /* Extra space to store the number of trees and ghosts in the send buffer */
    total_alloc += 2 * sizeof (t8_locidx_t);

    /* Debugging output about shipped trees and ghosts. */
    t8_debugf ("NT %i -- each %zdB\n", num_trees, sizeof (t8_ctree_struct_t));
    t8_debugf ("NGS %i - each %zdB\n", num_ghost_send, sizeof (t8_cghost_struct_t));
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
      buffer = (*send_buffer)[iproc - *send_first - flag] = T8_ALLOC (char, total_alloc);
    }
    else if (num_trees > 0 || num_ghost_send > 0) {
      *my_buffer = buffer = T8_ALLOC (char, total_alloc);
      *my_buffer_bytes = total_alloc;
    }
    else {
      my_buffer = NULL;
      buffer = NULL;
    }
    T8_ASSERT (num_trees + num_ghost_send == 0 || (!cmesh_from->set_partition && iproc == cmesh_from->mpirank)
               || t8_offset_sendsto (cmesh->mpirank, iproc, offset_from, offset_to));

    /* Copy all data to the send buffer */
    t8_cmesh_partition_copy_data (buffer, cmesh, cmesh_from, num_trees, attr_info_bytes, ghost_attr_info_bytes,
                                  ghost_neighbor_bytes, tree_neighbor_bytes, attr_bytes, &send_as_ghost, range_start,
                                  range_end, total_alloc, iproc);

    /* If we send to a remote process we post the MPI_Isend here */
    if (iproc != cmesh->mpirank) {
      if (num_trees + num_ghost_send > 0) {
        /* send buffer to remote process */
        t8_debugf ("Post send of %i trees/%zd bytes to %i\n",
                   *(t8_locidx_t *) (buffer + total_alloc - 2 * sizeof (t8_locidx_t)), total_alloc, iproc - flag);
        mpiret = sc_MPI_Isend (buffer, total_alloc, sc_MPI_BYTE, iproc, T8_MPI_PARTITION_CMESH, comm,
                               *requests + iproc - flag - *send_first);
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
    while (iproc < *send_last && !t8_offset_sendsto (cmesh_from->mpirank, iproc, offset_from, offset_to)) {
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
      t8_gloidx_t first_send;

      first_send = t8_offset_first (iproc, offset_to);
      while (first_send < t8_cmesh_get_num_local_trees (cmesh_from)
             && !t8_offset_sendstree (cmesh_from->mpirank, iproc, first_send, offset_from, offset_to)) {

        first_send++;
      }
      range_start = first_send - cmesh_from->first_tree;
      t8_debugf ("RS: %i\n", range_start);

      /* add number of trees in iproc to send to range_end */
      /* We have to be careful with locidx overflow when we go out of bounds
       * of our process */

      /* Set range end to the last tree of the receiver or to our last tree,
       * whichever is smaller */
      range_end = SC_MIN (t8_offset_last (iproc, offset_to) - cmesh_from->first_tree, cmesh_from->num_local_trees - 1);
      t8_debugf ("RE: %i\n", range_end);
      /* If the last tree was already present on the receiving process
       * we exclude it from sending */
      if (cmesh_from->mpirank != iproc && range_end == cmesh_from->num_local_trees - 1
          && !t8_offset_empty (iproc, offset_from)
          && t8_offset_first (iproc, offset_from) == cmesh_from->first_tree + cmesh_from->num_local_trees - 1) {
        range_end--;
        t8_debugf ("RE: %i\n", range_end);
      }
      if (range_end < range_start || range_start >= t8_cmesh_get_num_local_trees (cmesh_from)) {
        /* We do not send to this process and are finished */
        iproc = *send_last + 1;
        *send_last = -2;
      }
    }
  } /* sending loop ends here */
  T8_FREE (ghost_flag_send);
  sc_array_reset (&send_as_ghost);
  t8_debugf ("End send loop\n");
  return num_send_mpi;
}

/**
 * Receive a partition message from another process.
 * \param[in] cmesh The cmesh to receive the partition message for.
 * \param[in] comm The MPI communicator to use for the communication.
 * \param[in] proc_recv The rank of the process to receive the message from.
 * \param[in, out] status The MPI status object containing the message information.
 * \param[in] local_procid The local process id's of the processes that we receive from.
 * \param[in] recv_first The first process rank that we receive from.
 * \param[in, out] num_ghosts The number of ghosts that we received in this message.
 */
void
t8_cmesh_partition_receive_message (t8_cmesh_t cmesh, sc_MPI_Comm comm, const int proc_recv, sc_MPI_Status *status,
                                    const int *local_procid, const int recv_first, t8_locidx_t *num_ghosts)
{
  int mpiret;
  int recv_bytes;
  t8_part_tree_t recv_part;

  T8_ASSERT (proc_recv == status->MPI_SOURCE);
  T8_ASSERT (status->MPI_TAG == T8_MPI_PARTITION_CMESH);

  mpiret = sc_MPI_Get_count (status, sc_MPI_BYTE, &recv_bytes);
  SC_CHECK_MPI (mpiret);
  /* Allocate receive buffer */
  recv_part = t8_cmesh_trees_get_part (cmesh->trees, local_procid[proc_recv - recv_first]);
  /* take first tree of part and allocate recv_bytes */
  recv_part->first_tree = T8_ALLOC (char, recv_bytes);
  /* Receive message */
  mpiret = sc_MPI_Recv (recv_part->first_tree, recv_bytes, sc_MPI_BYTE, proc_recv, T8_MPI_PARTITION_CMESH, comm,
                        sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);
  /* Read num trees and num ghosts */
  recv_part->num_trees = *((t8_locidx_t *) (recv_part->first_tree + recv_bytes - 2 * sizeof (t8_locidx_t)));
  recv_part->num_ghosts = *((t8_locidx_t *) (recv_part->first_tree + recv_bytes - sizeof (t8_locidx_t)));
  *num_ghosts += recv_part->num_ghosts;

  t8_debugf ("Received %i trees/%i ghosts/%i bytes from %i to %i\n", recv_part->num_trees, recv_part->num_ghosts,
             recv_bytes, proc_recv, local_procid[proc_recv - recv_first]);
  /* If we are profiling, we count the number of trees and ghosts that
   * we received. */
  if (cmesh->profile != NULL && proc_recv != cmesh->mpirank) {
    cmesh->profile->partition_ghosts_recv += recv_part->num_ghosts;
    cmesh->profile->partition_trees_recv += recv_part->num_trees;
  }
}

/* Loop over all the processes that we receive trees from and get the
 * MPI messages. */
/* stores the number of received ghosts in num_ghosts */
/* fr and lr are only for debugging, see t8_cmesh_partition_debug_listprocs */
/* TODO: Remove the const qualifier at the cmesh_from parameter */
static void
t8_cmesh_partition_recvloop (t8_cmesh_t cmesh, const t8_cmesh *cmesh_from, const t8_gloidx_t *tree_offset,
                             char *my_buffer, size_t my_buffer_bytes, sc_MPI_Comm comm, [[maybe_unused]] int fr,
                             [[maybe_unused]] int lr)
{
  int num_receive, *local_procid; /* ranks of the processor from which we will receive */
  int mpiret, proc_recv, iproc;
  int recv_first, recv_last; /* ranks of the processor from which we will receive */
  int num_parts, myrank_part;
  t8_locidx_t num_trees, num_ghosts;
  const t8_gloidx_t *from_offsets;
  t8_part_tree_t recv_part;
  sc_MPI_Status status;

  num_trees = t8_offset_num_trees (cmesh->mpirank, tree_offset);
  /* Receive from other processes */
  if (cmesh_from->set_partition) {
    T8_ASSERT (cmesh_from->tree_offsets != NULL);
    from_offsets = t8_shmem_array_get_gloidx_array (cmesh_from->tree_offsets);
    t8_cmesh_partition_recvrange (cmesh, (t8_cmesh_t) cmesh_from, &recv_first, &recv_last);
#if T8_ENABLE_DEBUG
    if (recv_first <= recv_last) {
      T8_ASSERT (fr == recv_first);
      T8_ASSERT (lr == recv_last);
    }
#endif

    num_parts = t8_offset_range_send (recv_first, recv_last, cmesh->mpirank, from_offsets, tree_offset);
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
    if (t8_offset_sendsto (cmesh->mpirank, cmesh->mpirank, from_offsets, tree_offset)) {
      num_receive = num_parts - 1;
    }
    else {
      num_receive = num_parts;
    }
  }
  else {
    num_receive = 0;
  }

  num_receive = SC_MAX (num_receive, 0); /* num_receive could get negative, so we set it
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
    local_procid[iproc] = !t8_offset_sendsto (iproc + recv_first, cmesh->mpirank, from_offsets, tree_offset)
                            ? local_procid[iproc - 1]
                            : local_procid[iproc - 1] + 1;
  }
  if (cmesh_from->set_partition) {
    if (recv_first <= cmesh->mpirank && cmesh->mpirank <= recv_last
        && t8_offset_sendsto (cmesh->mpirank, cmesh->mpirank, from_offsets, tree_offset)) {
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

  T8_ASSERT (recv_first == -1 || !cmesh_from->set_partition || !t8_offset_empty (recv_first, from_offsets));

  /**************************************************/
  /*            Receive MPI messages                */
  /**************************************************/

  /****     Setup     ****/

  if (num_receive > 0) {

    /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for an message for each
     * sender and if there is one we receive it and remove the sender from
     * the list.
     * The last message can be received via probe */
    while (num_receive > 0) {
      t8_debugf ("Probing for %i messages.\n", num_receive);
      mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_PARTITION_CMESH, comm, &status);
      SC_CHECK_MPI (mpiret);
      num_receive--;
      /* There is a message to receive */
      proc_recv = status.MPI_SOURCE;
      /* TODO: assert that proc_recv is still contained in the list of receivers. */
      T8_ASSERT (status.MPI_TAG == T8_MPI_PARTITION_CMESH);
      T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last
                 && t8_offset_sendsto (proc_recv, cmesh->mpirank, from_offsets, tree_offset));
      t8_cmesh_partition_receive_message (cmesh, comm, proc_recv, &status, local_procid, recv_first, &num_ghosts);
    }
  }

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
    recv_part->num_trees = *(t8_locidx_t *) (recv_part->first_tree + my_buffer_bytes - 2 * sizeof (t8_locidx_t));
    recv_part->num_ghosts = *(t8_locidx_t *) (recv_part->first_tree + my_buffer_bytes - sizeof (t8_locidx_t));
    num_ghosts += recv_part->num_ghosts;

    t8_debugf ("Received %i trees/%i ghosts from myself to %i\n", recv_part->num_trees, recv_part->num_ghosts,
               myrank_part);
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
t8_cmesh_partition_debug_listprocs (const t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, sc_MPI_Comm comm, int *fs,
                                    int *ls, int *fr, int *lr)
{
  int mpiret, mpisize, mpirank, p;
  char out[BUFSIZ] = "";
  const t8_gloidx_t *from, *to;

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
      snprintf (out + strlen (out), BUFSIZ - strlen (out), "%i%c ", p, p == mpisize - 1 ? '!' : ',');
      *fs = SC_MIN (*fs, p);
      *ls = SC_MAX (*ls, p);
    }
  }
  t8_debugf ("I send to: %s\n", out);
  std::strcpy (out, " ");
  if (cmesh_from->set_partition) {
    for (p = 0; p < mpisize; p++) {
      if (t8_offset_sendsto (p, mpirank, from, to)) {
        snprintf (out + strlen (out), BUFSIZ - strlen (out), "%i%c ", p, p == mpisize - 1 ? '!' : ',');
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
t8_cmesh_partition_given (const t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from, const t8_gloidx_t *tree_offset,
                          sc_MPI_Comm comm)
{
  int send_first, send_last, num_request_alloc; /* ranks of the processor to which we will send */
  int iproc, num_send_mpi, mpiret;
  size_t my_buffer_bytes = -1;
  char **send_buffer = NULL, *my_buffer = NULL;

  int fs, ls;
  int fr = 0;
  int lr = 0;

  sc_MPI_Request *requests = NULL;
  t8_locidx_t num_ghosts, itree, num_trees;
  t8_part_tree_t recv_part;
  t8_ctree_t tree;

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
    t8_cmesh_partition_debug_listprocs (cmesh, (t8_cmesh_t) cmesh_from, comm, &fs, &ls, &fr, &lr);
  }

  /*********************************************/
  /*        Done with setup                    */
  /*********************************************/

  /* Send all trees and ghosts */
  num_send_mpi = t8_cmesh_partition_sendloop (cmesh, (t8_cmesh_t) cmesh_from, &num_request_alloc, &send_first,
                                              &send_last, &send_buffer, &my_buffer, &my_buffer_bytes, &requests, comm);
  T8_ASSERT (!cmesh_from->set_partition || send_first == -1 || send_first == fs);
  T8_ASSERT (!cmesh_from->set_partition || send_last == -2 || send_last == ls);

  /* receive all trees and ghosts */
  t8_cmesh_partition_recvloop (cmesh, cmesh_from, tree_offset, my_buffer, my_buffer_bytes, comm, fr, lr);
  if (num_send_mpi > 0) {
    mpiret = sc_MPI_Waitall (num_request_alloc, requests, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* Clean-up */
  for (iproc = 0; iproc < send_last - send_first + !(cmesh->mpirank >= send_first && cmesh->mpirank <= send_last);
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
    for (itree = recv_part->first_tree_id; itree < recv_part->first_tree_id + recv_part->num_trees; itree++) {
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
    /* Assign new local ids to the ghosts of this part, also set ghost_to_proc */
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
  t8_cmesh_t cmesh_from;
  t8_gloidx_t last_tree;
  const t8_gloidx_t *tree_offsets;

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
    const t8_scheme *scheme = cmesh->set_partition_scheme; /* The refinement scheme */
    /* Compute first and last tree index */
    T8_ASSERT (cmesh->tree_offsets == NULL);
    T8_ASSERT (scheme != NULL);
    t8_cmesh_uniform_bounds_for_irregular_refinement (cmesh_from, cmesh->set_partition_level, scheme,
                                                      &cmesh->first_tree, NULL, &last_tree, NULL,
                                                      &cmesh->first_tree_shared, comm);

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
    cmesh->num_local_trees = t8_offset_num_trees (cmesh->mpirank, tree_offsets);
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
  /* Deactivate the active tree. Tree related data (such as vertices) might have been moved by the new partition and
   * has to be loaded again if needed. */
  if (cmesh->geometry_handler != NULL) {
    cmesh->geometry_handler->deactivate_tree ();
  }
  /* If profiling is enabled, we measure the runtime of this routine. */
  if (cmesh->profile) {
    /* Runtime = current_time - start_time */
    cmesh->profile->partition_runtime = sc_MPI_Wtime () - cmesh->profile->partition_runtime;
  }

  t8_global_productionf ("Done cmesh partition\n");
}

void
t8_cmesh_offset_print ([[maybe_unused]] const t8_cmesh_t cmesh, [[maybe_unused]] sc_MPI_Comm comm)
{
#if T8_ENABLE_DEBUG
  int offset_isnew = 0;

  if (cmesh->set_partition) {
    if (cmesh->tree_offsets == NULL) {
      t8_cmesh_gather_treecount (cmesh, comm);
      offset_isnew = 1;
    }
    t8_offset_print (cmesh->tree_offsets, comm);
    T8_ASSERT (t8_offset_consistent (cmesh->mpisize, cmesh->tree_offsets, cmesh->num_trees));
    if (offset_isnew == 1) {
      t8_shmem_array_destroy (&cmesh->tree_offsets);
      T8_ASSERT (cmesh->tree_offsets == NULL);
    }
    else {
      t8_debugf ("Replicated cmesh with %lli trees.\n", (long long) cmesh->num_trees);
    }
  }
#endif
}

/* Create a partition that concentrates everything at a given proc */
t8_shmem_array_t
t8_cmesh_offset_concentrate (const int proc, sc_MPI_Comm comm, const t8_gloidx_t num_trees)
{
  int mpirank, mpiret, mpisize, iproc;
  t8_shmem_array_t shmem_array;
  t8_gloidx_t *offsets;
#if T8_ENABLE_DEBUG
  char out[BUFSIZ] = "";
#endif

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  shmem_array = t8_cmesh_alloc_offsets (mpisize, comm);
  if (t8_shmem_array_start_writing (shmem_array)) {
    offsets = t8_shmem_array_get_gloidx_array_for_writing (shmem_array);
    offsets[0] = 0;
    for (iproc = 1; iproc <= mpisize; iproc++) {
      if (iproc == proc + 1) {
        offsets[iproc] = num_trees;
      }
      else {
        offsets[iproc] = offsets[iproc - 1];
      }
#if T8_ENABLE_DEBUG
      snprintf (out + strlen (out), BUFSIZ - strlen (out), "%li,", offsets[iproc]);
#endif
    }
#if T8_ENABLE_DEBUG
    t8_debugf ("Partition with offsets:0,%s\n", out);
#endif
  }
  t8_shmem_array_end_writing (shmem_array);

  T8_ASSERT (t8_offset_consistent (mpisize, shmem_array, num_trees));
  return shmem_array;
}

/* Create a random partition */
/* if shared is nonzero than first trees can be shared */
t8_shmem_array_t
t8_cmesh_offset_random (sc_MPI_Comm comm, const t8_gloidx_t num_trees, const int shared, const unsigned seed)
{
  int iproc, mpisize, mpiret, random_number, mpirank;
  int first_shared;
  unsigned u_seed;
  t8_gloidx_t *offsets;
  t8_shmem_array_t shmem_array;
  t8_gloidx_t new_first;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  shmem_array = t8_cmesh_alloc_offsets (mpisize, comm);

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

  if (t8_shmem_array_start_writing (shmem_array)) {
    offsets = t8_shmem_array_get_gloidx_array_for_writing (shmem_array);

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
      if (shared && new_first < num_trees) { /* new first is num_trees, this process must be empty */
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
  }
  t8_shmem_array_end_writing (shmem_array);

  T8_ASSERT (t8_offset_consistent (mpisize, shmem_array, num_trees));
  return shmem_array;
}

t8_shmem_array_t
t8_cmesh_offset_percent (const t8_cmesh_t cmesh, sc_MPI_Comm comm, const int percent)
{
  t8_gloidx_t new_first_tree, old_first_tree;
  t8_locidx_t old_num_trees_pm1;
  t8_shmem_array_t partition_array;
  const t8_gloidx_t *old_partition;
  int mpirank, mpisize, mpiret;
  int created = 0;
#if T8_ENABLE_DEBUG
  int total = 0;
  int proc_perc = percent;
  sc_MPI_Allreduce (&proc_perc, &total, 1, sc_MPI_INT, sc_MPI_SUM, comm);
  T8_ASSERT (total == 100);
#endif

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
  old_num_trees_pm1 = t8_offset_num_trees (mpirank > 0 ? mpirank - 1 : 0, old_partition);
  /* Compute the new first local tree */
  new_first_tree = old_first_tree - old_num_trees_pm1 * percent / 100;
  /* Compute the new entry in the offset array.
   * If the old first tree was shared, then the new one will be as well
   * and if not it will not be shared. */
  if (mpirank != 0) {
    new_first_tree = t8_offset_first_tree_to_entry (new_first_tree, cmesh->first_tree_shared);
  }
  else {
    new_first_tree = 0;
  }
  t8_shmem_array_allgather (&new_first_tree, 1, T8_MPI_GLOIDX, partition_array, 1, T8_MPI_GLOIDX);
  if (t8_shmem_array_start_writing (partition_array)) {
    t8_shmem_array_set_gloidx (partition_array, mpisize, t8_cmesh_get_num_trees (cmesh));
  }
  t8_shmem_array_end_writing (partition_array);

  if (created) {
    /* We needed to create the old partition array and thus we clean it up
     * again. */
    t8_shmem_array_destroy (&cmesh->tree_offsets);
  }
  T8_ASSERT (t8_offset_consistent (mpisize, partition_array, t8_cmesh_get_num_trees (cmesh)));
  return partition_array;
}

/* Create a repartition array, where each process sends half of its
 * trees to the next process. The last process does not send any trees. */
/* TODO: This function was not tested with shared trees yet. */
t8_shmem_array_t
t8_cmesh_offset_half (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  return t8_cmesh_offset_percent (cmesh, comm, 50);
}
