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

/* Compute first and last process to which we will send data */
/* Returns the local tree_id of last local tree on send_first. */
/* TODO: This can propably be optimized since we changed to tree_offsets. using binary search or so */
static t8_locidx_t
t8_cmesh_partition_sendrange (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                              int *send_first, int *send_last)
{
  int                 iproc;
  t8_gloidx_t         first_tree = 0, last_tree, last_local_tree, ret;

  /* p_i is first process we send to if and only if
   *
   * tree_first[p_i]    <=    tree_first    <=    tree_last[pi]
   *      |                       ||                    |
   *  new_partition     cmesh_from->first_tree     new_partition
   */
  /* p_i is last process we send to if and only if
   *
   * tree_first[p_i]    <=    tree_last    <=    tree_last[pi]
   *        |                     ||                    |
   *  new_partition     cmesh_from->first_tree     new_partition
   *                       + num_local_trees-1
   */

  *send_first = -1;
  *send_last = -1;
  last_local_tree = cmesh_from->first_tree + cmesh_from->num_local_trees - 1;
  for (iproc = 0; iproc < cmesh->mpisize && last_tree > last_local_tree;
       iproc++) {
    /* first tree stores new first_tree of process iproc */
    /* last tree the new last tree of process iproc */
    first_tree = t8_glo_abs (cmesh->tree_offsets[iproc]);
    last_tree = t8_glo_abs (cmesh->tree_offsets[iproc+1]) - 1 +
        cmesh->tree_offsets[iproc + 1 ] < 0;

    if (*send_first == -1 &&
        first_tree <= cmesh_from->first_tree && cmesh_from->first_tree <=
        last_tree) {
      /* set send_first */
      *send_first = iproc;
      ret = last_tree - cmesh->first_tree;
    }
    if (*send_last == -1 &&
        first_tree <= last_local_tree && last_local_tree <= last_tree) {
      /* set send_last */
      *send_last = iproc;
    }
  }
  T8_ASSERT (*send_first >= 0);
  T8_ASSERT (*send_last >= 0);
  T8_ASSERT (ret >= 0 && ret <= cmesh->num_local_trees);
  T8_ASSERT (ret == (t8_locidx_t) ret);
  return (t8_locidx_t) ret;
}

static void
t8_cmesh_partition_recvrange (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                              int *recv_first, int *recv_last)
{
  /* p_i is first process we receive from if and only if
   *
   * tree_first[p_i]    <=    tree_first    <=    tree_last[pi]
   *      |                       ||                    |
   *  old_partition        cmesh->first_tree     old_partition
   */
  /* p_i is last process we send to if and only if
   *
   * tree_first[p_i]    <=    tree_last    <=    tree_last[pi]
   *        |                     ||                    |
   *  old_partition       cmesh->first_tree     old_partition
   *                      + num_local_trees-1
   */
  /* This is the reversed situation of send first/last, so we can
   * use this function with cmesh and cmesh_from changed. However, if we do
   * so we need the shared array on cmesh_from to exist. See also TODO comment
   * in partition_given. */

  T8_ASSERT (cmesh_from->tree_offsets != NULL);
  (void) t8_cmesh_partition_sendrange (cmesh_from, cmesh, recv_first,
                                       recv_last);
}

/* Determine whether a local tree or ghost should be send to process p.
 * This is the case if and only if:
 *  - tree is not already a ghost or tree on p
 * and
 *  - we are the smallest rank under all procs sending to p that
 *    has this tree as ghost or local tree. */
static int
t8_cmesh_send_ghost (t8_cmesh_t cmesh, int p, t8_locidx_t tree)
{
  return 0;
}

/* copy all tree/ghost/attribute data to the send buffer */
static void
t8_cmesh_partition_copy_data (char *send_buffer,
                              const struct t8_cmesh * cmesh_from,
                              t8_locidx_t num_trees, size_t attr_info_bytes,
                              size_t attr_bytes,
                              size_t ghost_neighbor_bytes,
                              size_t tree_neighbor_bytes,
                              sc_array_t * send_as_ghost,
                              t8_locidx_t send_first, t8_locidx_t send_last)
{
  t8_ctree_t          tree, tree_cpy;
  size_t              temp_offset_tree, temp_offset_att, iz, temp_offset;
  ssize_t             last_attribute_diff;
  t8_attribute_info_struct_t * attr_info;
  t8_locidx_t         num_ghost_send = send_as_ghost->elem_count, itree;
  t8_locidx_t        *face_neighbor, ghost_id;
  t8_gloidx_t        *face_neighbor_g, *face_neighbor_gnew, new_neighbor;
  int8_t             *ttf;
  t8_cghost_t         ghost, ghost_cpy;
  int                 iface;

  /* Copy all trees to the send buffer */
  tree = t8_cmesh_trees_get_tree_ext (cmesh_from->trees, send_first,
                                      &face_neighbor, &ttf);
  (void) memcpy (send_buffer, tree, num_trees * sizeof (t8_ctree_struct_t));
  /* Copy all face neighbor information to send_buffer */
  (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                 num_ghost_send * sizeof (t8_cghost_struct_t) +
                 ghost_neighbor_bytes, face_neighbor, tree_neighbor_bytes);
  /* Copy all attribute infos to send_buffer */
  (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                 num_ghost_send * sizeof (t8_cghost_struct_t) +
                 ghost_neighbor_bytes + tree_neighbor_bytes,
                 (char *)tree + tree->att_offset, attr_info_bytes);
  /* Copy all attribute data to send_buffer */
  (void) memcpy (send_buffer + num_trees * sizeof (t8_ctree_struct_t) +
                 num_ghost_send * sizeof (t8_cghost_struct_t) +
                 ghost_neighbor_bytes + tree_neighbor_bytes +
                 attr_info_bytes,
                 T8_TREE_ATTR(tree,T8_TREE_ATTR_INFO(tree, 0)), attr_bytes);
  /* Set new face_neighbor offsets */
  temp_offset = num_trees * sizeof (t8_ctree_struct_t) +
      num_ghost_send * sizeof (t8_cghost_struct_t) +
      ghost_neighbor_bytes; /* Computes the offset of the face neighbors of the new trees */
  temp_offset_att = attr_info_bytes; /* Compute the offset of the new attribute infos */
  temp_offset_tree = 0;
  attr_info = (t8_attribute_info_struct_t *)(send_buffer +
                  num_trees * sizeof (t8_ctree_struct_t) +
                  num_ghost_send * sizeof (t8_cghost_struct_t) +
                  ghost_neighbor_bytes + tree_neighbor_bytes);
  /* This has to be add to each attribute info offset to calculated the new offset */
  last_attribute_diff = attr_info_bytes - attr_info->attribute_offset;
  /* Set attribute offsets of trees and attribute data offsets of info objects */
  for (itree = send_first; itree <= send_last; itree++) {
    tree_cpy = (t8_ctree_t) (send_buffer + temp_offset_tree);
    /* new neighbor offset of tree */
    tree_cpy->neigh_offset = temp_offset - temp_offset_tree;
    /* compute neighbor offset for next tree */
    temp_offset += t8_eclass_num_faces[tree_cpy->eclass] *
        (sizeof (t8_locidx_t) + sizeof (int8_t))
        + ((4 - t8_eclass_num_faces[tree_cpy->eclass] *
        (sizeof (t8_locidx_t) + sizeof (int8_t)) % 4) %4);
    /* new attribute offset for tree */
    tree_cpy->att_offset = temp_offset_att - temp_offset_tree;
    /* set new attribtue data offsets */
    for (iz = 0;iz < tree_cpy->num_attributes;iz++) {
      attr_info++;
      attr_info->attribute_offset += last_attribute_diff;
    }
    temp_offset_att += sizeof (t8_attribute_info_struct_t);
    temp_offset_tree += sizeof (t8_ctree_t);
  }
  /* Copy all ghosts and set their face entries and offsets */
  /* Offset of ghost face_neighbor from first ghost */
  temp_offset = num_trees * sizeof (t8_ctree_struct_t) +
      num_ghost_send * sizeof (t8_cghost_struct_t);
  /* Offset of current ghost from first ghost */
  temp_offset_tree = 0;
  for (iz = 0;iz < send_as_ghost->elem_count;iz++) {
    ghost_id = *((t8_locidx_t *) sc_array_index (send_as_ghost,iz));
    ghost_cpy = (t8_cghost_t) (send_buffer +
                               num_trees * sizeof (t8_ctree_struct_t)
                               + iz * sizeof (t8_cghost_struct_t));
    ghost_cpy->neigh_offset = temp_offset - temp_offset_tree; /* New face neighbor offset */
    face_neighbor_gnew = (t8_gloidx_t *) T8_GHOST_FACE(ghost_cpy);
    if (ghost_id >= cmesh_from->num_local_trees) {
      /* The ghost that we send was a local ghost */
      ghost = t8_cmesh_trees_get_ghost_ext (cmesh_from->trees,
                                            ghost_id -
                                            cmesh_from->num_local_trees,
                                            &face_neighbor_g);
      ghost_cpy->eclass = ghost->eclass;
      ghost_cpy->treeid = ghost->treeid;
      /* Copt face_neighbor entries */
      memcpy (face_neighbor_gnew, face_neighbor_g,
              t8_eclass_num_faces[ghost_cpy->eclass] * sizeof(t8_gloidx_t));
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
      for (iface = 0; iface < t8_eclass_num_faces[ghost_cpy->eclass];iface++) {
        if (face_neighbor[iface] < 0) {
          /* TODO: think about this */
          new_neighbor = -1; /* boundary indicator */
        }
        else {
          /* Compute global index from local index */
          new_neighbor = face_neighbor[iface] < cmesh_from->num_local_trees ?
                face_neighbor[iface] + cmesh_from->first_tree :
                face_neighbor[iface] - cmesh_from->num_local_trees
                + cmesh_from->first_tree;
        }
        face_neighbor_gnew[iface] = new_neighbor;
      }
    }
    /* compute new offsets */
    temp_offset = t8_eclass_num_faces[ghost_cpy->eclass] * sizeof (t8_gloidx_t);
    temp_offset_tree += sizeof (t8_cghost_struct_t);
  }
}

static t8_gloidx_t
t8_glo_min (t8_gloidx_t A, t8_gloidx_t B)
{
  return A < B ? A : B;
}

static void
t8_cmesh_partition_given (t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from,
                          t8_gloidx_t * tree_offset)
{
  int                 send_first, send_last;    /* ranks of the processor to which we will send */
  int                 recv_first, recv_last, num_receive;    /* ranks of the processor from which we will receive */
  int                 iproc, iface, mpiret, proc_recv, recv_bytes;
  char               *send_buffer = NULL, *my_buffer = NULL;
  size_t              attr_bytes = 0, tree_neighbor_bytes,
      ghost_neighbor_bytes, total_alloc, attr_info_bytes;
  sc_MPI_Request     *requests;
  sc_MPI_Status      *status;
  t8_locidx_t         range_start, range_end, num_ghost_send, itree,
    num_trees, *face_neighbor;
  t8_part_tree_t      recv_part;
  int8_t             *ttf;
  t8_locidx_t         neighbor;
  t8_ctree_t          tree;
  sc_array_t          send_as_ghost;    /* Stores local id's of trees and ghosts that will be send as ghosts */
  sc_array_t          keep_as_ghost;    /* Store local id's of local trees and ghosts that will be ghost on this process */
  int8_t             *ghost_flag;       /* For each local tree and ghost set to 1 if it is in send_as_ghost
                                           and set to 3 if it is in keep_as_ghost */
  /* TODO: computing recv information needs the shared array of the old partition in cmesh_from,
   *       thus, two of those huge arrays need to exist at the same time.
   *       Using MPI-2.1 One-sided communication we could resolve this, since
   *       pi \in pj_recv if and only if pj \in pi_send. The latter can be computed
   *       using only the new partition array w/o communication and then via one-sided
   *       communication we can compute the number of receiving processes on each process, which
   *       should be enough to receive all messages in a while loop. */

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (cmesh_from != NULL);
  T8_ASSERT (cmesh_from->committed);

  ghost_flag = T8_ALLOC (int8_t, cmesh_from->num_local_trees +
                         cmesh_from->num_ghosts);
  sc_array_init (&send_as_ghost, sizeof (t8_locidx_t));
  sc_array_init (&keep_as_ghost, sizeof (t8_locidx_t));

  /* determine send and receive range. temp_tree is last local tree of send_first in new partition */
  t8_cmesh_partition_recvrange (cmesh, (t8_cmesh_t) cmesh_from, &recv_first,
                                &recv_last);
  range_end = t8_cmesh_partition_sendrange (cmesh, (t8_cmesh_t) cmesh_from,
                                            &send_first, &send_last);
  /* range_end stores (my rank) local tree_id of last tree on send_first */

  range_start = 0; /* Stores the first tree that was not send yet */
  requests = T8_ALLOC (sc_MPI_Request, send_last - send_first + 1);
  for (iproc = send_first; iproc <= send_last; iproc++) {
    attr_bytes = 0;
    tree_neighbor_bytes = 0;
    attr_info_bytes = 0;
    ghost_neighbor_bytes = 0;
    /* iproc != mpirank */

    memset (ghost_flag, 0,
            (cmesh_from->num_local_trees + cmesh_from->num_ghosts)
            * sizeof (int8_t));       /* Yes, i know that sizeof(int8_t) is always 1, so what?! */
    sc_array_truncate (&send_as_ghost);
    /* loop over all trees that will be send */
    for (itree = range_start; itree <= range_end; itree++) {
      tree= t8_cmesh_trees_get_tree_ext (cmesh_from->trees, itree,
                                         &face_neighbor, &ttf);
      /* Count the additional memory needed per tree from neighbors */
      tree_neighbor_bytes += t8_eclass_num_faces[tree->eclass] *
        (sizeof (*face_neighbor) + sizeof (*ttf));
      tree_neighbor_bytes += (4-tree_neighbor_bytes %4) %4; /* padding to make number of bytes per tree
                                                              a multiple of 4 */
      /*  Compute number of attribute bytes in this tree range.
       *       Not every tree has an attribute */
      attr_info_bytes += tree->num_attributes *
          sizeof (t8_attribute_info_struct_t);
      attr_bytes += t8_cmesh_trees_attribute_size (tree);

      /* loop over all faces of each tree to determine ghost to send */
      for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
        neighbor = face_neighbor[iface];
        if (neighbor >= 0 && neighbor < cmesh_from->num_local_trees &&
            (neighbor < range_start || neighbor > range_end)) {
          /* neighbor is a local tree or local ghost and will be ghost on iproc */
          if (ghost_flag[neighbor] == 0 &&
              t8_cmesh_send_ghost (cmesh, iproc, neighbor)) {
            /* we did not add this neighbor yet and it should be send to iproc */
            ghost_flag[neighbor] = 1;
            *((t8_locidx_t *) sc_array_push (&send_as_ghost)) = neighbor;
          }
        }
        else if (cmesh->first_tree <= neighbor + cmesh_from->first_tree
                 && neighbor + cmesh_from->first_tree <= cmesh->first_tree +
                 cmesh->num_local_trees - 1) {
          /* neighbor will be local on this process, thus tree will be ghost on this process */
          if (ghost_flag[itree] == 0) {
            ghost_flag[itree] = 3;
            *((t8_locidx_t *) sc_array_push (&keep_as_ghost)) = itree;
          }
        }
      }
    }
    /* loop over trees ends here */
    /* Calculate and allocate memory for send buffer */
    /**********************************************************/
    /*
     *      The data that we send has the layout
     *
     * | tree 0,1...| ghost 0,1...| face_neighbors0, tree_to_face0,...| ghost_neighbors 0,1,...| Attinfo00,01,...,Attinfoend | attdata00,01,... |
     *
     */
    attr_info_bytes += sizeof (t8_attribute_info_struct_t); /* We always put one attribute info struct at the end */
    num_trees = range_end - range_start + 1;
    num_ghost_send = send_as_ghost.elem_count;
    /* TODO: parse throgh send_as_ghost and add to ghost_neighbor_bytes */
    total_alloc = num_trees * sizeof (t8_ctree_struct_t) +
      num_ghost_send * sizeof (t8_cghost_struct_t) + ghost_neighbor_bytes +
        tree_neighbor_bytes + attr_info_bytes + attr_bytes;
    /* Extra space to store the number of trees and ghosts in the send buffer */
    total_alloc += 2 * sizeof (t8_locidx_t);
#if SC_ENABLE_USE_REALLOC
    send_buffer = SC_REALLOC (send_buffer, char, total_alloc);
#else
    T8_FREE (send_buffer);
    send_buffer = T8_ALLOC (char, total_alloc);
#endif
    /* Store number of trees and ghosts at the end of send buffer */
    *((t8_locidx_t *) send_buffer + total_alloc - 2 * sizeof(t8_locidx_t))
        = num_trees;
    *((t8_locidx_t *) send_buffer + total_alloc - sizeof(t8_locidx_t))
        = num_ghost_send;
    /* Copy all data to the send buffer */
    t8_cmesh_partition_copy_data (send_buffer, cmesh_from, num_trees,
                                  attr_info_bytes, attr_bytes,
                                  ghost_neighbor_bytes,
                                  tree_neighbor_bytes, &send_as_ghost,
                                  send_first, send_last);
    if (iproc == cmesh->mpirank) {
      /* TODO: */
      /* We do not need to send the trees since we already own them */
      /* and act as if we just received the send_buffer */
      my_buffer = T8_ALLOC (char, total_alloc);
      memcpy (my_buffer, send_buffer, total_alloc);
    }
    else {
      /* send buffer to remote process */
      mpiret = sc_MPI_Isend (send_buffer, total_alloc, sc_MPI_BYTE, iproc, 0,
                             cmesh->mpicomm, requests + iproc);
      SC_CHECK_MPI (mpiret);
    }
    /* compute new ranges here */
    /* If tree_offset of iproc + 1 is < 0 the last tree is shared and has to
     * be send to the next process as well */
    range_start = range_end + 1 - (tree_offset[iproc+1] < 0);
    /* add number of trees in iproc + 1 to send to range_end */
    /* We have to be careful with locidx overflow when we go out of bounds
     * of our process */
    range_end = t8_glo_min (range_end + t8_glo_abs (tree_offset[iproc + 1]) -
        t8_glo_abs (tree_offset[iproc]) + tree_offset[iproc + 1] < 0,
        cmesh_from->num_trees - 1);
  } /* sending loop ends here */
  T8_FREE (ghost_flag);
  sc_array_reset (&send_as_ghost);

  num_trees = t8_glo_abs (tree_offset[cmesh->mpirank + 1]) -
      t8_glo_abs (tree_offset[cmesh->mpirank])
      + tree_offset[cmesh->mpirank] < 0;
  /* Initialize trees structure with yet unknown number of ghosts */
  t8_cmesh_trees_init (&cmesh->trees, recv_last - recv_first + 1, num_trees,
                       0);
  /* Receive from other processes */
  /* Total number of processor from which we receive is one less if the current
   * rank is included */
  num_receive = recv_last + 1 - recv_first - (recv_first <= cmesh->mpirank
                                              && cmesh->mpirank <= recv_last);
  status = T8_ALLOC (sc_MPI_Status, num_receive);
  for (iproc = 0;iproc < num_receive;iproc++) {
    /* Blocking test for message. We do not know from which process or how big it is */
    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, 0, cmesh->mpicomm, status + iproc);
    SC_CHECK_MPI (mpiret);
    proc_recv = status[iproc].MPI_SOURCE;
    T8_ASSERT (status[iproc].MPI_TAG == 0);
    T8_ASSERT (recv_first <= proc_recv && proc_recv <= recv_last);
    mpiret = sc_MPI_Get_count (status + iproc, sc_MPI_BYTE, &recv_bytes);
    SC_CHECK_MPI (mpiret);
    /* Allocate receive buffer */
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, proc_recv - recv_first);
    /* take first tree of part and allocate recv_bytes */
    recv_part->first_tree = T8_ALLOC (char, recv_bytes);
    /* Receive message */
    mpiret = sc_MPI_Recv (recv_part->first_tree, recv_bytes, sc_MPI_BYTE,
                          proc_recv, 0, cmesh->mpicomm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    /* Read num trees and num ghosts */
    recv_part->num_trees = (t8_locidx_t) recv_part->first_tree[recv_bytes -
        2*sizeof(t8_locidx_t)];
    recv_part->num_ghosts = (t8_locidx_t) recv_part->first_tree[recv_bytes -
        sizeof(t8_locidx_t)];
    /* Free memory that was only used to store num_trees/ghosts */
    T8_FREE (recv_part->first_tree + recv_bytes - 2 * sizeof (t8_locidx_t));
  }
  /* Got trees and ghosts from myself */
  if (my_buffer != NULL) {
    recv_part = t8_cmesh_trees_get_part (cmesh->trees,
                                         cmesh->mpirank - recv_first);
    recv_part->first_tree = my_buffer;
    /* Read num trees and num ghosts */
    recv_part->num_trees = (t8_locidx_t) recv_part->first_tree[recv_bytes -
        2*sizeof(t8_locidx_t)];
    recv_part->num_ghosts = (t8_locidx_t) recv_part->first_tree[recv_bytes -
        sizeof(t8_locidx_t)];
    /* Free memory that was only used to store num_trees/ghosts */
    T8_FREE (recv_part->first_tree + recv_bytes - 2 * sizeof (t8_locidx_t));
  }
  /* set recv_part->first_tree_id/first_ghost_id */
  /* Use as temporary variables to store the first tree_id/ghost_id of the new parts */
  num_trees = 0;
  num_ghost_send = 0;
  for (iproc = 0;iproc < cmesh->trees->from_proc->elem_count;iproc++) {
    recv_part = t8_cmesh_trees_get_part (cmesh->trees, iproc);
    recv_part->first_tree_id = num_trees;
    recv_part->num_ghosts = num_ghost_send;
    num_trees += recv_part->num_trees;
    num_ghost_send += recv_part->num_ghosts;
  }
  /* TODO: set new local ids of face_neighbors */
}

void
t8_cmesh_partition (t8_cmesh_t cmesh)
{
  t8_cmesh_t          cmesh_from;
  t8_gloidx_t         last_tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (cmesh->set_from != NULL);
  T8_ASSERT (cmesh->set_from->committed);

  cmesh_from = (t8_cmesh_t) cmesh->set_from;
  /**********************************************/
  /*      Compute local number of trees         */
  /*         and trees per proc array           */
  /**********************************************/
  if (cmesh->set_level >= 0) {
    /* Compute first and last tree index */
    T8_ASSERT (cmesh->tree_offsets == NULL);
    t8_cmesh_uniform_bounds (cmesh_from, cmesh->set_level, &cmesh->first_tree,
                             NULL, &last_tree, NULL,
                             &cmesh->first_tree_shared);
    cmesh->num_local_trees = last_tree - cmesh->first_tree;
    /* To compute the tree_offsets correctly we have to invert the sign on the
     * first tree if last tree is shared, since we use it for MPI to allgather the tree_offsets */
    if (cmesh->first_tree_shared) {
      cmesh->first_tree = -cmesh->first_tree;
    }
    /* allocate and fill tree_offset array with number of trees per process */
    /* We have to be careful with shared memory in the case where cmesh comm is
     * duplicated, since we have to use the same communicator for allocating and freeing memory.
     * Thus this function must only be called after cmesh communicator was duplicated,
     * so we check whether mpisize has been set */
    T8_ASSERT (cmesh->mpisize > 0);
    cmesh->tree_offsets = (t8_gloidx_t *) SC_SHMEM_ALLOC (t8_locidx_t,
                                                          cmesh->mpisize + 1,
                                                          cmesh->mpicomm);
    sc_shmem_allgather (&cmesh->first_tree, 1, T8_MPI_GLOIDX,
                        cmesh->tree_offsets, 1, T8_MPI_LOCIDX,
                        cmesh->mpicomm);
    cmesh->tree_offsets[cmesh->mpisize] = cmesh_from->num_trees;
    /* tree_offsets was computed, reinvert the sign */
    if (cmesh->first_tree_shared) {
      T8_ASSERT (cmesh->first_tree <= 0);
      cmesh->first_tree = -cmesh->first_tree;
    }
  }
  else {
    /* We compute the partition after a given partition table in cmesh->tree_offsets */
    T8_ASSERT (cmesh->tree_offsets != NULL);
    cmesh->first_tree_shared = cmesh->tree_offsets[cmesh->mpirank] < 0;
    /* compute local first tree */
    cmesh->first_tree = t8_glo_abs (cmesh->tree_offsets[cmesh->mpirank]);
    /* compute local num trees */
    cmesh->num_local_trees = t8_glo_abs (cmesh->tree_offsets[cmesh->mpirank + 1])
        - cmesh->first_tree + cmesh->tree_offsets[cmesh->mpirank + 1] < 0;
  }
  /***************************************************/
  /*        Done with local num and tree_offset      */
  /***************************************************/
  t8_cmesh_partition_given (cmesh, cmesh->set_from, cmesh->tree_offsets);
}
