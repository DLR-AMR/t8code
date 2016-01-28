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

#if 0

/* Compute first and last process to which we will send data */
/* Returns the local tree_id of last local tree on send_first. */
/* TODO: This can propably be optimized since we changed to tree_offsets. using binary search or so */
static t8_locidx_t
t8_cmesh_partition_sendrange (t8_cmesh_t cmesh, t8_cmesh_t cmesh_from,
                              int *send_first, int *send_last)
{
  int                 iproc;
  int                 last_tree_shared;
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
    last_tree_shared = cmesh->tree_offsets[iproc] < 0;
    first_tree = t8_glo_abs (cmesh->tree_offsets[iproc]);
    last_tree = t8_glo_abs (cmesh->tree_offsets[iproc+1]) - 1 + last_tree_shared;

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

static void
t8_cmesh_partition_given (t8_cmesh_t cmesh, const struct t8_cmesh *cmesh_from,
                          t8_gloidx_t * tree_offset)
{
  int                 send_first, send_last;    /* ranks of the processor to which we will send */
  int                 recv_first, recv_last;    /* ranks of the processor from which we will receive */
  int                 iproc, iface, F;
  char               *send_buffer = NULL;
  size_t              attr_bytes =
    0, tree_additional_bytes, ghost_neighbor_bytes, total_alloc, count, iz;
  t8_attribute_info_struct_t *attr_info;
  t8_locidx_t         range_start, range_end, num_ghost_send, itree,
    num_trees, *face_neighbor;
  int8_t             *ttf;
  t8_locidx_t         neighbor;
  t8_ctree_t          tree, tree_cpy;
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

  range_start = 0;
  for (iproc = send_first; iproc <= send_last; iproc++) {
    attr_bytes = 0;
    tree_additional_bytes = 0;
    ghost_neighbor_bytes = 0;
    if (iproc == cmesh->mpirank) {
      /* We do not need to send the trees since we already own them */
    }
    else {
      /* iproc != mpirank */

      memset (ghost_flag, 0,
              (cmesh_from->num_local_trees + cmesh_from->num_ghosts)
              * sizeof (int8_t));       /* Yes, i know that sizeof(int8_t) is always 1, so what?! */
      sc_array_truncate (&send_as_ghost);
      /* loop over all trees that will be send */
      for (itree = 0; itree <= range_end; itree++) {
        tree = t8_cmesh_trees_get_tree (cmesh_from->trees, itree);
        /* Count the additional memory needed per tree from neighbors */
        tree_additional_bytes += t8_eclass_num_faces[tree->eclass] *
          (sizeof (*tree->face_neighbors) + sizeof (*tree->tree_to_face));        
        /*  Compute number of attribute bytes in this tree range.
         *       Not every tree has an attribute */
        if (tree->attributes != NULL) {
          tree_additional_bytes += tree->attributes->elem_count *
            sizeof (t8_attribute_info_struct_t);
          for (iz = 0; iz < tree->attributes->elem_count;iz++) {
            attr_info = (t8_attribute_info_struct_t *)
                sc_array_index (tree->attributes, iz);
            attr_bytes += attr_info->attribute_size;
          }
        }
        /* loop over all faces of each tree to determine ghost to send */
        for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
          neighbor = tree->face_neighbors[iface];
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
       * | tree 0,1...| ghost 0,1...| attribute 0,1...| (face_neighbors0, tree_to_face0,attribute_info00,01,...),...| ghost_neighbors 0,1,...|
       *
       *      The only thing that has to be allocated newly on receiver is metadata for the sc_array
       *      holding the attribute_info's.
       */
      num_trees = range_end - range_start + 1;
      num_ghost_send = send_as_ghost.elem_count;
      /* TODO: parse throgh send_as_ghost and add to ghost_neighbor_bytes */
      total_alloc = num_trees * sizeof (t8_ctree_struct_t) +
        num_ghost_send * sizeof (t8_cghost_struct_t) + attr_bytes +
        tree_additional_bytes + ghost_neighbor_bytes;
#if SC_ENABLE_USE_REALLOC
      send_buffer = SC_REALLOC (send_buffer, char, total_alloc);
#else
      T8_FREE (send_buffer);
      send_buffer = T8_ALLOC (char, total_alloc);
#endif
      count = 0; /* Counts the bytes from the beginning of send_buffer that we
                   already filled */
      ghost_neighbor_bytes = num_trees * sizeof (t8_ctree_struct_t) +
          num_ghost_send * sizeof (t8_cghost_struct_t) + attr_bytes +
          tree_additional_bytes; /* Offset of fisrt not yet copied ghost neighbors */
      tree_additional_bytes = ghost_neighbor_bytes - tree_additional_bytes;
      /* Offset of first, not yet copied tree face_neighbor/ttf/attribute */
      attr_bytes = tree_additional_bytes - attr_bytes; /* Offset of first non-set tree_attribute */
      for (itree = send_first; itree <= send_last; itree++) {
        /* copy trees to send_buffer here */
        tree = t8_cmesh_trees_get_tree (cmesh->trees, itree, face_neighbor,
                                        ttf);
        (void) memcpy (send_buffer + count, tree, sizeof (t8_ctree_struct_t));
        count += sizeof (t8_ctree_struct_t);
        tree_cpy = (t8_ctree_t) (send_buffer + count);
        F = t8_eclass_num_faces[tree->eclass];
        /* Copy face_neighbors */
        memcpy (send_buffer + tree_additional_bytes, face_neighbor,
                F * sizeof (t8_locidx_t));
        tree_cpy->f_neigh_offset = tree_additional_bytes;
        tree_additional_bytes += F * sizeof (t8_locidx_t);
        /* Copy tree_to_face */
        memccpy (send_buffer + tree_additional_bytes, ttf, F * sizeof (int8_t));
        tree_cpy->ttf_offset = tree_additional_bytes;
        tree_additional_bytes += F * sizeof (int8_t);
        /* Copy attributes (if exist) */
        if (tree->attributes != NULL) {
          for (iz = 0; iz < tree->attributes->elem_count;iz++) {
            attr_info = sc_array_index (tree->attributes, iz);
            t8_cmesh_trees_get_attribute (cmesh->trees, itree, a)
            memccpy ()
          }
        }
    }
    /* compute new ranges here */
  }                             /* sending loop ends here */
  T8_FREE (ghost_flag);
  sc_array_reset (&send_as_ghost);
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
                             &cmesh->last_tree_shared);
    cmesh->num_local_trees = last_tree - cmesh->first_tree;
    /* To compute the tree_offsets correctly we have to invert the sign on the
     * first tree if last tree is shared, since we use it for MPI to allgather the tree_offsets */
    if (cmesh->last_tree_shared) {
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
    if (cmesh->last_tree_shared) {
      T8_ASSERT (cmesh->first_tree <= 0);
      cmesh->first_tree = -cmesh->first_tree;
    }
  }
  else {
    /* We compute the partition after a given partition table in cmesh->tree_offsets */
    T8_ASSERT (cmesh->tree_offsets != NULL);
    cmesh->last_tree_shared = cmesh->tree_offsets[cmesh->mpirank] < 0;
    /* compute local first tree */
    cmesh->first_tree = cmesh->last_tree_shared ?
      -cmesh->tree_offsets[cmesh->mpirank] :
        cmesh->tree_offsets[cmesh->mpirank];
    /* compute local num trees */
    cmesh->num_local_trees = t8_glo_abs (cmesh->tree_offsets[cmesh->mpirank + 1])
        - cmesh->first_tree + cmesh->last_tree_shared;
  }
  /***************************************************/
  /*        Done with local num and tree_offset      */
  /***************************************************/
  t8_cmesh_partition_given (cmesh, cmesh->set_from, cmesh->tree_offsets);
}
#endif
