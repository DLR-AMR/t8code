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

/** \file t8_cmesh_commit.c
 *
 * TODO: document this file
 */

#include <t8_shmem.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"

struct ghost_facejoins_struct
{
  t8_gloidx_t         ghost_id; /* The id of the ghost */
  size_t              index;    /* The index in stash->joinfaces where this ghost is used */
  int8_t              flag;     /* First bit is true if we are certain if entry is a real ghost,
                                 * but may also be false in this case.
                                 * First but is guaranteed to be false if entry is only ghost of ghost.
                                 * Second bit is true if and only if this entry comes from the second
                                 * id member of joinfaces[index]. */
};

static int
t8_ghost_facejoins_compare (const void *fj1, const void *fj2)
{
  struct ghost_facejoins_struct *FJ1 = (struct ghost_facejoins_struct *) fj1;
  struct ghost_facejoins_struct *FJ2 = (struct ghost_facejoins_struct *) fj2;

  return FJ1->ghost_id < FJ2->ghost_id ? -1 : FJ1->ghost_id != FJ2->ghost_id;
}

/* requires the attributes array on the stash to be sorted */
static void
t8_cmesh_add_attributes (t8_cmesh_t cmesh, t8_stash_t stash)
{
  size_t              attr_bytes, attr_offset, num_attr, attr_index, si;
  int                 attr_id, key;
  t8_locidx_t         tree_id, temptree;

  attr_offset = 0;
  temptree = -1;
  /* TODO: optimize loop by looping only over local trees */
  for (si = 0; si < stash->attributes.elem_count; si++) {
    tree_id = t8_stash_get_attribute_tree_id (cmesh->stash, si) -
      cmesh->first_tree;
    T8_ASSERT (t8_stash_get_attribute_tree_id (cmesh->stash, si) -
               cmesh->first_tree == (t8_gloidx_t) tree_id);
    /* Only add local trees */
    if (0 <= tree_id && tree_id < cmesh->num_local_trees) {
      if (tree_id > temptree) {
        /* We are entering a new tree */
        temptree = tree_id;
        /* count the number of attributes at this tree... */
        num_attr = si;
        /* ...by counting after how many attributes we enter a new tree */
        while (temptree == tree_id
               && ++num_attr < stash->attributes.elem_count) {
          temptree =
            t8_stash_get_attribute_tree_id (cmesh->stash,
                                            num_attr) - cmesh->first_tree;
        }
        num_attr -= si;         /* now stores the number of attribute of the tree */
        /* initialize storage for tree attributes */
        t8_cmesh_trees_init_attributes (cmesh->trees, tree_id, num_attr);
        temptree = tree_id;     /* store the current tree_id in oldtree */
        attr_index = 0;
      }
      attr_bytes = t8_stash_get_attribute_size (cmesh->stash, si);
      key = t8_stash_get_attribute_key (cmesh->stash, si);
      attr_id = t8_stash_get_attribute_id (cmesh->stash, si);
      t8_cmesh_tree_add_attribute (cmesh->trees, 0, tree_id, attr_id, key,
                                   (char *)
                                   t8_stash_get_attribute (cmesh->stash,
                                                           si), attr_bytes,
                                   attr_offset, attr_index);
      attr_offset += attr_bytes;
      attr_index++;
    }
  }
}

static void
t8_cmesh_set_shmem_type (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);

  sc_shmem_set_type (cmesh->mpicomm, T8_SHMEM_BEST_TYPE);
}

/* TODO: set boundary face connections here.
 *       not trivial if replicated and not level 3 face_knowledg
 *       Edit: boundary face is default. If no face-connection is added then
 *             we assume a boundary face.
 * TODO: Implement a debug check for mesh consistency between processes.
 */
void
t8_cmesh_commit (t8_cmesh_t cmesh)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->mpicomm != sc_MPI_COMM_NULL);
  T8_ASSERT (!cmesh->committed);

  /* TODO: setup trees */
  if (!cmesh->set_partitioned) {
    if (cmesh->stash != NULL && cmesh->stash->classes.elem_count > 0) {
      t8_stash_t          stash = cmesh->stash;
      sc_array_t         *class_entries = &stash->classes;
      t8_stash_class_struct_t *entry;
      t8_topidx_t         num_trees = class_entries->elem_count, itree;
      size_t              si, attr_bytes;

      t8_cmesh_trees_init (&cmesh->trees, 1, num_trees, 0);
      /* compute size of attributes */
      attr_bytes = 0;
      for (si = 0; si < stash->attributes.elem_count; si++) {
        attr_bytes += t8_stash_get_attribute_size (stash, si);
      }
      cmesh->num_trees = cmesh->num_local_trees = num_trees;
      cmesh->first_tree = 0;
      t8_cmesh_trees_init_part (cmesh->trees, 0, 0, num_trees, 0, attr_bytes);
      /* set tree classes */
      for (itree = 0; itree < num_trees; itree++) {
        entry = (t8_stash_class_struct_t *)
          t8_sc_array_index_topidx (class_entries, itree);
        t8_cmesh_trees_add_tree (cmesh->trees, entry->id, 0, entry->eclass);
        cmesh->num_trees_per_eclass[entry->eclass]++;
      }
      /* set tree attributes */
      /* TODO: replace attribute sort by bucket sort into tree structs +
       *       sorting inside the tree structs by key.
       *       This will bring down the runtime of this step from O(nlog(n)) to
       *       O(nm) where m is the maximum number of attributes under all trees.
       */
      t8_stash_attribute_sort (cmesh->stash);
      t8_cmesh_add_attributes (cmesh, cmesh->stash);
    }
  }
  else {
    sc_array_t         *ghost_ids;
    size_t              joinfaces_it, attr_byte_count, iz;
    size_t             *ghost_offsets;
    ssize_t             ghost_ind, jz, class_end;
    t8_stash_joinface_struct_t *joinface;
    t8_gloidx_t         last_tree = cmesh->num_local_trees +
      cmesh->first_tree - 1, id1, id2;
    t8_stash_attribute_struct_t *attribute;
    t8_stash_class_struct_t *classentry;
    int                 id1_istree, id2_istree, F, face, face2, temp;
    t8_ctree_t          tree1;
    t8_cghost_t         ghost1;
    struct ghost_facejoins_struct *ghost_facejoin;
    int8_t              is_true_ghost;

    if (cmesh->face_knowledge != 3) {
      t8_global_errorf ("Expected a face knowledge of 3.\nAbort commit.");
      /* TODO: reset cmesh */
      return;
    }
    t8_cmesh_set_shmem_type (cmesh);    /* TODO: do we actually need the shared array? */
    t8_stash_class_sort (cmesh->stash);
//    t8_stash_joinface_sort (cmesh->stash); /* TODO: this is propably not usefull */
    t8_stash_attribute_sort (cmesh->stash);

    ghost_ids = sc_array_new (sizeof (struct ghost_facejoins_struct));
    /* Parse joinfaces array and save all indices of entries with ghosts involved */
    for (joinfaces_it = 0; joinfaces_it < cmesh->stash->joinfaces.elem_count;
         joinfaces_it++) {
      joinface =
        (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->stash->
                                                       joinfaces,
                                                       joinfaces_it);
      id1 = joinface->id1;
      id2 = joinface->id2;
      id2_istree = id2 <= last_tree && id2 >= cmesh->first_tree;
      id1_istree = id1 <= last_tree && id1 >= cmesh->first_tree;
      if (!id2_istree) {
        /* id2 is a ghost */
        ghost_facejoin =
          (struct ghost_facejoins_struct *) sc_array_push (ghost_ids);
        ghost_facejoin->ghost_id = id2;
        ghost_facejoin->index = joinfaces_it;
        ghost_facejoin->flag = id1_istree | 2;  /* The 2 specifies that we come from id2 entry */
      }
      if (!id1_istree) {
        /* id1 is a ghost */
        ghost_facejoin =
          (struct ghost_facejoins_struct *) sc_array_push (ghost_ids);
        ghost_facejoin->ghost_id = id1;
        ghost_facejoin->index = joinfaces_it;
        ghost_facejoin->flag = id2_istree;
      }
      /* For a true ghost at least one entry will have flag = 1 */
    }

    /* Added one ghostid entry for each face connection with a ghost and ghost of a ghost.
     * Need to count w/o duplicates.
     * Need to know which ones are proper ghosts */
    /* TODO: this step is the most expensive one in this algorihtm.
     * ~5-6 times slower than the sections above and below each.
     * Also this step scales with the number of non-local trees (not ghosts!)
     * which scales with the number of processes and approximates O(global trees).
     * Would be much faster if we only had ghosts here, since ghost = O(trees)
     */
    sc_array_sort (ghost_ids, t8_ghost_facejoins_compare);
    /* Count the ghosts without duplicates */
    /* TODO: we could also reduce the ghost_ids list and only keep entries with proper ghosts */
    id1 = -1;
    is_true_ghost = 0;

    for (cmesh->num_ghosts = 0, jz = -1, iz = 0; iz < ghost_ids->elem_count;
         iz++) {
      ghost_facejoin =
        (struct ghost_facejoins_struct *) sc_array_index (ghost_ids, iz);
      temp = ghost_facejoin->flag & 1;
      if (ghost_facejoin->ghost_id > id1) {
        id1 = ghost_facejoin->ghost_id;
        /* Check if the current global id belongs to a ghost */
        if (jz >= 0 && is_true_ghost == 1) {
          /* If so go back to first entry of this id and set flag to 1.
           * Count this entry as ghost. */
          ghost_facejoin =
            (struct ghost_facejoins_struct *) sc_array_index (ghost_ids, jz);
          ghost_facejoin->flag |= 1;
          cmesh->num_ghosts++;
        }
        jz = iz;                /* Store first index of new id */
        is_true_ghost = 0;
      }
      is_true_ghost |= temp;        /* is_true_ghost gets true if any flag's
                                                           first bit is true. */
    }
    /* This is to catch the last entry in the array */
    if (jz >= 0 && is_true_ghost == 1) {
      /* If so go back to first entry of this id and set flag to 1.
       * Count this entry as ghost. */
      ghost_facejoin =
        (struct ghost_facejoins_struct *) sc_array_index (ghost_ids, jz);
      ghost_facejoin->flag |= 1;
      cmesh->num_ghosts++;
    }
    /* Now that we know the number of ghosts, we can save the offsets of
     * the true ghosts into the ghost_ids array */
    ghost_offsets = T8_ALLOC (size_t, cmesh->num_ghosts);
    id1 = -1;
    for (iz = 0, jz = 0; iz < ghost_ids->elem_count;iz++) {
      ghost_facejoin =
        (struct ghost_facejoins_struct *) sc_array_index (ghost_ids, iz);
      if (ghost_facejoin->ghost_id > id1 && ghost_facejoin->flag & 1) {
        /* If we enter a new real ghost, store its offset */
        ghost_offsets[jz] = iz;
        id1 = ghost_facejoin->ghost_id;
        jz++;
      }
    }
    T8_ASSERT (jz == cmesh->num_ghosts);

    /* Count attribute bytes */
    attr_byte_count = 0;
    /* TODO: optimize: start in attributes array at position of the first tree,
     * resp. the first tree with attributes
     */
    for (iz = 0; iz < cmesh->stash->attributes.elem_count; iz++) {
      attribute =
        (t8_stash_attribute_struct_t *) sc_array_index (&cmesh->stash->
                                                        attributes, iz);
      if (cmesh->first_tree <= attribute->id && attribute->id <= last_tree) {
        /* TODO: check for duplicate attributes */
        attr_byte_count += attribute->attr_size;
      }
    }
    /********************************************************/
    /*      END COUNTING TREES/GHOSTS/ATTRIBUTES            */
    /********************************************************/
    /* Now that we know the number of trees/ghosts/and attribute_bytes we can
     * initialize the trees structure. */
    /* This even initializes the trees structure if there are neither trees
     * nor ghosts */
    t8_debugf ("Init trees with %li T, %li G\n",
               (long) cmesh->num_local_trees, (long) cmesh->num_ghosts);
    t8_cmesh_trees_init (&cmesh->trees, 1, cmesh->num_local_trees,
                         cmesh->num_ghosts);
    t8_cmesh_trees_init_part (cmesh->trees, 0, 0, cmesh->num_local_trees,
                              cmesh->num_ghosts, attr_byte_count);
    if (cmesh->num_local_trees > 0) {
      iz = (size_t) t8_stash_class_bsearch (cmesh->stash, cmesh->first_tree);
      class_end = (size_t) t8_stash_class_bsearch (cmesh->stash, last_tree);
      T8_ASSERT (t8_stash_class_bsearch (cmesh->stash, cmesh->first_tree) >=
                 0);
      T8_ASSERT (t8_stash_class_bsearch (cmesh->stash, last_tree) >= 0);
      T8_ASSERT ((t8_locidx_t) class_end - (t8_locidx_t) iz + 1 ==
                 cmesh->num_local_trees);
    }
    else {
      t8_debugf ("Empty partition.\n");
      iz = 0;
      class_end = -1;
    }
    /* TODO: optimize if non-hybrid mesh */
    /* loop over all local trees and set classes */
    for (; iz < (size_t) class_end + 1; iz++) {
      /* get class and tree id */
      classentry = (t8_stash_class_struct_t *)
        sc_array_index (&cmesh->stash->classes, iz);
      /* initialize tree */
      t8_cmesh_trees_add_tree (cmesh->trees,
                               classentry->id - cmesh->first_tree, 0,
                               classentry->eclass);
      cmesh->num_trees_per_eclass[classentry->eclass]++;
    }
    /* TODO: optimize if non-hybrid mesh */
    /* Iterate through ghosts and set classes */
    id1 = -1;
    for (iz = 0; iz < cmesh->num_ghosts; iz++) {
        ghost_facejoin = (struct ghost_facejoins_struct *)
          sc_array_index (ghost_ids, ghost_offsets[iz]);
      T8_ASSERT (ghost_facejoin->flag & 1);
      T8_ASSERT (ghost_facejoin->ghost_id > id1);
      id1 = ghost_facejoin->ghost_id;
      /* Get position of ghost in classes array */
      /* TODO: optimize so that we do not need this bsearch */
      ghost_ind = t8_stash_class_bsearch (cmesh->stash, id1);
      T8_ASSERT (ghost_ind >= 0);
      classentry = (t8_stash_class_struct_t *)
          sc_array_index_ssize_t (&cmesh->stash->classes, ghost_ind);
      t8_cmesh_trees_add_ghost (cmesh->trees, iz, id1, 0,
                                classentry->eclass);
      cmesh->trees->ghost_to_offset[iz] = iz;
    }
    /* We are done with stash->classes now  so we free memory.
     * Since the array is destroyed in stash_destroy we only reset it. */
    sc_array_reset (&cmesh->stash->classes);
    /* Go through all face_neighbour entries and parse every
     * local tree to local tree entry */
    /* TODO: could optimize by also using an questions array */
    for (iz = 0; iz < cmesh->stash->joinfaces.elem_count; iz++) {
      joinface = (t8_stash_joinface_struct_t *)
        sc_array_index (&cmesh->stash->joinfaces, iz);
      id1_istree = cmesh->first_tree <= joinface->id1 &&
        last_tree >= joinface->id1;
      id2_istree = cmesh->first_tree <= joinface->id2 &&
        last_tree >= joinface->id2;
      /* Both trees in the connection are local trees */
      if (id1_istree && id2_istree) {
        t8_cmesh_tree_set_join (cmesh->trees,
                                joinface->id1 - cmesh->first_tree,
                                joinface->id2 - cmesh->first_tree,
                                joinface->face1, joinface->face2,
                                joinface->orientation);
      }
    }
    /* Go through ghost_id array and set all facejoins with at least one ghost involved */
    id1 = -1;
    is_true_ghost = 0;
    for (iz = 0, jz = -1; iz < cmesh->num_ghosts; iz++) {
      temp = 0;
      ghost_facejoin = (struct ghost_facejoins_struct *)
        sc_array_index (ghost_ids, ghost_offsets[iz]);
      T8_ASSERT (ghost_facejoin->ghost_id > id1);
      T8_ASSERT (ghost_facejoin->flag & 1);
      ghost1 = t8_cmesh_trees_get_ghost (cmesh->trees, iz);
      joinface = (t8_stash_joinface_struct_t *)
          sc_array_index (&cmesh->stash->joinfaces, ghost_facejoin->index);
      /* Store the id of the neighbor */
      id2 = ghost_facejoin->flag & 2 ? joinface->id1 : joinface->id2;
      /* Store the face of the ghost */
      face = ghost_facejoin->flag & 2 ? joinface->face2 : joinface->face1;
      ghost1->neighbors[face] = id2;
      /* Check whether the second tree in face-connection is a local tree */
      /* If it is we have to add the connection to the local tree */
      id2_istree = cmesh->first_tree <= id2 && last_tree >= id2;
      if (id2_istree) {
        tree1 = t8_cmesh_trees_get_tree (cmesh->trees, id2 - cmesh->first_tree);
        /* Store the face of the tree */
        face2 =
          ghost_facejoin->flag & 2 ? joinface->face1 : joinface->face2;
        /* Set the local ghost id as neighbor */
        tree1->face_neighbors[face2] = jz + cmesh->first_tree;
        F = t8_eclass_num_faces[ghost1->eclass];
        tree1->tree_to_face[face2] = face * F + joinface->orientation;
      }
    }
    T8_FREE (ghost_offsets);
    sc_array_destroy (ghost_ids);
    /* Add attributes to the local trees */
    t8_cmesh_add_attributes (cmesh, cmesh->stash);
    /* compute global number of trees. id1 serves as buffer since
     * global number and local number have different datatypes */
    id1 = cmesh->num_local_trees;
    sc_MPI_Allreduce (&id1, &cmesh->num_trees, 1, T8_MPI_GLOIDX,
                      sc_MPI_SUM, cmesh->mpicomm);
  }

  cmesh->committed = 1;

  /* dup communicator if requested */
  if (cmesh->do_dup) {
    mpiret = sc_MPI_Comm_dup (cmesh->mpicomm, &comm_dup);
    SC_CHECK_MPI (mpiret);
    cmesh->mpicomm = comm_dup;
  }
  /* query communicator new */
  mpiret = sc_MPI_Comm_size (cmesh->mpicomm, &cmesh->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (cmesh->mpicomm, &cmesh->mpirank);
  SC_CHECK_MPI (mpiret);
  t8_stash_destroy (&cmesh->stash);
}
