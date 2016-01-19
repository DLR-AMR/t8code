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
  int                 mpiret, key;
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
      t8_locidx_t         tree_id, temptree;
      size_t              si, attr_bytes, attr_offset, num_attr, attr_index;
      int                 attr_id;

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
      attr_offset = 0;
      temptree = -1;
      for (si = 0; si < stash->attributes.elem_count; si++) {
        tree_id = t8_stash_get_attribute_tree_id (cmesh->stash, si) -
          cmesh->first_tree;
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
          num_attr -= si;       /* now stores the number of attribute of the tree */
          /* initialize storage for tree attributes */
          t8_cmesh_trees_init_attributes (cmesh->trees, tree_id, num_attr);
          temptree = tree_id;   /* store the current tree_id in oldtree */
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
  else {
    sc_array_t         *ghost_ids;
    size_t              joinfaces_it, attr_byte_count, iz, jz, class_end;
    ssize_t             ghost_ind;
    t8_stash_joinface_struct_t *joinface;
    t8_gloidx_t         last_tree = cmesh->num_local_trees +
      cmesh->first_tree - 1, id1, id2, *ghost;
    t8_stash_attribute_struct_t *attribute;
    t8_stash_class_struct_t *classentry;
    int                 id1_istree, id2_istree, F;
    t8_ctree_t          tree1;
    t8_cghost_t         ghost1;

    if (cmesh->face_knowledge != 3) {
      t8_global_errorf ("Expected a face knowledge of 3.\nAbort commit.");
      /* TODO: reset cmesh */
      return;
    }
    t8_cmesh_set_shmem_type (cmesh);    /* TODO: do we actually need the shared array? */
    t8_stash_class_sort (cmesh->stash);
    t8_stash_joinface_sort (cmesh->stash);
    t8_stash_attribute_sort (cmesh->stash);

    ghost_ids = sc_array_new (sizeof (t8_gloidx_t));
    for (joinfaces_it = 0; joinfaces_it < cmesh->stash->joinfaces.elem_count;
         joinfaces_it++) {
      joinface =
        (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->stash->
                                                       joinfaces,
                                                       joinfaces_it);
      id1 = joinface->id1;
      id2 = joinface->id2;
      if (cmesh->first_tree <= id1 && last_tree >= id1) {
        if (id2 > last_tree || id2 < cmesh->first_tree) {
          /* id2 is a ghost */
          ghost = (t8_gloidx_t *) sc_array_push (ghost_ids);
          *ghost = id2;
        }
      }
      else if (cmesh->first_tree <= id2 && last_tree >= id2) {
        /* id1 is a ghost */
        ghost = (t8_gloidx_t *) sc_array_push (ghost_ids);
        *ghost = id1;
      }
    }

    /* Added one ghostid entry for each face connection with a ghost.
     * Need to remove duplicates. */
    sc_array_sort (ghost_ids, t8_compare_gloidx);
    /* Count the ghosts without duplicates */
    /* Start with 1 if the ghost array is nonempty */
    /* Remove duplicate entries in ghost_id array */
    if (ghost_ids->elem_count > 0) {
      id1 = *((t8_gloidx_t *) sc_array_index (ghost_ids, 0));
      id2 = id1;
      jz = 0;
      iz = 0;
      while (iz < ghost_ids->elem_count) {
        for (; iz < ghost_ids->elem_count && id2 == id1; iz++) {
          id2 = *((t8_gloidx_t *) sc_array_index (ghost_ids, iz));
        }
        if (id2 != id1) {
          T8_ASSERT (id2 > id1);
          jz++;
          id1 = *((t8_gloidx_t *) sc_array_index (ghost_ids, jz)) = id2;
        }
      }
      sc_array_resize (ghost_ids, jz + 1);
    }
    cmesh->num_ghosts = ghost_ids->elem_count;

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
    /* Now that we know the number of trees/ghosts/and attribute_bytes we can
     * initialize the trees structure. */
    /* This even initializes the trees structure if there are neither trees
     * nor ghosts */
    t8_debugf ("Init trees with %li T, %li G\n", (long) cmesh->num_local_trees,
               (long) cmesh->num_ghosts);
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
      t8_debugf ("%i %i\n", (int) iz, (int) class_end);
      T8_ASSERT ((t8_locidx_t) class_end - (t8_locidx_t) iz + 1 ==
                 cmesh->num_local_trees);
    }
    else {
      t8_debugf ("Empty partition.\n");
      iz = 0;
      class_end = -1;
    }
    /* TODO: optimize if non-hybrid mesh */
    /* loop over all local trees */
    for (; iz < class_end + 1; iz++) {
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
    /* Iterate through ghosts */
    for (iz = 0; iz < ghost_ids->elem_count; iz++) {
      ghost = (t8_gloidx_t *) sc_array_index (ghost_ids, iz);
      /* Get position of ghost in classes array */
      ghost_ind = t8_stash_class_bsearch (cmesh->stash, *ghost);
      T8_ASSERT (ghost_ind >= 0);
      classentry = (t8_stash_class_struct_t *)
        sc_array_index_ssize_t (&cmesh->stash->classes, ghost_ind);
      t8_cmesh_trees_add_ghost (cmesh->trees, iz, *ghost, 0,
                                classentry->eclass);
      cmesh->trees->ghost_to_offset[iz] = iz;
    }
    /* We are done with stash->classes now  so we free memory.
     * Since the array is destroyed in stash_destroy we only reset it. */
    sc_array_reset (&cmesh->stash->classes);
    /* Go through all face_neighbour entries and parse every
     * entry that has a local tree or ghost involved */
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
      /* One tree is a local tree the other one a ghost */
      else if (id1_istree || id2_istree) {
        /* If second one is local tree we swap first and second one */
        if (id2_istree) {
          id1 = joinface->id1;
          joinface->id1 = joinface->id2;
          joinface->id2 = id1;
          F = joinface->face1;
          joinface->face1 = joinface->face2;
          joinface->face2 = F;
        }
        /* From here on the first tree is a local tree and the second one a ghost */
        ghost_ind = sc_array_bsearch (ghost_ids, &joinface->id2,
                                      t8_compare_gloidx);
        /* Search for local ghost id in array */
        t8_debugf
          ("Face-conn %lli <--> %lli via faces %i and %i. Ghost_index = %li\n",
           (long long) joinface->id1, (long long) joinface->id2,
           joinface->face1, joinface->face2, (long) ghost_ind);
        SC_CHECK_ABORTF (ghost_ind >= 0
                         && ghost_ind < cmesh->num_ghosts,
                         "Ghost tree (global id %lli) in face-connections that"
                         " is not in classes.\n", (long long) joinface->id2);
        tree1 =
          t8_cmesh_trees_get_tree (cmesh->trees,
                                   joinface->id1 - cmesh->first_tree);
        t8_debugf ("Set face %i of tree %i to G %zd\n", joinface->face1,
                   tree1->treeid, ghost_ind);
        ghost1 =
          t8_cmesh_trees_get_ghost (cmesh->trees, (t8_locidx_t) ghost_ind);
        t8_debugf ("Hello\n");
        tree1->face_neighbors[joinface->face1] = ghost_ind;
        F = t8_eclass_num_faces[ghost1->eclass];
        tree1->tree_to_face[joinface->face1] = joinface->face2 * F +
          joinface->orientation;
        ghost1->neighbors[joinface->face2] = joinface->id1;
      }
      else {
        /* Either both trees are ghosts or non is local tree or ghost */
        ghost1 = NULL;
        ghost_ind = sc_array_bsearch (ghost_ids, &joinface->id1,
                                      t8_compare_gloidx);
        if (ghost_ind >= 0) {
          /* First tree is local ghost */
          ghost1 = t8_cmesh_trees_get_ghost (cmesh->trees,
                                             (t8_locidx_t) ghost_ind);
          ghost1->neighbors[joinface->face1] = joinface->id2;
        }
        ghost_ind = sc_array_bsearch (ghost_ids, &joinface->id2,
                                      t8_compare_gloidx);
        if (ghost_ind >= 0) {
          /* Second tree is (maybe also) local ghost */
          ghost1 = t8_cmesh_trees_get_ghost (cmesh->trees,
                                             (t8_locidx_t) ghost_ind);
          ghost1->neighbors[joinface->face2] = joinface->id1;
        }
      }
    }
    sc_array_destroy (ghost_ids);

    //SC_ABORTF ("partitioned commit not implemented.%c", '\n');
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
