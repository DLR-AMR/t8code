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
#include <sc_flops.h>
#include <sc_statistics.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh_trees.h"

typedef struct ghost_facejoins_struct
{
  t8_gloidx_t         ghost_id; /* The id of the ghost */
  t8_locidx_t         local_id; /* The local id of the ghost */
} t8_ghost_facejoin_t;

static int
t8_ghost_facejoins_compare (const void *fj1, const void *fj2)
{
  struct ghost_facejoins_struct *FJ1 = (struct ghost_facejoins_struct *) fj1;
  struct ghost_facejoins_struct *FJ2 = (struct ghost_facejoins_struct *) fj2;

  return FJ1->ghost_id < FJ2->ghost_id ? -1 : FJ1->ghost_id != FJ2->ghost_id;
}

static int
t8_ghost_facejoin_equal (const void *v1, const void *v2, const void *u)
{
  return t8_ghost_facejoins_compare (v1, v2) == 0;
}

/* The hash value for a ghost is
 * global_id % num_local_trees
 *
 * This hash function gets as input the global ghost id
 * and as user data the number of local trees.
 */
static unsigned
t8_ghost_hash (const void *v, const void *u)
{
  t8_gloidx_t         ghost_id = *((t8_gloidx_t *) v);
  t8_locidx_t         num_local_trees = *((t8_locidx_t *) u);

  return ghost_id % num_local_trees;
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


  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[3];

  sc_flops_start (&fi);


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
    sc_hash_t          *ghost_ids;
    sc_mempool_t       *ghost_facejoin_mempool;
    struct ghost_facejoins_struct *ghost_facejoin, *temp_facejoin,
      **facejoin_pp;
    size_t              joinfaces_it, attr_byte_count, iz;
    t8_stash_joinface_struct_t *joinface;
    t8_gloidx_t         last_tree = cmesh->num_local_trees +
      cmesh->first_tree - 1, id1, id2;
    t8_locidx_t         temp_local_id;
    t8_stash_attribute_struct_t *attribute;
    t8_stash_class_struct_t *classentry;
    int                 id1_istree, id2_istree, F;
    t8_ctree_t          tree1, tree2;
    t8_cghost_t         ghost1, ghost2;


    sc_flops_snap (&fi, &snapshot);

    if (cmesh->face_knowledge != 3) {
      t8_global_errorf ("Expected a face knowledge of 3.\nAbort commit.");
      /* TODO: reset cmesh */
      return;
    }
    t8_cmesh_set_shmem_type (cmesh);    /* TODO: do we actually need the shared array? */
    t8_stash_class_sort (cmesh->stash); /* TODO: do we need this to be sorted */
//    t8_stash_joinface_sort (cmesh->stash); /* TODO: this is propably not usefull */
    t8_stash_attribute_sort (cmesh->stash);

    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[0], snapshot.iwtime, "cmesh_commit_sort");
    sc_flops_snap (&fi, &snapshot);

    ghost_facejoin_mempool = sc_mempool_new (sizeof (t8_ghost_facejoin_t));
    ghost_ids = sc_hash_new (t8_ghost_hash, t8_ghost_facejoin_equal,
                             &cmesh->num_local_trees, ghost_facejoin_mempool);
    cmesh->num_ghosts = 0;
    /* Parse joinfaces array and save all global id of local ghosts, and assign them a local id */
    for (joinfaces_it = 0; joinfaces_it < cmesh->stash->joinfaces.elem_count;
         joinfaces_it++) {
      joinface =
        (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->
                                                       stash->joinfaces,
                                                       joinfaces_it);
      id1 = joinface->id1;
      id2 = joinface->id2;
      id2_istree = id2 <= last_tree && id2 >= cmesh->first_tree;
      id1_istree = id1 <= last_tree && id1 >= cmesh->first_tree;
      if (id1_istree || id2_istree) {
        /* Only consider facejoins with local trees involved to geth the global
         * ids of all local ghosts. */
        if (!id2_istree) {
          /* id2 is a ghost */
          if (sc_hash_insert_unique (ghost_ids, &id2,
                                     (void ***) &facejoin_pp)) {
            /* If we did not already stored id2 in the hash we do so and assign the next local ghost id */
            ghost_facejoin = *facejoin_pp;
            ghost_facejoin->ghost_id = id2;
            ghost_facejoin->local_id = cmesh->num_ghosts++;
          }
        }
        if (!id1_istree) {
          /* id1 is a ghost */
          T8_ASSERT (id2_istree);
          if (sc_hash_insert_unique (ghost_ids, &id2,
                                     (void ***) &facejoin_pp)) {
            /* If we did not already stored id1 in the hash we do so and assign the next local ghost id */
            ghost_facejoin = *facejoin_pp;
            ghost_facejoin->ghost_id = id1;
            ghost_facejoin->local_id = cmesh->num_ghosts++;
          }
        }
      }
    }

    /* Count attribute bytes */
    attr_byte_count = 0;
    /* TODO: optimize: start in attributes array at position of the first tree,
     * resp. the first tree with attributes
     */
    for (iz = 0; iz < cmesh->stash->attributes.elem_count; iz++) {
      attribute =
        (t8_stash_attribute_struct_t *) sc_array_index (&cmesh->
                                                        stash->attributes,
                                                        iz);
      if (cmesh->first_tree <= attribute->id && attribute->id <= last_tree) {
        /* TODO: check for duplicate attributes */
        attr_byte_count += attribute->attr_size;
      }
    }

    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[1], snapshot.iwtime, "cmesh_commit_count");
    sc_flops_snap (&fi, &snapshot);
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

#ifdef T8_ENABLE_DEBUG
    if (cmesh->num_local_trees == 0) {
      t8_debugf ("Empty partition.\n");
    }
#endif
    /* TODO: optimize if non-hybrid mesh */
    /* Itterate through classes and add ghosts and trees */
    /* We need a temporary ghost_facejoin to check the hash for existing global ids */
    temp_facejoin = T8_ALLOC (struct ghost_facejoins_struct, 1);
    temp_facejoin->local_id = -1;
    for (iz = 0; iz < cmesh->stash->classes.elem_count; iz++) {
      /* get class and tree id */
      classentry = (t8_stash_class_struct_t *)
        sc_array_index (&cmesh->stash->classes, iz);
      if (cmesh->first_tree <= classentry->id && classentry->id <= last_tree) {
        /* initialize tree */
        t8_cmesh_trees_add_tree (cmesh->trees,
                                 classentry->id - cmesh->first_tree, 0,
                                 classentry->eclass);
        cmesh->num_trees_per_eclass[classentry->eclass]++;
      }
      else {
        ghost_facejoin->ghost_id = classentry->id;
        if (sc_hash_lookup (ghost_ids, temp_facejoin,
                            (void ***) &facejoin_pp)) {
          /* The classentry belongs to a local ghost */
          ghost_facejoin = *facejoin_pp;
          t8_cmesh_trees_add_ghost (cmesh->trees, ghost_facejoin->local_id,
                                    ghost_facejoin->ghost_id, 0,
                                    classentry->eclass);
        }
      }
    }
    /* We are done with stash->classes now  so we free memory.
     * Since the array is destroyed in stash_destroy we only reset it. */
    sc_array_reset (&cmesh->stash->classes);
    /* Go through all face_neighbour entries and parse every
     * important entry */
    for (iz = 0; iz < cmesh->stash->joinfaces.elem_count; iz++) {
      joinface = (t8_stash_joinface_struct_t *)
        sc_array_index (&cmesh->stash->joinfaces, iz);
      id1_istree = cmesh->first_tree <= joinface->id1 &&
        last_tree >= joinface->id1;
      id2_istree = cmesh->first_tree <= joinface->id2 &&
        last_tree >= joinface->id2;
      temp_facejoin->ghost_id = joinface->id1;
      id1 = joinface->id1;
      tree1 = NULL;
      ghost1 = NULL;
      /* There are the following cases:
       * Both trees are local trees.
       * One is a local tree and one a local ghost.
       * Both are local ghosts.
       * One is a local ghost and on neither ghost nor local tree.
       * For each of these cases we have to set the correct face connection.
       */
      /* TODO: This if else if else spaghetti code is quite ugly, think it over
       *       and find a better way */
      if (id1_istree) {
        /* First tree in the connection is a local tree */
        tree1 = t8_cmesh_trees_get_tree (cmesh->trees,
                                         joinface->id1 - cmesh->first_tree);
      }
      else if (sc_hash_lookup (ghost_ids, &temp_facejoin,
                               (void ***) &facejoin_pp)) {
        /* id1 is a local ghost */
        ghost_facejoin = *facejoin_pp;
        ghost1 = t8_cmesh_trees_get_ghost (cmesh->trees,
                                           ghost_facejoin->local_id);
        temp_local_id = ghost_facejoin->local_id;
      }
      id2 = temp_facejoin->ghost_id = joinface->id2;
      ghost2 = NULL;
      tree2 = NULL;
      if (id2_istree) {
        /* Second tree in the connection is a local tree */
        tree2 = t8_cmesh_trees_get_tree (cmesh->trees,
                                         joinface->id2 - cmesh->first_tree);
      }
      else if (sc_hash_lookup (ghost_ids, &temp_facejoin,
                               (void ***) &facejoin_pp)) {
        /* id1 is a local ghost */
        ghost_facejoin = *facejoin_pp;
        ghost2 = t8_cmesh_trees_get_ghost (cmesh->trees,
                                           ghost_facejoin->local_id);
      }
      if (tree1) {
        tree1->face_neighbors[joinface->face1] =
          tree2 ? id2 - cmesh->first_tree :
          ghost_facejoin->local_id + cmesh->num_local_trees;
        F =
          t8_eclass_num_faces[tree2 != NULL ? tree2->eclass : ghost2->eclass];
        tree1->tree_to_face[joinface->face1] =
          F * joinface->orientation + joinface->face2;
      }
      else if (ghost1) {
        ghost1->neighbors[joinface->face1] = id2;
      }
      else {
        T8_ASSERT (ghost2 != NULL);
        ghost2->neighbors[joinface->face2] = id1;
        /* Done with setting face join */
      }
      if (tree2) {
        tree2->face_neighbors[joinface->face2] =
          tree1 ? id1 - cmesh->first_tree :
          temp_local_id + cmesh->num_local_trees;
        F =
          t8_eclass_num_faces[tree1 != NULL ? tree1->eclass : ghost1->eclass];
        tree1->tree_to_face[joinface->face2] =
          F * joinface->orientation + joinface->face1;
      }
      else if (ghost2) {
        ghost2->neighbors[joinface->face2] = id2;
      }
      else {
        T8_ASSERT (ghost1 != NULL);
        /* Done with setting face join */
      }
    }
    sc_hash_destroy (ghost_ids);
    sc_mempool_destroy (ghost_facejoin_mempool);
    /* Add attributes to the local trees */
    t8_cmesh_add_attributes (cmesh, cmesh->stash);
    /* compute global number of trees. id1 serves as buffer since
     * global number and local number have different datatypes */
    id1 = cmesh->num_local_trees;
    sc_MPI_Allreduce (&id1, &cmesh->num_trees, 1, T8_MPI_GLOIDX,
                      sc_MPI_SUM, cmesh->mpicomm);

    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[2], snapshot.iwtime, "cmesh_commit_end");


    sc_stats_compute (sc_MPI_COMM_WORLD, 3, stats);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 3, stats, 1, 1);
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
