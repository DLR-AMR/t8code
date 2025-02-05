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

/** \file t8_cmesh_commit.cxx
 *
 * TODO: document this file
 */

#include <t8_data/t8_shmem.h>
#include <t8_cmesh.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_copy.h>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_handler.hxx>

typedef struct ghost_facejoins_struct
{
  t8_gloidx_t ghost_id; /* The id of the ghost */
  t8_locidx_t local_id; /* The local id of the ghost */
  t8_gloidx_t attr_id;  /* The current number of inserted ghost attributes */
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
 * This hash function gets as input a facejoins_struct
 * and as user data the number of hashes, which is the number of local trees if
 * nonzero and 10 otherwise.
 */
static unsigned
t8_ghost_hash (const void *v, const void *u)
{
  t8_gloidx_t ghost_id = ((t8_ghost_facejoin_t *) v)->ghost_id;
  t8_locidx_t num_hashes = *((t8_locidx_t *) u);

  return ghost_id % num_hashes;
}

static void
t8_cmesh_set_shmem_type (sc_MPI_Comm comm)
{
  T8_ASSERT (comm != sc_MPI_COMM_NULL);

  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
}

static void
t8_cmesh_add_attributes (const t8_cmesh_t cmesh, sc_hash_t *ghost_ids)
{
  t8_stash_attribute_struct_t *attribute;
  const t8_stash_t stash = cmesh->stash;
  t8_locidx_t ltree;
  size_t si, sj;
  t8_ghost_facejoin_t **facejoin_pp; /* used to lookup global ghost ids in the hash */

  t8_ghost_facejoin_t temp_facejoin = { 0, 0, 0 };

  t8_locidx_t ghosts_inserted = 0;
  ltree = -1;
  for (si = 0, sj = 0; si < stash->attributes.elem_count; si++, sj++) {
    attribute = (t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, si);
    if (attribute->id > ltree) {
      /* Enter a new tree */
      ltree = attribute->id;
      sj = 0;
    }
    if (cmesh->first_tree <= attribute->id && attribute->id < cmesh->first_tree + cmesh->num_local_trees) {
      /* attribute->id is a gloidx that is casted to a locidx here.
       * Should not cause problems, since mesh is replicated */
      T8_ASSERT (attribute->id - cmesh->first_tree == (t8_locidx_t) attribute->id - cmesh->first_tree);
      t8_cmesh_trees_add_attribute (cmesh->trees, 0, attribute, attribute->id - cmesh->first_tree, sj);
    }
    else {
      T8_ASSERT (ghost_ids != NULL);
      temp_facejoin.ghost_id = attribute->id;
      if (sc_hash_lookup (ghost_ids, &temp_facejoin, (void ***) &facejoin_pp)) {
        T8_ASSERT ((t8_locidx_t) sj == (t8_locidx_t) (*facejoin_pp)->attr_id);
        if (sj == 0) {
          ghosts_inserted++;
        }
        /* attribute is on a ghost tree */
        t8_cmesh_trees_add_ghost_attribute (cmesh->trees, attribute, (*facejoin_pp)->local_id, ghosts_inserted,
                                            (*facejoin_pp)->attr_id);
        (*facejoin_pp)->attr_id++;
      }
    }
  }
}

static void
t8_cmesh_commit_replicated_new (t8_cmesh_t cmesh)
{
  t8_stash_attribute_struct_t *attribute;
  t8_locidx_t *face_neigh, *face_neigh2;
  int8_t *ttf, *ttf2;
  t8_stash_joinface_struct_t *joinface;
  t8_ctree_t tree1;
  size_t si;

  if (cmesh->stash != NULL && cmesh->stash->classes.elem_count > 0) {
    t8_stash_t stash = cmesh->stash;
    sc_array_t *class_entries = &stash->classes;
    t8_stash_class_struct_t *entry;
    t8_locidx_t num_trees = class_entries->elem_count, ltree;

    t8_cmesh_trees_init (&cmesh->trees, 1, num_trees, 0);
    t8_cmesh_trees_start_part (cmesh->trees, 0, 0, num_trees, 0, 0, 1);
    /* set tree classes */
    for (ltree = 0; ltree < num_trees; ltree++) {
      entry = (t8_stash_class_struct_t *) t8_sc_array_index_locidx (class_entries, ltree);
      t8_cmesh_trees_add_tree (cmesh->trees, entry->id, 0, entry->eclass);
      cmesh->num_trees_per_eclass[entry->eclass]++;
      cmesh->num_local_trees_per_eclass[entry->eclass]++;
    }
    for (si = 0; si < stash->attributes.elem_count; si++) {
      attribute = (t8_stash_attribute_struct_t *) sc_array_index (&stash->attributes, si);
      /* attribute->id is a gloidx that is casted to a locidx here.
       * Should not cause problems, since mesh is replicated */
      T8_ASSERT (attribute->id == (t8_locidx_t) attribute->id);
      tree1 = t8_cmesh_trees_get_tree (cmesh->trees, attribute->id);
      tree1->num_attributes++;
      tree1->att_offset += attribute->attr_size;
    }
    /* Finish memory allocation of tree/ghost/face/attribute array
     * using the info calculated above */
    t8_cmesh_trees_finish_part (cmesh->trees, 0);
    cmesh->num_trees = cmesh->num_local_trees = num_trees;

    /* Add attributes */
    /* TODO: currently the attributes array still has to be sorted,
     *       find a way around it */
    t8_stash_attribute_sort (cmesh->stash);
    cmesh->num_trees = cmesh->num_local_trees = num_trees;
    cmesh->first_tree = 0;
    t8_cmesh_add_attributes (cmesh, NULL);

    /* Set all face connections */
    t8_cmesh_trees_set_all_boundary (cmesh, cmesh->trees);
    for (si = 0; si < cmesh->stash->joinfaces.elem_count; si++) {
      joinface = (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->stash->joinfaces, si);
      (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, joinface->id1, &face_neigh, &ttf);
      (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, joinface->id2, &face_neigh2, &ttf2);
      /* Set the correct face neighbors and tree_to_face values */
      face_neigh[joinface->face1] = (t8_locidx_t) joinface->id2;
      ttf[joinface->face1]
        = t8_cmesh_tree_to_face_encode (cmesh->dimension, (t8_locidx_t) joinface->face2, joinface->orientation);
      face_neigh2[joinface->face2] = (t8_locidx_t) joinface->id1;
      ttf2[joinface->face2]
        = t8_cmesh_tree_to_face_encode (cmesh->dimension, (t8_locidx_t) joinface->face1, joinface->orientation);
    }
  }
}

static void
t8_cmesh_commit_partitioned_new (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  /* Cmesh is partitioned and new */
  t8_stash_attribute_struct_t *attribute;
  t8_locidx_t *face_neigh, *face_neigh2;
  int8_t *ttf, *ttf2;
  t8_stash_joinface_struct_t *joinface;
  t8_ctree_t tree1, tree2;
  t8_cghost_t ghost1;
  int F;
  size_t si;

#if T8_ENABLE_DEBUG
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[3];
#endif

  sc_hash_t *ghost_ids;
  sc_mempool_t *ghost_facejoin_mempool;
  t8_ghost_facejoin_t *ghost_facejoin = NULL, *temp_facejoin, **facejoin_pp;
  size_t joinfaces_it, iz;
  t8_gloidx_t id1, id2;
  t8_locidx_t temp_local_id = 0;
  t8_locidx_t num_hashes;
  t8_gloidx_t *face_neigh_g, *face_neigh_g2;
  t8_stash_class_struct_t *classentry;
  int id1_istree, id2_istree;
#if T8_ENABLE_DEBUG
  t8_cghost_t ghost2;
#endif

#if T8_ENABLE_DEBUG
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
#endif

  T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));

  if (cmesh->face_knowledge != 3) {
    t8_global_errorf ("Expected a face knowledge of 3.\nAbort commit.");
    /* TODO: reset cmesh */
    return;
  }
  t8_shmem_init (comm);
  t8_cmesh_set_shmem_type (comm); /* TODO: do we actually need the shared array? */
  t8_stash_attribute_sort (cmesh->stash);

#if T8_ENABLE_DEBUG
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "cmesh_commit_sort");
  sc_flops_snap (&fi, &snapshot);
#endif

  if (cmesh->tree_offsets != NULL) {
    const t8_gloidx_t *tree_offsets = (t8_gloidx_t *) t8_shmem_array_get_gloidx_array (cmesh->tree_offsets);
    /* We partition using tree offsets */
    /* Get the first tree and whether it is shared */
    cmesh->first_tree = t8_offset_first (cmesh->mpirank, tree_offsets);
    cmesh->first_tree_shared = t8_shmem_array_get_gloidx (cmesh->tree_offsets, cmesh->mpirank) < 0;
    /* Get the number of local trees */
    cmesh->num_local_trees = t8_offset_num_trees (cmesh->mpirank, tree_offsets);
  }
  else {
    SC_CHECK_ABORT (cmesh->set_partition_level < 0,
                    "Do not use t8_cmesh_set_partition_uniform when creating a cmesh from stash!\n");
  }
  /* The first_tree and first_tree_shared entries must be set by now */
  T8_ASSERT (cmesh->first_tree >= 0);
  T8_ASSERT (cmesh->first_tree_shared >= 0);

  /* Compute the global id of the cmeshes last_tree.
   * This must happen after cmesh->first_tree is computed. */
  const t8_gloidx_t last_tree = cmesh->num_local_trees + cmesh->first_tree - 1;

  num_hashes = cmesh->num_local_trees > 0 ? cmesh->num_local_trees : 10;
  ghost_facejoin_mempool = sc_mempool_new (sizeof (t8_ghost_facejoin_t));
  ghost_ids = sc_hash_new (t8_ghost_hash, t8_ghost_facejoin_equal, &num_hashes, NULL);

  temp_facejoin = (t8_ghost_facejoin_t *) sc_mempool_alloc (ghost_facejoin_mempool);
  temp_facejoin->attr_id = 0;

  cmesh->num_ghosts = 0;
  /* Parse joinfaces array and save all global id of local ghosts, and assign them a local id */
  for (joinfaces_it = 0; joinfaces_it < cmesh->stash->joinfaces.elem_count; joinfaces_it++) {
    joinface = (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->stash->joinfaces, joinfaces_it);
    id1 = joinface->id1;
    id2 = joinface->id2;
    id2_istree = id2 <= last_tree && id2 >= cmesh->first_tree;
    id1_istree = id1 <= last_tree && id1 >= cmesh->first_tree;
    if (id1_istree || id2_istree) {
      /* Only consider facejoins with local trees involved to get the global
       * ids of all local ghosts. */
      if (!id2_istree) {
        /* id2 is a ghost */
        temp_facejoin->ghost_id = id2;
        if (sc_hash_insert_unique (ghost_ids, temp_facejoin, NULL)) {
          /* If we did not already stored id2 in the hash we do so and assign the next local ghost id */
          temp_facejoin->local_id = cmesh->num_ghosts++;
          temp_facejoin = (t8_ghost_facejoin_t *) sc_mempool_alloc (ghost_facejoin_mempool);
          temp_facejoin->attr_id = 0;
        }
      }
      if (!id1_istree) {
        /* id1 is a ghost */
        T8_ASSERT (id2_istree);
        temp_facejoin->ghost_id = id1;
        if (sc_hash_insert_unique (ghost_ids, temp_facejoin, NULL)) {
          /* If we did not already stored id1 in the hash we do so and assign the next local ghost id */
          temp_facejoin->local_id = cmesh->num_ghosts++;
          temp_facejoin = (t8_ghost_facejoin_t *) sc_mempool_alloc (ghost_facejoin_mempool);
          temp_facejoin->attr_id = 0;
        }
      }
    }
  }

#if T8_ENABLE_DEBUG
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[1], snapshot.iwtime, "cmesh_commit_count");
  sc_flops_snap (&fi, &snapshot);
#endif
  /********************************************************/
  /*             END COUNTING TREES/GHOSTS                */
  /********************************************************/
  /* Now that we know the number of trees and ghosts we can
   * start initializing the trees structure. */
  /* This even initializes the trees structure if there are neither trees
   * nor ghosts */
  t8_debugf ("Init trees with %li T, %li G\n", (long) cmesh->num_local_trees, (long) cmesh->num_ghosts);
  t8_cmesh_trees_init (&cmesh->trees, 1, cmesh->num_local_trees, cmesh->num_ghosts);
  t8_cmesh_trees_start_part (cmesh->trees, 0, 0, cmesh->num_local_trees, 0, cmesh->num_ghosts, 1);

#ifdef T8_ENABLE_DEBUG
  if (cmesh->num_local_trees == 0) {
    t8_debugf ("Empty partition.\n");
  }
#endif
  if (cmesh->num_local_trees != 0 || cmesh->num_ghosts != 0) {
    /* Only do something if the partition is not empty */
    /* TODO: optimize if non-hybrid mesh */
    /* Iterate through classes and add ghosts and trees */
    /* We need a temporary ghost_facejoin to check the hash for existing global ids */
    temp_facejoin->local_id = -10;
    for (iz = 0; iz < cmesh->stash->classes.elem_count; iz++) {
      /* get class and tree id */
      classentry = (t8_stash_class_struct_t *) sc_array_index (&cmesh->stash->classes, iz);
      temp_facejoin->ghost_id = classentry->id;
      if (cmesh->first_tree <= classentry->id && classentry->id <= last_tree) {
        /* initialize tree */
        t8_cmesh_trees_add_tree (cmesh->trees, classentry->id - cmesh->first_tree, 0, classentry->eclass);
        cmesh->num_local_trees_per_eclass[classentry->eclass]++;
      }
      else {
        if (sc_hash_lookup (ghost_ids, temp_facejoin, (void ***) &facejoin_pp)) {
          /* The classentry belongs to a local ghost */
          ghost_facejoin = *facejoin_pp;
          t8_cmesh_trees_add_ghost (cmesh->trees, ghost_facejoin->local_id, ghost_facejoin->ghost_id, 0,
                                    classentry->eclass, cmesh->num_local_trees);
        }
      }
    }
    /* We are done with stash->classes now  so we free memory.
     * Since the array is destroyed in stash_destroy we only reset it. */
    sc_array_reset (&cmesh->stash->classes);

    /* Parse through attributes to count the number of attributes per tree
     * and total size of attributes per tree */
    for (si = 0; si < cmesh->stash->attributes.elem_count; si++) {
      attribute = (t8_stash_attribute_struct_t *) sc_array_index (&cmesh->stash->attributes, si);
      temp_facejoin->ghost_id = attribute->id;
      if (cmesh->first_tree <= attribute->id && attribute->id < cmesh->first_tree + cmesh->num_local_trees) {
        /* attribute->id is a gloidx that is casted to a locidx here.
         * Should not cause problems, since mesh is replicated */
        T8_ASSERT (attribute->id - cmesh->first_tree == (t8_locidx_t) (attribute->id - cmesh->first_tree));
        tree1 = t8_cmesh_trees_get_tree (cmesh->trees, attribute->id - cmesh->first_tree);
        tree1->num_attributes++;
        /* temporarily misuse the att_offset to store the total attribute size, until t8_cmesh_trees_finish_part */
        tree1->att_offset += attribute->attr_size;
      }
      else if (sc_hash_lookup (ghost_ids, temp_facejoin, (void ***) &facejoin_pp)) {
        /* attribute is on a ghost tree */
        ghost1 = t8_cmesh_trees_get_ghost (cmesh->trees, (*facejoin_pp)->local_id);
        ghost1->num_attributes++;
        /* temporarily misuse the att_offset to store the total attribute size, until t8_cmesh_trees_finish_part */
        ghost1->att_offset += attribute->attr_size;
      }
    }
    t8_cmesh_trees_finish_part (cmesh->trees, 0);
    t8_cmesh_trees_set_all_boundary (cmesh, cmesh->trees);

    /* Go through all face_neighbour entries and parse every
     * important entry */
    for (iz = 0; iz < cmesh->stash->joinfaces.elem_count; iz++) {
      joinface = (t8_stash_joinface_struct_t *) sc_array_index (&cmesh->stash->joinfaces, iz);
      id1 = joinface->id1;
      id2 = joinface->id2;
      id1_istree = cmesh->first_tree <= id1 && last_tree >= id1;
      id2_istree = cmesh->first_tree <= id2 && last_tree >= id2;
      temp_facejoin->ghost_id = id1;
      tree1 = NULL;
#if T8_ENABLE_DEBUG
      ghost1 = NULL;
#endif
      /* There are the following cases:
       * Both trees are local trees.
       * One is a local tree and one a local ghost.
       * Both are local ghosts.
       * One is a local ghost and on neither ghost nor local tree.
       * For each of these cases we have to set the correct face connection.
       */
      /* TODO: This if else if else spaghetti code is quite ugly, think it over
       *       and find a better way */
      face_neigh = face_neigh2 = NULL;
      ttf = ttf2 = NULL;
      if (id1_istree) {
        /* First tree in the connection is a local tree */
        tree1 = t8_cmesh_trees_get_tree_ext (cmesh->trees, joinface->id1 - cmesh->first_tree, &face_neigh, &ttf);
      }
      else if (sc_hash_lookup (ghost_ids, temp_facejoin, (void ***) &facejoin_pp)) {
        /* id1 is a local ghost */
        ghost_facejoin = *facejoin_pp;
#if T8_ENABLE_DEBUG
        ghost1 = t8_cmesh_trees_get_ghost_ext (cmesh->trees, ghost_facejoin->local_id, &face_neigh_g, &ttf);
#else
        (void) t8_cmesh_trees_get_ghost_ext (cmesh->trees, ghost_facejoin->local_id, &face_neigh_g, &ttf);
#endif
        temp_local_id = ghost_facejoin->local_id;
      }
#if T8_ENABLE_DEBUG
      ghost2 = NULL;
#endif
      tree2 = NULL;
      temp_facejoin->ghost_id = joinface->id2;
      if (id2_istree) {
        /* Second tree in the connection is a local tree */
        tree2 = t8_cmesh_trees_get_tree_ext (cmesh->trees, joinface->id2 - cmesh->first_tree, &face_neigh2, &ttf2);
      }
      else if (sc_hash_lookup (ghost_ids, temp_facejoin, (void ***) &facejoin_pp)) {
        /* id1 is a local ghost */
        ghost_facejoin = *facejoin_pp;
#if T8_ENABLE_DEBUG
        ghost2 = t8_cmesh_trees_get_ghost_ext (cmesh->trees, ghost_facejoin->local_id, &face_neigh_g2, &ttf2);
#else
        (void) t8_cmesh_trees_get_ghost_ext (cmesh->trees, ghost_facejoin->local_id, &face_neigh_g2, &ttf2);
#endif
      }
      F = t8_eclass_max_num_faces[cmesh->dimension];
      if (ttf != NULL) {
        /* The first entry is either a tree or ghost */
        ttf[joinface->face1] = F * joinface->orientation + joinface->face2;
        if (tree1 != NULL) {
          /* First entry is a tree */
          T8_ASSERT (ghost2 != NULL || tree2 != NULL);
          face_neigh[joinface->face1]
            = tree2 ? id2 - cmesh->first_tree : ghost_facejoin->local_id + cmesh->num_local_trees;
        }
        else {
          /* First entry is a ghost */
          T8_ASSERT (ghost1 != NULL);
          face_neigh_g[joinface->face1] = id2;
        }
      }
      if (ttf2 != NULL) {
        /* The second entry is either a tree or a ghost */
        ttf2[joinface->face2] = F * joinface->orientation + joinface->face1;
        if (tree2 != NULL) {
          /* The second entry is a tree */
          T8_ASSERT (tree1 != NULL || ghost1 != NULL);
          face_neigh2[joinface->face2] = tree1 ? id1 - cmesh->first_tree : temp_local_id + cmesh->num_local_trees;
        }
        else {
          /* The second entry is a ghost */
          T8_ASSERT (ghost2 != NULL);
          face_neigh_g2[joinface->face2] = id1;
        }
      }
      /* Done with setting face join */
    }

    /* Add attributes */
    /* TODO: attribute is required to be sorted right now. Think this over
     *       if we cant work around it, at least use the sorted array above when
     *       counting the attributes per tree. */

    t8_stash_attribute_sort (cmesh->stash);
    t8_cmesh_add_attributes (cmesh, ghost_ids);

    /* compute global number of trees. id1 serves as buffer since
     * global number and local number have different datatypes */

  } /* End if nonempty partition */

  sc_mempool_free (ghost_facejoin_mempool, temp_facejoin);
  sc_hash_destroy (ghost_ids);
  sc_mempool_destroy (ghost_facejoin_mempool);

  id1 = cmesh->num_local_trees;
  /* We must not count shared trees. Thus, we subtract one if
   * the first tree is shared. However, we exclude the case where
   * this cmesh has no trees but nevertheless the first tree is marked
   * as shared. */
  if (cmesh->first_tree_shared && id1 > 0) {
    id1--;
  }
  sc_MPI_Allreduce (&id1, &cmesh->num_trees, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);

#if T8_ENABLE_DEBUG
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[2], snapshot.iwtime, "cmesh_commit_end");
  sc_stats_compute (comm, 3, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 3, stats, 1, 1);
#endif
}

void
t8_cmesh_commit_from_stash (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  T8_ASSERT (cmesh != NULL);

  if (cmesh->set_partition) {
    /* partitioned commit */
    t8_cmesh_commit_partitioned_new (cmesh, comm);
  }
  else {
    /* replicated commit */
    t8_cmesh_commit_replicated_new (cmesh);
  }
}

/* TODO: set boundary face connections here.
 *       not trivial if replicated and not level 3 face_knowledg
 *       Edit: boundary face is default. If no face-connection is added then
 *             we assume a boundary face.
 * TODO: Implement a debug check for mesh consistency between processes.
 */
/* TODO: split this up into smaller functions */
void
t8_cmesh_commit (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int mpiret;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (comm != sc_MPI_COMM_NULL);
  T8_ASSERT (!cmesh->committed);
  SC_CHECK_ABORT (0 <= cmesh->dimension && cmesh->dimension <= T8_ECLASS_MAX_DIM,
                  "Dimension of the cmesh is not set properly.\n");

  /* If profiling is enabled, we measure the runtime of  commit. */
  if (cmesh->profile != NULL) {
    cmesh->profile->commit_runtime = sc_MPI_Wtime ();
  }
  /* Get mpisize and rank */
  mpiret = sc_MPI_Comm_size (comm, &cmesh->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &cmesh->mpirank);
  SC_CHECK_MPI (mpiret);
  if (cmesh->set_from != NULL) {
    cmesh->dimension = cmesh->set_from->dimension;
    cmesh->allow_negative_volumes = cmesh->set_from->allow_negative_volumes;
    if (cmesh->face_knowledge == -1) {
      /* Keep the face knowledge of the from cmesh, if -1 was specified */
      cmesh->face_knowledge = cmesh->set_from->face_knowledge;
    }

    /* If present use the set geometry handler, otherwise take
     * over the handler from set_from. */
    if (cmesh->geometry_handler == NULL && cmesh->set_from->geometry_handler != NULL) {
      cmesh->geometry_handler = cmesh->set_from->geometry_handler;
      cmesh->geometry_handler->ref ();
    }

    if (cmesh->set_partition) {
      /* The cmesh should be partitioned */
      t8_cmesh_partition (cmesh, comm);
    }
    else {
      t8_cmesh_copy (cmesh, cmesh->set_from, comm);
    }
  } /* End set_from != NULL */
  else {
    t8_cmesh_commit_from_stash (cmesh, comm);
  }
  cmesh->committed = 1;

  /* Compute trees_per_eclass */
  t8_cmesh_gather_trees_per_eclass (cmesh, comm);

  if (cmesh->set_partition) {
    t8_cmesh_gather_treecount (cmesh, comm);
  }
  T8_ASSERT (cmesh->set_partition || cmesh->tree_offsets == NULL);

#if T8_ENABLE_DEBUG
  t8_debugf ("Cmesh is %spartitioned.\n", cmesh->set_partition ? "" : "not ");
  if (cmesh->set_partition) {
    t8_cmesh_offset_print (cmesh, comm);
  }
  //t8_cmesh_trees_print (cmesh, cmesh->trees);
#endif

  if (cmesh->set_from != NULL) {
    /* Unref set_from and set it to NULL */
    t8_cmesh_unref (&cmesh->set_from);
    cmesh->set_from = NULL;
  }
  if (cmesh->stash != NULL) {
    t8_stash_destroy (&cmesh->stash);
  }

  t8_debugf ("committed cmesh with %li local and %lli global trees and"
             " %li ghosts.\n",
             (long) cmesh->num_local_trees, (long long) cmesh->num_trees, (long) cmesh->num_ghosts);

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  #if T8_ENABLE_DEBUG
  if (cmesh->allow_negative_volumes == 0){
    T8_ASSERT (t8_cmesh_validate_geometry (cmesh) == 0);
  }
  else{
    int return_code = t8_cmesh_validate_geometry (cmesh);
    if (return_code == 1){
      t8_debugf("WARNING: Negative volume detected. Ignoring since negative volumes were allowed by the user.\n");
      T8_ASSERT (return_code < 2);
    }
  }
  #endif
  /* If profiling is enabled, we measure the runtime of  commit. */
  if (cmesh->profile != NULL) {
    cmesh->profile->commit_runtime = sc_MPI_Wtime () - cmesh->profile->commit_runtime;
    /* We also measure the number of shared trees,
     * it is the average over all first_tree_shared*mpisize values. */
    cmesh->profile->first_tree_shared = cmesh->first_tree_shared * cmesh->mpisize;
  }
}
