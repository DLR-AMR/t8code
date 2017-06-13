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

#include <sc_statistics.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_refcount.h>
#include <t8_shmem.h>
#ifdef T8_WITH_METIS
#include <metis.h>
#endif
#include "t8_cmesh_trees.h"

/** \file t8_cmesh.c
 *
 * TODO: document this file
 */

int
t8_cmesh_is_initialized (t8_cmesh_t cmesh)
{
  if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) &&
        !cmesh->committed)) {
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* TODO: check conditions that must always hold after init and before commit */
  if (0) {
    return 0;
  }
#endif

  return 1;
}

int
t8_cmesh_is_committed (t8_cmesh_t cmesh)
{
  if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) &&
        cmesh->committed)) {
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* TODO: check more conditions that must always hold after commit */
  if ((!t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees)) || 0) {
    return 0;
  }
#endif

  return 1;
}

#if 0
/* Compute a hash value for a ghost tree. */
/* deprecated */
static unsigned
t8_cmesh_ghost_hash_fn (const void *ghost, const void *data)
{
  t8_cmesh_t          cmesh;
  t8_cghost_t         G;

  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (cmesh->num_ghosts > 0);
  T8_ASSERT (cmesh->set_partition);
  T8_ASSERT (cmesh->num_local_trees > 0);

  G = (t8_cghost_t) ghost;
  /* TODO: is this a reasonable hash value? */
  return G->treeid % cmesh->num_ghosts;
}

static int
t8_cmesh_ghost_equal_fn (const void *ghost1, const void *ghost2,
                         const void *data)
{
  t8_cghost_t         G1, G2;

  G1 = (t8_cghost_t) ghost1;
  G2 = (t8_cghost_t) ghost2;

  return G1->treeid == G2->treeid;
}
#endif

/* Check whether a given communicator assigns the same rank and mpisize
 * as stored in a given cmesh. */
int
t8_cmesh_comm_is_valid (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 mpiret, mpisize, mpirank;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  if (mpisize != cmesh->mpisize || mpirank != cmesh->mpirank) {
    return 0;
  }
  return 1;
}

void
t8_cmesh_init (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (pcmesh != NULL);

  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  t8_refcount_init (&cmesh->rc);

  /* sensible (hard error) defaults */
  cmesh->set_refine_level = 0;  /*< sensible default TODO document */
  cmesh->set_partition_level = -1;
  cmesh->dimension = -1;        /*< ok; force user to select dimension */
  cmesh->mpirank = -1;
  cmesh->mpisize = -1;
  cmesh->face_knowledge = 3;    /*< sensible default TODO document */
  t8_stash_init (&cmesh->stash);

  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
}

#if 0
/* This function is not part of the interface. The number of trees is always clear
 * from the number of calls to t8_cmesh_set_tree_class.
 * It is set in t8_cmesh_commit */
/* TODO: rename num_trees to global_num_trees or num_gtrees etc.
 *       to always distinguish between local and global.
 *       Do this everywhere in the code.
 */
static void
t8_cmesh_set_num_trees (t8_cmesh_t cmesh, t8_gloidx_t num_trees)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  /* If the cmesh is entered as a partitioned cmesh,
   * this function sets the local number of trees;
   * (TODO ^^^ would require locidx -- better provide two arguments)
   * the global number then must have been set in cmesh_set_partition.
   * Otherwise the global number of trees is set here.
   * TODO: make this function behave consistently independent on prior
   *       calls to set_partition.
   *       We want the user to be free in the sequence of calls
   *       as much as possible.
   */
  if (cmesh->set_partition) {
    /* num_trees == 0 is allowed */
    T8_ASSERT (cmesh->num_trees > 0);
    T8_ASSERT (cmesh->num_local_trees == 0);
    cmesh->num_local_trees = num_trees;
  }
  else {
    /* num_trees == 0 is allowed */
    T8_ASSERT (cmesh->num_trees >= 0);
    cmesh->num_trees = cmesh->num_local_trees = num_trees;
  }
  /* As soon as we know the number of trees, we allocate
   * the ctree array.
   * TODO?
   */
}
#endif

void
t8_cmesh_set_derive (t8_cmesh_t cmesh, t8_cmesh_t set_from)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (set_from == NULL || t8_cmesh_is_committed (set_from));

  if (cmesh->set_from != NULL) {
    /* If we overwrite a previously set cmesh, then we unref it. */
    t8_cmesh_unref (&cmesh->set_from);
  }
  cmesh->set_from = set_from;

  if (set_from != NULL) {
    t8_cmesh_set_dimension (cmesh, set_from->dimension);
  }
}

#if 0
/* TODO: deprecated, remove */
void
t8_cmesh_set_partition (t8_cmesh_t cmesh, int set_partition,
                        int set_face_knowledge,
                        t8_gloidx_t first_local_tree,
                        t8_gloidx_t last_local_tree,
                        t8_gloidx_t * tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (0 <= set_face_knowledge && set_face_knowledge <= 3);
  /* TODO: allow -1 for set_face_knowledge to keep it unchanged?
   *      update: unchanged from what? face_knowledge is only important for the
   * information on the stash. When the cmesh is derived there is no
   * stash. A committed cmesh has always face_knowledge 3. */

  /* TODO: Careful with tese assumptions; allow the user maximum flexibility */
#if 0
  T8_ASSERT (cmesh->num_trees == 0);
  T8_ASSERT (cmesh->num_local_trees == 0);
  T8_ASSERT (cmesh->first_tree == 0);
#endif

  /* set cmesh->set_partition to 0 or 1 (no; we always treat nonzero as true) */
  cmesh->set_partition = set_partition;
  /* TODO: this is how to query boolean variables */
  if (set_partition) {
    cmesh->first_tree = first_local_tree;
    cmesh->num_local_trees = last_local_tree - first_local_tree + 1;
    /* Since num_local_trees is a locidx we have to check whether we did create an
     * overflow in the previous computation */
    T8_ASSERT (cmesh->num_local_trees ==
               last_local_tree - first_local_tree + 1);
    cmesh->face_knowledge = set_face_knowledge;
    /* Right now no other face_knowledge is supported */
    SC_CHECK_ABORTF (set_face_knowledge == 3, "Level %i of face knowledge"
                     "is not supported.\n", set_face_knowledge);
    cmesh->tree_offsets = tree_offsets;
  }
}
#endif

t8_shmem_array_t
t8_cmesh_alloc_offsets (int mpisize, sc_MPI_Comm comm)
{
  t8_shmem_array_t    offsets;
#ifdef T8_ENABLE_DEBUG
  int                 mpisize_debug, mpiret;
  mpiret = sc_MPI_Comm_size (comm, &mpisize_debug);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (mpisize == mpisize_debug);
  t8_debugf ("Allocating shared array with type %s\n",
             sc_shmem_type_to_string[sc_shmem_get_type (comm)]);
#endif

  t8_shmem_array_init (&offsets, sizeof (t8_gloidx_t), mpisize + 1, comm);
  return offsets;
}

void
t8_cmesh_set_partition_range (t8_cmesh_t cmesh, int set_face_knowledge,
                              t8_gloidx_t first_local_tree,
                              t8_gloidx_t last_local_tree)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  SC_CHECK_ABORT (set_face_knowledge == -1 || set_face_knowledge == 3,
                  "Face knowledge other than three is not implemented yet.");
  cmesh->face_knowledge = set_face_knowledge;
  cmesh->first_tree = first_local_tree;
  cmesh->num_local_trees = last_local_tree - first_local_tree + 1;
  cmesh->set_partition = 1;
  /* Overwrite previous partition settings */
  if (cmesh->tree_offsets != NULL) {
    t8_shmem_array_destroy (&cmesh->tree_offsets);
    cmesh->tree_offsets = NULL;
  }
  cmesh->set_partition_level = -1;
}

void
t8_cmesh_set_partition_offsets (t8_cmesh_t cmesh,
                                t8_shmem_array_t tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  if (cmesh->tree_offsets != NULL && cmesh->tree_offsets != tree_offsets) {
    /* We overwrite a previouly set offset array, so
     * we need to free its memory first. */
    t8_shmem_array_destroy (&cmesh->tree_offsets);
  }
  cmesh->tree_offsets = tree_offsets;
  cmesh->set_partition = 1;
  if (tree_offsets != NULL) {
    /* We overwrite any previously partition settings */
    cmesh->first_tree = -1;
    cmesh->num_local_trees = -1;
    cmesh->set_partition_level = -1;
  }
}

void
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, int element_level)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (element_level >= -1);

  cmesh->set_partition = 1;
  cmesh->set_partition_level = element_level;
  if (element_level >= 0) {
    /* We overwrite any previous partition settings */
    cmesh->first_tree = -1;
    cmesh->num_local_trees = -1;
    if (cmesh->tree_offsets != NULL) {
      t8_shmem_array_destroy (&cmesh->tree_offsets);
      cmesh->tree_offsets = NULL;
    }
  }
}

#if 0
/* No longer needed */
void
t8_cmesh_set_partition_from (t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from,
                             int level, t8_gloidx_t * tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (t8_cmesh_is_committed (cmesh_from));
  T8_ASSERT (cmesh_from->set_partition);

  cmesh->set_from = cmesh_from;
  cmesh->set_partition = 1;
  cmesh->face_knowledge = cmesh_from->face_knowledge;
  if (level >= 0) {
    cmesh->set_partition_level = level;
  }
  else {
    cmesh->tree_offsets = tree_offsets;
  }
  cmesh->from_method |= T8_CMESH_PARTITION;
}
#endif

void
t8_cmesh_set_refine (t8_cmesh_t cmesh, int level, t8_scheme_cxx_t * scheme)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (level >= 0);
  T8_ASSERT (scheme != NULL);

  cmesh->set_refine_level = level;
  cmesh->set_refine_scheme = scheme;
}

t8_gloidx_t
t8_cmesh_get_first_treeid (t8_cmesh_t cmesh)
{
  return cmesh->first_tree;
}

/* TODO: should get a gloidx?
 *       place after commit */
t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_locidx_t ltree_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (0 <= ltree_id && ltree_id < cmesh->num_local_trees);

  return t8_cmesh_trees_get_tree (cmesh->trees, ltree_id);
}

/* Returns the first local tree.
 * Returns NULL if there are no local trees. */
/* TODO: hide */
t8_ctree_t
t8_cmesh_get_first_tree (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->num_local_trees > 0 ? t8_cmesh_get_tree (cmesh, 0) : NULL;
}

/* returns the next local tree in the cmesh (by treeid)
 * after a given tree.
 * The given tree must be a valid and owned tree.
 * If the given tree is the last local tree, NULL is returned */
/* TODO: hide */
t8_ctree_t
t8_cmesh_get_next_tree (t8_cmesh_t cmesh, t8_ctree_t tree)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= tree->treeid && tree->treeid < cmesh->num_local_trees);
  T8_ASSERT (cmesh->committed);
  return tree->treeid <
    cmesh->num_local_trees -
    1 ? t8_cmesh_get_tree (cmesh, tree->treeid + 1) : NULL;
}

void
t8_cmesh_set_attribute (t8_cmesh_t cmesh, t8_gloidx_t gtree_id,
                        int package_id, int key, void *data, size_t data_size,
                        int data_persists)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, gtree_id, package_id, key, data_size,
                          data, data_persists);
}

void               *
t8_cmesh_get_attribute (t8_cmesh_t cmesh, int package_id, int key,
                        t8_locidx_t ltree_id)
{
  int                 is_ghost;

  T8_ASSERT (cmesh->committed);
  T8_ASSERT (0 <= ltree_id &&
             ltree_id < cmesh->num_ghosts + cmesh->num_local_trees);
  is_ghost = ltree_id >= cmesh->num_local_trees;
  if (is_ghost) {
    ltree_id = ltree_id - cmesh->num_local_trees;
  }
  return t8_cmesh_trees_get_attribute (cmesh->trees, ltree_id, package_id,
                                       key, NULL, is_ghost);
}

t8_shmem_array_t
t8_cmesh_get_partition_table (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (!cmesh->set_partition) {
    /* The mesh is not partitioned. We return NULL. */
    return NULL;
  }
  /* If the mesh is not stored, NULL is returned, otherwise the
   * partition array. */
  return cmesh->tree_offsets;
}

#if 0
/* Check whether a given tree_id belongs to a tree in the cmesh.
 * If partitioned only local trees are allowed.
 */
static int
t8_cmesh_tree_id_is_owned (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  T8_ASSERT (cmesh->committed);
  if (cmesh->set_partition) {
    return cmesh->first_tree <= tree_id
      && tree_id < cmesh->first_tree + cmesh->num_local_trees;
  }
  else {
    return 0 <= tree_id && tree_id < cmesh->num_trees;
  }
}

#endif

#if 0
/* Given a tree_id return the index of the specified tree in
 * cmesh's tree array
 */
static t8_locidx_t
t8_cmesh_tree_index (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  return cmesh->set_partition ? tree_id - cmesh->first_tree : tree_id;
}
#endif

void
t8_cmesh_set_dimension (t8_cmesh_t cmesh, int dim)
{
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  T8_ASSERT (0 <= dim && dim <= T8_ECLASS_MAX_DIM);

  cmesh->dimension = dim;
}

void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, t8_gloidx_t gtree_id,
                         t8_eclass_t tree_class)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (gtree_id >= 0);

  /* If we insert the first tree, set the dimension of the cmesh
   * to this tree's dimension. Otherwise check whether the dimension
   * of the tree to be inserted equals the dimension of the cmesh. */
  if (cmesh->dimension == -1) {
    cmesh->dimension = t8_eclass_to_dimension[tree_class];
  }
  else {
    /* TODO: This makes it illegal to set a tree to i.e. quad and change it
     *       to hex later. Even if we replace all trees with another dimension.
     *       We could move this check to commit. */
    /* TODO: If cmesh is partitioned and this part has no trees then the
     *       dimension remains unset forever. */
    T8_ASSERT (t8_eclass_to_dimension[tree_class] == cmesh->dimension);
  }

  t8_stash_add_class (cmesh->stash, gtree_id, tree_class);
#ifdef T8_ENABLE_DEBUG
  cmesh->inserted_trees++;
#endif
}

void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltree_id,
                            int package_id, int key,
                            double *vertices, int num_vertices)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (vertices != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, ltree_id, package_id, key,
                          3 * num_vertices * sizeof (double),
                          (void *) vertices, 1);
}

void
t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t gtree1, t8_gloidx_t gtree2,
                   int face1, int face2, int orientation)
{
  T8_ASSERT (0 <= orientation);

  t8_stash_add_facejoin (cmesh->stash, gtree1, gtree2, face1, face2,
                         orientation);
}

void
t8_cmesh_set_profiling (t8_cmesh_t cmesh, int set_profiling)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  if (set_profiling) {
    if (cmesh->profile == NULL) {
      /* Only do something if profiling is not enabled already */
      cmesh->profile = T8_ALLOC_ZERO (t8_cprofile_struct_t, 1);
    }
  }
  else {
    /* Free any profile that is already set */
    if (cmesh->profile != NULL) {
      T8_FREE (cmesh->profile);
    }
  }
}

/* returns true if cmesh_a equals cmesh_b */
int
t8_cmesh_is_equal (t8_cmesh_t cmesh_a, t8_cmesh_t cmesh_b)
/* TODO: rewrite */
{
  int                 is_equal;
  T8_ASSERT (cmesh_a != NULL && cmesh_b != NULL);

  if (cmesh_a == cmesh_b) {
    return 1;
  }
  /* check entries that are numbers */
  is_equal = cmesh_a->committed != cmesh_b->committed || cmesh_a->dimension !=
    cmesh_b->dimension ||
    cmesh_a->set_partition != cmesh_b->set_partition ||
    cmesh_a->mpirank != cmesh_b->mpirank ||
    cmesh_a->mpisize != cmesh_b->mpisize ||
    cmesh_a->num_trees != cmesh_b->num_trees ||
    cmesh_a->num_local_trees != cmesh_b->num_local_trees ||
    cmesh_a->num_ghosts != cmesh_b->num_ghosts ||
    cmesh_a->first_tree != cmesh_b->first_tree;
#if 0
  /* TODO: The inserted variables are counters that are only active if the
   * cmesh is committed from scratch. If a cmesh is commited via cmesh_copy,
   * then these counters are not active. So even for equal cmeshes
   * these counters must not store the same value. */
#ifdef T8_ENABLE_DEBUG
  is_equal = is_equal || cmesh_a->inserted_trees != cmesh_b->inserted_trees ||
    cmesh_a->inserted_ghosts != cmesh_b->inserted_ghosts;
#endif
#endif
  if (is_equal != 0) {
    return 0;
  }
  /* check arrays */
  is_equal = memcmp (cmesh_a->num_trees_per_eclass,
                     cmesh_b->num_trees_per_eclass,
                     T8_ECLASS_COUNT * sizeof (t8_gloidx_t));

  /* check tree_offsets */
  if (cmesh_a->tree_offsets != NULL) {
    if (cmesh_b->tree_offsets == NULL) {
      return 0;
    }
    else {
      is_equal = is_equal || !t8_shmem_array_is_equal (cmesh_a->tree_offsets,
                                                       cmesh_b->tree_offsets);
    }
  }
  if (is_equal != 0) {
    return 0;
  }
  /* check trees */
  if (cmesh_a->committed &&
      !t8_cmesh_trees_is_equal (cmesh_a, cmesh_a->trees, cmesh_b->trees)) {
    /* if we have committed check tree arrays */
    return 0;
  }
  else {
    if (!cmesh_a->committed &&
        !t8_stash_is_equal (cmesh_a->stash, cmesh_b->stash)) {
      /* if we have not committed check stash arrays */
      return 0;
    }
  }
  return 1;
}

#if 0
/* broadcast the tree attributes of a cmesh on root to all processors */
/* TODO: can we optimize it by just sending the memory of the mempools? */
static void
t8_cmesh_bcast_attributes (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int                 mpirank, mpisize, mpiret;
  int                 has_attr;
  t8_ctree_t          tree;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  for (tree = t8_cmesh_get_first_tree (cmesh_in); tree != NULL;
       tree = t8_cmesh_get_next_tree (cmesh_in, tree)) {
    if (mpirank == root && tree->attribute != NULL) {
      has_attr = 1;
    }
    else {
      has_attr = 0;
    }
    mpiret = sc_MPI_Bcast (&has_attr, 1, sc_MPI_INT, root, comm);
    SC_CHECK_MPI (mpiret);
    if (has_attr) {
      if (mpirank != root) {
        tree->attribute =
          sc_mempool_alloc (cmesh_in->tree_attributes_mem[tree->eclass]);
      }
      mpiret = sc_MPI_Bcast (tree->attribute,
                             t8_cmesh_get_attribute_size (cmesh_in,
                                                          tree->eclass),
                             sc_MPI_BYTE, root, comm);
      SC_CHECK_MPI (mpiret);
    }
  }
}
#endif

t8_cmesh_t
t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int                 mpirank, mpisize, mpiret;
  int                 iclass;

  struct
  {
    int                 dimension;
    t8_locidx_t         num_trees;
    t8_gloidx_t         num_trees_per_eclass[T8_ECLASS_COUNT];
    size_t              stash_elem_counts[3];
#ifdef T8_ENABLE_DEBUG
    t8_locidx_t         inserted_trees;
    sc_MPI_Comm         comm;
#endif
  } dimensions;

  /* TODO: BUG: running with two processes and a cmesh of one T8_ECLASS_LINE,
   *       the on both processes the face_neigbors and vertices arrays of
   *       the single tree point to the same physical memory.
   *       (face_neighbors on both processes are equal and vertices on both
   *        processes are equal)
   */
  /* TODO: Send the tree's vertices */

  /* TODO: rewrite */

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (0 <= root && root < mpisize);
  T8_ASSERT (mpirank == root || cmesh_in == NULL);
  T8_ASSERT (mpirank != root || cmesh_in != NULL);
  T8_ASSERT (mpirank != root || cmesh_in->set_partition == 0);
  /* The cmesh on the calling process must not be owned by something
   * else. */
  /* TODO: would it be useful to allow bcast even if the cmesh is referenced?
   * But then the bcasted version on other procs would have a different refcount
   * than the cmesh on the root */
  T8_ASSERT (mpirank != root || cmesh_in->rc.refcount == 1);

  /* At first we broadcast all meta information. */
  if (mpirank == root) {
    dimensions.dimension = cmesh_in->dimension;
    dimensions.num_trees = cmesh_in->num_trees;
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      dimensions.num_trees_per_eclass[iclass] =
        cmesh_in->num_trees_per_eclass[iclass];
    }
    dimensions.stash_elem_counts[0] = cmesh_in->stash->attributes.elem_count;
    dimensions.stash_elem_counts[1] = cmesh_in->stash->classes.elem_count;
    dimensions.stash_elem_counts[2] = cmesh_in->stash->joinfaces.elem_count;
#ifdef T8_ENABLE_DEBUG
    dimensions.inserted_trees = cmesh_in->inserted_trees;
#endif
  }
  /* TODO: we could optimize this by using IBcast */
  mpiret = sc_MPI_Bcast (&dimensions, sizeof (dimensions), sc_MPI_BYTE, root,
                         comm);
  SC_CHECK_MPI (mpiret);

  /* If not root store information in new cmesh and allocate memory for arrays. */
  if (mpirank != root) {
    t8_cmesh_init (&cmesh_in);
    cmesh_in->dimension = dimensions.dimension;
    cmesh_in->num_trees = dimensions.num_trees;
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      cmesh_in->num_trees_per_eclass[iclass] =
        dimensions.num_trees_per_eclass[iclass];
    }
#ifdef T8_ENABLE_DEBUG
    cmesh_in->inserted_trees = dimensions.inserted_trees;
    T8_ASSERT (dimensions.comm == comm);
#endif
  }
  /* broadcast all the stashed information about trees/neighbors/attributes */
  t8_stash_bcast (cmesh_in->stash, root, comm, dimensions.stash_elem_counts);
  return cmesh_in;
}

int
t8_cmesh_face_is_boundary (t8_cmesh_t cmesh, t8_locidx_t ltreeid, int face)
{
  int8_t             *ttf;
  int                 F;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (0 <= ltreeid
             && ltreeid < cmesh->num_ghosts + cmesh->num_local_trees);

  F = t8_eclass_max_num_faces[cmesh->dimension];
  if (ltreeid < cmesh->num_local_trees) {
    t8_locidx_t        *neighbors;
    /* The tree is a local tree */
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, ltreeid, &neighbors,
                                        &ttf);
    if (neighbors[face] == ltreeid && ttf[face] % F == face) {
      /* This face is a boundary face */
      return 1;
    }
  }
  else if (ltreeid < cmesh->num_local_trees + cmesh->num_ghosts) {
    /* The tree is a ghost */
    t8_cghost_t         ghost;
    t8_gloidx_t        *neighbors;
    ghost = t8_cmesh_trees_get_ghost_ext (cmesh->trees,
                                          ltreeid - cmesh->num_local_trees,
                                          &neighbors, &ttf);
    if (neighbors[face] == ghost->treeid && ttf[face] % F == face) {
      /* this is a boundary face */
      return 1;
    }
  }
  return 0;
}

#ifdef T8_WITH_METIS
void
t8_cmesh_reorder (t8_cmesh_t cmesh, sc_MPI_Comm comm, idx_t num_partitions)
{
  idx_t               ncon = 1, elemens;
  idx_t               volume, *partition, ipart, newpart;
  int                 num_faces, iface, count_face;
  idx_t              *xadj, *adjncy;
  idx_t               options[METIS_NOPTIONS];
  int                 success;
  t8_locidx_t        *new_number, itree, *tree_per_part_off, *tree_per_part;
  t8_locidx_t        *face_neighbor;
  t8_eclass_t         tree_class;

  t8_debugf ("Starting METIS\n");

  /* cmesh must be commited and not partitioned */
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (!cmesh->set_partition);

  elemens = cmesh->num_trees;
  T8_ASSERT ((t8_locidx_t) elemens == cmesh->num_trees);

  if (num_partitions == 1) {
    /* Only one partition, no repartioning is performed */
    return;
  }
  /* Count the number of tree-to-tree connections via a face */
  num_faces = 0;
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    for (iface = 0; iface < t8_eclass_num_faces[tree_class]; iface++) {
      if (!t8_cmesh_face_is_boundary (cmesh, itree, iface)) {
        num_faces++;
      }
    }
  }

  /* xadj and adjncy store the face-connections in a CSR format
   * xadj[treeid] = offset of the tree in adjncy
   * adjncy[xadj[treeid]]...adjncy[xadj[treeid]-1] are the trees with which
   * the tree has a face connection */
  xadj = T8_ALLOC_ZERO (idx_t, elemens + 1);
  adjncy = T8_ALLOC (idx_t, num_faces);

  /* fill xadj and adjncy arrays */
  for (itree = 0, count_face = 0; itree < cmesh->num_trees; itree++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor,
                                        NULL);
    xadj[itree + 1] = xadj[itree];
    for (iface = 0; iface < t8_eclass_num_faces[tree_class]; iface++) {
      if (!t8_cmesh_face_is_boundary (cmesh, itree, iface)) {
        adjncy[count_face++] = face_neighbor[iface];
        xadj[itree + 1]++;
      }
    }
  }

  T8_ASSERT (count_face == num_faces);

  /* partition stores the new partition number for each element */
  partition = T8_ALLOC (idx_t, elemens);
  /* Set the default Metis options */
  METIS_SetDefaultOptions (options);
  /* use c-style 0 indexed numbering */
  options[METIS_OPTION_NUMBERING] = 0;
  /* partition the elements in mpisize many partitions */
#if 1
  success =
    METIS_PartGraphRecursive (&elemens, &ncon, xadj, adjncy, NULL, NULL, NULL,
                              &num_partitions, NULL, NULL, options, &volume,
                              partition);
#else
  success =
    METIS_PartGraphKway (&elemens, &ncon, xadj, adjncy, NULL, NULL, NULL,
                         &idx_mpisize, NULL, NULL, options, &volume,
                         partition);
#endif

  T8_ASSERT (success == METIS_OK);
  /* memory to store the new treeid of a tree */
  new_number = T8_ALLOC (t8_locidx_t, cmesh->num_trees);
  /* Store the number of trees per partition */
  tree_per_part = T8_ALLOC_ZERO (t8_locidx_t, num_partitions);
  /* Store the treeid offset of each partition. */
  tree_per_part_off = T8_ALLOC_ZERO (t8_locidx_t, num_partitions + 1);
  tree_per_part_off[0] = 0;
  /* compute tree_per_part and prepare tree_per_part_off */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree_per_part[partition[itree]]++;
    tree_per_part_off[partition[itree] + 1]++;
  }
  /* compute tree_per_part_off */
  for (ipart = 1; ipart <= num_partitions; ipart++) {
    tree_per_part_off[ipart] += tree_per_part_off[ipart - 1];
  }
  /* Compute for each tree its new treeid */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    newpart = partition[itree];
    T8_ASSERT (tree_per_part[newpart] > 0);
    new_number[itree] =
      tree_per_part_off[newpart + 1] - tree_per_part[newpart];
    tree_per_part[newpart]--;
  }

  /* Reorder the trees */
  t8_cmesh_trees_reorder (cmesh, cmesh->trees, new_number);

  T8_FREE (partition);
  T8_FREE (xadj);
  T8_FREE (adjncy);
  T8_FREE (new_number);
  T8_FREE (tree_per_part);
  T8_FREE (tree_per_part_off);

  t8_debugf ("End METIS\n");
}
#endif

int
t8_cmesh_first_tree_is_shared (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->first_tree_shared;
}

t8_gloidx_t
t8_cmesh_get_num_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_trees;
}

t8_locidx_t
t8_cmesh_get_num_local_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->num_local_trees;
}

t8_locidx_t
t8_cmesh_get_num_ghosts (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->num_ghosts;
}

t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_locidx_t ltree_id)
{
  t8_ctree_t          tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  tree = t8_cmesh_get_tree (cmesh, ltree_id);
  return tree->eclass;
}

t8_eclass_t
t8_cmesh_get_ghost_class (t8_cmesh_t cmesh, t8_locidx_t lghost_id)
{
  t8_cghost_t         ghost;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (0 <= lghost_id && lghost_id < cmesh->num_ghosts);

  ghost = t8_cmesh_trees_get_ghost (cmesh->trees, lghost_id);
  return ghost->eclass;
}

t8_gloidx_t
t8_cmesh_get_global_id (t8_cmesh_t cmesh, t8_locidx_t local_id)
{
  T8_ASSERT (0 <= local_id && local_id <
             cmesh->num_ghosts + cmesh->num_local_trees);
  if (local_id < cmesh->num_local_trees) {
    return local_id + cmesh->first_tree;
  }
  else {
    return t8_cmesh_trees_get_ghost (cmesh->trees,
                                     local_id -
                                     cmesh->num_local_trees)->treeid;
  }
}

int
t8_cmesh_tree_is_local (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  if (0 <= ltreeid && ltreeid < t8_cmesh_get_num_local_trees (cmesh)) {
    return 1;
  }
  return 0;
}

t8_locidx_t
t8_cmesh_get_local_id (t8_cmesh_t cmesh, t8_gloidx_t global_id)
{
  t8_gloidx_t         temp_local_id;
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (0 <= global_id && global_id < cmesh->num_trees);

  if (!cmesh->set_partition) {
    /* If the cmesh is not partitioned the local id is the global id */
    return global_id;
  }
  temp_local_id = global_id - cmesh->first_tree;
  if (0 <= temp_local_id && temp_local_id < cmesh->num_local_trees) {
    /* The tree is a local tree */
    return temp_local_id;
  }
  else {
    /* The tree may be a ghost tree */
    return t8_cmesh_trees_get_ghost_local_id (cmesh->trees, global_id);
  }
}

void
t8_cmesh_print_profile (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (cmesh->profile != NULL) {
    /* Only print something if profiling is enabled */
    sc_statinfo_t       stats[T8_CPROFILE_NUM_STATS];
    t8_cprofile_t      *profile = cmesh->profile;

    /* Set the stats */
    sc_stats_set1 (&stats[0], profile->partition_trees_shipped,
                   "cmesh: Number of trees sent.");
    sc_stats_set1 (&stats[1],
                   profile->partition_ghosts_shipped,
                   "cmesh: Number of ghosts sent.");
    sc_stats_set1 (&stats[2], profile->partition_trees_recv,
                   "cmesh: Number of trees received.");
    sc_stats_set1 (&stats[3], profile->partition_ghosts_recv,
                   "cmesh: Number of ghosts received.");
    sc_stats_set1 (&stats[4], profile->partition_bytes_sent,
                   "cmesh: Number of bytes sent.");
    sc_stats_set1 (&stats[5], profile->partition_procs_sent,
                   "cmesh: Number of processes sent to.");
    sc_stats_set1 (&stats[6], profile->first_tree_shared,
                   "cmesh: First tree is shared.");
    sc_stats_set1 (&stats[7], profile->partition_runtime,
                   "cmesh: Partition runtime.");
    sc_stats_set1 (&stats[8], profile->commit_runtime,
                   "cmesh: Commit runtime.");
    /* compute stats */
    sc_stats_compute (sc_MPI_COMM_WORLD, T8_CPROFILE_NUM_STATS, stats);
    /* print stats */
    t8_logf (SC_LC_GLOBAL, SC_LP_STATISTICS, "Printing stats for cmesh.\n");
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS,
                    T8_CPROFILE_NUM_STATS, stats, 1, 1);
  }
}

void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                         t8_gloidx_t * first_local_tree,
                         t8_gloidx_t * child_in_tree_begin,
                         t8_gloidx_t * last_local_tree,
                         t8_gloidx_t * child_in_tree_end,
                         int8_t * first_tree_shared)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (level >= 0);

  *first_local_tree = 0;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = 0;
  }
  *last_local_tree = 0;
  if (child_in_tree_end != NULL) {
    *child_in_tree_end = 0;
  }

  if (cmesh->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0 || level == 0) {
    t8_gloidx_t         global_num_children;
    t8_gloidx_t         first_global_child;
    t8_gloidx_t         last_global_child;
    t8_gloidx_t         children_per_tree;
    t8_gloidx_t         prev_last_tree = -1;
    const uint64_t      one = 1;

    children_per_tree = one << cmesh->dimension * level;
    global_num_children = cmesh->num_trees * children_per_tree;

    if (cmesh->mpirank == 0) {
      first_global_child = 0;
      if (child_in_tree_begin != NULL) {
        *child_in_tree_begin = 0;
      }
    }
    else {
      /* The first global child of processor p
       * with P total processor is (the biggest int smaller than)
       * (total_num_children * p) / P
       * We cast to long double and double first to prevent integer overflow.
       */
      first_global_child =
        ((long double) global_num_children *
         cmesh->mpirank) / (double) cmesh->mpisize;
    }
    if (cmesh->mpirank != cmesh->mpisize - 1) {
      last_global_child =
        ((long double) global_num_children *
         (cmesh->mpirank + 1)) / (double) cmesh->mpisize;
    }
    else {
      last_global_child = global_num_children;
    }

    SC_CHECK_ABORT (first_global_child < last_global_child,
                    "forest does not support empty processes yet");
    T8_ASSERT (0 <= first_global_child
               && first_global_child <= global_num_children);
    T8_ASSERT (0 <= last_global_child
               && last_global_child <= global_num_children);
    *first_local_tree = first_global_child / children_per_tree;
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin =
        first_global_child - *first_local_tree * children_per_tree;
    }
    /* TODO: Just fixed this line from last_global_child -1 / cpt
     *       Why did we not notice this error before?
     *       Changed it back*/
    *last_local_tree = (last_global_child - 1) / children_per_tree;
    if (first_tree_shared != NULL) {
      prev_last_tree = (first_global_child - 1) / children_per_tree;
      T8_ASSERT (cmesh->mpirank > 0 || prev_last_tree <= 0);
      if (cmesh->mpirank > 0 && prev_last_tree == *first_local_tree &&
          first_global_child < last_global_child && last_global_child >= 0) {
        /* We exclude empty partitions here, by def their first_tree_shared flag is zero */
        /* We also exclude that the previous partition was empty at the beginning of the
         * partitions array */
        /* TODO: If empty partitions in the middle can occur then we have to think this over */
        *first_tree_shared = 1;
      }
      else {
        *first_tree_shared = 0;
      }
    }
    if (child_in_tree_end != NULL) {
      if (*last_local_tree > 0) {
        *child_in_tree_end =
          last_global_child - *last_local_tree * children_per_tree;
      }
      else {
        *child_in_tree_end = last_global_child;
      }
    }
    if (first_global_child >= last_global_child && cmesh->mpirank != 0) {
      /* This process is empty */
      *first_local_tree = prev_last_tree + 1;
    }
  }
  else {
    SC_ABORT ("Partition with level > 0 "
              "does not support pyramidal elements yet.");
  }
}

static void
t8_cmesh_reset (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  /* free tree_offset */
  if (cmesh->tree_offsets != NULL) {
#if T8_ENABLE_DEBUG
    sc_MPI_Comm         comm;
    /* Check whether a correct communicator was stored at tree_offsets.
     * This is useful for debugging. */
    comm = t8_shmem_array_get_comm (cmesh->tree_offsets);
    T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
#endif
    /* Destroy the shared memory array */
    t8_shmem_array_destroy (&cmesh->tree_offsets);
  }
  /*TODO: write this */
  if (!cmesh->committed) {
    t8_stash_destroy (&cmesh->stash);
    if (cmesh->set_from != NULL) {
      /* We unref our reference of set_from */
      t8_cmesh_unref (&cmesh->set_from);
    }
  }
  else {
    if (cmesh->trees != NULL) {
      t8_cmesh_trees_destroy (&cmesh->trees);
    }
    T8_ASSERT (cmesh->set_from == NULL);
  }
  if (cmesh->profile != NULL) {
    T8_FREE (cmesh->profile);
  }
  if (cmesh->set_refine_scheme != NULL) {
    t8_scheme_cxx_unref (&cmesh->set_refine_scheme);
  }

  T8_FREE (cmesh);
  *pcmesh = NULL;
}

void
t8_cmesh_ref (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  t8_refcount_ref (&cmesh->rc);
}

void
t8_cmesh_unref (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  if (t8_refcount_unref (&cmesh->rc)) {
    t8_cmesh_reset (pcmesh);
  }
}

void
t8_cmesh_destroy (t8_cmesh_t * pcmesh)
{
  T8_ASSERT (pcmesh != NULL && *pcmesh != NULL &&
             t8_refcount_is_last (&(*pcmesh)->rc));
  t8_cmesh_unref (pcmesh);
  T8_ASSERT (*pcmesh == NULL);
}

/* TODO: In p4est a tree edge is joined with itself to denote a domain boundary.
 *       Will we do it the same in t8code? This is not yet decided, however the
 *       function below stores these neighbourhood information in the cmesh. */
/* TODO: Eventually we may directly partition the mesh here */
/* Offset-1 is added to each tree_id, this is used in i.e. t8_cmesh_new_disjoint_bricks,
 * If offset is nonzero, then set_partition must be true and the cmesh is
 * partitioned and has all trees in conn as local trees.
 * The offsets on the different processes must add up! */
static              t8_cmesh_t
t8_cmesh_new_from_p4est_ext (void *conn, int dim,
                             sc_MPI_Comm comm, int set_partition,
                             t8_gloidx_t offset)
{
#define _T8_CMESH_P48_CONN(_ENTRY) \
  (dim == 2 ? ((p4est_connectivity_t *) conn)->_ENTRY \
            : ((p8est_connectivity_t *) conn)->_ENTRY)
  t8_cmesh_t          cmesh;
  t8_gloidx_t         ltree;
  p4est_topidx_t      treevertex;
  double              vertices[24];     /* Only 4 * 3 = 12 used in 2d */
  int                 num_tvertices;
  int                 num_faces;
  int                 ivertex, iface;
  int                 use_offset;
  int8_t              ttf;
  p4est_topidx_t      ttt;

  T8_ASSERT (dim == 2 || dim == 3);
  T8_ASSERT (dim == 3
             ||
             p4est_connectivity_is_valid ((p4est_connectivity_t *) (conn)));
  T8_ASSERT (dim == 2
             ||
             p8est_connectivity_is_valid ((p8est_connectivity_t *) (conn)));
  T8_ASSERT (offset == 0 || set_partition);
  if (offset) {
    offset--;
    use_offset = 1;
  }
  else {
    use_offset = 0;
  }
  T8_ASSERT (offset >= 0);
  /* TODO: Check offsets for consistency */
  num_tvertices = 1 << dim;     /*vertices per tree. 4 if dim = 2 and 8 if dim = 3. */
  num_faces = dim == 2 ? 4 : 6;
  /* basic setup */
  t8_cmesh_init (&cmesh);
  /* Add each tree to cmesh and get vertex information for each tree */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) {    /* loop over each tree */
    t8_cmesh_set_tree_class (cmesh, ltree + offset,
                             dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX);
    for (ivertex = 0; ivertex < num_tvertices; ivertex++) {     /* loop over each tree corner */
      treevertex =
        _T8_CMESH_P48_CONN (tree_to_vertex[num_tvertices * ltree + ivertex]);
      vertices[3 * ivertex] = _T8_CMESH_P48_CONN (vertices[3 * treevertex]);
      vertices[3 * ivertex + 1] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 1]);
      vertices[3 * ivertex + 2] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 2]);
    }
    t8_cmesh_set_tree_vertices (cmesh, ltree + offset, t8_get_package_id (),
                                0, vertices, num_tvertices);
  }
  /* get face neighbor information from conn and join faces in cmesh */
  for (ltree = 0; ltree < _T8_CMESH_P48_CONN (num_trees); ltree++) {    /* loop over each tree */
    for (iface = 0; iface < num_faces; iface++) {       /* loop over each face */
      ttf = _T8_CMESH_P48_CONN (tree_to_face[num_faces * ltree + iface]);
      ttt = _T8_CMESH_P48_CONN (tree_to_tree[num_faces * ltree + iface]);
      /* insert the face only if we did not insert it before */
      if (ltree < ttt || (ltree == ttt && iface < ttf % num_faces)) {
        t8_cmesh_set_join (cmesh, ltree + offset, ttt + offset, iface,
                           ttf % num_faces, ttf / num_faces);
      }
    }
  }
  if (set_partition) {
    /* TODO: a copy of this code exists below, make it a function */
    int                 mpirank, mpisize, mpiret;
    t8_gloidx_t         first_tree, last_tree, num_trees, num_local_trees;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    if (use_offset == 0) {
      /* The total number of trees is the number of trees in conn */
      num_trees = _T8_CMESH_P48_CONN (num_trees);
      /* First tree and last tree according to uniform level 0 partitioning */
      first_tree = (mpirank * num_trees) / mpisize;
      last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    }
    else {
      /* First_tree and last_tree are the first and last trees of conn plu the offset */
      num_local_trees = _T8_CMESH_P48_CONN (num_trees);
      first_tree = offset;
      last_tree = offset + num_local_trees - 1;
      /* The global number of trees is the sum over all numbers of trees
       * in conn on each process */
      sc_MPI_Allreduce (&num_local_trees, &num_trees, 1, T8_MPI_GLOIDX,
                        sc_MPI_SUM, comm);
      t8_debugf ("[H] Generating partitioned cmesh from connectivity\n"
                 "[H] Has %li global and %li local trees.\n", num_trees,
                 num_local_trees);
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
#undef _T8_CMESH_P48_CONN
}

t8_cmesh_t
t8_cmesh_new_from_p4est (p4est_connectivity_t * conn,
                         sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 2, comm, do_partition, 0);
}

t8_cmesh_t
t8_cmesh_new_from_p8est (p8est_connectivity_t * conn,
                         sc_MPI_Comm comm, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 3, comm, do_partition, 0);
}

static t8_cmesh_t
t8_cmesh_new_vertex (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[3] = {
    0, 0, 0
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_VERTEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 1);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_line (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[6] = {
    0, 0, 0,
    1, 0, 0
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 2);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[9] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 3);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    1, 1, 1,
    1, -1, -1,
    -1, 1, -1,
    -1, -1, 1
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 8);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_pyramid (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[15] = {
    -1, -1, 0, 1, -1, 0, -1, 1, 0, 1, 1, 0, 0, 0, sqrt (2)
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                              vertices, 15);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

static t8_cmesh_t
t8_cmesh_new_prism (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[18] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_from_class (t8_eclass_t eclass, sc_MPI_Comm comm)
{
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return t8_cmesh_new_vertex (comm);
    break;
  case T8_ECLASS_LINE:
    return t8_cmesh_new_line (comm);
    break;
  case T8_ECLASS_TRIANGLE:
    return t8_cmesh_new_tri (comm);
    break;
  case T8_ECLASS_QUAD:
    return t8_cmesh_new_quad (comm);
    break;
  case T8_ECLASS_TET:
    return t8_cmesh_new_tet (comm);
    break;
  case T8_ECLASS_HEX:
    return t8_cmesh_new_hex (comm);
    break;
  case T8_ECLASS_PYRAMID:
    return t8_cmesh_new_pyramid (comm);
    break;
  case T8_ECLASS_PRISM:
    return t8_cmesh_new_prism (comm);
    break;
  default:
    SC_ABORT ("Invalid eclass\n");
    return NULL;
  }
}

t8_cmesh_t
t8_cmesh_new_empty (sc_MPI_Comm comm, int do_partition)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* TODO: This is just a helper function that was needed when we changed the vertex interface
 *       to use attributes. Before we stored a list of vertex coordinates in the cmesh and each tree indexed into this list.
 *       Now each tree carries the coordinates of its vertices.
 *       This function translates from the first approached to the second
 *       and was introduced to avoid rewritting the already existing cmesh_new... functions below.
 *       It would be nice to eventually rewrite these functions correctly.
 */
static void
t8_cmesh_new_translate_vertices_to_attributes (t8_topidx_t *
                                               tvertices,
                                               double
                                               *vertices,
                                               double
                                               *attr_vertices,
                                               int num_vertices)
{
  int                 i;
  for (i = 0; i < num_vertices; i++) {
    attr_vertices[3 * i] = vertices[3 * tvertices[i]];
    attr_vertices[3 * i + 1] = vertices[3 * tvertices[i] + 1];
    attr_vertices[3 * i + 2] = vertices[3 * tvertices[i] + 2];
  }
}

/* The unit cube is constructed from trees of the same eclass.
 * For triangles the square is divided along the (0,0) -- (1,1) diagonal.
 * For prisms the front (y=0) and back (y=1) face are divided into triangles
 * as above.
 */
/* TODO: upgrade with int x,y,z for periodic faces */
t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_bcast,
                        int do_partition)
{
  t8_cmesh_t          cmesh;
  int                 num_trees_for_hypercube[T8_ECLASS_COUNT] = {
    1, 1, 1, 2, 1, 6, 2, 3
  };
  int                 i;
  t8_topidx_t         vertices[8];
  double              attr_vertices[24];
  int                 mpirank, mpiret;
  double              vertices_coords[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (!do_bcast || mpirank == 0) {
    t8_cmesh_init (&cmesh);
    for (i = 0; i < num_trees_for_hypercube[eclass]; i++) {
      t8_cmesh_set_tree_class (cmesh, i, eclass);
    }
    switch (eclass) {
    case T8_ECLASS_HEX:
      vertices[4] = 4;
      vertices[5] = 5;
      vertices[6] = 6;
      vertices[7] = 7;
    case T8_ECLASS_QUAD:
      vertices[3] = 3;
      vertices[2] = 2;
    case T8_ECLASS_LINE:
      vertices[1] = 1;
    case T8_ECLASS_VERTEX:
      vertices[0] = 0;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices,
                                                     t8_eclass_num_vertices
                                                     [eclass]);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices,
                                  t8_eclass_num_vertices[eclass]);
      break;
    case T8_ECLASS_PRISM:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 2;
      vertices[4] = 3;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      vertices[1] = 5;
      vertices[2] = 4;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      break;
    case T8_ECLASS_TRIANGLE:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 3);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 3);
      break;
    case T8_ECLASS_TET:
      t8_cmesh_set_join (cmesh, 0, 1, 2, 2, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 3, 2, 2, 0);
      t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 0);
      t8_cmesh_set_join (cmesh, 4, 5, 2, 2, 0);
      t8_cmesh_set_join (cmesh, 5, 0, 1, 1, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 5;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 5, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      break;
    case T8_ECLASS_PYRAMID:
      vertices[0] = 1;
      vertices[1] = 3;
      vertices[2] = 0;
      vertices[3] = 2;
      vertices[4] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      vertices[0] = 0;
      vertices[1] = 2;
      vertices[2] = 4;
      vertices[3] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      vertices[0] = 1;
      vertices[1] = 0;
      vertices[2] = 5;
      vertices[3] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      t8_cmesh_set_join (cmesh, 0, 1, 3, 2, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 0, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 0, 2, 0, 0);
      break;
    default:
      break;
    }
  }
  if (do_bcast) {
    if (mpirank != 0) {
      cmesh = NULL;
    }
    cmesh = t8_cmesh_bcast (cmesh, 0, comm);
  }

  if (do_partition) {
    int                 mpirank, mpisize, mpiret;
    int                 first_tree, last_tree, num_trees;
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    num_trees = num_trees_for_hypercube[eclass];
    first_tree = (mpirank * num_trees) / mpisize;
    last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic (sc_MPI_Comm comm, int dim)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         tree_class;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  T8_ASSERT (dim == 1 || dim == 2 || dim == 3);
  t8_cmesh_init (&cmesh);
  switch (dim) {
  case 1:
    tree_class = T8_ECLASS_LINE;
    break;
  case 2:
    tree_class = T8_ECLASS_QUAD;
    break;
  case 3:
    tree_class = T8_ECLASS_HEX;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  t8_cmesh_set_tree_class (cmesh, 0, tree_class);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                              vertices, 1 << dim);
  t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
  if (dim > 1) {
    t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  }
  if (dim == 3) {
    t8_cmesh_set_join (cmesh, 0, 0, 4, 5, 0);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_bigmesh (t8_eclass_t eclass, int num_trees, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  int                 i;

  t8_cmesh_init (&cmesh);
  for (i = 0; i < num_trees; i++) {
    t8_cmesh_set_tree_class (cmesh, i, eclass);
    if (cmesh->dimension > 0) {
      /* We join each tree with its successor along faces 0 and 1
       * to get a nontrivial connectivity */
      t8_cmesh_set_join (cmesh, i, (i + 1) % num_trees, 0, 1, 0);
    }
  }

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* On each process, create a num_x by num_y (by num_z) brick connectivity and
 * make a cmesh connectivity from the disjoint union of those.
 * Example: 2 processors,
 * On the first  num_x = 1, num_y = 1
 * On the second num_x = 2, num_y = 1
 *                            _
 * connectivity on first:    |_|
 *
 *                           _ _
 * connectivity on second:  |_|_|
 *
 *                     _    _ _
 * Leads to the cmesh |_|  |_|_|
 * which is partitioned accordingly.
 */
t8_cmesh_t
t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x, t8_gloidx_t num_y,
                              t8_gloidx_t num_z, int x_periodic,
                              int y_periodic, int z_periodic,
                              sc_MPI_Comm comm)
{
  p4est_connectivity_t *my_brick = NULL;        /* pre-initialized to prevent compiler warning */
  p8est_connectivity_t *my_brick_3d = NULL;
  t8_cmesh_t          cmesh;
  t8_gloidx_t         num_trees, offset;
  int                 dim;

  T8_ASSERT (num_x >= 0 && num_y >= 0 && num_z >= 0);
  /* Set the dimension to 3 if num_z > 0 and 2 otherwise. */
  if (num_z > 0) {
    dim = 3;
  }
  else {
    dim = 2;
  }
  num_trees = num_x * num_y;
  if (dim == 3) {
    num_trees *= num_z;
  }
  /* Create a p4est brick connectivity on the process with
   * num_x times num_y elements */
  if (num_trees > 0) {
    if (dim == 2) {
      my_brick = p4est_connectivity_new_brick (num_x, num_y, x_periodic,
                                               y_periodic);
    }
    else {
      my_brick_3d = p8est_connectivity_new_brick (num_x, num_y, num_z,
                                                  x_periodic, y_periodic,
                                                  z_periodic);
    }
  }
  else {
    num_x = num_y = num_z = 0;
    num_trees = 0;
    if (dim == 2) {
      my_brick = p4est_connectivity_new (0, 0, 0, 0);
    }
    else {
      my_brick_3d = p8est_connectivity_new (0, 0, 0, 0, 0, 0);
    }
  }

  /* Calculate the x and y offset of trees */
  sc_MPI_Scan (&num_trees, &offset, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  offset -= num_trees;
  t8_debugf ("[H] offset = %li\n", offset);

  if (dim == 2) {
    cmesh = t8_cmesh_new_from_p4est_ext ((void *) my_brick,
                                         dim, comm, 1, offset + 1);
    p4est_connectivity_destroy (my_brick);
  }
  else {
    cmesh = t8_cmesh_new_from_p4est_ext ((void *) my_brick_3d,
                                         dim, comm, 1, offset + 1);
    p8est_connectivity_destroy (my_brick_3d);
  }
  return cmesh;
}
