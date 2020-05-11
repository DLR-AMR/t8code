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
#include <t8_data/t8_shmem.h>
#include <t8_vec.h>
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

/* For a committed cmesh check whether the entries of num_trees_per_eclass
 * and num_local_trees_per_eclass are valid.
 * Thus, num_local_trees_per_eclass[i] <= num_trees_per_eclass[i]
 * and the sum of the local trees must match cmesh->num_local_trees
 * and the sum of the global trees must match cmesh->num_trees.
 *
 * Returns true, if everything is fine.
 */
#ifdef T8_ENABLE_DEBUG
static int
t8_cmesh_check_trees_per_eclass (t8_cmesh_t cmesh)
{
  int                 ieclass;
  t8_gloidx_t         glo_trees = 0;
  t8_locidx_t         lo_trees = 0;
  int                 ret = 0;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
    ret = ret && cmesh->num_local_trees_per_eclass[ieclass] <=
      cmesh->num_trees_per_eclass[ieclass];
    lo_trees += cmesh->num_local_trees_per_eclass[ieclass];
    glo_trees += cmesh->num_trees_per_eclass[ieclass];
  }
  return !ret && lo_trees == cmesh->num_local_trees
    && glo_trees == cmesh->num_trees;
}
#endif

int
t8_cmesh_is_committed (t8_cmesh_t cmesh)
{
  static int          is_checking = 0;

  /* We run into a stackoverflow if routines that we call here,
   * also call t8_cmesh_is_committed.
   * We prevent this with the static variable is_checking.
   * This variable lives beyond one execution of t8_cmesh_is_committed.
   * We use it as a form of lock to prevent entering an infinite recursion.
   */
  if (!is_checking) {
    is_checking = 1;

    if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) &&
          cmesh->committed)) {
      is_checking = 0;
      return 0;
    }

#ifdef T8_ENABLE_DEBUG
    /* TODO: check more conditions that must always hold after commit */
    if ((!t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees)) ||
        (!t8_cmesh_no_negative_volume (cmesh))
        || (!t8_cmesh_check_trees_per_eclass (cmesh))) {
      is_checking = 0;
      return 0;
    }
#endif
    is_checking = 0;
  }
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
  cmesh->first_tree = -1;
  cmesh->first_tree_shared = -1;
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
  if (first_local_tree < 0) {
    /* the first tree is shared */
    cmesh->first_tree = -first_local_tree - 1;
    cmesh->first_tree_shared = 1;
  }
  else {
    /* The first tree is not shared */
    cmesh->first_tree = first_local_tree;
    cmesh->first_tree_shared = 0;
  }
  cmesh->num_local_trees = last_local_tree - cmesh->first_tree + 1;
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
    cmesh->first_tree_shared = -1;
    cmesh->num_local_trees = -1;
    cmesh->set_partition_level = -1;
  }
}

void
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, int element_level,
                                t8_scheme_cxx_t * ts)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (element_level >= -1);
  T8_ASSERT (ts != NULL);

  cmesh->set_partition = 1;
  cmesh->set_partition_level = element_level;
  cmesh->set_partition_scheme = ts;
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

int
t8_cmesh_treeid_is_local_tree (const t8_cmesh_t cmesh,
                               const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return 0 <= ltreeid && ltreeid < t8_cmesh_get_num_local_trees (cmesh);
}

int
t8_cmesh_treeid_is_ghost (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const t8_locidx_t   num_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t   num_ghosts = t8_cmesh_get_num_ghosts (cmesh);

  return num_trees <= ltreeid && ltreeid < num_trees + num_ghosts;
}

t8_locidx_t
t8_cmesh_ltreeid_to_ghostid (const t8_cmesh_t cmesh,
                             const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_ghost (cmesh, ltreeid));

  return ltreeid - t8_cmesh_get_num_local_trees (cmesh);
}

/* TODO: should get a gloidx?
 *       place after commit */
t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_locidx_t ltree_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id));

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
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, tree->treeid));
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

double             *
t8_cmesh_get_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid));

  return (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0,
                                            ltreeid);
}

void               *
t8_cmesh_get_attribute (t8_cmesh_t cmesh, int package_id, int key,
                        t8_locidx_t ltree_id)
{
  int                 is_ghost;

  T8_ASSERT (cmesh->committed);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id)
             || t8_cmesh_treeid_is_ghost (cmesh, ltree_id));
  is_ghost = t8_cmesh_treeid_is_ghost (cmesh, ltree_id);

  if (is_ghost) {
    ltree_id = t8_cmesh_ltreeid_to_ghostid (cmesh, ltree_id);
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
static              t8_locidx_t
t8_cmesh_tree_index (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  return cmesh->set_partition ? tree_id - cmesh->first_tree : tree_id;
}
#endif

void
t8_cmesh_set_dimension (t8_cmesh_t cmesh, int dim)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
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

/* Compute erg = v_1 . v_2
 * the 3D scalar product.
 */
static double
t8_cmesh_tree_vertices_dot (double *v_1, double *v_2)
{
  double              erg = 0;
  int                 i;

  for (i = 0; i < 3; i++) {
    erg += v_1[i] * v_2[i];
  }
  return erg;
}

/* Compute erg = v_1 x v_2
 * the 3D cross product.
 */
static void
t8_cmesh_tree_vertices_cross (double *v_1, double *v_2, double *erg)
{
  int                 i;

  for (i = 0; i < 3; i++) {
    erg[i] = v_1[(i + 1) % 3] * v_2[(i + 2) % 3]
      - v_1[(i + 2) % 3] * v_2[(i + 1) % 3];
  }
}

/* Given a set of vertex coordinates for a tree of a given eclass.
 * Query whether the geometric volume of the tree with this coordinates
 * would be negative.
 * Returns true if a tree of the given eclass with the given vertex
 * coordinates does have negative volume.
 */
int
t8_cmesh_tree_vertices_negative_volume (t8_eclass_t eclass,
                                        double *vertices, int num_vertices)
{
  double              v_1[3], v_2[3], v_j[3], cross[3], sc_prod;
  int                 i, j;

  T8_ASSERT (num_vertices == t8_eclass_num_vertices[eclass]);

  if (t8_eclass_to_dimension[eclass] <= 2) {
    /* Only three dimensional eclass do have a volume */
    return 0;
  }

  T8_ASSERT (eclass == T8_ECLASS_TET || eclass == T8_ECLASS_HEX
             || eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_PYRAMID);
  T8_ASSERT (num_vertices >= 4);

  /*
   *      6 ______  7  For Hexes and pyramids, if the vertex 4 is below the 0-1-2-3 plane,
   *       /|     /     the volume is negative. This is the case if and only if
   *    4 /_____5/|     the scalar product of v_4 with the cross product of v_1 and v_2 is
   *      | | _ |_|     smaller 0:
   *      | 2   | / 3   < v_4, v_1 x v_2 > < 0
   *      |/____|/
   *     0      1
   *
   *
   *    For tets/prisms, if the vertex 3 is below/above the 0-1-2 plane, the volume
   *    is negative. This is the case if and only if
   *    the scalar product of v_3 with the cross product of v_1 and v_2 is
   *    greater 0:
   *
   *    < v_3, v_1 x v_2 > > 0
   *
   */

  /* build the vectors v_i as vertices_i - vertices_0 */

  if (eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM) {
    /* In the tet/prism case, the third vector is v_3 */
    j = 3;
  }
  else {
    /* For pyramids and Hexes, the third vector is v_4 */
    j = 4;
  }
  for (i = 0; i < 3; i++) {
    v_1[i] = vertices[3 + i] - vertices[i];
    v_2[i] = vertices[6 + i] - vertices[i];
    v_j[i] = vertices[3 * j + i] - vertices[i];
  }

  /* compute cross = v_1 x v_2 */
  t8_cmesh_tree_vertices_cross (v_1, v_2, cross);
  /* Compute sc_prod = <v_j, cross> */
  sc_prod = t8_cmesh_tree_vertices_dot (v_j, cross);

  T8_ASSERT (sc_prod != 0);
  return eclass == T8_ECLASS_TET ? sc_prod > 0 : sc_prod < 0;
}

#ifdef T8_ENABLE_DEBUG
/* After a cmesh is committed, check whether all trees in a cmesh do have positive volume.
 * Returns true if all trees have positive volume.
 */
int
t8_cmesh_no_negative_volume (t8_cmesh_t cmesh)
{
  t8_locidx_t         itree;
  double             *vertices;
  t8_eclass_t         eclass;
  int                 ret, res = 0;

  if (cmesh == NULL) {
    return 0;
  }
  /* Iterate over all trees, get their vertices and check the volume */
  for (itree = 0; itree < cmesh->num_local_trees; itree++) {
    vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
    ret = 1;
    if (vertices != NULL) {
      /* Vertices are set */
      eclass = t8_cmesh_get_tree_class (cmesh, itree);
      ret = t8_cmesh_tree_vertices_negative_volume (eclass, vertices,
                                                    t8_eclass_num_vertices
                                                    [eclass]);
      if (ret) {
        t8_debugf ("Detected negative volume in tree %li\n", (long) itree);
      }
      res |= ret;               /* res is true if one ret value is true */
    }
  }
  return !res;
}
#endif

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
  is_equal = is_equal || memcmp (cmesh_a->num_local_trees_per_eclass,
                                 cmesh_b->num_local_trees_per_eclass,
                                 T8_ECLASS_COUNT * sizeof (t8_locidx_t));

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
  t8_cmesh_t          cmesh_out;

  struct
  {
    t8_cmesh_struct_t   cmesh;
    t8_gloidx_t         num_trees_per_eclass[T8_ECLASS_COUNT];
    size_t              stash_elem_counts[3];
    int                 pre_commit;     /* True, if cmesh on root is not committed yet. */
#ifdef T8_ENABLE_DEBUG
    sc_MPI_Comm         comm;
#endif
  } meta_info;

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
    memcpy (&meta_info.cmesh, cmesh_in, sizeof (*cmesh_in));
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      meta_info.num_trees_per_eclass[iclass] =
        cmesh_in->num_trees_per_eclass[iclass];
      T8_ASSERT (cmesh_in->num_local_trees_per_eclass[iclass] ==
                 cmesh_in->num_trees_per_eclass[iclass]);
    }
    if (t8_cmesh_is_committed (cmesh_in)) {
      meta_info.pre_commit = 0;
    }
    else {
      meta_info.pre_commit = 1;
      meta_info.stash_elem_counts[0] = cmesh_in->stash->attributes.elem_count;
      meta_info.stash_elem_counts[1] = cmesh_in->stash->classes.elem_count;
      meta_info.stash_elem_counts[2] = cmesh_in->stash->joinfaces.elem_count;
    }
#ifdef T8_ENABLE_DEBUG
    meta_info.comm = comm;
#endif
    /* Root returns the input cmesh */
    cmesh_out = cmesh_in;
  }
  /* TODO: we could optimize this by using IBcast */
  mpiret = sc_MPI_Bcast (&meta_info, sizeof (meta_info), sc_MPI_BYTE, root,
                         comm);
  SC_CHECK_MPI (mpiret);

  /* If not root store information in new cmesh and allocate memory for arrays. */
  if (mpirank != root) {
    t8_cmesh_init (&cmesh_out);
    cmesh_out->dimension = meta_info.cmesh.dimension;
    cmesh_out->face_knowledge = meta_info.cmesh.face_knowledge;
    cmesh_out->set_partition = meta_info.cmesh.set_partition;
    cmesh_out->set_partition_level = meta_info.cmesh.set_partition_level;
    cmesh_out->set_refine_level = meta_info.cmesh.set_refine_level;
    cmesh_out->num_trees = meta_info.cmesh.num_trees;
    cmesh_out->num_local_trees = cmesh_out->num_trees;
    cmesh_out->first_tree = 0;
    cmesh_out->first_tree_shared = 0;
    cmesh_out->num_ghosts = 0;
    T8_ASSERT (cmesh_out->set_partition == 0);
    if (meta_info.cmesh.profile != NULL) {
      t8_cmesh_set_profiling (cmesh_in, 1);
    }
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      cmesh_out->num_trees_per_eclass[iclass] =
        meta_info.num_trees_per_eclass[iclass];
      cmesh_out->num_local_trees_per_eclass[iclass] =
        meta_info.num_trees_per_eclass[iclass];
    }
#ifdef T8_ENABLE_DEBUG
    T8_ASSERT (meta_info.comm == comm);
#endif
  }
  if (meta_info.pre_commit) {
    /* broadcast all the stashed information about trees/neighbors/attributes */
    t8_stash_bcast (cmesh_out->stash, root, comm,
                    meta_info.stash_elem_counts);
  }
  else {
    /* broadcast the stored information about the trees */
    t8_cmesh_trees_bcast (cmesh_out, root, comm);
    if (mpirank != root) {
      /* destroy stash and set to committed */
      t8_stash_destroy (&cmesh_out->stash);
      cmesh_out->committed = 1;
    }
  }

  cmesh_out->mpirank = mpirank;
  cmesh_out->mpisize = mpisize;
  /* Final checks */
#ifdef T8_ENABLE_DEBUG
  if (!meta_info.pre_commit) {
    T8_ASSERT (t8_cmesh_is_committed (cmesh_out));
    T8_ASSERT (t8_cmesh_comm_is_valid (cmesh_out, comm));
  }
#endif
  return cmesh_out;
}

#ifdef T8_WITH_METIS
void
t8_cmesh_reorder (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 mpisize, mpiret;
  idx_t               idx_mpisize;
  idx_t               ncon = 1, elemens;
  idx_t               volume, *partition, ipart, newpart;
  int                 num_faces, iface, count_face;
  idx_t              *xadj, *adjncy;
  int                 success;
  t8_locidx_t        *new_number, itree, *tree_per_part_off, *tree_per_part;
  t8_locidx_t        *face_neighbor;
  t8_locidx_t         neigh_id;
  t8_ctree_t          tree;

  /* cmesh must be commited and not partitioned */
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (!cmesh->set_partition);

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  idx_mpisize = mpisize;
  SC_CHECK_MPI (mpiret);

  elemens = cmesh->num_trees;
  T8_ASSERT ((t8_locidx_t) elemens == cmesh->num_trees);

  /* Count the number of tree-to-tree connections via a face */
  num_faces = 0;
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor,
                                        NULL);
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      if (face_neighbor[iface] >= 0)
        num_faces++;
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
    tree = t8_cmesh_get_tree (cmesh, itree);
    xadj[itree + 1] = xadj[itree];
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      if (face_neighbor[iface] >= 0) {
        adjncy[count_face++] = face_neighbor[iface];
        xadj[itree + 1]++;
      }
    }
  }

  /* partition stores the new partition number for each element */
  partition = T8_ALLOC (idx_t, elemens);
  /* partition the elements in mpisize many partitions */
  success =
    METIS_PartGraphRecursive (&elemens, &ncon, xadj, adjncy, NULL, NULL, NULL,
                              &idx_mpisize, NULL, NULL, NULL, &volume,
                              partition);
  T8_ASSERT (success == METIS_OK);
  /* memory to store the new treeid of a tree */
  new_number = T8_ALLOC (t8_locidx_t, cmesh->num_trees);
  /* Store the number of trees pinter partition */
  tree_per_part = T8_ALLOC_ZERO (t8_locidx_t, mpisize);
  /* Store the treeid offset of each partition. */
  tree_per_part_off = T8_ALLOC_ZERO (t8_locidx_t, mpisize + 1);
  tree_per_part_off[0] = 0;
  /* compute tree_per_part and prepare tree_per_part_off */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree_per_part[partition[itree]]++;
    tree_per_part_off[partition[itree] + 1]++;
  }
  /* compute tree_per_part_off */
  for (ipart = 1; ipart <= mpisize; ipart++) {
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
  /* Set for each tree its new treeid and the new ids of its neighbors */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor,
                                        NULL);
    tree->treeid = new_number[itree];
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      neigh_id = face_neighbor[iface];
      if (neigh_id >= 0) {
        face_neighbor[iface] = new_number[neigh_id];
      }
    }
  }
  T8_FREE (partition);
  T8_FREE (xadj);
  T8_FREE (adjncy);
  T8_FREE (new_number);
  T8_FREE (tree_per_part);
  T8_FREE (tree_per_part_off);
}
#endif

int
t8_cmesh_is_partitioned (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->set_partition != 0;
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

int
t8_cmesh_tree_face_is_boundary (const t8_cmesh_t cmesh,
                                const t8_locidx_t ltreeid, const int face)
{
  int8_t             *ttf;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  if (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid)) {
    /* The local tree id belongs to a tree */
    t8_locidx_t        *face_neighbor;
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, ltreeid, &face_neighbor,
                                        &ttf);

    if (face_neighbor[face] == ltreeid && ttf[face] == face) {
      /* The tree is connected to itself at the same face.
       * Thus this is a domain boundary */
      return 1;
    }
  }
  else {
    /* The local tree id belongs to a ghost */
    T8_ASSERT (t8_cmesh_treeid_is_ghost (cmesh, ltreeid));

    t8_gloidx_t        *face_neighbor;
    const t8_locidx_t   lghostid =
      t8_cmesh_ltreeid_to_ghostid (cmesh, ltreeid);
    (void) t8_cmesh_trees_get_ghost_ext (cmesh->trees, lghostid,
                                         &face_neighbor, &ttf);

    if (face_neighbor[face] == t8_cmesh_get_global_id (cmesh, ltreeid)
        && ttf[face] == face) {
      /* The ghost is connected to itself at the same face.
       * Thus this is a domain boundary */
      return 1;
    }
  }

  return 0;
}

t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_locidx_t ltree_id)
{
  t8_ctree_t          tree;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id));

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
  /* Check that we do not get wrong numbers when converting to locidx */
  T8_ASSERT ((t8_locidx_t) temp_local_id == temp_local_id);
  if (t8_cmesh_treeid_is_local_tree (cmesh, temp_local_id)) {
    /* The tree is a local tree */
    return temp_local_id;
  }
  else {
    /* The tree may be a ghost tree */
    return t8_cmesh_trees_get_ghost_local_id (cmesh->trees, global_id);
  }
}

/* Given a local tree id and a face number, get information about the face neighbor tree.
 * \param [in]      cmesh     The cmesh to be considered.
 * \param [in]      ltreeid   The local id of a tree or a ghost.
 * \param [in]      face      A face number of the tree/ghost.
 * \param [out]     dual_face If not NULL, the face number of the neighbor tree at this connection.
 * \param [out]     orientation If not NULL, the face orientation of the connection.
 * \return                    If non-negative: The local id of the neighbor tree or ghost.
 *                            If negative: There is no neighbor across this face. \a dual_face and
 *                            \a orientation remain unchanged.
 * \note If \a ltreeid is a ghost and it has a neighbor which is neither a local tree or ghost,
 *       then the return value will be negative.
 *       This, a negative return value does not necessarily mean that this is a domain boundary.
 *       To find out whether a tree is a domain boundary or not \see t8_cmesh_tree_face_is_boundary.
 */
t8_locidx_t
t8_cmesh_get_face_neighbor (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid,
                            const int face, int *dual_face, int *orientation)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid)
             || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  const int           is_ghost = t8_cmesh_treeid_is_ghost (cmesh, ltreeid);
  int8_t              ttf;
  t8_locidx_t         face_neigh;
  int                 dual_face_temp, orientation_temp;

  /* If this is a domain boundary, return -1 */
  if (t8_cmesh_tree_face_is_boundary (cmesh, ltreeid, face)) {
    return -1;
  }

  if (!is_ghost) {
    /* The local tree id belongs to a local tree (not a ghost) */
    /* Get the tree */
    const t8_ctree_t    tree = t8_cmesh_get_tree (cmesh, ltreeid);

#ifdef T8_ENABLE_DEBUG
    /* Get the eclass */
    t8_eclass_t         eclass = tree->eclass;
    /* Check that face is valid */
    T8_ASSERT (0 <= face && face < t8_eclass_num_faces[eclass]);
#endif

    /* Get the local id of the face neighbor */
    face_neigh = t8_cmesh_trees_get_face_neighbor_ext (tree, face, &ttf);
  }
  else {
    /* The local tree id belongs to a ghost */
    const t8_locidx_t   lghostid =
      ltreeid - t8_cmesh_get_num_local_trees (cmesh);
    /* Get the ghost */
    const t8_cghost_t   ghost =
      t8_cmesh_trees_get_ghost (cmesh->trees, lghostid);

    t8_gloidx_t         global_face_neigh;

#ifdef T8_ENABLE_DEBUG
    /* Get the eclass */
    t8_eclass_t         eclass = ghost->eclass;
    /* Check that face is valid */
    T8_ASSERT (0 <= face && face < t8_eclass_num_faces[eclass]);
#endif

    /* Get the global id of the face neighbor */
    global_face_neigh =
      t8_cmesh_trees_get_ghost_face_neighbor_ext (ghost, face, &ttf);
    /* Convert it into a local id */
    face_neigh = t8_cmesh_get_local_id (cmesh, global_face_neigh);

    /* TODO: Check whether this face is a boundary face */
    if (face_neigh < 0) {
      /* The neighbor is not local, return -1 */
      return -1;
    }
  }

  /* Decode the ttf information to get the orientation and the dual face */
  t8_cmesh_tree_to_face_decode (cmesh->dimension, ttf, &dual_face_temp,
                                &orientation_temp);
  if (dual_face != NULL) {
    *dual_face = dual_face_temp;
  }
  if (orientation != NULL) {
    *orientation = orientation_temp;
  }
  /* Return the face neighbor */
  return face_neigh;
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

  /* unref the refine scheme (if set) */
  if (cmesh->set_refine_scheme != NULL) {
    t8_scheme_cxx_unref (&cmesh->set_refine_scheme);
  }

  /* unref the partition scheme (if set) */
  if (cmesh->set_partition_scheme != NULL) {
    t8_scheme_cxx_unref (&cmesh->set_partition_scheme);
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
      t8_debugf ("Generating partitioned cmesh from connectivity\n"
                 "Has %li global and %li local trees.\n", num_trees,
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

static              t8_cmesh_t
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

/* Compute y = ax + b on an array of doubles, interpreting
 * each 3 as one vector x */
static void
t8_cmesh_coords_axb (const double *coords_in, double *coords_out,
                     int num_vertices, double alpha, const double b[3])
{
  int                 i;

  for (i = 0; i < num_vertices; i++) {
    t8_vec_axpyz (coords_in + i * 3, b, coords_out + i * 3, alpha);
  }
}

t8_cmesh_t
t8_cmesh_new_hypercube_hybrid (int dim, sc_MPI_Comm comm, int do_partition,
                               int periodic)
{
  int                 i;
  t8_cmesh_t          cmesh;
  t8_topidx_t         vertices[8];
  double              vertices_coords_temp[24];
  double              attr_vertices[24];
  double              null_vec[3] = { 0, 0, 0 };
  double              shift[7][3] = { {0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.5},
  {0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0.5}, {0, 0.5, 0.5}
  };
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

  t8_cmesh_init (&cmesh);
  /* This cmesh consists of 6 tets, 6 prisms and 3 hexes */
  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TET);
  }
  for (i = 6; i < 12; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }
  for (i = 12; i < 16; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_HEX);
  }

    /************************************/
  /*  The tetrahedra                  */
    /************************************/
  /* We place the tetrahedra at the origin of the unit cube.
   * They are essentially the tetrahedral hypercube scaled by 0.5 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       null_vec);
  t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 5, 0, 2, 1, 0);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 5;
  vertices[3] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                              attr_vertices, 4);
  vertices[1] = 3;
  vertices[2] = 1;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                              attr_vertices, 4);
  vertices[1] = 2;
  vertices[2] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                              attr_vertices, 4);
  vertices[1] = 6;
  vertices[2] = 2;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0,
                              attr_vertices, 4);
  vertices[1] = 4;
  vertices[2] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0,
                              attr_vertices, 4);
  vertices[1] = 5;
  vertices[2] = 4;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 5, t8_get_package_id (), 0,
                              attr_vertices, 4);

    /************************************/
  /*     The prisms                   */
    /************************************/
  /* We place the prism to the left, right, and top of the tetrahedra.
   * They are essentially the prism hypercube scaled by 0.5 and
   * shifted in 3 different direction. */
  /* trees 6 and 7 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[0]);
  vertices[0] = 0;
  vertices[1] = 6;
  vertices[2] = 4;
  vertices[3] = 1;
  vertices[4] = 7;
  vertices[5] = 5;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 6, t8_get_package_id (), 0,
                              attr_vertices, 6);
  vertices[1] = 2;
  vertices[2] = 6;
  vertices[4] = 3;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 7, t8_get_package_id (), 0,
                              attr_vertices, 6);

  t8_cmesh_set_join (cmesh, 6, 7, 2, 1, 0);
  /* trees 8 and 9 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[1]);
  vertices[0] = 0;
  vertices[1] = 5;
  vertices[2] = 1;
  vertices[3] = 2;
  vertices[4] = 7;
  vertices[5] = 3;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 8, t8_get_package_id (), 0,
                              attr_vertices, 6);
  vertices[1] = 4;
  vertices[2] = 5;
  vertices[4] = 6;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 9, t8_get_package_id (), 0,
                              attr_vertices, 6);
  t8_cmesh_set_join (cmesh, 8, 9, 2, 1, 0);
  /* trees 10 an 11 */
  t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                       shift[2]);
  vertices[0] = 0;
  vertices[1] = 1;
  vertices[2] = 3;
  vertices[3] = 4;
  vertices[4] = 5;
  vertices[5] = 7;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 10, t8_get_package_id (), 0,
                              attr_vertices, 6);
  vertices[1] = 3;
  vertices[2] = 2;
  vertices[4] = 7;
  vertices[5] = 6;
  t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                 vertices_coords_temp,
                                                 attr_vertices, 6);
  t8_cmesh_set_tree_vertices (cmesh, 11, t8_get_package_id (), 0,
                              attr_vertices, 6);
  t8_cmesh_set_join (cmesh, 10, 11, 1, 2, 0);

  /* Connect prisms and tets */
  t8_cmesh_set_join (cmesh, 0, 6, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 1, 7, 0, 3, 1);
  t8_cmesh_set_join (cmesh, 2, 8, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 3, 9, 0, 3, 1);
  t8_cmesh_set_join (cmesh, 4, 11, 0, 3, 0);
  t8_cmesh_set_join (cmesh, 5, 10, 0, 3, 1);

  /************************************/
  /*  The hexahedra                   */
  /************************************/

  for (i = 0; i < 8; i++) {
    vertices[i] = i;
  }

  for (i = 0; i < 4; i++) {
    t8_cmesh_coords_axb (vertices_coords, vertices_coords_temp, 8, 0.5,
                         shift[3 + i]);
    t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                   vertices_coords_temp,
                                                   attr_vertices, 8);
    t8_cmesh_set_tree_vertices (cmesh, 12 + i, t8_get_package_id (), 0,
                                attr_vertices, 8);
  }
  /* Join the hexes */
  t8_cmesh_set_join (cmesh, 12, 14, 5, 4, 0);
  t8_cmesh_set_join (cmesh, 13, 14, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 14, 15, 0, 1, 0);

  /* Join the prisms and hexes */
  t8_cmesh_set_join (cmesh, 6, 13, 0, 4, 1);
  t8_cmesh_set_join (cmesh, 7, 12, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 8, 12, 0, 0, 1);
  t8_cmesh_set_join (cmesh, 9, 15, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 10, 13, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 11, 15, 0, 2, 1);

  if (periodic) {
    /* Connect the sides of the cube to make it periodic */
    /* tets to prisms */
    t8_cmesh_set_join (cmesh, 0, 8, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 5, 9, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 3, 7, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 4, 6, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 1, 10, 3, 4, 0);
    t8_cmesh_set_join (cmesh, 2, 11, 3, 4, 0);
    /* prism to hex */
    t8_cmesh_set_join (cmesh, 6, 12, 1, 3, 0);
    t8_cmesh_set_join (cmesh, 9, 12, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 7, 13, 2, 5, 0);
    t8_cmesh_set_join (cmesh, 11, 13, 1, 1, 0);
    t8_cmesh_set_join (cmesh, 8, 15, 1, 5, 0);
    t8_cmesh_set_join (cmesh, 10, 15, 2, 3, 0);
    /* hex to hex */
    t8_cmesh_set_join (cmesh, 12, 14, 4, 5, 0);
    t8_cmesh_set_join (cmesh, 13, 14, 2, 3, 0);
    t8_cmesh_set_join (cmesh, 14, 15, 1, 0, 0);

  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* The unit cube is constructed from trees of the same eclass.
 * For triangles the square is divided along the (0,0) -- (1,1) diagonal.
 * For prisms the front (y=0) and back (y=1) face are divided into triangles
 * as above.
 */
/* TODO: upgrade with int x,y,z for periodic faces */
t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_bcast,
                        int do_partition, int periodic)
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

  SC_CHECK_ABORT (eclass != T8_ECLASS_PYRAMID || !periodic,
                  "The pyramid cube mesh cannot be periodic.");

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
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 4, 5, 0);
      }
    case T8_ECLASS_QUAD:
      vertices[3] = 3;
      vertices[2] = 2;
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
      }
    case T8_ECLASS_LINE:
      vertices[1] = 1;
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
      }
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
      vertices[2] = 3;
      vertices[3] = 4;
      vertices[4] = 5;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      vertices[1] = 3;
      vertices[2] = 2;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
        t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
        t8_cmesh_set_join (cmesh, 0, 0, 3, 4, 0);
        t8_cmesh_set_join (cmesh, 1, 1, 3, 4, 0);
      }
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
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
        t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
      }
      break;
    case T8_ECLASS_TET:
      t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 3, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 3, 4, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
      t8_cmesh_set_join (cmesh, 5, 0, 2, 1, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 3;
      vertices[2] = 1;
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
      vertices[1] = 6;
      vertices[2] = 2;
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
      vertices[1] = 5;
      vertices[2] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 5, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      if (periodic) {
        t8_cmesh_set_join (cmesh, 0, 4, 0, 3, 0);
        t8_cmesh_set_join (cmesh, 1, 3, 0, 3, 2);

        t8_cmesh_set_join (cmesh, 0, 2, 3, 0, 0);
        t8_cmesh_set_join (cmesh, 3, 5, 0, 3, 2);

        t8_cmesh_set_join (cmesh, 1, 5, 3, 0, 2);
        t8_cmesh_set_join (cmesh, 2, 4, 3, 0, 0);
      }
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
t8_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  double              vertices[12] = {
    0, 0, 0,
    0.2, 0, 0,
    0.6, 0, 0,
    1, 0, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_LINE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices + 3,
                              2);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices + 6,
                              2);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 0, 1, 0, 0);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_tri (sc_MPI_Comm comm)
{
  double              vertices[18] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 0,
    1, 1, 0,
    0, 1, 0
  };
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices + 9,
                              3);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 1, 0, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 1);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic_hybrid (sc_MPI_Comm comm)
{
  double              vertices[60] = {  /* Just all vertices of all trees. partly duplicated */
    0, 0, 0,                    /* tree 0, triangle */
    0.5, 0, 0,
    0.5, 0.5, 0,
    0, 0, 0,                    /* tree 1, triangle */
    0.5, 0.5, 0,
    0, 0.5, 0,
    0.5, 0, 0,                  /* tree 2, quad */
    1, 0, 0,
    0.5, 0.5, 0,
    1, 0.5, 0,
    0, 0.5, 0,                  /* tree 3, quad */
    0.5, 0.5, 0,
    0, 1, 0,
    0.5, 1, 0,
    0.5, 0.5, 0,                /* tree 4, triangle */
    1, 0.5, 0,
    1, 1, 0,
    0.5, 0.5, 0,                /* tree 5, triangle */
    1, 1, 0,
    0.5, 1, 0
  };
  t8_cmesh_t          cmesh;

  /*
   *  This is how the cmesh looks like. The numbers are the tree numbers:
   *
   *   +---+---+
   *   |   |5 /|
   *   | 3 | / |
   *   |   |/ 4|
   *   +---+---+
   *   |1 /|   |
   *   | / | 2 |
   *   |/0 |   |
   *   +---+---+
   */

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);

  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices + 9,
                              3);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                              vertices + 18, 4);
  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0,
                              vertices + 30, 4);
  t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0,
                              vertices + 42, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, t8_get_package_id (), 0,
                              vertices + 51, 3);

  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 2, 3, 0);

  t8_cmesh_set_join (cmesh, 1, 3, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);

  t8_cmesh_set_join (cmesh, 2, 4, 3, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 5, 2, 0, 1);

  t8_cmesh_set_join (cmesh, 3, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0);

  t8_cmesh_set_join (cmesh, 4, 5, 1, 2, 0);

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

t8_cmesh_t
t8_cmesh_new_line_zigzag (sc_MPI_Comm comm)
{
  int                 i;
  double              vertices[18] = { 1, 2, 0,
    2, 4, 1,
    1, 1, 2,
    2, 4, 1,
    1, 1, 2,
    3, 2, 5
  };
  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);
  for (i = 0; i < 3; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_LINE);
  }
  /*tree_num is joined with tree_num at face_num and face_num with orientation_num */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 0, 0, 0);

  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 2);
  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices + 6,
                              2);
  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                              vertices + 12, 2);

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism_cake (sc_MPI_Comm comm, int num_of_prisms)
{
  int                 i, j;
  /*num_of_prisms Prism a 6 vertices a 3 coords */
  /* TODO: This seems too be a lot of memory, can we also get by with only
     6 * 3 doubles? */
  double             *vertices = T8_ALLOC (double, num_of_prisms * 6 * 3);
  t8_cmesh_t          cmesh;
  double              degrees = 360. / num_of_prisms;

  T8_ASSERT (num_of_prisms > 2);

  for (i = 0; i < num_of_prisms; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * degrees * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * degrees * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] =
          cos ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] =
          sin ((i * degrees + degrees) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  t8_cmesh_init (&cmesh);
  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }

  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_join (cmesh, i, (i == (num_of_prisms - 1) ? 0 : i + 1), 1, 2,
                       0);
  }
  for (i = 0; i < num_of_prisms; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, t8_get_package_id (), 0,
                                vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  T8_FREE (vertices);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism_deformed (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[18] = { -1, -0.5, 0.25,
    1, 0, 0,
    1, 1, 0,
    0, 0, 0.75,
    1.25, 0, 1,
    2, 2, 2
  };
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/*rotates counterclockwise*/
static void
prism_rotate (double vertices[18], int rotation)
{
  double              helper[3] = { vertices[6], vertices[7], vertices[8] };
  int                 i, j;
  T8_ASSERT (3 > rotation && rotation > 0);
  for (i = 0; i < rotation; i++) {
    for (j = 8; j >= 0; j--) {
      vertices[j] = j >= 3 ? vertices[j - 3] : helper[j];
    }
    for (j = 0; j < 3; j++) {
      helper[j] = vertices[6 + j];
    }
  }
  for (i = 0; i < 3; i++) {
    helper[i] = vertices[15 + i];
  }
  for (i = 0; i < rotation; i++) {
    for (j = 17; j >= 9; j--) {
      vertices[j] = j >= 12 ? vertices[j - 3] : helper[j - 9];
    }
    for (j = 0; j < 3; j++) {
      helper[j] = vertices[15 + j];
    }
  }
}

t8_cmesh_t
t8_cmesh_new_prism_cake_funny_oriented (sc_MPI_Comm comm)
{
  int                 i, j;
  /*6 Prism a 6 vertices a 3 coords */
  double              vertices[108];
  t8_cmesh_t          cmesh;

  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  prism_rotate (vertices + 18, 2);
  prism_rotate (vertices + 36, 1);
  prism_rotate (vertices + 54, 1);
  prism_rotate (vertices + 72, 1);
  prism_rotate (vertices + 90, 2);

  t8_cmesh_init (&cmesh);
  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }

  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 3);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 0, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 0, 1, 1, 0);

  for (i = 0; i < 6; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, t8_get_package_id (), 0,
                                vertices + i * 18, 6);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism_geometry (sc_MPI_Comm comm)
{
  int                 i, j;
  /*8 Prism a 6 vertices a 3 coords */
  double              vertices[144];
  t8_cmesh_t          cmesh;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[i * 6 * 3 + j * 3] = 0;
        vertices[i * 6 * 3 + j * 3 + 1] = 0;
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 3 ? 1 : 0);
      }
      else if (j == 1 || j == 4) {
        vertices[i * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 4 ? 1 : 0);
      }
      else if (j == 2 || j == 5) {
        vertices[i * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 1] = sin ((i * 60 + 60) * M_PI / 180);
        vertices[i * 6 * 3 + j * 3 + 2] = (j == 5 ? 1 : 0);
      }
    }
  }
  for (i = 2; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      /*Get the edges at the unit circle */
      if (j == 0 || j == 3) {
        vertices[(i + 1) * 6 * 3 + j * 3] = 0;
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] = 0;
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 3 ? 2 : 1);
      }
      else if (j == 1 || j == 4) {
        vertices[(i + 1) * 6 * 3 + j * 3] = cos (i * 60 * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] = sin (i * 60 * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 4 ? 2 : 1);
      }
      else if (j == 2 || j == 5) {
        vertices[(i + 1) * 6 * 3 + j * 3] = cos ((i * 60 + 60) * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 1] =
          sin ((i * 60 + 60) * M_PI / 180);
        vertices[(i + 1) * 6 * 3 + j * 3 + 2] = (j == 5 ? 2 : 1);
      }
    }
  }
  vertices[126] = cos (300 * M_PI / 180);
  vertices[127] = sin (300 * M_PI / 180);
  vertices[128] = 1;
  vertices[129] = 1;
  vertices[130] = 0;
  vertices[131] = 1;
  vertices[132] = cos (300 * M_PI / 180) + 1;
  vertices[133] = sin (300 * M_PI / 180);
  vertices[134] = 1;
  vertices[135] = cos (300 * M_PI / 180);
  vertices[136] = sin (300 * M_PI / 180);
  vertices[137] = 2;
  vertices[138] = 1;
  vertices[139] = 0;
  vertices[140] = 2;
  vertices[141] = cos (300 * M_PI / 180) + 1;
  vertices[142] = sin (300 * M_PI / 180);
  vertices[143] = 2;
  prism_rotate (vertices + 18, 2);
  prism_rotate (vertices + 36, 1);
  prism_rotate (vertices + 72, 2);

  t8_cmesh_init (&cmesh);
  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_PRISM);
  }
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 3);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 4, 3, 3);
  t8_cmesh_set_join (cmesh, 3, 4, 2, 0, 3);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 5, 6, 2, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 0, 1, 0);

  for (i = 0; i < 8; i++) {
    t8_cmesh_set_tree_vertices (cmesh, i, t8_get_package_id (), 0,
                                vertices + i * 18, 6);
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

static void
t8_cmesh_translate_coordinates (const double *coords_in, double *coords_out,
                                int num_vertices, double translate[3])
{
  int                 i;

  for (i = 0; i < num_vertices; i++) {
    coords_out[3 * i] = coords_in[3 * i] + translate[0];
    coords_out[3 * i + 1] = coords_in[3 * i + 1] + translate[1];
    coords_out[3 * i + 2] = coords_in[3 * i + 2] + translate[2];
  }
}

/* Construct a tetrahedral cmesh that has all possible face to face
 * connections and orientations. */
t8_cmesh_t
t8_cmesh_new_tet_orientation_test (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  int                 i;

  double              vertices_coords[12] = {
    0, 0, 0,
    1, 0, 0,
    1, 0, 1,
    1, 1, 1
  };
  double              translated_coords[12];
  double              translate[3] = { 1, 0, 0 };
  const t8_gloidx_t   num_trees = 24;

  t8_cmesh_init (&cmesh);
  /* A tet has 4 faces and each face connection has 3 possible orientations,
   * we thus have (4+3+2+1)*3 = 30 possible face-to-face combinations.
   * We use a cmesh of 24 tetrahedron trees. */
  for (i = 0; i < num_trees; i++) {
    t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_TET);
  }
  /* face combinations:
   *  0 - 0 0 - 1 0 - 2 0 - 3
   *  1 - 1 1 - 2 1 - 3
   *  2 - 2 2 - 3
   *  3 - 3
   */
  /* i iterates over the orientations */
  for (i = 0; i < 3; i++) {
    /* Face 0 with face k */
    /* For trees 0 -> 1, 2 -> 3, 4 -> 5, ..., 22 -> 23 */
    t8_cmesh_set_join (cmesh, 8 * i, 8 * i + 1, 0, 0, i);
    t8_cmesh_set_join (cmesh, 8 * i + 2, 8 * i + 3, 0, 1, i);
    t8_cmesh_set_join (cmesh, 8 * i + 4, 8 * i + 5, 0, 2, i);
    t8_cmesh_set_join (cmesh, 8 * i + 6, 8 * i + 7, 0, 3, i);
    /* Each tree with an even number has face 0 connected */
    /* Trees 1,  9, 17 face 0
     * Trees 3, 11, 19 face 1
     * Trees 5, 13, 21 face 2
     * Trees 7, 15, 23 face 3 */

    /* Face 1 with face k */
    /* Connect face 1 of trees 0 -> 1, 2 -> 3, ..., 16->17 */
    t8_cmesh_set_join (cmesh, 6 * i, 6 * i + 1, 1, 1, i);
    t8_cmesh_set_join (cmesh, 6 * i + 2, 6 * i + 3, 1, 2, i);
    t8_cmesh_set_join (cmesh, 6 * i + 4, 6 * i + 5, 1, 3, i);
    /* Each tree with even number up to 16 has face 1 connected. */
    /* Trees 1,  7, 13 face 1
     * Trees 3,  9, 15 face 2
     * Trees 5, 11, 17 face 3
     */

    /* Face 2 with face k */
    /* Connect face 2 of trees 0 -> 1, 2 -> 3,...,10 -> 11 */
    t8_cmesh_set_join (cmesh, 4 * i, 4 * i + 12, 2, 2, i);
    t8_cmesh_set_join (cmesh, 4 * i + 2, 4 * i + 6, 2, 3, i);
    /* Each tree with even number up to 10 has face 2 connected */
    /* Trees  12, 16, 20 face 2
     * Trees   6, 10, 14 face 3
     */

    /* Face 3 with face k */
    /* Connect face 3 of tree 0 -> 1, 2 -> 3, 4 -> 5 */
    t8_cmesh_set_join (cmesh, 2 * i, 2 * i + 16, 3, 3, i);
    /* Trees  0,  2,  4 have face 3 connected */
    /* Trees 16, 18, 20 face 3 */
  }
  /* Set the coordinates. Each tet is just a translated version of
   * the root tet */
  for (i = 0; i < num_trees; i++) {
    translate[0] = (i & 1) + 2 * !!(i & 8);
    translate[1] = !!(i & 2) + 2 * !!(i & 16);
    translate[2] = !!(i & 4) + 2 * !!(i & 32);
    t8_debugf ("%i  %.0f %.0f %.0f\n", i, translate[0], translate[1],
               translate[2]);
    t8_cmesh_translate_coordinates (vertices_coords, translated_coords, 4,
                                    translate);
    t8_cmesh_set_tree_vertices (cmesh, i, t8_get_package_id (), 0,
                                translated_coords, 4);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hybrid_gate (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[32];
  int                 i;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_HEX);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 2, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 0);

  /* Tetrahedron 1 vertices */
  vertices[0] = 0.43;
  vertices[1] = 0;
  vertices[2] = 2;

  vertices[3] = 0;
  vertices[4] = 0;
  vertices[5] = 1;

  vertices[6] = 0.86;
  vertices[7] = -0.5;
  vertices[8] = 1;

  vertices[9] = 0.86;
  vertices[10] = 0.5;
  vertices[11] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);

  /* Tetrahedron 2 vertices */
  for (i = 0; i < 3; i++) {
    vertices[i] = vertices[i] + (i == 0 ? 1 + 0.86 : 0);
    vertices[3 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
    vertices[9 + i] = vertices[9 + i] + (i == 0 ? 1 : 0);
  }
  vertices[6] = 1 + 2 * 0.86;
  vertices[7] = 0;
  vertices[8] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices, 4);

  /* Prism 1 vertices */

  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;

  vertices[3] = 0.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices, 6);

  /* Prism 2 vertices */

  for (i = 0; i < 3; i++) {
    vertices[3 + i] = vertices[i] + (i == 0 ? 1 + 2 * 0.86 : 0);
    vertices[6 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
  }

  vertices[0] = 0.86 + 1;
  vertices[1] = -0.5;
  vertices[2] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0, vertices, 6);

  /* Hex coordinates */
  vertices[0] = 0.86;
  vertices[1] = -0.5;
  vertices[2] = 0;

  vertices[3] = 1.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  vertices[9] = 1.86;
  vertices[10] = 0.5;
  vertices[11] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 4; i++) {
    vertices[12 + 3 * i] = vertices[3 * i];
    vertices[12 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[12 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }

  t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0, vertices, 8);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hybrid_gate_deformed (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  double              vertices[32];
  int                 i;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TET);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_HEX);
  t8_cmesh_set_join (cmesh, 0, 2, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 4, 0);
  t8_cmesh_set_join (cmesh, 2, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 1, 1, 0);

  /* Tetrahedron 1 vertices */
  vertices[0] = 1;
  vertices[1] = -1;
  vertices[2] = 2.7;

  vertices[3] = 0;
  vertices[4] = -0.5;
  vertices[5] = 2;

  vertices[6] = 0.86;
  vertices[7] = -0.5;
  vertices[8] = 1;

  vertices[9] = 0.86;
  vertices[10] = 0.5;
  vertices[11] = 1;

  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);

  /* Tetrahedron 2 vertices */
  for (i = 0; i < 3; i++) {
    vertices[i] = vertices[i] + (i == 0 ? 1 + 0.86 : 0);
    vertices[3 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
    vertices[9 + i] = vertices[9 + i] + (i == 0 ? 1 : 0);
  }
  vertices[0] = 1.7;
  vertices[1] = 0.3;
  vertices[2] = 2.5;

  vertices[6] = 1 + 2 * 0.86;
  vertices[7] = 0;
  vertices[8] = 1.2;

  vertices[3] = 1.5;
  vertices[4] = -0.2;
  vertices[5] = 0.8;

  t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0, vertices, 4);

  /* Prism 1 vertices */

  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;

  vertices[3] = 0.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[2] = 0.2;
  vertices[9] = 0;
  vertices[10] = -0.5;
  vertices[11] = 2;
  vertices[3] = 0.9;
  vertices[4] = -0.7;
  vertices[5] = 0.3;

  t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0, vertices, 6);

  /* Prism 2 vertices */

  for (i = 0; i < 3; i++) {
    vertices[3 + i] = vertices[i] + (i == 0 ? 1 + 2 * 0.86 : 0);
    vertices[6 + i] = vertices[6 + i] + (i == 0 ? 1 : 0);
  }

  vertices[0] = 0.86 + 1;
  vertices[1] = -0.5;
  vertices[2] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 3; i++) {
    vertices[9 + 3 * i] = vertices[3 * i];
    vertices[9 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[9 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[6] = 2;
  vertices[7] = 0.2;
  vertices[8] = -0.3;

  vertices[9] = 1.5;
  vertices[10] = -0.2;
  vertices[11] = 0.8;
  t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0, vertices, 6);

  /* Hex coordinates */
  vertices[0] = 0.9;
  vertices[1] = -0.7;
  vertices[2] = 0.3;

  vertices[3] = 1.86;
  vertices[4] = -0.5;
  vertices[5] = 0;

  vertices[6] = 0.86;
  vertices[7] = 0.5;
  vertices[8] = 0;

  vertices[9] = 1.86;
  vertices[10] = 0.5;
  vertices[11] = 0;

  /* Translate +1 in z-axis for the upper vertices */
  for (i = 0; i < 4; i++) {
    vertices[12 + 3 * i] = vertices[3 * i];
    vertices[12 + 3 * i + 1] = vertices[3 * i + 1];
    vertices[12 + 3 * i + 2] = vertices[3 * i + 2] + 1;
  }
  vertices[9] = 2;
  vertices[10] = 0.2;
  vertices[11] = -0.3;

  vertices[12] = 0.86;
  vertices[13] = -0.5;
  vertices[14] = 1;

  vertices[15] = 1.5;
  vertices[16] = -0.2;
  vertices[17] = 0.8;

  t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0, vertices, 8);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
