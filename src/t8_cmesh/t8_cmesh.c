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
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_refcount.h>
#include <t8_data/t8_shmem.h>
#include <t8_vec.h>
#include <t8_eclass.h>
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
  if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) && !cmesh->committed)) {
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
  int ieclass;
  t8_gloidx_t glo_trees = 0;
  t8_locidx_t lo_trees = 0;
  int ret = 0;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  for (ieclass = 0; ieclass < T8_ECLASS_COUNT; ieclass++) {
    ret = ret && cmesh->num_local_trees_per_eclass[ieclass] <= cmesh->num_trees_per_eclass[ieclass];
    lo_trees += cmesh->num_local_trees_per_eclass[ieclass];
    glo_trees += cmesh->num_trees_per_eclass[ieclass];
  }
  return !ret && lo_trees == cmesh->num_local_trees && glo_trees == cmesh->num_trees;
}
#endif

int
t8_cmesh_is_committed (t8_cmesh_t cmesh)
{
  static int is_checking = 0;

  /* We run into a stackoverflow if routines that we call here,
   * also call t8_cmesh_is_committed.
   * We prevent this with the static variable is_checking.
   * This variable lives beyond one execution of t8_cmesh_is_committed.
   * We use it as a form of lock to prevent entering an infinite recursion.
   */
  /* TODO: This is_checking is not thread safe. If two threads call cmesh routines
   *       that call t8_cmesh_is_committed, only one of them will correctly check the cmesh. */
  if (!is_checking) {
    is_checking = 1;

    if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) && cmesh->committed)) {
      is_checking = 0;
      return 0;
    }

#ifdef T8_ENABLE_DEBUG
    /* TODO: check more conditions that must always hold after commit */
    if ((!t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) || (!t8_cmesh_no_negative_volume (cmesh))
        || (!t8_cmesh_check_trees_per_eclass (cmesh))) {
      is_checking = 0;
      return 0;
    }
    if (t8_cmesh_get_num_local_trees (cmesh) > 0 && t8_cmesh_is_empty (cmesh)) {
      is_checking = 0;
      return 0;
    }
#endif
    is_checking = 0;
  }
  return 1;
}

/* Check whether a given communicator assigns the same rank and mpisize
 * as stored in a given cmesh. */
int
t8_cmesh_comm_is_valid (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int mpiret, mpisize, mpirank;

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
t8_cmesh_init (t8_cmesh_t *pcmesh)
{
  t8_cmesh_t cmesh;
  T8_ASSERT (pcmesh != NULL);

  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  t8_refcount_init (&cmesh->rc);

  /* sensible (hard error) defaults */
  cmesh->set_refine_level = 0; /*< sensible default TODO document */
  cmesh->set_partition_level = -1;
  cmesh->dimension = -1; /*< ok; force user to select dimension */
  cmesh->mpirank = -1;
  cmesh->mpisize = -1;
  cmesh->first_tree = -1;
  cmesh->first_tree_shared = -1;
  cmesh->face_knowledge = 3; /*< sensible default TODO document */
  t8_stash_init (&cmesh->stash);
  /* Set the geometry handler to NULL.
   * It will get initialized either when a geometry is registered
   * or when the cmesh gets committed. */
  cmesh->geometry_handler = NULL;

  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
}

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

t8_shmem_array_t
t8_cmesh_alloc_offsets (int mpisize, sc_MPI_Comm comm)
{
  t8_shmem_array_t offsets;
#ifdef T8_ENABLE_DEBUG
  int mpisize_debug, mpiret;
  mpiret = sc_MPI_Comm_size (comm, &mpisize_debug);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (mpisize == mpisize_debug);
#endif

  t8_shmem_array_init (&offsets, sizeof (t8_gloidx_t), mpisize + 1, comm);
  t8_debugf ("Allocating shared array with type %s\n", sc_shmem_type_to_string[sc_shmem_get_type (comm)]);
  return offsets;
}

void
t8_cmesh_set_partition_range (t8_cmesh_t cmesh, int set_face_knowledge, t8_gloidx_t first_local_tree,
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
t8_cmesh_set_partition_offsets (t8_cmesh_t cmesh, t8_shmem_array_t tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  if (cmesh->tree_offsets != NULL && cmesh->tree_offsets != tree_offsets) {
    /* We overwrite a previously set offset array, so
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
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, int element_level, t8_scheme_cxx_t *ts)
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

void
t8_cmesh_set_refine (t8_cmesh_t cmesh, int level, t8_scheme_cxx_t *scheme)
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
t8_cmesh_treeid_is_local_tree (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return 0 <= ltreeid && ltreeid < t8_cmesh_get_num_local_trees (cmesh);
}

int
t8_cmesh_treeid_is_ghost (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghosts = t8_cmesh_get_num_ghosts (cmesh);

  return num_trees <= ltreeid && ltreeid < num_trees + num_ghosts;
}

t8_locidx_t
t8_cmesh_ltreeid_to_ghostid (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
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
  return tree->treeid < cmesh->num_local_trees - 1 ? t8_cmesh_get_tree (cmesh, tree->treeid + 1) : NULL;
}

void
t8_cmesh_set_attribute (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key, void *data, size_t data_size,
                        int data_persists)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, gtree_id, package_id, key, data_size, data, !data_persists);
}

void
t8_cmesh_set_attribute_string (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key, const char *string)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  /* The size is the string's length + the terminating '\0' */
  size_t size = strlen (string) + 1;
  /* Add the string as an attribute. */
  t8_cmesh_set_attribute (cmesh, gtree_id, package_id, key, (void *) string, size, 0);
}

double *
t8_cmesh_get_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid));

  return (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, ltreeid);
}

void *
t8_cmesh_get_attribute (t8_cmesh_t cmesh, int package_id, int key, t8_locidx_t ltree_id)
{
  int is_ghost;

  T8_ASSERT (cmesh->committed);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id) || t8_cmesh_treeid_is_ghost (cmesh, ltree_id));
  is_ghost = t8_cmesh_treeid_is_ghost (cmesh, ltree_id);

  if (is_ghost) {
    ltree_id = t8_cmesh_ltreeid_to_ghostid (cmesh, ltree_id);
  }
  return t8_cmesh_trees_get_attribute (cmesh->trees, ltree_id, package_id, key, NULL, is_ghost);
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

void
t8_cmesh_set_dimension (t8_cmesh_t cmesh, int dim)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (0 <= dim && dim <= T8_ECLASS_MAX_DIM);

  cmesh->dimension = dim;
}

void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, t8_eclass_t tree_class)
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

/* Given a set of vertex coordinates for a tree of a given eclass.
 * Query whether the geometric volume of the tree with this coordinates
 * would be negative.
 * Returns true if a tree of the given eclass with the given vertex
 * coordinates does have negative volume.
 */
int
t8_cmesh_tree_vertices_negative_volume (t8_eclass_t eclass, double *vertices, int num_vertices)
{
  double v_1[3], v_2[3], v_j[3], cross[3], sc_prod;
  int i, j;

  T8_ASSERT (num_vertices == t8_eclass_num_vertices[eclass]);

  if (t8_eclass_to_dimension[eclass] <= 2) {
    /* Only three dimensional eclass do have a volume */
    return 0;
  }

  T8_ASSERT (eclass == T8_ECLASS_TET || eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_PRISM
             || eclass == T8_ECLASS_PYRAMID);
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
  t8_vec_cross (v_1, v_2, cross);
  /* Compute sc_prod = <v_j, cross> */
  sc_prod = t8_vec_dot (v_j, cross);

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
  t8_locidx_t itree;
  double *vertices;
  t8_eclass_t eclass;
  int ret, res = 0;

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
      ret = t8_cmesh_tree_vertices_negative_volume (eclass, vertices, t8_eclass_num_vertices[eclass]);
      if (ret) {
        t8_debugf ("Detected negative volume in tree %li\n", (long) itree);
      }
      res |= ret; /* res is true if one ret value is true */
    }
  }
  return !res;
}
#endif

void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltree_id, double *vertices, int num_vertices)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (vertices != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, ltree_id, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY,
                          3 * num_vertices * sizeof (double), (void *) vertices, 1);
}

void
t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t gtree1, t8_gloidx_t gtree2, int face1, int face2, int orientation)
{
  T8_ASSERT (0 <= orientation);

  t8_stash_add_facejoin (cmesh->stash, gtree1, gtree2, face1, face2, orientation);
}

/* Allocate a cmesh profile if not yet present and set default
 * values. */
static void
t8_cmesh_init_profile (t8_cmesh_t cmesh)
{
  if (cmesh->profile == NULL) {
    /* Allocate new profile if it is not enabled already */
    cmesh->profile = T8_ALLOC_ZERO (t8_cprofile_struct_t, 1);
  }
  /* Set default values */
  cmesh->profile->commit_runtime = 0;
  cmesh->profile->first_tree_shared = -1; /* invalid until commit */
  cmesh->profile->geometry_evaluate_runtime = 0;
  cmesh->profile->geometry_evaluate_num_calls = 0;
  cmesh->profile->partition_bytes_sent = 0;
  cmesh->profile->partition_ghosts_recv = 0;
  cmesh->profile->partition_ghosts_shipped = 0;
  cmesh->profile->partition_procs_sent = 0;
  cmesh->profile->partition_runtime = 0;
  cmesh->profile->partition_trees_recv = 0;
  cmesh->profile->partition_trees_shipped = 0;
}

void
t8_cmesh_set_profiling (t8_cmesh_t cmesh, int set_profiling)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  if (set_profiling) {
    t8_cmesh_init_profile (cmesh);
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
  int is_equal;
  T8_ASSERT (cmesh_a != NULL && cmesh_b != NULL);

  if (cmesh_a == cmesh_b) {
    return 1;
  }
  /* check entries that are numbers */
  is_equal = cmesh_a->committed != cmesh_b->committed || cmesh_a->dimension != cmesh_b->dimension
             || cmesh_a->set_partition != cmesh_b->set_partition || cmesh_a->mpirank != cmesh_b->mpirank
             || cmesh_a->mpisize != cmesh_b->mpisize || cmesh_a->num_trees != cmesh_b->num_trees
             || cmesh_a->num_local_trees != cmesh_b->num_local_trees || cmesh_a->num_ghosts != cmesh_b->num_ghosts
             || cmesh_a->first_tree != cmesh_b->first_tree;
  if (is_equal != 0) {
    return 0;
  }
  /* check arrays */
  is_equal
    = memcmp (cmesh_a->num_trees_per_eclass, cmesh_b->num_trees_per_eclass, T8_ECLASS_COUNT * sizeof (t8_gloidx_t));
  is_equal = is_equal
             || memcmp (cmesh_a->num_local_trees_per_eclass, cmesh_b->num_local_trees_per_eclass,
                        T8_ECLASS_COUNT * sizeof (t8_locidx_t));

  /* check tree_offsets */
  if (cmesh_a->tree_offsets != NULL) {
    if (cmesh_b->tree_offsets == NULL) {
      return 0;
    }
    else {
      is_equal = is_equal || !t8_shmem_array_is_equal (cmesh_a->tree_offsets, cmesh_b->tree_offsets);
    }
  }
  if (is_equal != 0) {
    return 0;
  }
  /* check trees */
  if (cmesh_a->committed && !t8_cmesh_trees_is_equal (cmesh_a, cmesh_a->trees, cmesh_b->trees)) {
    /* if we have committed check tree arrays */
    return 0;
  }
  else {
    if (!cmesh_a->committed && !t8_stash_is_equal (cmesh_a->stash, cmesh_b->stash)) {
      /* if we have not committed check stash arrays */
      return 0;
    }
  }
  return 1;
}

int
t8_cmesh_is_empty (t8_cmesh_t cmesh)
{
  return cmesh->num_trees == 0;
}

t8_cmesh_t
t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int mpirank, mpisize, mpiret;
  int iclass;
  t8_cmesh_t cmesh_out = NULL; /* NULL initializer prevents compiler warning. */

  struct
  {
    t8_cmesh_struct_t cmesh;
    t8_gloidx_t num_trees_per_eclass[T8_ECLASS_COUNT];
    size_t stash_elem_counts[3];
    int pre_commit; /* True, if cmesh on root is not committed yet. */
#ifdef T8_ENABLE_DEBUG
    sc_MPI_Comm comm;
#endif
  } meta_info;

  /* TODO: BUG: running with two processes and a cmesh of one T8_ECLASS_LINE,
   *       then on both processes the face_neighbors and vertices arrays of
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
    /* Check whether geometries are set. If so, abort.
     * We cannot broadcast the geometries, since they are pointers to derived 
     * classes that we cannot know of on the receiving process.
     * Geometries must therefore be added after broadcasting. */
    SC_CHECK_ABORT (cmesh_in->geometry_handler == NULL,
                    "Error: Broadcasting a cmesh with registered geometries is not possible.\n"
                    "We recommend to broadcast first and register the geometries after.\n");
    memcpy (&meta_info.cmesh, cmesh_in, sizeof (*cmesh_in));
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      meta_info.num_trees_per_eclass[iclass] = cmesh_in->num_trees_per_eclass[iclass];
      T8_ASSERT (cmesh_in->num_local_trees_per_eclass[iclass] == cmesh_in->num_trees_per_eclass[iclass]);
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

    /* Root returns the input cmesh */
    cmesh_out = cmesh_in;
  }
  /* TODO: we could optimize this by using IBcast */
  mpiret = sc_MPI_Bcast (&meta_info, sizeof (meta_info), sc_MPI_BYTE, root, comm);

  SC_CHECK_MPI (mpiret);
#ifdef T8_ENABLE_DEBUG
  mpiret = sc_MPI_Comm_dup (comm, &(meta_info.comm));
  SC_CHECK_MPI (mpiret);
#endif

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
      cmesh_out->num_trees_per_eclass[iclass] = meta_info.num_trees_per_eclass[iclass];
      cmesh_out->num_local_trees_per_eclass[iclass] = meta_info.num_trees_per_eclass[iclass];
    }
#ifdef T8_ENABLE_DEBUG
    int result;
    mpiret = sc_MPI_Comm_compare (comm, meta_info.comm, &result);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (result == sc_MPI_CONGRUENT);
#endif
  }
  if (meta_info.pre_commit) {
    /* broadcast all the stashed information about trees/neighbors/attributes */
    t8_stash_bcast (cmesh_out->stash, root, comm, meta_info.stash_elem_counts);
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
  mpiret = sc_MPI_Comm_free (&meta_info.comm);
  SC_CHECK_MPI (mpiret);
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
  int mpisize, mpiret;
  idx_t idx_mpisize;
  idx_t ncon = 1, elements;
  idx_t volume, *partition, ipart, newpart;
  int num_faces, iface, count_face;
  idx_t *xadj, *adjncy;
  int success;
  t8_locidx_t *new_number, itree, *tree_per_part_off, *tree_per_part;
  t8_locidx_t *face_neighbor;
  t8_locidx_t neigh_id;
  t8_ctree_t tree;

  /* cmesh must be committed and not partitioned */
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (!cmesh->set_partition);

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  idx_mpisize = mpisize;
  SC_CHECK_MPI (mpiret);

  elements = cmesh->num_trees;
  T8_ASSERT ((t8_locidx_t) elements == cmesh->num_trees);

  /* Count the number of tree-to-tree connections via a face */
  num_faces = 0;
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor, NULL);
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      if (face_neighbor[iface] >= 0)
        num_faces++;
    }
  }

  /* xadj and adjncy store the face-connections in a CSR format
   * xadj[treeid] = offset of the tree in adjncy
   * adjncy[xadj[treeid]]...adjncy[xadj[treeid]-1] are the trees with which
   * the tree has a face connection */
  xadj = T8_ALLOC_ZERO (idx_t, elements + 1);
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
  partition = T8_ALLOC (idx_t, elements);
  /* partition the elements in mpisize many partitions */
  success = METIS_PartGraphRecursive (&elements, &ncon, xadj, adjncy, NULL, NULL, NULL, &idx_mpisize, NULL, NULL, NULL,
                                      &volume, partition);
  T8_ASSERT (success == METIS_OK);
  /* memory to store the new treeid of a tree */
  new_number = T8_ALLOC (t8_locidx_t, cmesh->num_trees);
  /* Store the number of trees pointer partition */
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
    new_number[itree] = tree_per_part_off[newpart + 1] - tree_per_part[newpart];
    tree_per_part[newpart]--;
  }
  /* Set for each tree its new treeid and the new ids of its neighbors */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor, NULL);
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
t8_cmesh_tree_face_is_boundary (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int face)
{
  int8_t *ttf;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  if (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid)) {
    /* The local tree id belongs to a tree */
    t8_locidx_t *face_neighbor;
    (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, ltreeid, &face_neighbor, &ttf);

    if (face_neighbor[face] == ltreeid && ttf[face] == face) {
      /* The tree is connected to itself at the same face.
       * Thus this is a domain boundary */
      return 1;
    }
  }
  else {
    /* The local tree id belongs to a ghost */
    T8_ASSERT (t8_cmesh_treeid_is_ghost (cmesh, ltreeid));

    t8_gloidx_t *face_neighbor;
    const t8_locidx_t lghostid = t8_cmesh_ltreeid_to_ghostid (cmesh, ltreeid);
    (void) t8_cmesh_trees_get_ghost_ext (cmesh->trees, lghostid, &face_neighbor, &ttf);

    if (face_neighbor[face] == t8_cmesh_get_global_id (cmesh, ltreeid) && ttf[face] == face) {
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
  t8_ctree_t tree;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id));

  tree = t8_cmesh_get_tree (cmesh, ltree_id);
  return tree->eclass;
}

t8_eclass_t
t8_cmesh_get_ghost_class (t8_cmesh_t cmesh, t8_locidx_t lghost_id)
{
  t8_cghost_t ghost;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (0 <= lghost_id && lghost_id < cmesh->num_ghosts);

  ghost = t8_cmesh_trees_get_ghost (cmesh->trees, lghost_id);
  return ghost->eclass;
}

t8_gloidx_t
t8_cmesh_get_global_id (t8_cmesh_t cmesh, t8_locidx_t local_id)
{
  T8_ASSERT (0 <= local_id && local_id < cmesh->num_ghosts + cmesh->num_local_trees);
  if (local_id < cmesh->num_local_trees) {
    return local_id + cmesh->first_tree;
  }
  else {
    return t8_cmesh_trees_get_ghost (cmesh->trees, local_id - cmesh->num_local_trees)->treeid;
  }
}

t8_locidx_t
t8_cmesh_get_local_id (t8_cmesh_t cmesh, t8_gloidx_t global_id)
{
  t8_gloidx_t temp_local_id;
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
t8_cmesh_get_face_neighbor (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int face, int *dual_face,
                            int *orientation)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid) || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  const int is_ghost = t8_cmesh_treeid_is_ghost (cmesh, ltreeid);
  int8_t ttf;
  t8_locidx_t face_neigh;
  int dual_face_temp, orientation_temp;

  /* If this is a domain boundary, return -1 */
  if (t8_cmesh_tree_face_is_boundary (cmesh, ltreeid, face)) {
    return -1;
  }

  if (!is_ghost) {
    /* The local tree id belongs to a local tree (not a ghost) */
    /* Get the tree */
    const t8_ctree_t tree = t8_cmesh_get_tree (cmesh, ltreeid);

#ifdef T8_ENABLE_DEBUG
    /* Get the eclass */
    t8_eclass_t eclass = tree->eclass;
    /* Check that face is valid */
    T8_ASSERT (0 <= face && face < t8_eclass_num_faces[eclass]);
#endif

    /* Get the local id of the face neighbor */
    face_neigh = t8_cmesh_trees_get_face_neighbor_ext (tree, face, &ttf);
  }
  else {
    /* The local tree id belongs to a ghost */
    const t8_locidx_t lghostid = ltreeid - t8_cmesh_get_num_local_trees (cmesh);
    /* Get the ghost */
    const t8_cghost_t ghost = t8_cmesh_trees_get_ghost (cmesh->trees, lghostid);

    t8_gloidx_t global_face_neigh;

#ifdef T8_ENABLE_DEBUG
    /* Get the eclass */
    t8_eclass_t eclass = ghost->eclass;
    /* Check that face is valid */
    T8_ASSERT (0 <= face && face < t8_eclass_num_faces[eclass]);
#endif

    /* Get the global id of the face neighbor */
    global_face_neigh = t8_cmesh_trees_get_ghost_face_neighbor_ext (ghost, face, &ttf);
    /* Convert it into a local id */
    face_neigh = t8_cmesh_get_local_id (cmesh, global_face_neigh);

    /* TODO: Check whether this face is a boundary face */
    if (face_neigh < 0) {
      /* The neighbor is not local, return -1 */
      return -1;
    }
  }

  /* Decode the ttf information to get the orientation and the dual face */
  t8_cmesh_tree_to_face_decode (cmesh->dimension, ttf, &dual_face_temp, &orientation_temp);
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
    sc_statinfo_t stats[T8_CPROFILE_NUM_STATS];
    t8_cprofile_t *profile = cmesh->profile;

    /* Set the stats */
    sc_stats_set1 (&stats[0], profile->partition_trees_shipped, "cmesh: Number of trees sent.");
    sc_stats_set1 (&stats[1], profile->partition_ghosts_shipped, "cmesh: Number of ghosts sent.");
    sc_stats_set1 (&stats[2], profile->partition_trees_recv, "cmesh: Number of trees received.");
    sc_stats_set1 (&stats[3], profile->partition_ghosts_recv, "cmesh: Number of ghosts received.");
    sc_stats_set1 (&stats[4], profile->partition_bytes_sent, "cmesh: Number of bytes sent.");
    sc_stats_set1 (&stats[5], profile->partition_procs_sent, "cmesh: Number of processes sent to.");
    sc_stats_set1 (&stats[6], profile->first_tree_shared, "cmesh: First tree is shared.");
    sc_stats_set1 (&stats[7], profile->partition_runtime, "cmesh: Partition runtime.");
    sc_stats_set1 (&stats[8], profile->commit_runtime, "cmesh: Commit runtime.");
    sc_stats_set1 (&stats[9], profile->geometry_evaluate_num_calls, "cmesh: Number of geometry evaluations.");
    sc_stats_set1 (&stats[10], profile->geometry_evaluate_runtime, "cmesh: Accumulated geometry evaluation runtime.");
    /* compute stats */
    sc_stats_compute (sc_MPI_COMM_WORLD, T8_CPROFILE_NUM_STATS, stats);
    /* print stats */
    t8_logf (SC_LC_GLOBAL, SC_LP_STATISTICS, "Printing stats for cmesh.\n");
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, T8_CPROFILE_NUM_STATS, stats, 1, 1);
  }
}

static void
t8_cmesh_reset (t8_cmesh_t *pcmesh)
{
  t8_cmesh_t cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  /* free tree_offset */
  if (cmesh->tree_offsets != NULL) {
#if T8_ENABLE_DEBUG
    sc_MPI_Comm comm;
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

  /* Unref the geometry handler. */
  if (cmesh->geometry_handler != NULL) {
    t8_geom_handler_unref (&cmesh->geometry_handler);
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
t8_cmesh_unref (t8_cmesh_t *pcmesh)
{
  t8_cmesh_t cmesh;
  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  if (t8_refcount_unref (&cmesh->rc)) {
    t8_cmesh_reset (pcmesh);
  }
}

void
t8_cmesh_destroy (t8_cmesh_t *pcmesh)
{
  T8_ASSERT (pcmesh != NULL && *pcmesh != NULL && t8_refcount_is_last (&(*pcmesh)->rc));
  t8_cmesh_unref (pcmesh);
  T8_ASSERT (*pcmesh == NULL);
}

void
t8_cmesh_translate_coordinates (const double *coords_in, double *coords_out, int num_vertices, double translate[3])
{
  int i;

  for (i = 0; i < num_vertices; i++) {
    coords_out[3 * i] = coords_in[3 * i] + translate[0];
    coords_out[3 * i + 1] = coords_in[3 * i + 1] + translate[1];
    coords_out[3 * i + 2] = coords_in[3 * i + 2] + translate[2];
  }
}

/* TODO: This is just a helper function that was needed when we changed the vertex interface
 *       to use attributes. Before we stored a list of vertex coordinates in the cmesh and each tree indexed into this list.
 *       Now each tree carries the coordinates of its vertices.
 *       This function translates from the first approached to the second
 *       and was introduced to avoid rewriting the already existing cmesh_new... functions below.
 *       It would be nice to eventually rewrite these functions correctly.
 */
void
t8_cmesh_new_translate_vertices_to_attributes (const t8_locidx_t *tvertices, const double *vertices,
                                               double *attr_vertices, const int num_vertices)
{
  int i;
  for (i = 0; i < num_vertices; i++) {
    attr_vertices[3 * i] = vertices[3 * tvertices[i]];
    attr_vertices[3 * i + 1] = vertices[3 * tvertices[i] + 1];
    attr_vertices[3 * i + 2] = vertices[3 * tvertices[i] + 2];
  }
}

/* Compute y = ax + b on an array of doubles, interpreting
 * each 3 as one vector x */
void
t8_cmesh_coords_axb (const double *coords_in, double *coords_out, int num_vertices, double alpha, const double b[3])
{
  int i;

  for (i = 0; i < num_vertices; i++) {
    t8_vec_axpyz (coords_in + i * 3, b, coords_out + i * 3, alpha);
  }
}

#ifdef T8_ENABLE_DEBUG
/**
 * \warning This function is only available in debug-modus and should only 
 * be used in debug-modus.
 * 
 * Prints the vertices of each local tree. 
 * 
 * \param[in] cmesh   source-cmesh, which trees get printed.
 */
static void
t8_cmesh_print_local_trees (const t8_cmesh_t cmesh)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_vertices = t8_eclass_num_vertices[tree_class];
    const t8_gloidx_t gtree = t8_cmesh_get_global_id (cmesh, itree);
    for (int ivertex = 0; ivertex < num_vertices; ivertex++) {
      const int vert_x = 3 * ivertex;
      const int vert_y = 3 * ivertex + 1;
      const int vert_z = 3 * ivertex + 2;
      t8_debugf ("[Global_tree: %li, local_tree: %i, vertex: %i] eclass: %s\t %f, %f, %f\n", gtree, itree, ivertex,
                 t8_eclass_to_string[tree_class], vertices[vert_x], vertices[vert_y], vertices[vert_z]);
    }
    t8_debugf ("\n");
  }
}
#endif

void
t8_cmesh_debug_print_trees (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
#ifdef T8_ENABLE_DEBUG
  /* This function is probably rather slow, linear in the number of processes and therefore
   * only available if the debug-modus is enabled. */
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  if (t8_cmesh_is_partitioned (cmesh)) {
    /* The cmesh is partitioned */
    int rank;
    int size;
    int mpiret;
    mpiret = sc_MPI_Comm_rank (comm, &rank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &size);
    SC_CHECK_MPI (mpiret);
    const int send_rank = (rank + 1) % size;
    const int recv_rank = rank - 1;
    /* Print the local trees in process-order. */
    if (rank != 0) {
      int source_rank;
      /* Send the rank as an additional check */
      mpiret = sc_MPI_Recv (&source_rank, 1, sc_MPI_INT, recv_rank, 0, comm, sc_MPI_STATUS_IGNORE);
      SC_CHECK_MPI (mpiret);
      T8_ASSERT (source_rank == (rank - 1));
    }
    t8_cmesh_print_local_trees (cmesh);
    if (rank != (size - 1)) {
      mpiret = sc_MPI_Send (&rank, 1, sc_MPI_INT, send_rank, 0, comm);
      SC_CHECK_MPI (mpiret);
    }
  }
  else {
    /* The cmesh is not partitioned, only one rank prints the trees. */
    if (cmesh->mpirank == 0) {
      t8_cmesh_print_local_trees (cmesh);
    }
  }

#else
  t8_global_errorf ("Do not call t8_cmesh_debug_print_trees if t8code is not compiled with --enable-debug.\n");
#endif /* T8_ENABLE_DEBUG */
}
