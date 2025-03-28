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
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_refcount.h>
#include <t8_data/t8_shmem.h>
#include <t8_types/t8_vec.h>
#include <t8_eclass.h>
#include "t8_cmesh_types.h"
#ifdef T8_WITH_METIS
#include <metis.h>

#endif
#include "t8_cmesh_trees.h"

#include <t8_element.h>

/** \file t8_cmesh.cxx
 *  This file collects all general cmesh routines that need c++ compilation.
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
t8_cmesh_is_committed (const t8_cmesh_t cmesh)
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
    if ((!t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) || (!t8_cmesh_check_trees_per_eclass (cmesh))) {
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

#if T8_ENABLE_DEBUG
int
t8_cmesh_validate_geometry (const t8_cmesh_t cmesh)
{
  /* After a cmesh is committed, check whether all trees in a cmesh are compatible
 * with their geometry and if they have positive volume.
 * Returns true if all trees are valid. Returns also true if no geometries are
 * registered yet, since the validity computation depends on the used geometry.
 */

  /* Geometry handler is not constructed yet */
  if (cmesh->geometry_handler == NULL) {
    return true;
  }
  if (cmesh == NULL) {
    return true;
  }
  if (cmesh->geometry_handler->get_num_geometries () > 0) {
    /* Iterate over all trees, get their vertices and check the volume */
    for (t8_locidx_t itree = 0; itree < cmesh->num_local_trees; itree++) {
      /* Check if tree and geometry are compatible. */
      const int geometry_compatible
        = cmesh->geometry_handler->tree_compatible_with_geom (cmesh, t8_cmesh_get_global_id (cmesh, itree));
      if (!geometry_compatible) {
        t8_debugf ("Detected incompatible geometry for tree %li\n", (long) itree);
        return false;
      }
      if (geometry_compatible) {
        /* Check for negative volume. This only makes sense if the geometry is valid for the tree. */
        const int negative_volume
          = cmesh->geometry_handler->tree_negative_volume (cmesh, t8_cmesh_get_global_id (cmesh, itree));
        if (negative_volume) {
          t8_debugf ("Detected negative volume in tree %li\n", (long) itree);
          return false;
        }
      }
    }
  }
  return true;
}
#endif /* T8_ENABLE_DEBUG */

/* Check whether a given communicator assigns the same rank and mpisize
 * as stored in a given cmesh. */
int
t8_cmesh_comm_is_valid (const t8_cmesh_t cmesh, sc_MPI_Comm comm)
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

t8_cmesh_t
t8_cmesh_new ()
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  return cmesh;
}

void
t8_cmesh_set_derive (const t8_cmesh_t cmesh, const t8_cmesh_t set_from)
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
    SC_CHECK_ABORT (cmesh->stash->attributes.elem_count == 0,
                    "ERROR: Cannot add attributes to cmesh when deriving from another cmesh.\n");
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
t8_cmesh_set_partition_range (t8_cmesh_t cmesh, const int set_face_knowledge, const t8_gloidx_t first_local_tree,
                              const t8_gloidx_t last_local_tree)
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
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, const int element_level, const t8_scheme *scheme)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (element_level >= -1);
  T8_ASSERT (scheme != NULL);

  cmesh->set_partition = 1;
  cmesh->set_partition_level = element_level;
  cmesh->set_partition_scheme = scheme;
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

t8_gloidx_t
t8_cmesh_get_first_treeid (const t8_cmesh_t cmesh)
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
t8_cmesh_get_tree (const t8_cmesh_t cmesh, const t8_locidx_t ltree_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id));

  return t8_cmesh_trees_get_tree (cmesh->trees, ltree_id);
}

/* Returns the first local tree.
 * Returns NULL if there are no local trees. */
/* TODO: hide */
t8_ctree_t
t8_cmesh_get_first_tree (const t8_cmesh_t cmesh)
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
t8_cmesh_get_next_tree (const t8_cmesh_t cmesh, const t8_ctree_t tree)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (tree != NULL);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, tree->treeid));
  T8_ASSERT (cmesh->committed);
  return tree->treeid < cmesh->num_local_trees - 1 ? t8_cmesh_get_tree (cmesh, tree->treeid + 1) : NULL;
}

void
t8_cmesh_set_attribute (t8_cmesh_t cmesh, const t8_gloidx_t gtree_id, const int package_id, const int key,
                        void *const data, const size_t data_size, const int data_persists)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  SC_CHECK_ABORT (cmesh->set_from == NULL, "ERROR: Cannot add attributes to cmesh when deriving from another cmesh.\n");

  t8_stash_add_attribute (cmesh->stash, gtree_id, package_id, key, data_size, data, !data_persists);
}

void
t8_cmesh_set_attribute_string (t8_cmesh_t cmesh, const t8_gloidx_t gtree_id, const int package_id, const int key,
                               const char *string)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  /* The size is the string's length + the terminating '\0' */
  size_t size = strlen (string) + 1;
  /* Add the string as an attribute. */
  t8_cmesh_set_attribute (cmesh, gtree_id, package_id, key, (void *) string, size, 0);
}

void
t8_cmesh_set_attribute_gloidx_array (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key,
                                     const t8_gloidx_t *data, size_t data_count, int data_persists)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  const size_t data_size = data_count * sizeof (*data);
  t8_cmesh_set_attribute (cmesh, gtree_id, package_id, key, (void *) data, data_size, data_persists);
}

double *
t8_cmesh_get_tree_vertices (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid) || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  return (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, ltreeid);
}

void *
t8_cmesh_get_attribute (const t8_cmesh_t cmesh, const int package_id, const int key, const t8_locidx_t ltree_id)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id) || t8_cmesh_treeid_is_ghost (cmesh, ltree_id));
  const int is_ghost = t8_cmesh_treeid_is_ghost (cmesh, ltree_id);

  return t8_cmesh_trees_get_attribute (
    cmesh->trees, is_ghost ? t8_cmesh_ltreeid_to_ghostid (cmesh, ltree_id) : ltree_id, package_id, key, NULL, is_ghost);
}

t8_gloidx_t *
t8_cmesh_get_attribute_gloidx_array (const t8_cmesh_t cmesh, const int package_id, const int key,
                                     const t8_locidx_t ltree_id,
                                     [[maybe_unused]] const size_t data_count)  //TODO: remove data_count
{
  return (t8_gloidx_t *) t8_cmesh_get_attribute (cmesh, package_id, key, ltree_id);
}

t8_shmem_array_t
t8_cmesh_get_partition_table (const t8_cmesh_t cmesh)
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
t8_cmesh_set_dimension (t8_cmesh_t cmesh, const int dim)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (0 <= dim && dim <= T8_ECLASS_MAX_DIM);

  cmesh->dimension = dim;
}

int
t8_cmesh_get_dimension (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (0 <= cmesh->dimension && cmesh->dimension <= T8_ECLASS_MAX_DIM);

  return cmesh->dimension;
}

void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, const t8_gloidx_t gtree_id, const t8_eclass_t tree_class)
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
t8_cmesh_tree_vertices_negative_volume (const t8_eclass_t eclass, const double *vertices, const int num_vertices)
{
  T8_ASSERT (num_vertices == t8_eclass_num_vertices[eclass]);

  /* Points and lines do not have a volume orientation. */
  if (t8_eclass_to_dimension[eclass] < 2) {
    return 0;
  }

  T8_ASSERT (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TET
             || eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_PYRAMID);

  /* Skip negative volume check (orientation of face normal) of 2D elements
   * when z-coordinates are not (almost) zero. */
  if (t8_eclass_to_dimension[eclass] < 3) {
    for (int ivert = 0; ivert < num_vertices; ivert++) {
      const double z_coordinate = vertices[3 * ivert + 2];
      if (std::abs (z_coordinate) > 10 * T8_PRECISION_EPS) {
        return false;
      }
    }
  }

  /*
   *      z             For 2D meshes we enforce the right-hand-rule in terms
   *      |             of node ordering. The volume is defined by the parallelepiped
   *      | 2- - -(3)   spanned by the vectors between nodes 0:1 and 0:2 as well as the
   *      |/____ /      unit vector in z-direction. This definition works for both triangles and quads.
   *      0     1
   *
   *      6 ______  7   For Hexes and pyramids, if the vertex 4 is below the 0-1-2-3 plane,
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

  /* Build the vectors v_i as vertices_i - vertices_0. */
  double v_1[3], v_2[3], v_j[3], cross[3], sc_prod;

  if (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_QUAD) {
    for (int i = 0; i < 3; i++) {
      v_1[i] = vertices[3 + i] - vertices[i];
      v_2[i] = vertices[6 + i] - vertices[i];
    }

    /* Unit vector in z-direction. */
    v_j[0] = 0.0;
    v_j[1] = 0.0;
    v_j[2] = 1.0;

    /* Compute cross = v_1 x v_2. */
    t8_cross_3D (v_1, v_2, cross);
    /* Compute sc_prod = <v_j, cross>. */
    sc_prod = t8_dot (v_j, cross);

    T8_ASSERT (sc_prod != 0);
    return sc_prod < 0;
  }

  int j;
  if (eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM) {
    /* In the tet/prism case, the third vector is v_3 */
    j = 3;
  }
  else {
    /* For pyramids and Hexes, the third vector is v_4 */
    j = 4;
  }
  for (int i = 0; i < 3; i++) {
    v_1[i] = vertices[3 + i] - vertices[i];
    v_2[i] = vertices[6 + i] - vertices[i];
    v_j[i] = vertices[3 * j + i] - vertices[i];
  }
  /* compute cross = v_1 x v_2 */
  t8_cross_3D (v_1, v_2, cross);
  /* Compute sc_prod = <v_j, cross> */
  sc_prod = t8_dot (v_j, cross);

  T8_ASSERT (sc_prod != 0);
  return eclass == T8_ECLASS_TET ? sc_prod > 0 : sc_prod < 0;
}

void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, const t8_gloidx_t gtree_id, const double *vertices,
                            const int num_vertices)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (vertices != NULL);
  T8_ASSERT (!cmesh->committed);

  const size_t data_size = 3 * num_vertices * sizeof (double);

  t8_cmesh_set_attribute (cmesh, gtree_id, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, (void *) vertices,
                          data_size, 0);
}

void
t8_cmesh_set_join (t8_cmesh_t cmesh, const t8_gloidx_t gtree1, const t8_gloidx_t gtree2, const int face1,
                   const int face2, const int orientation)
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
t8_cmesh_set_profiling (t8_cmesh_t cmesh, const int set_profiling)
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
t8_cmesh_is_equal (const t8_cmesh_t cmesh_a, const t8_cmesh_t cmesh_b)
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
t8_cmesh_is_empty (const t8_cmesh_t cmesh)
{
  return cmesh->num_trees == 0;
}

t8_cmesh_t
t8_cmesh_bcast (const t8_cmesh_t cmesh_in, const int root, sc_MPI_Comm comm)
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
    if (cmesh_in->geometry_handler != NULL) {
      SC_CHECK_ABORT (cmesh_in->geometry_handler->get_num_geometries () == 0,
                      "Error: Broadcasting a cmesh with registered geometries is not possible.\n"
                      "We recommend to broadcast first and register the geometries after.\n");
    }
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
t8_cmesh_is_partitioned (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->set_partition != 0;
}

t8_gloidx_t
t8_cmesh_get_num_trees (const t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_trees;
}

t8_locidx_t
t8_cmesh_get_num_local_trees (const t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->num_local_trees;
}

t8_locidx_t
t8_cmesh_get_num_ghosts (const t8_cmesh_t cmesh)
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
t8_cmesh_get_tree_class (const t8_cmesh_t cmesh, const t8_locidx_t ltree_id)
{
  t8_ctree_t tree;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltree_id));

  tree = t8_cmesh_get_tree (cmesh, ltree_id);
  return tree->eclass;
}

t8_eclass_t
t8_cmesh_get_ghost_class (const t8_cmesh_t cmesh, const t8_locidx_t lghost_id)
{
  t8_cghost_t ghost;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (0 <= lghost_id && lghost_id < cmesh->num_ghosts);

  ghost = t8_cmesh_trees_get_ghost (cmesh->trees, lghost_id);
  return ghost->eclass;
}

t8_gloidx_t
t8_cmesh_get_global_id (const t8_cmesh_t cmesh, const t8_locidx_t local_id)
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
t8_cmesh_get_local_id (const t8_cmesh_t cmesh, const t8_gloidx_t global_id)
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
t8_cmesh_print_profile (const t8_cmesh_t cmesh)
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
    if (t8_cmesh_is_committed (cmesh)) {
      comm = t8_shmem_array_get_comm (cmesh->tree_offsets);
      T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
    }
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

  if (cmesh->geometry_handler != NULL) {
    cmesh->geometry_handler->unref ();
    cmesh->geometry_handler = NULL;
  }

  /* unref the partition scheme (if set) */
  if (cmesh->set_partition_scheme != NULL) {
    cmesh->set_partition_scheme->unref ();
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
t8_cmesh_translate_coordinates (const double *coords_in, double *coords_out, const int num_vertices,
                                const double translate[3])
{
  for (int ivertex = 0; ivertex < num_vertices; ivertex++) {
    coords_out[3 * ivertex] = coords_in[3 * ivertex] + translate[0];
    coords_out[3 * ivertex + 1] = coords_in[3 * ivertex + 1] + translate[1];
    coords_out[3 * ivertex + 2] = coords_in[3 * ivertex + 2] + translate[2];
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
    t8_axpyz (coords_in + i * 3, b, coords_out + i * 3, alpha);
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
t8_cmesh_debug_print_trees ([[maybe_unused]] const t8_cmesh_t cmesh, [[maybe_unused]] sc_MPI_Comm comm)
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

void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, const int level, const t8_scheme *scheme, t8_gloidx_t *first_local_tree,
                         t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                         int8_t *first_tree_shared)
{
  int is_empty;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (level >= 0);
  T8_ASSERT (scheme != NULL);

  *first_local_tree = 0;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = 0;
  }
  *last_local_tree = 0;
  if (child_in_tree_end != NULL) {
    *child_in_tree_end = 0;
  }

  t8_gloidx_t global_num_children;
  t8_gloidx_t first_global_child;
  t8_gloidx_t child_in_tree_begin_temp;
  t8_gloidx_t last_global_child;
  t8_gloidx_t children_per_tree = 0;
#ifdef T8_ENABLE_DEBUG
  t8_gloidx_t prev_last_tree = -1;
#endif
  int tree_class;

  /* Compute the number of children on level in each tree */
  global_num_children = 0;
  for (tree_class = T8_ECLASS_ZERO; tree_class < T8_ECLASS_COUNT; ++tree_class) {
    /* We iterate over each element class and get the number of children for this
     * tree class.
     */
    if (cmesh->num_trees_per_eclass[tree_class] > 0) {
      children_per_tree = scheme->count_leaves_from_root (static_cast<t8_eclass_t> (tree_class), level);
      T8_ASSERT (children_per_tree >= 0);
      global_num_children += cmesh->num_trees_per_eclass[tree_class] * children_per_tree;
    }
  }
  T8_ASSERT (children_per_tree != 0);

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
    first_global_child = ((long double) global_num_children * cmesh->mpirank) / (double) cmesh->mpisize;
  }
  if (cmesh->mpirank != cmesh->mpisize - 1) {
    last_global_child = ((long double) global_num_children * (cmesh->mpirank + 1)) / (double) cmesh->mpisize;
  }
  else {
    last_global_child = global_num_children;
  }

  T8_ASSERT (0 <= first_global_child && first_global_child <= global_num_children);
  T8_ASSERT (0 <= last_global_child && last_global_child <= global_num_children);

  *first_local_tree = first_global_child / children_per_tree;
  child_in_tree_begin_temp = first_global_child - *first_local_tree * children_per_tree;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = child_in_tree_begin_temp;
  }

  *last_local_tree = (last_global_child - 1) / children_per_tree;

  is_empty = *first_local_tree >= *last_local_tree && first_global_child >= last_global_child;
  if (first_tree_shared != NULL) {
#ifdef T8_ENABLE_DEBUG
    prev_last_tree = (first_global_child - 1) / children_per_tree;
    T8_ASSERT (cmesh->mpirank > 0 || prev_last_tree <= 0);
#endif
    if (!is_empty && cmesh->mpirank > 0 && child_in_tree_begin_temp > 0) {
      /* We exclude empty partitions here, by def their first_tree_shared flag is zero */
      /* We also exclude that the previous partition was empty at the beginning of the
       * partitions array */
      /* We also exclude the case that we have the first global element but
       * are not rank 0. */
      *first_tree_shared = 1;
    }
    else {
      *first_tree_shared = 0;
    }
  }
  if (child_in_tree_end != NULL) {
    if (*last_local_tree > 0) {
      *child_in_tree_end = last_global_child - *last_local_tree * children_per_tree;
    }
    else {
      *child_in_tree_end = last_global_child;
    }
  }
  if (is_empty) {
    /* This process is empty */
    /* We now set the first local tree to the first local tree on the
     * next nonempty rank, and the last local tree to first - 1 */
    *first_local_tree = last_global_child / children_per_tree;
    if (first_global_child % children_per_tree != 0) {
      /* The next nonempty process shares this tree. */
      (*first_local_tree)++;
    }

    *last_local_tree = *first_local_tree - 1;
  }
}

/* #################### new for hybrid stuff starts here.  #################### */

/* This struct is used to temporarily store some context values
 * when we determine the partition bounds of processes in
 * t8_cmesh_determine_partition. */
typedef struct
{
  int num_procs;                   /* Number of MPI processes. */
  int process_offset;              /* The first process that this process sends to. */
  t8_gloidx_t global_num_elements; /* The global number of elements. */
} t8_cmesh_partition_query_t;

/* Helper function to set the correct return value in the case that 
 * the uniform partition is empty.
 * We set all values to -1 and first_tree_shared to 0. */
static void
t8_cmesh_uniform_set_return_parameters_to_empty (t8_gloidx_t *first_local_tree, t8_gloidx_t *child_in_tree_begin,
                                                 t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                                                 int8_t *first_tree_shared)
{
  *first_local_tree = *last_local_tree = -1;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = -1;
  }
  if (child_in_tree_end != NULL) {
    *child_in_tree_end = -1;
  }
  if (first_tree_shared != NULL) {
    *first_tree_shared = 0;
  }
}

/* Helper function to compute (A*B)/C for large integers, when A*B might not fit in a
 * 64-bit int anymore. Uses floor(floor(A/C)*B + ((A%C)/C)*B) instead. */
static inline t8_gloidx_t
t8_A_times_B_over_C_gloidx (t8_gloidx_t A, t8_gloidx_t B, t8_gloidx_t C)
{
  const t8_gloidx_t a_over_c = A / C;
  const t8_gloidx_t a_o_c_times_b = a_over_c * B;

  /* We check whether computing A/C * B will cause an overflow.
   * This can be achieved by checking if dividing the result by B
   * yields A/C again. */
  T8_ASSERT (a_over_c == 0 || a_o_c_times_b / a_over_c == B);

  return (t8_gloidx_t) (a_o_c_times_b + (((long double) (A % C)) / C) * B);
}

/* Version of t8_A_times_B_over_C where A is an integer.
 * We reuse the gloidx version but check before whether an int fits into
 * a t8_gloidx_t.
 * This function can also be used if B or C are ints. */
static inline t8_gloidx_t
t8_A_times_B_over_C_intA (int A, t8_gloidx_t B, t8_gloidx_t C)
{
  T8_ASSERT (sizeof (int) <= sizeof (t8_gloidx_t));
  return t8_A_times_B_over_C_gloidx (A, B, C);
}

/* Compute and return the first global element index of a process in the 
 * uniform partition.
 * The first index of a process 0 <= p < P among E elements is
 * floor ((p * E) / P)
 */
static inline t8_gloidx_t
t8_cmesh_get_first_element_of_process (int process, int mpisize, t8_gloidx_t global_num_elements)
{
  return t8_A_times_B_over_C_intA (process, global_num_elements, mpisize);
}

void
t8_cmesh_uniform_bounds_equal_element_count (t8_cmesh_t cmesh, int level, t8_scheme *tree_scheme,
                                             t8_gloidx_t *first_local_tree, t8_gloidx_t *child_in_tree_begin,
                                             t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                                             int8_t *first_tree_shared)
{
  int is_empty;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (level >= 0);
  T8_ASSERT (tree_scheme != NULL);

  *first_local_tree = 0;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = 0;
  }
  *last_local_tree = 0;
  if (child_in_tree_end != NULL) {
    *child_in_tree_end = 0;
  }

  t8_gloidx_t global_num_children;
  t8_gloidx_t first_global_child;
  t8_gloidx_t child_in_tree_begin_temp;
  t8_gloidx_t last_global_child;
  t8_gloidx_t children_per_tree = 0;
#ifdef T8_ENABLE_DEBUG
  t8_gloidx_t prev_last_tree = -1;
#endif
  int tree_class;

  /* Compute the number of children on level in each tree */
  global_num_children = 0;
#ifdef T8_ENABLE_DEBUG
  t8_gloidx_t children_per_tree_check = 0;
  for (tree_class = T8_ECLASS_ZERO; tree_class < T8_ECLASS_COUNT; ++tree_class) {
    /* Get the number of children on level of the first tree class that is used*/
    if (cmesh->num_trees_per_eclass[tree_class] > 0) {
      children_per_tree_check = tree_scheme->count_leaves_from_root ((t8_eclass_t) tree_class, level);
      break;
    }
  }
  /* Compare it to the other used tree classes*/
  for (tree_class = T8_ECLASS_ZERO; tree_class < T8_ECLASS_COUNT; ++tree_class) {
    if (cmesh->num_trees_per_eclass[tree_class] > 0) {
      T8_ASSERT (children_per_tree_check == tree_scheme->count_leaves_from_root ((t8_eclass_t) tree_class, level));
    }
  }

#endif
  for (tree_class = T8_ECLASS_ZERO; tree_class < T8_ECLASS_COUNT; ++tree_class) {
    /* We iterate over each element class and get the number of children for this
     * tree class.
     */
    if (cmesh->num_trees_per_eclass[tree_class] > 0) {
      T8_ASSERT (tree_scheme != NULL);
      children_per_tree = tree_scheme->count_leaves_from_root ((t8_eclass_t) tree_class, level);
      T8_ASSERT (children_per_tree >= 0);
      global_num_children += cmesh->num_trees_per_eclass[tree_class] * children_per_tree;
    }
  }
  T8_ASSERT (children_per_tree != 0);

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
    first_global_child = ((long double) global_num_children * cmesh->mpirank) / (double) cmesh->mpisize;
  }
  if (cmesh->mpirank != cmesh->mpisize - 1) {
    last_global_child = ((long double) global_num_children * (cmesh->mpirank + 1)) / (double) cmesh->mpisize;
  }
  else {
    last_global_child = global_num_children;
  }

  T8_ASSERT (0 <= first_global_child && first_global_child <= global_num_children);
  T8_ASSERT (0 <= last_global_child && last_global_child <= global_num_children);

  *first_local_tree = first_global_child / children_per_tree;
  child_in_tree_begin_temp = first_global_child - *first_local_tree * children_per_tree;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = child_in_tree_begin_temp;
  }

  *last_local_tree = (last_global_child - 1) / children_per_tree;

  is_empty = *first_local_tree >= *last_local_tree && first_global_child >= last_global_child;
  if (first_tree_shared != NULL) {
#ifdef T8_ENABLE_DEBUG
    prev_last_tree = (first_global_child - 1) / children_per_tree;
    T8_ASSERT (cmesh->mpirank > 0 || prev_last_tree <= 0);
#endif
    if (!is_empty && cmesh->mpirank > 0 && child_in_tree_begin_temp > 0) {
      /* We exclude empty partitions here, by def their first_tree_shared flag is zero */
      /* We also exclude that the previous partition was empty at the beginning of the
       * partitions array */
      /* We also exclude the case that we have the first global element but
       * are not rank 0. */
      *first_tree_shared = 1;
    }
    else {
      *first_tree_shared = 0;
    }
  }
  if (child_in_tree_end != NULL) {
    if (*last_local_tree > 0) {
      *child_in_tree_end = last_global_child - *last_local_tree * children_per_tree;
    }
    else {
      *child_in_tree_end = last_global_child;
    }
  }
  if (is_empty) {
    /* This process is empty */
    /* We now set the first local tree to the first local tree on the
     * next nonempty rank, and the last local tree to first - 1 */
    *first_local_tree = last_global_child / children_per_tree;
    if (first_global_child % children_per_tree != 0) {
      /* The next nonempty process shares this tree. */
      (*first_local_tree)++;
    }

    *last_local_tree = *first_local_tree - 1;
  }
}

/* Given an array (partition) storing the global index of the 
 * first element of each (pure) local tree and
 * a (pure) local tree index, return the first process that
 * will have elements of this tree in a uniform partition. 
 * The data pointer must point to a valid t8_cmesh_partition_query_t,
 * storing the number of processes and the global number of elements. 
 * 
 * This function is used standalone and as callback of sc_array_split. */
static size_t
t8_cmesh_determine_partition (sc_array_t *first_element_tree, size_t pure_local_tree, void *data)
{
  T8_ASSERT (data != NULL);
  T8_ASSERT (first_element_tree != NULL);

  const t8_cmesh_partition_query_t *query_data = (const t8_cmesh_partition_query_t *) data;
  size_t first_proc_adjusted;
  const t8_gloidx_t element_index = *(t8_gloidx_t *) sc_array_index (first_element_tree, pure_local_tree);
  const t8_gloidx_t mirror_element_index = query_data->global_num_elements - element_index - 1;
  int first_proc_rank;

  if (element_index == query_data->global_num_elements) {
    return query_data->num_procs - query_data->process_offset;
  }
  else {
    first_proc_rank
      = query_data->num_procs - 1
        - t8_A_times_B_over_C_intA (query_data->num_procs, mirror_element_index, query_data->global_num_elements);

    /* If this is called by array split, we need to adjust for relative processes starting
     * add the first process. */
    first_proc_adjusted = first_proc_rank - query_data->process_offset;
  }

  /* Safety checks */
#ifdef T8_ENABLE_DEBUG
  T8_ASSERT (0 <= first_proc_rank && (int) first_proc_rank < query_data->num_procs);
  //T8_ASSERT (0 <= first_proc_adjusted);
  /* Check that the element lies in the partition of the computed proc. */
  T8_ASSERT (
    t8_cmesh_get_first_element_of_process (first_proc_rank, query_data->num_procs, query_data->global_num_elements)
    <= element_index);
  if ((int) first_proc_rank != query_data->num_procs - 1) {
    T8_ASSERT (t8_cmesh_get_first_element_of_process (first_proc_rank + 1, query_data->num_procs,
                                                      query_data->global_num_elements)
               > element_index);
  }
  if (element_index == query_data->global_num_elements) {
    T8_ASSERT ((int) first_proc_rank == query_data->num_procs);
  }
#endif
  return first_proc_adjusted;
}
void
t8_cmesh_uniform_bounds_from_unpartioned (t8_cmesh_t cmesh, const t8_gloidx_t local_num_children, const int level,
                                          const t8_scheme *scheme, t8_gloidx_t *first_local_tree,
                                          t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree,
                                          t8_gloidx_t *child_in_tree_end, int8_t *first_tree_shared)
{
  const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (cmesh);
  t8_debugf ("Cmesh is not partitioned.\n");
  /* Compute the first and last element of this process. Then loop over
     * all trees to find the trees in which these are contained.
     * We cast to long double and double to prevent overflow. */
  /* cmesh is replicated, therefore the computation of local_num_children equals the global number of children*/
  const t8_gloidx_t first_child
    = t8_cmesh_get_first_element_of_process (cmesh->mpirank, cmesh->mpisize, local_num_children);
  const t8_gloidx_t last_child
    = t8_cmesh_get_first_element_of_process (cmesh->mpirank + 1, cmesh->mpisize, local_num_children) - 1;
  /* Can't we optimize this linear loop by using a binary search?
     * -> No, we cannot. Since we need in any case compute the t8_element_count_leaves_from_root
     *    for each tree.
     */
  t8_gloidx_t current_tree_element_offset = 0;
  for (t8_gloidx_t igtree = 0; igtree < num_trees; ++igtree) {
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, (t8_locidx_t) igtree);
    /* TODO: We can optimize by buffering the elem_in_tree value. Thus, if 
         the computation is expensive (may be for non-morton-type schemes),
         we do it only once. */
    const t8_gloidx_t elem_in_tree = scheme->count_leaves_from_root (tree_class, level);
    /* Check if the first element is on the current tree */
    if (current_tree_element_offset <= first_child && first_child < current_tree_element_offset + elem_in_tree) {
      if (child_in_tree_begin != NULL) {
        *child_in_tree_begin = first_child - current_tree_element_offset;
      }
      *first_local_tree = igtree;
      /* If our first element is not the very first element in the tree, we share
         * this tree with the previous process. */
      if (first_tree_shared != NULL) {
        *first_tree_shared = current_tree_element_offset < first_child ? 1 : 0;
      }
    }
    if (last_child < first_child) {
      /* This process is empty. We needed to identify the 'first_local_tree' since it is the 
         * first local tree of the next process, which we store in this case.
         */
      t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                       child_in_tree_end, first_tree_shared);
      *first_local_tree = 0;
      return;
    }
    /* Check if the last element is on the current tree */
    if (current_tree_element_offset <= last_child && last_child < current_tree_element_offset + elem_in_tree) {
      if (child_in_tree_end != NULL) {
        *child_in_tree_end = last_child - current_tree_element_offset + 1;
      }
      *last_local_tree = igtree;

      /* We have found the last tree and can immediately return,
         * ending the for loop. */
      T8_ASSERT (*first_local_tree <= *last_local_tree);
      T8_ASSERT (0 <= *first_local_tree && *first_local_tree < num_trees);
      T8_ASSERT (0 <= *last_local_tree && *last_local_tree < num_trees);
      return;
    }
    current_tree_element_offset += elem_in_tree;
  }
  /* If we reach this part, we do not have any trees - the cmesh is empty */
  T8_ASSERT (num_trees == 0);
  t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                   child_in_tree_end, first_tree_shared);
  return;
}

void
t8_cmesh_bounds_send_start (t8_cmesh_t cmesh, const int proc_is_empty, const t8_locidx_t first_puretree_of_current_proc,
                            const int first_tree_shared_shift, const t8_gloidx_t iproc, sc_array *send_requests,
                            sc_array *send_buffer, int *current_pos_in_send_buffer,
                            const t8_gloidx_t first_element_in_tree_index_of_current_proc,
                            t8_gloidx_t *first_local_tree, int8_t *first_tree_shared, t8_gloidx_t *child_in_tree_begin,
                            bool *expect_start_message, t8_cmesh_partition_query_t *data, sc_MPI_Comm comm
#ifdef T8_ENABLE_DEBUG
                            ,
                            int *num_received_start_messages, int *num_message_sent
#endif
)
{
  const t8_gloidx_t global_id_of_first_tree
    = proc_is_empty ? -1 : first_puretree_of_current_proc + t8_cmesh_get_first_treeid (cmesh) + first_tree_shared_shift;

  if (iproc != cmesh->mpirank) {
    const int num_entries = 2;
    t8_gloidx_t *message;
    /* Allocate a new request */
    sc_MPI_Request *request = (sc_MPI_Request *) sc_array_push (send_requests);

    /* Grow total send buffer */
    sc_array_push_count (send_buffer, num_entries);
    /* Set message to current position in the buffer. */
    message = (t8_gloidx_t *) sc_array_index_int (send_buffer, *current_pos_in_send_buffer);
    *current_pos_in_send_buffer += num_entries;
    /* Send the global id of the first tree of this process */
    message[0] = global_id_of_first_tree;
    /* The index in the tree is the index of the element minus the offset of the tree. */
    message[1] = first_element_in_tree_index_of_current_proc;

    const int mpiret
      = sc_MPI_Isend (message, num_entries, T8_MPI_GLOIDX, iproc, T8_MPI_CMESH_UNIFORM_BOUNDS_START, comm, request);
    SC_CHECK_MPI (mpiret);
#ifdef T8_ENABLE_DEBUG
    (*num_message_sent)++;
#endif
    t8_debugf ("Sending start message (%li, %li) to %li (global num el %li)\n", message[0], message[1], iproc,
               data->global_num_elements);

    T8_ASSERT (proc_is_empty || (0 <= message[0] && message[0] < t8_cmesh_get_num_trees (cmesh)));
    T8_ASSERT (proc_is_empty || (0 <= message[1] && message[1] < data->global_num_elements));
    T8_ASSERT (!(proc_is_empty && message[0] != -1));
    T8_ASSERT (!(proc_is_empty && message[1] != -1));
  }
  else { /* We are the current proc, so we just copy the data. */
    *first_local_tree = global_id_of_first_tree;
    if (first_element_in_tree_index_of_current_proc > 0) {
      /* The first tree is shared */
      if (first_tree_shared != NULL) {
        *first_tree_shared = 1;
      }
    }
    else {
      if (first_tree_shared != NULL) {
        *first_tree_shared = 0;
      }
    }
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin = first_element_in_tree_index_of_current_proc;
    }
    /* We do not expect this message from another proc */
    *expect_start_message = false;
#ifdef T8_ENABLE_DEBUG
    (*num_received_start_messages)++;
#endif
  }
}

void
t8_cmesh_uniform_bounds_from_partition (t8_cmesh_t cmesh, t8_gloidx_t local_num_children, const int level,
                                        const t8_scheme *scheme, t8_gloidx_t *first_local_tree,
                                        t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree,
                                        t8_gloidx_t *child_in_tree_end, int8_t *first_tree_shared, sc_MPI_Comm comm)
{
#ifdef T8_ENABLE_DEBUG
  int num_received_start_messages = 0;
  int num_received_end_messages = 0;
  int num_message_sent = 0;
  int num_message_recv = 0;
#endif
  t8_shmem_array_t offset_array;
  t8_shmem_array_init (&offset_array, sizeof (t8_gloidx_t), cmesh->mpisize + 1, comm);
  /* Fill the offset array for each process with the global index of its first element in
   * the uniform partition.
   * (0, l_n_c_0, l_n_c_0 + l_n_c_1, l_n_c_0 + l_n_c_1 + l_n_c_2, ...) */
  t8_shmem_prefix (&local_num_children, offset_array, 1, T8_MPI_GLOIDX, sc_MPI_SUM);

  t8_cmesh_partition_query_t data;
  data.num_procs = cmesh->mpisize;
  /* Get global number of elements */
  data.global_num_elements = t8_shmem_array_get_gloidx (offset_array, cmesh->mpisize);
  SC_CHECK_ABORTF (0 <= data.global_num_elements && data.global_num_elements < T8_GLOIDX_MAX,
                   "Overflow in number of elements.\n");

  /* Compute number of non-shared-trees and the local index of the first non-shared-tree */
  const int first_tree_shared_shift = cmesh->first_tree_shared ? 1 : 0;
  const t8_locidx_t pure_local_trees = cmesh->num_local_trees - first_tree_shared_shift;
  sc_array send_requests;
  bool expect_start_message = true;
  bool expect_end_message = true;
  if (pure_local_trees > 0) {
    /* Compute which trees and elements to send to which process.
     * We skip empty processes. */
    t8_locidx_t igtree = first_tree_shared_shift;
    sc_array_t first_element_tree;
    sc_array_init_size (&first_element_tree, sizeof (t8_gloidx_t), pure_local_trees + 1);

    /* Set the first entry of first_element_tree to the global index of
     * the first element of our first pure local tree. */
    t8_gloidx_t *elem_index_pointer = (t8_gloidx_t *) sc_array_index_int (&first_element_tree, 0);
    *elem_index_pointer = t8_shmem_array_get_gloidx (offset_array, cmesh->mpirank);

    /* Compute the first element in every pure local tree.
     * This array stores for each tree the global element index offset. 
     * Example: 2 local trees, each has 8 Elements. First element index 12: | 12 | 20 | 28 | */
    for (t8_locidx_t itree = 0; itree < pure_local_trees; ++itree, ++igtree) {
      const t8_gloidx_t *first_element_of_tree = (const t8_gloidx_t *) sc_array_index_int (&first_element_tree, itree);
      t8_gloidx_t *first_element_of_next_tree = (t8_gloidx_t *) sc_array_index_int (&first_element_tree, itree + 1);
      const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, igtree);
      /* Set the first element of the next tree by adding the number of element of the current tree. */
      *first_element_of_next_tree = *first_element_of_tree + scheme->count_leaves_from_root (tree_class, level);
    }

    /* Check that the last tree + 1 stores as offset the start of the next process */
    T8_ASSERT (*((t8_gloidx_t *) t8_sc_array_index_locidx (&first_element_tree, pure_local_trees))
               == t8_shmem_array_get_gloidx (offset_array, cmesh->mpirank + 1));

    /* Compute the range of mpiranks to minimize the loop over processes */

    /* Process offset must be 0 when calling t8_cmesh_determine_partition
     * the first times out of sc_array_split. */
    data.process_offset = 0;
    /* Compute the process that will own the first element of our first tree. */

    const t8_gloidx_t send_first = (t8_gloidx_t) t8_cmesh_determine_partition (&first_element_tree, 0, &data);

    /* Compute the process that may own the last element of our first tree. */
    t8_gloidx_t send_last = (t8_gloidx_t) t8_cmesh_determine_partition (&first_element_tree, pure_local_trees, &data);
    if (send_last >= cmesh->mpisize) {
      /* t8_cmesh_determine_partition will return mpisize if we plug in the number of global
       * trees as tree index, which happens on the last process that has data.
       * We need to correct by subtracting 1. */
      send_last = cmesh->mpisize - 1;
    }
    const t8_gloidx_t num_procs_we_send_to = send_last - send_first + 1;
    sc_array_t offset_partition;
    sc_array_init_size (&offset_partition, sizeof (size_t), num_procs_we_send_to);
    /* In array split the 'types' that we compute are the processes we send to,
     * but offset by the first process, so that they start at 0 and end at num_procs_we_send_to. */
    data.process_offset = send_first;
    /* For each process that we send to find the first tree whose first element
     * belongs to this process.
     * These tree indices will be stored in offset_partition. */
    sc_array_split (&first_element_tree, &offset_partition, (size_t) num_procs_we_send_to + 1,
                    t8_cmesh_determine_partition, (void *) &data);

    /* We know: Lowest process and highest process we need to send trees to
       Tree0 Tree1         TreeN
       | ---- | --- | .... | --- |
       send_first                 send_last

       We need the information: Given a process in range, in which tree do its elements start
       and in which tree do they end.
       Note: If proc i ends in tree j, then proc i+1 may start in tree j or tree j+1.

       offset_partition:
       At position i for process send_first + i the first local tree whose first element
       belongs to send_first+i.
       If no such tree exists, then the index of the previous process is stored.

       Examples: 
       A                                     B                           C
       Proc 0 needs tree 0 and 1             Proc 0 needs tree 0 and 1      Proc 0 needs tree 0 and 1
       Proc 1 needs tree 1 and 2             Proc 1 needs tree 2            Proc 1 needs tree 1
       Proc 2 needs tree 1 and 2
       | 0 | 2 | 3 |                         | 0 | 2 | 3 |                  | 0 | 0 | 2 | 3 |

       Need to identify the first tree we send to proc i:
       If a tree with the first element on proc i exists, then the first tree is either this tree
       or the tree before. We can determine this by comparing the element_offset of the proc with the
       element offset of the tree.
       If no tree exists, then the process only requires elements from one tree (otherwise it would have the
       first element of its second tree).
       This tree must then be the last tree of the previous process.

       Need to identify last tree that a process has elements of.
       Compute last element of process i via ((cmesh->mpirank+1)*data.global_num_elements)/data.num_procs - 1
       Get offset of first tree of process i+1: first_element_tree[offset_partition[i+1]]
       If this index is <= to last element of process i then the last tree of i is the first tree of i+1,
       otherwise this index is > and the last tree of i is the tree before the first tree of i+1.

     */

    /* Initialize the array of MPI Requests */
    sc_array_init (&send_requests, sizeof (sc_MPI_Request));

    /* Initialize the send buffer for all messages.
     * Since we use MPI_Isend, we need to keep the buffers alive until
     * we post the MPI_Wait. Since we first send all messages, we need to
     * buffer them all. */
    sc_array send_buffer;
    sc_array_init (&send_buffer, sizeof (t8_gloidx_t));
    int current_pos_in_send_buffer = 0;

    /* Iterate over offset_partition to find boundaries
     * and send the MPI messages. */

    const t8_gloidx_t last_el_index_of_last_tree
      = *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, pure_local_trees) - 1;
    for (t8_gloidx_t iproc = send_first; iproc <= send_last; iproc++) {

      const t8_gloidx_t first_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (iproc, data.num_procs, data.global_num_elements);
      const t8_gloidx_t last_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (iproc + 1, data.num_procs, data.global_num_elements) - 1;
      const int proc_is_empty = last_element_index_of_current_proc < first_element_index_of_current_proc;
      bool send_start_message = true;
      bool send_end_message = true;

      t8_locidx_t first_puretree_of_current_proc = -1;
      t8_locidx_t last_puretree_of_current_proc = -1;
      t8_gloidx_t first_element_in_tree_index_of_current_proc = -1;
      t8_gloidx_t last_element_in_tree_index_of_current_proc = -2;

      if (!proc_is_empty) {
        /* This process' partition is not empty. */
        const t8_locidx_t possibly_first_puretree_of_current_proc
          = *((size_t *) sc_array_index_int (&offset_partition, iproc - send_first));
        const t8_locidx_t possibly_first_puretree_of_next_proc
          = *((size_t *) sc_array_index_int (&offset_partition, iproc + 1 - send_first));
        const t8_gloidx_t first_el_index_of_first_tree
          = *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, possibly_first_puretree_of_current_proc);
        if (first_element_index_of_current_proc >= last_el_index_of_last_tree + 1) {
          /* We do not send to this process at all. Its first element belongs 
           * to the next process. */
          send_start_message = send_end_message = false;
          first_puretree_of_current_proc = -1;
          last_puretree_of_current_proc = -1;
        }
        if (send_start_message) {
          /* Compute the first tree of this proc and whether we need to send a start message to it. */
          if (first_el_index_of_first_tree > first_element_index_of_current_proc) {
            /* The first element of this proc does lie on possibly_first_puretree_of_current_proc - 1.
             * We check whether we own this tree and if not we do not send anything. */
            if (possibly_first_puretree_of_current_proc - 1 < 0) {
              /* We do not send any Start message to the proc. */
              send_start_message = false;
              first_puretree_of_current_proc = -1;
            }
            else {
              /* The first element of this proc lies on the previous tree. */
              first_puretree_of_current_proc = possibly_first_puretree_of_current_proc - 1;
            }
          }
          else {
            first_puretree_of_current_proc = possibly_first_puretree_of_current_proc;
          }
        }
        /* Compute the last tree of this proc and whether we need to send an end message to it. */
        /* We know 
         * 
         * possibly_first_puretree_of_next_proc - The tree whose first element lies on the next process.
         * 
         * 
         * If the next process is empty, then possibly_first_puretree_of_next_proc = possibly_first_puretree_of_current_proc
         * and this is our last tree.
         */
        if (send_end_message) {
          if (last_el_index_of_last_tree < last_element_index_of_current_proc) {
            /* The last element of this proc does not lie in our partition.
             * We do not need to send any End information to this process. */
            send_end_message = false;
            last_puretree_of_current_proc = -1;
          }
          else {
            if (iproc == send_last) {
              /* The very last process must have our last tree as its last tree. */
              last_puretree_of_current_proc = pure_local_trees - 1;
            }
            else if (possibly_first_puretree_of_next_proc == possibly_first_puretree_of_current_proc) {
              /* The next process is empty. This can only happen if
               *   each process gets 0 or 1 element.
               * Hence, the current process has only one element and it must lie on the first tree. */
              last_puretree_of_current_proc = first_puretree_of_current_proc;
            }
            else {
              last_puretree_of_current_proc = possibly_first_puretree_of_next_proc - 1;
            }
          }
        }
#ifdef T8_ENABLE_DEBUG
        /* Check that the trees have valid values. */
        if (send_start_message) {
          T8_ASSERT (0 <= first_puretree_of_current_proc && first_puretree_of_current_proc < pure_local_trees);
        }
        else {
          T8_ASSERT (first_puretree_of_current_proc == -1);
        }

        if (send_end_message) {
          T8_ASSERT (0 <= last_puretree_of_current_proc && last_puretree_of_current_proc < pure_local_trees);
        }
        else {
          T8_ASSERT (last_puretree_of_current_proc == -1);
        }
        if (first_puretree_of_current_proc != -1 && last_puretree_of_current_proc != -1) {
          T8_ASSERT (first_puretree_of_current_proc <= last_puretree_of_current_proc);
        }
#endif
        if (send_start_message) {
          /* Compute the index int8_cmesh_uniform_bounds_equal_element_countside the tree of the first element. */
          const t8_gloidx_t first_el_of_first_tree
            = *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, first_puretree_of_current_proc);
          first_element_in_tree_index_of_current_proc = first_element_index_of_current_proc - first_el_of_first_tree;
        }
        if (send_end_message) {
          /* Compute the index inside the tree of the last element. */
          const t8_gloidx_t first_el_of_last_tree
            = *(t8_gloidx_t *) t8_sc_array_index_locidx (&first_element_tree, last_puretree_of_current_proc);
          last_element_in_tree_index_of_current_proc = last_element_index_of_current_proc - first_el_of_last_tree;
        }
      }
      else {
        /* This process is empty. */
        send_start_message = send_end_message = false;
        first_puretree_of_current_proc = -1;
        last_puretree_of_current_proc = -1;
      }

      /*
       *
       *  MPI Communication starts here
       * 
       */

      /* Post the start message, we send the first tree id of the process
       * and (if desired) the index of the first element in this tree. */
      if (send_start_message) {
        t8_cmesh_bounds_send_start (cmesh, proc_is_empty, first_puretree_of_current_proc, first_tree_shared_shift,
                                    iproc, &send_requests, &send_buffer, &current_pos_in_send_buffer,
                                    first_element_in_tree_index_of_current_proc, first_local_tree, first_tree_shared,
                                    child_in_tree_begin, &expect_start_message, &data, comm
#ifdef T8_ENABLE_DEBUG
                                    ,
                                    &num_received_start_messages, &num_message_sent

#endif
        );
      } /* End sending of start message */

      t8_debugf ("End message: %i\n", send_end_message);
      /* Post the end message, we send the last tree id of the process
       * and the index of the last element in this tree. */
      if (send_end_message) {
        const t8_gloidx_t global_id_of_last_tree
          = proc_is_empty ? -1
                          : last_puretree_of_current_proc + t8_cmesh_get_first_treeid (cmesh) + first_tree_shared_shift;
        if (iproc != cmesh->mpirank) {
          const int num_entries = 2;
          t8_gloidx_t *message;
          /* Allocate a new request */
          sc_MPI_Request *request = (sc_MPI_Request *) sc_array_push (&send_requests);

          /* Grow total send buffer */
          sc_array_push_count (&send_buffer, num_entries);
          /* Set message to current position in the buffer. */
          message = (t8_gloidx_t *) sc_array_index_int (&send_buffer, current_pos_in_send_buffer);
          current_pos_in_send_buffer += num_entries;

          message[0] = global_id_of_last_tree;
          /* The index in the tree is the index of the element minus the offset of the tree. */
          message[1] = last_element_in_tree_index_of_current_proc;
          const int mpiret
            = sc_MPI_Isend (message, num_entries, T8_MPI_GLOIDX, iproc, T8_MPI_CMESH_UNIFORM_BOUNDS_END, comm, request);
          SC_CHECK_MPI (mpiret);
#ifdef T8_ENABLE_DEBUG
          num_message_sent++;
#endif
          t8_debugf ("Sending end message (%li, %li) to %li\n", message[0], message[1], iproc);
          T8_ASSERT (proc_is_empty || (0 <= message[0] && message[0] < t8_cmesh_get_num_trees (cmesh)));
          T8_ASSERT (proc_is_empty || (0 <= message[1] && message[1] < data.global_num_elements));
          T8_ASSERT (!(proc_is_empty && message[0] != -1));
          T8_ASSERT (!(proc_is_empty && message[1] != -2));
        }
        else { /* We are the current proc, so we just copy the data. */
          *last_local_tree = global_id_of_last_tree;
          if (child_in_tree_end != NULL) {
            *child_in_tree_end = last_element_in_tree_index_of_current_proc + 1;
          }
          /* We do not expect this message from another proc */
          expect_end_message = false;
#ifdef T8_ENABLE_DEBUG
          num_received_end_messages++;
#endif
        }
      } /* End sending of end message */
    }   /* End loop over processes */
    sc_array_reset (&first_element_tree);
    sc_array_reset (&offset_partition);
    sc_array_reset (&send_buffer);
  } /* if (pure_local_trees > 0) */

  const t8_gloidx_t first_element_index_of_current_proc
    = t8_cmesh_get_first_element_of_process (cmesh->mpirank, data.num_procs, data.global_num_elements);
  const t8_gloidx_t last_element_index_of_current_proc
    = t8_cmesh_get_first_element_of_process (cmesh->mpirank + 1, data.num_procs, data.global_num_elements) - 1;

  if (first_element_index_of_current_proc > last_element_index_of_current_proc) {
    /* We do not expect a start/end message if this proc is empty. */
    expect_start_message = false;
    expect_end_message = false;
    *first_local_tree = 0;
    *last_local_tree = -1;
    if (first_tree_shared != NULL) {
      *first_tree_shared = 0;
    }
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin = -1;
    }
    if (child_in_tree_end != NULL) {
      *child_in_tree_end = -1;
    }
#ifdef T8_ENABLE_DEBUG
    num_received_end_messages++;
    num_received_start_messages++;
#endif
  }

  /* Post the receives. */
  if (expect_start_message) {
    t8_debugf ("expect start message\n");
    const int num_entries = 2;
    t8_gloidx_t *message = T8_ALLOC (t8_gloidx_t, num_entries);
    const int mpiret = sc_MPI_Recv (message, num_entries, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE,
                                    T8_MPI_CMESH_UNIFORM_BOUNDS_START, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);

#ifdef T8_ENABLE_DEBUG
    num_message_recv++;
    num_received_start_messages++;
#endif
    /* Copy the received data to output parameters */
    *first_local_tree = message[0];
    T8_ASSERT (*first_local_tree == -1
               || (0 <= *first_local_tree && *first_local_tree < t8_cmesh_get_num_trees (cmesh)));
    if (message[1] > 0) {
      /* The first tree is shared */
      if (first_tree_shared != NULL) {
        *first_tree_shared = 1;
      }
    }
    else {
      if (first_tree_shared != NULL) {
        *first_tree_shared = 0;
      }
    }
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin = message[1];
      T8_ASSERT (*child_in_tree_begin == -1
                 || (0 <= *child_in_tree_begin && *child_in_tree_begin < data.global_num_elements));
    }
    T8_FREE (message);
  } /* End receiving start message */
  if (expect_end_message) {
    const int num_entries = 2;
    t8_gloidx_t *message = T8_ALLOC (t8_gloidx_t, num_entries);
    const int mpiret = sc_MPI_Recv (message, num_entries, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE,
                                    T8_MPI_CMESH_UNIFORM_BOUNDS_END, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
#ifdef T8_ENABLE_DEBUG
    num_received_end_messages++;
    num_message_recv++;
#endif
    t8_debugf ("Receiving end message (%li, %li) global num %li\n", message[0], message[1], data.global_num_elements);

    /* Copy the received data to output parameters */
    *last_local_tree = message[0];
    T8_ASSERT (*last_local_tree == -1 || (0 <= *last_local_tree && *last_local_tree <= t8_cmesh_get_num_trees (cmesh)));
    if (child_in_tree_end != NULL) {
      *child_in_tree_end = message[1] + 1;
      T8_ASSERT (*child_in_tree_end == -1
                 || (0 <= *child_in_tree_end && *child_in_tree_end <= data.global_num_elements));
    }
    T8_FREE (message);
  } /* End receiving end message */
  if (child_in_tree_begin != NULL && child_in_tree_end != NULL) {
    /* Check for empty partition */
    if (*child_in_tree_end < 0 || *child_in_tree_begin < 0) {
      /* This partition is empty */
      t8_debugf ("[H] EMPTY hybrid %li %li %li %li\n", *first_local_tree, *last_local_tree, *child_in_tree_begin,
                 *child_in_tree_end);

      t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                       child_in_tree_end, first_tree_shared);
    }
  }
  t8_gloidx_t num_messages_sent = 0;
  int mpiret = sc_MPI_Waitall (num_messages_sent, (sc_MPI_Request *) send_requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (num_received_start_messages == 1);
  T8_ASSERT (num_received_end_messages == 1);

  if (pure_local_trees > 0) {
    sc_array_reset (&send_requests);
  }

  t8_shmem_array_destroy (&offset_array);

#ifdef T8_ENABLE_DEBUG
  int total_num_sent = 0;
  int total_num_recv = 0;

  mpiret = sc_MPI_Allreduce (&num_message_sent, &total_num_sent, 1, sc_MPI_INT, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (&num_message_recv, &total_num_recv, 1, sc_MPI_INT, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (total_num_sent == total_num_recv);
#endif

  t8_debugf ("Done with t8_cmesh_uniform_bounds_hybrid.\n");
  return;
}

/* TODO: Shared trees, binary search in offset-array to avoid recv_any,
 * use partition_given to partition the cmesh*/
void
t8_cmesh_uniform_bounds_for_irregular_refinement (const t8_cmesh_t cmesh, const int level, const t8_scheme *scheme,
                                                  t8_gloidx_t *first_local_tree, t8_gloidx_t *child_in_tree_begin,
                                                  t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                                                  int8_t *first_tree_shared, sc_MPI_Comm comm)
{
  T8_ASSERT (cmesh != NULL);
  /* TODO: Clean up size_t and gloidx_t data types, ensure that each variables has the 
   *          matching type. */

  t8_debugf ("Into t8_cmesh_uniform_bounds_hybrid.\n");

  if (t8_cmesh_is_empty (cmesh)) {
    t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                     child_in_tree_end, first_tree_shared);
    return;
  }

#ifdef T8_ENABLE_DEBUG
  {
    /* Check that comm matches mpirank and size stored in cmesh */
    int mpirank, mpisize;
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (mpirank == cmesh->mpirank);
    T8_ASSERT (mpisize == cmesh->mpisize);
  }
#endif
  t8_gloidx_t local_num_children = 0;
  /*Compute number of local elements. Ignore shared trees */
  for (int ieclass = T8_ECLASS_ZERO; ieclass < T8_ECLASS_COUNT; ieclass++) {
    local_num_children
      += cmesh->num_local_trees_per_eclass[ieclass] * scheme->count_leaves_from_root ((t8_eclass_t) ieclass, level);
  }

  /* Do not consider shared trees */
  if (cmesh->first_tree_shared && cmesh->set_partition) {
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, 0);
    local_num_children -= scheme->count_leaves_from_root (tree_class, level);
  }

  /* 
   *
   *  Cmesh is not partitioned
   * 
   */
  /* If the initial cmesh is not partitioned, every process knows "everything" and we do not
   * need any communication.*/
  if (!cmesh->set_partition) {
    t8_cmesh_uniform_bounds_from_unpartioned (cmesh, local_num_children, level, scheme, first_local_tree,
                                              child_in_tree_begin, last_local_tree, child_in_tree_end,
                                              first_tree_shared);
    return;
  }
  else {
    /* 
    *
    *  Cmesh is partitioned
    * 
    */
    t8_cmesh_uniform_bounds_from_partition (cmesh, local_num_children, level, scheme, first_local_tree,
                                            child_in_tree_begin, last_local_tree, child_in_tree_end, first_tree_shared,
                                            comm);
    return;
  }
}
