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

#include <t8_refcount.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

/** \file t8_cmesh.h
 *
 * TODO: document this file
 */

/* *INDENT-OFF* */
static int
t8_cmesh_tree_id_is_owned (t8_cmesh_t cmesh, t8_topidx_t tree_id);

static t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_topidx_t tree_id);
/* *INDENT-ON* */

/* Compute a hash value for a ghost tree. */
static unsigned
t8_cmesh_ghost_hash_fn (const void *ghost, const void *data)
{
  t8_cmesh_t          cmesh;
  t8_cghost_t         G;

  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (cmesh->num_ghosts > 0);
  T8_ASSERT (cmesh->set_partitioned);
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

void
t8_cmesh_init (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (pcmesh != NULL);

  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  t8_refcount_init (&cmesh->rc);

  /* sensible (hard error) defaults */
  cmesh->dimension = -1;
  cmesh->mpicomm = sc_MPI_COMM_WORLD;
  cmesh->mpirank = -1;
  cmesh->mpisize = -1;
}

void
t8_cmesh_set_mpicomm (t8_cmesh_t cmesh, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount > 0);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->mpicomm == sc_MPI_COMM_WORLD);

  T8_ASSERT (mpicomm != sc_MPI_COMM_NULL);

  cmesh->mpicomm = mpicomm;
  cmesh->do_dup = do_dup;
}

sc_MPI_Comm
t8_cmesh_get_mpicomm (t8_cmesh_t cmesh, int *do_dup)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount > 0);
  T8_ASSERT (cmesh->mpicomm != sc_MPI_COMM_NULL);

  *do_dup = cmesh->do_dup;
  return cmesh->mpicomm;
}

void
t8_cmesh_set_partitioned (t8_cmesh_t cmesh, int set_partitioned,
                          t8_topidx_t num_global_trees,
                          t8_topidx_t first_local_tree,
                          t8_topidx_t num_ghosts)
{
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->set_partitioned == 0);
  T8_ASSERT (cmesh->num_trees == 0);
  T8_ASSERT (cmesh->num_local_trees == 0);
  T8_ASSERT (cmesh->first_tree == 0);

  if ((cmesh->set_partitioned = set_partitioned) == 0) {
    /* The mesh is replicated, and this function just serves
     * as set_num_trees.
     * first_local_tree and num_ghosts are ignored. */
    t8_cmesh_set_num_trees (cmesh, num_global_trees);
    return;
  }
  else {
    cmesh->num_trees = num_global_trees;
    cmesh->first_tree = first_local_tree;
    cmesh->num_ghosts = num_ghosts;
    cmesh->ghosts = sc_hash_array_new (sizeof (t8_cghost_struct_t),
                                       t8_cmesh_ghost_hash_fn,
                                       t8_cmesh_ghost_equal_fn,
                                       (void *) cmesh);
  }
}

/* Return a pointer to the ctree of a given tree_id. */
static              t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  t8_topidx_t         index;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));

  index = cmesh->set_partitioned ? tree_id - cmesh->first_tree : tree_id;
  return (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, index);
}

/** Set the number of corners for a cmesh.
 * If there are no corners (num_corners = 0) this function does not need to
 * be called.
 * It is not allowed to call this function after \see t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     num_corners The number of corners to be set.
 */
static void
t8_cmesh_set_num_corners (t8_cmesh_t cmesh, t8_topidx_t num_corners)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (num_corners > 0);
  T8_ASSERT (cmesh->num_corners == 0);

  cmesh->num_corners = num_corners;
  if (cmesh->set_partitioned == 0) {
    cmesh->num_local_corners = num_corners;
  }
}

/** Set the number of local corners for a partitioned cmesh.
*  If there are no local corners (num_local_corners = 0) this function does not need to
 * be called.
 * It is not allowed to call this function after \see t8_cmesh_commit.
 * Must be called after \ref t8_cmesh_set_num_corners.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     num_local_corners The number of local corners to be set.
 */
static void
t8_cmesh_set_num_local_corners (t8_cmesh_t cmesh,
                                t8_topidx_t num_local_corners)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (num_local_corners > 0);
  T8_ASSERT (cmesh->num_corners >= num_local_corners);
  T8_ASSERT (cmesh->set_partitioned);

  cmesh->num_local_corners = num_local_corners;
}

void
t8_cmesh_set_attribute_sizes (t8_cmesh_t cmesh, size_t attr_sizes[],
                              int num_sizes)
{
  int                 iclass;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (num_sizes == T8_ECLASS_LAST);

  for (iclass = 0; iclass < num_sizes; iclass++) {
    cmesh->tree_attributes_mem[iclass] = sc_mempool_new (attr_sizes[iclass]);
  }
}

static              size_t
t8_cmesh_get_attribute_size (t8_cmesh_t cmesh, t8_eclass_t eclass)
{
  T8_ASSERT (cmesh->tree_attributes_mem[eclass] != NULL);
  return cmesh->tree_attributes_mem[eclass]->elem_size;
}

void
t8_cmesh_set_attribute_size_single (t8_cmesh_t cmesh, size_t attr_size,
                                    t8_eclass_t tree_class)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (T8_ECLASS_FIRST <= tree_class && tree_class < T8_ECLASS_LAST);

  cmesh->tree_attributes_mem[(int) tree_class] = sc_mempool_new (attr_size);
}

void
t8_cmesh_tree_set_attribute (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                             void *attribute)
{
  t8_ctree_t          tree;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));

  tree = t8_cmesh_get_tree (cmesh, tree_id);
  T8_ASSERT (tree->eclass != T8_ECLASS_LAST);
  T8_ASSERT (tree->attribute == NULL);
  T8_ASSERT (t8_cmesh_get_attribute_size (cmesh, tree->eclass) > 0);

  tree->attribute =
    sc_mempool_alloc (cmesh->tree_attributes_mem[tree->eclass]);
  memcpy (tree->attribute, attribute,
          t8_cmesh_get_attribute_size (cmesh, tree->eclass));
}

void
t8_cmesh_set_num_trees (t8_cmesh_t cmesh, t8_topidx_t num_trees)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  /* If the cmesh is entered as a partitioned cmesh,
   * this function sets the local number of trees;
   * the global number then must have been set in cmesh_set_partitioned.
   * Otherwise the global number of trees is set here.
   */
  if (cmesh->set_partitioned) {
    /* num_trees == 0 is allowed */
    T8_ASSERT (cmesh->num_trees > 0);
    T8_ASSERT (cmesh->num_local_trees == 0);
    cmesh->num_local_trees = num_trees;
  }
  else {
    /* num_trees == 0 is not allowed */
    T8_ASSERT (num_trees > 0);
    T8_ASSERT (cmesh->num_trees == 0);
    cmesh->num_trees = cmesh->num_local_trees = num_trees;
  }
  /* As soon as we know the number of trees, we allocate
   * the ctree array.
   */
  cmesh->ctrees = sc_array_new_size (sizeof (t8_ctree_struct_t), num_trees);
}

/* Check whether a given tree_id belongs to a tree in the cmesh.
 * If partitioned only local trees are allowed.
 */
static int
t8_cmesh_tree_id_is_owned (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  if (cmesh->set_partitioned) {
    return cmesh->first_tree <= tree_id
      && tree_id < cmesh->first_tree + cmesh->num_local_trees;
  }
  else {
    return 0 <= tree_id && tree_id < cmesh->num_trees;
  }
}

/* Given a tree_id return the index of the specified tree in
 * cmesh's tree array
 */
static              t8_topidx_t
t8_cmesh_tree_index (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  return cmesh->set_partitioned ? tree_id - cmesh->first_tree : tree_id;
}

void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                         t8_eclass_t tree_class)
{
  t8_ctree_t          tree;
  int                 i, num_neighbors;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));

  /* If we insert the first tree, set the dimension of the cmesh
   * to this tree's dimension. Otherwise check whether the dimension
   * of the tree to be inserted equals the dimension of the cmesh. */
  if (cmesh->dimension == -1) {
    cmesh->dimension = t8_eclass_to_dimension[tree_class];
  }
  else {
    T8_ASSERT (t8_eclass_to_dimension[tree_class] == cmesh->dimension);
  }
  cmesh->num_trees_per_eclass[tree_class]++;

  tree = t8_cmesh_get_tree (cmesh, tree_id);

  tree->eclass = tree_class;
  tree->treeid = tree_id;
  num_neighbors = t8_eclass_num_faces[tree_class];
  /* Allocate neighbors and set entries to invalid values. */
  tree->face_neighbors =
    T8_ALLOC (t8_ctree_fneighbor_struct_t, num_neighbors);
  for (i = 0; i < num_neighbors; i++) {
    tree->face_neighbors[i].treeid = -1;
    tree->face_neighbors[i].tree_to_face = -1;
    tree->face_neighbors[i].is_owned = -1;
  }
  tree->corners = NULL;
  tree->vertices = NULL;
  tree->attribute = NULL;
#ifdef T8_ENABLE_DEBUG
  cmesh->inserted_trees++;
#endif
}

/** Set the corners of a tree in the cmesh.
 * It is not allowed to call this function after \see t8_cmesh_commit.
 * The eclass of the tree has to be set before calling this function.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree_id      The global number of the tree.
 * \param [in]     corners     An array of as many corner indices as the tree
 *                              has corners.
 * \param [in]     num_corners The number of corners in \a corners. Must
 *                              match the number of corners of the tree.
 */
static void
t8_cmesh_set_tree_corners (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                           t8_topidx_t * corners, t8_topidx_t num_corners)
{
  t8_ctree_t          tree;
  int                 vi;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));
  T8_ASSERT (corners != NULL);

  tree = t8_cmesh_get_tree (cmesh, tree_id);
  T8_ASSERT (tree->eclass != T8_ECLASS_LAST);
  T8_ASSERT (num_corners == t8_eclass_num_vertices[tree->eclass]);
  T8_ASSERT (tree->corners == NULL);

  tree->corners = T8_ALLOC (t8_topidx_t, num_corners);

#ifdef T8_ENABLE_DEBUG
  for (vi = 0; vi < num_corners; vi++) {
    T8_ASSERT (0 <= corners[vi] && corners[vi] < cmesh->num_corners);
    tree->corners[vi] = corners[vi];
  }
#else
  memcpy (tree->corners, corners, num_corners * sizeof (*corners));
#endif
}

void
t8_cmesh_set_attribute_to_vertices (t8_cmesh_t cmesh)
{
  size_t              attribute_sizes[T8_ECLASS_LAST];
  int                 iclass;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  for (iclass = T8_ECLASS_FIRST; iclass < T8_ECLASS_LAST; iclass++) {
    attribute_sizes[iclass] =
      3 * sizeof (double) * t8_eclass_num_vertices[iclass];
  }
  t8_cmesh_set_attribute_sizes (cmesh, attribute_sizes, T8_ECLASS_LAST);
}

void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                            double *vertices, t8_topidx_t num_vertices)
{
  t8_ctree_t          tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));
  T8_ASSERT (vertices != NULL);

  tree = t8_cmesh_get_tree (cmesh, tree_id);
  T8_ASSERT (tree->eclass != T8_ECLASS_LAST);
  T8_ASSERT (num_vertices * 3 * sizeof (double) ==
             t8_cmesh_get_attribute_size (cmesh, tree->eclass));
  T8_ASSERT (tree->attribute == NULL);

  t8_cmesh_tree_set_attribute (cmesh, tree_id, (void *) vertices);
}

void
t8_cmesh_set_ghost (t8_cmesh_t cmesh, t8_topidx_t ghost_id,
                    t8_eclass_t ghost_eclass)
{
  t8_cghost_t         Ghost;
  int                 i;
  void               *check_ret;

  T8_ASSERT (cmesh->set_partitioned);
  T8_ASSERT (0 <= ghost_id && ghost_id < cmesh->num_ghosts);
  /* If we insert the very first tree, set the dimension of the cmesh
   * to this tree's dimension. Otherwise check whether the dimension
   * of the tree to be inserted equals the dimension of the cmesh. */
  if (cmesh->dimension == -1) {
    cmesh->dimension = t8_eclass_to_dimension[ghost_eclass];
  }
  else {
    T8_ASSERT (t8_eclass_to_dimension[ghost_eclass] == cmesh->dimension);
  }
  Ghost = T8_ALLOC (t8_cghost_struct_t, 1);
  Ghost->eclass = ghost_eclass;
  Ghost->treeid = ghost_id;
  Ghost->owning_proc = -1;
  Ghost->local_neighbors = T8_ALLOC (t8_topidx_t,
                                     t8_eclass_num_faces[ghost_eclass]);
  for (i = 0; i < t8_eclass_num_faces[ghost_eclass]; i++) {
    Ghost->local_neighbors[i] = -1;
  }
  check_ret = sc_hash_array_insert_unique (cmesh->ghosts, Ghost, NULL);
  if (check_ret == NULL) {
    SC_ABORTF ("Ghost tree %i inserted twice.", ghost_id);
  }
#ifdef T8_ENABLE_DEBUG
  cmesh->inserted_ghosts++;
#endif
}

void
t8_cmesh_join_faces (t8_cmesh_t cmesh, t8_topidx_t tree1, t8_topidx_t tree2,
                     int face1, int face2, int orientation)
{
  t8_ctree_t          T1, T2;

  T8_ASSERT (0 <= tree1 && tree1 < cmesh->num_trees);
  T8_ASSERT (0 <= tree2 && tree2 < cmesh->num_trees);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree1)
             || t8_cmesh_tree_id_is_owned (cmesh, tree2));      /* At least one of the trees
                                                                 * must belong to this process. */
  T8_ASSERT (0 <= orientation);

  if (t8_cmesh_tree_id_is_owned (cmesh, tree1)
      || t8_cmesh_tree_id_is_owned (cmesh, tree2))
    /* Both trees belong to this process. */
  {
    T1 = t8_cmesh_get_tree (cmesh, tree1);
    T2 = t8_cmesh_get_tree (cmesh, tree2);
    /* Check if the trees were added to cmesh before. */
    T8_ASSERT (T1->treeid == tree1 && T2->treeid == tree2);
    T8_ASSERT (0 <= face1 && face1 < t8_eclass_num_faces[T1->eclass]);
    T8_ASSERT (0 <= face2 && face2 < t8_eclass_num_faces[T2->eclass]);
    /* Check if both faces are of the same type (i.e. do not join a triangle and a square) */
    T8_ASSERT (t8_eclass_face_types[T1->eclass][face1] ==
               t8_eclass_face_types[T2->eclass][face2]);
    /* Check if both faces are not joined already */
    T8_ASSERT (T1->face_neighbors[face1].treeid == -1);
    T8_ASSERT (T2->face_neighbors[face2].treeid == -1);
    /* Compute the tree_to_face index according to the given orientation. */
    T1->face_neighbors[face1].is_owned = 1;
    T1->face_neighbors[face1].treeid = tree2;
    T1->face_neighbors[face1].tree_to_face = orientation *
      t8_eclass_num_faces[T2->eclass] + face2;
    T2->face_neighbors[face2].is_owned = 1;
    T2->face_neighbors[face2].treeid = tree1;
    T2->face_neighbors[face2].tree_to_face = orientation *
      t8_eclass_num_faces[T1->eclass] + face1;
  }
  else
    /* One of the trees is not owned by this process. */
  {
    t8_topidx_t         ghost_id;
    t8_topidx_t         owned_id;
#ifdef T8_ENABLE_DEBUG
    int                 owned_face;
#endif
    int                 ghost_face;
    t8_cghost_t         Ghost;
    t8_eclass_t         ghost_eclass;
    size_t              pos;

    /* Find out which one is owned and which one not. */
    if (t8_cmesh_tree_id_is_owned (cmesh, tree1)) {
      owned_id = tree1;
      ghost_id = tree2;
#ifdef T8_ENABLE_DEBUG
      owned_face = face1;
#endif
      ghost_face = face2;
    }
    else {
      T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree2));
      owned_id = tree2;
      ghost_id = tree1;
#ifdef T8_ENABLE_DEBUG
      owned_face = face2;
#endif
      ghost_face = face1;
    }
    T1 = t8_cmesh_get_tree (cmesh, owned_id);

    Ghost = T8_ALLOC (t8_cghost_struct_t, 1);
    Ghost->treeid = ghost_id;
    if (sc_hash_array_insert_unique (cmesh->ghosts, Ghost, &pos) == NULL)
      /* The ghost already exists in the array and we only need to add data to it. */
    {
      T8_FREE (Ghost);
      Ghost = (t8_cghost_t) sc_array_index (&cmesh->ghosts->a, pos);
      ghost_eclass = Ghost->eclass;
      T8_ASSERT (Ghost->treeid == ghost_id);
    }
    else {
      SC_ABORTF ("The ghost tree %i was not found in the coarse mesh.",
                 ghost_id);
    }
    /* Check if both faces are of the same type (i.e. do not join a triangle and a square) */
    T8_ASSERT (t8_eclass_face_types[T1->eclass][owned_face] ==
               t8_eclass_face_types[ghost_eclass][ghost_face]);
    /* Check if both faces are not joined already */
    T8_ASSERT (T1->face_neighbors[face1].treeid == -1);
    T8_ASSERT (Ghost->local_neighbors[ghost_face] == -1);

    Ghost->local_neighbors[ghost_face] = owned_id;
    /* Compute the tree_to_face index according to the given orientation. */
    T1->face_neighbors[ghost_face].is_owned = 0;
    T1->face_neighbors[ghost_face].treeid = ghost_id;
    T1->face_neighbors[ghost_face].tree_to_face = orientation *
      t8_eclass_num_faces[ghost_eclass] + ghost_face;
  }
}

/* compare two arrays of face_neighbors for equality */
static int
t8_cmesh_face_n_is_equal (t8_ctree_fneighbor_struct_t * face_a,
                          t8_ctree_fneighbor_struct_t * face_b, int num_neigh)
{
  int                 iface;

  for (iface = 0; iface < num_neigh; iface++) {
    if (face_a[iface].is_owned != face_b[iface].is_owned ||
        face_a[iface].treeid != face_b[iface].treeid ||
        face_a[iface].tree_to_face != face_b[iface].tree_to_face) {
      return 0;
    }
  }
  return 1;
}

static int
t8_cmesh_ctree_is_equal (t8_ctree_t tree_a, t8_ctree_t tree_b)
{
  int                 is_equal;
  T8_ASSERT (tree_a != NULL && tree_b != NULL);

  is_equal = tree_a->treeid != tree_b->treeid ||
    tree_a->eclass != tree_b->eclass;
  if (is_equal != 0) {
    return 0;
  }
  if (tree_a->corners != NULL) {
    if (tree_b->corners == NULL) {
      return 0;
    }
    else {
      if (memcmp (tree_a->corners, tree_b->corners,
                  t8_eclass_num_vertices[tree_a->eclass] *
                  sizeof (t8_topidx_t))) {
        return 0;
      }
    }
  }
  if (!t8_cmesh_face_n_is_equal
      (tree_a->face_neighbors, tree_b->face_neighbors,
       t8_eclass_num_faces[tree_a->eclass])) {
    return 0;
  }

  /* TODO check attributes */
  return 1;
}

/* returns true if cmesh_a equals cmesh_b */
int
t8_cmesh_is_equal (t8_cmesh_t cmesh_a, t8_cmesh_t cmesh_b)
{
  int                 is_equal;
  int                 iclass;
  t8_topidx_t         itree;
  T8_ASSERT (cmesh_a != NULL && cmesh_b != NULL);

  if (cmesh_a == cmesh_b) {
    return 1;
  }
  /* check entries that are numbers */
  is_equal = cmesh_a->committed != cmesh_b->committed || cmesh_a->dimension !=
    cmesh_b->dimension || cmesh_a->do_dup != cmesh_b->do_dup ||
    cmesh_a->set_partitioned != cmesh_b->set_partitioned ||
    cmesh_a->mpirank != cmesh_b->mpirank ||
    cmesh_a->mpisize != cmesh_b->mpisize ||
    cmesh_a->num_corners != cmesh_b->num_corners ||
    cmesh_a->num_local_corners != cmesh_b->num_local_corners ||
    cmesh_a->num_trees != cmesh_b->num_trees ||
    cmesh_a->num_local_trees != cmesh_b->num_local_trees ||
    cmesh_a->num_ghosts != cmesh_b->num_ghosts ||
    cmesh_a->first_tree != cmesh_b->first_tree;
#ifdef T8_ENABLE_DEBUG
  is_equal = is_equal || cmesh_a->inserted_trees != cmesh_b->inserted_trees ||
    cmesh_a->inserted_ghosts != cmesh_b->inserted_ghosts;
#endif
  if (cmesh_a->do_dup == 0) {
    is_equal = is_equal || cmesh_a->mpicomm != cmesh_b->mpicomm;
  }
  if (is_equal != 0) {
    return 0;
  }
  /* check arrays */
  is_equal = memcmp (cmesh_a->num_trees_per_eclass,
                     cmesh_b->num_trees_per_eclass,
                     T8_ECLASS_LAST * sizeof (t8_topidx_t));
  /* check attribute sizes */
  for (iclass = 0; iclass < T8_ECLASS_LAST; iclass++) {
    if (cmesh_a->tree_attributes_mem[iclass] != NULL) {
      if (cmesh_b->tree_attributes_mem[iclass] == NULL) {
        return 0;
      }
      else {
        is_equal = is_equal
          || t8_cmesh_get_attribute_size (cmesh_a,
                                          (t8_eclass_t ) iclass) !=
          t8_cmesh_get_attribute_size (cmesh_b, (t8_eclass_t ) iclass)
          || sc_mempool_memory_used (cmesh_a->tree_attributes_mem[iclass]) !=
          sc_mempool_memory_used (cmesh_b->tree_attributes_mem[iclass]);
      }
    }
  }
  /* check tree_offsets */
  if (cmesh_a->tree_offsets != NULL) {
    if (cmesh_b->tree_offsets == NULL) {
      return 0;
    }
    else {
      is_equal = is_equal || memcmp (cmesh_a->tree_offsets,
                                     cmesh_b->tree_offsets,
                                     cmesh_a->mpisize * sizeof (t8_topidx_t));
    }
  }
  if (is_equal != 0) {
    return 0;
  }
  /* check trees */
  for (itree = 0; itree < cmesh_a->num_trees; itree++) {
    if (!t8_cmesh_ctree_is_equal (t8_cmesh_get_tree (cmesh_a, itree),
                                  t8_cmesh_get_tree (cmesh_b, itree))) {
      return 0;
    }
  }
  /* check ghosts */
  if (cmesh_a->num_ghosts > 0 &&
      !sc_array_is_equal (&cmesh_a->ghosts->a, &cmesh_b->ghosts->a)) {
    return 0;
  }
  return 1;
}

t8_cmesh_t
t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int                 mpirank, mpisize, mpiret;
  int                 iclass;
  t8_ctree_t          tree;
  t8_topidx_t         itree;
  t8_topidx_t         count_face;
  t8_topidx_t         num_neighbors;
  t8_ctree_fneighbor_struct_t *fneighbors;

  struct
  {
    int                 dimension;
    int                 do_dup;
    t8_topidx_t         num_trees;
    t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST];
    size_t              tree_attribute_size[T8_ECLASS_LAST];
#ifdef T8_ENABLE_DEBUG
    t8_topidx_t         inserted_trees;
#endif
  } dimensions;

  /* TODO: BUG: running with two processes and a cmesh of one T8_ECLASS_LINE,
   *       the on both processes the face_neigbors and vertices arrays of
   *       the single tree point to the same physical memory.
   *       (face_neighbors on both processes are equal and vertices on both
   *        processes are equal)
   */
  /* TODO: Send the tree's vertices */

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (0 <= root && root < mpisize);
  T8_ASSERT (mpirank == root || cmesh_in == NULL);
  T8_ASSERT (mpirank != root || cmesh_in != NULL);
  T8_ASSERT (mpirank != root || cmesh_in->mpicomm == comm);
  T8_ASSERT (mpirank != root || cmesh_in->set_partitioned == 0);
  /* The cmesh on the calling process must not be owned by something
   * else. */
  T8_ASSERT (mpirank != root || cmesh_in->rc.refcount == 1);

  /* At first we broadcast all information needed to allocate the tree
   * arrays. */
  if (mpirank == root) {
    /* TODO: vertices are missing yet */
    dimensions.dimension = cmesh_in->dimension;
    dimensions.do_dup = cmesh_in->do_dup;
    dimensions.num_trees = cmesh_in->num_trees;
    for (iclass = 0; iclass < T8_ECLASS_LAST; iclass++) {
      dimensions.num_trees_per_eclass[iclass] =
        cmesh_in->num_trees_per_eclass[iclass];
      dimensions.tree_attribute_size[iclass] =
        t8_cmesh_get_attribute_size (cmesh_in, (t8_eclass_t ) iclass);
    }
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
    cmesh_in->mpicomm = comm;
    cmesh_in->dimension = dimensions.dimension;
    cmesh_in->do_dup = dimensions.do_dup;
    /* set tree num and allocate trees */
    t8_cmesh_set_num_trees (cmesh_in, dimensions.num_trees);
    for (iclass = 0; iclass < T8_ECLASS_LAST; iclass++) {
      cmesh_in->num_trees_per_eclass[iclass] =
        dimensions.num_trees_per_eclass[iclass];
    }
    t8_cmesh_set_attribute_sizes (cmesh_in, dimensions.tree_attribute_size,
                                  T8_ECLASS_LAST);
#ifdef T8_ENABLE_DEBUG
    cmesh_in->inserted_trees = dimensions.inserted_trees;
#endif
  }
  /* broadcast all the trees */
  /* TODO: this step relies on the sc_array implementation.
   *       can we do it differently ? */
  mpiret = sc_MPI_Bcast (cmesh_in->ctrees->array,
                         cmesh_in->num_trees * sizeof (t8_ctree_struct_t),
                         sc_MPI_BYTE, root, comm);
  SC_CHECK_MPI (mpiret);
  if (mpirank != root) {
    /* iterate through trees and allocate neighbor and corner arrays.
     * We cannot do it before because we have to know a tree's eclass for this */
    for (itree = 0; itree < cmesh_in->num_trees; itree++) {
      tree = t8_cmesh_get_tree (cmesh_in, itree);
      tree->face_neighbors = T8_ALLOC (t8_ctree_fneighbor_struct_t,
                                       t8_eclass_num_faces[tree->eclass]);
      tree->vertices = T8_ALLOC (t8_topidx_t,
                                 t8_eclass_num_vertices[tree->eclass]);
      tree->corners = NULL;
      if (tree->attribute != NULL) {
        tree->attribute =
          T8_ALLOC (char, cmesh_in->tree_attribute_size[tree->eclass]);
      }
    }
  }
  /* broadcast attributes */
  for (itree = 0; itree < cmesh_in->num_trees; itree++) {
    tree = t8_cmesh_get_tree (cmesh_in, itree);
    if (tree->attribute != NULL) {
      sc_MPI_Bcast (tree->attribute,
                    cmesh_in->tree_attribute_size[tree->eclass], sc_MPI_BYTE,
                    root, comm);
    }
  }
  /* Since broadcasting one big data set instead of several small ones is much
   * faster, we collect all corner and face neighbor information in arrays and
   * broadcast those.
   */
  /* count total number of face neighbors */
  num_neighbors = 0;
  for (iclass = 0; iclass < T8_ECLASS_LAST; iclass++) {
    num_neighbors += cmesh_in->num_trees_per_eclass[iclass] *
      t8_eclass_num_faces[iclass];
  }
  fneighbors = T8_ALLOC (t8_ctree_fneighbor_struct_t, num_neighbors);
  /* fill face_neighbor arrays on root */
  if (mpirank == 0) {
    count_face = 0;
    for (itree = 0; itree < cmesh_in->num_trees; itree++) {
      tree = t8_cmesh_get_tree (cmesh_in, itree);
      memcpy (fneighbors + count_face, tree->face_neighbors,
              t8_eclass_num_faces[tree->eclass] *
              sizeof (t8_ctree_fneighbor_struct_t));
      count_face += t8_eclass_num_faces[tree->eclass];
    }
    T8_ASSERT (count_face == num_neighbors);
  }

  mpiret = sc_MPI_Bcast (fneighbors, num_neighbors *
                         sizeof (t8_ctree_fneighbor_struct_t),
                         sc_MPI_BYTE, root, comm);
  SC_CHECK_MPI (mpiret);

  if (mpirank != 0) {
    count_face = 0;
    for (itree = 0; itree < cmesh_in->num_trees; itree++) {
      tree = t8_cmesh_get_tree (cmesh_in, itree);
      memcpy (tree->face_neighbors, fneighbors + count_face,
              t8_eclass_num_faces[tree->eclass] *
              sizeof (t8_ctree_fneighbor_struct_t));
      count_face += t8_eclass_num_faces[tree->eclass];
    }
    T8_ASSERT (count_face == num_neighbors);
  }
  T8_FREE (fneighbors);
  return cmesh_in;
}

void
t8_cmesh_commit (t8_cmesh_t cmesh)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->mpicomm != sc_MPI_COMM_NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (cmesh->num_trees > 0);

  T8_ASSERT (cmesh->num_trees == cmesh->inserted_trees);
  T8_ASSERT (cmesh->num_ghosts == cmesh->inserted_ghosts);
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
}

t8_topidx_t
t8_cmesh_get_num_corners (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_corners;
}

t8_topidx_t
t8_cmesh_get_num_vertices (t8_cmesh_t cmesh)
{
  int                 iclass;
  t8_topidx_t         num_vertices = 0;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  for (iclass = T8_ECLASS_FIRST; iclass < T8_ECLASS_LAST; iclass++) {
    num_vertices += t8_eclass_num_vertices[iclass] *
      cmesh->num_trees_per_eclass[iclass];
  }
  return num_vertices;
}

t8_topidx_t
t8_cmesh_get_num_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_trees;
}

t8_topidx_t
t8_cmesh_get_local_num_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  if (cmesh->set_partitioned) {
    return cmesh->num_local_trees;
  }
  else {
    return cmesh->num_trees;
  }
}

t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  t8_ctree_t          tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree_id));

  tree = t8_cmesh_get_tree (cmesh, tree_id);
  return tree->eclass;
}

void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                         t8_topidx_t * first_local_tree,
                         t8_gloidx_t *
                         child_in_tree_begin,
                         t8_topidx_t * last_local_tree,
                         t8_gloidx_t * child_in_tree_end)
{
  *first_local_tree = 0;
  *child_in_tree_begin = 0;
  *last_local_tree = 0;
  *child_in_tree_end = 0;

  if (cmesh->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0) {
    t8_gloidx_t         global_num_children;
    t8_gloidx_t         first_global_child;
    t8_gloidx_t         last_global_child;
    t8_gloidx_t         children_per_tree;
    const t8_gloidx_t   one = 1;

    children_per_tree = one << cmesh->dimension * level;
    global_num_children = cmesh->num_trees * children_per_tree;

    if (cmesh->mpirank == 0) {
      *child_in_tree_begin = first_global_child = 0;
    }
    else {
      /* The first global children of processor p
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
    T8_ASSERT (0 <= first_global_child
               && first_global_child <= global_num_children);
    T8_ASSERT (0 <= last_global_child
               && last_global_child <= global_num_children);
    *first_local_tree = first_global_child / children_per_tree;
    *child_in_tree_begin =
      first_global_child - *first_local_tree * children_per_tree;
    if (first_global_child < last_global_child) {
      *last_local_tree = (last_global_child - 1) / children_per_tree;
    }
    else {
      /* empty processor */
      *last_local_tree = *first_local_tree;
    }
    if (*last_local_tree > 0) {
      *child_in_tree_end =
        last_global_child - *last_local_tree * children_per_tree;
    }
    else {
      *child_in_tree_end = last_global_child;
    }
  }
  else {
    SC_ABORT ("Partition does not support pyramidal elements yet.");
  }
}

static void
t8_cmesh_reset (t8_cmesh_t * pcmesh)
{
  int                 mpiret;
  int                 iclass;
  t8_cmesh_t          cmesh;
  t8_ctree_t          treeit;
  t8_cghost_t         ghostit;
  t8_topidx_t         ti;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  if (cmesh->do_dup && cmesh->committed) {
    mpiret = sc_MPI_Comm_free (&cmesh->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  /* free attributes */
  for (iclass = T8_ECLASS_FIRST; iclass < T8_ECLASS_LAST; iclass++) {
    sc_mempool_destroy (cmesh->tree_attributes_mem[iclass]);
  }
  /* free trees */
  if (cmesh->ctrees != NULL) {
    for (ti = 0; ti < cmesh->num_local_trees; ti++) {
      treeit = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, ti);
      T8_FREE (treeit->face_neighbors);
      T8_FREE (treeit->corners);
      T8_FREE (treeit->vertices);
    }
    sc_array_destroy (cmesh->ctrees);
  }
  /* free tree_offset */
  if (cmesh->tree_offsets != NULL) {
    T8_FREE (cmesh->tree_offsets);
  }
  /* free ghosts */
  if (cmesh->ghosts != NULL) {
    for (ti = 0; ti < cmesh->num_ghosts; ti++) {
      ghostit =
        (t8_cghost_t) t8_sc_array_index_topidx (&cmesh->ghosts->a, ti);
      if (ghostit != NULL) {
        T8_FREE (ghostit->local_neighbors);
        T8_FREE (ghostit);
      }
    }
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

static void
t8_cmesh_new_translate_vertices_to_attributes (t8_topidx_t * tvertices,
                                               double *vertices,
                                               double *attr_vertices,
                                               int num_vertices)
{
  int                 i;

  for (i = 0; i < num_vertices; i++) {
    attr_vertices[3 * i] = vertices[tvertices[i]];
    attr_vertices[3 * i + 1] = vertices[tvertices[i] + 1];
    attr_vertices[3 * i + 2] = vertices[tvertices[i] + 2];
  }
}

t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[9] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_attribute_to_vertices (cmesh);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    1, 1, 1,
    1, -1, -1,
    -1, 1, -1,
    -1, -1, 1
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_attribute_to_vertices (cmesh);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_attribute_to_vertices (cmesh);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm, int do_dup)
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
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_attribute_to_vertices (cmesh);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

/* The unit cube is constructed from trees of the same eclass.
 * For triangles the square is divided along the (0,0) -- (1,1) diagonal.
 * For prisms the front (y=0) and back (y=1) face are divided into triangles
 * as above.
 */
t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_dup,
                        int do_bcast)
{
  t8_cmesh_t          cmesh;
  int                 num_trees_for_hypercube[T8_ECLASS_LAST] =
    { 1, 1, 1, 2, 1, 6, 2, 3 };
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
    t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
    t8_cmesh_set_num_trees (cmesh, num_trees_for_hypercube[eclass]);
    t8_cmesh_set_attribute_to_vertices (cmesh);
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
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices,
                                  t8_eclass_num_vertices[eclass]);
      break;
    case T8_ECLASS_PRISM:
      t8_cmesh_join_faces (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 2;
      vertices[4] = 3;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 6);
      vertices[1] = 5;
      vertices[2] = 4;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 6);
      break;
    case T8_ECLASS_TRIANGLE:
      t8_cmesh_join_faces (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 3);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 3);
      break;
    case T8_ECLASS_TET:
      t8_cmesh_join_faces (cmesh, 0, 1, 1, 2, 0);
      t8_cmesh_join_faces (cmesh, 1, 2, 1, 2, 0);
      t8_cmesh_join_faces (cmesh, 2, 3, 1, 2, 0);
      t8_cmesh_join_faces (cmesh, 3, 4, 1, 2, 0);
      t8_cmesh_join_faces (cmesh, 4, 5, 1, 2, 0);
      t8_cmesh_join_faces (cmesh, 5, 0, 1, 2, 0);
      vertices[0] = 0;
      vertices[3] = 7;
      vertices[1] = 5;
      vertices[2] = 1;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, attr_vertices, 4);
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 1, attr_vertices, 4);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 2, attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 3, attr_vertices, 4);
      vertices[1] = 6;
      vertices[2] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 4, attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 5;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 5, attr_vertices, 4);
      break;
    case T8_ECLASS_PYRAMID:
      t8_cmesh_join_faces (cmesh, 0, 1, 3, 2, 0);
      t8_cmesh_join_faces (cmesh, 1, 2, 0, 0, 0);
      t8_cmesh_join_faces (cmesh, 2, 0, 2, 0, 1);
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

  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic (sc_MPI_Comm comm, int do_dup, int dim)
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
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_attribute_to_vertices (cmesh);
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
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 1 << dim);
  t8_cmesh_join_faces (cmesh, 0, 0, 0, 1, 0);
  if (dim > 1) {
    t8_cmesh_join_faces (cmesh, 0, 0, 2, 3, 0);
  }
  if (dim == 3) {
    t8_cmesh_join_faces (cmesh, 0, 0, 4, 5, 0);
  }
  t8_cmesh_commit (cmesh);
  return cmesh;
}
