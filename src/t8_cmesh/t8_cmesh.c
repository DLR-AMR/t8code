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

/** \file t8_cmesh.h
 *
 * TODO: document this file
 */

/** This structure hold the connectivity data of the coarse mesh.
 *  It can either be replicated, then each process stores a copy of the whole
 *  mesh, or partitioned. In the latter case, each process only stores a local
 *  portion of the mesh plus information about ghost elements.
 *
 *  The coarse mesh is a collection of coarse trees that can be identified
 *  along faces.
 *  The array ctrees stores these coarse trees sorted by their (global) tree_id.
 *  If the mesh if partitioned it is partitioned according to an (possible only
 *  virtually existing) underlying fine mesh. Therefore the ctrees array can
 *  store duplicated trees on different processes, if each of these processes
 *  owns elements of the same tree in the fine mesh.
 *
 *  Each tree stores information about its face-neighbours in an array of
 *  \ref t8_ctree_fneighbor. \see t8_ctree_fneighbor
 *
 *  If partitioned the ghost trees are stored in a hash table that is backed up
 *  by an array. The hash value of a ghost tree is its tree_id modulo the number
 *  of ghosts on this process.
 */
typedef struct t8_cmesh
{
  /* TODO: make the comments more legible */
  int                 committed;
  int                 dimension; /**< The dimension of the cmesh. It is set when the first tree is inserted. */
  int                 do_dup;   /**< Communicator shall be duped. */
  int                 set_partitioned; /**< If nonzero the cmesh is partitioned.
                                            If zero each process has the whole cmesh. */
  sc_MPI_Comm         mpicomm;  /**< MPI communicator to use. */
  int                 mpirank;  /**< Number of this MPI process. */
  int                 mpisize;  /**< Number of MPI processes. */
  t8_refcount_t       rc; /**< The reference count of the cmesh. */
  t8_topidx_t         num_vertices; /**< The global number of vertices. */
  t8_topidx_t         num_trees;  /**< The global number of trees */
  t8_topidx_t         num_local_trees; /**< If partitioned the number of trees on this process. Otherwise the global number of trees. */
  t8_topidx_t         num_ghosts; /**< If partitioned the number of neighbor trees
                                    owned by different processes. */
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST]; /**< After commit the number of
                                                                 trees for each eclass. */
  sc_array_t         *ctrees; /**< An array of all trees in the cmesh. */
  sc_hash_array_t    *ghosts; /**< The trees that do not belong to this process
                                   but are a face-neighbor of at least one local tree. */
  t8_topidx_t         first_tree; /**< The global index of the first full tree
                                       on this process. Zero if the cmesh is not partitioned. -1 if this processor is empty. */
  t8_topidx_t        *tree_offsets; /**< If partitioned the global number of the
                                         first full tree of each process. */
#ifdef T8_ENABLE_DEBUG
  t8_topidx_t         inserted_trees; /**< Count the number of inserted trees to
                                           check at commit if it equals the total number. */
  t8_topidx_t         inserted_ghosts; /**< Count the number of inserted ghosts to
                                           check at commit if it equals the total number. */
#endif
  /* TODO: make tree_offsets shared array as soon as libsc is updated */
}
t8_cmesh_struct_t;

/** This structure holds the data of a face-neighbor of a tree.
 * The tree_to_face index is computed as follows.
 * Let F be the number of faces of the neighbor tree, then
 * ttf % F is the face number and ttf / F is the orientation.
 * The orientation is determined as follows.  Let my_face and other_face
 * be the two face numbers of the connecting trees.  Then the first
 * face corner of the lower of my_face and other_face connects to a face
 * corner in the higher of my_face and other_face.  The face
 * orientation is defined as the number of this corner.
 * If my_face == other_face, treating
 * either of both faces as the lower one leads to the same result.
 */
/* TODO: This last statement about the same result has to be checked!
 *       It depends on the numbering of the faces as soon as different element
 *       types occur */
typedef struct t8_ctree_fneighbor
{
  t8_topidx_t         treeid; /**< The global number of this neighbor. */
  int                 is_owned; /**< Nonzero if the neighbor belongs to this process. */
  int8_t              tree_to_face;     /* TODO: think of an encoding and document */
}
t8_ctree_fneighbor_struct_t;

typedef struct t8_cghost
{
  t8_topidx_t         treeid; /**< The global number of this ghost. */
  t8_eclass_t         eclass; /**< The eclass of this ghost. */
  int                 owning_proc; /**< The number of the owning process. */
  t8_topidx_t        *local_neighbors; /** Neighbors of this ghost that
                                           are owned by this process. */
}
t8_cghost_struct_t;

typedef struct t8_ctree
{
  t8_topidx_t         treeid; /**< The global number of this tree. */
  t8_eclass_t         eclass; /**< The eclass of this tree. */
  t8_ctree_fneighbor_struct_t *face_neighbors; /**< Information about the face neighbors of this tree. */
}
t8_ctree_struct_t;

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

void
t8_cmesh_set_num_vertices (t8_cmesh_t cmesh, t8_topidx_t num_vertices)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);
  T8_ASSERT (num_vertices > 0);
  T8_ASSERT (cmesh->num_vertices == 0);

  cmesh->num_vertices = num_vertices;
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
static t8_topidx_t
t8_cmesh_tree_index (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  return cmesh->set_partitioned ? tree_id - cmesh->first_tree : tree_id;
}

void
t8_cmesh_set_tree (t8_cmesh_t cmesh, t8_topidx_t tree_id,
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

  tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees,
                                                t8_cmesh_tree_index (cmesh,
                                                                     tree_id));

  tree->eclass = tree_class;
  tree->treeid = tree_id;
  num_neighbors = t8_eclass_num_faces[tree_class];
  /* Allocate neighbors and set entries to invalid values. */
  tree->face_neighbors =
    T8_ALLOC (t8_ctree_fneighbor_struct_t, num_neighbors);
  for (i = 0; i < num_neighbors; i++) {
    tree->face_neighbors[i].treeid = -1;
    tree->face_neighbors[i].tree_to_face = -1;
  }
#ifdef T8_ENABLE_DEBUG
  cmesh->inserted_trees++;
#endif
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
    T1 = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, tree1);
    T2 = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, tree2);
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
    int                 owned_face;
    int                 ghost_face;
    t8_cghost_t         Ghost;
    t8_eclass_t         ghost_eclass;
    size_t              pos;

    /* Find out which one is owned and which one not. */
    if (t8_cmesh_tree_id_is_owned (cmesh, tree1)) {
      owned_id = tree1;
      ghost_id = tree2;
      owned_face = face1;
      ghost_face = face2;
    }
    else {
      T8_ASSERT (t8_cmesh_tree_id_is_owned (cmesh, tree2));
      owned_id = tree2;
      ghost_id = tree1;
      owned_face = face2;
      ghost_face = face1;
    }
    T1 = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, owned_id);

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
t8_cmesh_get_num_vertices (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_vertices;
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

  tree = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees,
                                                t8_cmesh_tree_index (cmesh,
                                                                     tree_id));
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
  /* free trees */
  if (cmesh->ctrees != NULL) {
    for (ti = 0; ti < cmesh->num_local_trees; ti++) {
      treeit = (t8_ctree_t) t8_sc_array_index_topidx (cmesh->ctrees, ti);
      T8_FREE (treeit->face_neighbors);
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

t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
  t8_cmesh_set_tree (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  int                 num_trees_for_hypercube[T8_ECLASS_LAST] =
    { 1, 1, 1, 2, 1, 6, 2, 3 };
  int                 i;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, num_trees_for_hypercube[eclass]);
  for (i = 0; i < num_trees_for_hypercube[eclass]; i++) {
    t8_cmesh_set_tree (cmesh, i, eclass);
  }
  switch (eclass) {
  case T8_ECLASS_PRISM:
  case T8_ECLASS_TRIANGLE:
    t8_cmesh_join_faces (cmesh, 0, 1, 1, 2, 0);
    break;
  case T8_ECLASS_TET:
    t8_cmesh_join_faces (cmesh, 0, 1, 1, 2, 0);
    t8_cmesh_join_faces (cmesh, 1, 2, 1, 2, 0);
    t8_cmesh_join_faces (cmesh, 2, 3, 1, 2, 0);
    t8_cmesh_join_faces (cmesh, 3, 4, 1, 2, 0);
    t8_cmesh_join_faces (cmesh, 4, 5, 1, 2, 0);
    t8_cmesh_join_faces (cmesh, 5, 0, 1, 2, 0);
    break;
  case T8_ECLASS_PYRAMID:
    t8_cmesh_join_faces (cmesh, 0, 1, 3, 2, 0);
    t8_cmesh_join_faces (cmesh, 1, 2, 0, 0, 0);
    t8_cmesh_join_faces (cmesh, 2, 0, 2, 0, 1);
    break;
  default:
    break;
  }

  t8_cmesh_commit (cmesh);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic (sc_MPI_Comm comm, int do_dup, int dim)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         tree_class;

  T8_ASSERT (dim == 1 || dim == 2 || dim == 3);
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  t8_cmesh_set_num_trees (cmesh, 1);
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

  t8_cmesh_set_tree (cmesh, 0, tree_class);
  /* TODO: if orientation is specified, check whether 0 is the correct choice here */
  t8_cmesh_join_faces (cmesh, 0, 0, 0, 1, 0);
  t8_cmesh_join_faces (cmesh, 0, 0, 1, 0, 0);
  if (dim > 1) {
    t8_cmesh_join_faces (cmesh, 0, 0, 2, 3, 0);
    t8_cmesh_join_faces (cmesh, 0, 0, 3, 2, 0);
  }
  if (dim == 3) {
    t8_cmesh_join_faces (cmesh, 0, 0, 4, 5, 0);
    t8_cmesh_join_faces (cmesh, 0, 0, 5, 4, 0);
  }
  t8_cmesh_commit (cmesh);
  return cmesh;
}
