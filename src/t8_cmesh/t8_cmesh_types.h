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

#ifndef T8_CMESH_TYPES_H
#define T8_CMESH_TYPES_H

#include <t8.h>
#include <t8_refcount.h>

/** \file t8_cmesh_types.h
 * We define here the datatypes needed for internal cmesh routines.
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
  t8_topidx_t         num_corners; /**< The global number of corners that help define the topology. Is allowed to be zero if topology and geometry are equal. */
  t8_topidx_t         num_local_corners; /**< If partitioned the local number of corners. Otherwise the global number of corners. */
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
  size_t              tree_attribute_size[T8_ECLASS_LAST];
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
  /* TODO: write a macro instead of is_owned */
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
  t8_topidx_t        *corners; /**< The corner indices of this tree's corners. Can be NULL if \a cmesh.num_corners is 0. */
  t8_topidx_t        *vertices; /**< The vertex indices of this tree's corners. This defines an embedding of the tree into \f$R^3$\f. */
  t8_ctree_fneighbor_struct_t *face_neighbors; /**< Information about the face neighbors of this tree. */
  void               *attribute;
}
t8_ctree_struct_t;

#endif /* !T8_CMESH_TYPES_H */
