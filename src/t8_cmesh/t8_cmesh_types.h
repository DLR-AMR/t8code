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
#include <t8_cmesh/t8_cmesh_part_tree.h>

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
  t8_gloidx_t         num_trees;   /**< The global number of trees */
  t8_topidx_t         num_local_trees; /**< If partitioned the number of trees on this process. Otherwise the global number of trees. */
  t8_topidx_t         num_ghosts; /**< If partitioned the number of neighbor trees
                                    owned by different processes. */
  t8_gloidx_t         num_trees_per_eclass[T8_ECLASS_LAST]; /**< After commit the number of
                                                                 trees for each eclass. */

  t8_part_tree_t     *trees_ghosts;
  int                 num_parts; /** Number of entries in \a trees_ghosts */
  size_t             *attribute_size; /** If attributes are used, for each tree the size of its attribute */
  size_t             *attribute_offset; /* TODO: document, for each tree offset into the part array */
#if 0
  sc_array_t         *ctrees; /**< An array of all trees in the cmesh. */
  sc_array_t         *ghosts; /**< The trees that do not belong to this process
                                   but are a face-neighbor of at least one local tree. */
#endif
  t8_gloidx_t         first_tree; /**< The global index of the first local tree
                                       on this process. Zero if the cmesh is not partitioned. -1 if this processor is empty. */
  t8_topidx_t        *tree_per_proc; /**< If partitioned twice the number of local
                                          trees on each process plus one if the last tree of the respective
                                          process is the first tree of the next process */
  sc_mempool_t       *tree_attributes_mem[T8_ECLASS_LAST]; /**< For each eclass we can specify an
                                         attribute size and attach attributes of this size to each trees */
#ifdef T8_ENABLE_DEBUG
  t8_topidx_t         inserted_trees; /**< Count the number of inserted trees to
                                           check at commit if it equals the total number. */
  t8_topidx_t         inserted_ghosts; /**< Count the number of inserted ghosts to
                                           check at commit if it equals the total number. */
#endif
  /* TODO: make tree_offsets shared array as soon as libsc is updated */
}
t8_cmesh_struct_t;

typedef struct t8_cghost
{
  t8_gloidx_t         treeid; /**< The global number of this ghost. */
  t8_eclass_t         eclass; /**< The eclass of this ghost. */
  t8_gloidx_t        *neighbors; /**< Global id's of all neighbors of this ghost */
}
t8_cghost_struct_t;

/** This structure holds the data of a local tree including the information
 * about face neighbors. For those
 * the tree_to_face index is computed as follows.
 * Let F be the number of faces of the neighbor tree, then
 * ttf % F is the face number and ttf / F is the orientation.
 * The orientation is determined as follows.  Let my_face and other_face
 * be the two face numbers of the connecting trees.
 * We chose a master_face from them as follows: Either both trees have the same
 * element class, then the face with the lower face number is the master_face or
 * the trees belong to different classes in which case the face belonging to the
 * tree with the lower class according to the ordering
 * triangle < square,
 * hex < tet < prism < pyramid,
 * is the master_face.
 * Then the first face corner of the master_face connects to a face
 * corner in the other face.  The face
 * orientation is defined as the number of this corner.
 * If the classes are equal and my_face == other_face, treating
 * either of both faces as the master_face leads to the same result.
 */
typedef struct t8_ctree
{
  t8_topidx_t         treeid; /**< The local number of this tree. */
  /* TODO: The local id of a tree should be clear from context, the entry can
   *       be optimized out. */
  t8_eclass_t         eclass; /**< The eclass of this tree. */
  t8_topidx_t        *face_neighbors; /**< For each face the local index of the face neighbor
                                          of this tree at the face. Indices greater than
                                          the number of local trees refer to ghosts. */
  int8_t             *tree_to_face; /**< For each face the encoding of the face neighbor orientation. */
}
t8_ctree_struct_t;

/* TODO: document */
typedef struct t8_part_tree
{
  char              *first_tree;
  t8_topidx_t        num_trees;
  t8_topidx_t        num_ghosts;
#if 0
  /* TODO: Do we need this? */
  size_t             num_bytes_for_attributes;
#endif
}
t8_part_tree_struct_t;

#endif /* !T8_CMESH_TYPES_H */
