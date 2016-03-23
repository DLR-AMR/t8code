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
#include "t8_cmesh_stash.h"

/** \file t8_cmesh_types.h
 * We define here the datatypes needed for internal cmesh routines.
 */

typedef struct t8_stash *t8_stash_t;
typedef struct t8_part_tree *t8_part_tree_t;
typedef struct t8_cmesh_trees *t8_cmesh_trees_t;

/** This structure holds the connectivity data of the coarse mesh.
 *  It can either be replicated, then each process stores a copy of the whole
 *  mesh, or partitioned. In the latter case, each process only stores a local
 *  portion of the mesh plus information about ghost elements.
 *
 *  The coarse mesh is a collection of coarse trees that can be identified
 *  along faces.
 *  TODO: this description is outdated. rewrite it.
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
  int                 face_knowledge;  /**< If partitioned the level of face knowledge that is expected. \ref t8_mesh_set_partioned */
  int8_t              set_level;       /**< Non-negative if the cmesh should be partition from an already existing cmesh
                                         with an assumes \a level uniform mesh underneath */
  struct t8_cmesh    *set_from; /**< If this cmesh is a modified (i.e. partitioned) version
                                        of another cmesh, we store a pointer to this other cmesh here. */
  sc_MPI_Comm         mpicomm;  /**< MPI communicator to use. */
  int                 mpirank;  /**< Number of this MPI process. */
  int                 mpisize;  /**< Number of MPI processes. */
  t8_refcount_t       rc; /**< The reference count of the cmesh. */
  t8_gloidx_t         num_trees;   /**< The global number of trees */
  t8_locidx_t         num_local_trees; /**< If partitioned the number of trees on this process. Otherwise the global number of trees. */
  t8_locidx_t         num_ghosts; /**< If partitioned the number of neighbor trees
                                    owned by different processes. */
  /* TODO: wouldnt a local num_trees_per_eclass be better? */
  t8_gloidx_t         num_trees_per_eclass[T8_ECLASS_LAST]; /**< After commit the number of
                                                                 trees for each eclass. */

  t8_cmesh_trees_t    trees; /**< structure that holds all local trees and ghosts */

  t8_gloidx_t         first_tree; /**< The global index of the first local tree
                                       on this process. Zero if the cmesh is not partitioned. -1 if this processor is empty. */
  int8_t             first_tree_shared; /**< If partitioned true if the first tree on this process is also the last tree on the next process.
                                             Always zero if num_local_trees = 0 */
  /* TODO: deprecated, replaced by offset */
  t8_gloidx_t        *tree_offsets;  /**< If partitioned for each process the global index of its first local tree
                                        or -(first local tree) - 1
                                        if the first tree on that process is shared.
                                        Since this is very memory consuming we only fill it when needed. */
#ifdef T8_ENABLE_DEBUG
  t8_topidx_t         inserted_trees; /**< Count the number of inserted trees to
                                           check at commit if it equals the total number. */
  t8_topidx_t         inserted_ghosts; /**< Count the number of inserted ghosts to
                                           check at commit if it equals the total number. */
#endif
  t8_stash_t          stash; /**< Used as temporary storage for the trees before commit. */
  /* TODO: make tree_offsets shared array as soon as libsc is updated */
}
t8_cmesh_struct_t;

typedef struct t8_cghost
{
  t8_gloidx_t         treeid; /**< The global number of this ghost. */
  t8_eclass_t         eclass; /**< The eclass of this ghost. */
  size_t              neigh_offset; /* TODO: document */
}
t8_cghost_struct_t;

/** This structure holds the data of a local tree including the information
 * about face neighbors. For those
 * the tree_to_face index is computed as follows.
 * Let F be the maximal number of faces of any eclass of the cmesh's dimension, then
 * ttf % F is the face number and ttf / F is the orientation. (\ref t8_eclass_max_num_faces)
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
  size_t              neigh_offset;  /**< Adding this offset to the adress of the tree
                                       yield the array of face_neighbor entries */
  size_t              att_offset;    /**< Adding this offset to the adress of the tree
                                       yield the array of attribute_info entries */
  int                 num_attributes; /**< The number of attributes at this tree */
}
t8_ctree_struct_t;

/** This structure hold the information associated to an attribute of a tree.
 *  The attributes of each are stored in a key-value storage, where the key consists
 *  of the two entries (package_id,key) both being integers.
 *  The package_id serves to identify the application layer that added the attribute
 *  and the key identifies the attribute whithin that application layer.
 *
 *  All attribute info objects of one tree are stored in an array and adding
 *  a tree's att_offset entry to the tree's adress yields this array.
 *  The attributes themselfes are stored in an array directly behind the array of
 *  the attribute infos.
 */
typedef struct t8_attribute_info
{
  int       package_id; /**< The identifier of the application layer that added this attribute */
  int       key; /**< The (tree unique) key of the attribute whithin this AL. */
  size_t    attribute_offset; /**< The offset of the attribute data from the first
                    attribute info of the tree.
                    (Thus, the attribute is stored at adress tree + tree->att_offset + attribute_offset) */
  /* TODO: eventually remove the size */
  size_t    attribute_size; /**< The size in bytes of the attribute */
} t8_attribute_info_struct_t;

/* TODO: document, process is a bad naming, since it does not refer to MPI ranks here */
typedef struct t8_cmesh_trees
{
  sc_array_t         *from_proc;        /* array of t8_part_tree, one for each process */
  int                *tree_to_proc;     /* for each tree its process */
  int                *ghost_to_proc;    /* for each ghost its process */
} t8_cmesh_trees_struct_t;

/* TODO: document */
typedef struct t8_part_tree
{
  char               *first_tree;       /* Stores the trees, the ghosts and the attributes.
                                           The last 2*sizeof(t8_topidx) bytes store num_trees and num_ghosts */
  t8_locidx_t         first_tree_id;    /* local tree_id of the first tree. -1 if num_trees = 0 */
  t8_locidx_t         first_ghost_id;   /* TODO: document. -1 if num_ghost=0 */
  t8_locidx_t         num_trees;
  t8_locidx_t         num_ghosts;
}
t8_part_tree_struct_t;

#endif /* !T8_CMESH_TYPES_H */
