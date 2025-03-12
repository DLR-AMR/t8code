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
#include <t8_data/t8_shmem.h>
#include <t8_geometry/t8_geometry.h>
#include "t8_cmesh_stash.h"
#include "t8_element.h"

/** \file t8_cmesh_types.h
 * We define here the datatypes needed for internal cmesh routines.
 */

typedef struct t8_part_tree *t8_part_tree_t;
typedef struct t8_cmesh_trees *t8_cmesh_trees_t;
typedef struct t8_cprofile t8_cprofile_t; /* Defined below */

/* TODO: no longer needed.
 *       User may use set_derived_from, then set_from is non-NULL.
 *       When refinement is desired, level shall be set positive.
 *       When partition is desired, we expect set_partition to be true.
 *       Any combination can be called.
 *       In t8_commit we check whether the parameters are consistent,
 *       and whether a combination is currently supported,
 *       and provide informative error messages both in debug and non-debug.
 */

/* Definitions for attribute identifiers that are reserved for a special purpose. 
 * T8_CMESH_NEXT_POSSIBLE_KEY is the first unused key, hence it can be repurposed for different attributes.*/
#define T8_CMESH_VERTICES_ATTRIBUTE_KEY 0            /* Used to store vertex coordinates. */
#define T8_CMESH_GEOMETRY_ATTRIBUTE_KEY 1            /* Used to store the name of a tree's geometry. */
#define T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY 2            /* Used to store which edge is linked to which geometry */
#define T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY 3 /* Used to store edge parameters */
#define T8_CMESH_CAD_FACE_ATTRIBUTE_KEY \
  T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY \
  +T8_ECLASS_MAX_EDGES /* Used to store which face is linked to which surface */
#define T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY \
  T8_CMESH_CAD_FACE_ATTRIBUTE_KEY + 1 /* Used to store face parameters */
#define T8_CMESH_LAGRANGE_POLY_DEGREE_KEY T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + T8_ECLASS_MAX_FACES
#define T8_CMESH_NEXT_POSSIBLE_KEY \
  T8_CMESH_LAGRANGE_POLY_DEGREE_KEY + 1 /* The next free value for a t8code attribute key */

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
  int committed; /**< Flag that specifies whether the cmesh is committed or not. \ref t8_cmesh_commit */
  int dimension; /**< The dimension of the cmesh. It is set when the first tree is inserted. */

  int set_partition;  /**< If nonzero the cmesh is partitioned.
                                            If zero each process has the whole cmesh. */
  int face_knowledge; /**< If partitioned the level of face knowledge that is expected. \ref t8_cmesh_set_partitioned;
                            see \ref t8_cmesh_set_partition.
*/

  const t8_scheme_c *set_partition_scheme; /**< If the cmesh is to be partitioned according to a uniform level,
                                                the scheme that describes the refinement pattern. See \ref t8_cmesh_set_partition. */
  int8_t set_partition_level;  /**< Non-negative if the cmesh should be partitioned from an already existing cmesh
                                         with an assumed \a level uniform mesh underneath. */
  struct t8_cmesh *set_from;   /**< If this cmesh shall be derived from an
                                  existing cmesh by copy or more elaborate
                                  modification, we store a pointer to this
                                  other cmesh here. */
  int mpirank;                 /**< Number of this MPI process. */
  int mpisize;                 /**< Number of MPI processes. */
  t8_refcount_t rc;            /**< The reference count of the cmesh. */
  t8_gloidx_t num_trees;       /**< The global number of trees */
  t8_locidx_t num_local_trees; /**< If partitioned the number of trees on this process. 
                                    Otherwise the global number of trees. */
  t8_locidx_t num_ghosts;      /**< If partitioned the number of neighbor trees
                                    owned by different processes. */
  /* TODO: wouldnt a local num_trees_per_eclass be better?
   *       only as an additional info. we need the global count. i.e. for forest_maxlevel computation.
   */
  t8_locidx_t num_local_trees_per_eclass[T8_ECLASS_COUNT]; /**< After commit the number of local
                                                                 trees for each eclass.
                                                                 Stores the same entries as \a num_trees_per_eclass,
                                                                 if the cmesh is replicated. */
  t8_gloidx_t num_trees_per_eclass[T8_ECLASS_COUNT];       /**< After commit the number of global
                                                                 trees for each eclass. */

  t8_cmesh_trees_t trees; /**< structure that holds all local trees and ghosts */

  t8_gloidx_t first_tree;   /**< The global index of the first local tree on this process. 
                                       Zero if the cmesh is not partitioned. -1 if this processor is empty.
                                       See also https://github.com/DLR-AMR/t8code/wiki/Tree-indexing */
  int8_t first_tree_shared; /**< If partitioned true if the first tree on this process is also the last tree 
                                  on the next process. Always zero if num_local_trees = 0 */

  t8_shmem_array_t tree_offsets; /**< If partitioned for each process the global index of its first local tree
                                        or -(first local tree) - 1
                                        if the first tree on that process is shared.
                                        Since this is very memory consuming we only fill it when needed. */

  t8_geometry_handler_c *geometry_handler; /**< Handles all geometries that are used by trees in this cmesh. */
  int negative_volume_check;               /**< Check cmesh for negative volumes during commit. */

#ifdef T8_ENABLE_DEBUG
  t8_locidx_t inserted_trees;  /**< Count the number of inserted trees to
                                           check at commit if it equals the total number. */
  t8_locidx_t inserted_ghosts; /**< Count the number of inserted ghosts to
                                           check at commit if it equals the total number. */
#endif
  t8_stash_t stash;       /**< Used as temporary storage for the trees before commit. */
  t8_cprofile_t *profile; /**< Used to measure runtimes and statistics of the cmesh algorithms. */
} t8_cmesh_struct_t;

/* TODO: cghost could be the same type as ctree.
 *       treeid for ghosts then negative?!
 *       1. typedef cghost ctree
 *       2. completely replace
 */
/* TODO: document */
typedef struct t8_cghost
{
  t8_gloidx_t treeid;  /**< The global number of this ghost. */
  t8_eclass_t eclass;  /**< The eclass of this ghost. */
  size_t neigh_offset; /** Offset to the array of face neighbors of this ghost.
                                        This count has to be added to the address of the ghost to get its face neighbors. */
  size_t att_offset;   /**< Adding this offset to the address of the ghost
                                       yields the array of attribute_info entries */
  /* TODO: Could be a size_t */
  int num_attributes; /**< The number of attributes at this ghost */
} t8_cghost_struct_t;

/** This structure holds the data of a local tree including the information
 * about face neighbors. For those
 * the tree_to_face index is computed as follows.
 * Let F be the maximal number of faces of any eclass of the cmesh's dimension, then
 * ttf % F is the face number and ttf / F is the orientation. (\ref t8_eclass_max_num_faces)
 * The orientation is determined as follows.  Let my_face and other_face
 * be the two face numbers of the connecting trees.
 * We chose a main_face from them as follows: Either both trees have the same
 * element class, then the face with the lower face number is the main_face or
 * the trees belong to different classes in which case the face belonging to the
 * tree with the lower class according to the ordering
 * triangle < square,
 * hex < tet < prism < pyramid,
 * is the main_face.
 * Then face corner 0 of the main_face connects to a face
 * corner k in the other face.  The face orientation is defined as the number k.
 * If the classes are equal and my_face == other_face, treating
 * either of both faces as the main_face leads to the same result.
 * See https://arxiv.org/pdf/1611.02929.pdf for more details.
 */
typedef struct t8_ctree
{
  t8_locidx_t treeid; /**< The local number of this tree. */
  /* TODO: The local id of a tree should be clear from context, the entry can
   *       be optimized out. */
  t8_eclass_t eclass;  /**< The eclass of this tree. */
  size_t neigh_offset; /**< Adding this offset to the address of the tree
                                       yields the array of face_neighbor entries */
  size_t att_offset;   /**< Adding this offset to the address of the tree
                                       yields the array of attribute_info entries */
  /* TODO: Could be a size_t */
  int num_attributes; /**< The number of attributes at this tree */
} t8_ctree_struct_t;

/** This structure holds the information associated to an attribute of a tree.
 *  The attributes of each are stored in a key-value storage, where the key consists
 *  of the two entries (package_id,key) both being integers.
 *  The package_id serves to identify the application layer that added the attribute
 *  and the key identifies the attribute within that application layer.
 *
 *  All attribute info objects of one tree are stored in an array and adding
 *  a tree's att_offset entry to the tree's address yields this array.
 *  The attributes themselves are stored in an array directly behind the array of
 *  the attribute infos.
 */
typedef struct t8_attribute_info
{
  int package_id;
  /**< The identifier of the application layer that added this attribute */
  int key;
  /**< The (tree unique) key of the attribute within this AL. */
  size_t attribute_offset;
  /**< The offset of the attribute data from the first
                    attribute info of the tree.
                    (Thus, the attribute is stored at address tree + tree->att_offset + attribute_offset) */
  /* TODO: eventually remove the size */
  size_t attribute_size;
  /**< The size in bytes of the attribute */
} t8_attribute_info_struct_t;

/* TODO: document, process is a bad naming, since it does not refer to MPI ranks here */
typedef struct t8_cmesh_trees
{
  sc_array_t *from_proc;                 /* array of t8_part_tree, one for each process */
  int *tree_to_proc;                     /* for each tree its process */
  int *ghost_to_proc;                    /* for each ghost its process */
  sc_hash_t *ghost_globalid_to_local_id; /* A hash table storing the map
                                                           global_id -> local_id for the ghost trees.
                                                           The local_id is the local ghost id starting at num_local_trees  */
  sc_mempool_t *global_local_mempool;    /* Memory pool for the entries in the hash table */
} t8_cmesh_trees_struct_t;

/* TODO: document */
typedef struct t8_part_tree
{
  char *first_tree;           /* Stores the trees, the ghosts and the attributes.
                                           The last 2*sizeof(t8_locidx) bytes store num_trees and num_ghosts */
  t8_locidx_t first_tree_id;  /* local tree_id of the first tree. -1 if num_trees = 0 */
  t8_locidx_t first_ghost_id; /* TODO: document. -1 if num_ghost=0, 0 for the first part, (not num_local_trees!)
                                           0 <= first_ghost_id < num_ghosts */
  t8_locidx_t num_trees;
  t8_locidx_t num_ghosts;
} t8_part_tree_struct_t;

/* TODO: Extend this structure with meaningful entries.
 *       Maybe the number of shipped trees per process is useful?
 */
/** This struct is used to profile cmesh algorithms.
 * The cmesh struct stores a pointer to a profile struct, and if
 * it is nonzero, various runtimes and data measurements are stored here.
 * \see t8_cmesh_set_profiling and \see t8_cmesh_print_profile
 */
typedef struct t8_cprofile
{
  t8_locidx_t partition_trees_shipped;  /**< The number of trees this process has
                                                 sent to other in the last partition call. */
  t8_locidx_t partition_ghosts_shipped; /**< The number of ghosts this process has
                                                 sent to other in the last partition call. */
  t8_locidx_t partition_trees_recv;     /**< The number of trees this process has
                                                 received from other in the last partition call. */
  t8_locidx_t partition_ghosts_recv;    /**< The number of ghosts this process has
                                                 received from other in the last partition call. */
  size_t partition_bytes_sent;          /**< The total number of bytes sent to other processes in the
                                                last partition call. */
  int partition_procs_sent;             /**< The number of different processes this process has send
                                           local trees or ghosts to in the last partition call. */
  int first_tree_shared;                /**< 1 if this processes' first tree is shared. 0 if not. */
  double partition_runtime;             /**< The runtime of  the last call to \a t8_cmesh_partition. */
  double commit_runtime;                /**< The runtime of the last call to \a t8_cmesh_commit. */
  double geometry_evaluate_num_calls;   /**< The number of calls to \a t8_geometry_evaluate. */
  double geometry_evaluate_runtime;     /**< The accumulated runtime of calls to \a t8_geometry_evaluate. */
} t8_cprofile_struct_t;

/** The number of entries in a cprofile struct */
#define T8_CPROFILE_NUM_STATS 11

#endif /* !T8_CMESH_TYPES_H */
