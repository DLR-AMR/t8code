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

#ifndef T8_FOREST_TYPES_H
#define T8_FOREST_TYPES_H

/** \file t8_forest_types.h
 * We define here the datatypes needed for internal forest routines.
 */

#include <t8.h>
#include <t8_refcount.h>
#include <t8_cmesh.h>
#include <t8_element.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>

typedef struct t8_profile t8_profile_t; /* Defined below */
typedef struct t8_forest_ghost *t8_forest_ghost_t;      /* Defined below */

typedef enum t8_forest_from
{
  T8_FOREST_FROM_FIRST,
  T8_FOREST_FROM_COPY = T8_FOREST_FROM_FIRST,
  T8_FOREST_FROM_ADAPT,
  T8_FOREST_FROM_PARTITION,
  T8_FOREST_FROM_LAST
}
t8_forest_from_t;

/** This structure is private to the implementation. */
typedef struct t8_forest
{
  t8_refcount_t       rc;               /**< Reference counter. */

  int                 set_level;        /**< Level to use in new construction. */
  int                 set_for_coarsening;       /**< Change partition to allow
                                                     for one round of coarsening */

  sc_MPI_Comm         mpicomm;          /**< MPI communicator to use. */
  t8_cmesh_t          cmesh;            /**< Coarse mesh to use. */
  //t8_scheme_t        *scheme;        /**< Scheme for element types. */
  t8_scheme_cxx_t    *scheme_cxx;        /**< Scheme for element types. */
  int                 do_dup;           /**< Communicator shall be duped. */
  int                 dimension;        /**< Dimension inferred from \b cmesh. */

  t8_forest_t         set_from;         /**< Temporarily store source forest. */
  t8_forest_from_t    from_method;      /**< Method to derive from \b set_from. */
  t8_forest_replace_t set_replace_fn;   /**< Replace function. Called when \b from_method
                                             is set to T8_FOREST_FROM_ADAPT. */
  t8_forest_adapt_t   set_adapt_fn;     /**< refinement and coarsen function. Called when \b from_method
                                             is set to T8_FOREST_FROM_ADAPT. */
  int                 set_adapt_recursive; /**< Flag to decide whether coarsen and refine
                                                are carried out recursive */
  int                 do_ghost;         /**< If True, a ghost layer will be created when the forest is committed. */
  t8_ghost_type_t     ghost_type;       /**< If a ghost layer will be created, the type of neighbors that count as ghost. */
  void               *user_data;        /**< Pointer for arbitrary user data. \see t8_forest_set_user_data. */
  int                 committed;        /**< \ref t8_forest_commit called? */
  int                 mpisize;          /**< Number of MPI processes. */
  int                 mpirank;          /**< Number of this MPI process. */

  t8_gloidx_t         first_local_tree;
  t8_gloidx_t         last_local_tree;
  t8_gloidx_t         global_num_trees; /**< The total number of global trees */
  sc_array_t         *trees;
  t8_forest_ghost_t   ghosts;           /**< The ghost elements. \see t8_forest_ghost.h */
  t8_shmem_array_t    element_offsets; /**< If partitioned, for each process the global index
                                            of its first element. Since it is memory consuming,
                                            it is usually only constructed when needed and otherwise unallocated. */
  t8_shmem_array_t    global_first_desc; /**< If partitioned, for each process the linear id (at maxlevel) of its
                                              first element's first descendant.
                                             \ref t8_element_set_linear_id. Stores 0 for empty processes.
                                            Since it is memory consuming,
                                            it is usually only constructed when needed and otherwise unallocated. */
  t8_shmem_array_t    tree_offsets; /**<  If partitioned for each process the global index of its first local tree
                                          or -(first local tree) - 1
                                          if the first tree on that process is shared.
                                          Since this is memory consuming we only construct it when needed.
                                          This array follows the same logic as \a tree_offsets in \a t8_cmesh_t */

  t8_locidx_t         local_num_elements;  /**< Number of elements on this processor. */
  t8_gloidx_t         global_num_elements; /**< Number of elements on all processors. */
  t8_profile_t       *profile; /**< If not NULL, runtimes and statistics about forest_commit are stored here. */

}
t8_forest_struct_t;

/** The t8 tree datatype */
typedef struct t8_tree
{
  sc_array_t          elements;              /**< locally stored elements */
  t8_eclass_t         eclass;                /**< The element class of this tree */
  /* TODO: We will need the *_desc variables later for shure. */
  t8_element_t       *first_desc,            /**< first local descendant */
                     *last_desc;             /**< last local descendant */
  t8_locidx_t         elements_offset;      /**< cumulative sum over earlier
                                                  trees on this processor
                                                  (locals only) */
}
t8_tree_struct_t;

/** This struct is used to profile forest algorithms.
 * The forest struct stores a pointer to a profile struct, and if
 * it is nonzero, various runtimes and data measurements are stored here.
 * \see t8_cmesh_set_profiling and \see t8_cmesh_print_profile
 */

/** The number of statistics collected by a profile struct. */
#define T8_PROFILE_NUM_STATS 6
typedef struct t8_profile
{
  t8_locidx_t         partition_elements_shipped; /**< The number of elements this process has
                                                  sent to other in the last partition call. */
  t8_locidx_t         partition_elements_recv; /**< The number of elements this process has
                                                  received from other in the last partition call. */
  size_t              partition_bytes_sent; /**< The total number of bytes sent to other processes in the
                                                 last partition call. */
  int                 partition_procs_sent;  /**< The number of different processes this process has send
                                            local elements to in the last partition call. */
  double              partition_runtime; /**< The runtime of  the last call to \a t8_cmesh_partition. */
  double              commit_runtime; /**< The runtim of the last call to \a t8_cmesh_commit. */

}
t8_profile_struct_t;

/* TODO: document */
typedef struct t8_forest_ghost
{
  t8_refcount_t       rc; /**< The reference counter. */

  t8_ghost_type_t     ghost_type;   /**< Describes which neighbors are considered ghosts. */
  sc_array_t         *ghost_trees;      /* ghost tree data:
                                           global_id.
                                           eclass.
                                           elements. In linear id order */
  sc_hash_t          *global_tree_to_ghost_tree;        /* Indexes into ghost_trees.
                                                           Given a global tree id I give the index
                                                           i such that the tree is in ghost_trees[i]
                                                         */
  sc_hash_t          *process_offsets;  /* Given a process, return the first ghost tree and
                                           whithin it the first element of that process. */
#if 0
  /* TODO: obsolete by remote_processes below. */
  sc_array_t         *processes;        /* ranks of the processes */
#endif
  sc_hash_array_t    *remote_ghosts;    /* array of local trees that have ghost elements for another process.
                                           for each tree an array of t8_element_t * pointing to the local ghost elements.
                                           It is a hash table, hashed with the rank of a remote process.
                                           Sorted within each process by linear id.
                                         */
  sc_array_t         *remote_processes; /* The ranks of the processes for which local elements are ghost.
                                           Array of int's. */

  sc_mempool_t       *glo_tree_mempool;
  sc_mempool_t       *proc_offset_mempool;
} t8_forest_ghost_struct_t;

#endif /* ! T8_FOREST_TYPES_H! */
