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

/** \file t8_forest_ghost.cxx
 * Implements functions declared in \ref t8_forest_ghost.h.
 */

#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_data/t8_containers.h>
#include <sc_statistics.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/**
 * The information for a remote process, what data we have to send to them.
 */
using t8_ghost_mpi_send_info_t = struct
{
  int recv_rank;           /**< The rank to which we send. */
  size_t num_bytes;        /**< The number of bytes that we send. */
  sc_MPI_Request *request; /**< Communication request, not owned by this struct. */
  char *buffer;            /**< The send buffer. */
};

/** 
 * The information stored for the ghost trees 
 */
using t8_ghost_tree_t = struct
{
  t8_gloidx_t global_id;       /**< Global id of the tree */
  t8_locidx_t element_offset;  /**< The count of all ghost elements in all smaller ghost trees */
  t8_element_array_t elements; /**< The ghost elements of that tree */
  t8_eclass_t eclass;          /**< The trees element class */
};

/**
 * The data structure stored in the global_tree_to_ghost_tree hash table. 
 */
using t8_ghost_gtree_hash_t = struct
{
  t8_gloidx_t global_id; /**< global tree id */
  size_t index;          /**< the index of that global tree in the ghost_trees array. */
};

/**
 * The data structure stored in the process_offsets array. 
 */
using t8_ghost_process_hash_t = struct
{
  int mpirank;              /**< rank of the process */
  t8_locidx_t ghost_offset; /**< The number of ghost elements for all previous ranks */
  size_t tree_index;        /**< index of first ghost tree of this process in ghost_trees */
  size_t first_element;     /**< the index of the first element in the elements array of the ghost tree. */
};

/**
 * The information stored for the remote trees.
 * Each remote process stores an array of these 
 */
using t8_ghost_remote_tree_t = struct
{
  t8_gloidx_t global_id;       /**< global id of the tree */
  int mpirank;                 /**< The mpirank of the remote process */
  t8_element_array_t elements; /**< The remote elements of that tree */
  sc_array_t element_indices;  /**< The (tree) local indices of the ghost elements. */
  t8_eclass_t eclass;          /**< The trees element class */
};

/**
 * This struct stores information about the data that the current process needs from a specific remote_process
 * as ghost data, such as the number of remote elements and the remote trees.
*/
using t8_ghost_remote_t = struct
{
  int remote_rank;          /**< The rank of the remote process */
  t8_locidx_t num_elements; /**< The number of remote elements for this process */
  sc_array_t remote_trees;  /**< Array of the remote trees of this process */
};

/** 
 * The hash function for the global tree hash.
 * As hash value we just return the global tree id. 
 * 
 * \param[in] ghost_gtree_hash  Global tree hash.
 * \param[in] data              Unused dummy Argument to allow passing this function to sc_hash_new. 
 * 
 * \return The global tree id.
 */
static unsigned
t8_ghost_gtree_hash_function (const void *ghost_gtree_hash, [[maybe_unused]] const void *data)
{
  const t8_ghost_gtree_hash_t *object = (const t8_ghost_gtree_hash_t *) ghost_gtree_hash;

  return (unsigned) object->global_id;
}

/**
 * The equal function for two global tree hash objects.
 * Two t8_ghost_gtree_hash_t are considered equal if their global tree ids are the same.
 * 
 * \param[in] ghost_gtreea  Global tree hash object A
 * \param[in] ghost_gtreeb  Global tree hash object B
 * \param[in] user          Unused dummy Argument to allow passing this function to sc_hash_new. 
 * 
 * \return 1 if equal, 0 otherwise.
 */
static int
t8_ghost_gtree_equal_function (const void *ghost_gtreea, const void *ghost_gtreeb, [[maybe_unused]] const void *user)
{
  const t8_ghost_gtree_hash_t *objecta = (const t8_ghost_gtree_hash_t *) ghost_gtreea;
  const t8_ghost_gtree_hash_t *objectb = (const t8_ghost_gtree_hash_t *) ghost_gtreeb;

  /* return true if and only if the global_ids are the same */
  return objecta->global_id == objectb->global_id;
}

/** 
 * The hash value for an entry of the process_offsets hash is the processes mpirank. 
 * 
 * \param[in] process_data    Process data as void pointer.
 * \param[in] user_data       Unused dummy Argument to allow passing this function to sc_hash_new. 
 * 
 * \return The process' MPI rank.
 * 
 */
static unsigned
t8_ghost_process_hash_function (const void *process_data, [[maybe_unused]] const void *user_data)
{
  const t8_ghost_process_hash_t *process = (const t8_ghost_process_hash_t *) process_data;

  return process->mpirank;
}

/** 
 * The equal function for the process_offsets array.
 * Two entries are the same if their mpiranks are equal. 
 * 
 * \param[in] process_dataa   Process offset array A
 * \param[in] process_datab   Process offset array B
 * \param[in] user            Unused dummy Argument to allow passing this function to sc_hash_new. 
 * 
 * \return 1 if equal, 0 if not.
 */
static int
t8_ghost_process_equal_function (const void *process_dataa, const void *process_datab,
                                 [[maybe_unused]] const void *user)
{
  const t8_ghost_process_hash_t *processa = (const t8_ghost_process_hash_t *) process_dataa;
  const t8_ghost_process_hash_t *processb = (const t8_ghost_process_hash_t *) process_datab;

  return processa->mpirank == processb->mpirank;
}

/** 
 * The hash function for the remote_ghosts hash table.
 * The hash value for an mpirank is just the rank.
 * 
 * \param[in] remote_data The remote ghost data as void pointer.
 * \param[in] user_data   Unused dummy Argument to allow passing this function to sc_hash_new.
 * 
 * \return The mpi rank.
 */
static unsigned
t8_ghost_remote_hash_function (const void *remote_data, [[maybe_unused]] const void *user_data)
{
  const t8_ghost_remote_t *remote = (const t8_ghost_remote_t *) remote_data;

  return remote->remote_rank;
}

/**
 * The equal function for the remote hash table.
 * Two entries are the same if they have the same rank. 
 * 
 * \param[in] remote_dataa  Remote hash table A.
 * \param[in] remote_datab  Remote hash table B.
 * \param[in] user          Unused dummy Argument to allow passing this function to sc_hash_new. 
 * 
 * \return 1 if the two have the same rank, 0 if not.
 */
static int
t8_ghost_remote_equal_function (const void *remote_dataa, const void *remote_datab, [[maybe_unused]] const void *user)
{
  const t8_ghost_remote_t *remotea = (const t8_ghost_remote_t *) remote_dataa;
  const t8_ghost_remote_t *remoteb = (const t8_ghost_remote_t *) remote_datab;

  return remotea->remote_rank == remoteb->remote_rank;
}

/** 
 * This struct is used during a ghost data exchange.
 * Since we use asynchronous communication, we store the
 * send buffers and mpi requests until we end the communication.
 */
using t8_ghost_data_exchange_t = struct
{
  int num_remotes;               /**< The number of processes, we send to. */
  char **send_buffers;           /**< For each remote the send buffer. */
  sc_MPI_Request *send_requests; /**< For each process we send to, the MPI request used. */
  sc_MPI_Request *recv_requests; /**< For each process we receive from, the MPI request used. */
};

void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type)
{
  t8_forest_ghost_t ghost;

  /* We currently only support face-neighbor ghosts */
  T8_ASSERT (ghost_type == T8_GHOST_FACES);

  /* Allocate memory for ghost */
  ghost = *pghost = T8_ALLOC_ZERO (t8_forest_ghost_struct_t, 1);
  /* initialize the reference counter */
  t8_refcount_init (&ghost->rc);
  /* Set the ghost type */
  ghost->ghost_type = ghost_type;

  /* Allocate the trees array */
  ghost->ghost_trees = sc_array_new (sizeof (t8_ghost_tree_t));

  /* initialize the global_tree_to_ghost_tree hash table */
  ghost->glo_tree_mempool = sc_mempool_new (sizeof (t8_ghost_gtree_hash_t));
  ghost->global_tree_to_ghost_tree
    = sc_hash_new (t8_ghost_gtree_hash_function, t8_ghost_gtree_equal_function, nullptr, nullptr);

  /* initialize the process_offset hash table */
  ghost->proc_offset_mempool = sc_mempool_new (sizeof (t8_ghost_process_hash_t));
  ghost->process_offsets
    = sc_hash_new (t8_ghost_process_hash_function, t8_ghost_process_equal_function, nullptr, nullptr);
  /* initialize the remote ghosts hash table */
  ghost->remote_ghosts = sc_hash_array_new (sizeof (t8_ghost_remote_t), t8_ghost_remote_hash_function,
                                            t8_ghost_remote_equal_function, nullptr);
  /* initialize the remote processes array */
  ghost->remote_processes = sc_array_new (sizeof (int));
}

/** 
 *  Return the remote struct of a given remote rank 
 *  
 *  \param[in] forest A committed forest.
 *  \param[in] remote The rank of the remote process.
 * 
 * \return The remote struct of the rank as \see t8_ghost_remote_t pointer.
 */
static t8_ghost_remote_t *
t8_forest_ghost_get_remote (t8_forest_t forest, int remote)
{
  t8_ghost_remote_t remote_search;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  size_t index;

  T8_ASSERT (t8_forest_is_committed (forest));

  remote_search.remote_rank = remote;
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (forest->ghosts->remote_ghosts, &remote_search, &index);
  T8_ASSERT (ret);
  return (t8_ghost_remote_t *) sc_array_index (&forest->ghosts->remote_ghosts->a, index);
}

/** 
 *  Return a remote processes info about the stored ghost elements 
 *  
 *  \param[in] forest A committed forest.
 *  \param[in] remote The rank of the remote.
 * 
 *  \return The remote process info about the stored ghost elements, as \see t8_ghost_process_hash_t.
 */
static t8_ghost_process_hash_t *
t8_forest_ghost_get_proc_info (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t proc_hash_search, **pproc_hash_found, *proc_hash_found;
#if T8_ENABLE_DEBUG
  int ret;
#endif

  T8_ASSERT (t8_forest_is_committed (forest));

  proc_hash_search.mpirank = remote;
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_lookup (forest->ghosts->process_offsets, &proc_hash_search, (void ***) &pproc_hash_found);
  T8_ASSERT (ret);
  proc_hash_found = *pproc_hash_found;
  T8_ASSERT (proc_hash_found->mpirank == remote);
  return proc_hash_found;
}

/* Return the number of trees in a ghost */
t8_locidx_t
t8_forest_ghost_num_trees (const t8_forest_t forest)
{
  if (forest->ghosts == nullptr) {
    return 0;
  }
  T8_ASSERT (forest->ghosts != NULL);
  if (forest->ghosts->num_ghosts_elements <= 0) {
    return 0;
  }
  T8_ASSERT (forest->ghosts->ghost_trees != NULL);

  return forest->ghosts->ghost_trees->elem_count;
}

/** 
 * Given an index into the ghost_trees array return the ghost tree 
 * 
 * \param[in] forest      A committed forest.
 * \param[in] lghost_tree Index of the tree within the ghost_trees array.
 * 
 * \return The ghost tree.
 */
static t8_ghost_tree_t *
t8_forest_ghost_get_tree (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  t8_forest_ghost_t ghost;

  T8_ASSERT (t8_forest_is_committed (forest));
  ghost = forest->ghosts;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->ghost_trees != NULL);
  T8_ASSERT (0 <= lghost_tree && lghost_tree < t8_forest_ghost_num_trees (forest));

  ghost_tree = (t8_ghost_tree_t *) t8_sc_array_index_locidx (ghost->ghost_trees, lghost_tree);
  return ghost_tree;
}

t8_locidx_t
t8_forest_ghost_get_tree_element_offset (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return t8_forest_ghost_get_tree (forest, lghost_tree)->element_offset;
}

/** Given an index in the ghost_tree array, return this tree's number of elements */
t8_locidx_t
t8_forest_ghost_tree_num_leaf_elements (t8_forest_t forest, t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return t8_element_array_get_count (&ghost_tree->elements);
}

t8_element_array_t *
t8_forest_ghost_get_tree_leaf_elements (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  return &t8_forest_ghost_get_tree (forest, lghost_tree)->elements;
}

t8_locidx_t
t8_forest_ghost_get_ghost_treeid (t8_forest_t forest, t8_gloidx_t gtreeid)
{
  t8_ghost_gtree_hash_t query, *found, **pfound;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  query.global_id = gtreeid;
  if (sc_hash_lookup (forest->ghosts->global_tree_to_ghost_tree, &query, (void ***) &pfound)) {
    /* The tree was found */
    found = *pfound;
    return found->index;
  }
  else {
    /* The tree was not found */
    return -1;
  }
}

/* Given an index in the ghost_tree array, return this tree's element class */
t8_eclass_t
t8_forest_ghost_get_tree_class (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->eclass;
}

/* Given an index in the ghost_tree array, return this tree's global id */
t8_gloidx_t
t8_forest_ghost_get_global_treeid (const t8_forest_t forest, const t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->global_id;
}

/* Given an index into the ghost_trees array and for that tree an element index,
   return the corresponding element. */
t8_element_t *
t8_forest_ghost_get_leaf_element (t8_forest_t forest, t8_locidx_t lghost_tree, t8_locidx_t lelement)
{
  t8_ghost_tree_t *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  T8_ASSERT (0 <= lelement && lelement < t8_forest_ghost_tree_num_leaf_elements (forest, lghost_tree));
  /* TODO: In future, make return type const (and offer additional mutable version) and call t8_element_array_index_locidx (the const version). */
  return t8_element_array_index_locidx_mutable (&ghost_tree->elements, lelement);
}

/** Initialize a t8_ghost_remote_tree_t.
 * 
 *  \param[in]  forest    The forest.
 *  \param[in]  gtreeid   The global ID of the remote tree.
 *  \param[in]  remote_rank The rank of the reomte process holding the tree.
 *  \param[in]  tree_class  The eclass of the remote tree.
 *  \param[in, out] remote_tree A pointer to the t8_ghost_remote_tree_t to be initialized. 
 *                              Has to be non-NULL on input. On output, it is initialized.
*/
static void
t8_ghost_init_remote_tree (t8_forest_t forest, t8_gloidx_t gtreeid, int remote_rank, t8_eclass_t tree_class,
                           t8_ghost_remote_tree_t *remote_tree)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_locidx_t local_treeid;

  T8_ASSERT (remote_tree != NULL);

  local_treeid = gtreeid - t8_forest_get_first_local_tree_id (forest);
  /* Set the entries of the new remote tree */
  remote_tree->global_id = gtreeid;
  remote_tree->mpirank = remote_rank;
  remote_tree->eclass = t8_forest_get_eclass (forest, local_treeid);
  /* Initialize the array to store the element */
  t8_element_array_init (&remote_tree->elements, scheme, tree_class);
  /* Initialize the array to store the element indices. */
  sc_array_init (&remote_tree->element_indices, sizeof (t8_locidx_t));
}

/** 
 * Add a new element to the remote hash table (if not already in it).
 * Must be called for elements in linear order
 * element_index is the tree local index of this element.
 * 
 * \param[in] forest          The forest.
 * \param[in] ghost           The ghost structure.
 * \param[in] remote_rank     The remote rank.
 * \param[in] ltreeid         Local id of the tree within the forest.
 * \param[in] elem            The element to be added.
 * \param[in] element_index   The element's tree-local id.
 */
static void
t8_ghost_add_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int remote_rank, t8_locidx_t ltreeid,
                     const t8_element_t *elem, t8_locidx_t element_index)
{
  t8_ghost_remote_t remote_entry_lookup, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_element_t *elem_copy;
  sc_array_t *remote_array;
  size_t index, element_count;
  int *remote_process_entry;
  int level, copy_level = 0;

  /* Get the tree's element class and the scheme */
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_gloidx_t gtreeid = t8_forest_get_first_local_tree_id (forest) + ltreeid;

  /* Check whether the remote_rank is already present in the remote ghosts
   * array. */
  remote_entry_lookup.remote_rank = remote_rank;
  /* clang-format off */
  remote_entry = (t8_ghost_remote_t *) sc_hash_array_insert_unique (ghost->remote_ghosts, (void *) &remote_entry_lookup,
                                                                    &index);
  /* clang-format on */
  if (remote_entry != nullptr) {
    /* The remote rank was not in the array and was inserted now */
    remote_entry->remote_rank = remote_rank;
    remote_entry->num_elements = 0;
    /* Initialize the tree array of the new entry */
    sc_array_init_size (&remote_entry->remote_trees, sizeof (t8_ghost_remote_tree_t), 1);
    /* Get a pointer to the new entry */
    remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees, 0);
    /* initialize the remote_tree */
    t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, tree_class, remote_tree);
    /* Since the rank is a new remote rank, we also add it to the remote ranks array */
    remote_process_entry = (int *) sc_array_push (ghost->remote_processes);
    *remote_process_entry = remote_rank;
  }
  else {
    /* The remote rank already is contained in the remotes array at position index. */
    remote_array = &ghost->remote_ghosts->a;
    remote_entry = (t8_ghost_remote_t *) sc_array_index (remote_array, index);
    T8_ASSERT (remote_entry->remote_rank == remote_rank);
    /* Check whether the tree has already an entry for this process.
     * Since we only add in local tree order the current tree is either
     * the last entry or does not have an entry yet. */
    remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees,
                                                             remote_entry->remote_trees.elem_count - 1);
    if (remote_tree->global_id != gtreeid) {
      /* The tree does not exist in the array. We thus need to add it and
       * initialize it. */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_push (&remote_entry->remote_trees);
      t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, tree_class, remote_tree);
    }
  }
  /* remote_tree now points to a valid entry for the tree.
   * We can add a copy of the element to the elements array
   * if it does not exist already. If it exists it is the last entry in the array. */
#if T8_ENABLE_DEBUG
  {
    /* debugging assertion that the element is really not contained already */
    int ielem;
    int elem_count = t8_element_array_get_count (&remote_tree->elements);
    for (ielem = 0; ielem < elem_count - 1; ielem++) {
      const t8_element_t *test_el = t8_element_array_index_int (&remote_tree->elements, ielem);
      SC_CHECK_ABORTF (!scheme->element_is_equal (tree_class, test_el, elem),
                       "Local element %i already in remote ghosts at pos %i\n", element_index, ielem);
    }
  }
#endif
  elem_copy = nullptr;
  level = scheme->element_get_level (tree_class, elem);
  element_count = t8_element_array_get_count (&remote_tree->elements);
  if (element_count > 0) {
    elem_copy = t8_element_array_index_locidx_mutable (&remote_tree->elements, element_count - 1);
    copy_level = scheme->element_get_level (tree_class, elem_copy);
  }
  /* Check if the element was not contained in the array.
   * If so, we add a copy of elem to the array.
   * Otherwise, we do nothing. */
  if (elem_copy == nullptr || level != copy_level
      || scheme->element_get_linear_id (tree_class, elem_copy, copy_level)
           != scheme->element_get_linear_id (tree_class, elem, level)) {
    /* Add the element */
    elem_copy = t8_element_array_push (&remote_tree->elements);
    scheme->element_copy (tree_class, elem, elem_copy);
    /* Add the index of the element */
    *(t8_locidx_t *) sc_array_push (&remote_tree->element_indices) = element_index;
    remote_entry->num_elements++;
  }
}

/**
 * This struct stores the ghost boundary data of a forest.
 */
using t8_forest_ghost_boundary_data_t = struct
{

  sc_array_t bounds_per_level; /**< For each level from the nca to the parent of the current element
                                 *   we store for each face the lower and upper bounds of the owners at
                                 *   this face. We also store bounds for the element's owners.
                                 *   Each entry is an array of 2 * (max_num_faces + 1) integers,
      |                          *   face_0 low | face_0 high | ... | face_n low | face_n high | owner low | owner high | 
                                 */
  sc_array_t face_owners;      /**< Temporary storage for all owners at a leaf's face. */
  const t8_scheme *scheme;     /**< The scheme. */
  t8_gloidx_t gtreeid;         /**< The global tree id. */

  int level_nca;      /**< The refinement level of the root element in the search.
                                     At position element_level - level_nca in bounds_per_level are the bounds
                                     for the parent of element. */
  int max_num_faces;  /**< The maximum number of faces. */
  t8_eclass_t eclass; /**< The element class. */
#if T8_ENABLE_DEBUG
  t8_locidx_t left_out; /**< Count the elements for which we skip the search */
#endif
};

/**
 *  This function is used as callback search function within \ref t8_forest_search to check whether the neighbors of
 *  the current element are on another rank. If so, add the element to the ghost structures.
 *  
 *  \param[in] forest           The forest.
 *  \param[in] ltreeid          Local ID of the ghost tree.
 *  \param[in] element          The current element.
 *  \param[in] is_leaf          Switch indicating whether \a element is a leaf.
 *  \param[in] leaves           The array of leaves (used for debug purposes only).
 *  \param[in] tree_leaf_index  If the element is a leaf, its tree-local id.
 * 
 *  \return 0 if the element and its face neighbors are completely owned by the current rank; 1 otherwise.
 */
static int
t8_forest_ghost_search_boundary (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 const int is_leaf, [[maybe_unused]] const t8_element_array_t *leaves,
                                 const t8_locidx_t tree_leaf_index)
{
  t8_forest_ghost_boundary_data_t *data = (t8_forest_ghost_boundary_data_t *) t8_forest_get_user_data (forest);
  int num_faces, iface, faces_totally_owned, level;
  int parent_face;
  int lower, upper, *bounds, *new_bounds, parent_lower, parent_upper;
  int el_lower, el_upper;
  int element_is_owned, iproc, remote_rank;

  /* First part: the search enters a new tree, we need to reset the user_data */
  if (t8_forest_global_tree_id (forest, ltreeid) != data->gtreeid) {
    int max_num_faces;
    /* The search has entered a new tree, store its eclass and element scheme */
    data->gtreeid = t8_forest_global_tree_id (forest, ltreeid);
    data->eclass = t8_forest_get_eclass (forest, ltreeid);
    data->scheme = t8_forest_get_scheme (forest);
    data->level_nca = data->scheme->element_get_level (data->eclass, element);
    data->max_num_faces = data->scheme->element_get_max_num_faces (data->eclass, element);
    max_num_faces = data->max_num_faces;
    sc_array_reset (&data->bounds_per_level);
    sc_array_init_size (&data->bounds_per_level, 2 * (max_num_faces + 1) * sizeof (int), 1);
    /* Set the (imaginary) owner bounds for the parent of the root element */
    bounds = (int *) sc_array_index (&data->bounds_per_level, 0);
    for (iface = 0; iface < max_num_faces + 1; iface++) {
      bounds[iface * 2] = 0;
      bounds[iface * 2 + 1] = forest->mpisize - 1;
    }
    /* TODO: compute bounds */
  }

  /* The level of the current element */
  level = data->scheme->element_get_level (data->eclass, element);
  /* Get a pointer to the owner at face bounds of this element, if there doesnt exist
   * an entry for this in the bounds_per_level array yet, we allocate it */
  T8_ASSERT (level >= data->level_nca);
  if (data->bounds_per_level.elem_count <= (size_t) level - data->level_nca + 1) {
    T8_ASSERT (data->bounds_per_level.elem_count == (size_t) level - data->level_nca + 1);
    new_bounds = (int *) sc_array_push (&data->bounds_per_level);
  }
  else {
    new_bounds = (int *) sc_array_index (&data->bounds_per_level, level - data->level_nca + 1);
  }

  /* Get a pointer to the owner bounds of the parent */
  bounds = (int *) sc_array_index (&data->bounds_per_level, level - data->level_nca);
  /* Get bounds for the element's parent's owners */
  parent_lower = bounds[2 * data->max_num_faces];
  parent_upper = bounds[2 * data->max_num_faces + 1];
  /* Temporarily store them to serve as bounds for this element's owners */
  el_lower = parent_lower;
  el_upper = parent_upper;
  /* Compute bounds for the element's owners */
  t8_forest_element_owners_bounds (forest, data->gtreeid, element, data->eclass, &el_lower, &el_upper);
  /* Set these as the new bounds */
  new_bounds[2 * data->max_num_faces] = el_lower;
  new_bounds[2 * data->max_num_faces + 1] = el_upper;
  element_is_owned = (el_lower == el_upper);
  num_faces = data->scheme->element_get_num_faces (data->eclass, element);
  faces_totally_owned = 1;

  /* TODO: we may not carry on with the face computations if the element is not
   *       totally owned and immediately return 1. However, how do we set the bounds for
   *       the face owners then?
   */
  for (iface = 0; iface < num_faces; iface++) {
    /* Compute the face number of the parent to reuse the bounds */
    parent_face = data->scheme->element_face_get_parent_face (data->eclass, element, iface);
    if (parent_face >= 0) {
      /* This face was also a face of the parent, we reuse the computed bounds */
      lower = bounds[parent_face * 2];
      upper = bounds[parent_face * 2 + 1];
    }
    else {
      /* this is an inner face, thus the face owners must be owners of the parent element */
      lower = parent_lower;
      upper = parent_upper;
    }

    if (!is_leaf) {
      /* The element is not a leaf, we compute bounds for the face neighbor owners,
       * if all face neighbors are owned by this rank, and the element is completely
       * owned, then we do not continue the search. */
      /* Compute the owners of the neighbor at this face of the element */
      t8_forest_element_owners_at_neigh_face_bounds (forest, ltreeid, element, iface, &lower, &upper);
      /* Store the new bounds at the entry for this element */
      new_bounds[iface * 2] = lower;
      new_bounds[iface * 2 + 1] = upper;
      if (lower != upper or lower != forest->mpirank) {
        faces_totally_owned = 0;
      }
    }
    else {
      /* The element is a leaf, we compute all of its face neighbor owners
       * and add the element as a remote element to all of them. */
      sc_array_resize (&data->face_owners, 2);
      /* The first and second entry in the face_owners array serve as lower and upper bound */
      *(int *) sc_array_index (&data->face_owners, 0) = lower;
      *(int *) sc_array_index (&data->face_owners, 1) = upper;
      t8_forest_element_owners_at_neigh_face (forest, ltreeid, element, iface, &data->face_owners);
      /*TODO: add as remotes */
      for (iproc = 0; iproc < (int) data->face_owners.elem_count; iproc++) {
        remote_rank = *(int *) sc_array_index (&data->face_owners, iproc);
        if (remote_rank != forest->mpirank) {
          t8_ghost_add_remote (forest, forest->ghosts, remote_rank, ltreeid, element, tree_leaf_index);
        }
      }
    }
  } /* end face loop */
  if (faces_totally_owned && element_is_owned) {
    /* The element only has local descendants and all of its face neighbors are local as well. 
     * We do not continue the search */
#if T8_ENABLE_DEBUG
    if (tree_leaf_index < 0) {
      data->left_out += t8_element_array_get_count (leaves);
    }
#endif
    return 0;
  }
  /* Continue the top-down search if this element or its face neighbors are not completely owned by the rank. */
  return 1;
}

/** 
 * Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 * 
 * \param[in,out] forest  the forest.
 */
static void
t8_forest_ghost_fill_remote_v3 (t8_forest_t forest)
{
  t8_forest_ghost_boundary_data_t data;
  void *store_user_data = nullptr;

  /* Start with invalid entries in the user data.
   * These are set in t8_forest_ghost_search_boundary each time a new tree is entered */
  data.eclass = T8_ECLASS_COUNT;
  data.gtreeid = -1;
  data.scheme = nullptr;
#if T8_ENABLE_DEBUG
  data.left_out = 0;
#endif
  sc_array_init (&data.face_owners, sizeof (int));
  /* This is a dummy init, since we call sc_array_reset in ghost_search_boundary
   * and we should not call sc_array_reset on a non-initialized array */
  sc_array_init (&data.bounds_per_level, 1);
  /* Store any user data that may reside on the forest */
  store_user_data = t8_forest_get_user_data (forest);
  /* Set the user data for the search routine */
  t8_forest_set_user_data (forest, &data);
  /* Loop over the trees of the forest */
  t8_forest_search (forest, t8_forest_ghost_search_boundary, nullptr, nullptr);

  /* Reset the user data from before search */
  t8_forest_set_user_data (forest, store_user_data);

  /* Reset the data arrays */
  sc_array_reset (&data.face_owners);
  sc_array_reset (&data.bounds_per_level);
#if T8_ENABLE_DEBUG
#endif
}

/** 
 * Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 * If ghost_method is 0, then we assume a balanced forest and
 * construct the remote processes by looking at the half neighbors of an element.
 * Otherwise, we use the owners_at_face method.
 * 
 * \param[in] forest        The forest.
 * \param[in] ghost         The forest's ghost.
 * \param[in] ghost_method  Switch indicating the ghost type.
 */
static void
t8_forest_ghost_fill_remote (t8_forest_t forest, t8_forest_ghost_t ghost, int ghost_method)
{
  t8_element_t **half_neighbors = nullptr;
  t8_locidx_t num_local_trees, num_tree_elems;
  t8_locidx_t itree, ielem;
  t8_tree_t tree;
  t8_eclass_t last_class;
  t8_gloidx_t neighbor_tree;

  int iface, num_faces;
  int num_face_children, max_num_face_children = 0;
  int ichild, owner;
  sc_array_t owners, tree_owners;
  int is_atom;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  last_class = T8_ECLASS_COUNT;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  if (ghost_method != 0) {
    sc_array_init (&owners, sizeof (int));
    sc_array_init (&tree_owners, sizeof (int));
  }

  /* Loop over the trees of the forest */
  for (itree = 0; itree < num_local_trees; itree++) {
    /* Get a pointer to the tree, the class of the tree, the
     * scheme associated to the class and the number of elements in this tree. */
    tree = t8_forest_get_tree (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);

    /* Loop over the elements of this tree */
    num_tree_elems = t8_forest_get_tree_leaf_element_count (tree);
    for (ielem = 0; ielem < num_tree_elems; ielem++) {
      /* Get the element of the tree */
      const t8_element_t *elem = t8_forest_get_tree_leaf_element (tree, ielem);
      num_faces = scheme->element_get_num_faces (tree_class, elem);
      if (scheme->element_get_level (tree_class, elem) == scheme->get_maxlevel (tree_class)) {
        /* flag to decide whether this element is at the maximum level */
        is_atom = 1;
      }
      else {
        is_atom = 0;
      }
      for (iface = 0; iface < num_faces; iface++) {
        /* TODO: Check whether the neighbor element is inside the forest,
         *       if not then do not compute the half_neighbors.
         *       This will save computing time. Needs an "element is in forest" function
         *       Currently we perform this check in the half_neighbors function. */

        /* Get the element class of the neighbor tree */
        const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest, itree, elem, iface);
        if (ghost_method == 0) {
          /* Use half neighbors */
          /* Get the number of face children of the element at this face */
          num_face_children = scheme->element_get_num_face_children (tree_class, elem, iface);
          /* regrow the half_neighbors array if necessary.
           * We also need to reallocate it, if the element class of the neighbor
           * changes */
          if (max_num_face_children < num_face_children || last_class != neigh_class) {
            half_neighbors = T8_ALLOC (t8_element_t *, num_face_children);
            /* Allocate memory for the half size face neighbors */
            scheme->element_new (neigh_class, num_face_children, half_neighbors);
            max_num_face_children = num_face_children;
            last_class = neigh_class;
          }
          if (!is_atom) {
            /* Construct each half size neighbor */
            neighbor_tree = t8_forest_element_half_face_neighbors (forest, itree, elem, half_neighbors, neigh_class,
                                                                   iface, num_face_children, nullptr);
          }
          else {
            int dummy_neigh_face;
            /* This element has maximum level, we only construct its neighbor */
            neighbor_tree = t8_forest_element_face_neighbor (forest, itree, elem, half_neighbors[0], neigh_class, iface,
                                                             &dummy_neigh_face);
          }
          if (neighbor_tree >= 0) {
            /* If there exist face neighbor elements (we are not at a domain boundary */
            /* Find the owner process of each face_child */
            for (ichild = 0; ichild < num_face_children; ichild++) {
              /* find the owner */
              owner = t8_forest_element_find_owner (forest, neighbor_tree, half_neighbors[ichild], neigh_class);
              T8_ASSERT (0 <= owner && owner < forest->mpisize);
              if (owner != forest->mpirank) {
                /* Add the element as a remote element */
                t8_ghost_add_remote (forest, ghost, owner, itree, elem, ielem);
              }
            }
          }
          scheme->element_destroy (neigh_class, num_face_children, half_neighbors);
          T8_FREE (half_neighbors);
        } /* end ghost_method 0 */
        else {
          size_t iowner;
          /* Construct the owners at the face of the neighbor element */
          t8_forest_element_owners_at_neigh_face (forest, itree, elem, iface, &owners);
          /* Iterate over all owners and if any is not the current process,
           * add this element as remote */
          for (iowner = 0; iowner < owners.elem_count; iowner++) {
            owner = *(int *) sc_array_index (&owners, iowner);
            T8_ASSERT (0 <= owner && owner < forest->mpisize);
            if (owner != forest->mpirank) {
              /* Add the element as a remote element */
              t8_ghost_add_remote (forest, ghost, owner, itree, elem, ielem);
            }
          }
          sc_array_truncate (&owners);
        }
      } /* end face loop */
    }   /* end element loop */
  }     /* end tree loop */

  if (forest->profile != nullptr) {
    /* If profiling is enabled, we count the number of remote processes. */
    forest->profile->ghosts_remotes = ghost->remote_processes->elem_count;
  }
  /* Clean-up memory */
  if (ghost_method != 0) {
    sc_array_reset (&owners);
    sc_array_reset (&tree_owners);
  }
}

/** 
 * Begin sending the ghost elements from the remote ranks
 * using non-blocking communication.
 * Afterwards,
 *  t8_forest_ghost_send_end
 * must be called to end the communication.
 * 
 * \param[in]  forest     The forest.
 * \param[in]  ghost      The forest's ghost.
 * \param[out] requests   The send requests as an array of pointers to sc_MPI_Requests.
 * 
 * \return An array of mpi_send_info_t, holding one entry for each remote rank.
 */
static t8_ghost_mpi_send_info_t *
t8_forest_ghost_send_start (t8_forest_t forest, t8_forest_ghost_t ghost, sc_MPI_Request **requests)
{
  int proc_index, remote_rank;
  int num_remotes;
  size_t remote_index;
  t8_ghost_remote_t *remote_entry;
  sc_array_t *remote_trees;
  t8_ghost_remote_tree_t *remote_tree = nullptr;
  t8_ghost_mpi_send_info_t *send_info, *current_send_info;
  char *current_buffer;
  size_t bytes_written, element_bytes, element_count, element_size;
#if T8_ENABLE_DEBUG
  size_t acc_el_count = 0;
#endif
  int mpiret;

  /* Allocate a send_buffer for each remote rank */
  num_remotes = ghost->remote_processes->elem_count;
  send_info = T8_ALLOC (t8_ghost_mpi_send_info_t, num_remotes);
  *requests = T8_ALLOC (sc_MPI_Request, num_remotes);

  /* Loop over all remote processes */
  for (proc_index = 0; proc_index < (int) ghost->remote_processes->elem_count; proc_index++) {
    current_send_info = send_info + proc_index;
    /* Get the rank of the current remote process. */
    remote_rank = *(int *) sc_array_index_int (ghost->remote_processes, proc_index);
    t8_debugf ("Filling send buffer for process %i\n", remote_rank);
    /* initialize the send_info for the current rank */
    current_send_info->recv_rank = remote_rank;
    current_send_info->num_bytes = 0;
    current_send_info->request = *requests + proc_index;
    /* Lookup the ghost elements for the first tree of this remote */
    remote_entry = t8_forest_ghost_get_remote (forest, remote_rank);
    T8_ASSERT (remote_entry->remote_rank == remote_rank);
    /* Loop over all trees of the remote rank and count the bytes */
    /* At first we store the number of remote trees in the buffer */
    current_send_info->num_bytes += sizeof (size_t);
    /* add padding before the eclass */
    current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
    /* TODO: put this in a function */
    /* TODO: Use remote_entry to count the number of bytes while inserting
     *        the remote ghosts. */
    remote_trees = &remote_entry->remote_trees;
    for (remote_index = 0; remote_index < remote_trees->elem_count; remote_index++) {
      /* Get the next remote tree. */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (remote_trees, remote_index);
      /* We will store the global tree id, the element class and the list
       * of elements in the send_buffer. */
      current_send_info->num_bytes += sizeof (t8_gloidx_t);
      /* add padding before the eclass */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += sizeof (t8_eclass_t);
      /* add padding before the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_count = t8_element_array_get_count (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* We will store the number of elements */
      current_send_info->num_bytes += sizeof (size_t);
      /* add padding before the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += element_bytes;
      /* add padding after the elements */
      current_send_info->num_bytes += T8_ADD_PADDING (current_send_info->num_bytes);
    }

    /* We now now the number of bytes for our send_buffer and thus allocate it. */
    current_send_info->buffer = T8_ALLOC_ZERO (char, current_send_info->num_bytes);

    /* We iterate through the tree again and store the tree info and the elements into the send_buffer. */
    current_buffer = current_send_info->buffer;
    bytes_written = 0;
    /* Start with the number of remote trees in the buffer */
    memcpy (current_buffer + bytes_written, &remote_trees->elem_count, sizeof (size_t));
    bytes_written += sizeof (size_t);
    bytes_written += T8_ADD_PADDING (bytes_written);
#if T8_ENABLE_DEBUG
    acc_el_count = 0;
#endif
    for (remote_index = 0; remote_index < remote_trees->elem_count; remote_index++) {
      /* Get a pointer to the tree */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (remote_trees, remote_index);
      T8_ASSERT (remote_tree->mpirank == remote_rank);

      /* Copy the global tree id */
      memcpy (current_buffer + bytes_written, &remote_tree->global_id, sizeof (t8_gloidx_t));
      bytes_written += sizeof (t8_gloidx_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Copy the trees element class */
      memcpy (current_buffer + bytes_written, &remote_tree->eclass, sizeof (t8_eclass_t));
      bytes_written += sizeof (t8_eclass_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Store the number of elements in the buffer */
      element_count = t8_element_array_get_count (&remote_tree->elements);
      memcpy (current_buffer + bytes_written, &element_count, sizeof (size_t));
      bytes_written += sizeof (size_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* Copy the elements into the send buffer */
      memcpy (current_buffer + bytes_written, t8_element_array_get_data (&remote_tree->elements), element_bytes);
      bytes_written += element_bytes;
      /* add padding after the elements */
      bytes_written += T8_ADD_PADDING (bytes_written);

      /* Add to the counter of remote elements. */
      ghost->num_remote_elements += element_count;
#if T8_ENABLE_DEBUG
      acc_el_count += element_count;
#endif
    } /* End tree loop */

    T8_ASSERT (bytes_written == current_send_info->num_bytes);
    /* We can now post the MPI_Isend for the remote process */
    mpiret = sc_MPI_Isend (current_buffer, bytes_written, sc_MPI_BYTE, remote_rank, T8_MPI_GHOST_FOREST,
                           forest->mpicomm, *requests + proc_index);
    SC_CHECK_MPI (mpiret);
  } /* end process loop */
  return send_info;
}

/**
 *  End the communication of the ghost element sends and receives.
 * 
 *  \param[in] forest          A committed forest. 
 *  \param[in] ghost           The ghost data.
 *  \param[in, out] send_info  The send information array. On output, its memory is freed.
 *  \param[in, out] requests   The array of MPI requests. On output, its memory is freed.
 * 
 */
static void
t8_forest_ghost_send_end ([[maybe_unused]] t8_forest_t forest, t8_forest_ghost_t ghost,
                          t8_ghost_mpi_send_info_t *send_info, sc_MPI_Request *requests)
{
  int num_remotes;
  int proc_pos, mpiret;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (ghost != NULL);

  /* Get the number of remote processes */
  num_remotes = ghost->remote_processes->elem_count;

  /* We wait for all communication to end. */
  mpiret = sc_MPI_Waitall (num_remotes, requests, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* Clean-up */
  for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
    T8_FREE (send_info[proc_pos].buffer);
  }
  T8_FREE (send_info);
  T8_FREE (requests);
}

/**
 * Receive a single message from a remote process, after the message was successfully probed.
 * Returns the allocated receive buffer and the number of bytes received 
 * 
 * \param[in] recv_rank   The receiving rank.
 * \param[in] comm        The MPI communicator.
 * \param[in] status      The sc_MPI_Status.
 * \param[in] recv_bytes  The number of bytes to be received.
 * 
 * \return The allocated receive buffer and the number of bytes received.
 */
static char *
t8_forest_ghost_receive_message (int recv_rank, sc_MPI_Comm comm, sc_MPI_Status status, int *recv_bytes)
{
  char *recv_buffer;
  int mpiret;

  T8_ASSERT (recv_rank == status.MPI_SOURCE);
  T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);

  /* Get the number of bytes in the message */
  mpiret = sc_MPI_Get_count (&status, sc_MPI_BYTE, recv_bytes);

  /* Allocate receive buffer */
  recv_buffer = T8_ALLOC_ZERO (char, *recv_bytes);
  /* receive the message */
  mpiret
    = sc_MPI_Recv (recv_buffer, *recv_bytes, sc_MPI_BYTE, recv_rank, T8_MPI_GHOST_FOREST, comm, sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);

  return recv_buffer;
}

/**
 * Parse a message from a remote process and correctly include the received
 * elements in the ghost structure.
 * The message looks like:
 * num_trees | pad | treeid 0 | pad | eclass 0 | pad | num_elems 0 | pad | elements | pad | treeid 1 | ...
 *  size_t   |     |t8_gloidx |     |t8_eclass |     | size_t      |     | t8_element_t |
 *
 * pad is paddind, see T8_ADD_PADDING
 *
 * current_element_offset is updated in each step to store the element offset
 * of the next ghost tree to be inserted.
 * When called with the first message, current_element_offset must be set to 0.
 *
 * Currently we expect that the messages arrive in order of the sender's rank. 
 * 
 * \param[in] forest                      The forest.
 * \param[in] ghost                       The forest's ghost.
 * \param[in, out] current_element_offset The current element offset. Has to be zero on input and is updated in each step.
 * \param[in] recv_rank                   The receiving rank.
 * \param[in] recv_buffer                 The receive buffer.
 * \param[in] recv_bytes                  The number of bytes received.
 */
static void
t8_forest_ghost_parse_received_message (t8_forest_t forest, t8_forest_ghost_t ghost,
                                        t8_locidx_t *current_element_offset, int recv_rank, char *recv_buffer,
                                        int recv_bytes)
{
  size_t bytes_read, first_tree_index = 0, first_element_index = 0;
  t8_locidx_t num_trees, itree;
  t8_gloidx_t global_id;
  t8_eclass_t eclass;
  size_t num_elements, old_elem_count, ghosts_offset;
  t8_ghost_gtree_hash_t *tree_hash, **pfound_tree, *found_tree;
  t8_ghost_tree_t *ghost_tree;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_element_t *element_insert;
  t8_ghost_process_hash_t *process_hash;
#if T8_ENABLE_DEBUG
  int added_process;
#endif

  bytes_read = 0;
  /* read the number of trees */
  num_trees = *(size_t *) recv_buffer;
  bytes_read += sizeof (size_t);
  bytes_read += T8_ADD_PADDING (bytes_read);

  t8_debugf ("Received %li trees from %i (%i bytes)\n", (long) num_trees, recv_rank, recv_bytes);

  /* Count the total number of ghosts that we receive from this rank */
  ghosts_offset = ghost->num_ghosts_elements;
  for (itree = 0; itree < num_trees; itree++) {
    /* Get tree id */
    /* search if tree was inserted */
    /* if not: new entry, add elements */
    /* if yes: add the elements to the end of the tree's element array. */

    /* read the global id of this tree. */
    global_id = *(t8_gloidx_t *) (recv_buffer + bytes_read);
    bytes_read += sizeof (t8_gloidx_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* read the element class of the tree */
    eclass = *(t8_eclass_t *) (recv_buffer + bytes_read);
    bytes_read += sizeof (t8_eclass_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* read the number of elements sent */
    num_elements = *(size_t *) (recv_buffer + bytes_read);

    /* Add to the counter of ghost elements. */
    ghost->num_ghosts_elements += num_elements;

    bytes_read += sizeof (size_t);
    bytes_read += T8_ADD_PADDING (bytes_read);
    /* Search for the tree in the ghost_trees array */
    tree_hash = (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    tree_hash->global_id = global_id;

    /* Get the scheme for this tree */
    if (sc_hash_insert_unique (ghost->global_tree_to_ghost_tree, tree_hash, (void ***) &pfound_tree)) {
      /* The tree was not stored already, tree_hash is now an entry in the hash table. */
      /* If the tree was not contained, it is the newest tree in the array and
       * thus has as index the number of currently inserted trees. */
      tree_hash->index = ghost->ghost_trees->elem_count;
      found_tree = tree_hash;
      /* We grow the array by one and initialize the entry */
      ghost_tree = (t8_ghost_tree_t *) sc_array_push (ghost->ghost_trees);
      ghost_tree->global_id = global_id;
      ghost_tree->eclass = eclass;
      /* Initialize the element array */
      t8_element_array_init_size (&ghost_tree->elements, scheme, eclass, num_elements);
      /* pointer to where the elements are to be inserted */
      element_insert = t8_element_array_get_data_mutable (&ghost_tree->elements);
      /* Compute the element offset of this new tree by adding the offset
       * of the previous tree to the element count of the previous tree. */
      ghost_tree->element_offset = *current_element_offset;
      /* Allocate a new tree_hash for the next search */
      old_elem_count = 0;
      tree_hash = (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    }
    else {
      /* The entry was found in the trees array */
      found_tree = *pfound_tree;
      T8_ASSERT (found_tree->global_id == global_id);
      /* Get a pointer to the tree */
      ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees, found_tree->index);
      T8_ASSERT (ghost_tree->eclass == eclass);
      T8_ASSERT (ghost_tree->global_id == global_id);
      T8_ASSERT (ghost_tree->elements.tree_class == eclass);

      old_elem_count = t8_element_array_get_count (&ghost_tree->elements);

      /* Grow the elements array of the tree to fit the new elements */
      t8_element_array_resize (&ghost_tree->elements, old_elem_count + num_elements);
      /* Get a pointer to where the new elements are to be inserted */
      element_insert = t8_element_array_index_locidx_mutable (&ghost_tree->elements, old_elem_count);
    }

    if (itree == 0) {
      /* We store the index of the first tree and the first element of this
       * rank */
      first_tree_index = found_tree->index;
      first_element_index = old_elem_count;
    }
    /* Insert the new elements */
    memcpy (element_insert, recv_buffer + bytes_read, num_elements * scheme->get_element_size (eclass));

    bytes_read += num_elements * scheme->get_element_size (eclass);
    bytes_read += T8_ADD_PADDING (bytes_read);
    *current_element_offset += num_elements;
  }
  T8_ASSERT (bytes_read == (size_t) recv_bytes);
  T8_FREE (recv_buffer);

  /* At last we add the receiving rank to the ghosts process_offset hash table */
  process_hash = (t8_ghost_process_hash_t *) sc_mempool_alloc (ghost->proc_offset_mempool);
  process_hash->mpirank = recv_rank;
  process_hash->tree_index = first_tree_index;
  process_hash->first_element = first_element_index;
  process_hash->ghost_offset = ghosts_offset;
  /* Insert this rank into the hash table. We assert if the rank was not already contained. */
#if T8_ENABLE_DEBUG
  added_process =
#else
  (void)
#endif
    sc_hash_insert_unique (ghost->process_offsets, process_hash, nullptr);
  T8_ASSERT (added_process);
}

/**
 * In forest_ghost_receive we need a lookup table to give us the position
 * of a process in the ghost->remote_processes array, given the rank of a process.
 * We implement this via a hash table with the following struct as entry. 
 */
using t8_recv_list_entry_t = struct t8_recv_list_entry_struct
{
  int rank;                    /**< The rank of this process */
  int pos_in_remote_processes; /**< The position of this process in the remote_processes array */
};

/** We hash these entries by their rank. */
unsigned
t8_recv_list_entry_hash (const void *v1, [[maybe_unused]] const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;

  return e1->rank;
}

/** Two entries are considered equal if they have the same rank. */
int
t8_recv_list_entry_equal (const void *v1, const void *v2, [[maybe_unused]] const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;
  const t8_recv_list_entry_t *e2 = (const t8_recv_list_entry_t *) v2;

  return e1->rank == e2->rank;
}

/**
 * Probe for all incoming messages from the remote ranks and receive them.
 * We receive the message in the order in which they arrive. To achieve this,
 * we have to use polling. 
 * 
 * \param[in] forest  The forest.
 * \param[in] ghost   The forest's ghost.
 */
static void
t8_forest_ghost_receive (t8_forest_t forest, t8_forest_ghost_t ghost)
{
  int num_remotes;
  int proc_pos;
  int recv_rank;
  int mpiret;
  sc_MPI_Comm comm;
  sc_MPI_Status status;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (ghost != NULL);

  comm = forest->mpicomm;
  /* Get the number of remote processes */
  num_remotes = ghost->remote_processes->elem_count;

  if (num_remotes == 0) {
    /* There is nothing to do */
    return;
  }

  {
    /**
     * This code receives the message in order of their arrival.
     * This is effective in terms of runtime, but makes it more difficult
     * to store the received data, since the data has to be stored in order of
     * ascending ranks.
     * We include the received data into the ghost structure in order of the
     * ranks of the receivers and we do this as soon as the message from
     * the next rank that we can include was received.
     */
    char **buffer;
    int *recv_bytes;
    int received_messages = 0;
    int *received_flag;
    int last_rank_parsed = -1, parse_it;
#undef T8_POLLING /* activates polling for Mpi messages, currently used for testing */
#ifdef T8_POLLING
    sc_link_t *proc_it, *prev;
    int iprobe_flag;
    sc_list_t *receivers;
#else
    t8_recv_list_entry_t **pfound, *found;
#if T8_ENABLE_DEBUG
    int ret;
#endif
    sc_hash_t *recv_list_entries_hash;
#endif
    t8_recv_list_entry_t recv_list_entry, *recv_list_entries;
    t8_locidx_t current_element_offset = 0;

    buffer = T8_ALLOC (char *, num_remotes);
    recv_bytes = T8_ALLOC (int, num_remotes);
    received_flag = T8_ALLOC_ZERO (int, num_remotes);
    recv_list_entries = T8_ALLOC (t8_recv_list_entry_t, num_remotes);

    /* Sort the array of remote processes, such that the ranks are in
     * ascending order. */
    sc_array_sort (ghost->remote_processes, sc_int_compare);

    /* We build a hash table of all ranks from which we receive and their position
     * in the remote_processes array. */
#ifdef T8_POLLING /* polling */
    receivers = sc_list_new (NULL);
#else
    recv_list_entries_hash = sc_hash_new (t8_recv_list_entry_hash, t8_recv_list_entry_equal, nullptr, nullptr);
#endif
    for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
      recv_list_entries[proc_pos].rank = *(int *) sc_array_index_int (ghost->remote_processes, proc_pos);
      recv_list_entries[proc_pos].pos_in_remote_processes = proc_pos;
#ifndef T8_POLLING
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_insert_unique (recv_list_entries_hash, recv_list_entries + proc_pos, nullptr);
      T8_ASSERT (ret == 1);
#else /* polling */
      sc_list_append (receivers, recv_list_entries + proc_pos);
#endif
    }

    /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for a message for each
     * sender and if there is one we receive it and remove the sender from the list.
     * The last message can be received via probe */
#ifdef T8_POLLING
    while (received_messages < num_remotes - 1) {
      /* TODO: This part of the code using polling and IProbe to receive the
       *       messages. We replaced with a non-polling version that uses the
       *       blocking Probe. */
      iprobe_flag = 0;
      prev = NULL; /* ensure that if the first receive entry is matched first,
                                it is removed properly. */
      for (proc_it = receivers->first; proc_it != NULL && iprobe_flag == 0;) {
        /* pointer to the rank of a receiver */
        recv_rank = ((t8_recv_list_entry_t *) proc_it->data)->rank;
        proc_pos = ((t8_recv_list_entry_t *) proc_it->data)->pos_in_remote_processes;

        mpiret = sc_MPI_Iprobe (recv_rank, T8_MPI_GHOST_FOREST, comm, &iprobe_flag, &status);
        SC_CHECK_MPI (mpiret);
#else
    while (received_messages < num_remotes) {
      /* blocking probe for a message. */
      mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_GHOST_FOREST, comm, &status);
      SC_CHECK_MPI (mpiret);
#endif
#ifdef T8_POLLING
        if (iprobe_flag == 0) {
          /* There is no message to receive, we continue */
          prev = proc_it;
          proc_it = proc_it->next;
        }
        else {
#else
      /* There is a message to receive, we receive it. */
      recv_rank = status.MPI_SOURCE;
      /* Get the position of this rank in the remote processes array */
      recv_list_entry.rank = recv_rank;
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (recv_list_entries_hash, &recv_list_entry, (void ***) &pfound);
      T8_ASSERT (ret != 0);
      found = *pfound;
      proc_pos = found->pos_in_remote_processes;
#endif
          T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);
          buffer[proc_pos] = t8_forest_ghost_receive_message (recv_rank, comm, status, recv_bytes + proc_pos);
          /* mark this entry as received. */
          T8_ASSERT (received_flag[proc_pos] == 0);
          received_flag[proc_pos] = 1;
          received_messages++;
          /* Parse all messages that we can parse now.
           * We have to parse the messages in order of their rank. */
          T8_ASSERT (last_rank_parsed < proc_pos);
          /* For all ranks that we haven't parsed yet, but can be parsed in order */
          for (parse_it = last_rank_parsed + 1; parse_it < num_remotes && received_flag[parse_it] == 1; parse_it++) {
            recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
            t8_forest_ghost_parse_received_message (forest, ghost, &current_element_offset, recv_rank, buffer[parse_it],
                                                    recv_bytes[parse_it]);
            last_rank_parsed++;
          }

#ifdef T8_POLLING /* polling */
          /* Remove the process from the list of receivers. */
          proc_it = proc_it->next;
          sc_list_remove (receivers, prev);
        }
      } /* end for */
#endif
    } /* end while */

#ifdef T8_POLLING
    /* polling */
    T8_ASSERT (receivers->elem_count == 1);
    /* Get the last rank from which we didnt receive yet */
    recv_list_entry = *(t8_recv_list_entry_t *) sc_list_pop (receivers);
    recv_rank = recv_list_entry.rank;
    proc_pos = recv_list_entry.pos_in_remote_processes;
    /* destroy the list */
    sc_list_destroy (receivers);
    /* Blocking probe for the last message */
    mpiret = sc_MPI_Probe (recv_rank, T8_MPI_GHOST_FOREST, comm, &status);
    SC_CHECK_MPI (mpiret);
    /* Receive the message */
    T8_ASSERT (received_flag[proc_pos] == 0);
    buffer[proc_pos] = t8_forest_ghost_receive_message (recv_rank, comm, status, recv_bytes + proc_pos);
    received_flag[proc_pos] = 1;
    received_messages++;
    T8_ASSERT (received_messages == num_remotes);
    /* parse all messages that are left */
    /* For all ranks that we haven't parsed yet, but can be parsed in order */
    for (parse_it = last_rank_parsed + 1; parse_it < num_remotes && received_flag[parse_it] == 1; parse_it++) {
      recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
      t8_forest_ghost_parse_received_message (forest, ghost, &current_element_offset, recv_rank, buffer[parse_it],
                                              recv_bytes[parse_it]);
      last_rank_parsed++;
    }
#endif
#if T8_ENABLE_DEBUG
    for (parse_it = 0; parse_it < num_remotes; parse_it++) {
      T8_ASSERT (received_flag[parse_it] == 1);
    }
#endif
    T8_ASSERT (last_rank_parsed == num_remotes - 1);

    /* clean-up */
#ifndef T8_POLLING
    sc_hash_destroy (recv_list_entries_hash);
#endif
    T8_FREE (buffer);
    T8_FREE (received_flag);
    T8_FREE (recv_list_entries);
    T8_FREE (recv_bytes);
  }
}

/**
 * Create one layer of ghost elements, following the algorithm
 * in: p4est: Scalable Algorithms For Parallel Adaptive
 *     Mesh Refinement On Forests of Octrees
 *     C. Burstedde, L. C. Wilcox, O. Ghattas
 * for unbalanced_version = 0 (balanced forest only) or
 *     Recursive algorithms for distributed forests of octrees
 *     T. Isaac, C. Burstedde, L. C. Wilcox and O. Ghattas
 * for unbalanced_version = 1 (also unbalanced forests possible).
 *
 * version 3 with top-down search
 * for unbalanced_version = -1
 */
void
t8_forest_ghost_create_ext (t8_forest_t forest, int unbalanced_version)
{
  t8_forest_ghost_t ghost = nullptr;
  t8_ghost_mpi_send_info_t *send_info;
  sc_MPI_Request *requests;
  int create_tree_array = 0, create_gfirst_desc_array = 0;
  int create_element_array = 0;

  T8_ASSERT (t8_forest_is_committed (forest));

  t8_productionf ("Into t8_forest_ghost with %i local elements.\n", t8_forest_get_local_num_leaf_elements (forest));

  /* In parallel, check forest for deleted elements. The ghost algorithm currently
  * does not work on forests with deleted elements.
  * See also: https://github.com/DLR-AMR/t8code/issues/825
  * See also the test case: TODO Add a test case that currently fails. */
  SC_CHECK_ABORT (
    !forest->incomplete_trees || forest->mpisize == 1,
    "ERROR: Cannot compute ghost layer for forest with deleted elements (incomplete trees/holes in the mesh).\n");

  if (forest->profile != nullptr) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start ghost at %f  %f\n", sc_MPI_Wtime (), forest->profile->ghost_runtime);
  }

  if (forest->element_offsets == nullptr) {
    /* create element offset array if not done already */
    create_element_array = 1;
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == nullptr) {
    /* Create tree offset array if not done already */
    create_tree_array = 1;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == nullptr) {
    /* Create global first desc array if not done already */
    create_gfirst_desc_array = 1;
    t8_forest_partition_create_first_desc (forest);
  }

  if (t8_forest_get_local_num_leaf_elements (forest) > 0) {
    if (forest->ghost_type == T8_GHOST_NONE) {
      t8_debugf ("WARNING: Trying to construct ghosts with ghost_type NONE. "
                 "Ghost layer is not constructed.\n");
      return;
    }
    /* Currently we only support face ghosts */
    T8_ASSERT (forest->ghost_type == T8_GHOST_FACES);

    /* Initialize the ghost structure */
    t8_forest_ghost_init (&forest->ghosts, forest->ghost_type);
    ghost = forest->ghosts;

    if (unbalanced_version == -1) {
      t8_forest_ghost_fill_remote_v3 (forest);
    }
    else {
      /* Construct the remote elements and processes. */
      t8_forest_ghost_fill_remote (forest, ghost, unbalanced_version != 0);
    }

    /* Start sending the remote elements */
    send_info = t8_forest_ghost_send_start (forest, ghost, &requests);

    /* Receive the ghost elements from the remote processes */
    t8_forest_ghost_receive (forest, ghost);

    /* End sending the remote elements */
    t8_forest_ghost_send_end (forest, ghost, send_info, requests);
  }

  if (forest->profile != nullptr) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime += sc_MPI_Wtime ();
    /* We also store the number of ghosts and remotes */
    if (ghost != nullptr) {
      forest->profile->ghosts_received = ghost->num_ghosts_elements;
      forest->profile->ghosts_shipped = ghost->num_remote_elements;
    }
    else {
      forest->profile->ghosts_received = 0;
      forest->profile->ghosts_shipped = 0;
    }
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End ghost at %f  %f\n", sc_MPI_Wtime (), forest->profile->ghost_runtime);
  }

  if (create_element_array) {
    /* Free the offset memory, if created */
    t8_shmem_array_destroy (&forest->element_offsets);
  }
  if (create_tree_array) {
    /* Free the offset memory, if created */
    t8_shmem_array_destroy (&forest->tree_offsets);
  }
  if (create_gfirst_desc_array) {
    /* Free the offset memory, if created */
    t8_shmem_array_destroy (&forest->global_first_desc);
  }

  t8_productionf ("Done t8_forest_ghost with %i local elements and %i"
                  " ghost elements.\n",
                  t8_forest_get_local_num_leaf_elements (forest), t8_forest_get_num_ghosts (forest));
}

void
t8_forest_ghost_create (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->mpisize > 1) {
    /* call unbalanced version of ghost algorithm */
    t8_forest_ghost_create_ext (forest, 1);
  }
}

void
t8_forest_ghost_create_balanced_only (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->mpisize > 1) {
    /* TODO: assert that forest is balanced */
    /* Call balanced version of ghost algorithm */
    t8_forest_ghost_create_ext (forest, 0);
  }
}

void
t8_forest_ghost_create_topdown (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  t8_forest_ghost_create_ext (forest, -1);
}

/** Return the array of remote ranks.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in,out] num_remotes On output the number of remote ranks is stored here.
 * \return              The array of remote ranks in ascending order.
 */
int *
t8_forest_ghost_get_remotes (t8_forest_t forest, int *num_remotes)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->ghosts == nullptr) {
    *num_remotes = 0;
    return nullptr;
  }
  T8_ASSERT (forest->ghosts != NULL);

  *num_remotes = forest->ghosts->remote_processes->elem_count;
  return (int *) forest->ghosts->remote_processes->array;
}

/** Return the first local ghost tree of a remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The ghost tree id of the first ghost tree that stores ghost
 *                      elements of \a remote.
 */
t8_locidx_t
t8_forest_ghost_remote_first_tree (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t *proc_entry;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  proc_entry = t8_forest_ghost_get_proc_info (forest, remote);
  T8_ASSERT (proc_entry->mpirank == remote);
  return proc_entry->tree_index;
}

/** Return the local index of the first ghost element that belongs to a given remote rank.
 * \param [in] forest   A forest with constructed ghost layer.
 * \param [in] remote   A remote rank of the ghost layer in \a forest.
 * \return              The index i in the ghost elements of the first element of rank \a remote
 */
t8_locidx_t
t8_forest_ghost_remote_first_elem (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t *proc_entry;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->ghosts != NULL);

  proc_entry = t8_forest_ghost_get_proc_info (forest, remote);
  T8_ASSERT (proc_entry->mpirank == remote);

  return proc_entry->ghost_offset;
}

/**
 * Fill the send buffer for a ghost data exchange for on remote rank.
 * 
 * \param[in]   forest        The forest.
 * \param[in]   remote        The remote rank to send to.
 * \param[out]  pbuffer       The send buffer, allocated within this function.
 * \param[in]   element_data  The element data.
 * 
 * \return The number of bytes in the buffer. 
*/
static size_t
t8_forest_ghost_exchange_fill_send_buffer (t8_forest_t forest, int remote, char **pbuffer, sc_array_t *element_data)
{
  char *buffer;
  t8_ghost_remote_t lookup_rank, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_forest_ghost_t ghost;
  size_t index, element_index, data_size;
  size_t elements_inserted, byte_count;
  t8_tree_t local_tree;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  t8_locidx_t itree, ielement, element_pos;
  t8_locidx_t ltreeid;
  size_t elem_count;

  ghost = forest->ghosts;
  data_size = element_data->elem_size;
  elements_inserted = 0;
  lookup_rank.remote_rank = remote;

  /* Lookup the remote entry of this remote process */
#if T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (ghost->remote_ghosts, &lookup_rank, &index);
  T8_ASSERT (ret != 0);
  remote_entry = (t8_ghost_remote_t *) sc_array_index (&ghost->remote_ghosts->a, index);
  T8_ASSERT (remote_entry->remote_rank == remote);

  /* allocate memory for the send buffer */
  byte_count = data_size * remote_entry->num_elements;
  buffer = *pbuffer = T8_ALLOC (char, byte_count);

  /* We now iterate over the remote trees and their elements to find the
   * local element indices of the remote elements */
  for (itree = 0; itree < (t8_locidx_t) remote_entry->remote_trees.elem_count; itree++) {
    /* tree loop */
    remote_tree = (t8_ghost_remote_tree_t *) t8_sc_array_index_locidx (&remote_entry->remote_trees, itree);
    /* Get the local id of this tree */
    /* TODO: Why does remote_tree store the global id? could be local instead */
    ltreeid = t8_forest_get_local_id (forest, remote_tree->global_id);
    /* Get a pointer to the forest tree */
    local_tree = t8_forest_get_tree (forest, ltreeid);
    elem_count = t8_element_array_get_count (&remote_tree->elements);
    for (ielement = 0; ielement < (t8_locidx_t) elem_count; ielement++) {
      /* element loop */
      /* Get the index of this remote element in its local tree */
      element_pos = *(t8_locidx_t *) t8_sc_array_index_locidx (&remote_tree->element_indices, ielement);
      T8_ASSERT (0 <= element_pos);
      /* Compute the index of this element in the element_data array */
      element_index = local_tree->elements_offset + element_pos;
      /* Copy the data of this element from the element_data array to the send buffer */
      memcpy (buffer + elements_inserted * data_size, sc_array_index (element_data, element_index), data_size);
      elements_inserted++;
    }
  }
  return byte_count;
}

/**
 *  Begin the ghost data exchange communication.
 * 
 *  \param[in] forest         A committed forest.
 *  \param[in] element_data   The element data array.
 * 
 *  \return The ghost data exchange type, as pointer to \see t8_ghost_data_exchange_t.
 */
static t8_ghost_data_exchange_t *
t8_forest_ghost_exchange_begin (t8_forest_t forest, sc_array_t *element_data)
{
  t8_ghost_data_exchange_t *data_exchange;
  t8_forest_ghost_t ghost;
  size_t bytes_to_send, ghost_start;
  int iremote, remote_rank;
  int mpiret, recv_rank, bytes_recv;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  char **send_buffers;
  t8_ghost_process_hash_t lookup_proc, *process_entry, **pfound;
  t8_locidx_t remote_offset, next_offset;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (element_data != NULL);
  T8_ASSERT (forest->ghosts != NULL);

  ghost = forest->ghosts;

  /* Allocate the new exchange context */
  data_exchange = T8_ALLOC (t8_ghost_data_exchange_t, 1);
  /* The number of processes we need to send to */
  data_exchange->num_remotes = ghost->remote_processes->elem_count;
  /* Allocate MPI requests */
  data_exchange->send_requests = T8_ALLOC (sc_MPI_Request, data_exchange->num_remotes);
  data_exchange->recv_requests = T8_ALLOC (sc_MPI_Request, data_exchange->num_remotes);
  /* Allocate pointers to send buffers */
  send_buffers = data_exchange->send_buffers = T8_ALLOC (char *, data_exchange->num_remotes);

  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* Iterate over all remote processes and fill their send buffers */
    remote_rank = *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    /* Fill the send buffers and compute the number of bytes to send */
    bytes_to_send
      = t8_forest_ghost_exchange_fill_send_buffer (forest, remote_rank, send_buffers + iremote, element_data);

    /* Post the asynchronuos send */
    mpiret = sc_MPI_Isend (send_buffers[iremote], bytes_to_send, sc_MPI_BYTE, remote_rank, T8_MPI_GHOST_EXC_FOREST,
                           forest->mpicomm, data_exchange->send_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }

  /* The index in element_data at which the ghost elements start */
  ghost_start = t8_forest_get_local_num_leaf_elements (forest);
  /* Receive the incoming messages */
  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* We need to compute the offset in element_data to which we can receive the message */
    /* Search for this processes' entry in the ghost struct */
    recv_rank = *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    lookup_proc.mpirank = recv_rank;
#if T8_ENABLE_DEBUG
    ret =
#else
    (void)
#endif
      sc_hash_lookup (ghost->process_offsets, &lookup_proc, (void ***) &pfound);
    T8_ASSERT (ret);
    process_entry = *pfound;
    /* In process_entry we stored the offset of this ranks ghosts under all
     * ghosts. Thus in element_data we look at the position
     *  ghost_start + offset
     */
    remote_offset = process_entry->ghost_offset;
    /* Compute the offset of the next remote rank */
    if (iremote + 1 < data_exchange->num_remotes) {
      lookup_proc.mpirank = *(int *) sc_array_index_int (ghost->remote_processes, iremote + 1);
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &lookup_proc, (void ***) &pfound);
      T8_ASSERT (ret);
      process_entry = *pfound;
      next_offset = process_entry->ghost_offset;
    }
    else {
      /* We are the last rank, the next offset is the total number of ghosts */
      next_offset = ghost->num_ghosts_elements;
    }
    /* Calculate the number of bytes to receive */
    bytes_recv = (next_offset - remote_offset) * element_data->elem_size;
    /* receive the message */
    mpiret = sc_MPI_Irecv (sc_array_index (element_data, ghost_start + remote_offset), bytes_recv, sc_MPI_BYTE,
                           recv_rank, T8_MPI_GHOST_EXC_FOREST, forest->mpicomm, data_exchange->recv_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }
  return data_exchange;
}

/**
 *  Free the memory of the ghost data exchange.
 * 
 *  After all associated communication has terminated, this function is used to free the 
 *  allocated memory of the send and receive buffers.
 * 
 *  \param[in] data_exchange  The ghost data exchange type to free memory for.
*/
static void
t8_forest_ghost_exchange_end (t8_ghost_data_exchange_t *data_exchange)
{
  int iproc;

  T8_ASSERT (data_exchange != NULL);
  /* Wait for all communications to end */
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->recv_requests, sc_MPI_STATUSES_IGNORE);
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->send_requests, sc_MPI_STATUSES_IGNORE);

  /* Free the send buffers */
  for (iproc = 0; iproc < data_exchange->num_remotes; iproc++) {
    T8_FREE (data_exchange->send_buffers[iproc]);
  }
  T8_FREE (data_exchange->send_buffers);
  /* free requests */
  T8_FREE (data_exchange->send_requests);
  T8_FREE (data_exchange->recv_requests);
  T8_FREE (data_exchange);
}

void
t8_forest_ghost_exchange_data (t8_forest_t forest, sc_array_t *element_data)
{
  t8_ghost_data_exchange_t *data_exchange;

  t8_debugf ("Entering ghost_exchange_data\n");
  T8_ASSERT (t8_forest_is_committed (forest));

  if (forest->ghosts == nullptr) {
    /* This process has no ghosts */
    return;
  }

  T8_ASSERT (forest->ghosts != NULL);
  T8_ASSERT (element_data != NULL);
  T8_ASSERT ((t8_locidx_t) element_data->elem_count
             == t8_forest_get_local_num_leaf_elements (forest) + t8_forest_get_num_ghosts (forest));

  data_exchange = t8_forest_ghost_exchange_begin (forest, element_data);
  if (forest->profile != nullptr) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime = -sc_MPI_Wtime ();
  }
  t8_forest_ghost_exchange_end (data_exchange);
  if (forest->profile != nullptr) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime += sc_MPI_Wtime ();
  }
  t8_debugf ("Finished ghost_exchange_data\n");
}

/* Print a forest ghost structure */
void
t8_forest_ghost_print (t8_forest_t forest)
{
  t8_forest_ghost_t ghost;
  t8_ghost_remote_t *remote_found;
  t8_ghost_remote_tree_t *remote_tree;
  t8_ghost_process_hash_t proc_hash, **pfound, *found;
  size_t iremote, itree;
#if T8_ENABLE_DEBUG
  int ret;
#endif
  int remote_rank;
  char remote_buffer[BUFSIZ] = "";
  char buffer[BUFSIZ] = "";

  if (forest->ghosts == nullptr) {
    return;
  }
  T8_ASSERT (forest->ghosts != NULL);
  ghost = forest->ghosts;
  snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer), "\tRemotes:\n");
  snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), "\tReceived:\n");

  if (ghost->num_ghosts_elements > 0) {
    for (iremote = 0; iremote < ghost->remote_processes->elem_count; iremote++) {
      /* Get the rank of the remote process */
      remote_rank = *(int *) sc_array_index (ghost->remote_processes, iremote);
      /* Get this remote's entry */
      remote_found = t8_forest_ghost_get_remote (forest, remote_rank);
      /* investigate the entry of this remote process */
      snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer), "\t[Rank %i] (%li trees):\n",
                remote_found->remote_rank, remote_found->remote_trees.elem_count);
      for (itree = 0; itree < remote_found->remote_trees.elem_count; itree++) {
        remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_found->remote_trees, itree);
        snprintf (remote_buffer + strlen (remote_buffer), BUFSIZ - strlen (remote_buffer),
                  "\t\t[id: %lli, class: %s, #elem: %li]\n", (long long) remote_tree->global_id,
                  t8_eclass_to_string[remote_tree->eclass], (long) t8_element_array_get_count (&remote_tree->elements));
      }

      /* Investigate the elements that we received from this process */
      proc_hash.mpirank = remote_rank;
      /* look up this rank in the hash table */
#if T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &proc_hash, (void ***) &pfound);

      T8_ASSERT (ret);
      found = *pfound;
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                "\t[Rank %i] First tree: %li\n\t\t First element: %li\n", remote_rank, (long) found->tree_index,
                (long) found->first_element);
    }
  }
  t8_debugf ("Ghost structure:\n%s\n%s\n", remote_buffer, buffer);
}

/**
 * Completely destroy a ghost structure 
 * 
 * \param[in,out] pghost  The ghost structure to be destroyed.
 */
static void
t8_forest_ghost_reset (t8_forest_ghost_t *pghost)
{
  t8_forest_ghost_t ghost;
  size_t it, it_trees;
  t8_ghost_tree_t *ghost_tree;
  t8_ghost_remote_t *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->rc.refcount == 0);

  /* Clean-up the arrays */
  for (it_trees = 0; it_trees < ghost->ghost_trees->elem_count; it_trees++) {
    ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees, it_trees);
    t8_element_array_reset (&ghost_tree->elements);
  }

  sc_array_destroy (ghost->ghost_trees);
  sc_array_destroy (ghost->remote_processes);
  /* Clean-up the hashtables */
  sc_hash_destroy (ghost->global_tree_to_ghost_tree);
  sc_hash_destroy (ghost->process_offsets);
  /* Clean-up the remote ghost entries */
  for (it = 0; it < ghost->remote_ghosts->a.elem_count; it++) {
    remote_entry = (t8_ghost_remote_t *) sc_array_index (&ghost->remote_ghosts->a, it);
    for (it_trees = 0; it_trees < remote_entry->remote_trees.elem_count; it_trees++) {
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (&remote_entry->remote_trees, it_trees);
      t8_element_array_reset (&remote_tree->elements);
      sc_array_reset (&remote_tree->element_indices);
    }
    sc_array_reset (&remote_entry->remote_trees);
  }
  sc_hash_array_destroy (ghost->remote_ghosts);

  /* Clean-up the memory pools for the data inside
   * the hash tables */
  sc_mempool_destroy (ghost->glo_tree_mempool);
  sc_mempool_destroy (ghost->proc_offset_mempool);

  /* Free the ghost */
  T8_FREE (ghost);
  pghost = nullptr;
}

void
t8_forest_ghost_ref (t8_forest_ghost_t ghost)
{
  T8_ASSERT (ghost != NULL);

  t8_refcount_ref (&ghost->rc);
}

void
t8_forest_ghost_unref (t8_forest_ghost_t *pghost)
{
  t8_forest_ghost_t ghost;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);

  if (t8_refcount_unref (&ghost->rc)) {
    t8_forest_ghost_reset (pghost);
  }
}

void
t8_forest_ghost_destroy (t8_forest_ghost_t *pghost)
{
  T8_ASSERT (pghost != NULL && *pghost != NULL && t8_refcount_is_last (&(*pghost)->rc));
  t8_forest_ghost_unref (pghost);
  T8_ASSERT (*pghost == NULL);
}

T8_EXTERN_C_END ();
