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

#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_element_cxx.hxx>
#include <t8_data/t8_containers.h>
#include <sc_statistics.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* The information for a remote process, what data
 * we have to send to them.
 */
typedef struct
{
  int                 recv_rank;        /* The rank to which we send. */
  size_t              num_bytes;        /* The number of bytes that we send. */
  sc_MPI_Request     *request;  /* Commuication request, not owned by this struct. */
  char               *buffer;   /* The send buffer. */
} t8_ghost_mpi_send_info_t;

/* The information stored for the ghost trees */
typedef struct
{
  t8_gloidx_t         global_id;        /* global id of the tree */
  t8_locidx_t         element_offset;   /* The count of all ghost elements in all smaller ghost trees */
  t8_element_array_t  elements; /* The ghost elements of that tree */
  t8_eclass_t         eclass;   /* The trees element class */
} t8_ghost_tree_t;

/* The data structure stored in the global_tree_to_ghost_tree hash table. */
typedef struct
{
  t8_gloidx_t         global_id;        /* global tree id */
  size_t              index;    /* the index of that global tree in the ghost_trees array. */
} t8_ghost_gtree_hash_t;

/* The data structure stored in the process_offsets array. */
typedef struct
{
  int                 mpirank;  /* rank of the process */
  t8_locidx_t         ghost_offset;     /* The number of ghost elements for all previous ranks */
  size_t              tree_index;       /* index of first ghost tree of this process in ghost_trees */
  size_t              first_element;    /* the index of the first element in the elements array
                                           of the ghost tree. */
} t8_ghost_process_hash_t;

/* The information stored for the remote trees.
 * Each remote process stores an array of these */
typedef struct
{
  t8_gloidx_t         global_id;        /* global id of the tree */
  int                 mpirank;  /* The mpirank of the remote process */
  t8_element_array_t  elements; /* The remote elements of that tree */
  sc_array_t          element_indices;  /* The (tree) local indices of the ghost elements. */
  t8_eclass_t         eclass;   /* The trees element class */
} t8_ghost_remote_tree_t;

typedef struct
{
  int                 remote_rank;      /* The rank of the remote process */
  t8_locidx_t         num_elements;     /* The number of remote elements for this process */
  sc_array_t          remote_trees;     /* Array of the remote trees of this process */
} t8_ghost_remote_t;

#if 0
/* Compare two ghost_tree entries. We need this function to sort the
 * ghost_trees array by global_id. */
static int
t8_ghost_tree_compare (const void *tree_a, const void *tree_b)
{
  const t8_ghost_tree_t *A = (const t8_ghost_tree_t *) tree_a;
  const t8_ghost_tree_t *B = (const t8_ghost_tree_t *) tree_b;

  if (A->global_id < B->global_id) {
    return -1;
  }
  return A->global_id != B->global_id;
}
#endif

/* The hash function for the global tree hash.
 * As hash value we just return the global tree id. */
static unsigned
t8_ghost_gtree_hash_function (const void *ghost_gtree_hash, const void *data)
{
  const t8_ghost_gtree_hash_t *object =
    (const t8_ghost_gtree_hash_t *) ghost_gtree_hash;

  return (unsigned) object->global_id;
}

/* The equal function for two global tree hash objects.
 * Two t8_ghost_gtree_hash_t are considered equal if theit global
 * tree ids are the same.
 */
static int
t8_ghost_gtree_equal_function (const void *ghost_gtreea,
                               const void *ghost_gtreeb, const void *user)
{
  const t8_ghost_gtree_hash_t *objecta =
    (const t8_ghost_gtree_hash_t *) ghost_gtreea;
  const t8_ghost_gtree_hash_t *objectb =
    (const t8_ghost_gtree_hash_t *) ghost_gtreeb;

  /* return true if and only if the global_ids are the same */
  return objecta->global_id == objectb->global_id;
}

/* The hash value for an entry of the process_offsets hash is the
 * processes mpirank. */
static unsigned
t8_ghost_process_hash_function (const void *process_data,
                                const void *user_data)
{
  const t8_ghost_process_hash_t *process =
    (const t8_ghost_process_hash_t *) process_data;

  return process->mpirank;
}

/* The equal function for the process_offsets array.
 * Two entries are the same if their mpiranks are equal. */
static int
t8_ghost_process_equal_function (const void *process_dataa,
                                 const void *process_datab, const void *user)
{
  const t8_ghost_process_hash_t *processa =
    (const t8_ghost_process_hash_t *) process_dataa;
  const t8_ghost_process_hash_t *processb =
    (const t8_ghost_process_hash_t *) process_datab;

  return processa->mpirank == processb->mpirank;
}

/* The hash funtion for the remote_ghosts hash table.
 * The hash value for an mpirank is just the rank */
static unsigned
t8_ghost_remote_hash_function (const void *remote_data, const void *user_data)
{
  const t8_ghost_remote_t *remote = (const t8_ghost_remote_t *) remote_data;

  return remote->remote_rank;
}

/* The equal function for the remote hash table.
 * Two entries are the same if they have the same rank. */
static int
t8_ghost_remote_equal_function (const void *remote_dataa,
                                const void *remote_datab, const void *user)
{
  const t8_ghost_remote_t *remotea = (const t8_ghost_remote_t *) remote_dataa;
  const t8_ghost_remote_t *remoteb = (const t8_ghost_remote_t *) remote_datab;

  return remotea->remote_rank == remoteb->remote_rank;
}

/** This struct is used during a ghost data exchange.
 * Since we use asynchronuous communication, we store the
 * send buffers and mpi requests until we end the communication.
 */
typedef struct
{
  int                 num_remotes;
                    /** The number of processes, we send to */
  char              **send_buffers;
                      /** For each remote the send buffer */
  sc_MPI_Request     *send_requests;
                           /** For each process we send to, the MPI request used */
  sc_MPI_Request     *recv_requests;
                           /** For each process we receive from, the MPI request used */
} t8_ghost_data_exchange_t;

void
t8_forest_ghost_init (t8_forest_ghost_t *pghost, t8_ghost_type_t ghost_type)
{
  t8_forest_ghost_t   ghost;

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
  ghost->global_tree_to_ghost_tree =
    sc_hash_new (t8_ghost_gtree_hash_function, t8_ghost_gtree_equal_function,
                 NULL, NULL);

  /* initialize the process_offset hash table */
  ghost->proc_offset_mempool =
    sc_mempool_new (sizeof (t8_ghost_process_hash_t));
  ghost->process_offsets =
    sc_hash_new (t8_ghost_process_hash_function,
                 t8_ghost_process_equal_function, NULL, NULL);
  /* initialize the remote ghosts hash table */
  ghost->remote_ghosts =
    sc_hash_array_new (sizeof (t8_ghost_remote_t),
                       t8_ghost_remote_hash_function,
                       t8_ghost_remote_equal_function, NULL);
  /* initialize the remote processes array */
  ghost->remote_processes = sc_array_new (sizeof (int));
}

/* Return the remote struct of a given remote rank */
static t8_ghost_remote_t *
t8_forest_ghost_get_remote (t8_forest_t forest, int remote)
{
  t8_ghost_remote_t   remote_search;
#ifdef T8_ENABLE_DEBUG
  int                 ret;
#endif
  size_t              index;

  T8_ASSERT (t8_forest_is_committed (forest));

  remote_search.remote_rank = remote;
#ifdef T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (forest->ghosts->remote_ghosts, &remote_search,
                          &index);
  T8_ASSERT (ret);
  return (t8_ghost_remote_t *)
    sc_array_index (&forest->ghosts->remote_ghosts->a, index);
}

/* Return a remote processes info about the stored ghost elements */
static t8_ghost_process_hash_t *
t8_forest_ghost_get_proc_info (t8_forest_t forest, int remote)
{
  t8_ghost_process_hash_t proc_hash_search, **pproc_hash_found,
    *proc_hash_found;
#ifdef T8_ENABLE_DEBUG
  int                 ret;
#endif

  T8_ASSERT (t8_forest_is_committed (forest));

  proc_hash_search.mpirank = remote;
#ifdef T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_lookup (forest->ghosts->process_offsets, &proc_hash_search,
                    (void ***) &pproc_hash_found);
  T8_ASSERT (ret);
  proc_hash_found = *pproc_hash_found;
  T8_ASSERT (proc_hash_found->mpirank == remote);
  return proc_hash_found;
}

/* return the number of trees in a ghost */
t8_locidx_t
t8_forest_ghost_num_trees (t8_forest_t forest)
{
  if (forest->ghosts == NULL) {
    return 0;
  }
  T8_ASSERT (forest->ghosts != NULL);
  if (forest->ghosts->num_ghosts_elements <= 0) {
    return 0;
  }
  T8_ASSERT (forest->ghosts->ghost_trees != NULL);

  return forest->ghosts->ghost_trees->elem_count;
}

/* Given an index into the ghost_trees array return the ghost tree */
static t8_ghost_tree_t *
t8_forest_ghost_get_tree (t8_forest_t forest, t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t    *ghost_tree;
  t8_forest_ghost_t   ghost;

  T8_ASSERT (t8_forest_is_committed (forest));
  ghost = forest->ghosts;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->ghost_trees != NULL);
  T8_ASSERT (0 <= lghost_tree &&
             lghost_tree < t8_forest_ghost_num_trees (forest));

  ghost_tree =
    (t8_ghost_tree_t *) t8_sc_array_index_locidx (ghost->ghost_trees,
                                                  lghost_tree);
  return ghost_tree;
}

t8_locidx_t
t8_forest_ghost_get_tree_element_offset (t8_forest_t
                                         forest, t8_locidx_t lghost_tree)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return t8_forest_ghost_get_tree (forest, lghost_tree)->element_offset;
}

/* Given an index in the ghost_tree array, return this tree's number of elements */
t8_locidx_t
t8_forest_ghost_tree_num_elements (t8_forest_t forest,
                                   t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t    *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return t8_element_array_get_count (&ghost_tree->elements);
}

t8_element_array_t *
t8_forest_ghost_get_tree_elements (t8_forest_t forest,
                                   t8_locidx_t lghost_tree)
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
  if (sc_hash_lookup (forest->ghosts->global_tree_to_ghost_tree, &query,
                      (void ***) &pfound)) {
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
t8_forest_ghost_get_tree_class (t8_forest_t forest, t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t    *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->eclass;
}

/* Given an index in the ghost_tree array, return this tree's global id */
t8_gloidx_t
t8_forest_ghost_get_global_treeid (t8_forest_t forest,
                                   t8_locidx_t lghost_tree)
{
  t8_ghost_tree_t    *ghost_tree;
  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  return ghost_tree->global_id;
}

/* Given an index into the ghost_trees array and for that tree an element index,
 * return the corresponding element. */
t8_element_t       *
t8_forest_ghost_get_element (t8_forest_t forest, t8_locidx_t lghost_tree,
                             t8_locidx_t lelement)
{
  t8_ghost_tree_t    *ghost_tree;

  T8_ASSERT (t8_forest_is_committed (forest));

  ghost_tree = t8_forest_ghost_get_tree (forest, lghost_tree);
  T8_ASSERT (0 <= lelement &&
             lelement < t8_forest_ghost_tree_num_elements (forest,
                                                           lghost_tree));
  return t8_element_array_index_locidx (&ghost_tree->elements, lelement);
}

/* Initialize a t8_ghost_remote_tree_t */
static void
t8_ghost_init_remote_tree (t8_forest_t forest, t8_gloidx_t gtreeid,
                           int remote_rank,
                           t8_eclass_t eclass,
                           t8_ghost_remote_tree_t *remote_tree)
{
  t8_eclass_scheme_c *ts;
  t8_locidx_t         local_treeid;

  T8_ASSERT (remote_tree != NULL);

  ts = t8_forest_get_eclass_scheme (forest, eclass);
  local_treeid = gtreeid - t8_forest_get_first_local_tree_id (forest);
  /* Set the entries of the new remote tree */
  remote_tree->global_id = gtreeid;
  remote_tree->mpirank = remote_rank;
  remote_tree->eclass = t8_forest_get_eclass (forest, local_treeid);
  /* Initialize the array to store the element */
  t8_element_array_init (&remote_tree->elements, ts);
  /* Initialize the array to store the element indices. */
  sc_array_init (&remote_tree->element_indices, sizeof (t8_locidx_t));
}

/* Add a new element to the remote hash table (if not already in it).
 * Must be called for elements in linear order
 * element_index is the tree local index of this element */
static void
t8_ghost_add_remote (t8_forest_t forest, t8_forest_ghost_t ghost,
                     int remote_rank, t8_locidx_t ltreeid,
                     const t8_element_t *elem, t8_locidx_t element_index)
{
  t8_ghost_remote_t   remote_entry_lookup, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_element_t       *elem_copy;
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass;
  sc_array_t         *remote_array;
  size_t              index, element_count;
  t8_gloidx_t         gtreeid;
  int                *remote_process_entry;
  int                 level, copy_level = 0;

  /* Get the tree's element class and the scheme */
  eclass = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  gtreeid = t8_forest_get_first_local_tree_id (forest) + ltreeid;

  /* Check whether the remote_rank is already present in the remote ghosts
   * array. */
  remote_entry_lookup.remote_rank = remote_rank;
  remote_entry = (t8_ghost_remote_t *)
    sc_hash_array_insert_unique (ghost->remote_ghosts,
                                 (void *) &remote_entry_lookup, &index);
  if (remote_entry != NULL) {
    /* The remote rank was not in the array and was inserted now */
    remote_entry->remote_rank = remote_rank;
    remote_entry->num_elements = 0;
    /* Initialize the tree array of the new entry */
    sc_array_init_size (&remote_entry->remote_trees,
                        sizeof (t8_ghost_remote_tree_t), 1);
    /* Get a pointer to the new entry */
    remote_tree = (t8_ghost_remote_tree_t *)
      sc_array_index (&remote_entry->remote_trees, 0);
    /* initialize the remote_tree */
    t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, eclass,
                               remote_tree);
    /* Since the rank is a new remote rank, we also add it to the
     * remote ranks array */
    remote_process_entry = (int *) sc_array_push (ghost->remote_processes);
    *remote_process_entry = remote_rank;
  }
  else {
    /* The remote rank alrady is contained in the remotes array at
     * position index. */
    remote_array = &ghost->remote_ghosts->a;
    remote_entry = (t8_ghost_remote_t *) sc_array_index (remote_array, index);
    T8_ASSERT (remote_entry->remote_rank == remote_rank);
    /* Check whether the tree has already an entry for this process.
     * Since we only add in local tree order the current tree is either
     * the last entry or does not have an entry yet. */
    remote_tree = (t8_ghost_remote_tree_t *)
      sc_array_index (&remote_entry->remote_trees,
                      remote_entry->remote_trees.elem_count - 1);
    if (remote_tree->global_id != gtreeid) {
      /* The tree does not exist in the array. We thus need to add it and
       * initialize it. */
      remote_tree = (t8_ghost_remote_tree_t *)
        sc_array_push (&remote_entry->remote_trees);
      t8_ghost_init_remote_tree (forest, gtreeid, remote_rank, eclass,
                                 remote_tree);
    }
  }
  /* remote_tree now points to a valid entry for the tree.
   * We can add a copy of the element to the elements array
   * if it does not exist already. If it exists it is the last entry in
   * the array. */
#ifdef T8_ENABLE_DEBUG
  {
    /* debugging assertion that the element is really not contained already */
    int                 ielem;
    t8_element_t       *test_el;
    int                 elem_count =
      t8_element_array_get_count (&remote_tree->elements);
    for (ielem = 0; ielem < elem_count - 1; ielem++) {
      test_el = t8_element_array_index_int (&remote_tree->elements, ielem);
      SC_CHECK_ABORTF (ts->t8_element_compare (test_el, elem),
                       "Local element %i already in remote ghosts at pos %i\n",
                       element_index, ielem);
    }
  }
#endif
  elem_copy = NULL;
  level = ts->t8_element_level (elem);
  element_count = t8_element_array_get_count (&remote_tree->elements);
  if (element_count > 0) {
    elem_copy =
      t8_element_array_index_locidx (&remote_tree->elements,
                                     element_count - 1);
    copy_level = ts->t8_element_level (elem_copy);
  }
  /* Check if the element was not contained in the array.
   * If so, we add a copy of elem to the array.
   * Otherwise, we do nothing. */
  if (elem_copy == NULL ||
      level != copy_level ||
      ts->t8_element_get_linear_id (elem_copy, copy_level) !=
      ts->t8_element_get_linear_id (elem, level)) {
    /* Add the element */
    elem_copy = t8_element_array_push (&remote_tree->elements);
    ts->t8_element_copy (elem, elem_copy);
    /* Add the index of the element */
    *(t8_locidx_t *) sc_array_push (&remote_tree->element_indices) =
      element_index;
    remote_entry->num_elements++;
  }
}

#if 0
/* In ghost version 3, the remote elements are not added in their linear order
 * to the ghost struct, and same elements may be added more than once.
 * Since we need to call t8_ghost_add_remote in linear order and only once per element,
 * we temporally store the element indices in a hash_array,
 * after all remote elements of a tree are parsed, we sort the hash_array and
 * call t8_ghost_add_remote for each entry. */

typedef struct
{
  t8_locidx_t         element_index;    /* The tree local index of this element */
  sc_array_t          remote_ranks;     /* All ranks that this element is a remote of */
} t8_forest_ghost_rem_el_index_t;

/* hash function usese element_index as hash */
static unsigned
t8_forest_ghost_rem_el_index_hash (const void *index, const void *u)
{
  return ((t8_forest_ghost_rem_el_index_t *) index)->element_index;
}

/* comparison function for sorting uses element_index for equality */
static int
t8_forest_ghost_rem_el_index_compare (const void *indexa, const void *indexb)
{
  t8_forest_ghost_rem_el_index_t *Ia =
    (t8_forest_ghost_rem_el_index_t *) indexa;
  t8_forest_ghost_rem_el_index_t *Ib =
    (t8_forest_ghost_rem_el_index_t *) indexb;

  return Ia->element_index - Ib->element_index;
}

/* equal function uses element_index for equality */
static int
t8_forest_ghost_rem_el_index_eq (const void *indexa, const void *indexb,
                                 const void *u)
{
  return !t8_forest_ghost_rem_el_index_compare (indexa, indexb);
}

/* Sort all added remote indices, parse them and add the elements as ghosts.
 * This function also truncates the hash_array. */
static void
t8_forest_ghost_add_remote_indices (t8_forest_t forest,
                                    t8_forest_ghost_t ghost,
                                    t8_locidx_t ltreeid,
                                    sc_hash_array_t * rem_el_indices)
{
  sc_array_t         *indices_view;
  size_t              iremote, ielem;
  t8_forest_ghost_rem_el_index_t *index;
  t8_element_t       *element;
  int                 remote_rank;

  /* Get a pointer to the underlying array */
  indices_view = &rem_el_indices->a;
  /* sort the array */
  sc_array_sort (indices_view, t8_forest_ghost_rem_el_index_compare);

  /* iterate through the element and add as remotes */
  for (ielem = 0; ielem < indices_view->elem_count; ielem++) {
    /* Get the element's index */
    index =
      (t8_forest_ghost_rem_el_index_t *) sc_array_index (indices_view, ielem);
    /* Get a pointer to the element */
    element =
      t8_forest_get_element_in_tree (forest, ltreeid, index->element_index);
    /* parse all ranks that this element is a remote of and add the elemetn as
     * remote ghost */
    for (iremote = 0; iremote < index->remote_ranks.elem_count; iremote++) {
      remote_rank = *(int *) sc_array_index (&index->remote_ranks, iremote);
      t8_ghost_add_remote (forest, ghost, remote_rank, ltreeid, element,
                           index->element_index);
    }
    /* Clean-up the memory for the remote ranks */
    sc_array_reset (&index->remote_ranks);
  }
  /* Clean the hash_array */
  sc_hash_array_truncate (rem_el_indices);
}

/* Add an entry to the hash_array of element indices.
 * If this element was already considered for the remote then the element
 * index is not added */
static void
t8_forest_ghost_add_remote_index (sc_hash_array_t * rem_el_indices,
                                  t8_locidx_t element_index, int remote_rank)
{
  t8_forest_ghost_rem_el_index_t index_search, *index_found;
  size_t              position, iremote;
  int                 check_rank;

  index_search.element_index = element_index;
  /* Try to insert this entry. If this element index already has an entry,
   * then position is set to the array position of the contained entry and
   * index_found is NULL.
   * otherwise to the position of the newly inserted entry, and
   * index_found points to this entry, */
  index_found = (t8_forest_ghost_rem_el_index_t *)
    sc_hash_array_insert_unique (rem_el_indices, &index_search, &position);
  if (index_found != NULL) {
    /* This is a new index and we need to initialize it first */
    index_found->element_index = element_index;
    sc_array_init (&index_found->remote_ranks, sizeof (int));
    /* we add the remote_rank to this entry */
    *(int *) sc_array_push (&index_found->remote_ranks) = remote_rank;
    return;
  }
  else {
    /* The entry was already contained in the hash_array */
    /* Get a pointer to the found entry */
    index_found =
      (t8_forest_ghost_rem_el_index_t *) sc_array_index (&rem_el_indices->a,
                                                         position);
    T8_ASSERT (index_found->element_index == element_index);
    /* Search whether the remote_rank was already added as a remote rank of
     * this element */
    /* TODO: The number of remote ranks of an element is, in theory, not bounded,
     *       thus this search can be very expensive. Is this a problem in praxis? */
    for (iremote = 0; iremote < index_found->remote_ranks.elem_count;
         iremote++) {
      check_rank =
        *(int *) sc_array_index (&index_found->remote_ranks, iremote);
      if (check_rank == remote_rank) {
        /* The rank was already stored as a remote for this element,
         * we can abort here */
        return;
      }
    }
    /* The remote_rank was not added as a rank for this element, we add it */
    *(int *) sc_array_push (&index_found->remote_ranks) = remote_rank;
  }
}

typedef int8_t      t8_element_face_flag_t;

typedef struct
{
  sc_hash_array_t    *rem_el_indices;   /* The indices of the so far added remote elements of this tree */
  int                 face_owner_low, face_owner_high;  /* The lowest and highest owner
                                                           at the neighbor face for the parent element */
  int                 neighbor_unique_owner;    /* If non-negative, then the owner of the neighbor face
                                                   is unique and this stores its rank */
  t8_eclass_t         eclass;
} t8_forest_ghost_iterate_face_data_t;

/* As soon as the ghost top-down search finds an element whose descendants at a face f
 * are all owned by the current rank, we start a face iteration for this element and face.
 * Thus, all leafs of the element that lie on the face are traversed.
 * This happens in a top-down search, and on each intermediate level this callback
 * function is called. If it returns true, the top-down iteration continues.
 * In this callback, we compute the owners at the neighbor face and if
 *  a) all neighbor leafs are owned by this rank,
 *     then we do not need to consider the leafs as remote elements at face f,
 *     thus we return false and the top-down search stops.
 *  b) the neighbor leafs are owned by ranks different to the current rank.
 *     In this case, we continue the search.
 *  c) the element is a leaf,
 *     then we compute its owners and add it as remote to all owners that are not the current rank.
 */
static int
t8_forest_ghost_iterate_face_add_remote (t8_forest_t forest,
                                         t8_locidx_t ltreeid,
                                         const t8_element_t *element,
                                         int face, void *user_data,
                                         t8_locidx_t leaf_index)
{
  t8_forest_ghost_iterate_face_data_t *data;
  sc_array_t          owners_at_face;   /* TODO: we could also at an sc_array to data and reuse it everytime */
  int                 lower, upper;
  int                 remote_rank;
  size_t              iown;

  data = (t8_forest_ghost_iterate_face_data_t *) user_data;
  lower = data->face_owner_low;
  upper = data->face_owner_high;

  if (leaf_index >= 0) {
    /* The element is a leaf, we compute its neighbor owners and add it
     * as a remote to the ghost struct */
    if (data->neighbor_unique_owner < 0) {
      /* There may be more than one owner processes of the neighbor face,
       * since the owners of this leaf have not been computed, we compute them
       * now. */
      sc_array_init_size (&owners_at_face, sizeof (int), 2);
      /* Set the lower and upper bound of the face owners */
      *(int *) sc_array_index (&owners_at_face, 0) = lower;
      *(int *) sc_array_index (&owners_at_face, 1) = upper;
      t8_forest_element_owners_at_neigh_face (forest, ltreeid, element,
                                              face, &owners_at_face);
      /* parse through all owners of the neighbor element and add this leaf
       * as a remote for each */
      for (iown = 0; iown < owners_at_face.elem_count; iown++) {
        remote_rank = *(int *) sc_array_index (&owners_at_face, iown);
        /* TODO: Does this still work, even though the remotes are not added in
         *       order? */
        if (remote_rank != forest->mpirank) {
          t8_forest_ghost_add_remote_index (data->rem_el_indices, leaf_index,
                                            remote_rank);
        }
      }
      sc_array_reset (&owners_at_face);
    }
    else {
      /* The owner is unique, add the leaf as remote to this owner */
      t8_forest_ghost_add_remote_index (data->rem_el_indices, leaf_index,
                                        data->neighbor_unique_owner);
    }
    return 0;                   /* return value is ignored for leafs */
  }

  if (lower > upper) {
    /* This face does not have any neighbors (domain boundary) */
    /* Do not continue recursion */
    SC_ABORT_NOT_REACHED ();    /* TODO: remove if this never hits. We just keep it in as a correctness check */
    return 0;
  }

  /* we now compute the bounds for owners at the neighbor face for this element */
  t8_forest_element_owners_at_neigh_face_bounds (forest, ltreeid, element,
                                                 face, &lower, &upper);

#if 0
  /* TODO: We used this to store the lower and upper bounds for the next
   * level, however if the top-down search enters a new sibling, this information
   * is not correct any longer. Now, we do not save new lower and upper bounds at
   * all. We could store for each level the bounds, memory O (maxlevel) */
  /* Set the new lower and upper bound for the owners */
  data->face_owner_low = *(int *) sc_array_index (&owners_at_face, 0);
  if (owners_at_face.elem_count >= 2) {
    data->face_owner_high = *(int *)
      sc_array_index (&owners_at_face, owners_at_face.elem_count - 1);
  }
  else {
    data->face_owner_high = data->face_owner_low;
  }
#endif
  if (lower > upper) {
    /* there is no neighbor */
    return 0;
  }
  if (lower == upper) {
    /* There is only one owner of the neighbor element at the face.
     * We check if it is the current rank, if so, there is nothing left to do. */
    if (lower == forest->mpirank) {
      /* do not continue recursion */
      return 0;
    }
    /* There is exactly one owner across the face and it is not the current rank */
    data->neighbor_unique_owner = lower;
    return 1;
  }
  else {
    /* We cannot say whether the owner of the neighbor face is unique and
     * different from the current rank */
    data->neighbor_unique_owner = -1;
  }
  return 1;
}
#endif

typedef struct
{
  sc_array_t          bounds_per_level; /* For each level from the nca to the parent of the current element
                                           we store for each face the lower and upper bounds of the owners at
                                           this face. We also store bounds for the element's owners.
                                           Each entry is an array of 2 * (max_num_faces + 1) integers,
                                           | face_0 low | face_0 high | ... | face_n low | face_n high | owner low | owner high | */
  sc_array_t          face_owners;      /* Temporary storage for all owners at a leaf's face */
  t8_eclass_scheme_c *ts;
  t8_gloidx_t         gtreeid;
  int                 level_nca;        /* The refinement level of the root element in the search.
                                           At position element_level - level_nca in bounds_per_level are the bounds
                                           for the parent of element. */
  int                 max_num_faces;
  t8_eclass_t         eclass;
#ifdef T8_ENABLE_DEBUG
  t8_locidx_t         left_out; /* Count the elements for which we skip the search */
#endif
} t8_forest_ghost_boundary_data_t;

static int
t8_forest_ghost_search_boundary (t8_forest_t forest, t8_locidx_t ltreeid,
                                 const t8_element_t *element,
                                 const int is_leaf,
                                 t8_element_array_t *leafs,
                                 t8_locidx_t tree_leaf_index, void *query,
                                 size_t query_index)
{
  t8_forest_ghost_boundary_data_t *data =
    (t8_forest_ghost_boundary_data_t *) t8_forest_get_user_data (forest);
  int                 num_faces, iface, faces_totally_owned, level;
  int                 parent_face;
  int                 lower, upper, *bounds, *new_bounds, parent_lower,
    parent_upper;
  int                 el_lower, el_upper;
  int                 element_is_owned, iproc, remote_rank;

  /* First part: the search enters a new tree, we need to reset the user_data */
  if (t8_forest_global_tree_id (forest, ltreeid) != data->gtreeid) {
    int                 max_num_faces;
    /* The search has entered a new tree, store its eclass and element scheme */
    data->gtreeid = t8_forest_global_tree_id (forest, ltreeid);
    data->eclass = t8_forest_get_eclass (forest, ltreeid);
    data->ts = t8_forest_get_eclass_scheme (forest, data->eclass);
    data->level_nca = data->ts->t8_element_level (element);
    data->max_num_faces = data->ts->t8_element_max_num_faces (element);
    max_num_faces = data->max_num_faces;
    sc_array_reset (&data->bounds_per_level);
    sc_array_init_size (&data->bounds_per_level,
                        2 * (max_num_faces + 1) * sizeof (int), 1);
    /* Set the (imaginary) owner bounds for the parent of the root element */
    bounds = (int *) sc_array_index (&data->bounds_per_level, 0);
    for (iface = 0; iface < max_num_faces + 1; iface++) {
      bounds[iface * 2] = 0;
      bounds[iface * 2 + 1] = forest->mpisize - 1;
    }
    /* TODO: compute bounds */
  }

  /* The level of the current element */
  level = data->ts->t8_element_level (element);
  /* Get a pointer to the owner at face bounds of this element, if there doesnt exist
   * an entry for this in the bounds_per_level array yet, we allocate it */
  T8_ASSERT (level >= data->level_nca);
  if (data->bounds_per_level.elem_count <=
      (size_t) level - data->level_nca + 1) {
    T8_ASSERT (data->bounds_per_level.elem_count ==
               (size_t) level - data->level_nca + 1);
    new_bounds = (int *) sc_array_push (&data->bounds_per_level);
  }
  else {
    new_bounds = (int *) sc_array_index (&data->bounds_per_level,
                                         level - data->level_nca + 1);
  }

  /* Get a pointer to the owner bounds of the parent */
  bounds =
    (int *) sc_array_index (&data->bounds_per_level, level - data->level_nca);
  /* Get bounds for the element's parent's owners */
  parent_lower = bounds[2 * data->max_num_faces];
  parent_upper = bounds[2 * data->max_num_faces + 1];
  /* Temporarily store them to serve as bounds for this element's owners */
  el_lower = parent_lower;
  el_upper = parent_upper;
  /* Compute bounds for the element's owners */
  t8_forest_element_owners_bounds (forest, data->gtreeid, element,
                                   data->eclass, &el_lower, &el_upper);
  /* Set these as the new bounds */
  new_bounds[2 * data->max_num_faces] = el_lower;
  new_bounds[2 * data->max_num_faces + 1] = el_upper;
  element_is_owned = (el_lower == el_upper);
  num_faces = data->ts->t8_element_num_faces (element);
  faces_totally_owned = 1;

  /* TODO: we may not carry on with the face computations if the element is not
   *       totally onwed and immediately return 1. However, how do we set the bounds for
   *       the face owners then?
   */
  for (iface = 0; iface < num_faces; iface++) {
    /* Compute the face number of the parent to reuse the bounds */
    parent_face = data->ts->t8_element_face_parent_face (element, iface);
    if (parent_face >= 0) {
      /* This face was also a face of the parent, we reuse the computed bounds */
      lower = bounds[parent_face * 2];
      upper = bounds[parent_face * 2 + 1];
    }
    else {
      /* this is an inner face, thus the face owners must be owners of the
       * parent element */
      lower = parent_lower;
      upper = parent_upper;
    }

    if (!is_leaf) {
      /* The element is not a leaf, we compute bounds for the face neighbor owners,
       * if all face neighbors are owned by this rank, and the element is completely
       * owned, then we do not continue the search. */
      /* Compute the owners of the neighbor at this face of the element */
      t8_forest_element_owners_at_neigh_face_bounds (forest, ltreeid, element,
                                                     iface, &lower, &upper);
      /* Store the new bounds at the entry for this element */
      new_bounds[iface * 2] = lower;
      new_bounds[iface * 2 + 1] = upper;
      if (lower == upper && lower == forest->mpirank) {
        /* All neighbor leafs at this face are owned by the current rank */
        faces_totally_owned = faces_totally_owned && 1;
      }
      else {
        faces_totally_owned = 0;
      }
    }
    else {
      /* The element is a leaf, we compute all of its face neighbor owners
       * and add the element as a remote element to all of them. */
      sc_array_resize (&data->face_owners, 2);
      /* The first and second entry in the face_owners array serve as lower
       * and upper bound */
      *(int *) sc_array_index (&data->face_owners, 0) = lower;
      *(int *) sc_array_index (&data->face_owners, 1) = upper;
      t8_forest_element_owners_at_neigh_face (forest, ltreeid, element,
                                              iface, &data->face_owners);
      /*TODO: add as remotes */
      for (iproc = 0; iproc < (int) data->face_owners.elem_count; iproc++) {
        remote_rank = *(int *) sc_array_index (&data->face_owners, iproc);
        if (remote_rank != forest->mpirank) {
          t8_ghost_add_remote (forest, forest->ghosts, remote_rank, ltreeid,
                               element, tree_leaf_index);
        }
      }
    }
  }                             /* end face loop */
#if 0
  /* TODO: can we remove this code? */
  if (element_is_owned || face_totally_owned) {
    /* Either all descendants of element are owned by the current rank
     * or all of its leafs at the face are. */
    face_it_data.face_owner_low = 0;
    face_it_data.face_owner_high = forest->mpisize - 1;
    face_it_data.eclass = data->eclass;
    face_it_data.neighbor_unique_owner = -1;
    face_it_data.rem_el_indices = data->rem_el_indices;
    t8_forest_iterate_faces (forest, ltreeid, element, iface, leafs,
                             &face_it_data, correct_tree_leaf_index,
                             t8_forest_ghost_iterate_face_add_remote);
  }
#endif
  if (faces_totally_owned && element_is_owned) {
    /* The element only has local descendants and all of its face neighbors
     * are local as well. We do not continue the search */
#ifdef T8_ENABLE_DEBUG
    if (tree_leaf_index < 0) {
      data->left_out += t8_element_array_get_count (leafs);
    }
#endif
    return 0;
  }
  /* Continue the top-down search if this element or its face neighbors are
   * not completely owned by the rank. */
  return 1;
}

/* Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 */
static void
t8_forest_ghost_fill_remote_v3 (t8_forest_t forest)
{
  t8_forest_ghost_boundary_data_t data;
  void               *store_user_data = NULL;

  /* Start with invalid entries in the user data.
   * These are set in t8_forest_ghost_search_boundary each time
   * a new tree is entered */
  data.eclass = T8_ECLASS_COUNT;
  data.gtreeid = -1;
  data.ts = NULL;
#ifdef T8_ENABLE_DEBUG
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
  t8_forest_search (forest, t8_forest_ghost_search_boundary, NULL, NULL);

  /* Reset the user data from before search */
  t8_forest_set_user_data (forest, store_user_data);

  /* Reset the data arrays */
  sc_array_reset (&data.face_owners);
  sc_array_reset (&data.bounds_per_level);
#ifdef T8_ENABLE_DEBUG
#endif
}

/* Fill the remote ghosts of a ghost structure.
 * We iterate through all elements and check if their neighbors
 * lie on remote processes. If so, we add the element to the
 * remote_ghosts array of ghost.
 * We also fill the remote_processes here.
 * If ghost_method is 0, then we assume a balanced forest and
 * construct the remote processes by looking at the half neighbors of an element.
 * Otherwise, we use the owners_at_face method.
 */
static void
t8_forest_ghost_fill_remote (t8_forest_t forest, t8_forest_ghost_t ghost,
                             int ghost_method)
{
  t8_element_t       *elem, **half_neighbors = NULL;
  t8_locidx_t         num_local_trees, num_tree_elems;
  t8_locidx_t         itree, ielem;
  t8_tree_t           tree;
  t8_eclass_t         tree_class, neigh_class, last_class;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_scheme_c *ts, *neigh_scheme = NULL, *prev_neigh_scheme = NULL;

  int                 iface, num_faces;
  int                 num_face_children, max_num_face_children = 0;
  int                 ichild, owner;
  sc_array_t          owners, tree_owners;
  int                 is_atom;

  last_class = T8_ECLASS_COUNT;
  num_local_trees = t8_forest_get_num_local_trees (forest);
  if (ghost_method != 0) {
    sc_array_init (&owners, sizeof (int));
    sc_array_init (&tree_owners, sizeof (int));
  }

  /* Loop over the trees of the forest */
  for (itree = 0; itree < num_local_trees; itree++) {
    /* Get a pointer to the tree, the class of the tree, the
     * scheme associated to the class and the number of elements in
     * this tree. */
    tree = t8_forest_get_tree (forest, itree);
    tree_class = t8_forest_get_tree_class (forest, itree);
    ts = t8_forest_get_eclass_scheme (forest, tree_class);

    /* Loop over the elements of this tree */
    num_tree_elems = t8_forest_get_tree_element_count (tree);
    for (ielem = 0; ielem < num_tree_elems; ielem++) {
      /* Get the element of the tree */
      elem = t8_forest_get_tree_element (tree, ielem);
      num_faces = ts->t8_element_num_faces (elem);
      if (ts->t8_element_level (elem) == ts->t8_element_maxlevel ()) {
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
        neigh_class =
          t8_forest_element_neighbor_eclass (forest, itree, elem, iface);
        neigh_scheme = t8_forest_get_eclass_scheme (forest, neigh_class);
        if (ghost_method == 0) {
          /* Use half neighbors */
          /* Get the number of face children of the element at this face */
          num_face_children = ts->t8_element_num_face_children (elem, iface);
          /* regrow the half_neighbors array if neccessary.
           * We also need to reallocate it, if the element class of the neighbor
           * changes */
          if (max_num_face_children < num_face_children ||
              last_class != neigh_class) {
            if (max_num_face_children > 0) {
              /* Clean-up memory */
              prev_neigh_scheme->t8_element_destroy (max_num_face_children,
                                                     half_neighbors);
              T8_FREE (half_neighbors);
            }
            half_neighbors = T8_ALLOC (t8_element_t *, num_face_children);
            /* Allocate memory for the half size face neighbors */
            neigh_scheme->t8_element_new (num_face_children, half_neighbors);
            max_num_face_children = num_face_children;
            last_class = neigh_class;
            prev_neigh_scheme = neigh_scheme;
          }
          if (!is_atom) {
            /* Construct each half size neighbor */
            neighbor_tree =
              t8_forest_element_half_face_neighbors (forest, itree, elem,
                                                     half_neighbors,
                                                     neigh_scheme, iface,
                                                     num_face_children, NULL);
          }
          else {
            int                 dummy_neigh_face;
            /* This element has maximum level, we only construct its neighbor */
            neighbor_tree =
              t8_forest_element_face_neighbor (forest, itree, elem,
                                               half_neighbors[0],
                                               neigh_scheme, iface,
                                               &dummy_neigh_face);
          }
          if (neighbor_tree >= 0) {
            /* If there exist face neighbor elements (we are not at a domain boundary */
            /* Find the owner process of each face_child */
            for (ichild = 0; ichild < num_face_children; ichild++) {
              /* find the owner */
              owner =
                t8_forest_element_find_owner (forest, neighbor_tree,
                                              half_neighbors[ichild],
                                              neigh_class);
              T8_ASSERT (0 <= owner && owner < forest->mpisize);
              if (owner != forest->mpirank) {
                /* Add the element as a remote element */
                t8_ghost_add_remote (forest, ghost, owner, itree, elem,
                                     ielem);
              }
            }
          }
        }                       /* end ghost_method 0 */
        else {
          size_t              iowner;
          /* Construc the owners at the face of the neighbor element */
          t8_forest_element_owners_at_neigh_face (forest, itree, elem, iface,
                                                  &owners);
          T8_ASSERT (owners.elem_count >= 0);
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
      }                         /* end face loop */
    }                           /* end element loop */
  }                             /* end tree loop */

  if (forest->profile != NULL) {
    /* If profiling is enabled, we count the number of remote processes. */
    forest->profile->ghosts_remotes = ghost->remote_processes->elem_count;
  }
  /* Clean-up memory */
  if (ghost_method == 0) {
    if (half_neighbors != NULL) {
      neigh_scheme->t8_element_destroy (max_num_face_children,
                                        half_neighbors);
      T8_FREE (half_neighbors);
    }
  }
  else {
    sc_array_reset (&owners);
    sc_array_reset (&tree_owners);
  }
}

/* Begin sending the ghost elements from the remote ranks
 * using non-blocking communication.
 * Afterward
 *  t8_forest_ghost_send_end
 * must be called to end the communication.
 * Returns an array of mpi_send_info_t, one for each remote rank.
 */
static t8_ghost_mpi_send_info_t *
t8_forest_ghost_send_start (t8_forest_t forest, t8_forest_ghost_t ghost,
                            sc_MPI_Request ** requests)
{
  int                 proc_index, remote_rank;
  int                 num_remotes;
  size_t              remote_index;
  t8_ghost_remote_t  *remote_entry;
  sc_array_t         *remote_trees;
  t8_ghost_remote_tree_t *remote_tree = NULL;
  t8_ghost_mpi_send_info_t *send_info, *current_send_info;
  char               *current_buffer;
  size_t              bytes_written, element_bytes, element_count,
    element_size;
#ifdef T8_ENABLE_DEBUG
  size_t              acc_el_count = 0;
#endif
  int                 mpiret;

  /* Allocate a send_buffer for each remote rank */
  num_remotes = ghost->remote_processes->elem_count;
  send_info = T8_ALLOC (t8_ghost_mpi_send_info_t, num_remotes);
  *requests = T8_ALLOC (sc_MPI_Request, num_remotes);

  /* Loop over all remote processes */
  for (proc_index = 0; proc_index < (int) ghost->remote_processes->elem_count;
       proc_index++) {
    current_send_info = send_info + proc_index;
    /* Get the rank of the current remote process. */
    remote_rank = *(int *) sc_array_index_int (ghost->remote_processes,
                                               proc_index);
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
    current_send_info->num_bytes +=
      T8_ADD_PADDING (current_send_info->num_bytes);
    /* TODO: put this in a funtion */
    /* TODO: Use remote_entry to count the number of bytes while inserting
     *        the remote ghosts. */
    remote_trees = &remote_entry->remote_trees;
    for (remote_index = 0; remote_index < remote_trees->elem_count;
         remote_index++) {
      /* Get the next remote tree. */
      remote_tree = (t8_ghost_remote_tree_t *) sc_array_index (remote_trees,
                                                               remote_index);
      /* We will store the global tree id, the element class and the list
       * of elements in the send_buffer. */
      current_send_info->num_bytes += sizeof (t8_gloidx_t);
      /* add padding before the eclass */
      current_send_info->num_bytes +=
        T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += sizeof (t8_eclass_t);
      /* add padding before the elements */
      current_send_info->num_bytes +=
        T8_ADD_PADDING (current_send_info->num_bytes);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_count = t8_element_array_get_count (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* We will store the number of elements */
      current_send_info->num_bytes += sizeof (size_t);
      /* add padding before the elements */
      current_send_info->num_bytes +=
        T8_ADD_PADDING (current_send_info->num_bytes);
      current_send_info->num_bytes += element_bytes;
      /* add padding after the elements */
      current_send_info->num_bytes +=
        T8_ADD_PADDING (current_send_info->num_bytes);
    }

    /* We now now the number of bytes for our send_buffer and thus
     * allocate it. */
    current_send_info->buffer = T8_ALLOC_ZERO (char,
                                               current_send_info->num_bytes);

    /* We iterate through the tree again and store the tree info and the elements
     * into the send_buffer. */
    current_buffer = current_send_info->buffer;
    bytes_written = 0;
    /* Start with the number of remote trees in the buffer */
    memcpy (current_buffer + bytes_written, &remote_trees->elem_count,
            sizeof (size_t));
    bytes_written += sizeof (size_t);
    bytes_written += T8_ADD_PADDING (bytes_written);
#ifdef T8_ENABLE_DEBUG
    acc_el_count = 0;
#endif
    for (remote_index = 0; remote_index < remote_trees->elem_count;
         remote_index++) {
      /* Get a pointer to the tree */
      remote_tree =
        (t8_ghost_remote_tree_t *) sc_array_index (remote_trees,
                                                   remote_index);
      T8_ASSERT (remote_tree->mpirank == remote_rank);

      /* Copy the global tree id */
      memcpy (current_buffer + bytes_written, &remote_tree->global_id,
              sizeof (t8_gloidx_t));
      bytes_written += sizeof (t8_gloidx_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Copy the trees element class */
      memcpy (current_buffer + bytes_written, &remote_tree->eclass,
              sizeof (t8_eclass_t));
      bytes_written += sizeof (t8_eclass_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* Store the number of elements in the buffer */
      element_count = t8_element_array_get_count (&remote_tree->elements);
      memcpy (current_buffer + bytes_written, &element_count,
              sizeof (size_t));
      bytes_written += sizeof (size_t);
      bytes_written += T8_ADD_PADDING (bytes_written);
      /* The byte count of the elements */
      element_size = t8_element_array_get_size (&remote_tree->elements);
      element_bytes = element_size * element_count;
      /* Copy the elements into the send buffer */
      memcpy (current_buffer + bytes_written,
              t8_element_array_get_data (&remote_tree->elements),
              element_bytes);
      bytes_written += element_bytes;
      /* add padding after the elements */
      bytes_written += T8_ADD_PADDING (bytes_written);

      /* Add to the counter of remote elements. */
      ghost->num_remote_elements += element_count;
#ifdef T8_ENABLE_DEBUG
      acc_el_count += element_count;
#endif
    }                           /* End tree loop */

    T8_ASSERT (bytes_written == current_send_info->num_bytes);
    /* We can now post the MPI_Isend for the remote process */
    mpiret =
      sc_MPI_Isend (current_buffer, bytes_written, sc_MPI_BYTE, remote_rank,
                    T8_MPI_GHOST_FOREST, forest->mpicomm,
                    *requests + proc_index);
    SC_CHECK_MPI (mpiret);
  }                             /* end process loop */
  return send_info;
}

static void
t8_forest_ghost_send_end (t8_forest_t forest, t8_forest_ghost_t ghost,
                          t8_ghost_mpi_send_info_t *send_info,
                          sc_MPI_Request * requests)
{
  int                 num_remotes;
  int                 proc_pos, mpiret;

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

/* Receive a single message from a remote process, after the message was
 * successfully probed.
 * Returns the allocated receive buffer and the number of bytes received */
static char        *
t8_forest_ghost_receive_message (int recv_rank, sc_MPI_Comm comm,
                                 sc_MPI_Status status, int *recv_bytes)
{
  char               *recv_buffer;
  int                 mpiret;

  T8_ASSERT (recv_rank == status.MPI_SOURCE);
  T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);

  /* Get the number of bytes in the message */
  mpiret = sc_MPI_Get_count (&status, sc_MPI_BYTE, recv_bytes);

  /* Allocate receive buffer */
  recv_buffer = T8_ALLOC_ZERO (char, *recv_bytes);
  /* receive the message */
  mpiret = sc_MPI_Recv (recv_buffer, *recv_bytes, sc_MPI_BYTE, recv_rank,
                        T8_MPI_GHOST_FOREST, comm, sc_MPI_STATUS_IGNORE);
  SC_CHECK_MPI (mpiret);

  return recv_buffer;
}

/* Parse a message from a remote process and correctly include the received
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
 */
/* Currently we expect that the messages arrive in order of the sender's rank. */
static void
t8_forest_ghost_parse_received_message (t8_forest_t forest,
                                        t8_forest_ghost_t ghost,
                                        t8_locidx_t *current_element_offset,
                                        int recv_rank, char *recv_buffer,
                                        int recv_bytes)
{
  size_t              bytes_read, first_tree_index = 0, first_element_index =
    0;
  t8_locidx_t         num_trees, itree;
  t8_gloidx_t         global_id;
  t8_eclass_t         eclass;
  size_t              num_elements, old_elem_count, ghosts_offset;
  t8_ghost_gtree_hash_t *tree_hash, **pfound_tree, *found_tree;
  t8_ghost_tree_t    *ghost_tree;
  t8_eclass_scheme_c *ts;
  t8_element_t       *element_insert;
  t8_ghost_process_hash_t *process_hash;
#ifdef T8_ENABLE_DEBUG
  int                 added_process;
#endif

  bytes_read = 0;
  /* read the number of trees */
  num_trees = *(size_t *) recv_buffer;
  bytes_read += sizeof (size_t);
  bytes_read += T8_ADD_PADDING (bytes_read);

  t8_debugf ("Received %li trees from %i (%i bytes)\n",
             (long) num_trees, recv_rank, recv_bytes);

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
    tree_hash =
      (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    tree_hash->global_id = global_id;

    /* Get the element scheme for this tree */
    ts = t8_forest_get_eclass_scheme (forest, eclass);
    if (sc_hash_insert_unique (ghost->global_tree_to_ghost_tree, tree_hash,
                               (void ***) &pfound_tree)) {
      /* The tree was not stored already, tree_hash is now an entry in the hash table. */
      /* If the tree was not contained, it is the newest tree in the array and
       * thus has as index the number of currently inserted trees. */
      tree_hash->index = ghost->ghost_trees->elem_count;
      found_tree = tree_hash;
      /* We grow the array by one and initilize the entry */
      ghost_tree = (t8_ghost_tree_t *) sc_array_push (ghost->ghost_trees);
      ghost_tree->global_id = global_id;
      ghost_tree->eclass = eclass;
      /* Initialize the element array */
      t8_element_array_init_size (&ghost_tree->elements, ts, num_elements);
      /* pointer to where the elements are to be inserted */
      element_insert = t8_element_array_get_data (&ghost_tree->elements);
      /* Compute the element offset of this new tree by adding the offset
       * of the previous tree to the element count of the previous tree. */
      ghost_tree->element_offset = *current_element_offset;
      /* Allocate a new tree_hash for the next search */
      old_elem_count = 0;
      tree_hash =
        (t8_ghost_gtree_hash_t *) sc_mempool_alloc (ghost->glo_tree_mempool);
    }
    else {
      /* The entry was found in the trees array */
      found_tree = *pfound_tree;
      T8_ASSERT (found_tree->global_id == global_id);
      /* Get a pointer to the tree */
      ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees,
                                                       found_tree->index);
      T8_ASSERT (ghost_tree->eclass == eclass);
      T8_ASSERT (ghost_tree->global_id == global_id);
      T8_ASSERT (ghost_tree->elements.scheme == ts);

      old_elem_count = t8_element_array_get_count (&ghost_tree->elements);

      /* Grow the elements array of the tree to fit the new elements */
      t8_element_array_resize (&ghost_tree->elements,
                               old_elem_count + num_elements);
      /* Get a pointer to where the new elements are to be inserted */
      element_insert = t8_element_array_index_locidx (&ghost_tree->elements,
                                                      old_elem_count);
    }
    if (itree == 0) {
      /* We store the index of the first tree and the first element of this
       * rank */
      first_tree_index = found_tree->index;
      first_element_index = old_elem_count;
    }
    /* Insert the new elements */
    memcpy (element_insert, recv_buffer + bytes_read,
            num_elements * ts->t8_element_size ());

    bytes_read += num_elements * ts->t8_element_size ();
    bytes_read += T8_ADD_PADDING (bytes_read);
    *current_element_offset += num_elements;
  }
  T8_ASSERT (bytes_read == (size_t) recv_bytes);
  T8_FREE (recv_buffer);

  /* At last we add the receiving rank to the ghosts process_offset hash table */
  process_hash =
    (t8_ghost_process_hash_t *) sc_mempool_alloc (ghost->proc_offset_mempool);
  process_hash->mpirank = recv_rank;
  process_hash->tree_index = first_tree_index;
  process_hash->first_element = first_element_index;
  process_hash->ghost_offset = ghosts_offset;
  /* Insert this rank into the hash table. We assert if the rank was not already
   * contained. */
#ifdef T8_ENABLE_DEBUG
  added_process =
#else
  (void)
#endif
    sc_hash_insert_unique (ghost->process_offsets, process_hash, NULL);
  T8_ASSERT (added_process);
}

/* In forest_ghost_receive we need a lookup table to give us the position
 * of a process in the ghost->remote_processes array, given the rank of
 * a process. We implement this via a hash table with the following struct
 * as entry. */
typedef struct t8_recv_list_entry_struct
{
  int                 rank;     /* The rank of this process */
  int                 pos_in_remote_processes;  /* The position of this process in the remote_processes array */
} t8_recv_list_entry_t;

/* We hash these entries by their rank */
unsigned
t8_recv_list_entry_hash (const void *v1, const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;

  return e1->rank;
}

/* two entries are considered equal if they have the same rank. */
int
t8_recv_list_entry_equal (const void *v1, const void *v2, const void *u)
{
  const t8_recv_list_entry_t *e1 = (const t8_recv_list_entry_t *) v1;
  const t8_recv_list_entry_t *e2 = (const t8_recv_list_entry_t *) v2;

  return e1->rank == e2->rank;
}

/* Probe for all incoming messages from the remote ranks and receive them.
 * We receive the message in the order in which they arrive. To achieve this,
 * we have to use polling. */
static void
t8_forest_ghost_receive (t8_forest_t forest, t8_forest_ghost_t ghost)
{
  int                 num_remotes;
  int                 proc_pos;
  int                 recv_rank;
  int                 mpiret;
  sc_MPI_Comm         comm;
  sc_MPI_Status       status;

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
    /*       This code receives the message in order of their arrival.
     *       This is effective in terms of runtime, but makes it more difficult
     *       to store the received data, since the data has to be stored in order of
     *       ascending ranks.
     *       We include the received data into the ghost structure in order of the
     *       ranks of the receivers and we do this as soon as the message from
     *       the next rank that we can include was received. */
    char              **buffer;
    int                *recv_bytes;
    int                 received_messages = 0;
    int                *received_flag;
    int                 last_rank_parsed = -1, parse_it;
#undef T8_POLLING               /* activates polling for Mpi messages, currently used for testing */
#ifdef T8_POLLING
    sc_link_t          *proc_it, *prev;
    int                 iprobe_flag;
    sc_list_t          *receivers;
#else
    t8_recv_list_entry_t **pfound, *found;
#ifdef T8_ENABLE_DEBUG
    int                 ret;
#endif
    sc_hash_t          *recv_list_entries_hash;
#endif
    t8_recv_list_entry_t recv_list_entry, *recv_list_entries;
    t8_locidx_t         current_element_offset = 0;

    buffer = T8_ALLOC (char *, num_remotes);
    recv_bytes = T8_ALLOC (int, num_remotes);
    received_flag = T8_ALLOC_ZERO (int, num_remotes);
    recv_list_entries = T8_ALLOC (t8_recv_list_entry_t, num_remotes);

    /* Sort the array of remote processes, such that the ranks are in
     * ascending order. */
    sc_array_sort (ghost->remote_processes, sc_int_compare);

    /* We build a hash table of all ranks from which we receive and their position
     * in the remote_processes array. */
#ifdef T8_POLLING               /* polling */
    receivers = sc_list_new (NULL);
#else
    recv_list_entries_hash = sc_hash_new (t8_recv_list_entry_hash,
                                          t8_recv_list_entry_equal, NULL,
                                          NULL);
#endif
    for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
      recv_list_entries[proc_pos].rank =
        *(int *) sc_array_index_int (ghost->remote_processes, proc_pos);
      recv_list_entries[proc_pos].pos_in_remote_processes = proc_pos;
#ifndef T8_POLLING
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_insert_unique (recv_list_entries_hash,
                               recv_list_entries + proc_pos, NULL);
      T8_ASSERT (ret == 1);
#else /* polling */
      sc_list_append (receivers, recv_list_entries + proc_pos);
#endif
    }

  /****     Actual communication    ****/

    /* Until there is only one sender left we iprobe for a message for each
     * sender and if there is one we receive it and remove the sender from
     * the list.
     * The last message can be received via probe */
#ifdef T8_POLLING
    while (received_messages < num_remotes - 1) {
      /* TODO: This part of the code using polling and IProbe to receive the
       *       messages. We replaced with a non-polling version that uses the
       *       blocking Probe. */
      iprobe_flag = 0;
      prev = NULL;              /* ensure that if the first receive entry is matched first,
                                   it is removed properly. */
      for (proc_it = receivers->first; proc_it != NULL && iprobe_flag == 0;) {
        /* pointer to the rank of a receiver */
        recv_rank = ((t8_recv_list_entry_t *) proc_it->data)->rank;
        proc_pos =
          ((t8_recv_list_entry_t *) proc_it->data)->pos_in_remote_processes;

        mpiret = sc_MPI_Iprobe (recv_rank, T8_MPI_GHOST_FOREST, comm,
                                &iprobe_flag, &status);
        SC_CHECK_MPI (mpiret);
#else
    while (received_messages < num_remotes) {
      /* blocking probe for a message. */
      mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_GHOST_FOREST, comm,
                             &status);
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
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (recv_list_entries_hash, &recv_list_entry,
                        (void ***) &pfound);
      T8_ASSERT (ret != 0);
      found = *pfound;
      proc_pos = found->pos_in_remote_processes;
#endif
      T8_ASSERT (status.MPI_TAG == T8_MPI_GHOST_FOREST);
      buffer[proc_pos] =
        t8_forest_ghost_receive_message (recv_rank, comm, status,
                                         recv_bytes + proc_pos);
      /* mark this entry as received. */
      T8_ASSERT (received_flag[proc_pos] == 0);
      received_flag[proc_pos] = 1;
      received_messages++;
      /* Parse all messages that we can parse now.
       * We have to parse the messages in order of their rank. */
      T8_ASSERT (last_rank_parsed < proc_pos);
      /* For all ranks that we haven't parsed yet, but can be parsed in order */
      for (parse_it = last_rank_parsed + 1; parse_it < num_remotes &&
           received_flag[parse_it] == 1; parse_it++) {
        recv_rank =
          *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
        t8_forest_ghost_parse_received_message (forest, ghost,
                                                &current_element_offset,
                                                recv_rank, buffer[parse_it],
                                                recv_bytes[parse_it]);
        last_rank_parsed++;
      }

#ifdef T8_POLLING               /* polling */
      /* Remove the process from the list of receivers. */
      proc_it = proc_it->next;
      sc_list_remove (receivers, prev);
    }
  }                             /* end for */
#endif
#if 0
  {                             /* this is for indent */
    {
#endif
    }                           /* end while */
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
    buffer[proc_pos] = t8_forest_ghost_receive_message (recv_rank, comm,
                                                        status,
                                                        recv_bytes +
                                                        proc_pos);
    received_flag[proc_pos] = 1;
    received_messages++;
    T8_ASSERT (received_messages == num_remotes);
    /* parse all messages that are left */
    /* For all ranks that we haven't parsed yet, but can be parsed in order */
    for (parse_it = last_rank_parsed + 1; parse_it < num_remotes &&
         received_flag[parse_it] == 1; parse_it++) {
      recv_rank =
        *(int *) sc_array_index_int (ghost->remote_processes, parse_it);
      t8_forest_ghost_parse_received_message (forest, ghost,
                                              &current_element_offset,
                                              recv_rank, buffer[parse_it],
                                              recv_bytes[parse_it]);
      last_rank_parsed++;
    }
#endif
#ifdef T8_ENABLE_DEBUG
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

#if 0
  /* Receive the message in order of the sender's rank,
   * this is a non-optimized version of the code below.
   * slow but simple. */

  /* Sort the array of remote processes, such that the ranks are in
   * ascending order. */
  sc_array_sort (ghost->remote_processes, sc_int_compare);

  for (proc_pos = 0; proc_pos < num_remotes; proc_pos++) {
    recv_rank =
      (int *) sc_array_index_int (ghost->remote_processes, proc_pos);
    /* blocking probe for a message. */
    mpiret = sc_MPI_Probe (*recv_rank, T8_MPI_GHOST_FOREST, comm, &status);
    SC_CHECK_MPI (mpiret);
    /* receive message */
    buffer =
      t8_forest_ghost_receive_message (*recv_rank, comm, status, &recv_bytes);
    t8_forest_ghost_parse_received_message (forest, ghost, *recv_rank, buffer,
                                            recv_bytes);
  }
#endif
}

/* Create one layer of ghost elements, following the algorithm
 * in: p4est: Scalable Algorithms For Parallel Adaptive
 *     Mesh Refinement On Forests of Octrees
 *     C. Burstedde, L. C. Wilcox, O. Ghattas
 * for unbalanced_version = 0 (balanced forest only) or
 *     Recursive algorithms for distributed forests of octrees
 *     T. Isaac, C. Burstedde, L. C. Wilcox and O. Ghattas
 * for unbalanced_version = 1 (also unbalanced forests possible).
 *
 * verion 3 with top-down search
 * for unbalanced_version = -1
 */
void
t8_forest_ghost_create_ext (t8_forest_t forest, int unbalanced_version)
{
  t8_forest_ghost_t   ghost = NULL;
  t8_ghost_mpi_send_info_t *send_info;
  sc_MPI_Request     *requests;
  int                 create_tree_array = 0, create_gfirst_desc_array = 0;
  int                 create_element_array = 0;

  T8_ASSERT (t8_forest_is_committed (forest));

  t8_global_productionf ("Into t8_forest_ghost with %i local elements.\n",
                         t8_forest_get_local_num_elements (forest));

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start ghost at %f  %f\n", sc_MPI_Wtime (),
                           forest->profile->ghost_runtime);
  }

  if (forest->element_offsets == NULL) {
    /* create element offset array if not done already */
    create_element_array = 1;
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == NULL) {
    /* Create tree offset array if not done already */
    create_tree_array = 1;
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    /* Create global first desc array if not done already */
    create_gfirst_desc_array = 1;
    t8_forest_partition_create_first_desc (forest);
  }

  if (t8_forest_get_local_num_elements (forest) > 0) {
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

    /* Reveive the ghost elements from the remote processes */
    t8_forest_ghost_receive (forest, ghost);

    /* End sending the remote elements */
    t8_forest_ghost_send_end (forest, ghost, send_info, requests);

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

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of ghost_create */
    forest->profile->ghost_runtime += sc_MPI_Wtime ();
    /* We also store the number of ghosts and remotes */
    if (ghost != NULL) {
      forest->profile->ghosts_received = ghost->num_ghosts_elements;
      forest->profile->ghosts_shipped = ghost->num_remote_elements;
    }
    else {
      forest->profile->ghosts_received = 0;
      forest->profile->ghosts_shipped = 0;
    }
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End ghost at %f  %f\n", sc_MPI_Wtime (),
                           forest->profile->ghost_runtime);
  }

  t8_global_productionf ("Done t8_forest_ghost with %i local elements and %i"
                         " ghost elements.\n",
                         t8_forest_get_local_num_elements (forest),
                         t8_forest_get_num_ghosts (forest));
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
int                *
t8_forest_ghost_get_remotes (t8_forest_t forest, int *num_remotes)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->ghosts == NULL) {
    *num_remotes = 0;
    return NULL;
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

/* Fill the send buffer for a ghost data exchange for on remote rank.
 * returns the number of bytes in the buffer. */
static size_t
t8_forest_ghost_exchange_fill_send_buffer (t8_forest_t forest, int remote,
                                           char **pbuffer,
                                           sc_array_t *element_data)
{
  char               *buffer;
  t8_ghost_remote_t   lookup_rank, *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;
  t8_forest_ghost_t   ghost;
  size_t              index, element_index, data_size;
  size_t              elements_inserted, byte_count;
  t8_tree_t           local_tree;
#ifdef T8_ENABLE_DEBUG
  int                 ret;
#endif
  t8_locidx_t         itree, ielement, element_pos;
  t8_locidx_t         ltreeid;
  size_t              elem_count;

  ghost = forest->ghosts;
  data_size = element_data->elem_size;
  elements_inserted = 0;
  lookup_rank.remote_rank = remote;

  /* Lookup the remote entry of this remote process */
#ifdef T8_ENABLE_DEBUG
  ret =
#else
  (void)
#endif
    sc_hash_array_lookup (ghost->remote_ghosts, &lookup_rank, &index);
  T8_ASSERT (ret != 0);
  remote_entry =
    (t8_ghost_remote_t *) sc_array_index (&ghost->remote_ghosts->a, index);
  T8_ASSERT (remote_entry->remote_rank == remote);

  /* allocate memory for the send buffer */
  byte_count = data_size * remote_entry->num_elements;
  buffer = *pbuffer = T8_ALLOC (char, byte_count);

  /* We now iterate over the remote trees and their elements to find the
   * local element indices of the remote elements */
  for (itree = 0; itree < (t8_locidx_t) remote_entry->remote_trees.elem_count;
       itree++) {
    /* tree loop */
    remote_tree = (t8_ghost_remote_tree_t *)
      t8_sc_array_index_locidx (&remote_entry->remote_trees, itree);
    /* Get the local id of this tree */
    /* TODO: Why does remote_tree store the global id? could be local instead */
    ltreeid = t8_forest_get_local_id (forest, remote_tree->global_id);
    /* Get a pointer to the forest tree */
    local_tree = t8_forest_get_tree (forest, ltreeid);
    elem_count = t8_element_array_get_count (&remote_tree->elements);
    for (ielement = 0; ielement < (t8_locidx_t) elem_count; ielement++) {
      /* element loop */
      /* Get the index of this remote element in its local tree */
      element_pos = *(t8_locidx_t *)
        t8_sc_array_index_locidx (&remote_tree->element_indices, ielement);
      T8_ASSERT (0 <= element_pos);
      /* Compute the index of this element in the element_data array */
      element_index = local_tree->elements_offset + element_pos;
      /* Copy the data of this element from the element_data array to the send buffer */
      memcpy (buffer + elements_inserted * data_size,
              sc_array_index (element_data, element_index), data_size);
      elements_inserted++;
    }
  }
  return byte_count;
}

static t8_ghost_data_exchange_t *
t8_forest_ghost_exchange_begin (t8_forest_t forest, sc_array_t *element_data)
{
  t8_ghost_data_exchange_t *data_exchange;
  t8_forest_ghost_t   ghost;
  size_t              bytes_to_send, ghost_start;
  int                 iremote, remote_rank;
  int                 mpiret, recv_rank, bytes_recv;
#ifdef T8_ENABLE_DEBUG
  int                 ret;
#endif
  char              **send_buffers;
  t8_ghost_process_hash_t lookup_proc, *process_entry, **pfound;
  t8_locidx_t         remote_offset, next_offset;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (element_data != NULL);
  T8_ASSERT (forest->ghosts != NULL);

  ghost = forest->ghosts;

  /* Allocate the new exchange context */
  data_exchange = T8_ALLOC (t8_ghost_data_exchange_t, 1);
  /* The number of processes we need to send to */
  data_exchange->num_remotes = ghost->remote_processes->elem_count;
  /* Allocate MPI requests */
  data_exchange->send_requests = T8_ALLOC (sc_MPI_Request,
                                           data_exchange->num_remotes);
  data_exchange->recv_requests = T8_ALLOC (sc_MPI_Request,
                                           data_exchange->num_remotes);
  /* Allocate pointers to send buffers */
  send_buffers = data_exchange->send_buffers =
    T8_ALLOC (char *, data_exchange->num_remotes);

  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* Iterate over all remote processes and fill their send buffers */
    remote_rank =
      *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    /* Fill the send buffers and compute the number of bytes to send */
    bytes_to_send =
      t8_forest_ghost_exchange_fill_send_buffer (forest, remote_rank,
                                                 send_buffers + iremote,
                                                 element_data);

    /* Post the asynchronuos send */
    mpiret = sc_MPI_Isend (send_buffers[iremote], bytes_to_send, sc_MPI_BYTE,
                           remote_rank, T8_MPI_GHOST_EXC_FOREST,
                           forest->mpicomm,
                           data_exchange->send_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }

  /* The index in element_data at which the ghost elements start */
  ghost_start = t8_forest_get_local_num_elements (forest);
  /* Receive the incoming messages */
#if 0
  while (received_messages < data_exchange->num_remotes) {
    /* Blocking test for incoming message */
    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, T8_MPI_GHOST_EXC_FOREST,
                           forest->mpicomm, &recv_status);
    SC_CHECK_MPI (mpiret);
    recv_rank = recv_status.MPI_SOURCE;
    /* Get the number of bytes to receive */
    mpiret = sc_MPI_Get_count (&recv_status, sc_MPI_BYTE, &bytes_recv);
    SC_CHECK_MPI (mpiret);

    /* We need to compute the offset in element_data to which we can receive the message */
    /* Search for this process' entry in the ghost struct */
    lookup_proc.mpirank = recv_rank;
    ret =
      sc_hash_lookup (ghost->process_offsets, &lookup_proc,
                      (void ***) &pfound);
    T8_ASSERT (ret);
    process_entry = *pfound;
    /* In process_entry we stored the offset of this ranks ghosts under all
     * ghosts. Thus in element_data we look at the position
     *  ghost_start + offset
     */
    /* receive the message */
    sc_MPI_Recv (sc_array_index
                 (element_data, ghost_start + process_entry->ghost_offset),
                 bytes_recv, sc_MPI_BYTE, recv_rank, T8_MPI_GHOST_EXC_FOREST,
                 forest->mpicomm, sc_MPI_STATUS_IGNORE);
    received_messages++;
  }
#endif
  for (iremote = 0; iremote < data_exchange->num_remotes; iremote++) {
    /* We need to compute the offset in element_data to which we can receive the message */
    /* Search for this processes' entry in the ghost struct */
    recv_rank =
      *(int *) sc_array_index_int (ghost->remote_processes, iremote);
    lookup_proc.mpirank = recv_rank;
#ifdef T8_ENABLE_DEBUG
    ret =
#else
    (void)
#endif
      sc_hash_lookup (ghost->process_offsets, &lookup_proc,
                      (void ***) &pfound);
    T8_ASSERT (ret);
    process_entry = *pfound;
    /* In process_entry we stored the offset of this ranks ghosts under all
     * ghosts. Thus in element_data we look at the position
     *  ghost_start + offset
     */
    remote_offset = process_entry->ghost_offset;
    /* Compute the offset of the next remote rank */
    if (iremote + 1 < data_exchange->num_remotes) {
      lookup_proc.mpirank =
        *(int *) sc_array_index_int (ghost->remote_processes, iremote + 1);
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &lookup_proc,
                        (void ***) &pfound);
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
    mpiret =
      sc_MPI_Irecv (sc_array_index
                    (element_data, ghost_start + remote_offset), bytes_recv,
                    sc_MPI_BYTE, recv_rank, T8_MPI_GHOST_EXC_FOREST,
                    forest->mpicomm, data_exchange->recv_requests + iremote);
    SC_CHECK_MPI (mpiret);
  }
  return data_exchange;
}

static void
t8_forest_ghost_exchange_end (t8_ghost_data_exchange_t *data_exchange)
{
  int                 iproc;

  T8_ASSERT (data_exchange != NULL);
  /* Wait for all communications to end */
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->recv_requests,
                  sc_MPI_STATUSES_IGNORE);
  sc_MPI_Waitall (data_exchange->num_remotes, data_exchange->send_requests,
                  sc_MPI_STATUSES_IGNORE);

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

  if (forest->ghosts == NULL) {
    /* This process has no ghosts */
    return;
  }

  T8_ASSERT (forest->ghosts != NULL);
  T8_ASSERT (element_data != NULL);
  T8_ASSERT ((t8_locidx_t) element_data->elem_count ==
             t8_forest_get_local_num_elements (forest)
             + t8_forest_get_num_ghosts (forest));

  data_exchange = t8_forest_ghost_exchange_begin (forest, element_data);
  if (forest->profile != NULL) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime = -sc_MPI_Wtime ();
  }
  t8_forest_ghost_exchange_end (data_exchange);
  if (forest->profile != NULL) {
    /* Measure the time for ghost_exchange_end */
    forest->profile->ghost_waittime += sc_MPI_Wtime ();
  }
  t8_debugf ("Finished ghost_exchange_data\n");
}

/* Print a forest ghost structure */
void
t8_forest_ghost_print (t8_forest_t forest)
{
  t8_forest_ghost_t   ghost;
  t8_ghost_remote_t  *remote_found;
  t8_ghost_remote_tree_t *remote_tree;
  t8_ghost_process_hash_t proc_hash, **pfound, *found;
  size_t              iremote, itree;
#ifdef T8_ENABLE_DEBUG
  int                 ret;
#endif
  int                 remote_rank;
  char                remote_buffer[BUFSIZ] = "";
  char                buffer[BUFSIZ] = "";

  if (forest->ghosts == NULL) {
    return;
  }
  T8_ASSERT (forest->ghosts != NULL);
  ghost = forest->ghosts;
  snprintf (remote_buffer + strlen (remote_buffer),
            BUFSIZ - strlen (remote_buffer), "\tRemotes:\n");
  snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
            "\tReceived:\n");

  if (ghost->num_ghosts_elements > 0) {
    for (iremote = 0; iremote < ghost->remote_processes->elem_count;
         iremote++) {
      /* Get the rank of the remote process */
      remote_rank =
        *(int *) sc_array_index (ghost->remote_processes, iremote);
      /* Get this remote's entry */
      remote_found = t8_forest_ghost_get_remote (forest, remote_rank);
      /* investigate the entry of this remote process */
      snprintf (remote_buffer + strlen (remote_buffer),
                BUFSIZ - strlen (remote_buffer), "\t[Rank %i] (%li trees):\n",
                remote_found->remote_rank,
                remote_found->remote_trees.elem_count);
      for (itree = 0; itree < remote_found->remote_trees.elem_count; itree++) {
        remote_tree = (t8_ghost_remote_tree_t *)
          sc_array_index (&remote_found->remote_trees, itree);
        snprintf (remote_buffer + strlen (remote_buffer),
                  BUFSIZ - strlen (remote_buffer),
                  "\t\t[id: %lli, class: %s, #elem: %li]\n",
                  (long long) remote_tree->global_id,
                  t8_eclass_to_string[remote_tree->eclass],
                  (long) t8_element_array_get_count (&remote_tree->elements));
      }

      /* Investigate the elements that we received from this process */
      proc_hash.mpirank = remote_rank;
      /* look up this rank in the hash table */
#ifdef T8_ENABLE_DEBUG
      ret =
#else
      (void)
#endif
        sc_hash_lookup (ghost->process_offsets, &proc_hash,
                        (void ***) &pfound);

      T8_ASSERT (ret);
      found = *pfound;
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer),
                "\t[Rank %i] First tree: %li\n\t\t First element: %li\n",
                remote_rank,
                (long) found->tree_index, (long) found->first_element);
    }
  }
  t8_debugf ("Ghost structure:\n%s\n%s\n", remote_buffer, buffer);
}

/* Completely destroy a ghost structure */
static void
t8_forest_ghost_reset (t8_forest_ghost_t *pghost)
{
  t8_forest_ghost_t   ghost;
  size_t              it, it_trees;
  t8_ghost_tree_t    *ghost_tree;
  t8_ghost_remote_t  *remote_entry;
  t8_ghost_remote_tree_t *remote_tree;

  T8_ASSERT (pghost != NULL);
  ghost = *pghost;
  T8_ASSERT (ghost != NULL);
  T8_ASSERT (ghost->rc.refcount == 0);

  /* Clean-up the arrays */
  for (it_trees = 0; it_trees < ghost->ghost_trees->elem_count; it_trees++) {
    ghost_tree = (t8_ghost_tree_t *) sc_array_index (ghost->ghost_trees,
                                                     it_trees);
    t8_element_array_reset (&ghost_tree->elements);
  }

  sc_array_destroy (ghost->ghost_trees);
  sc_array_destroy (ghost->remote_processes);
  /* Clean-up the hashtables */
  sc_hash_destroy (ghost->global_tree_to_ghost_tree);
  sc_hash_destroy (ghost->process_offsets);
  /* Clean-up the remote ghost entries */
  for (it = 0; it < ghost->remote_ghosts->a.elem_count; it++) {
    remote_entry = (t8_ghost_remote_t *)
      sc_array_index (&ghost->remote_ghosts->a, it);
    for (it_trees = 0; it_trees < remote_entry->remote_trees.elem_count;
         it_trees++) {
      remote_tree = (t8_ghost_remote_tree_t *)
        sc_array_index (&remote_entry->remote_trees, it_trees);
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
  pghost = NULL;
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
  t8_forest_ghost_t   ghost;

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
  T8_ASSERT (pghost != NULL && *pghost != NULL &&
             t8_refcount_is_last (&(*pghost)->rc));
  t8_forest_ghost_unref (pghost);
  T8_ASSERT (*pghost == NULL);
}

T8_EXTERN_C_END ();
