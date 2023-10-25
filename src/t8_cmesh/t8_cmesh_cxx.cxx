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

#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include "t8_cmesh_types.h"
#include <t8_data/t8_shmem.h>

/** \file t8_cmesh_cxx.cxx
 *  This file collects all general cmesh routines that need c++ compilation.
 *  Particularly those functions that use the element interface from \ref t8_element_cxx.hxx.
 *
 * TODO: document this file
 */

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
  t8_debugf ("[D] proc is empty\n");
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
#ifdef T8_ENABLE_DEBUG
  {
    /* We check whether computing A/C * B will cause an overflow.
     * This can be achieved by checking if dividing the result by B
     * yields A/C again. */
    const t8_gloidx_t a_over_c = A / C;
    const t8_gloidx_t a_o_c_times_b = a_over_c * B;
    T8_ASSERT (a_over_c == 0 || a_o_c_times_b / a_over_c == B);
  }
#endif
  return (t8_gloidx_t) ((A / C) * B + (((long double) (A % C)) / C) * B);
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
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level, t8_scheme_cxx_t *ts, t8_gloidx_t *first_local_tree,
                         t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                         int8_t *first_tree_shared)
{
  int is_empty;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (level >= 0);
  T8_ASSERT (ts != NULL);

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
  t8_eclass_scheme_c *tree_scheme;

  /* Compute the number of children on level in each tree */
  global_num_children = 0;
  for (tree_class = T8_ECLASS_ZERO; tree_class < T8_ECLASS_COUNT; ++tree_class) {
    /* We iterate over each element class and get the number of children for this
     * tree class.
     */
    if (cmesh->num_trees_per_eclass[tree_class] > 0) {
      tree_scheme = ts->eclass_schemes[tree_class];
      T8_ASSERT (tree_scheme != NULL);
      children_per_tree = tree_scheme->t8_element_count_leafs_from_root (level);
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
 * storing the number of processe and the global number of elements. 
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

  t8_debugf ("[H] tree %li el %li mirror %li\n", pure_local_tree, element_index, mirror_element_index);
  if (element_index == query_data->global_num_elements) {
    return query_data->num_procs - query_data->process_offset;
  }
  else {
    first_proc_rank
      = query_data->num_procs - 1
        - t8_A_times_B_over_C_intA (query_data->num_procs, mirror_element_index, query_data->global_num_elements);

    //       (int)(((long double)mirror_element_index/query_data->global_num_elements)*query_data->num_procs);

    /* If this is called by array split, we need to adjust for relative processes starting
     * add the first process. */
    first_proc_adjusted = first_proc_rank - query_data->process_offset;
  }
  t8_debugf ("[H] ptree %zd, first element %li, on proc %zd\n", pure_local_tree, element_index, first_proc_adjusted);

  /* Safety checks */
  T8_ASSERT (0 <= first_proc_rank && (int) first_proc_rank < query_data->num_procs);
  T8_ASSERT (0 <= first_proc_adjusted);
  /* Check that the element lies in the partition of the computed proc. */
  T8_ASSERT (
    t8_cmesh_get_first_element_of_process (first_proc_rank, query_data->num_procs, query_data->global_num_elements)
    <= element_index);
#ifdef T8_ENABLE_DEBUG
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

/* TODO: Empty processes, shared trees, binary search in offset-array to avoid recv_any,
 * use partition_given to partition the cmesh*/
void
t8_cmesh_uniform_bounds_hybrid (t8_cmesh_t cmesh, int level, t8_scheme_cxx_t *scheme, t8_gloidx_t *first_local_tree,
                                t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree,
                                t8_gloidx_t *child_in_tree_end, int8_t *first_tree_shared, sc_MPI_Comm comm)
{
  t8_gloidx_t local_num_children = 0;
  t8_gloidx_t *elem_index_pointer;
  int send_first, send_last, send_first_nonempty, num_procs_we_send_to = 0;
  t8_cmesh_partition_query_t data;
  t8_shmem_array_t offset_array;
  sc_array_t offset_partition, first_element_tree;
  int mpiret, num_messages_sent = 0;
  int expect_start_message = 1;
  int expect_end_message = 1;
  sc_array send_requests, send_buffer;
  int current_pos_in_send_buffer = 0;
  t8_eclass_scheme_c *tree_scheme;
  T8_ASSERT (cmesh != NULL);
  const int first_tree_shared_shift = cmesh->first_tree_shared ? 1 : 0;
#ifdef T8_ENABLE_DEBUG
  int num_received_start_messages = 0;
  int num_received_end_messages = 0;
#endif

  /* TODO: Clean up size_t and gloidx_t data types, ensure that each variables has the 
   *          matching type. */

  t8_debugf ("Into t8_cmesh_uniform_bounds_hybrid.\n");

  if (t8_cmesh_is_empty (cmesh)) {
    t8_debugf ("[D] cmesh_uniform_bounds_hybrid is empty\n");
    t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                     child_in_tree_end, first_tree_shared);
    return;
  }

#ifdef T8_ENABLE_DEBUG
  {
    /* Check that comm matches mpirank and size stored in cmesh */
    int mpirank, mpisize;
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (mpirank == cmesh->mpirank);
    T8_ASSERT (mpisize == cmesh->mpisize);
  }
#endif
  /*Compute number of local elements. Ignore shared trees */
  for (int ieclass = T8_ECLASS_ZERO; ieclass < T8_ECLASS_COUNT; ieclass++) {
    tree_scheme = scheme->eclass_schemes[ieclass];
    local_num_children
      += cmesh->num_local_trees_per_eclass[ieclass] * tree_scheme->t8_element_count_leafs_from_root (level);
  }

  /* Do not consider shared trees */
  if (cmesh->first_tree_shared && cmesh->set_partition) {
    const int ieclass = t8_cmesh_get_tree_class (cmesh, 0);
    tree_scheme = scheme->eclass_schemes[ieclass];
    local_num_children -= tree_scheme->t8_element_count_leafs_from_root (level);
  }

  /* 
   *
   *  Cmesh is not partitioned
   * 
   */
  /* If the initial cmesh is not partitioned, every process knows "everything" and we do not
   * need any communication.*/
  if (!cmesh->set_partition) {
    t8_gloidx_t current_tree_element_offset;
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
    t8_debugf ("[H] hybrid fc = %li, lc = %li\n", first_child, last_child);

    /* Can't we optimize this linear loop by using a binary search?
     * -> No, we cannot. Since we need in any case compute the t8_element_count_leafs_from_root
     *    for each tree.
     */
    current_tree_element_offset = 0;
    for (t8_gloidx_t igtree = 0; igtree < num_trees; ++igtree) {
      const int ieclass = t8_cmesh_get_tree_class (cmesh, (t8_locidx_t) igtree);
      tree_scheme = scheme->eclass_schemes[ieclass];
      /* TODO: We can optimize by buffering the elem_in_tree value. Thus, if 
         the computation is expensive (may be for non-morton-type schemes),
         we do it only once. */
      const t8_gloidx_t elem_in_tree = tree_scheme->t8_element_count_leafs_from_root (level);
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
        t8_debugf ("[D] lc < fc\n");
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
    t8_debugf ("[D] cmesh is empty\n");
    t8_cmesh_uniform_set_return_parameters_to_empty (first_local_tree, child_in_tree_begin, last_local_tree,
                                                     child_in_tree_end, first_tree_shared);
    return;
  }

  /* 
   *
   *  Cmesh is partitioned
   * 
   */

  t8_shmem_array_init (&offset_array, sizeof (t8_gloidx_t), cmesh->mpisize + 1, comm);
  /* Fill the offset array for each process the global index of its first element in
   * the uniform partition.
   * (0, l_n_c_0, l_n_c_0 + l_n_c_1, l_n_c_0 + l_n_c_1 + l_n_c_2, ...) */
  t8_shmem_prefix (&local_num_children, offset_array, 1, T8_MPI_GLOIDX, sc_MPI_SUM);

  data.num_procs = cmesh->mpisize;
  /* Get global number of elements */
  data.global_num_elements = t8_shmem_array_get_gloidx (offset_array, cmesh->mpisize);
  SC_CHECK_ABORTF (0 <= data.global_num_elements && data.global_num_elements < T8_GLOIDX_MAX,
                   "Overflow in number of elements.\n");
  t8_debugf ("[H] global num %li\n", data.global_num_elements);

  /*Compute number of non-shared-trees and the local index of the first non-shared-tree */
  const t8_locidx_t pure_local_trees = cmesh->num_local_trees - first_tree_shared_shift;

  if (pure_local_trees > 0) {
    /* Compute which trees and elements to send to which process.
     * We skip empty processes. */
    t8_locidx_t igtree = first_tree_shared_shift;

    sc_array_init_size (&first_element_tree, sizeof (t8_gloidx_t), pure_local_trees + 1);

    /* Set the first entry of first_element_tree to the global index of
     * the first element of our first pure local tree. */
    elem_index_pointer = (t8_gloidx_t *) sc_array_index_int (&first_element_tree, 0);
    *elem_index_pointer = t8_shmem_array_get_gloidx (offset_array, cmesh->mpirank);

    /* Compute the first element in every pure local tree.
     * This array stores for each tree the global element index offset. 
     * Example: 2 local trees, each has 8 Elements. First element index 12: | 12 | 20 | 28 | */
    for (t8_locidx_t itree = 0; itree < pure_local_trees; ++itree, ++igtree) {
      const t8_gloidx_t *first_element_of_tree = (const t8_gloidx_t *) sc_array_index_int (&first_element_tree, itree);
      t8_gloidx_t *first_element_of_next_tree = (t8_gloidx_t *) sc_array_index_int (&first_element_tree, itree + 1);
      const int ieclass = t8_cmesh_get_tree_class (cmesh, igtree);
      tree_scheme = scheme->eclass_schemes[ieclass];
      /* Set the first element of the next tree by adding the number of element of the current tree. */
      *first_element_of_next_tree = *first_element_of_tree + tree_scheme->t8_element_count_leafs_from_root (level);
    }

    /* Check that the last tree + 1 stores as offset the start of the next process */
    T8_ASSERT (*((t8_gloidx_t *) t8_sc_array_index_locidx (&first_element_tree, pure_local_trees))
               == t8_shmem_array_get_gloidx (offset_array, cmesh->mpirank + 1));

    /* Compute the range of mpiranks to minimize the loop over processes */

    /* Process offset must be 0 when calling t8_cmesh_determine_partition
     * the first times out of sc_array_split. */
    data.process_offset = 0;
    /* Compute the process that will own the first element of our first tree. */
    send_first_nonempty = (t8_gloidx_t) t8_cmesh_determine_partition (&first_element_tree, 0, &data);
    if (!cmesh->first_tree_shared && t8_cmesh_get_global_id (cmesh, 0) == 0) {
      /* If our first tree is 0 and not shared, then we need to send to all processes
       * below the start process. Even if their partition is empty. */
      t8_debugf("LLL: Would have set send_first to 0, instead: %i", send_first_nonempty);
      send_first = send_first_nonempty;
    }
    else {
      send_first = send_first_nonempty;
    }

    /* Compute the process that may own the last element of our first tree. */
    send_last = (t8_gloidx_t) t8_cmesh_determine_partition (&first_element_tree, pure_local_trees, &data);
    if (send_last >= cmesh->mpisize) {
      /* t8_cmesh_determine_partition will return mpisize if we plug in the number of global
       * trees as tree index, which happens on the last process that has data.
       * We need to correct by subtracting 1. */
      send_last = cmesh->mpisize - 1;
    }
    num_procs_we_send_to = send_last - send_first + 1;

    sc_array_init_size (&offset_partition, sizeof (size_t), num_procs_we_send_to);
    /* In array split the 'types' that we compute are the processes we send to,
     * but offset by the first process, so that they start at 0 and end at num_procs_we_send_to. */
    data.process_offset = send_first;
    /* For each process that we send to find the first tree whose first element
     * belongs to this process.
     * These tree indices will be stored in offset_partition. */
    sc_array_split (&first_element_tree, &offset_partition, (size_t) num_procs_we_send_to + 1,
                    t8_cmesh_determine_partition, (void *) &data);
    // TODO: Is determination of last proc we send to wrong?

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
    sc_array_init (&send_buffer, sizeof (t8_gloidx_t));
    current_pos_in_send_buffer = 0;

    /* Iterate over offset_partition to find boundaries
     * and send the MPI messages. */
    t8_debugf ("[H] sf: %i  sl: %i\n", send_first, send_last);

    const t8_gloidx_t last_el_index_of_last_tree
      = *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, pure_local_trees) - 1;
    for (int iproc = send_first; iproc <= send_last; iproc++) {

      const t8_gloidx_t first_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (iproc, data.num_procs, data.global_num_elements);
      const t8_gloidx_t last_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (iproc + 1, data.num_procs, data.global_num_elements) - 1;
      const int proc_is_empty = last_element_index_of_current_proc < first_element_index_of_current_proc;
      int send_start_message = 1;
      int send_end_message = 1;

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
          send_start_message = send_end_message = 0;
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
              send_start_message = 0;
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
            send_end_message = 0;
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
        t8_debugf ("[H] Set tree for %i: %i to %i, (#pure = %i)\n", iproc, first_puretree_of_current_proc,
                   last_puretree_of_current_proc, pure_local_trees);
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
          /* Compute the index inside the tree of the first element. */
          const t8_gloidx_t first_el_of_first_tree
            = *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, first_puretree_of_current_proc);
          first_element_in_tree_index_of_current_proc = first_element_index_of_current_proc - first_el_of_first_tree;
        }
        if (send_end_message) {
          /* Compute the index inside the tree of the last element. */
          const t8_gloidx_t first_el_of_last_tree
            = *(t8_gloidx_t *) t8_sc_array_index_locidx (&first_element_tree, last_puretree_of_current_proc);
          last_element_in_tree_index_of_current_proc = last_element_index_of_current_proc - first_el_of_last_tree;
          t8_debugf ("[H] Computing last element as %li - %li = %li\n\n", last_element_index_of_current_proc,
                     first_el_of_last_tree, last_element_in_tree_index_of_current_proc);
          t8_debugf ("[H] Last pure tree = %i, first element in it %li\n", last_puretree_of_current_proc,
                     first_el_of_last_tree);
        }
      }

      /*
       *
       *  MPI Communication starts here
       * 
       */

      /* Post the start message, we send the first tree id of the process
       * and (if desired) the index of the first element in this tree. */
      if (send_start_message) {
        const t8_gloidx_t global_id_of_first_tree
          = proc_is_empty
              ? -1
              : first_puretree_of_current_proc + t8_cmesh_get_first_treeid (cmesh) + first_tree_shared_shift;

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
          /* Send the global id of the first tree of this process */
          message[0] = global_id_of_first_tree;
          /* The index in the tree is the index of the element minus the offset of the tree. */
          message[1] = first_element_in_tree_index_of_current_proc;
          t8_debugf ("[H]  first e index %li, first in tree %li\n", first_element_index_of_current_proc,
                     first_puretree_of_current_proc >= 0
                       ? *(t8_gloidx_t *) t8_sc_array_index_gloidx (&first_element_tree, first_puretree_of_current_proc)
                       : -1);
          mpiret = sc_MPI_Isend (message, num_entries, T8_MPI_GLOIDX, iproc, T8_MPI_CMESH_UNIFORM_BOUNDS_START, comm,
                                 request);
          SC_CHECK_MPI (mpiret);
          t8_debugf ("Sending start message (%li, %li) to %i (global num el %li)\n", message[0], message[1], iproc,
                     data.global_num_elements);

          T8_ASSERT (proc_is_empty || (0 <= message[0] && message[0] < t8_cmesh_get_num_trees (cmesh)));
          T8_ASSERT (proc_is_empty || (0 <= message[1] && message[1] < data.global_num_elements));
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
          expect_start_message = 0;
#ifdef T8_ENABLE_DEBUG
          num_received_start_messages++;
#endif

          t8_debugf ("[H] Copied first tree %li element %li to self\n", global_id_of_first_tree,
                     first_element_in_tree_index_of_current_proc);
        }
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
          mpiret
            = sc_MPI_Isend (message, num_entries, T8_MPI_GLOIDX, iproc, T8_MPI_CMESH_UNIFORM_BOUNDS_END, comm, request);
          SC_CHECK_MPI (mpiret);
          t8_debugf ("Sending end message (%li, %li) to %i\n", message[0], message[1], iproc);
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
          expect_end_message = 0;
          t8_debugf ("[H] Copied last tree %li element %li to self\n", global_id_of_last_tree,
                     last_element_in_tree_index_of_current_proc + 1);
#ifdef T8_ENABLE_DEBUG
          num_received_end_messages++;
#endif
        }
      } /* End sending of end message */
    }   /* End loop over processes */
  }     /* if (pure_local_trees > 0) */

  const t8_gloidx_t first_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (cmesh->mpirank, data.num_procs, data.global_num_elements);
  const t8_gloidx_t last_element_index_of_current_proc
        = t8_cmesh_get_first_element_of_process (cmesh->mpirank + 1, data.num_procs, data.global_num_elements) - 1;

  if(first_element_index_of_current_proc > last_element_index_of_current_proc){
    expect_start_message = 0;
    expect_end_message = 0;
    *first_local_tree = 0;
    *last_local_tree = -1;
    *first_tree_shared = 0;
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin = -1;
    }
    if (child_in_tree_end != NULL) {
      *child_in_tree_end = - 1;
    }
    num_received_end_messages++;
    num_received_start_messages++;
  }

  /* Post the receives. */
  if (expect_start_message) {
    t8_debugf("expect start message\n");
    const int num_entries = 2;
    t8_gloidx_t *message = T8_ALLOC (t8_gloidx_t, num_entries);
    mpiret = sc_MPI_Recv (message, num_entries, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE, T8_MPI_CMESH_UNIFORM_BOUNDS_START,
                          comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);

#ifdef T8_ENABLE_DEBUG
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
      t8_debugf ("[H] Received cit begin = %li\n", *child_in_tree_begin);
      T8_ASSERT (*child_in_tree_begin == -1
                 || (0 <= *child_in_tree_begin && *child_in_tree_begin < data.global_num_elements));
    }
    T8_FREE (message);
  } /* End receiving start message */
  if (expect_end_message) {
    const int num_entries = 2;
    t8_gloidx_t *message = T8_ALLOC (t8_gloidx_t, num_entries);
    mpiret = sc_MPI_Recv (message, num_entries, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE, T8_MPI_CMESH_UNIFORM_BOUNDS_END, comm,
                          sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
#ifdef T8_ENABLE_DEBUG
    num_received_end_messages++;
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

  mpiret = sc_MPI_Waitall (num_messages_sent, (sc_MPI_Request *) send_requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (num_received_start_messages == 1);
  T8_ASSERT (num_received_end_messages == 1);

  if (pure_local_trees > 0) {
    sc_array_reset (&first_element_tree);
    sc_array_reset (&offset_partition);
    sc_array_reset (&send_buffer);
    sc_array_reset (&send_requests);
  }

  t8_shmem_array_destroy (&offset_array);

  t8_debugf ("Done with t8_cmesh_uniform_bounds_hybrid.\n");
  return;
}
