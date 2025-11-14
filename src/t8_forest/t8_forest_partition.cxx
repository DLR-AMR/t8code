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

#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_schemes/t8_scheme.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/**
 * For each tree that we send elements from to other processes,
 * we send the information stored in this struct to the other process
 */
typedef struct
{
  t8_gloidx_t gtree_id;     /**< The global id of that tree. TODO: we could optimize this out */
  t8_eclass_t eclass;       /**< The element class of that tree. */
  t8_locidx_t num_elements; /**< The number of elements from this tree that were sent. */
} t8_forest_partition_tree_info_t;

/* Given the element offset array and a rank, return the first local element id of this rank */
static t8_gloidx_t
t8_forest_partition_first_element (const t8_gloidx_t *offset, int rank)
{
  return offset[rank];
}

/* Given the element offset array and a rank, return the last local element id of this rank */
static t8_gloidx_t
t8_forest_partition_last_element (const t8_gloidx_t *offset, int rank)
{
  return offset[rank + 1] - 1;
}

/* Query whether a given process is assigned no elements in an offset array */
static int
t8_forest_partition_empty (const t8_gloidx_t *offset, int rank)
{
  if (t8_forest_partition_first_element (offset, rank) >= t8_forest_partition_first_element (offset, rank + 1)) {
    return 1;
  }
  return 0;
}

/* Compute the global index of the first local element.
 * This function is collective. */
static t8_gloidx_t
t8_forest_compute_first_local_element_id (t8_forest_t forest)
{
  t8_gloidx_t first_element, local_num_elements;
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Convert local_num_leaf_elements to t8_gloidx_t */
  local_num_elements = forest->local_num_leaf_elements;
  /* MPI Scan over local_num_elements lead the global index of the first local element */
  sc_MPI_Scan (&local_num_elements, &first_element, 1, T8_MPI_GLOIDX, sc_MPI_SUM, forest->mpicomm);
  /* MPI_Scan is inklusive, thus it counts our own data.
   * Therefore, we have to subtract it again */
  first_element -= local_num_elements;

  return first_element;
}

/* For a committed forest create the array of element_offsets
 * and store it in forest->element_offsets
 */
void
t8_forest_partition_create_offsets (t8_forest_t forest)
{
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->forest_offsets_runtime = -sc_MPI_Wtime ();
  }

  sc_MPI_Comm comm;
  t8_gloidx_t first_local_element;

  T8_ASSERT (t8_forest_is_committed (forest));

  T8_ASSERT (forest->element_offsets == NULL);
  t8_debugf ("Building offsets for forest %p\n", (void *) forest);
  comm = forest->mpicomm;
  /* Set the shmem array type of comm */
  t8_shmem_init (comm);
  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  /* Initialize the offset array as a shmem array
   * holding mpisize+1 many t8_gloidx_t */
  t8_shmem_array_init (&forest->element_offsets, sizeof (t8_gloidx_t), forest->mpisize + 1, comm);
  /* Calculate the global index of the first local element */
  first_local_element = t8_forest_compute_first_local_element_id (forest);

  /* Collect all first global indices in the array */
  t8_shmem_array_allgather (&first_local_element, 1, T8_MPI_GLOIDX, forest->element_offsets, 1, T8_MPI_GLOIDX);
  if (t8_shmem_array_start_writing (forest->element_offsets)) {
    t8_shmem_array_set_gloidx (forest->element_offsets, forest->mpisize, forest->global_num_leaf_elements);
  }
  t8_shmem_array_end_writing (forest->element_offsets);
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->forest_offsets_runtime += sc_MPI_Wtime ();
  }
}

#if T8_ENABLE_DEBUG
/* Test if all first descendants of the elements in the first tree have
 * a greater or equal linear id than the stored first descendant. */
static void
t8_forest_partition_test_desc (t8_forest_t forest)
{
  t8_element_t *elem_desc;
  t8_linearidx_t first_desc_id;
  t8_locidx_t ielem;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_tree_t tree;
  int level;

  if (t8_forest_get_num_local_trees (forest) == 0) {
    /* This forest is empty, nothing to do */
    return;
  }

  tree = t8_forest_get_tree (forest, 0);
  const t8_eclass_t tree_class = tree->eclass;
  /* Get the first descendant id of this rank */
  first_desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, forest->mpirank);
  scheme->element_new (tree_class, 1, &elem_desc);
  for (ielem = 0; ielem < t8_forest_get_tree_leaf_element_count (tree); ielem++) {
    /* Iterate over elems, for each one create the first descendant and check
     * its linear id versus the linear id of first_desc. */
    const t8_element_t *element = t8_element_array_index_locidx (&tree->leaf_elements, ielem);
    scheme->element_get_first_descendant (tree_class, element, elem_desc, forest->maxlevel);
    level = scheme->element_get_level (tree_class, elem_desc);
    T8_ASSERT (level == scheme->element_get_level (tree_class, elem_desc));
    T8_ASSERT (level == forest->maxlevel);
    T8_ASSERT (scheme->element_get_linear_id (tree_class, elem_desc, level) >= first_desc_id);
  }
  scheme->element_destroy (tree_class, 1, &elem_desc);
}
#endif

void
t8_forest_partition_test_boundary_element ([[maybe_unused]] const t8_forest_t forest)
{
#if T8_ENABLE_DEBUG
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->global_first_desc != NULL);

  if (forest->mpisize == 1) {
    return;
  }
  /* Get the local number of elements of the first tree of process rank + 1.
   * If a rank contains no trees, its local_tree_num_elements is set to 0. */
  int mpirank_from;
  int mpirank_to;
  int mpiret;
  sc_MPI_Request request;
  sc_MPI_Status status;
  t8_locidx_t local_tree_num_elements;
  t8_locidx_t local_tree_num_elements_my;
  if (t8_forest_get_num_local_trees (forest) > 0) {
    local_tree_num_elements_my = t8_forest_get_tree_num_leaf_elements (forest, 0);
  }
  else {
    local_tree_num_elements_my = 0;
  }
  if (forest->mpirank == 0) {
    mpirank_from = forest->mpirank + 1;
    mpirank_to = forest->mpisize - 1;
  }
  else if (forest->mpirank == forest->mpisize - 1) {
    mpirank_from = 0;
    mpirank_to = forest->mpirank - 1;
  }
  else {
    mpirank_from = forest->mpirank + 1;
    mpirank_to = forest->mpirank - 1;
  }
  mpiret = sc_MPI_Irecv (&local_tree_num_elements, 1, T8_MPI_LOCIDX, mpirank_from, 0, forest->mpicomm, &request);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Send (&local_tree_num_elements_my, 1, T8_MPI_LOCIDX, mpirank_to, 0, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Wait (&request, &status);
  SC_CHECK_MPI (mpiret);

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  if (t8_forest_get_local_num_leaf_elements (forest) == 0) {
    /* This forest is empty, nothing to do */
    return;
  }
  if (forest->mpirank == forest->mpisize - 1) {
    /* The last process can only share a tree with process rank-1.
     * However, this is already tested by process rank-1. */
    return;
  }
  const t8_gloidx_t global_tree_id = t8_shmem_array_get_gloidx (forest->tree_offsets, forest->mpirank + 1);
  if (global_tree_id >= 0) {
    /* The first tree on process rank+1 is not shared with current rank, nothing to do */
    return;
  }
  /* The first tree on process rank+1 may be shared but empty.
   * Thus, the first descendant id of rank+1 is not of the first local tree. */
  if (local_tree_num_elements == 0) {
    /* check if first not shared tree of process rank+1 contains elements */
    return;
  }

  /* Get last element of current rank and its last descendant id */
  t8_locidx_t itree = num_local_trees - 1;
  while (t8_forest_get_tree_num_leaf_elements (forest, itree) < 1) {
    itree--;
    T8_ASSERT (itree > -1);
  }
  const t8_tree_t tree = t8_forest_get_tree (forest, itree);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = tree->eclass;
  t8_element_t *element_last_desc;
  scheme->element_new (tree_class, 1, &element_last_desc);
  /* last element of current rank */
  const t8_element_t *element_last
    = t8_forest_get_leaf_element_in_tree (forest, itree, t8_forest_get_tree_leaf_element_count (tree) - 1);
  T8_ASSERT (scheme->element_is_valid (tree_class, element_last));
  /* last and finest possiple element of current rank */
  scheme->element_get_last_descendant (tree_class, element_last, element_last_desc, forest->maxlevel);
  T8_ASSERT (scheme->element_is_valid (tree_class, element_last_desc));
  const int level = scheme->element_get_level (tree_class, element_last_desc);
  T8_ASSERT (level == scheme->element_get_level (tree_class, element_last_desc));
  T8_ASSERT (level == forest->maxlevel);
  const t8_linearidx_t last_desc_id = scheme->element_get_linear_id (tree_class, element_last_desc, level);
  /* Get the first descendant id of rank+1 */
  const t8_linearidx_t first_desc_id
    = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, forest->mpirank + 1);
  /* The following inequality must apply, if our last element is on the same tree :
   * last_desc_id of last element of rank < first_desc_id of first element of rank+1 */
  /** TODO: This assertion might still be wrong, when our last element is the last element of the tree*/
  T8_ASSERT (itree < num_local_trees - 1 || last_desc_id < first_desc_id);
  /* clean up */
  scheme->element_destroy (tree_class, 1, &element_last_desc);
#endif
}

void
t8_forest_partition_create_first_desc (t8_forest_t forest)
{
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->first_descendant_runtime = -sc_MPI_Wtime ();
  }
  sc_MPI_Comm comm;
  t8_linearidx_t local_first_desc;
  t8_element_t *first_desc = NULL;

  T8_ASSERT (t8_forest_is_committed (forest));

  T8_ASSERT (forest->global_first_desc == NULL);
  t8_debugf ("Building global first descendants for forest %p\n", (void *) forest);
  comm = forest->mpicomm;

  if (forest->global_first_desc == NULL) {
    /* Set the shmem array type of comm */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    /* Initialize the offset array as a shmem array
     * holding mpisize+1 many t8_linearidx_t to store the elements linear ids */
    t8_shmem_array_init (&forest->global_first_desc, sizeof (t8_linearidx_t), forest->mpisize, comm);
  }
  /* Assert whether the array is set up correctly */
  T8_ASSERT (t8_shmem_array_get_elem_count (forest->global_first_desc) == (size_t) forest->mpisize);
  T8_ASSERT (t8_shmem_array_get_elem_size (forest->global_first_desc) == sizeof (t8_linearidx_t));
  T8_ASSERT (t8_shmem_array_get_comm (forest->global_first_desc) == comm);
  if (forest->local_num_leaf_elements <= 0) {
    /* This process is empty, we store 0 in the array */
    local_first_desc = 0;
  }
  else {
    const t8_element_t *first_element = NULL;
    /* Get a pointer to the first local element. */
    if (forest->incomplete_trees) {
      for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
        if (t8_forest_get_tree_num_leaf_elements (forest, itree) > 0) {
          first_element = t8_forest_get_leaf_element_in_tree (forest, itree, 0);
          break;
        }
      }
    }
    else {
      first_element = t8_forest_get_leaf_element_in_tree (forest, 0, 0);
    }
    /* This process is not empty, the element was found, so we compute its first descendant. */
    if (first_element != NULL) {
      /* Get the eclass_scheme of the element. */
      const t8_scheme *scheme = t8_forest_get_scheme (forest);
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);
      scheme->element_new (tree_class, 1, &first_desc);
      scheme->element_get_first_descendant (tree_class, first_element, first_desc, forest->maxlevel);
      /* Compute the linear id of the descendant. */
      local_first_desc = scheme->element_get_linear_id (tree_class, first_desc, forest->maxlevel);
      scheme->element_destroy (tree_class, 1, &first_desc);
    }
  }
  /* Collect all first global indices in the array */
#if T8_ENABLE_DEBUG
#ifdef SC_ENABLE_MPI
  {
    /* We assert that we use the correct data size in the allgather call. */
    int mpiret, size;
    mpiret = MPI_Type_size (T8_MPI_LINEARIDX, &size);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (size == sizeof (t8_linearidx_t));
  }
#endif
#endif
  t8_shmem_array_allgather (&local_first_desc, 1, T8_MPI_LINEARIDX, forest->global_first_desc, 1, T8_MPI_LINEARIDX);
#if T8_ENABLE_DEBUG
  {
    int iproc;
    char buffer[BUFSIZ] = {};
    t8_linearidx_t desc_id;

    for (iproc = 0; iproc < forest->mpisize; iproc++) {
      desc_id = *(t8_linearidx_t *) t8_shmem_array_index (forest->global_first_desc, iproc);
      snprintf (buffer + strlen (buffer), BUFSIZ - strlen (buffer), " %llu,", (long long unsigned) desc_id);
    }
  }
  t8_forest_partition_test_desc (forest);
#endif
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->first_descendant_runtime += sc_MPI_Wtime ();
  }
}

void
t8_forest_partition_create_tree_offsets (t8_forest_t forest)
{

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->cmesh_offsets_runtime = -sc_MPI_Wtime ();
  }

  t8_gloidx_t tree_offset;
  sc_MPI_Comm comm;
  int is_empty, has_empty;

  T8_ASSERT (t8_forest_is_committed (forest));

  comm = forest->mpicomm;

  /* Calculate this process's tree offset */
  tree_offset = t8_forest_first_tree_shared (forest) ? -forest->first_local_tree - 1 : forest->first_local_tree;
  if (t8_forest_get_local_num_leaf_elements (forest) <= 0) {
    /* This forest is empty */
    is_empty = 1;
    /* Set the global number of trees as offset (temporarily) */
    tree_offset = forest->global_num_trees;
  }
  else {
    is_empty = 0;
  }

  if (forest->tree_offsets == NULL) {
    /* Set the shmem array type of comm */
    t8_shmem_init (comm);
    t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    /* Only allocate the shmem array, if it is not already allocated */
    t8_shmem_array_init (&forest->tree_offsets, sizeof (t8_gloidx_t), forest->mpisize + 1, comm);
  }
  /* Assert whether the array is set up correctly */
  T8_ASSERT (t8_shmem_array_get_elem_count (forest->tree_offsets) == (size_t) forest->mpisize + 1);
  T8_ASSERT (t8_shmem_array_get_elem_size (forest->tree_offsets) == sizeof (t8_gloidx_t));
  T8_ASSERT (t8_shmem_array_get_comm (forest->tree_offsets) == comm);
  /* gather all tree offsets from all processes */
  t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX, forest->tree_offsets, 1, T8_MPI_GLOIDX);

  /* Store the global number of trees at the entry mpisize in the array */
  if (t8_shmem_array_start_writing (forest->tree_offsets)) {
    t8_shmem_array_set_gloidx (forest->tree_offsets, forest->mpisize, forest->global_num_trees);
  }
  t8_shmem_array_end_writing (forest->tree_offsets);

  /* Communicate whether we have empty processes */
  sc_MPI_Allreduce (&is_empty, &has_empty, 1, sc_MPI_INT, sc_MPI_LOR, forest->mpicomm);

  if (has_empty) {
    int next_nonempty;
    /* there exist empty ranks, we have to recalculate the offset.
     * Each empty rank stores the offset of the next nonempty rank */
    if (is_empty) {
      const t8_gloidx_t *tree_offset_array = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
      /* Find the next rank that is not empty */
      next_nonempty = forest->mpirank + 1;
      while (next_nonempty < forest->mpisize && tree_offset_array[next_nonempty] >= forest->global_num_trees) {
        next_nonempty++;
      }
      /* Set the tree offset to the first nonshared tree of the next rank */
      tree_offset = t8_offset_first (next_nonempty, tree_offset_array);
      if (tree_offset_array[next_nonempty] < 0) {
        tree_offset++;
      }
    }
    /* Communicate the new tree offsets */
    t8_shmem_array_allgather (&tree_offset, 1, T8_MPI_GLOIDX, forest->tree_offsets, 1, T8_MPI_GLOIDX);
  }
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->cmesh_offsets_runtime += sc_MPI_Wtime ();
  }
}

// Compute forest->element_offsets according to the weight function, if provided
static void
t8_forest_partition_compute_new_offset (t8_forest_t forest, t8_weight_fcn_t *weight_fcn)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->element_offsets == NULL);

  t8_forest_t forest_from = forest->set_from;
  sc_MPI_Comm comm = forest_from->mpicomm;
  int const mpirank = forest_from->mpirank;
  int const mpisize = forest_from->mpisize;
  t8_gloidx_t const global_num_leaf_elements = forest_from->global_num_leaf_elements;

  /* Initialize the shmem array */
  t8_shmem_init (comm);
  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  t8_shmem_array_init (&forest->element_offsets, sizeof (t8_gloidx_t), mpisize + 1, comm);

  if (global_num_leaf_elements == 0) {
    if (t8_shmem_array_start_writing (forest->element_offsets)) {
      t8_gloidx_t *element_offsets = t8_shmem_array_get_gloidx_array_for_writing (forest->element_offsets);
      std::fill_n (element_offsets, mpisize + 1, t8_gloidx_t { 0 });
    }
    t8_shmem_array_end_writing (forest->element_offsets);
    return;
  }

  if (not weight_fcn) {
    if (t8_shmem_array_start_writing (forest->element_offsets)) {
      t8_gloidx_t *element_offsets = t8_shmem_array_get_gloidx_array_for_writing (forest->element_offsets);
      for (int i = 0; i < mpisize; ++i) {
        element_offsets[i] = std::floor (static_cast<double> (global_num_leaf_elements) * i / mpisize);
      }
      element_offsets[mpisize] = global_num_leaf_elements;
    }
    t8_shmem_array_end_writing (forest->element_offsets);
    return;
  }

  double const partition_weight = [&] () {  // sum of the weights on the local partition
    double local_sum = 0.;
    for (t8_locidx_t ltreeid = 0; ltreeid < t8_forest_get_num_local_trees (forest_from); ++ltreeid) {
      for (t8_locidx_t ielm = 0; ielm < t8_forest_get_tree_num_leaf_elements (forest_from, ltreeid); ++ielm) {
        local_sum += weight_fcn (forest_from, ltreeid, ielm);
      }
    }
    return local_sum;
  }();

  double const partition_weight_offset = [&] () {  // partial sum of the partition weights, excluding the local rank
    double local_offset = 0.;
    double local_partition_weight = partition_weight;  // because MPI does not like const variables
    sc_MPI_Exscan (&local_partition_weight, &local_offset, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
    return mpirank > 0 ? local_offset : 0;  // because the result of MPI_Exscan is undefined on rank 0
  }();

  double const forest_weight = [&] () {  // complete sum of the partition weights
    double total_weight = partition_weight_offset + partition_weight;
    sc_MPI_Bcast (&total_weight, 1, sc_MPI_DOUBLE, mpisize - 1, comm);
    return total_weight;
  }();

  // The [rank_begin, rank_end) slice of the new offsets land in the local partition
  int const rank_begin = std::ceil (mpisize * partition_weight_offset / forest_weight);
  int const rank_end = std::ceil (mpisize * (partition_weight_offset + partition_weight) / forest_weight);
  std::vector<t8_gloidx_t> local_offsets ( rank_end - rank_begin, 0 );

  double accumulated_weight = partition_weight_offset;
  t8_gloidx_t global_elm_idx = t8_forest_get_first_local_leaf_element_id (forest_from);
  int i = rank_begin;

  for (t8_locidx_t ltreeid = 0; ltreeid < t8_forest_get_num_local_trees (forest_from); ++ltreeid) {
    for (t8_locidx_t ielm = 0; ielm < t8_forest_get_tree_num_leaf_elements (forest_from, ltreeid); ++ielm) {
      T8_ASSERT (0 <= global_elm_idx && global_elm_idx < global_num_leaf_elements);
      accumulated_weight += weight_fcn (forest_from, ltreeid, ielm);
      while (accumulated_weight
             > forest_weight * i / mpisize) {  // there may be empty partitions, hence while and not if
        T8_ASSERT (rank_begin <= i && i < rank_end);
        local_offsets[i - rank_begin] = global_elm_idx;
        ++i;
      }
      ++global_elm_idx;
    }
  }
  T8_ASSERT (i == rank_end);  // i.e. local_offsets has been filled properly

  t8_shmem_array_allgatherv (local_offsets.data (), local_offsets.size (), T8_MPI_GLOIDX, forest->element_offsets,
                             T8_MPI_GLOIDX, comm);

  if (t8_shmem_array_start_writing (forest->element_offsets)) {
    t8_gloidx_t *element_offsets = t8_shmem_array_get_gloidx_array_for_writing (forest->element_offsets);
    element_offsets[0] = 0;
    element_offsets[mpisize] = global_num_leaf_elements;
  }
  t8_shmem_array_end_writing (forest->element_offsets);
}

/* Find the owner of a given element.
 */
static int
t8_forest_partition_owner_of_element (const int mpisize, const int mpirank, const t8_gloidx_t gelement,
                                      const t8_gloidx_t *offset)
{
  /* Tree offsets are stored similar enough that we can exploit their function */
  /* In the element offset logic, an element cannot be owned by more than one
   * process, thus any owner must be the unique owner. */
  return t8_offset_any_owner_of_tree_ext (mpisize, mpirank, gelement, offset);
}

/* Compute the first and last rank that we need to receive elements from */
static void
t8_forest_partition_recvrange (t8_forest_t forest, int *recv_first, int *recv_last)
{
  t8_gloidx_t first_element, last_element;

  /* Get the old element offset array */
  const t8_gloidx_t *offset_old = t8_shmem_array_get_gloidx_array (forest->set_from->element_offsets);
  /* Get the new element offset array */
  const t8_gloidx_t *offset_new = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  /* Compute new first and last element on this process from offset array */
  first_element = t8_forest_partition_first_element (offset_new, forest->mpirank);
  last_element = t8_forest_partition_last_element (offset_new, forest->mpirank);
  if (last_element < first_element) {
    /* There are no elements on this process and thus we cannot receive anything */
    *recv_first = 0;
    *recv_last = -1;
    return;
  }
  /* Calculate the first and last process we receive from */
  *recv_first = t8_forest_partition_owner_of_element (forest->mpisize, forest->mpirank, first_element, offset_old);
  *recv_last = t8_forest_partition_owner_of_element (forest->mpisize, forest->mpirank, last_element, offset_old);
}

/* Compute the first and last rank that we need to send elements to */
static void
t8_forest_partition_sendrange (t8_forest_t forest, int *send_first, int *send_last)
{
  t8_gloidx_t first_element, last_element;

  t8_debugf ("Calculate sendrange\n");
  if (forest->set_from->local_num_leaf_elements == 0) {
    /* There are no elements to send */
    *send_first = 0;
    *send_last = -1;
    return;
  }
  /* Get the old element offset array */
  const t8_gloidx_t *offset_old = t8_shmem_array_get_gloidx_array (forest->set_from->element_offsets);
  t8_debugf ("Partition forest from:\n");
  t8_offset_print (forest->set_from->element_offsets, forest->mpicomm);
  /* Get the new element offset array */
  const t8_gloidx_t *offset_new = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  t8_debugf ("Partition forest to:\n");
  t8_offset_print (forest->element_offsets, forest->mpicomm);
  /* Compute old first and last element on this process from offset array */
  first_element = t8_forest_partition_first_element (offset_old, forest->mpirank);
  last_element = t8_forest_partition_last_element (offset_old, forest->mpirank);
  /* Calculate the first and last process we send to */
  *send_first = t8_forest_partition_owner_of_element (forest->mpisize, forest->mpirank, first_element, offset_new);
  *send_last = t8_forest_partition_owner_of_element (forest->mpisize, forest->mpirank, last_element, offset_new);
}

/* Given a tree and its local id, the first and last element id that we need to send to a proc
 * and the first tree we send elements from,
 * calculate the first and last element of this tree that we need to send.
 * The returned element indices are local to the tree.
 * Returns true, if the last element that we send is also the last element of the tree
 */
static int
t8_forest_partition_tree_first_last_el (t8_tree_t tree, t8_locidx_t tree_id, t8_locidx_t first_element_send,
                                        t8_locidx_t last_element_send, t8_locidx_t current_tree,
                                        t8_locidx_t *first_tree_el, t8_locidx_t *last_tree_el)
{
  size_t num_elements;
  if (tree_id == current_tree) {
    /* For the first tree, the first element of that tree that we send is the
     * first element that we send */
    *first_tree_el = first_element_send - tree->elements_offset;
  }
  else {
    *first_tree_el = 0;
  }
  num_elements = t8_element_array_get_count (&tree->leaf_elements);
  if (tree->elements_offset + num_elements > (size_t) last_element_send) {
    /* For the last tree, the last element we send is the overall
     * last element that we send */
    *last_tree_el = last_element_send - tree->elements_offset;
  }
  else {
    /* Else it is the last element of this tree */
    *last_tree_el = num_elements - 1;
  }
  /* return true if the last element that we send is also the last
   * element of the tree */
  if (last_element_send - tree->elements_offset == (t8_locidx_t) num_elements - 1) {
    return 1;
  }
  return 0;
}

/* Fill the send buffers for one send operation.
 * \param [in]  forest_from     The original forest
 * \param [in]  send_buffer     Unallocated send_buffer
 * \param [out] buffer_alloc    The number of bytes in the send buffer
 * \param [in]  current_tree    On input the id of the first tree that we need to send
 *                              elements from. On output the id of the next tree that
 *                              we would send elements from to the next process.
 * \param [in]  first_element_send The local id of the first element that we need to send.
 * \param [in]  last_element_send The local id of the last element that we need to send.
 */
/* The send buffer will look like this:
 *
 * | number of trees | padding | tree_1 info | ... | tree_n info | tree_1 elements | ... | tree_n elements |
 */
/* If send_data is true, data must be an array of length forest_from->num_local elements
 * and instead of shipping the elements of forest_from, we ship the data entries. */
static void
t8_forest_partition_fill_buffer (t8_forest_t forest_from, char **send_buffer, int *buffer_alloc,
                                 t8_locidx_t *current_tree, t8_locidx_t first_element_send,
                                 t8_locidx_t last_element_send)
{
  t8_locidx_t num_elements_send;
  t8_tree_t tree;
  t8_locidx_t current_element, tree_id, num_trees_send;
  t8_locidx_t first_tree_element, last_tree_element;
  int element_alloc, byte_alloc, tree_info_pos, element_pos;
  int last_element_is_last_tree_element = 0;
  t8_forest_partition_tree_info_t *tree_info;
  t8_locidx_t *pnum_trees_send;
  size_t elem_size;

  current_element = first_element_send;
  tree_id = *current_tree;
  element_alloc = 0;
  num_trees_send = 0;
  /* At first we calculate the number of bytes that fit in the buffer */
  while (current_element <= last_element_send) {
    /* Get the first tree that we send elements from */
    tree = t8_forest_get_tree (forest_from, tree_id);
    /* clang-format off */
    last_element_is_last_tree_element = t8_forest_partition_tree_first_last_el (tree, tree_id, first_element_send,
                                                                                last_element_send, *current_tree,
                                                                                &first_tree_element,
                                                                                &last_tree_element);
    /* clang-format on */
    /* We now know how many elements this tree will send */
    num_elements_send = last_tree_element - first_tree_element + 1;
    T8_ASSERT (num_elements_send >= 0);
    elem_size = t8_element_array_get_size (&tree->leaf_elements);
    element_alloc += num_elements_send * elem_size;
    current_element += num_elements_send;
    num_trees_send++;
    tree_id++;
  }
  /* We calculate the total number of bytes that we need to allocate and allocate the buffer */
  /* The buffer consists of the number of trees, ... */
  byte_alloc = sizeof (t8_locidx_t);
  /* padding, ... */
  byte_alloc += T8_ADD_PADDING (byte_alloc);
  /* Store the position of the first tree info struct in the buffer */
  tree_info_pos = byte_alloc;
  /* an info struct for each tree, ... */
  byte_alloc += num_trees_send * sizeof (t8_forest_partition_tree_info_t);
  /* Store the position of the first element in the buffer */
  element_pos = byte_alloc;
  /* and the bytes for each tree's elements */
  byte_alloc += element_alloc;
  /* Note, that we do not add padding after the info structs and
   * each tree's elements, since these are multiples of structs and
   * structs are padded correctly */
  /* We allocate the buffer */
  *send_buffer = T8_ALLOC (char, byte_alloc);
  /* We store the number of trees at first in the send buffer */
  pnum_trees_send = (t8_locidx_t *) *send_buffer;
  *pnum_trees_send = num_trees_send;
  for (tree_id = 0; tree_id < num_trees_send; tree_id++) {
    /* Get the first tree that we send elements from */
    tree = t8_forest_get_tree (forest_from, tree_id + *current_tree);
    (void) t8_forest_partition_tree_first_last_el (tree, tree_id + *current_tree, first_element_send, last_element_send,
                                                   *current_tree, &first_tree_element, &last_tree_element);
    /* We now know how many elements this tree will send */

    num_elements_send = last_tree_element - first_tree_element + 1;

    T8_ASSERT (num_elements_send >= 0);
    /* Get the tree info struct for this tree and fill it */
    tree_info = (t8_forest_partition_tree_info_t *) (*send_buffer + tree_info_pos);
    tree_info->eclass = tree->eclass;
    tree_info->gtree_id = tree_id + *current_tree + forest_from->first_local_tree;
    tree_info->num_elements = num_elements_send;
    tree_info_pos += sizeof (t8_forest_partition_tree_info_t);
    /* We can now fill the send buffer with all elements of that tree */
    if (num_elements_send > 0) {
      const t8_element_t *pfirst_element = t8_element_array_index_locidx (&tree->leaf_elements, first_tree_element);
      elem_size = t8_element_array_get_size (&tree->leaf_elements);
      memcpy (*send_buffer + element_pos, (const void *) pfirst_element, num_elements_send * elem_size);
      element_pos += num_elements_send * elem_size;
    }
  }
  *current_tree += num_trees_send - 1 + last_element_is_last_tree_element;
  *buffer_alloc = byte_alloc;
  t8_debugf ("Post send of %i trees\n", num_trees_send);
}

/* Fill the send buffers for one send operation in send_data mode.
 * \param [in]  forest_from     The original forest
 * \param [in]  send_buffer     Unallocated send_buffer
 * \param [out] buffer_alloc    The number of bytes in the send buffer
 * \param [in]  first_element_send The local id of the first element that we need to send.
 * \param [in]  last_element_send The local id of the last element that we need to send.
 */
static void
t8_forest_partition_fill_buffer_data ([[maybe_unused]] t8_forest_t forest_from, char **send_buffer, int *buffer_alloc,
                                      t8_locidx_t first_element_send, t8_locidx_t last_element_send,
                                      const sc_array_t *data)
{
  void *data_entry;

  /* Check dimensions of data. */
  T8_ASSERT (data != NULL);
  T8_ASSERT (data->elem_count == (size_t) forest_from->local_num_leaf_elements);

  /* Calculate the byte count */
  *buffer_alloc = (last_element_send - first_element_send + 1) * data->elem_size;

  /* Allocate a multiple of the padding size capable of holding all bytes that will be sent. */
  const int internal_buffer_alloc = *buffer_alloc + T8_ADD_PADDING (*buffer_alloc);

  /* Must be a multiple of T8_PADDING_SIZE. */
  T8_ASSERT (T8_ADD_PADDING (internal_buffer_alloc) == 0);

  /* A negative amount of bytes ought to be send cannot exist and the amount of bytes that will be allocated
   * has to be at least equal to or greater than the amount of bytes that will be sent. */
  T8_ASSERT (*buffer_alloc >= 0 && *buffer_alloc <= internal_buffer_alloc);

  /* Allocate the send buffer. */
  *send_buffer = T8_ALLOC (char, internal_buffer_alloc);

  /* Copy the data to the send_buffer. */
  data_entry = t8_sc_array_index_locidx ((sc_array_t *) data, first_element_send);
  memcpy (*send_buffer, data_entry, *buffer_alloc);
}

/* Carry out all sending of elements */
/* If send_data is true, the elements are not send but element data
 * stored in an sc_array of length forest->set_from->num_local_elements.
 * Returns true if we sent to ourselves. */
static int
t8_forest_partition_sendloop (t8_forest_t forest, const int send_first, const int send_last, sc_MPI_Request **requests,
                              int *num_request_alloc, char ***send_buffer, const int send_data,
                              const sc_array_t *data_in, size_t *byte_to_self)
{
  int iproc, mpiret;
  t8_gloidx_t gfirst_element_send, glast_element_send;
  t8_gloidx_t gfirst_local_element;
  t8_locidx_t first_element_send, last_element_send;
  t8_locidx_t current_tree;
  t8_locidx_t num_elements_send;
  t8_forest_t forest_from;
  char **buffer;
  int buffer_alloc;
  sc_MPI_Comm comm;
  int to_self = 0;

  t8_debugf ("Start send loop\n");
  /* If send_data is false, the forest must not be committed but initialized.
   * If send_data is true, the forest must be committed */
  T8_ASSERT (send_data || t8_forest_is_initialized (forest));
  T8_ASSERT (!send_data || t8_forest_is_committed (forest));
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  /* If send data is true, data_in must be non-zero and of length num_local_elements */
  T8_ASSERT (!send_data || data_in != NULL);
  T8_ASSERT (!send_data || data_in->elem_count == (size_t) forest_from->local_num_leaf_elements);

  comm = forest->mpicomm;
  /* Determine the number of requests for MPI communication. */
  *num_request_alloc = send_last - send_first + 1;
  if (*num_request_alloc < 0) {
    /* If there are no processes to send to, this value could get negative */
    *num_request_alloc = 0;
    T8_ASSERT (send_last - send_first + 1 == 0);
  }
  *requests = T8_ALLOC (sc_MPI_Request, *num_request_alloc);

  /* Allocate memory for pointers to the send buffers */
  /* We allocate zero in order to set unused pointers to NULL so that we can pass them to free */
  *send_buffer = T8_ALLOC_ZERO (char *, send_last - send_first + 1);

  /* Get the new and old offset array */
  const t8_gloidx_t *offset_to = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const t8_gloidx_t *offset_from = t8_shmem_array_get_gloidx_array (forest_from->element_offsets);

  /* Compute the global id of the current first local element */
  gfirst_local_element = offset_from[forest->mpirank];
  /* loop over all processes that we send to */
  for (iproc = send_first; iproc <= send_last; iproc++) {
    /* At first, we compute the local index of the first and last element
     * that we send to proc */
    if (iproc == send_first) {
      /* If this is the first process we send to, the first element we send is
       * our very first element */
      first_element_send = 0;
      current_tree = 0;
    }
    else {
      /* Otherwise, the first element we send is the new first element on the process */
      gfirst_element_send = offset_to[iproc];
      first_element_send = gfirst_element_send - gfirst_local_element;
      /* assert for overflow error */
      T8_ASSERT ((t8_gloidx_t) first_element_send == gfirst_element_send - gfirst_local_element);
    }
    if (iproc == send_last) {
      /* To the last process we send all our remaining elements */
      last_element_send = forest_from->local_num_leaf_elements - 1;
    }
    else {
      /* Otherwise, the last element we send to proc is the last element on proc in the new partition. */
      glast_element_send = offset_to[iproc + 1] - 1;
      last_element_send = glast_element_send - gfirst_local_element;
    }
    num_elements_send = last_element_send - first_element_send + 1;
    if (num_elements_send < 0) {
      num_elements_send = 0;
    }
    /* We now know the local indices of the first and last element that we send to proc. */
    buffer = *send_buffer + iproc - send_first;
    if (num_elements_send > 0) {
      if (iproc == forest->mpirank) {
        to_self = 1;
      }
      if (!send_data) {
        /* Fill the buffer with the elements and calculate the next tree from which to send elements */
        t8_forest_partition_fill_buffer (forest_from, buffer, &buffer_alloc, &current_tree, first_element_send,
                                         last_element_send);
      }
      else {
        T8_ASSERT (send_data);
        /* We are in send data mode. Fill the send buffer with the data */
        t8_forest_partition_fill_buffer_data (forest_from, buffer, &buffer_alloc, first_element_send, last_element_send,
                                              data_in);
      }
      /* Post the MPI Send.
       * TODO: This will also send to ourselves if proc==mpirank */
      if (iproc != forest->mpirank) {
        t8_debugf ("Post send of %li elements (%i bytes) to process %i\n", (long) num_elements_send, buffer_alloc,
                   iproc);
        mpiret = sc_MPI_Isend (*buffer, buffer_alloc, sc_MPI_BYTE, iproc, T8_MPI_PARTITION_FOREST, comm,
                               *requests + iproc - send_first);
        SC_CHECK_MPI (mpiret);
      }
      else {
        *byte_to_self = buffer_alloc;
        *(*requests + iproc - send_first) = sc_MPI_REQUEST_NULL;
      }
      if (!send_data && forest->profile != NULL) {
        if (iproc != forest->mpirank) {
          /* If profiling is enabled we count the number of elements sent to other processes */
          forest->profile->partition_elements_shipped += num_elements_send;
          /* The number of procs we send to */
          forest->profile->partition_procs_sent += 1;
          /* The number of bytes that we send */
          forest->profile->partition_bytes_sent += buffer_alloc;
        }
      }
    }
    else {
      /* We do not send any elements to iproc (iproc is empty in new partition) */
      /* Set the request to NULL, such that it is ignored when we wait for the requests to complete */
      *(*requests + iproc - send_first) = sc_MPI_REQUEST_NULL;
    }
  }
  t8_debugf ("End send loop\n");
  return to_self;
}

/* Receive a message in data sending mode, send in sendloop.
 * \param [in]  forest      The new forest.
 * \param [in]  comm        The MPI communicator.
 * \param [in]  proc        The rank from which we receive.
 * \param [in]  status      MPI status with which we probed for the message.
 * \param [in,out] last_loc_elem_recv On input the local index of the last element
 *                          that was received by this rank. Updated on output.
 * \param [out] data_out    The received data.
 * \param [in]  sent_to_self If proc equals the rank of this process, the message
 *                          should be passed as this parameter.
 * \param [in]  byte_to_self If proc equals the rank of this process, the number of
 *                          bytes in the message.
 * It is important, that we receive the messages in order to properly fill the
 * data_out array.
 */
static void
t8_forest_partition_recv_message_data (t8_forest_t forest, sc_MPI_Comm comm, int proc, sc_MPI_Status *status,
                                       t8_locidx_t *last_loc_elem_recvd, sc_array_t *data_out, char *sent_to_self,
                                       size_t byte_to_self)
{
  int mpiret, recv_bytes;
  char *recv_buffer;
  size_t data_offset;

  /* data_out must have the correct dimensions */
  T8_ASSERT (data_out != NULL);
  T8_ASSERT (data_out->elem_count == (size_t) forest->local_num_leaf_elements);

  /* TODO: The next part is duplicated in t8_forest_partition_recv_message.
   *       Put duplicated code in function */
  /* further assertions */

  if (proc != forest->mpirank) {
    T8_ASSERT (proc == status->MPI_SOURCE);
    T8_ASSERT (status->MPI_TAG == T8_MPI_PARTITION_FOREST);

    /* Get the number of bytes to receive */
    mpiret = sc_MPI_Get_count (status, sc_MPI_BYTE, &recv_bytes);
    SC_CHECK_MPI (mpiret);
    t8_debugf ("Receiving message of %i bytes from process %i\n", recv_bytes, proc);
    /* allocate the receive buffer */
    recv_buffer = T8_ALLOC (char, recv_bytes);
    /* receive the message */
    mpiret
      = sc_MPI_Recv (recv_buffer, recv_bytes, sc_MPI_BYTE, proc, T8_MPI_PARTITION_FOREST, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  else {
    recv_buffer = sent_to_self;
    recv_bytes = byte_to_self;
  }

  /* Compute the place where to insert the data */
  data_offset = data_out->elem_size * *last_loc_elem_recvd;
  /* Copy the data */
  memcpy (data_out->array + data_offset, recv_buffer, recv_bytes);

  /* update the last element received */
  T8_ASSERT (recv_bytes % data_out->elem_size == 0);
  *last_loc_elem_recvd += recv_bytes / data_out->elem_size;

  if (proc != forest->mpirank) {
    /* free the receive buffer */
    T8_FREE (recv_buffer);
  }
}

/* Receive a message send in sendloop to this rank.
 * \param [in]  forest      The new forest.
 * \param [in]  comm        The MPI communicator.
 * \param [in]  proc        The rank from which we receive.
 * \param [in]  status      MPI status with which we probed for the message.
 * \param [in]  prev_recvd  The count of messages that we already received.
 * \param [in]  sent_to_self If proc equals the rank of this process, the message
 *                          should be passed as this parameter.
 * \param [in]  byte_to_self If proc equals the rank of this process, the number of
 *                          bytes in the message.
 * It is important, that we receive the messages in order to properly fill the forest->trees array.
 */
static void
t8_forest_partition_recv_message (t8_forest_t forest, sc_MPI_Comm comm, int proc, sc_MPI_Status *status, int prev_recvd,
                                  char *sent_to_self, size_t byte_to_self)
{
  int mpiret;
  int recv_bytes;
  char *recv_buffer;
  t8_locidx_t num_trees, itree;
  t8_locidx_t num_elements_recv;
  t8_locidx_t old_num_elements, new_num_elements;
  size_t tree_cursor, element_cursor;
  t8_forest_partition_tree_info_t *tree_info;
  t8_tree_t tree, last_tree;
  size_t element_size {};
  const t8_scheme *scheme = t8_forest_get_scheme (forest->set_from);

  if (proc != forest->mpirank) {
    T8_ASSERT (proc == status->MPI_SOURCE);
    T8_ASSERT (status->MPI_TAG == T8_MPI_PARTITION_FOREST);
    /* Get the number of bytes to receive */
    mpiret = sc_MPI_Get_count (status, sc_MPI_BYTE, &recv_bytes);
    SC_CHECK_MPI (mpiret);
  }
  else {
    recv_bytes = byte_to_self;
  }
  t8_debugf ("Receiving message of %i bytes from process %i\n", recv_bytes, proc);

  if (proc != forest->mpirank) {
    /* allocate the receive buffer */
    recv_buffer = T8_ALLOC (char, recv_bytes);
    /* clang-format off */
    /* receive the message */
    mpiret = sc_MPI_Recv (recv_buffer, recv_bytes, sc_MPI_BYTE, proc, T8_MPI_PARTITION_FOREST,
                          comm, sc_MPI_STATUS_IGNORE);
    /* clang-format on */
    SC_CHECK_MPI (mpiret);
  }
  else {
    recv_buffer = sent_to_self;
    recv_bytes = byte_to_self;
  }
  /* Read the number of trees, it is the first locidx_t in recv_buffer */
  num_trees = *(t8_locidx_t *) recv_buffer;
  /* Set the tree cursor to the first tree info entry in recv_buffer */
  tree_cursor = sizeof (t8_locidx_t) + T8_ADD_PADDING (sizeof (t8_locidx_t));
  /* Set the element cursor to the first element entry in recv_buffer */
  element_cursor = tree_cursor + num_trees * sizeof (t8_forest_partition_tree_info_t);
  /* Get the information for the first tree */
  tree_info = (t8_forest_partition_tree_info_t *) (recv_buffer + tree_cursor);
  if (prev_recvd == 0) {
    /* This is the first tree ever that we receive */
    /* We set the forests first local tree id */
    forest->first_local_tree = tree_info->gtree_id;
    /* In last_local_tree we keep track of the latest tree we received */
    forest->last_local_tree = tree_info->gtree_id - 1;
  }
  num_elements_recv = 0;
  for (itree = 0; itree < num_trees; itree++) {
    num_elements_recv += tree_info->num_elements;
    T8_ASSERT (tree_info->gtree_id >= forest->last_local_tree);
    if (tree_info->gtree_id > forest->last_local_tree) {
      /* We will insert a new tree in the forest */
      tree = (t8_tree_t) sc_array_push (forest->trees);
      tree->eclass = tree_info->eclass;
      /* Calculate the element offset of the new tree */
      if (forest->last_local_tree >= forest->first_local_tree) {
        /* If there is a previous tree, we read it */
        T8_ASSERT (forest->trees->elem_count >= 2); /* We added one tree and the current tree */
        last_tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, forest->trees->elem_count - 2);
        /* The element offset is the offset of the previous tree plus the number of elements in the previous tree */
        tree->elements_offset = last_tree->elements_offset + t8_forest_get_tree_leaf_element_count (last_tree);
      }
      else {
        /* This is the first tree, the element offset is thus zero */
        tree->elements_offset = 0;
      }
      /* Done calculating the element offset */
      /* Get the size of an element of the tree */
      element_size = scheme->get_element_size (tree->eclass);
      /* initialize the elements array and copy the elements from the receive buffer */
      T8_ASSERT (element_cursor + tree_info->num_elements * element_size <= (size_t) recv_bytes);
      t8_element_array_init_copy (&tree->leaf_elements, scheme, tree->eclass,
                                  (t8_element_t *) (recv_buffer + element_cursor), tree_info->num_elements);
    }
    else {
      T8_ASSERT (itree == 0); /* This situation only happens for the first tree */
      /* The tree is already present in the forest and we need to add elements to it */
      T8_ASSERT (forest->last_local_tree == tree_info->gtree_id);
      /* Get a pointer to the tree */
      tree = t8_forest_get_tree (forest, forest->last_local_tree - forest->first_local_tree);
      /* assert for correctness */
      T8_ASSERT (tree->eclass == tree_info->eclass);
      /* Get the old number of elements in the tree and calculate the new number */
      old_num_elements = t8_forest_get_tree_leaf_element_count (tree);
      new_num_elements = old_num_elements + tree_info->num_elements;
      /* Enlarge the elements array */
      t8_element_array_resize (&tree->leaf_elements, new_num_elements);
      if (tree_info->num_elements > 0) {
        t8_element_t *first_new_element
          = t8_element_array_index_locidx_mutable (&tree->leaf_elements, old_num_elements);
        /* Get the size of an element of the tree */
        element_size = scheme->get_element_size (tree->eclass);
        T8_ASSERT (element_size == t8_element_array_get_size (&tree->leaf_elements));
        /* Copy the elements from the receive buffer to the elements array */
        memcpy ((void *) first_new_element, recv_buffer + element_cursor, tree_info->num_elements * element_size);
      }
    }

    /* compute the new number of local elements */
    forest->local_num_leaf_elements += tree_info->num_elements;
    /* Set the new last local tree */
    forest->last_local_tree = tree_info->gtree_id;
    /* advance the element cursor */
    element_cursor += element_size * tree_info->num_elements;
    /* Advance to the next tree_info entry in the recv buffer */
    tree_cursor += sizeof (t8_forest_partition_tree_info_t);
    tree_info += 1;
  }

  if (proc != forest->mpirank) {
    T8_FREE (recv_buffer);
  }
  if (forest->profile != NULL) {
    if (proc != forest->mpirank) {
      /* If profiling is enabled we count the number of elements received from other processes */
      forest->profile->partition_elements_recv += num_elements_recv;
    }
  }
}

/* Receive the elements from all processes, we receive from.
 * The message are received in order of the sending rank,
 * since then we can easily build up the new trees array.
 */
static void
t8_forest_partition_recvloop (t8_forest_t forest, int recv_first, int recv_last, const int recv_data,
                              sc_array_t *data_out, char *sent_to_self, size_t byte_to_self)
{
  int iproc, prev_recvd;
  t8_locidx_t last_received_local_element = 0;
  t8_forest_t forest_from;
  int mpiret;
  sc_MPI_Comm comm;
  sc_MPI_Status status;

  /* Initial checks and inits */
  T8_ASSERT (recv_data || t8_forest_is_initialized (forest));
  T8_ASSERT (!recv_data || t8_forest_is_committed (forest));
  T8_ASSERT (!recv_data || data_out->elem_count == (size_t) forest->local_num_leaf_elements);
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  const t8_gloidx_t *offset_from = t8_shmem_array_get_gloidx_array (forest_from->element_offsets);
  comm = forest->mpicomm;

  /****     Actual communication    ****/

  /* In order of their ranks, receive the trees and elements from the other processes. */

  prev_recvd = 0;
  if (!recv_data) {
    forest->local_num_leaf_elements = 0;
  }
  for (iproc = recv_first; iproc <= recv_last; iproc++) {
    if (!t8_forest_partition_empty (offset_from, iproc)) {
      /* We receive from each nonempty rank between recv_first and recv_last */
      if (iproc != forest->mpirank) {
        /* Probe for the message */
        mpiret = sc_MPI_Probe (iproc, T8_MPI_PARTITION_FOREST, comm, &status);
        SC_CHECK_MPI (mpiret);
        /* Consistency checks */
        T8_ASSERT (iproc == status.MPI_SOURCE);
        T8_ASSERT (status.MPI_TAG == T8_MPI_PARTITION_FOREST);
      }
      /* Receive the actual message */
      if (!recv_data) {
        t8_forest_partition_recv_message (forest, comm, iproc, &status, prev_recvd, sent_to_self, byte_to_self);
      }
      else {
        T8_ASSERT (data_out != NULL);
        T8_ASSERT (data_out->elem_count == (size_t) forest->local_num_leaf_elements);
        t8_forest_partition_recv_message_data (forest, comm, iproc, &status, &last_received_local_element, data_out,
                                               sent_to_self, byte_to_self);
      }
      prev_recvd++;
    }
  }
}

/* Partition a forest from forest->set_from and the element offsets set in forest->element_offsets
 */
static void
t8_forest_partition_given (t8_forest_t forest, const int send_data, const sc_array_t *data_in, sc_array_t *data_out)
{
  int send_first, send_last, recv_first, recv_last;
  sc_MPI_Request *requests = NULL;
  int num_request_alloc; /* The count of elements in the request array */
  char **send_buffer, *sent_to_self;
  int mpiret, i, to_self;
  t8_locidx_t num_new_elements;
  size_t byte_to_self = 0;

  t8_debugf ("Start partition_given\n");
  T8_ASSERT (send_data || t8_forest_is_initialized (forest));
  T8_ASSERT (!send_data || t8_forest_is_committed (forest));
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (t8_forest_is_committed (forest->set_from));
  /* Compute the first and last rank that we send to */
  t8_forest_partition_sendrange (forest, &send_first, &send_last);
  t8_debugf ("send_first = %i\n", send_first);
  t8_debugf ("send_last = %i\n", send_last);

  /* Send all elements to other ranks */
  to_self = t8_forest_partition_sendloop (forest, send_first, send_last, &requests, &num_request_alloc, &send_buffer,
                                          send_data, data_in, &byte_to_self);
  if (to_self) {
    /* We have sent data to ourselves. */
    sent_to_self = *(send_buffer + forest->mpirank - send_first);
  }
  else {
    sent_to_self = NULL;
  }

  /* Compute the number of new elements on this forest */
  if (!send_data) {
    num_new_elements = t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank + 1)
                       - t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank);
  }
  else {
    num_new_elements = t8_forest_get_local_num_leaf_elements (forest);
  }

  if (num_new_elements > 0) {
    /* Receive all element from other ranks */
    t8_forest_partition_recvrange (forest, &recv_first, &recv_last);
    t8_forest_partition_recvloop (forest, recv_first, recv_last, send_data, data_out, sent_to_self, byte_to_self);
  }
  else if (!send_data) {
    /* This forest is empty, set first and last local tree such
     * that t8_forest_get_num_local_trees return 0 */
    forest->first_local_tree = 0;
    forest->last_local_tree = -1;
    forest->local_num_leaf_elements = 0;
  }
  /* Wait for all sends to complete */
  if (num_request_alloc > 0) {
    mpiret = sc_MPI_Waitall (num_request_alloc, requests, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  T8_FREE (requests);
  for (i = 0; i < num_request_alloc; i++) {
    T8_FREE (send_buffer[i]);
  }
  T8_FREE (send_buffer);

  t8_debugf ("Done partition_given\n");
}

/* Populate a forest with the partitioned elements of forest->set_from.
 * Currently the elements are distributed evenly (each element has the same weight).
 */
void
t8_forest_partition (t8_forest_t forest, t8_weight_fcn_t *weight_callback)
{
  t8_forest_t forest_from;
  int create_offset_from = 0;

  t8_global_productionf ("Enter  forest partition.\n");
  t8_log_indent_push ();
  T8_ASSERT (t8_forest_is_initialized (forest));
  forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->partition_runtime = sc_MPI_Wtime ();

    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start partition %f %f\n", sc_MPI_Wtime (), forest->profile->partition_runtime);
  }

  if (forest_from->element_offsets == NULL) {
    /* We create the partition table of forest_from */
    create_offset_from = 1;
    t8_forest_partition_create_offsets (forest_from);
  }
  /* TODO: if offsets already exist on forest_from, check it for consistency */

  /* We now calculate the new element offsets */
  t8_forest_partition_compute_new_offset (forest, weight_callback);
  t8_forest_partition_given (forest, 0, NULL, NULL);

  T8_ASSERT ((size_t) t8_forest_get_num_local_trees (forest_from) == forest_from->trees->elem_count);
  T8_ASSERT ((size_t) t8_forest_get_num_local_trees (forest) == forest->trees->elem_count);

  if (create_offset_from) {
    /* Delete the offset memory that we allocated */
    t8_shmem_array_destroy (&forest_from->element_offsets);
  }

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->partition_runtime = sc_MPI_Wtime () - forest->profile->partition_runtime;

    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End partition %f %f\n", sc_MPI_Wtime (), forest->profile->partition_runtime);
  }

  t8_log_indent_pop ();
  t8_global_productionf ("Done forest partition.\n");
}

void
t8_forest_partition_data (t8_forest_t forest_from, t8_forest_t forest_to, const sc_array_t *data_in,
                          sc_array_t *data_out)
{
  t8_forest_t save_set_from;

  t8_global_productionf ("Enter forest partition data.\n");
  t8_log_indent_push ();

  /* Assertions */
  T8_ASSERT (t8_forest_is_committed (forest_from));
  T8_ASSERT (t8_forest_is_committed (forest_to));
  T8_ASSERT (data_in != NULL && data_out != NULL);
  T8_ASSERT (data_in->elem_size == data_out->elem_size);

  /* data_in must have length of forest_from number of elements.
   * data_out length of forest_to number of elements */
  T8_ASSERT (data_in->elem_count == (size_t) forest_from->local_num_leaf_elements);
  T8_ASSERT (data_out->elem_count == (size_t) forest_to->local_num_leaf_elements);

  /* Create partition tables if not existent yet */
  if (forest_from->element_offsets == NULL) {
    /* We create the partition table of forest_from */
    t8_forest_partition_create_offsets (forest_from);
  }

  if (forest_to->element_offsets == NULL) {
    /* We create the partition table of forest_to */
    t8_forest_partition_create_offsets (forest_to);
  }

  /* perform the actual partitioning */
  save_set_from = forest_to->set_from;
  forest_to->set_from = forest_from;
  t8_forest_partition_given (forest_to, 1, data_in, data_out);
  forest_to->set_from = save_set_from;

  t8_log_indent_pop ();
  t8_global_productionf ("Done forest partition data.\n");
}

T8_EXTERN_C_END ();
