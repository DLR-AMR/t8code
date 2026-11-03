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

/** \file t8_forest_partition.cxx
 * Implements functions declared in \ref t8_forest_partition.h.
 */

#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition_for_coarsening.h>
#include <t8_forest/t8_forest_splitters.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_schemes/t8_scheme.hxx>
#include <algorithm>
#include <cmath>
#include <vector>

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

t8_forest_t
t8_forest_new_gather (const t8_forest_t forest_from, const int gather_rank)
{
  // Declare and initialize forest to be returned.
  t8_forest_t forest_gather;
  t8_forest_init (&forest_gather);

  // Set partition (includes also sanity checks)
  t8_forest_set_partition (forest_gather, forest_from, 0);

  // Determine global ID of the element that will be the first local one:
  // To gather all elements on gather_rank, the first local element is set to
  // zero for all processes up to gather_rank and to number of elements for the
  // remaining ones.
  const t8_gloidx_t first_local_element
    = (forest_from->mpirank <= gather_rank) ? 0 : forest_from->global_num_leaf_elements;

  // Set the partition offset accordingly.
  t8_forest_set_partition_offset (forest_gather, first_local_element);

  // Commit the forest.
  t8_forest_commit (forest_gather);

  // Return gathered forest.
  return forest_gather;
}

void
t8_forest_set_partition_offset (t8_forest_t forest, const t8_gloidx_t first_global_element)
{
  // Set flag indicating a manual partition.
  forest->set_partition_offset = 1;

  // Set global ID of first local element.
  forest->set_first_global_element = first_global_element;
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
  SC_CHECK_ABORT (t8_shmem_init (comm) > 0, "Error in shared memory setup. Could not partition forest.");
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
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  if (!splitters.has_elements) {
    /* This forest is empty, nothing to do */
    return;
  }
  /* The last descendant of the last element of this rank lies before the first element of the next rank */
  t8_locidx_t itree = t8_forest_get_num_local_trees (forest) - 1;
  while (t8_forest_get_tree_num_leaf_elements (forest, itree) < 1) {
    itree--;
    T8_ASSERT (itree > -1);
  }
  const t8_tree_t tree = t8_forest_get_tree (forest, itree);
  const t8_element_t *element_last
    = t8_forest_get_leaf_element_in_tree (forest, itree, t8_forest_get_tree_leaf_element_count (tree) - 1);
  T8_ASSERT (t8_forest_get_scheme (forest)->element_is_valid (tree->eclass, element_last));
  t8_forest_position first, last;
  t8_forest_element_positions (forest, t8_forest_global_tree_id (forest, itree), tree->eclass, element_last, &first,
                               &last);
  T8_ASSERT (splitters.first_position <= first && last < splitters.next_position);
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
    SC_CHECK_ABORT (t8_shmem_init (comm) > 0, "Error in shared memory setup. Could not partition forest.");
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
    SC_CHECK_ABORT (t8_shmem_init (comm) > 0, "Error in shared memory setup. Could not partition forest.");
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

void
t8_forest_set_partition_weight_function (t8_forest_t forest, t8_weight_fcn_t *weight_fcn)
{
  forest->weight_function = weight_fcn;
}

/** The range of global element indices [first, next) of a rank in the new partition. */
typedef struct
{
  int rank;          /**< The rank. */
  t8_gloidx_t first; /**< The global index of the first element of the rank. */
  t8_gloidx_t next;  /**< The global index of the first element after the rank. */
} t8_forest_partition_range_t;

/** The global index of the first element of a rank in the uniform partition, in which every rank gets the same number
 * of elements up to one.
 * \param [in] forest_from  A committed forest.
 * \param [in] rank         A rank, or the number of ranks for the end of the last rank.
 * \return The global index of the first element of \a rank, the global number of elements for the end.
 */
static t8_gloidx_t
t8_forest_partition_uniform_first (const t8_forest_t forest_from, const int rank)
{
  const t8_gloidx_t global_num_leaf_elements = forest_from->global_num_leaf_elements;
  if (rank == forest_from->mpisize) {
    return global_num_leaf_elements;
  }
  return std::floor (static_cast<double> (global_num_leaf_elements) * rank / forest_from->mpisize);
}

/** The global index of the first element of this rank in the new partition, before the correction for coarsening.
 * Without a weight function, it is computed (or the offset set with \ref t8_forest_set_partition_offset). With a
 * weight function, the rank owning the element in the old partition computes it and sends it here.
 * \param [in] forest  The forest to partition, with forest->set_from committed.
 * \return The global index of the first element of this rank.
 * \note This function is collective.
 */
static t8_gloidx_t
t8_forest_partition_new_first (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->set_from != NULL);

  // Preparation: Define frequently used forest properties as own local variables.
  t8_forest_t forest_from = forest->set_from;
  sc_MPI_Comm comm = forest_from->mpicomm;
  const int mpirank = forest_from->mpirank;
  const int mpisize = forest_from->mpisize;
  const t8_gloidx_t global_num_leaf_elements = forest_from->global_num_leaf_elements;
  t8_weight_fcn_t *weight_fcn = forest->weight_function;

  // If the global number of elements is zero, all element offsets are too.
  if (global_num_leaf_elements == 0) {
    return 0;
  }
  // A manually set partition.
  if (forest->set_partition_offset != 0) {
    return forest->set_first_global_element;
  }
  // Without weights, the offsets are computed solely based on the number of elements.
  if (weight_fcn == nullptr) {
    return t8_forest_partition_uniform_first (forest_from, mpirank);
  }

  // Else: Weighted load balancing:
  // ------------------------------
  // Rank p > 0 starts at the first element whose prefix sum of the weights exceeds p / mpisize of the total weight,
  // and the owner of this element in the old partition computes it. For the ranks to agree on that owner, the end of
  // the prefix sums of a rank must equal the start of the next one, which MPI_Exscan computes in another order. So the
  // weights are rounded to multiples of a power of two (the unit) that is small against the total weight, which makes
  // all their sums exact.

  // The sum of the weights of the local elements
  const auto local_sum = [&] (const auto &weight) {
    double sum = 0.;
    for (t8_locidx_t ltreeid = 0; ltreeid < t8_forest_get_num_local_trees (forest_from); ++ltreeid) {
      for (t8_locidx_t ielm = 0; ielm < t8_forest_get_tree_num_leaf_elements (forest_from, ltreeid); ++ielm) {
        sum += weight (ltreeid, ielm);
      }
    }
    return sum;
  };
  const auto user_weight = [&] (const t8_locidx_t ltreeid, const t8_locidx_t ielm) {
    const double weight = weight_fcn (forest_from, ltreeid, ielm);
    T8_ASSERT (weight >= 0 && std::isfinite (weight));
    return weight;
  };

  // The largest sum of the weights of a rank, which is the same on all ranks, since the maximum is exact.
  double local_user_weight = local_sum (user_weight);
  double max_user_weight;
  int mpiret = sc_MPI_Allreduce (&local_user_weight, &max_user_weight, 1, sc_MPI_DOUBLE, sc_MPI_MAX, comm);
  SC_CHECK_MPI (mpiret);
  if (!(max_user_weight > 0)) {
    // Without any weight, the elements are distributed evenly.
    return t8_forest_partition_uniform_first (forest_from, mpirank);
  }
  // The total weight is at most mpisize * max_user_weight < 2^exponent = 2^51 unit. The rounded weights differ by at
  // most unit / 2 from the weights, so for less than 2^52 elements they add up to less than 2^52 unit, and all their
  // sums are exact.
  int exponent;
  std::frexp (mpisize * max_user_weight, &exponent);
  const double unit = std::ldexp (1., exponent - 51);
  const auto rounded_weight = [&] (const t8_locidx_t ltreeid, const t8_locidx_t ielm) {
    return std::round (user_weight (ltreeid, ielm) / unit) * unit;
  };

  // The sum of the rounded weights of this rank, of the ranks before it, and of all ranks
  double partition_weight = local_sum (rounded_weight);
  double partition_weight_offset = 0;
  mpiret = sc_MPI_Exscan (&partition_weight, &partition_weight_offset, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  if (mpirank == 0) {
    // The result of MPI_Exscan is undefined on rank 0
    partition_weight_offset = 0;
  }
  double forest_weight;
  mpiret = sc_MPI_Allreduce (&partition_weight, &forest_weight, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  if (!(forest_weight > 0)) {
    return t8_forest_partition_uniform_first (forest_from, mpirank);
  }

  // The ideal start of the weight of a rank, computed in the same way on all ranks
  const auto ideal_start = [&] (const int rank) { return forest_weight * rank / mpisize; };
  // The first rank p > 0 whose ideal start is not below a weight, or mpisize
  const auto first_rank_from = [&] (const double weight) {
    int lo = 1, hi = mpisize;
    while (lo < hi) {
      const int mid = lo + (hi - lo) / 2;
      if (ideal_start (mid) < weight) {
        lo = mid + 1;
      }
      else {
        hi = mid;
      }
    }
    return lo;
  };
  // This rank computes the first elements of the ranks whose ideal start lies in the range of its prefix sums. The end
  // of this range is the start of the range of the next rank, since the sums are exact.
  const int rank_begin = first_rank_from (partition_weight_offset);
  const int rank_end = first_rank_from (partition_weight_offset + partition_weight);
  std::vector<t8_gloidx_t> local_offsets (rank_end - rank_begin, 0);

  // Prepare computation of local element offsets.
  double accumulated_weight = partition_weight_offset;
  t8_gloidx_t global_elm_idx = t8_forest_get_first_local_leaf_element_id (forest_from);
  int iproc = rank_begin;

  // Loop over local trees and their elements:
  // - Add element weight to accumulated sum.
  // - If accumulated weight is above ideal bound for proc iproc, set local offset and go to next process.
  for (t8_locidx_t ltreeid = 0; ltreeid < t8_forest_get_num_local_trees (forest_from); ++ltreeid) {
    for (t8_locidx_t ielm = 0; ielm < t8_forest_get_tree_num_leaf_elements (forest_from, ltreeid); ++ielm) {
      T8_ASSERT (0 <= global_elm_idx && global_elm_idx < global_num_leaf_elements);
      accumulated_weight += rounded_weight (ltreeid, ielm);
      // there may be empty partitions, hence while and not if
      while (iproc < rank_end && accumulated_weight > ideal_start (iproc)) {
        local_offsets[iproc - rank_begin] = global_elm_idx;
        ++iproc;
      }
      ++global_elm_idx;
    }
  }
  T8_ASSERT (accumulated_weight == partition_weight_offset + partition_weight);
  T8_ASSERT (iproc == rank_end);  // i.e. local_offsets has been filled properly

  // Send each offset to its rank, which receives exactly one. The first rank starts at the first element, even if the
  // first elements have no weight.
  std::vector<sc_MPI_Request> requests (local_offsets.size ());
  for (int ioffset = 0; ioffset < (int) local_offsets.size (); ioffset++) {
    mpiret = sc_MPI_Isend (&local_offsets[ioffset], 1, T8_MPI_GLOIDX, rank_begin + ioffset, T8_MPI_PARTITION_OFFSETS,
                           comm, &requests[ioffset]);
    SC_CHECK_MPI (mpiret);
  }
  t8_gloidx_t first = 0;
  if (mpirank > 0) {
    mpiret
      = sc_MPI_Recv (&first, 1, T8_MPI_GLOIDX, sc_MPI_ANY_SOURCE, T8_MPI_PARTITION_OFFSETS, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  return first;
}

/** The global index of the first element of the next rank in the new partition, which that rank sends here.
 * \param [in] forest  The forest to partition, with forest->set_from committed.
 * \param [in] first   The global index of the first element of this rank in the new partition.
 * \return The global index of the first element of the next rank, the global number of elements for the last rank.
 * \note This function is collective.
 */
static t8_gloidx_t
t8_forest_partition_new_next (t8_forest_t forest, t8_gloidx_t first)
{
  const t8_forest_t forest_from = forest->set_from;
  const int mpirank = forest_from->mpirank;
  sc_MPI_Request request = sc_MPI_REQUEST_NULL;
  int mpiret;
  if (mpirank > 0) {
    mpiret
      = sc_MPI_Isend (&first, 1, T8_MPI_GLOIDX, mpirank - 1, T8_MPI_PARTITION_NEXT, forest_from->mpicomm, &request);
    SC_CHECK_MPI (mpiret);
  }
  t8_gloidx_t next = forest_from->global_num_leaf_elements;
  if (mpirank < forest_from->mpisize - 1) {
    mpiret = sc_MPI_Recv (&next, 1, T8_MPI_GLOIDX, mpirank + 1, T8_MPI_PARTITION_NEXT, forest_from->mpicomm,
                          sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Waitall (1, &request, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (first <= next);
  return next;
}

/** The ranges of the new partition that intersect the range of this rank in the old partition, sorted by rank: each
 * rank routes its new range to the ranks owning its elements in the old partition (\ref t8_forest_route).
 * \param [in] forest_from  The forest with the old partition.
 * \param [in] first        The global index of the first element of this rank in the new partition.
 * \param [in] next         The global index of the first element after this rank in the new partition.
 * \return The ranks to send elements to, with their new ranges.
 * \note This function is collective.
 */
static std::vector<t8_forest_partition_range_t>
t8_forest_partition_send_plan (const t8_forest_t forest_from, const t8_gloidx_t first, const t8_gloidx_t next)
{
  std::vector<t8_forest_partition_range_t> ranges;
  if (first < next) {
    ranges.push_back ({ forest_from->mpirank, first, next });
  }
  t8_forest_route (forest_from, ranges,
                   [] (const t8_forest_partition_range_t &range, const t8_forest_bisection &level) {
                     return t8_forest_route_halves_elements (level, range.first, range.next - 1, false);
                   });
  std::sort (
    ranges.begin (), ranges.end (),
    [] (const t8_forest_partition_range_t &a, const t8_forest_partition_range_t &b) { return a.rank < b.rank; });
  return ranges;
}

/** The ranges of the uniform partition that intersect the range of this rank in the old partition, sorted by rank.
 * Their first elements are computed, so no communication is needed.
 * \param [in] forest_from  The forest with the old partition.
 * \return The ranks to send elements to, with their new ranges.
 */
static std::vector<t8_forest_partition_range_t>
t8_forest_partition_send_plan_uniform (const t8_forest_t forest_from)
{
  std::vector<t8_forest_partition_range_t> ranges;
  const t8_gloidx_t old_first = t8_forest_get_first_local_leaf_element_id (forest_from);
  const t8_gloidx_t old_next = old_first + forest_from->local_num_leaf_elements;
  if (old_first == old_next) {
    return ranges;
  }
  /* The new owner of the first old element is the last rank whose first element is not after it */
  int lower = 0, upper = forest_from->mpisize;
  while (upper - lower > 1) {
    const int mid = lower + (upper - lower) / 2;
    if (t8_forest_partition_uniform_first (forest_from, mid) <= old_first) {
      lower = mid;
    }
    else {
      upper = mid;
    }
  }
  for (int rank = lower;
       rank < forest_from->mpisize && t8_forest_partition_uniform_first (forest_from, rank) < old_next; rank++) {
    const t8_gloidx_t first = t8_forest_partition_uniform_first (forest_from, rank);
    const t8_gloidx_t next = t8_forest_partition_uniform_first (forest_from, rank + 1);
    if (first < next) {
      ranges.push_back ({ rank, first, next });
    }
  }
  return ranges;
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

/** Send the elements of this rank in the old partition, or their data, to their ranks in the new partition.
 * \param [in] forest          The new forest, with forest->set_from the old one.
 * \param [in] plan            The ranks to send to with their new ranges, sorted by rank.
 * \param [in] send_data       If true, the entries of \a data_in are sent instead of the elements.
 * \param [in] data_in         With \a send_data, one entry per element of forest->set_from.
 * \param [in,out] requests    The requests of the sends are appended.
 * \param [in,out] buffers     The send buffers are appended, to be freed after the sends completed.
 * \param [out] to_self        The buffer for this rank itself, or NULL.
 * \param [out] bytes_to_self  The size of the buffer for this rank.
 */
static void
t8_forest_partition_send (t8_forest_t forest, const std::vector<t8_forest_partition_range_t> &plan, const int send_data,
                          const sc_array_t *data_in, std::vector<sc_MPI_Request> &requests,
                          std::vector<char *> &buffers, char **to_self, size_t *bytes_to_self)
{
  const t8_forest_t forest_from = forest->set_from;
  const t8_gloidx_t old_first = t8_forest_get_first_local_leaf_element_id (forest_from);
  const t8_gloidx_t old_next = old_first + forest_from->local_num_leaf_elements;
  t8_locidx_t current_tree = 0;
  *to_self = NULL;
  *bytes_to_self = 0;
  for (const t8_forest_partition_range_t &target : plan) {
    /* The local indices of the first and last element of the target */
    const t8_locidx_t first_send = SC_MAX (target.first, old_first) - old_first;
    const t8_locidx_t last_send = SC_MIN (target.next, old_next) - 1 - old_first;
    T8_ASSERT (first_send <= last_send);
    char *buffer;
    int buffer_alloc;
    if (!send_data) {
      t8_forest_partition_fill_buffer (forest_from, &buffer, &buffer_alloc, &current_tree, first_send, last_send);
    }
    else {
      t8_forest_partition_fill_buffer_data (forest_from, &buffer, &buffer_alloc, first_send, last_send, data_in);
    }
    buffers.push_back (buffer);
    if (target.rank == forest->mpirank) {
      *to_self = buffer;
      *bytes_to_self = buffer_alloc;
      continue;
    }
    t8_debugf ("Post send of %li elements (%i bytes) to process %i\n", (long) (last_send - first_send + 1),
               buffer_alloc, target.rank);
    requests.push_back (sc_MPI_REQUEST_NULL);
    const int mpiret
      = sc_MPI_Isend (buffer, buffer_alloc, sc_MPI_BYTE, target.rank,
                      send_data ? T8_MPI_PARTITION_DATA : T8_MPI_PARTITION_FOREST, forest->mpicomm, &requests.back ());
    SC_CHECK_MPI (mpiret);
    if (!send_data && forest->profile != NULL) {
      /* If profiling is enabled we count the elements, processes and bytes that we send to other processes */
      forest->profile->partition_elements_shipped += last_send - first_send + 1;
      forest->profile->partition_procs_sent += 1;
      forest->profile->partition_bytes_sent += buffer_alloc;
    }
  }
}

/** The number of elements in a message of the partition.
 * \param [in] buffer     The message.
 * \param [in] bytes      The size of the message.
 * \param [in] data_size  If the message contains data, the size of an entry, otherwise 0.
 * \return The number of elements in the message.
 */
static t8_locidx_t
t8_forest_partition_message_num_elements (const char *buffer, const size_t bytes, const size_t data_size)
{
  if (data_size > 0) {
    T8_ASSERT (bytes % data_size == 0);
    return bytes / data_size;
  }
  /* The message starts with the number of trees, followed by the tree info structs */
  const t8_locidx_t num_trees = *(const t8_locidx_t *) buffer;
  const t8_forest_partition_tree_info_t *tree_info
    = (const t8_forest_partition_tree_info_t *) (buffer + sizeof (t8_locidx_t) + T8_ADD_PADDING (sizeof (t8_locidx_t)));
  t8_locidx_t num_elements = 0;
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    num_elements += tree_info[itree].num_elements;
  }
  return num_elements;
}

/** Insert the data of a message received in \ref t8_forest_partition_receive.
 * \param [in]  forest      The new forest.
 * \param [in]  buffer      The message.
 * \param [in]  bytes       The size of the message.
 * \param [in,out] last_loc_elem_recvd On input the number of elements whose data was received before. Updated on
 *                          output.
 * \param [out] data_out    The received data.
 * The messages must be inserted in the order of the sending ranks.
 */
static void
t8_forest_partition_insert_message_data ([[maybe_unused]] t8_forest_t forest, const char *buffer, const size_t bytes,
                                         t8_locidx_t *last_loc_elem_recvd, sc_array_t *data_out)
{
  /* data_out must have the correct dimensions */
  T8_ASSERT (data_out != NULL);
  T8_ASSERT (data_out->elem_count == (size_t) forest->local_num_leaf_elements);
  T8_ASSERT (bytes % data_out->elem_size == 0);
  /* Copy the data behind the data received before */
  memcpy (data_out->array + data_out->elem_size * *last_loc_elem_recvd, buffer, bytes);
  *last_loc_elem_recvd += bytes / data_out->elem_size;
}

/** Insert the trees and elements of a message received in \ref t8_forest_partition_receive.
 * \param [in]  forest      The new forest.
 * \param [in]  proc        The rank that sent the message.
 * \param [in]  recv_buffer The message.
 * \param [in]  recv_bytes  The size of the message.
 * \param [in]  prev_recvd  The number of messages inserted before.
 * The messages must be inserted in the order of the sending ranks, to properly fill the forest->trees array.
 */
static void
t8_forest_partition_insert_message (t8_forest_t forest, const int proc, const char *recv_buffer,
                                    [[maybe_unused]] const size_t recv_bytes, const int prev_recvd)
{
  t8_locidx_t num_trees, itree;
  t8_locidx_t num_elements_recv;
  t8_locidx_t old_num_elements, new_num_elements;
  size_t tree_cursor, element_cursor;
  const t8_forest_partition_tree_info_t *tree_info;
  t8_tree_t tree, last_tree;
  size_t element_size {};
  const t8_scheme *scheme = t8_forest_get_scheme (forest->set_from);

  t8_debugf ("Inserting message of %zu bytes from process %i\n", recv_bytes, proc);
  /* Read the number of trees, it is the first locidx_t in recv_buffer */
  num_trees = *(const t8_locidx_t *) recv_buffer;
  /* Set the tree cursor to the first tree info entry in recv_buffer */
  tree_cursor = sizeof (t8_locidx_t) + T8_ADD_PADDING (sizeof (t8_locidx_t));
  /* Set the element cursor to the first element entry in recv_buffer */
  element_cursor = tree_cursor + num_trees * sizeof (t8_forest_partition_tree_info_t);
  /* Get the information for the first tree */
  tree_info = (const t8_forest_partition_tree_info_t *) (recv_buffer + tree_cursor);
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
      T8_ASSERT (element_cursor + tree_info->num_elements * element_size <= recv_bytes);
      t8_element_array_init_copy (&tree->leaf_elements, scheme, tree->eclass,
                                  (const t8_element_t *) (recv_buffer + element_cursor), tree_info->num_elements);
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

  if (forest->profile != NULL && proc != forest->mpirank) {
    /* If profiling is enabled we count the number of elements received from other processes */
    forest->profile->partition_elements_recv += num_elements_recv;
  }
}

/** Receive the elements of this rank in the new partition, or their data, and insert them in the order of the
 * sending ranks. These are the owners of the elements in the old partition. They are not known here, but the number
 * of elements is, so we receive from any rank until we have all elements. No message of another partition is
 * received: before sending the elements of the next partition, a rank commits a forest after this partition, and the
 * splitters computed in the commit need all ranks (the partition of data uses its own tag, and every data message
 * answers the range that this rank routed in \ref t8_forest_partition_send_plan).
 * \param [in] forest            The new forest.
 * \param [in] num_new_elements  The number of elements of this rank in the new partition.
 * \param [in] recv_data         If true, data is received instead of elements.
 * \param [out] data_out         With \a recv_data, one entry per element of \a forest.
 * \param [in] to_self           The message of this rank to itself, or NULL.
 * \param [in] bytes_to_self     The size of the message of this rank to itself.
 */
static void
t8_forest_partition_receive (t8_forest_t forest, const t8_locidx_t num_new_elements, const int recv_data,
                             sc_array_t *data_out, char *to_self, const size_t bytes_to_self)
{
  /* A received message */
  struct message
  {
    int proc;     /* The sender. */
    char *buffer; /* The message. */
    size_t bytes; /* The size of the message. */
  };
  const size_t data_size = recv_data ? data_out->elem_size : 0;
  std::vector<message> messages;
  t8_locidx_t num_received = 0;
  if (to_self != NULL) {
    messages.push_back ({ forest->mpirank, to_self, bytes_to_self });
    num_received += t8_forest_partition_message_num_elements (to_self, bytes_to_self, data_size);
  }
  while (num_received < num_new_elements) {
    sc_MPI_Status status;
    const int tag = recv_data ? T8_MPI_PARTITION_DATA : T8_MPI_PARTITION_FOREST;
    int mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, tag, forest->mpicomm, &status);
    SC_CHECK_MPI (mpiret);
    int recv_bytes;
    mpiret = sc_MPI_Get_count (&status, sc_MPI_BYTE, &recv_bytes);
    SC_CHECK_MPI (mpiret);
    char *buffer = T8_ALLOC (char, recv_bytes);
    mpiret
      = sc_MPI_Recv (buffer, recv_bytes, sc_MPI_BYTE, status.MPI_SOURCE, tag, forest->mpicomm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    t8_debugf ("Received message of %i bytes from process %i\n", recv_bytes, status.MPI_SOURCE);
    num_received += t8_forest_partition_message_num_elements (buffer, recv_bytes, data_size);
    messages.push_back ({ status.MPI_SOURCE, buffer, (size_t) recv_bytes });
  }
  T8_ASSERT (num_received == num_new_elements);

  /* The senders owned consecutive ranges of elements in the order of their ranks */
  std::sort (messages.begin (), messages.end (), [] (const message &a, const message &b) { return a.proc < b.proc; });
  if (!recv_data) {
    forest->local_num_leaf_elements = 0;
  }
  t8_locidx_t last_received_local_element = 0;
  for (size_t imessage = 0; imessage < messages.size (); imessage++) {
    const message &received = messages[imessage];
    if (!recv_data) {
      t8_forest_partition_insert_message (forest, received.proc, received.buffer, received.bytes, imessage);
    }
    else {
      t8_forest_partition_insert_message_data (forest, received.buffer, received.bytes, &last_received_local_element,
                                               data_out);
    }
    if (received.proc != forest->mpirank) {
      T8_FREE (received.buffer);
    }
  }
}

/** Partition a forest from forest->set_from, or data on it, given the range of this rank in the new partition.
 * \param [in,out] forest  The new forest, with forest->set_from the old one. Committed if \a send_data.
 * \param [in] first       The global index of the first element of this rank in the new partition.
 * \param [in] next        The global index of the first element after this rank in the new partition.
 * \param [in] send_data   If true, the entries of \a data_in are partitioned to \a data_out instead of the elements.
 * \param [in] data_in     With \a send_data, one entry per element of forest->set_from.
 * \param [out] data_out   With \a send_data, one entry per element of \a forest.
 * \param [in] uniform     If true, the new partition is the uniform one, whose ranges are computed instead of routed.
 * \note This function is collective.
 */
static void
t8_forest_partition_given (t8_forest_t forest, const t8_gloidx_t first, const t8_gloidx_t next, const int send_data,
                           const sc_array_t *data_in, sc_array_t *data_out, const bool uniform)
{
  t8_debugf ("Start partition_given\n");
  T8_ASSERT (send_data || t8_forest_is_initialized (forest));
  T8_ASSERT (!send_data || t8_forest_is_committed (forest));
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (t8_forest_is_committed (forest->set_from));
  T8_ASSERT (first <= next);

  /* Send the elements to their ranks in the new partition */
  const std::vector<t8_forest_partition_range_t> plan
    = uniform ? t8_forest_partition_send_plan_uniform (forest->set_from)
              : t8_forest_partition_send_plan (forest->set_from, first, next);
  std::vector<sc_MPI_Request> requests;
  std::vector<char *> buffers;
  char *to_self;
  size_t bytes_to_self;
  t8_forest_partition_send (forest, plan, send_data, data_in, requests, buffers, &to_self, &bytes_to_self);

  /* Receive the elements of this rank */
  const t8_locidx_t num_new_elements = next - first;
  T8_ASSERT (!send_data || num_new_elements == t8_forest_get_local_num_leaf_elements (forest));
  if (num_new_elements > 0) {
    t8_forest_partition_receive (forest, num_new_elements, send_data, data_out, to_self, bytes_to_self);
  }
  else if (!send_data) {
    /* This forest is empty, set first and last local tree such
     * that t8_forest_get_num_local_trees return 0 */
    forest->first_local_tree = 0;
    forest->last_local_tree = -1;
    forest->local_num_leaf_elements = 0;
  }
  /* Wait for all sends to complete */
  const int mpiret = sc_MPI_Waitall (requests.size (), requests.data (), sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  for (char *buffer : buffers) {
    T8_FREE (buffer);
  }
  t8_debugf ("Done partition_given\n");
}

/* Populate a forest with the partitioned elements of forest->set_from.
 * Without a weight function, the elements are distributed evenly (each element has the same weight).
 */
void
t8_forest_partition (t8_forest_t forest)
{
  t8_global_productionf ("Enter  forest partition.\n");
  t8_log_indent_push ();
  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (t8_forest_is_committed (forest->set_from));

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of partition */
    forest->profile->partition_runtime = sc_MPI_Wtime ();

    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start partition %f %f\n", sc_MPI_Wtime (), forest->profile->partition_runtime);
  }

  /* The range of this rank in the new partition. With partition for coarsening, its first element is moved behind a
   * family that would be split. The ranges of the uniform partition are known for all ranks. */
  const bool uniform
    = forest->weight_function == nullptr && forest->set_partition_offset == 0 && forest->set_for_coarsening == 0;
  t8_gloidx_t first = t8_forest_partition_new_first (forest);
  if (forest->set_for_coarsening != 0) {
    first = t8_forest_pfc_correct_first (forest, first);
  }
  const t8_gloidx_t next = uniform ? t8_forest_partition_uniform_first (forest->set_from, forest->set_from->mpirank + 1)
                                   : t8_forest_partition_new_next (forest, first);
  t8_forest_partition_given (forest, first, next, 0, NULL, NULL, uniform);

  T8_ASSERT ((size_t) t8_forest_get_num_local_trees (forest->set_from) == forest->set_from->trees->elem_count);
  T8_ASSERT ((size_t) t8_forest_get_num_local_trees (forest) == forest->trees->elem_count);

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

  /* perform the actual partitioning */
  const t8_gloidx_t first = t8_forest_get_first_local_leaf_element_id (forest_to);
  save_set_from = forest_to->set_from;
  forest_to->set_from = forest_from;
  t8_forest_partition_given (forest_to, first, first + forest_to->local_num_leaf_elements, 1, data_in, data_out, false);
  forest_to->set_from = save_set_from;

  t8_log_indent_pop ();
  t8_global_productionf ("Done forest partition data.\n");
}

T8_EXTERN_C_END ();
