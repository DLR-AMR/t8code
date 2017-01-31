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

#include <sc_statistics.h>
#include <t8_refcount.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest_vtk.h>
#include <t8_cmesh/t8_cmesh_offset.h>

void
t8_forest_init (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);

  forest = *pforest = T8_ALLOC_ZERO (t8_forest_struct_t, 1);
  t8_refcount_init (&forest->rc);

  /* sensible (hard error) defaults */
  forest->mpicomm = sc_MPI_COMM_NULL;
  forest->dimension = -1;
  forest->from_method = T8_FOREST_FROM_LAST;

  forest->mpisize = -1;
  forest->mpirank = -1;
  forest->first_local_tree = -1;
  forest->global_num_elements = -1;
  forest->set_adapt_recursive = -1;
}

int
t8_forest_is_initialized (t8_forest_t forest)
{
  if (!(forest != NULL && t8_refcount_is_active (&forest->rc) &&
        !forest->committed)) {
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* TODO: check conditions that must always hold after init and before commit */
  if (0) {
    return 0;
  }
#endif

  return 1;
}

int
t8_forest_is_committed (t8_forest_t forest)
{
  if (!(forest != NULL && t8_refcount_is_active (&forest->rc)
        && forest->committed)) {
    return 0;
  }
#ifdef T8_ENABLE_DEBUG
  /* TODO: check more conditions that must always hold after commit */
  if (0) {
    return 0;
  }
#endif
  return 1;
}

static void
t8_forest_set_mpicomm (t8_forest_t forest, sc_MPI_Comm mpicomm, int do_dup)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (mpicomm != sc_MPI_COMM_NULL);

  forest->mpicomm = mpicomm;
  forest->do_dup = do_dup;
}

/* TODO: Change forest mpi logic */
void
t8_forest_set_cmesh (t8_forest_t forest, t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 do_dup;
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (cmesh != NULL);

  if (forest->cmesh != NULL) {
    t8_cmesh_unref (&forest->cmesh);
  }
  if (cmesh != NULL) {
    t8_cmesh_ref (cmesh);
    T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  }
  forest->cmesh = cmesh;
  do_dup = 0;
  t8_forest_set_mpicomm (forest, comm, do_dup);
}

void
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_t * scheme)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (scheme != NULL);

  forest->scheme = scheme;
}

void
t8_forest_set_level (t8_forest_t forest, int level)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);

  T8_ASSERT (0 <= level);

  forest->set_level = level;
}

void
t8_forest_set_copy (t8_forest_t forest, const t8_forest_t set_from)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_COPY;
}

void
t8_forest_set_partition (t8_forest_t forest, const t8_forest_t set_from,
                         int set_for_coarsening)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_for_coarsening = set_for_coarsening;

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_PARTITION;
}

void
t8_forest_set_adapt (t8_forest_t forest, const t8_forest_t set_from,
                     t8_forest_adapt_t adapt_fn,
                     t8_forest_replace_t replace_fn, int recursive)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme == NULL);
  T8_ASSERT (forest->set_from == NULL);
  T8_ASSERT (forest->set_adapt_fn == NULL);
  T8_ASSERT (forest->set_adapt_recursive == -1);

  forest->set_adapt_fn = adapt_fn;
  forest->set_replace_fn = replace_fn;
  forest->set_adapt_recursive = recursive != 0;
  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_ADAPT;
}

void
t8_forest_set_user_data (t8_forest_t forest, void *data)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  forest->user_data = data;
}

void *
t8_forest_get_user_data (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  return forest->user_data;
}

void
t8_forest_comm_global_num_elements (t8_forest_t forest)
{
  int                 mpiret;
  t8_gloidx_t         local_num_el;
  t8_gloidx_t         global_num_el;

  local_num_el = (t8_gloidx_t) forest->local_num_elements;
  mpiret = sc_MPI_Allreduce (&local_num_el, &global_num_el, 1,
                             T8_MPI_GLOIDX, sc_MPI_SUM, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  forest->global_num_elements = global_num_el;
}

/* For each tree in a forest compute its first and last descendant */
static void
t8_forest_compute_desc (t8_forest_t forest)
{
  t8_locidx_t         itree_id, num_trees;
  t8_tree_t           itree;
  t8_eclass_scheme_t *ts;
  t8_element_t       *element;

  T8_ASSERT (forest != NULL);
  /* Iterate over all trees */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree_id = 0; itree_id < num_trees; itree_id++) {
    /* get a pointer to the tree */
    itree = t8_forest_get_tree (forest, itree_id);
    /* get a pointer to the first element of itree */
    element = (t8_element_t *) t8_sc_array_index_locidx (&itree->elements, 0);
    /* get the eclass scheme associated to tree */
    ts = forest->scheme->eclass_schemes[itree->eclass];
    /* get memory for the trees first descendant */
    t8_element_new (ts, 1, &itree->first_desc);
    /* calculate the first descendant of the first element */
    t8_element_first_descendant (ts, element, itree->first_desc);
    /* get a pointer to the last element of itree */
    element = (t8_element_t *)
      t8_sc_array_index_locidx (&itree->elements,
                                itree->elements.elem_count - 1);
    /* get memory for the trees first descendant */
    t8_element_new (ts, 1, &itree->last_desc);
    /* calculate the last descendant of the first element */
    t8_element_last_descendant (ts, element, itree->last_desc);
  }
}

/* Create the elements on this process given a uniform partition
 * of the coarse mesh. */
static void
t8_forest_populate (t8_forest_t forest)
{
  t8_gloidx_t         child_in_tree_begin;
  t8_gloidx_t         child_in_tree_end;
  t8_locidx_t         count_elements;
  t8_locidx_t         num_tree_elements;
  t8_locidx_t         num_local_trees;
  t8_gloidx_t         jt, first_ctree;
  t8_gloidx_t         start, end, et;
  t8_tree_t           tree;
  t8_element_t       *element, *element_succ;
  sc_array_t         *telements;
  t8_eclass_t         tree_class;
  t8_eclass_scheme_t *eclass_scheme;
  t8_gloidx_t         cmesh_first_tree, cmesh_last_tree;

  /* TODO: create trees and quadrants according to uniform refinement */
  t8_cmesh_uniform_bounds (forest->cmesh, forest->set_level,
                           &forest->first_local_tree, &child_in_tree_begin,
                           &forest->last_local_tree, &child_in_tree_end,
                           NULL);

  cmesh_first_tree = t8_cmesh_get_first_treeid (forest->cmesh);
  cmesh_last_tree = cmesh_first_tree +
    t8_cmesh_get_num_local_trees (forest->cmesh) - 1;
  SC_CHECK_ABORT (forest->first_local_tree >= cmesh_first_tree
                  && forest->last_local_tree <= cmesh_last_tree,
                  "cmesh partition does not match the planned forest partition");

  forest->global_num_elements = forest->local_num_elements = 0;
  /* create only the non-empty tree objects */
  if (forest->first_local_tree >= forest->last_local_tree
      && child_in_tree_begin >= child_in_tree_end) {
    /* This processor is empty
     * we still set the tree array to store 0 as the number of trees here */
    forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
    count_elements = 0;
  }
  else {
    /* for each tree, allocate elements */
    num_local_trees = forest->last_local_tree - forest->first_local_tree + 1;
    forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
    sc_array_resize (forest->trees, num_local_trees);
    first_ctree = t8_cmesh_get_first_treeid (forest->cmesh);
    for (jt = forest->first_local_tree, count_elements = 0;
         jt <= forest->last_local_tree; jt++) {
      tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees,
                                                   jt -
                                                   forest->first_local_tree);
      tree_class = tree->eclass = t8_cmesh_get_tree_class (forest->cmesh,
                                                           jt - first_ctree);
      tree->elements_offset = count_elements;
      eclass_scheme = forest->scheme->eclass_schemes[tree_class];
      T8_ASSERT (eclass_scheme != NULL);
      telements = &tree->elements;
      /* calculate first and last element on this tree */
      start = (jt == forest->first_local_tree) ? child_in_tree_begin : 0;
      end = (jt == forest->last_local_tree) ? child_in_tree_end :
        t8_eclass_count_leaf (tree_class, forest->set_level);
      num_tree_elements = end - start;
      T8_ASSERT (num_tree_elements > 0);
      /* Allocate elements for this processor. */
      sc_array_init_size (telements, t8_element_size (eclass_scheme),
                          num_tree_elements);
      element = (t8_element_t *) t8_sc_array_index_locidx (telements, 0);
      eclass_scheme->elem_set_linear_id (element, forest->set_level, start);
      count_elements++;
      for (et = start + 1; et < end; et++, count_elements++) {
        element_succ =
          (t8_element_t *) t8_sc_array_index_locidx (telements, et - start);
        eclass_scheme->elem_successor (element, element_succ,
                                       forest->set_level);
        /* TODO: process elements here */
        element = element_succ;
      }
    }
  }
  forest->local_num_elements = count_elements;
  /* TODO: if no tree has pyramid type we can optimize this to
   * global_num_elements = global_num_trees * 2^(dim*level)
   */
  t8_forest_comm_global_num_elements (forest);
  /* TODO: figure out global_first_position, global_first_quadrant without comm */
}

/* return nonzero if the first tree of a forest is shared with a smaller
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
static int
t8_forest_first_tree_shared (t8_forest_t forest)
{
  t8_tree_t           first_tree;
  t8_element_t       *first_desc, *first_element;
  t8_eclass_t         eclass;
  t8_eclass_scheme_t *ts;
  int                 ret;

  T8_ASSERT (forest != NULL);
  if (forest->trees == NULL
      || forest->first_local_tree > forest->last_local_tree) {
    /* This forest is empty and therefore the first tree is not shared */
    return 0;
  }
  /* Get a pointer to the first tree */
  first_tree = (t8_tree_t) sc_array_index (forest->trees, 0);
  /* Get the eclass scheme of the first tree */
  eclass = first_tree->eclass;
  /* Get the eclass scheme of the first tree */
  ts = forest->scheme->eclass_schemes[eclass];
  /* Calculate the first possible descendant of the first tree */
  /* we do this by first creating a level 0 child of the tree, then
   * calculating its first descendant */
  t8_element_new (ts, 1, &first_element);
  t8_element_set_linear_id (ts, first_element, 0, 0);
  t8_element_new (ts, 1, &first_desc);
  t8_element_first_descendant (ts, first_element, first_desc);
  /* We can now check whether the first possible descendant matches the
   * first local descendant */
  ret = t8_element_compare (ts, first_desc, first_tree->first_desc);
  t8_element_destroy (ts, 1, &first_element);
  t8_element_destroy (ts, 1, &first_desc);
  /* If the descendants are the same then ret is zero and we return false.
   * We return true otherwise */
  return ret;
}

/* Allocate memory for trees and set their values as in from.
 * For each tree allocate enough element memory to fit the elements of from.
 * If copy_elements is true, copy the elements of from into the element memory.
 */
static void
t8_forest_copy_trees (t8_forest_t forest, t8_forest_t from, int copy_elements)
{
  t8_tree_t           tree, fromtree;
  t8_gloidx_t         num_tree_elements;
  t8_locidx_t         jt, number_of_trees;
  t8_eclass_scheme_t *eclass_scheme;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (from != NULL);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (from->committed);

  number_of_trees = from->trees->elem_count;
  forest->trees =
    sc_array_new_size (sizeof (t8_tree_struct_t), number_of_trees);
  sc_array_copy (forest->trees, from->trees);
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt);
    fromtree = (t8_tree_t) t8_sc_array_index_locidx (from->trees, jt);
    tree->eclass = fromtree->eclass;
    eclass_scheme = forest->scheme->eclass_schemes[tree->eclass];
    num_tree_elements = fromtree->elements.elem_count;
    sc_array_init_size (&tree->elements, t8_element_size (eclass_scheme),
                        num_tree_elements);
    /* TODO: replace with t8_elem_copy (not existing yet), in order to
     * eventually copy additional pointer data stored in the elements? */
    if (copy_elements) {
      sc_array_copy (&tree->elements, &fromtree->elements);
      tree->elements_offset = fromtree->elements_offset;
    }
    else {
      sc_array_truncate (&tree->elements);
    }
  }
  forest->first_local_tree = from->first_local_tree;
  forest->last_local_tree = from->last_local_tree;
  if (copy_elements) {
    forest->local_num_elements = from->local_num_elements;
    forest->global_num_elements = from->global_num_elements;
  }
  else {
    forest->local_num_elements = 0;
    forest->global_num_elements = 0;
  }
}

void
t8_forest_commit (t8_forest_t forest)
{
  int                 mpiret;
  sc_MPI_Comm         comm_dup;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of commit */
    forest->profile->commit_runtime = sc_MPI_Wtime ();
  }

  if (forest->set_from == NULL) {
    T8_ASSERT (forest->mpicomm != sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh != NULL);
    T8_ASSERT (forest->scheme != NULL);
    T8_ASSERT (forest->from_method == T8_FOREST_FROM_LAST);

    /* dup communicator if requested */
    if (forest->do_dup) {
      mpiret = sc_MPI_Comm_dup (forest->mpicomm, &comm_dup);
      SC_CHECK_MPI (mpiret);
      forest->mpicomm = comm_dup;
    }

    /* Set mpirank and mpisize */
    mpiret = sc_MPI_Comm_size (forest->mpicomm, &forest->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (forest->mpicomm, &forest->mpirank);
    SC_CHECK_MPI (mpiret);
    /* populate a new forest with tree and quadrant objects */
    t8_forest_populate (forest);
    forest->global_num_trees = t8_cmesh_get_num_trees (forest->cmesh);
  }
  else {
    T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh == NULL);
    T8_ASSERT (forest->scheme == NULL);
    T8_ASSERT (!forest->do_dup);
    T8_ASSERT (forest->from_method >= T8_FOREST_FROM_FIRST &&
               forest->from_method < T8_FOREST_FROM_LAST);

    /* TODO: optimize all this when forest->set_from has reference count one */

    /* we must prevent the case that set_from frees the source communicator */
    if (!forest->set_from->do_dup) {
      forest->mpicomm = forest->set_from->mpicomm;
    }
    else {
      mpiret = sc_MPI_Comm_dup (forest->set_from->mpicomm, &forest->mpicomm);
      SC_CHECK_MPI (mpiret);
    }
    forest->do_dup = forest->set_from->do_dup;

    /* Set mpirank and mpisize */
    mpiret = sc_MPI_Comm_size (forest->mpicomm, &forest->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (forest->mpicomm, &forest->mpirank);
    SC_CHECK_MPI (mpiret);

    /* increase reference count of cmesh and scheme from the input forest */
    t8_cmesh_ref (forest->cmesh = forest->set_from->cmesh);
    t8_scheme_ref (forest->scheme = forest->set_from->scheme);
    /* set the dimension, cmesh and scheme from the old forest */
    forest->dimension = forest->set_from->dimension;
    forest->cmesh = forest->set_from->cmesh;
    forest->scheme = forest->set_from->scheme;
    forest->global_num_trees = forest->set_from->global_num_trees;

    /* TODO: currently we can only handle copy, adapt, and partition */
    /* T8_ASSERT (forest->from_method == T8_FOREST_FROM_COPY); */
    if (forest->from_method == T8_FOREST_FROM_ADAPT) {
      if (forest->set_adapt_fn != NULL) {
        t8_forest_copy_trees (forest, forest->set_from, 0);
        t8_forest_adapt (forest);
      }
    }
    else if (forest->from_method == T8_FOREST_FROM_PARTITION) {
      forest->global_num_elements = forest->set_from->global_num_elements;
      /* Initialize the trees array of the forest */
      forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
      /* partition the forest */
      t8_forest_partition (forest);
    }

    /* decrease reference count of input forest, possibly destroying it */
    t8_forest_unref (&forest->set_from);
  }
  /* Compute first and last descendant for each tree */
  t8_forest_compute_desc (forest);

  /* we do not need the set parameters anymore */
  forest->set_level = 0;
  forest->set_for_coarsening = 0;
  forest->set_from = NULL;
  forest->committed = 1;
  t8_debugf ("Committed forest with %li local elements and %lli "
             "global elements.\n\tTree range ist from %lli to %lli.\n",
             (long) forest->local_num_elements,
             (long long) forest->global_num_elements,
             (long long) forest->first_local_tree,
             (long long) forest->last_local_tree);
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of commit */
    forest->profile->commit_runtime = sc_MPI_Wtime () -
      forest->profile->commit_runtime;
  }
}

t8_locidx_t
t8_forest_get_num_element (t8_forest_t forest)
{
  return forest->local_num_elements;
}

/* Currently this function is not used */
#if 0
static t8_element_t *
t8_forest_get_first_element (t8_forest_t forest)
{
  t8_tree_t           tree;

  if (forest->trees == NULL || forest->trees->elem_count == 0) {
    return NULL;
  }
  tree = t8_forest_get_tree (forest, 0);
  return (t8_element_t *) sc_array_index (&tree->elements, 0);
}
#endif

/* Compute the offset array for a partition cmesh that should match the
 * forest's partition.
 */
static t8_shmem_array_t
t8_forest_compute_cmesh_offset (t8_forest_t forest, sc_MPI_Comm comm)
{
  t8_shmem_array_t    offset;
  t8_gloidx_t         local_offset;
  int                 first_tree_shared;

  /* initialize the shared memory array */
  t8_shmem_array_init (&offset, sizeof (t8_gloidx_t), forest->mpisize + 1,
                       comm);
  /* Compute whether our first local tree is shared with a smaller rank */
  first_tree_shared = t8_forest_first_tree_shared (forest);
  /* Calculate our entry in the offset array */
  local_offset = t8_offset_first_tree_to_entry (forest->first_local_tree,
                                                first_tree_shared);
  /* allgather the local entries of the offset array */
  t8_shmem_array_allgather (&local_offset, 1, T8_MPI_GLOIDX, offset, 1,
                            T8_MPI_GLOIDX);
  /* Set the last entry of the offset array to the global number of trees */
  t8_shmem_array_set_gloidx (offset, forest->mpisize,
                             forest->global_num_trees);
  return offset;
}

void
t8_forest_partition_cmesh (t8_forest_t forest, sc_MPI_Comm comm,
                           int set_profiling)
{
  t8_cmesh_t          cmesh_partition;

  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, forest->cmesh);
  /* set partition range of new cmesh according to forest trees */
  t8_cmesh_set_partition_offsets (cmesh_partition,
                                  t8_forest_compute_cmesh_offset (forest,
                                                                  comm));
  /* Set the profiling of the cmesh */
  t8_cmesh_set_profiling (cmesh_partition, set_profiling);
  /* Commit the new cmesh */
  t8_cmesh_commit (cmesh_partition, comm);
  /* unref the old one and set the new cmesh as the cmesh of the forest */
  t8_cmesh_unref (&forest->cmesh);
  forest->cmesh = cmesh_partition;
}

t8_locidx_t
t8_forest_get_num_local_trees (t8_forest_t forest)
{
  t8_locidx_t         num_trees;

  num_trees = forest->last_local_tree - forest->first_local_tree + 1;
  /* assert for possible overflow */
  T8_ASSERT ((t8_gloidx_t) num_trees == forest->last_local_tree
             - forest->first_local_tree + 1);
  if (num_trees < 0) {
    /* Set number of trees to zero if there are none */
    num_trees = 0;
  }
  return num_trees;
}

/* TODO: We use this function in forest_partition when the
 * forest is only partially committed. Thus, we cannot check whether the
 * forest is committed here. */
t8_tree_t
t8_forest_get_tree (t8_forest_t forest, t8_locidx_t ltree_id)
{
  T8_ASSERT (forest->trees != NULL);
  T8_ASSERT (0 <= ltree_id
             && ltree_id < (t8_locidx_t) forest->trees->elem_count);
  return (t8_tree_t) t8_sc_array_index_locidx (forest->trees, ltree_id);
}

t8_cmesh_t
t8_forest_get_cmesh (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return forest->cmesh;
}

t8_locidx_t
t8_forest_get_tree_element_count (t8_tree_t tree)
{
  t8_locidx_t         element_count;

  T8_ASSERT (tree != NULL);
  element_count = tree->elements.elem_count;
  T8_ASSERT ((size_t) element_count == tree->elements.elem_count);
  return element_count;
}

/* Return the global index of the first local element */
t8_gloidx_t
t8_forest_get_first_local_element_id (t8_forest_t forest)
{
  t8_gloidx_t         first_element, local_num_elements;
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Convert local_num_elements to t8_gloidx_t */
  local_num_elements = forest->local_num_elements;
  /* MPI Scan over local_num_elements lead the global index of the first
   * local element */
  sc_MPI_Scan (&local_num_elements, &first_element, 1, T8_MPI_GLOIDX,
               sc_MPI_SUM, forest->mpicomm);
  /* MPI_Scan is inklusive, thus it counts our own data.
   * Therefore, we have to subtract it again */
  first_element -= local_num_elements;

  return first_element;
}

t8_eclass_t
t8_forest_get_eclass (t8_forest_t forest, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return t8_forest_get_tree (forest, ltreeid)->eclass;
}

t8_locidx_t
t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest, t8_locidx_t ltreeid)
{
  t8_gloidx_t         cmesh_gfirst;

  T8_ASSERT (forest->cmesh != NULL);

  cmesh_gfirst = t8_cmesh_get_first_treeid (forest->cmesh);
  return forest->first_local_tree - cmesh_gfirst + ltreeid;
}

void
t8_forest_set_profiling (t8_forest_t forest, int set_profiling)
{
  T8_ASSERT (t8_forest_is_initialized (forest));

  if (set_profiling) {
    if (forest->profile == NULL) {
      /* Only do something if profiling is not enabled already */
      forest->profile = T8_ALLOC_ZERO (t8_profile_struct_t, 1);
    }
  }
  else {
    /* Free any profile that is already set */
    if (forest->profile != NULL) {
      T8_FREE (forest->profile);
    }
  }
}

void
t8_forest_print_profile (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    /* Only print something if profiling is enabled */
    sc_statinfo_t       stats[T8_PROFILE_NUM_STATS];
    t8_profile_t       *profile = forest->profile;

    /* Set the stats */
    sc_stats_set1 (&stats[0], profile->partition_elements_shipped,
                   "forest: Number of elements sent.");
    sc_stats_set1 (&stats[1], profile->partition_elements_recv,
                   "forest: Number of elements received.");
    sc_stats_set1 (&stats[2], profile->partition_bytes_sent,
                   "forest: Number of bytes sent.");
    sc_stats_set1 (&stats[3], profile->partition_procs_sent,
                   "forest: Number of processes sent to.");
    sc_stats_set1 (&stats[4], profile->partition_runtime,
                   "forest: Partition runtime.");
    sc_stats_set1 (&stats[5], profile->commit_runtime,
                   "forest: Commit runtime.");
    /* compute stats */
    sc_stats_compute (sc_MPI_COMM_WORLD, T8_PROFILE_NUM_STATS, stats);
    /* print stats */
    t8_logf (SC_LC_GLOBAL, SC_LP_STATISTICS, "Printing stats for forest.\n");
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS,
                    T8_PROFILE_NUM_STATS, stats, 1, 1);
  }
}

void
t8_forest_write_vtk (t8_forest_t forest, const char *filename)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);

  t8_forest_vtk_write_file (forest, filename, 1);
}

/* Given function values at the four edge points of a unit square and
 * a point within that square, interpolate the function value at this point.
 * \param [in]    vertex  An array of size at least dim given the coordinates of the vertex to interpolate
 * \param [in]    corner_values An array of size 2^dim * 3, giving for each corner (in zorder) of
 *                        the unit square/cube its function values in 3D space.
 * \param [out]   evaluated_function An array of size 3, on output the function values
 *                        at \a vertex are stored here.
 */
static void
t8_forest_bilinear_interpolation (const double *vertex,
                                  const double *corner_values,
                                  int dim, double *evaluated_function)
{
  int                 i;
  double              temp[3] = { 0 };

  for (i = 0; i < dim; i++) {
    temp[i] = corner_values[0 * 3 + i] * (1 - vertex[0]) * (1 - vertex[1])      /* x=0 y=0 */
      +corner_values[1 * 3 + i] * vertex[0] * (1 - vertex[1])   /* x=1 y=0 */
      +corner_values[2 * 3 + i] * (1 - vertex[0]) * vertex[1]   /* x=0 y=1 */
      +corner_values[3 * 3 + i] * vertex[0] * vertex[1];        /* x=1 y=1 */
    if (dim == 3) {
      temp[i] *= (1 - vertex[2]);
      temp[i] += (corner_values[4 * 3 + i] * (1 - vertex[0]) * (1 - vertex[1])  /* x=0 y=0 z=1 */
                  +corner_values[5 * 3 + i] * vertex[0] * (1 - vertex[1])       /* x=1 y=0 z=1 */
                  +corner_values[6 * 3 + i] * (1 - vertex[0]) * vertex[1]       /* x=0 y=1 z=1 */
                  +corner_values[7 * 3 + i] * vertex[0] * vertex[1])    /* x=1 y=1 z=1 */
        *vertex[2];
    }
    evaluated_function[i] = temp[i];
  }
}

/* given an element in a coarse tree, the corner coordinates of the coarse tree
 * and a corner number of the element compute the coordinates of that corner
 * within the coarse tree.
 */
void
t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id,
                              t8_element_t * element, const double *vertices,
                              int corner_number, double *coordinates)
{
  int                 corner_coords[3], i;
  double              vertex_coords[3];
  t8_scheme_t        *scheme;
  t8_eclass_t         eclass;
  double              len;
  int                 dim;

  T8_ASSERT (forest != NULL);
  scheme = forest->scheme;
  T8_ASSERT (scheme != NULL);
  eclass = t8_forest_get_tree (forest, ltree_id)->eclass;
  T8_ASSERT (eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET
             || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX);

  dim = t8_eclass_to_dimension[eclass];
  len = 1. / t8_element_root_len (scheme->eclass_schemes[eclass], element);
  t8_element_vertex_coords (scheme->eclass_schemes[eclass], element,
                            corner_number, corner_coords);
  switch (eclass) {
  case T8_ECLASS_TRIANGLE:
    corner_coords[2] = 0;
  case T8_ECLASS_TET:
    for (i = 0; i < 3; i++) {
      coordinates[i] =
        len * (vertices[3 + i] - vertices[i]) * corner_coords[0] +
        (dim ==
         3 ? len * (vertices[9 + i] -
                    vertices[6 + i]) * corner_coords[1] : 0.) +
        len * (vertices[6 + i] - vertices[3 + i]) * corner_coords[dim - 1]
        + vertices[i];
    }
    break;
  case T8_ECLASS_QUAD:
    corner_coords[2] = 0;
  case T8_ECLASS_HEX:
    /* Store the coordinates of the corner scaled to the unit square/cube */
    for (i = 0; i < 3; i++) {
      vertex_coords[i] = len * corner_coords[i];
    }
    t8_forest_bilinear_interpolation ((const double *) vertex_coords,
                                      vertices, dim, coordinates);
    break;
  default:
    SC_ABORT ("Forest coordinate computation is supported only for "
              "triangles/tets/quads/hexes.");
  }
  return;
}

/* Iterate through all the trees and free the element memory as well as
 * the tree memory.
 */
static void
t8_forest_free_trees (t8_forest_t forest)
{
  t8_tree_t           tree;
  t8_locidx_t         jt, number_of_trees;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->committed);

  number_of_trees = forest->trees->elem_count;
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt);
    sc_array_reset (&tree->elements);
  }
  sc_array_destroy (forest->trees);
}

/* Completely destroy a forest and unreference all structs that the
 * forest has taken ownership on.
 */
static void
t8_forest_reset (t8_forest_t * pforest)
{
  int                 mpiret;
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount == 0);

  if (!forest->committed) {
    if (forest->set_from != NULL) {
      /* in this case we have taken ownership and not released it yet */
      t8_forest_unref (&forest->set_from);
    }
  }
  else {
    T8_ASSERT (forest->set_from == NULL);
  }

  /* undup communicator if necessary */
  if (forest->committed) {
    if (forest->do_dup) {
      mpiret = sc_MPI_Comm_free (&forest->mpicomm);
      SC_CHECK_MPI (mpiret);
    }
    t8_forest_free_trees (forest);
  }

  /* we have taken ownership on calling t8_forest_set_* */
  if (forest->scheme != NULL) {
    t8_scheme_unref (&forest->scheme);
  }
  if (forest->cmesh != NULL) {
    t8_cmesh_unref (&forest->cmesh);
  }

  /* free the memory of the offset array */
  if (forest->element_offsets != NULL) {
    t8_shmem_array_destroy (&forest->element_offsets);
  }
  if (forest->profile != NULL) {
    T8_FREE (forest->profile);
  }
  T8_FREE (forest);
  *pforest = NULL;
}

void
t8_forest_ref (t8_forest_t forest)
{
  T8_ASSERT (forest != NULL);
  t8_refcount_ref (&forest->rc);
}

void
t8_forest_unref (t8_forest_t * pforest)
{
  t8_forest_t         forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest != NULL);

  if (t8_refcount_unref (&forest->rc)) {
    t8_forest_reset (pforest);
  }
}
