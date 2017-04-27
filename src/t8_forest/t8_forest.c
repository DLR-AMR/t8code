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
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest_vtk.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_trees.h>

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
    T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
  }
  forest->cmesh = cmesh;
  do_dup = 0;
  t8_forest_set_mpicomm (forest, comm, do_dup);
}

void
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_cxx_t * scheme)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->scheme_cxx == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (scheme != NULL);

  forest->scheme_cxx = scheme;
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
  T8_ASSERT (forest->scheme_cxx == NULL);
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
  T8_ASSERT (forest->scheme_cxx == NULL);
  T8_ASSERT (forest->set_from == NULL);

  T8_ASSERT (set_from != NULL);

  forest->set_for_coarsening = set_for_coarsening;

  forest->set_from = set_from;
  forest->from_method = T8_FOREST_FROM_PARTITION;
}

void
t8_forest_set_ghost (t8_forest_t forest, int do_ghost,
                     t8_ghost_type_t ghost_type)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  /* We currently only support face ghosts */
  SC_CHECK_ABORT (do_ghost == 0 || ghost_type == T8_GHOST_FACES,
                  "Ghost neighbors other than face-neighbors are not supported.\n");

  if (ghost_type == T8_GHOST_NONE) {
    /* none type disables ghost */
    forest->do_ghost = 0;
  }
  else {
    forest->do_ghost = (do_ghost != 0); /* True if and only if do_ghost != 0 */
  }
  forest->ghost_type = ghost_type;
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
  T8_ASSERT (forest->scheme_cxx == NULL);
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
  T8_ASSERT (!t8_forest_is_committed (forest));
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
    T8_ASSERT (forest->scheme_cxx != NULL);
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
  else {                        /* set_from != NULL */
    T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh == NULL);
    T8_ASSERT (forest->scheme_cxx == NULL);
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
    t8_scheme_cxx_ref (forest->scheme_cxx = forest->set_from->scheme_cxx);
    /* set the dimension, cmesh and scheme from the old forest */
    forest->dimension = forest->set_from->dimension;
    forest->cmesh = forest->set_from->cmesh;
    forest->scheme_cxx = forest->set_from->scheme_cxx;
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
  }                             /* end set_from != NULL */
  /* Compute the element offset of the trees */
  t8_forest_compute_elements_offset (forest);
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

  /* From here on, the forest passes the t8_forest_is_committed check */
  /* re-partition the cmesh */
  if (forest->cmesh->set_partition &&
      forest->from_method == T8_FOREST_FROM_PARTITION) {
    t8_forest_partition_cmesh (forest, forest->mpicomm,
                               forest->profile != NULL);
  }
  /* Construct a ghost layer, if desired */
  if (forest->do_ghost) {
    /* TODO: ghost type */
    t8_forest_ghost_create (forest);
  }

  forest->do_ghost = 0;
}

t8_locidx_t
t8_forest_get_num_element (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->local_num_elements;
}

t8_gloidx_t
t8_forest_get_global_num_elements (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->global_num_elements;
}

t8_locidx_t
t8_forest_get_num_ghosts (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Return the number of ghost elements, or 0 if no ghost structure
   * existst. */
  if (forest->ghosts == NULL) {
    return 0;
  }
  return forest->ghosts->num_ghosts_elements;
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
static              t8_shmem_array_t
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

  t8_debugf ("Partitioning cmesh according to forest\n");

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
  /* set the new cmesh as the cmesh of the forest */
  forest->cmesh = cmesh_partition;
  t8_debugf ("Done partitioning cmesh\n");
}

sc_MPI_Comm
t8_forest_get_mpicomm (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return forest->mpicomm;
}

t8_gloidx_t
t8_forest_get_first_local_tree_id (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->first_local_tree;
}

t8_locidx_t
t8_forest_get_num_ghost_trees (t8_forest_t forest)
{
  if (forest->ghosts != NULL) {
    return t8_forest_ghost_num_trees (forest);
  }
  else {
    return 0;
  }
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

t8_gloidx_t
t8_forest_get_num_global_trees (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->global_num_trees;
}

t8_gloidx_t
t8_forest_global_tree_id (t8_forest_t forest, t8_locidx_t ltreeid)
{
  t8_locidx_t         num_local_trees;
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest)
             + t8_forest_ghost_num_trees (forest));

  num_local_trees = t8_forest_get_num_local_trees (forest);
  if (ltreeid < num_local_trees) {
    /* The tree is a local tree */
    return ltreeid + forest->first_local_tree;
  }
  else {
    T8_ASSERT (forest->ghosts != NULL);
    /* Return the global id of the ghost tree */
    return t8_forest_ghost_get_global_treeid (forest,
                                              ltreeid - num_local_trees);
  }
}

/* TODO: We use this function in forest_partition when the
 * forest is only partially committed. Thus, we cannot check whether the
 * forest is committed here. */
t8_tree_t
t8_forest_get_tree (t8_forest_t forest, t8_locidx_t ltree_id)
{
  T8_ASSERT (forest->trees != NULL);
  T8_ASSERT (0 <= ltree_id
             && ltree_id < t8_forest_get_num_local_trees (forest));
  return (t8_tree_t) t8_sc_array_index_locidx (forest->trees, ltree_id);
}

t8_cmesh_t
t8_forest_get_cmesh (t8_forest_t forest)
{
  return forest->cmesh;
}

/* Compare function for the binary search in t8_forest_get_element.
 * Given a local element id and tree, this function returns 0
 * if the  element is inside the tree, -1 if it is inside a tree with
 * bigger local tree id and +1 if the element is inside a tree with
 * smaller local tree id.
 */
static int
t8_forest_compare_elem_tree (const void *lelement_id, const void *ltree)
{
  t8_locidx_t         leid = *(const t8_locidx_t *) lelement_id;
  const t8_tree_t     tree = (const t8_tree_t) ltree;

  if (tree->elements_offset > leid) {
    /* We have to look further to the left */
    return -1;
  }
  else if (tree->elements_offset + tree->elements.elem_count > leid) {
    /* We have found the tree */
    return 0;
  }
  else {
    /* We have to look further right */
    return 1;
  }
}

t8_element_t       *
t8_forest_get_element (t8_forest_t forest, t8_locidx_t lelement_id,
                       t8_locidx_t * ltreeid)
{
  t8_tree_t           tree;
  t8_locidx_t         ltree;
#ifdef T8_ENABLE_DEBUG
  t8_locidx_t         ltreedebug;
#endif

  T8_ASSERT (lelement_id >= 0);
  if (lelement_id >= t8_forest_get_num_element (forest)) {
    return NULL;
  }
  /* We optimized the binary search out by using sc_bsearch,
   * but keep it in for debugging. We check whether the hand-written
   * binary search matches the sc_array_bsearch. */
#ifdef T8_ENABLE_DEBUG
  {
    t8_locidx_t         ltree_a, ltree_b;
    ltree_a = 0;
    ltree_b = t8_forest_get_num_local_trees (forest);
    ltreedebug = (ltree_a + ltree_b) / 2;
    while (ltree_a < ltree_b) {
      ltreedebug = (ltree_a + ltree_b) / 2;
      tree = t8_forest_get_tree (forest, ltreedebug);
      if (tree->elements_offset > lelement_id) {
        /* We have to look further to the left */
        ltree_b = ltreedebug;
      }
      else if (tree->elements_offset + tree->elements.elem_count >
               lelement_id) {
        /* We have found the tree */
        ltree_a = ltree_b;
      }
      else {
        /* We have to look further right */
        ltree_a = ltreedebug;
      }
    }
  }
#endif
  ltree =
    sc_array_bsearch (forest->trees, &lelement_id,
                      t8_forest_compare_elem_tree);
  T8_ASSERT (ltreedebug == ltree);
  if (ltreeid != NULL) {
    *ltreeid = ltree;
  }

  /* The tree that contains the element is now local tree ltree.
   * Or the element is not a local element. */
  tree = t8_forest_get_tree (forest, ltree);
  if (tree->elements_offset <= lelement_id && lelement_id <
      tree->elements_offset + tree->elements.elem_count) {
    return (t8_element_t *)
      t8_sc_array_index_locidx (&tree->elements,
                                lelement_id - tree->elements_offset);
  }
  /* The element was not found.
   * This case is covered by the first if and should therefore
   * never happen. */
  SC_ABORT_NOT_REACHED ();
  return NULL;
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

t8_eclass_t
t8_forest_get_tree_class (t8_forest_t forest, t8_locidx_t ltreeid)
{
  t8_locidx_t         num_local_trees =
    t8_forest_get_num_local_trees (forest);

  T8_ASSERT (0 <= ltreeid
             && ltreeid <
             num_local_trees + t8_forest_get_num_ghost_trees (forest));
  if (ltreeid < num_local_trees) {
    /* The id belongs to a local tree */
    return t8_forest_get_tree (forest, ltreeid)->eclass;
  }
  else {
    /* The id belongs to a ghost tree */
    return t8_forest_ghost_get_tree_class (forest, ltreeid - num_local_trees);
  }
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

t8_scheme_cxx_t *
t8_forest_get_scheme (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->scheme_cxx != NULL);

  return forest->scheme_cxx;
}

t8_eclass_scheme_c *
t8_forest_get_eclass_scheme (t8_forest_t forest, t8_eclass_t eclass)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->scheme_cxx != NULL);
  T8_ASSERT (eclass != T8_ECLASS_COUNT);

  return forest->scheme_cxx->eclass_schemes[eclass];
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
  t8_locidx_t         num_local_trees;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->cmesh != NULL);
  num_local_trees = t8_forest_get_num_local_trees (forest);
  T8_ASSERT (0 <= ltreeid
             && ltreeid < num_local_trees
             + t8_forest_ghost_num_trees (forest));

  t8_debugf ("[H] lookup local id %i\n", ltreeid);
  if (ltreeid < num_local_trees) {
    /* This a local tree and not a ghost */

    cmesh_gfirst = t8_cmesh_get_first_treeid (forest->cmesh);
    return forest->first_local_tree - cmesh_gfirst + ltreeid;
  }
  else {
    /* This is a ghost */
    t8_gloidx_t         globalid;
    t8_locidx_t         cmesh_local_id;
    /* Compute the global id of this ghost tree */
    globalid = t8_forest_ghost_get_global_treeid (forest,
                                                  ltreeid - num_local_trees);
    /* Compute the cmesh local id of the ghost */
    cmesh_local_id = t8_cmesh_get_local_id (forest->cmesh, globalid);
    t8_debugf ("[H] is ghost with global id %li local id %i %p\n", globalid,
               cmesh_local_id,
               forest->cmesh->trees->ghost_globalid_to_local_id);
    /* is < 0 if this ghost does not exist */
    T8_ASSERT (cmesh_local_id >= 0);
    return cmesh_local_id;
  }
}

t8_locidx_t
t8_forest_cmesh_ltreeid_to_ltreeid (t8_forest_t forest, t8_locidx_t lctreeid)
{
  t8_locidx_t         ltreeid;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->cmesh != NULL);

  ltreeid = t8_cmesh_get_first_treeid (forest->cmesh) -
    t8_forest_get_first_local_tree_id (forest) + lctreeid;
  if (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest)) {
    /* The tree is a forest local tree */
    return ltreeid;
  }
  else {
    /* The tree is not forest local */
    return -1;
  }
}

t8_ctree_t
t8_forest_get_coarse_tree_ext (t8_forest_t forest,
                               t8_locidx_t ltreeid,
                               t8_locidx_t ** face_neigh, int8_t ** ttf)
{
  t8_locidx_t         lctreeid;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* Compute the coarse tree's local id */
  lctreeid = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);

  return t8_cmesh_trees_get_tree_ext (forest->cmesh->trees, lctreeid,
                                      face_neigh, ttf);
}

t8_ctree_t
t8_forest_get_coarse_tree (t8_forest_t forest, t8_locidx_t ltreeid)
{
  return t8_forest_get_coarse_tree_ext (forest, ltreeid, NULL, NULL);
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
    sc_stats_set1 (&stats[4], profile->ghosts_shipped,
                   "forest: Number of ghost elements sent.");
    sc_stats_set1 (&stats[5], profile->ghosts_received,
                   "forest: Number of ghost elements received.");
    sc_stats_set1 (&stats[6], profile->ghosts_remotes,
                   "forest: Number of processes we sent ghosts to/received from.");
    sc_stats_set1 (&stats[7], profile->partition_runtime,
                   "forest: Partition runtime.");
    sc_stats_set1 (&stats[8], profile->commit_runtime,
                   "forest: Commit runtime.");
    sc_stats_set1 (&stats[9], profile->ghost_runtime,
                   "forest: Ghost runtime.");
    /* compute stats */
    sc_stats_compute (sc_MPI_COMM_WORLD, T8_PROFILE_NUM_STATS, stats);
    /* print stats */
    t8_logf (SC_LC_GLOBAL, SC_LP_STATISTICS, "Printing stats for forest.\n");
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS,
                    T8_PROFILE_NUM_STATS, stats, 1, 1);
  }
}

void
t8_forest_compute_elements_offset (t8_forest_t forest)
{
  t8_locidx_t         itree, num_trees;
  t8_locidx_t         current_offset;
  t8_tree_t           tree;

  T8_ASSERT (t8_forest_is_initialized (forest));

  /* Get the number of local trees */
  num_trees = t8_forest_get_num_local_trees (forest);
  current_offset = 0;
  /* Iterate through all trees, sum up the element counts and set it as
   * the element_offsets */
  for (itree = 0; itree < num_trees; itree++) {
    tree = t8_forest_get_tree (forest, itree);
    tree->elements_offset = current_offset;
    current_offset += t8_forest_get_tree_element_count (tree);
  }
  /* At the end, we counted all elements */
  T8_ASSERT (current_offset == forest->local_num_elements);
}

void
t8_forest_write_vtk (t8_forest_t forest, const char *filename)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);

  t8_forest_vtk_write_file (forest, filename, 1, 1, 1, 1);
}

t8_forest_t
t8_forest_new_uniform (t8_cmesh_t cmesh, t8_scheme_cxx_t * scheme,
                       int level, int do_face_ghost, sc_MPI_Comm comm)
{
  t8_forest_t         forest;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (scheme != NULL);
  T8_ASSERT (0 <= level);

  /* Initialize the forest */
  t8_forest_init (&forest);
  /* Set the cmesh, scheme and level */
  t8_forest_set_cmesh (forest, cmesh, comm);
  t8_forest_set_scheme (forest, scheme);
  t8_forest_set_level (forest, level);
  if (do_face_ghost) {
    t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  }
  /* commit the forest */
  t8_forest_commit (forest);

  return forest;
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
  if (forest->scheme_cxx != NULL) {
    t8_scheme_cxx_unref (&forest->scheme_cxx);
  }
  if (forest->cmesh != NULL) {
    t8_cmesh_unref (&forest->cmesh);
  }

  /* free the memory of the offset array */
  if (forest->element_offsets != NULL) {
    t8_shmem_array_destroy (&forest->element_offsets);
  }
  /* free the memory of the global_first_desc array */
  if (forest->global_first_desc != NULL) {
    t8_shmem_array_destroy (&forest->global_first_desc);
  }
  /* free the memory of the tree_offsets array */
  if (forest->tree_offsets != NULL) {
    t8_shmem_array_destroy (&forest->tree_offsets);
  }
  if (forest->profile != NULL) {
    T8_FREE (forest->profile);
  }
  /* Dereference the ghost layer if it exists */
  if (forest->ghosts != NULL) {
    t8_forest_ghost_unref (&forest->ghosts);
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
