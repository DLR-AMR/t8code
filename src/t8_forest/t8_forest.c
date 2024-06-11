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
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_transition.h>
#include <t8_forest/t8_forest_vtk.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_element_c_interface.h>

void
t8_forest_init (t8_forest_t *pforest)
{
  t8_forest_t forest;

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
  forest->set_balance = -1;
  forest->maxlevel_existing = -1;
  forest->stats_computed = 0;
  forest->incomplete_trees = -1;
  forest->set_subelements = 0;
  forest->is_transitioned = 0;
}

int
t8_forest_is_initialized (t8_forest_t forest)
{
  return forest != NULL && t8_refcount_is_active (&forest->rc) && !forest->committed;
}

/** Check whether at least one eclass scheme of forest supports transitioning
 * \param [in] forest           A forest
 * \return                      True if at least one eclass scheme in forest has an implementation for subelements
 */
int
t8_forest_supports_transitioning (t8_forest_t forest)
{
  int                 supports_transition = 1;
  int                 supports_transition_all_procs = 0;        /* Result over all procs */
  int                 int_eclass;
  int                 mpiret;
  t8_eclass_scheme_c *tscheme;

  /* Iterate over all eclasses */
  for (int_eclass = (int) T8_ECLASS_ZERO; int_eclass < (int) T8_ECLASS_COUNT;
       int_eclass++) {
    /* If the forest has trees of the current eclass, check if elements of this
     * eclass supports transitioning. */
    if (forest->cmesh->num_local_trees_per_eclass[int_eclass] > 0) {
      tscheme = forest->scheme_cxx->eclass_schemes[int_eclass];
      supports_transition = supports_transition
        && t8_element_scheme_supports_transitioning (tscheme);
    }
  }
  /* Combine the process-local results via a logic or and distribute the
   * result over all procs (in the communicator).*/
  mpiret =
    sc_MPI_Allreduce (&supports_transition, &supports_transition_all_procs, 1,
                      sc_MPI_INT, sc_MPI_LAND, forest->mpicomm);
  SC_CHECK_MPI (mpiret);

  return supports_transition_all_procs;
}

int
t8_forest_is_committed (const t8_forest_t forest)
{
  return forest != NULL && t8_refcount_is_active (&forest->rc) && forest->committed;
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
  int do_dup;
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
t8_forest_set_scheme (t8_forest_t forest, t8_scheme_cxx_t *scheme)
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
  /* Set the from_method to COPY. This overwrites any previous setting of ADAPT, PARTITION, or BALANCE */
  forest->from_method = T8_FOREST_FROM_COPY;

  /* Overwrite any previous setting */
  forest->set_adapt_fn = NULL;
  forest->set_adapt_recursive = -1;
  forest->set_balance = -1;
  forest->set_for_coarsening = -1;
}

void
t8_forest_set_partition (t8_forest_t forest, const t8_forest_t set_from, int set_for_coarsening)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme_cxx == NULL);

  forest->set_for_coarsening = set_for_coarsening;

  if (set_from != NULL) {
    /* If set_from = NULL, we assume a previous forest_from was set */
    forest->set_from = set_from;
  }
  /* Add PARTITION to the from_method.
   * This overwrites T8_FOREST_FROM_COPY */
  if (forest->from_method == T8_FOREST_FROM_LAST) {
    forest->from_method = T8_FOREST_FROM_PARTITION;
  }
  else {
    forest->from_method |= T8_FOREST_FROM_PARTITION;
  }
}

void
t8_forest_set_balance (t8_forest_t forest, const t8_forest_t set_from, int no_repartition)
{
  T8_ASSERT (t8_forest_is_initialized (forest));

  if (no_repartition) {
    /* We do not repartition the forest during balance */
    forest->set_balance = T8_FOREST_BALANCE_NO_REPART;
  }
  else {
    /* Repartition during balance */
    forest->set_balance = T8_FOREST_BALANCE_REPART;
  }

  if (set_from != NULL) {
    /* If set_from = NULL, we assume a previous forest_from was set */
    forest->set_from = set_from;
  }

  /* Add BALANCE to the from_method.
   * This overwrites T8_FOREST_FROM_COPY */
  if (forest->from_method == T8_FOREST_FROM_LAST) {
    forest->from_method = T8_FOREST_FROM_BALANCE;
  }
  else {
    forest->from_method |= T8_FOREST_FROM_BALANCE;
  }
}

void
t8_forest_set_transition (t8_forest_t forest, const t8_forest_t set_from,
                          int set_transition_with_balance)
{
  T8_ASSERT (t8_forest_is_initialized (forest));

  if (set_transition_with_balance) {
    t8_debugf("-------------set transition with balance -------------\n");
    /* balance with repartition */
    t8_forest_set_balance (forest, set_from, 0);
  }


  if (set_from != NULL) {
    /* Note that it is possible to apply transitioning to a forest without transition implementation.
     * In this case, the transition refine routine will return 0, keeping the forest unchanged. 
     * Nevertheless, we assert here in this case. */
    T8_ASSERT (t8_forest_supports_transitioning (set_from));
    /* If set_from = NULL, we assume a previous forest_from was set */
    forest->set_from = set_from;
  }
  else {
    T8_ASSERT (forest->set_from != NULL);
    T8_ASSERT (t8_forest_supports_transitioning (forest->set_from));
  }

  /* Add SUBELEMENTS to the from_method.
   * This overwrites T8_FOREST_FROM_COPY */
  if (forest->from_method == T8_FOREST_FROM_LAST) {
    forest->from_method = T8_FOREST_FROM_TRANSITION;
  }
  else {
    forest->from_method |= T8_FOREST_FROM_TRANSITION;
  }

  /* set the forests subelement flag, which is for example used by the LFN routine */
  forest->set_subelements = 1;
  // forest->is_transitioned = 1;
}

void
t8_forest_set_ghost_ext (t8_forest_t forest, int do_ghost, t8_ghost_type_t ghost_type, int ghost_version)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  /* We currently only support face ghosts */
  SC_CHECK_ABORT (do_ghost == 0 || ghost_type == T8_GHOST_FACES,
                  "Ghost neighbors other than face-neighbors are not supported.\n");
  SC_CHECK_ABORT (1 <= ghost_version && ghost_version <= 3, "Invalid choice for ghost version. Choose 1, 2, or 3.\n");

  if (ghost_type == T8_GHOST_NONE) {
    /* none type disables ghost */
    forest->do_ghost = 0;
  }
  else {
    forest->do_ghost = (do_ghost != 0); /* True if and only if do_ghost != 0 */
  }
  if (forest->do_ghost) {
    forest->ghost_type = ghost_type;
    forest->ghost_algorithm = ghost_version;
  }
}

void
t8_forest_set_ghost (t8_forest_t forest, int do_ghost, t8_ghost_type_t ghost_type)
{
  /* Use ghost version 3, top-down search and for unbalanced forests. */
  t8_forest_set_ghost_ext (forest, do_ghost, ghost_type, 3);
}

void
t8_forest_set_adapt (t8_forest_t forest, const t8_forest_t set_from, t8_forest_adapt_t adapt_fn, int recursive)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
  T8_ASSERT (forest->cmesh == NULL);
  T8_ASSERT (forest->scheme_cxx == NULL);
  T8_ASSERT (forest->set_adapt_fn == NULL);
  T8_ASSERT (forest->set_adapt_recursive == -1);

  forest->set_adapt_fn = adapt_fn;
  forest->set_adapt_recursive = recursive != 0;

  if (set_from != NULL) {
    /* If set_from = NULL, we assume a previous forest_from was set */
    forest->set_from = set_from;
  }

  /* Add ADAPT to the from_method.
   * This overwrites T8_FOREST_FROM_COPY */

  if (forest->from_method == T8_FOREST_FROM_LAST) {
    forest->from_method = T8_FOREST_FROM_ADAPT;
  }
  else {
    forest->from_method |= T8_FOREST_FROM_ADAPT;
  }
}

void
t8_forest_set_user_data (t8_forest_t forest, void *data)
{
  /* TODO: we need a is_initialized and may be committed check! */
  forest->user_data = data;
}

void *
t8_forest_get_user_data (const t8_forest_t forest)
{
  return forest->user_data;
}

void
t8_forest_set_user_function (t8_forest_t forest, t8_generic_function_pointer function)
{
  T8_ASSERT (t8_forest_is_initialized (forest) || t8_forest_is_committed (forest));
  forest->user_function = function;
}

t8_generic_function_pointer
t8_forest_get_user_function (t8_forest_t forest)
{
  //T8_ASSERT (t8_forest_is_initialized (forest) || t8_forest_is_committed (forest));
  return forest->user_function;
}

void
t8_forest_comm_global_num_elements (t8_forest_t forest)
{
  int mpiret;
  t8_gloidx_t local_num_el;
  t8_gloidx_t global_num_el;

  local_num_el = (t8_gloidx_t) forest->local_num_elements;
  mpiret = sc_MPI_Allreduce (&local_num_el, &global_num_el, 1, T8_MPI_GLOIDX, sc_MPI_SUM, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  forest->global_num_elements = global_num_el;
}

/** Adapt callback function to refine every element in the forest.
 * It is merely used to build a new forest with pyramids. 
 * 
 * \param [in] forest       the forest to which the new elements belong
 * \param [in] forest_from  the forest that is adapted.
 * \param [in] which_tree   the local tree containing \a elements
 * \param [in] lelement_id  the local element id in \a forest_old in the tree of the current element
 * \param [in] ts           the eclass scheme of the tree
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements the number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 * \return                  Always return 1, to refine every element
 */
static int
t8_forest_refine_everything (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                             const int num_elements, t8_element_t *elements[])
{

  return 1;
}

/**
 * Check if any tree in a forest refines irregularly.
 * An irregular refining tree is a tree with an element that does not
 * refine into 2^dim children. For example the default implementation
 * of pyramids. 
 * \note This function is MPI collective
 * 
 * \param[in] forest    The forest to check
 * \return          non-zero if any tree refines irregular
 */
static int
t8_forest_refines_irregular (t8_forest_t forest)
{
  int irregular = 0;
  int irregular_all_procs = 0; /* Result over all procs */
  int int_eclass;
  int mpiret;
  t8_eclass_scheme_c *tscheme;
  /* Iterate over all eclasses */
  for (int_eclass = (int) T8_ECLASS_ZERO; int_eclass < (int) T8_ECLASS_COUNT; int_eclass++) {
    /* If the forest has trees of the current eclass, check if elements of this eclass refine irregular. */
    if (forest->cmesh->num_local_trees_per_eclass[int_eclass] > 0) {
      tscheme = t8_forest_get_eclass_scheme_before_commit (forest, (t8_eclass_t) int_eclass);
      irregular = irregular || t8_element_refines_irregular (tscheme);
    }
  }
  /* Combine the process-local results via a logic or and distribute the result over all procs (in the communicator).*/
  mpiret = sc_MPI_Allreduce (&irregular, &irregular_all_procs, 1, sc_MPI_INT, sc_MPI_LOR, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  return irregular_all_procs;
}

/** Algorithm to populate a forest, if any tree refines irregularly.
 * Create the elements on this process given a uniform partition
 * of the coarse mesh. We can not use the function t8_forest_populate, because
 * it assumes a regular refinement for all trees.
 * \param[in] forest  The forest to populate
*/
static void
t8_forest_populate_irregular (t8_forest_t forest)
{
  t8_forest_t forest_zero;
  t8_forest_t forest_tmp;
  t8_forest_t forest_tmp_partition;
  t8_cmesh_ref (forest->cmesh);
  t8_scheme_cxx_ref (forest->scheme_cxx);
  /* We start with a level 0 uniform refinement */
  t8_forest_init (&forest_zero);
  t8_forest_set_level (forest_zero, 0);
  t8_forest_set_cmesh (forest_zero, forest->cmesh, forest->mpicomm);
  t8_forest_set_scheme (forest_zero, forest->scheme_cxx);

  t8_forest_commit (forest_zero);
  /* Up to the specified level we refine every element. */
  for (int i = 1; i <= forest->set_level; i++) {
    t8_forest_init (&forest_tmp);
    t8_forest_set_level (forest_tmp, i);
    t8_forest_set_adapt (forest_tmp, forest_zero, t8_forest_refine_everything, 0);
    t8_forest_commit (forest_tmp);
    /* Partition the forest to even the load */
    t8_forest_init (&forest_tmp_partition);
    t8_forest_set_partition (forest_tmp_partition, forest_tmp, 0);
    t8_forest_commit (forest_tmp_partition);

    forest_zero = forest_tmp_partition;
  }

  /* Copy all elements over to the original forest. */
  t8_forest_copy_trees (forest, forest_zero, 1);

  t8_forest_unref (&forest_tmp_partition);
}

void
t8_forest_commit (t8_forest_t forest)
{
  int mpiret;
  int partitioned = 0;
  sc_MPI_Comm comm_dup;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (!forest->committed);
  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of commit */
    forest->profile->commit_runtime = sc_MPI_Wtime ();
  }

  if (forest->set_from == NULL) {
    /* This forest is constructed solely from its cmesh as a uniform
     * forest */
    T8_ASSERT (forest->mpicomm != sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh != NULL);
    T8_ASSERT (forest->scheme_cxx != NULL);
    T8_ASSERT (forest->from_method == T8_FOREST_FROM_LAST);
    T8_ASSERT (forest->incomplete_trees == -1);
    T8_ASSERT (forest->set_subelements == 0);
    T8_ASSERT (forest->is_transitioned == 0);
    /* dup communicator if requested */
    if (forest->do_dup) {
      mpiret = sc_MPI_Comm_dup (forest->mpicomm, &comm_dup);
      SC_CHECK_MPI (mpiret);
      forest->mpicomm = comm_dup;
    }
    forest->dimension = forest->cmesh->dimension;
    /* Set mpirank and mpisize */
    mpiret = sc_MPI_Comm_size (forest->mpicomm, &forest->mpisize);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Comm_rank (forest->mpicomm, &forest->mpirank);
    SC_CHECK_MPI (mpiret);
    /* Compute the maximum allowed refinement level */
    t8_forest_compute_maxlevel (forest);
    T8_ASSERT (forest->set_level <= forest->maxlevel);

    /* populate a new forest with tree and quadrant objects */
    if (t8_forest_refines_irregular (forest) && forest->set_level > 0) {
      /* On root level we will also use the normal algorithm */
      t8_forest_populate_irregular (forest);
    }
    else {

      t8_forest_populate (forest);
    }
    forest->global_num_trees = t8_cmesh_get_num_trees (forest->cmesh);
    forest->incomplete_trees = 0;
    forest->is_transitioned = 0;
  }
  else { 
                                   /* set_from != NULL */
    t8_forest_t forest_from = forest->set_from; /* temporarily store set_from, since we may overwrite it */

    T8_ASSERT (forest->mpicomm == sc_MPI_COMM_NULL);
    T8_ASSERT (forest->cmesh == NULL);
    T8_ASSERT (forest->scheme_cxx == NULL);
    T8_ASSERT (!forest->do_dup);
    T8_ASSERT (forest->from_method >= T8_FOREST_FROM_FIRST && forest->from_method < T8_FOREST_FROM_LAST);
    T8_ASSERT (forest->set_from->incomplete_trees > -1);

    /* TODO: optimize all this when forest->set_from has reference count one */
    /* TODO: Get rid of duping the communicator */
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
    t8_cmesh_ref (forest->set_from->cmesh);
    t8_scheme_cxx_ref (forest->set_from->scheme_cxx);
    /* set the dimension, cmesh and scheme from the old forest */
    forest->dimension = forest->set_from->dimension;
    forest->cmesh = forest->set_from->cmesh;
    forest->scheme_cxx = forest->set_from->scheme_cxx;
    forest->global_num_trees = forest->set_from->global_num_trees;

    /* Compute the maximum allowed refinement level */
    t8_forest_compute_maxlevel (forest);
    if (forest->from_method == T8_FOREST_FROM_COPY) {
      SC_CHECK_ABORT (forest->set_from != NULL, "No forest to copy from was specified.");
      t8_forest_copy_trees (forest, forest->set_from, 1);
    }
    /* TODO: currently we can only handle copy, adapt, partition, and balance */

    /* T8_ASSERT (forest->from_method == T8_FOREST_FROM_COPY); */
    if (forest->from_method & T8_FOREST_FROM_ADAPT) {
      SC_CHECK_ABORT (forest->set_adapt_fn != NULL, "No adapt function specified");

      forest->from_method -= T8_FOREST_FROM_ADAPT;
      if (forest->from_method > 0) {
        /* The forest should also be partitioned/balanced/transitioned.
         * We first untransition the forest, then adapt the forest, then balance and then partition and then possibly transition the forest again*/
        if (forest->set_from->is_transitioned){
          T8_ASSERT(forest->is_transitioned == 0);
          t8_forest_untransition(forest);
        }
        
        t8_forest_t forest_adapt;

        t8_forest_init (&forest_adapt);       
        /* forest_adapt should not change ownership of forest->set_from */
        if (forest_from == forest->set_from) {
            t8_forest_ref (forest->set_from);
         }
        /* set user data of forest to forest_adapt */
        t8_forest_set_user_data (forest_adapt, t8_forest_get_user_data (forest));
        /* Construct an intermediate, adapted forest */
        t8_forest_set_adapt (forest_adapt, forest->set_from, forest->set_adapt_fn, forest->set_adapt_recursive);
        /* Set profiling if enabled */
        t8_forest_set_profiling (forest_adapt, forest->profile != NULL);
        t8_forest_commit (forest_adapt);
        /* The new forest will be partitioned/balanced from forest_adapt */
        forest->set_from = forest_adapt;
        /* Set the user data of forest_from to forest_adapt */
        t8_forest_set_user_data (forest_adapt, t8_forest_get_user_data (forest_from));
        /* If profiling is enabled copy the runtime of adapt. */
        if (forest->profile != NULL) {
          forest->profile->adapt_runtime = forest_adapt->profile->adapt_runtime;
        }
      }
      else {
        /* This forest should only be adapted */
        t8_forest_copy_trees (forest, forest->set_from, 0);
        
        t8_forest_adapt (forest);        
      }
    }
    if (forest->from_method & T8_FOREST_FROM_PARTITION) {
      partitioned = 1;
      /* Partition this forest */
      forest->from_method -= T8_FOREST_FROM_PARTITION;

      if (forest->from_method > 0) {
        /* The forest should also be balanced/transitioned after partition */
        t8_forest_t forest_partition;

        t8_forest_init (&forest_partition);
        if (forest_from == forest->set_from) {
          /* forest_partition should not change ownership of forest->set_from */
          t8_forest_ref (forest->set_from);
        }
        t8_forest_set_partition (forest_partition, forest->set_from, forest->set_for_coarsening);
        /* activate profiling, if this forest has profiling */
        t8_forest_set_profiling (forest_partition, forest->profile != NULL);
        /* Commit the partitioned forest */
        t8_forest_commit (forest_partition);
        forest->set_from = forest_partition;
        if (forest->profile != NULL) {
          forest->profile->partition_bytes_sent = forest_partition->profile->partition_bytes_sent;
          forest->profile->partition_elements_recv = forest_partition->profile->partition_elements_recv;
          forest->profile->partition_elements_shipped = forest_partition->profile->partition_elements_shipped;
          forest->profile->partition_procs_sent = forest_partition->profile->partition_procs_sent;
          forest->profile->partition_runtime = forest_partition->profile->partition_runtime;
        }
      }
      else {
        forest->incomplete_trees = forest->set_from->incomplete_trees;
        /* Partitioning is the last routine, no balance was set */
        forest->global_num_elements = forest->set_from->global_num_elements;
        /* Initialize the trees array of the forest */
        forest->trees = sc_array_new (sizeof (t8_tree_struct_t));
        /* partition the forest */
        t8_forest_partition (forest);

      }
    }
    if (forest->from_method & T8_FOREST_FROM_BALANCE) {
      /* balance the forest */
      forest->from_method -= T8_FOREST_FROM_BALANCE;
      if (forest->from_method > 0) {
        int flag_rep;
        if (forest->set_balance == T8_FOREST_BALANCE_NO_REPART) {
          /* balance without repartition */
          flag_rep = 1;
        }
        else {
          /* balance with repartition */
          flag_rep = 0;
        }
        /* in this case, we will transition after balancing */
        t8_forest_t forest_balance;
        t8_forest_init (&forest_balance);
        /* forest_balance should not change ownership of forest->set_from */
         if (forest_from == forest->set_from) {
            t8_forest_ref (forest->set_from);
         }
        t8_forest_set_balance (forest_balance, forest->set_from, flag_rep);
        /* Set profiling if enabled */
        t8_forest_set_profiling (forest_balance, forest->profile != NULL);
        t8_forest_commit (forest_balance);
        /* The new forest will be partitioned/transitioned from forest_balance */
        forest->set_from = forest_balance;
      }
      else{
      /* Only execute t8_balance.*/
      T8_ASSERT (forest->from_method == 0);

      /* This forest should only be balanced */
      if (forest->set_balance == T8_FOREST_BALANCE_NO_REPART) {
        /* balance without repartition */
        t8_forest_balance (forest, 0);
      }
      else {
        /* balance with repartition */
        t8_forest_balance (forest, 1);
      }
      }
    }
    if (forest->from_method & T8_FOREST_FROM_TRANSITION) {
      forest->from_method -= T8_FOREST_FROM_TRANSITION;
      /* this is the last from method that we execute,
       * nothing should be left todo */
      T8_ASSERT (forest->from_method == 0);
      /* use subelements */
      t8_forest_transition (forest);
      forest->is_transitioned = 1;
    }
    if (forest_from != forest->set_from) {
      /* decrease reference count of intermediate input forest, possibly destroying it */
      t8_forest_unref (&forest->set_from);
    }
    /* reset forest->set_from */
    forest->set_from = forest_from;

    /* decrease reference count of input forest, possibly destroying it */
    t8_forest_unref (&forest->set_from);
  } /* end set_from != NULL */
  /* Compute the element offset of the trees */
  t8_forest_compute_elements_offset (forest);

  /* Compute first and last descendant for each tree */
  t8_forest_compute_desc (forest);

  /* we do not need the set parameters anymore */
  forest->set_level = 0;
  forest->set_for_coarsening = 0;
  forest->set_from = NULL;
  forest->committed = 1;
  t8_debugf("Forest is transitioned %i \n", forest->is_transitioned);
  t8_debugf ("Committed forest with %li local elements and %lli "
             "global elements.\n\tTree range is from %lli to %lli.\n",
             (long) forest->local_num_elements, (long long) forest->global_num_elements,
             (long long) forest->first_local_tree, (long long) forest->last_local_tree);

  if (forest->tree_offsets == NULL) {
    /* Compute the tree offset array */
    t8_forest_partition_create_tree_offsets (forest);

  }

  if (forest->element_offsets == NULL) {
    /* Compute element offsets */
    t8_forest_partition_create_offsets (forest);

  }
  if (forest->global_first_desc == NULL) {
    /* Compute global first desc array */
    t8_forest_partition_create_first_desc (forest);
  }

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of commit */
    forest->profile->commit_runtime = sc_MPI_Wtime () - forest->profile->commit_runtime;
  }

  /* From here on, the forest passes the t8_forest_is_committed check */

  /* re-partition the cmesh */
  if (forest->cmesh->set_partition && partitioned) {
    t8_forest_partition_cmesh (forest, forest->mpicomm, forest->profile != NULL);
  }

  if (forest->mpisize > 1) {
    /* Construct a ghost layer, if desired */
    if (forest->do_ghost) {
      /* TODO: ghost type */
      switch (forest->ghost_algorithm) {
      case 1:
        t8_forest_ghost_create_balanced_only (forest);
        break;
      case 2:
        t8_forest_ghost_create (forest);
        break;
      case 3:
        t8_forest_ghost_create_topdown (forest);
        break;
      default:
        SC_ABORT ("Invalid choice of ghost algorithm");
      }
    }
    forest->do_ghost = 0;
  }
#ifdef T8_ENABLE_DEBUG
  t8_forest_partition_test_boundary_element (forest);
#endif
}

t8_locidx_t
t8_forest_get_local_num_elements (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->local_num_elements;
}

t8_gloidx_t
t8_forest_get_global_num_elements (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->global_num_elements;
}

t8_gloidx_t
t8_forest_get_global_num_subelements (t8_forest_t forest)
{
  T8_ASSERT (forest->global_num_subelements <= forest->global_num_elements);
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->global_num_subelements;
}

void
t8_forest_comm_global_num_subelements (t8_forest_t forest)
{
  int                 mpiret;
  t8_gloidx_t         local_num_subel;
  t8_gloidx_t         global_num_subel;

  local_num_subel = (t8_gloidx_t) forest->local_num_subelements;
  mpiret = sc_MPI_Allreduce (&local_num_subel, &global_num_subel, 1,
                             T8_MPI_GLOIDX, sc_MPI_SUM, forest->mpicomm);
  SC_CHECK_MPI (mpiret);
  forest->global_num_subelements = global_num_subel;
}

t8_locidx_t
t8_forest_get_num_ghosts (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Return the number of ghost elements, or 0 if no ghost structure exists. */
  if (forest->ghosts == NULL) {
    return 0;
  }
  return forest->ghosts->num_ghosts_elements;
}

/* Currently this function is not used */

/* Compute the offset array for a partition cmesh that should match the
 * forest's partition.
 */
static t8_shmem_array_t
t8_forest_compute_cmesh_offset (t8_forest_t forest, sc_MPI_Comm comm)
{
  t8_shmem_array_t offset;

  if (forest->tree_offsets == NULL) {
    /* Create the tree offsets if necessary */
    t8_forest_partition_create_tree_offsets (forest);
  }
  /* initialize the shared memory array */
  t8_shmem_array_init (&offset, sizeof (t8_gloidx_t), forest->mpisize + 1, comm);

  /* Copy the contents */
  t8_shmem_array_copy (offset, forest->tree_offsets);

  return offset;
}

void
t8_forest_partition_cmesh (t8_forest_t forest, sc_MPI_Comm comm, int set_profiling)
{
  t8_cmesh_t cmesh_partition;
  t8_shmem_array_t offsets;

  t8_debugf ("Partitioning cmesh according to forest\n");

  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, forest->cmesh);
  /* set partition range of new cmesh according to forest trees */
  if (forest->tree_offsets == NULL) {
    t8_forest_partition_create_tree_offsets (forest);
  }
  offsets = t8_forest_compute_cmesh_offset (forest, comm);

  t8_cmesh_set_partition_offsets (cmesh_partition, offsets);
  /* Set the profiling of the cmesh */
  t8_cmesh_set_profiling (cmesh_partition, set_profiling);
  /* Commit the new cmesh */
  t8_cmesh_commit (cmesh_partition, comm);
  /* set the new cmesh as the cmesh of the forest */
  forest->cmesh = cmesh_partition;
  t8_debugf ("Done partitioning cmesh\n");
}

sc_MPI_Comm
t8_forest_get_mpicomm (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return forest->mpicomm;
}

t8_gloidx_t
t8_forest_get_first_local_tree_id (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->first_local_tree;
}

t8_locidx_t
t8_forest_get_num_ghost_trees (const t8_forest_t forest)
{
  if (forest->ghosts != NULL) {
    return t8_forest_ghost_num_trees (forest);
  }
  else {
    return 0;
  }
}

t8_locidx_t
t8_forest_get_num_local_trees (const t8_forest_t forest)
{
  t8_locidx_t num_trees;

  num_trees = forest->last_local_tree - forest->first_local_tree + 1;
  /* assert for possible overflow */
  T8_ASSERT ((t8_gloidx_t) num_trees == forest->last_local_tree - forest->first_local_tree + 1);
  if (num_trees < 0) {
    /* Set number of trees to zero if there are none */
    num_trees = 0;
  }
  return num_trees;
}

t8_gloidx_t
t8_forest_get_num_global_trees (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->global_num_trees;
}

t8_gloidx_t
t8_forest_global_tree_id (const t8_forest_t forest, const t8_locidx_t ltreeid)
{
  t8_locidx_t num_local_trees;
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest) + t8_forest_ghost_num_trees (forest));

  num_local_trees = t8_forest_get_num_local_trees (forest);
  if (ltreeid < num_local_trees) {
    /* The tree is a local tree */
    return ltreeid + forest->first_local_tree;
  }
  else {
    T8_ASSERT (forest->ghosts != NULL);
    /* Return the global id of the ghost tree */
    return t8_forest_ghost_get_global_treeid (forest, ltreeid - num_local_trees);
  }
}

/* TODO: We use this function in forest_partition when the
 * forest is only partially committed. Thus, we cannot check whether the
 * forest is committed here. */
t8_tree_t
t8_forest_get_tree (const t8_forest_t forest, const t8_locidx_t ltree_id)
{
  T8_ASSERT (forest->trees != NULL);
  T8_ASSERT (0 <= ltree_id && ltree_id < t8_forest_get_num_local_trees (forest));
  return (t8_tree_t) t8_sc_array_index_locidx (forest->trees, ltree_id);
}

double *
t8_forest_get_tree_vertices (t8_forest_t forest, t8_locidx_t ltreeid)
{
  return t8_cmesh_get_tree_vertices (forest->cmesh, t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid));
}

t8_element_array_t *
t8_forest_tree_get_leaves (const t8_forest_t forest, const t8_locidx_t ltree_id)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltree_id && ltree_id < t8_forest_get_num_local_trees (forest));

  return &t8_forest_get_tree (forest, ltree_id)->elements;
}

t8_cmesh_t
t8_forest_get_cmesh (const t8_forest_t forest)
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
  t8_locidx_t leid = *(const t8_locidx_t *) lelement_id;
  const t8_tree_t tree = (const t8_tree_t) ltree;

  if (tree->elements_offset > leid) {
    /* We have to look further to the left */
    return -1;
  }
  else if (tree->elements_offset + (t8_locidx_t) t8_element_array_get_count (&tree->elements) > leid) {
    /* We have found the tree */
    return 0;
  }
  else {
    /* We have to look further right */
    return 1;
  }
}

t8_element_t *
t8_forest_get_element (t8_forest_t forest, t8_locidx_t lelement_id, t8_locidx_t *ltreeid)
{
  t8_tree_t tree;
  t8_locidx_t ltree;
#ifdef T8_ENABLE_DEBUG
  t8_locidx_t ltreedebug;
#endif

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (lelement_id >= 0);
  if (lelement_id >= t8_forest_get_local_num_elements (forest)) {
    return NULL;
  }
  /* We optimized the binary search out by using sc_bsearch,
   * but keep it in for debugging. We check whether the hand-written
   * binary search matches the sc_array_bsearch. */
#ifdef T8_ENABLE_DEBUG
  {
    t8_locidx_t ltree_a, ltree_b;
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
      else if (tree->elements_offset + (t8_locidx_t) t8_element_array_get_count (&tree->elements) > lelement_id) {
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
  ltree = sc_array_bsearch (forest->trees, &lelement_id, t8_forest_compare_elem_tree);
  T8_ASSERT (ltreedebug == ltree);
  if (ltreeid != NULL) {
    *ltreeid = ltree;
  }

  /* The tree that contains the element is now local tree ltree.
   * Or the element is not a local element. */
  tree = t8_forest_get_tree (forest, ltree);
  if (tree->elements_offset <= lelement_id
      && lelement_id < tree->elements_offset + (t8_locidx_t) t8_element_array_get_count (&tree->elements)) {
    return t8_element_array_index_locidx (&tree->elements, lelement_id - tree->elements_offset);
  }
  /* The element was not found.
   * This case is covered by the first if and should therefore never happen. */
  SC_ABORT_NOT_REACHED ();
  return NULL;
}

const t8_element_t *
t8_forest_get_element_in_tree (t8_forest_t forest, t8_locidx_t ltreeid, t8_locidx_t leid_in_tree)
{
  t8_tree_t tree;
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  tree = t8_forest_get_tree (forest, ltreeid);
  return t8_forest_get_tree_element (tree, leid_in_tree);
}

t8_locidx_t
t8_forest_get_tree_element_offset (const t8_forest_t forest, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return t8_forest_get_tree (forest, ltreeid)->elements_offset;
}

t8_locidx_t
t8_forest_get_tree_element_count (t8_tree_t tree)
{
  t8_locidx_t element_count;

  T8_ASSERT (tree != NULL);
  element_count = t8_element_array_get_count (&tree->elements);
  /* check for type conversion errors */
  T8_ASSERT ((size_t) element_count == t8_element_array_get_count (&tree->elements));
  return element_count;
}

t8_locidx_t
t8_forest_get_tree_num_elements (t8_forest_t forest, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  return t8_forest_get_tree_element_count (t8_forest_get_tree (forest, ltreeid));
}

t8_eclass_t
t8_forest_get_tree_class (const t8_forest_t forest, const t8_locidx_t ltreeid)
{
  t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  T8_ASSERT (0 <= ltreeid && ltreeid < num_local_trees + t8_forest_get_num_ghost_trees (forest));
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
  T8_ASSERT (t8_forest_is_committed (forest));

  if (forest->element_offsets != NULL) {
    return t8_shmem_array_get_gloidx (forest->element_offsets, forest->mpirank);
  }
  return -1;
}

t8_scheme_cxx_t *
t8_forest_get_scheme (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->scheme_cxx != NULL);

  return forest->scheme_cxx;
}

t8_eclass_scheme_c *
t8_forest_get_eclass_scheme (const t8_forest_t forest, const t8_eclass_t eclass)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->scheme_cxx != NULL);
  T8_ASSERT (eclass != T8_ECLASS_COUNT);

  return forest->scheme_cxx->eclass_schemes[eclass];
}

t8_eclass_scheme_c *
t8_forest_get_eclass_scheme_before_commit (t8_forest_t forest, t8_eclass_t eclass)
{
  T8_ASSERT (t8_forest_is_initialized (forest));
  T8_ASSERT (forest->scheme_cxx != NULL);
  T8_ASSERT (eclass != T8_ECLASS_COUNT);

  return forest->scheme_cxx->eclass_schemes[eclass];
}

t8_eclass_t
t8_forest_get_eclass (const t8_forest_t forest, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  return t8_forest_get_tree (forest, ltreeid)->eclass;
}

t8_locidx_t
t8_forest_get_local_id (const t8_forest_t forest, const t8_gloidx_t gtreeid)
{
  t8_gloidx_t ltreeid;
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_forest_get_num_global_trees (forest));

  /* If the tree is local then its local id is the global id minus the
   * first global tree id on this forest. If this number is not in the
   * range of local tree ids then the tree is not local. */
  /* we use a gloidx for ltreeid to prevent overflow and false positives */
  ltreeid = gtreeid - t8_forest_get_first_local_tree_id (forest);
  /* Check if this tree is a local tree and if so return its local id */
  if (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest)) {
    return ltreeid;
  }
  else {
    return -1;
  }
}
t8_locidx_t
t8_forest_get_local_or_ghost_id (const t8_forest_t forest, const t8_gloidx_t gtreeid)
{
  t8_gloidx_t ltreeid;
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_forest_get_num_global_trees (forest));

  /* If the tree is local then its local id is the global id minus the
   * first global tree id on this forest. If this number is not in the
   * range of local tree ids then the tree is not local. */
  /* we use a gloidx for ltreeid to prevent overflow and false positives */
  ltreeid = gtreeid - t8_forest_get_first_local_tree_id (forest);
  /* Check if this tree is a local tree and if so return its local id */
  if (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest)) {
    return ltreeid;
  }
  else {
    t8_locidx_t ghost_id = t8_forest_ghost_get_ghost_treeid (forest, gtreeid);
    if (ghost_id >= 0)
      return t8_forest_get_num_local_trees (forest) + ghost_id;
    return -1;
  }
}

t8_locidx_t
t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest, t8_locidx_t ltreeid)
{
  t8_gloidx_t cmesh_gfirst;
  t8_locidx_t num_local_trees;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->cmesh != NULL);
  num_local_trees = t8_forest_get_num_local_trees (forest);
  T8_ASSERT (0 <= ltreeid && ltreeid < num_local_trees + t8_forest_ghost_num_trees (forest));

  if (ltreeid < num_local_trees) {
    /* This a local tree and not a ghost */

    cmesh_gfirst = t8_cmesh_get_first_treeid (forest->cmesh);
    return forest->first_local_tree - cmesh_gfirst + ltreeid;
  }
  else {
    /* This is a ghost */
    t8_gloidx_t globalid;
    t8_locidx_t cmesh_local_id;
    /* Compute the global id of this ghost tree */
    globalid = t8_forest_ghost_get_global_treeid (forest, ltreeid - num_local_trees);
    /* Compute the cmesh local id of the ghost */
    cmesh_local_id = t8_cmesh_get_local_id (forest->cmesh, globalid);
    /* is < 0 if this ghost does not exist */
    T8_ASSERT (cmesh_local_id >= 0);
    return cmesh_local_id;
  }
}

t8_locidx_t
t8_forest_cmesh_ltreeid_to_ltreeid (t8_forest_t forest, t8_locidx_t lctreeid)
{
  t8_locidx_t ltreeid;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->cmesh != NULL);

  ltreeid = t8_cmesh_get_first_treeid (forest->cmesh) - t8_forest_get_first_local_tree_id (forest) + lctreeid;
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
t8_forest_get_coarse_tree_ext (t8_forest_t forest, t8_locidx_t ltreeid, t8_locidx_t **face_neigh, int8_t **ttf)
{
  t8_locidx_t lctreeid;

  T8_ASSERT (t8_forest_is_committed (forest));

  /* Compute the coarse tree's local id */
  lctreeid = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltreeid);

  return t8_cmesh_trees_get_tree_ext (forest->cmesh->trees, lctreeid, face_neigh, ttf);
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
t8_forest_compute_profile (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    /* Only print something if profiling is enabled */
    t8_profile_t *profile = forest->profile;

    /* Set the stats */
    sc_stats_set1 (&forest->stats[0], profile->partition_elements_shipped, "forest: Number of elements sent.");
    sc_stats_set1 (&forest->stats[1], profile->partition_elements_recv, "forest: Number of elements received.");
    sc_stats_set1 (&forest->stats[2], profile->partition_bytes_sent, "forest: Number of bytes sent.");
    sc_stats_set1 (&forest->stats[3], profile->partition_procs_sent, "forest: Number of processes sent to.");
    sc_stats_set1 (&forest->stats[4], profile->ghosts_shipped, "forest: Number of ghost elements sent.");
    sc_stats_set1 (&forest->stats[5], profile->ghosts_received, "forest: Number of ghost elements received.");
    sc_stats_set1 (&forest->stats[6], profile->ghosts_remotes,
                   "forest: Number of processes we sent ghosts to/received from.");
    sc_stats_set1 (&forest->stats[7], profile->adapt_runtime, "forest: Adapt runtime.");
    sc_stats_set1 (&forest->stats[8], profile->partition_runtime, "forest: Partition runtime.");
    sc_stats_set1 (&forest->stats[9], profile->commit_runtime, "forest: Commit runtime.");
    sc_stats_set1 (&forest->stats[10], profile->ghost_runtime, "forest: Ghost runtime.");
    sc_stats_set1 (&forest->stats[11], profile->ghost_waittime, "forest: Ghost waittime.");
    sc_stats_set1 (&forest->stats[12], profile->balance_runtime, "forest: Balance runtime.");
    sc_stats_set1 (&forest->stats[13], profile->balance_rounds, "forest: Balance rounds.");
    /* compute stats */
    sc_stats_compute (sc_MPI_COMM_WORLD, T8_PROFILE_NUM_STATS, forest->stats);
    forest->stats_computed = 1;
  }
}

void
t8_forest_print_profile (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    /* Compute the stats if not already computed. */
    if (!forest->stats_computed) {
      t8_forest_compute_profile (forest);
    }
    /* print stats */
    t8_logf (SC_LC_GLOBAL, SC_LP_STATISTICS, "Printing stats for forest.\n");
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, T8_PROFILE_NUM_STATS, forest->stats, 1, 1);
  }
}

const sc_statinfo_t *
t8_forest_profile_get_adapt_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[7];
}

const sc_statinfo_t *
t8_forest_profile_get_ghost_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[10];
}

const sc_statinfo_t *
t8_forest_profile_get_partition_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[8];
}

const sc_statinfo_t *
t8_forest_profile_get_commit_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[9];
}

const sc_statinfo_t *
t8_forest_profile_get_balance_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[12];
}

const sc_statinfo_t *
t8_forest_profile_get_balance_rounds_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->profile != NULL);
  T8_ASSERT (forest->stats_computed);
  return &forest->stats[13];
}

double
t8_forest_profile_get_adapt_time (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    return forest->profile->adapt_runtime;
  }
  return 0;
}

double
t8_forest_profile_get_partition_time (t8_forest_t forest, int *procs_sent)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    *procs_sent = forest->profile->partition_procs_sent;
    return forest->profile->partition_runtime;
  }
  return 0;
}

double
t8_forest_profile_get_balance_time (t8_forest_t forest, int *balance_rounds)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    *balance_rounds = forest->profile->balance_rounds;
    return forest->profile->balance_runtime;
  }
  return 0;
}

double
t8_forest_profile_get_ghost_time (t8_forest_t forest, t8_locidx_t *ghosts_sent)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    *ghosts_sent = forest->profile->ghosts_shipped;
    return forest->profile->ghost_runtime;
  }
  *ghosts_sent = 0;
  return 0;
}

double
t8_forest_profile_get_ghostexchange_waittime (t8_forest_t forest)
{

  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    return forest->profile->ghost_waittime;
  }
  return 0;
}

double
t8_forest_profile_get_balance (t8_forest_t forest, int *balance_rounds)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  if (forest->profile != NULL) {
    *balance_rounds = forest->profile->balance_rounds;
    return forest->profile->balance_runtime;
  }
  return 0;
}

void
t8_forest_compute_elements_offset (t8_forest_t forest)
{
 
  t8_locidx_t itree, num_trees;
  t8_locidx_t current_offset;
  t8_tree_t tree;

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

int
t8_forest_write_vtk_ext (t8_forest_t forest, const char *fileprefix, const int write_treeid, const int write_mpirank,
                         const int write_level, const int write_element_id, const int write_ghosts,
                         const int write_curved, int do_not_use_API, const int num_data, t8_vtk_data_field_t *data)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);

#if T8_WITH_VTK
  if (do_not_use_API && write_curved) {
    t8_errorf ("WARNING: Export of curved elements not yet available with the inbuild function. "
               "Using the VTK API instead.\n");
    do_not_use_API = 0;
  }
#else
  /* We are not linked against the VTK library, so
   * we do not use the API by default.
   */
  if (write_curved) {
    t8_errorf ("WARNING: Export of curved elements not yet available with the inbuild function. "
               "Please link to VTK.\n"
               "Using the inbuild function to write out uncurved elements instead.\n");
  }
  do_not_use_API = 1;
#endif
  if (!do_not_use_API) {
    return t8_forest_vtk_write_file_via_API (forest, fileprefix, write_treeid, write_mpirank, write_level,
                                             write_element_id, write_ghosts, write_curved, num_data, data);
  }
  else {
    T8_ASSERT (!write_curved);
    return t8_forest_vtk_write_file (forest, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                     write_ghosts, num_data, data);
  }
}

int
t8_forest_write_vtk (t8_forest_t forest, const char *fileprefix)
{
  return t8_forest_write_vtk_ext (forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, 0, NULL);
}

t8_forest_t
t8_forest_new_uniform (t8_cmesh_t cmesh, t8_scheme_cxx_t *scheme, const int level, const int do_face_ghost,
                       sc_MPI_Comm comm)
{
  t8_forest_t forest;

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (scheme != NULL);
  T8_ASSERT (0 <= level);

  /* Initialize the forest */
  t8_forest_init (&forest);
  T8_ASSERT(t8_forest_is_initialized(forest));
  /* Set the cmesh, scheme and level */
  t8_forest_set_cmesh (forest, cmesh, comm);
  t8_forest_set_scheme (forest, scheme);
  t8_forest_set_level (forest, level);
  if (do_face_ghost) {
    t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  }
  /* commit the forest */
  t8_forest_commit (forest);
  t8_global_productionf ("Constructed uniform forest with %lli global elements.\n",
                         (long long) forest->global_num_elements);

  return forest;
}

t8_forest_t
t8_forest_new_adapt (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int recursive, int do_face_ghost,
                     void *user_data)
{
  t8_forest_t forest;

  t8_forest_init (&forest);
  t8_forest_set_adapt (forest, forest_from, adapt_fn, recursive);
  t8_forest_set_ghost (forest, do_face_ghost, T8_GHOST_FACES);
  if (user_data != NULL) {
    t8_forest_set_user_data (forest, user_data);
  }
  t8_forest_commit (forest);
  return forest;
}

/* Iterate through all the trees and free the element memory as well as
 * the tree memory.
 */
static void
t8_forest_free_trees (t8_forest_t forest)
{
  t8_tree_t tree;
  t8_locidx_t jt, number_of_trees;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->committed);

  number_of_trees = forest->trees->elem_count;
  for (jt = 0; jt < number_of_trees; jt++) {
    tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, jt);

    if (t8_forest_get_tree_element_count (tree) < 1) {
      /* if local tree is empty */
      T8_ASSERT (forest->incomplete_trees);
      continue;
    }
    t8_element_array_reset (&tree->elements);
    /* destroy first and last descendant */
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, jt);
    const t8_eclass_scheme_c *scheme = forest->scheme_cxx->eclass_schemes[eclass];
    t8_element_destroy (scheme, 1, &tree->first_desc);
    t8_element_destroy (scheme, 1, &tree->last_desc);
  }
  sc_array_destroy (forest->trees);
}

/* Completely destroy a forest and unreference all structs that the
 * forest has taken ownership on.
 */
static void
t8_forest_reset (t8_forest_t *pforest)
{
  int mpiret;
  t8_forest_t forest;

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

  /* Destroy the ghost layer if it exists */
  if (forest->ghosts != NULL) {
    t8_forest_ghost_unref (&forest->ghosts);
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
t8_forest_unref (t8_forest_t *pforest)
{
  t8_forest_t forest;

  T8_ASSERT (pforest != NULL);
  forest = *pforest;
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest != NULL);
  if (t8_refcount_unref (&forest->rc)) {
      t8_forest_reset (pforest);
  }
}