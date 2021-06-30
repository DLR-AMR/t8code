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
#include <t8_forest/t8_forest_subelements.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This is the adapt function called during one round of balance.
 * We refine an element if it has any face neighbor with a level larger
 * than the element's level + 1.
 */
/* TODO: We currently do not adapt recursively since some functions such
 * as half neighbor computation require the forest to be committed. Thus,
 * we pass forest_from as a parameter. But doing so is not valid anymore
 * if we refine recursively. */

static int
t8_forest_subelements_adapt (t8_forest_t forest, t8_forest_t forest_from,
                         t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         int num_elements, t8_element_t * elements[])
{
  /* NOTE do something */
  return 2;
}

/* Collective function to compute the maximum occurring refinement level in a forest */
static void
t8_forest_compute_max_element_level (t8_forest_t forest)
{
  t8_locidx_t         ielement, elem_in_tree;
  t8_locidx_t         itree, num_trees;
  t8_element_t       *elem;
  t8_eclass_scheme_c *scheme;
  int                 local_max_level = 0, elem_level;

  /* Iterate over all local trees and all local elements and comupte the maximum occurring level */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_trees; itree++) {
    elem_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    scheme =
      t8_forest_get_eclass_scheme (forest,
                                   t8_forest_get_tree_class (forest, itree));
    for (ielement = 0; ielement < elem_in_tree; ielement++) {
      /* Get the element and compute its level */
      elem = t8_forest_get_element_in_tree (forest, itree, ielement);
      elem_level = scheme->t8_element_level (elem);
      local_max_level = SC_MAX (local_max_level, elem_level);
    }
  }
  /* Communicate the local maximum levels */
  sc_MPI_Allreduce (&local_max_level, &forest->maxlevel_existing, 1,
                    sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
}

void
t8_forest_subelements (t8_forest_t forest)
{
  
  t8_forest_t         forest_temp, forest_from, forest_partition;
  int                 done = 0, done_global = 0;
  int                 count = 0, num_stats, i;
  double              ada_time, ghost_time, part_time;
  sc_statinfo_t      *adap_stats, *ghost_stats, *partition_stats;

  t8_global_productionf
    ("Into t8_forest_balance with %lli global elements.\n",
     (long long) t8_forest_get_global_num_elements (forest->set_from));
  t8_log_indent_push ();

  /* Set default value to prevent compiler warning */
  adap_stats = ghost_stats = partition_stats = NULL;

  if (forest->profile != NULL) {
    /* Profiling is enable, so we measure the runtime of balance */
    forest->profile->balance_runtime = -sc_MPI_Wtime ();
    /* We store the individual adapt, ghost, and partition runtimes */
    num_stats = 5;
    adap_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats);
    ghost_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats);
  }

  /* Compute the maximum occurring refinement level in the forest */
  t8_forest_compute_max_element_level (forest->set_from);
  t8_global_productionf ("Computed maximum occurring level:\t%i\n",
                         forest->set_from->maxlevel_existing);
  /* Use set_from as the first forest to adapt */
  forest_from = forest->set_from;
  /* This function is reference neutral regarding forest_from */
  t8_forest_ref (forest_from);

  if (forest->set_from->ghosts == NULL) {
    forest->set_from->ghost_type = T8_GHOST_FACES;
    t8_forest_ghost_create_topdown (forest->set_from);
  }
  
  done = 1;

  T8_ASSERT (forest_from->maxlevel_existing >= 0);
  /* Initialize the temp forest to be adapted from forest_from */
  t8_forest_init (&forest_temp);
  /* Update the maximum occurring level */
  forest_temp->maxlevel_existing = forest_from->maxlevel_existing;
  /* Adapt the forest */
  t8_forest_set_adapt (forest_temp, forest_from, t8_forest_subelements_adapt, 0);
  forest_temp->t8code_data = &done;
  /* If profiling is enabled, measure ghost/adapt rumtimes */
  if (forest->profile != NULL) {
    t8_forest_set_profiling (forest_temp, 1);
  }
  t8_global_productionf ("Profiling: %i\n", forest->profile != NULL);
  /* Adapt the forest */
  t8_forest_commit (forest_temp);
  /* Store the runtimes of adapt and ghost */
  if (forest->profile != NULL) {
    while (count >= num_stats - 2) {
      /* re-allocate memory for stats */
      num_stats += 2;
      adap_stats = T8_REALLOC (adap_stats, sc_statinfo_t, num_stats);
      ghost_stats = T8_REALLOC (ghost_stats, sc_statinfo_t, num_stats);
    }
    sc_stats_set1 (&adap_stats[count], forest_temp->profile->adapt_runtime,
                   "forest balance: Adapt time");
  }

  /* Compute the logical and of all process local done values, if this results
   * in 1 then all processes are finished */
  sc_MPI_Allreduce (&done, &done_global, 1, sc_MPI_INT, sc_MPI_LAND,
                    forest->mpicomm);

  /* Adapt forest_temp in the next round */
  forest_from = forest_temp;
  count++;

  /* Forest_temp is now balanced, we copy its trees and elements to forest */
  t8_forest_copy_trees (forest, forest_temp, 1);
  /* TODO: Also copy ghost elements if ghost creation is set */

  t8_log_indent_pop ();
  t8_global_productionf
    ("Done t8_forest_balance with %lli global elements.\n",
     (long long) t8_forest_get_global_num_elements (forest_temp));
  t8_debugf ("t8_forest_balance needed %i rounds.\n", count);
  /* clean-up */
  t8_forest_unref (&forest_temp);

  if (forest->profile != NULL) {
    /* Profiling is enabled, so we measure the runtime of balance. */
    forest->profile->balance_runtime += sc_MPI_Wtime ();
    forest->profile->balance_rounds = count;
    /* Print the runtime of adapt/ghost/partition */
    /* Compute the overall runtime and store in last entry */
    ada_time = ghost_time = part_time = 0;
    t8_debugf ("ada stats %f\n", adap_stats[count].sum_values);
    for (i = 0; i < count; i++) {
      ada_time += adap_stats[i].sum_values;
      ghost_time += ghost_stats[i].sum_values;
    }
    sc_stats_set1 (&adap_stats[count], ada_time,
                   "forest balance: Total adapt time");
    sc_stats_set1 (&ghost_stats[count], ghost_time,
                   "forest balance: Total ghost time");

    /* Compute and print the intermediate stats */
    sc_stats_compute (forest->mpicomm, count + 1, adap_stats);
    sc_stats_compute (forest->mpicomm, count + 1, ghost_stats);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count + 1,
                    adap_stats, 1, 1);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count + 1,
                    ghost_stats, 1, 1);
    T8_FREE (adap_stats);
    T8_FREE (ghost_stats);
  }
}

/* Check whether all hanging nodes are eliminated. */
int
t8_forest_subelements_used (t8_forest_t forest)
{
  /* NOTE do something */
  return 1;
}

T8_EXTERN_C_END ();
