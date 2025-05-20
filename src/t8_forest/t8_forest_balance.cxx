/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_c_interface.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_schemes/t8_scheme.hxx>

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
t8_forest_balance_adapt (t8_forest_t forest, t8_forest_t forest_from, const t8_locidx_t ltree_id,
                         const t8_eclass_t tree_class, [[maybe_unused]] const t8_locidx_t lelement_id,
                         const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                         [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int *pdone, iface, num_faces, num_half_neighbors, ineigh;
  t8_gloidx_t neighbor_tree;
  t8_eclass_t neigh_class;
  const t8_element_t *element = elements[0];
  t8_element_t **half_neighbors;

  /* We only need to check an element, if its level is smaller then the maximum
   * level in the forest minus 2.
   * Otherwise there cannot exist neighbors of level greater than the element's level plus one.
   * The variable maxlevel_existing is only set if we enter this function from t8_forest_balance.
   * If we enter from the check function is_balanced, then it may not be set.
   */

  if (forest_from->maxlevel_existing <= 0
      || scheme->element_get_level (tree_class, element) <= forest_from->maxlevel_existing - 2) {

    pdone = (int *) forest->t8code_data;

    num_faces = scheme->element_get_num_faces (tree_class, element);
    for (iface = 0; iface < num_faces; iface++) {
      /* Get the element class and scheme of the face neighbor */
      neigh_class = t8_forest_element_neighbor_eclass (forest_from, ltree_id, element, iface);
      /* Allocate memory for the number of half face neighbors */
      num_half_neighbors = scheme->element_get_num_face_children (tree_class, element, iface);
      half_neighbors = T8_ALLOC (t8_element_t *, num_half_neighbors);
      scheme->element_new (neigh_class, num_half_neighbors, half_neighbors);
      /* Compute the half face neighbors of element at this face */
      neighbor_tree = t8_forest_element_half_face_neighbors (forest_from, ltree_id, element, half_neighbors,
                                                             neigh_class, iface, num_half_neighbors, NULL);
      if (neighbor_tree >= 0) {
        /* The face neighbors do exist, check for each one, whether it has
         * local or ghost leaf descendants in the forest.
         * If so, the element will be refined. */
        for (ineigh = 0; ineigh < num_half_neighbors; ineigh++) {
          if (t8_forest_element_has_leaf_desc (forest_from, neighbor_tree, half_neighbors[ineigh], neigh_class)) {
            /* This element should be refined */
            *pdone = 0;
            /* clean-up */
            scheme->element_destroy (neigh_class, num_half_neighbors, half_neighbors);
            T8_FREE (half_neighbors);
            return 1;
          }
        }
      }
      /* clean-up */
      scheme->element_destroy (neigh_class, num_half_neighbors, half_neighbors);
      T8_FREE (half_neighbors);
    }
  }

  return 0;
}

/* Collective function to compute the maximum occurring refinement level in a forest */
static void
t8_forest_compute_max_element_level (t8_forest_t forest)
{
  t8_locidx_t ielement, elem_in_tree;
  t8_locidx_t itree, num_trees;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  int local_max_level = 0;

  /* Iterate over all local trees and all local elements and comupte the maximum occurring level */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_trees; itree++) {
    elem_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (ielement = 0; ielement < elem_in_tree; ielement++) {
      /* Get the element and compute its level */
      const t8_element_t *elem = t8_forest_get_element_in_tree (forest, itree, ielement);
      const int elem_level = scheme->element_get_level (tree_class, elem);
      local_max_level = SC_MAX (local_max_level, elem_level);
    }
  }
  /* Communicate the local maximum levels */
  sc_MPI_Allreduce (&local_max_level, &forest->maxlevel_existing, 1, sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
}

void
t8_forest_balance (t8_forest_t forest, int repartition)
{
  t8_forest_t forest_temp, forest_from, forest_partition;
  int done = 0, done_global = 0;
  int count_rounds = 0;
  /* The following variables are only required if profiling is
   * enabled. */
  int num_stats_allocated, istats;
  const int stat_alloc_chunk_size = 10; /* How many stats we allocate (and add if we need more) */
  int count_adapt_stats = 0, count_ghost_stats = 0;
  int count_partition_stats = 0;
  double ada_time, ghost_time, part_time;
  sc_statinfo_t *adap_stats, *ghost_stats, *partition_stats;
  int create_ghost_definition = 0; /* flag if create ghost_definition */

  t8_global_productionf ("Into t8_forest_balance with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_elements (forest->set_from));
  t8_log_indent_push ();

  /* Set default value to prevent compiler warning */
  adap_stats = ghost_stats = partition_stats = NULL;

  if (forest->profile != NULL) {
    /* Profiling is enable, so we measure the runtime of balance */
    forest->profile->balance_runtime = -sc_MPI_Wtime ();
    /* We store the individual adapt, ghost, and partition runtimes */
    /* We reserve memory for stat_alloc_chunk_size - 1 many balance rounds
     * (the extra entry is required for the total sum).
     * We will grow the statistics arrays dynamically if more rounds are required. */
    T8_ASSERT (stat_alloc_chunk_size > 0);
    num_stats_allocated = stat_alloc_chunk_size;
    adap_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    ghost_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    if (repartition) {
      partition_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    }
  }

  /* Compute the maximum occurring refinement level in the forest */
  t8_forest_compute_max_element_level (forest->set_from);
  t8_global_productionf ("Computed maximum occurring level:\t%i\n", forest->set_from->maxlevel_existing);
  /* Use set_from as the first forest to adapt */
  forest_from = forest->set_from;
  /* This function is reference neutral regarding forest_from */
  t8_forest_ref (forest_from);

  /* if the set_from forest of the current forest has no ghost layer computed,
   * compute a ghost layer for the set_from forest */
  if (forest->set_from->ghosts == NULL) {
    /* If the forest does not yet have a ghost_definition */
    if (forest->set_from->ghost_definition == NULL) {
      /* create a ghost_definition of type face with top-down-search */
      forest->set_from->ghost_definition = t8_forest_ghost_definition_face_new (3);
      create_ghost_definition = 1;
    }
    /* compute ghost layer for set_from forest */
    t8_forest_ghost_create_topdown (forest->set_from);
    if (create_ghost_definition) { /* if a ghost_definition has been created, it will be deleted here */
      t8_forest_ghost_definition_unref (&(forest->set_from->ghost_definition));
      forest->set_from->ghost_definition = NULL;
    }
  }

  while (!done_global) {
    done = 1;

    T8_ASSERT (forest_from->maxlevel_existing >= 0);
    /* Initialize the temp forest to be adapted from forest_from */
    t8_forest_init (&forest_temp);
    /* Update the maximum occurring level */
    forest_temp->maxlevel_existing = forest_from->maxlevel_existing;
    /* Adapt the forest */
    t8_forest_set_adapt (forest_temp, forest_from, t8_forest_balance_adapt, 0);
    if (!repartition) {
      t8_forest_set_ghost (forest_temp, 1, T8_GHOST_FACES);
    }
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
      if (count_rounds > num_stats_allocated - 2) {
        T8_ASSERT (count_adapt_stats <= count_rounds);
        T8_ASSERT (count_ghost_stats <= count_rounds);
        T8_ASSERT (count_partition_stats <= count_rounds);
        /* re-allocate memory for stats */
        num_stats_allocated += stat_alloc_chunk_size;
        adap_stats = T8_REALLOC (adap_stats, sc_statinfo_t, num_stats_allocated);
        ghost_stats = T8_REALLOC (ghost_stats, sc_statinfo_t, num_stats_allocated);
        if (repartition) {
          partition_stats = T8_REALLOC (partition_stats, sc_statinfo_t, num_stats_allocated);
        }
      }
      sc_stats_set1 (&adap_stats[count_adapt_stats], forest_temp->profile->adapt_runtime, "forest balance: Adapt time");
      count_adapt_stats++;
      if (!repartition) {
        sc_stats_set1 (&ghost_stats[count_ghost_stats], forest_temp->profile->ghost_runtime,
                       "forest balance: Ghost time");
        count_ghost_stats++;
      }
    }

    /* Compute the logical and of all process local done values, if this results
     * in 1 then all processes are finished */
    sc_MPI_Allreduce (&done, &done_global, 1, sc_MPI_INT, sc_MPI_LAND, forest->mpicomm);

    if (repartition && !done_global) {
      /* If repartitioning is used, we partition the forest */
      t8_forest_init (&forest_partition);
      /* Update the maximum occurring level */
      forest_partition->maxlevel_existing = forest_temp->maxlevel_existing;
      t8_forest_set_partition (forest_partition, forest_temp, 0);
      t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
      /* If profiling is enabled, measure partition rumtimes */
      if (forest->profile != NULL) {
        t8_forest_set_profiling (forest_partition, 1);
      }
      t8_forest_commit (forest_partition);

      /* Store the runtimes of partition */
      if (forest->profile != NULL) {
        sc_stats_set1 (&partition_stats[count_partition_stats], forest_partition->profile->partition_runtime,
                       "forest balance: Partition time");
        count_partition_stats++;
        sc_stats_set1 (&ghost_stats[count_ghost_stats], forest_partition->profile->ghost_runtime,
                       "forest balance: Ghost time");
        count_ghost_stats++;
      }

      forest_temp = forest_partition;
      forest_partition = NULL;
    }
    /* Adapt forest_temp in the next round */
    forest_from = forest_temp;
    count_rounds++;
  }

  T8_ASSERT (t8_forest_is_balanced (forest_temp));
  /* Forest_temp is now balanced, we copy its trees and elements to forest */
  t8_forest_copy_trees (forest, forest_temp, 1);
  /* TODO: Also copy ghost elements if ghost creation is set */

  t8_log_indent_pop ();
  t8_global_productionf ("Done t8_forest_balance with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_elements (forest_temp));
  t8_debugf ("t8_forest_balance needed %i rounds.\n", count_rounds);
  /* clean-up */
  t8_forest_unref (&forest_temp);

  if (forest->profile != NULL) {
    /* Profiling is enabled, so we measure the runtime of balance. */
    forest->profile->balance_runtime += sc_MPI_Wtime ();
    forest->profile->balance_rounds = count_rounds;
    /* Print the runtime of adapt/ghost/partition */
    /* Compute the overall runtime and store in last entry */
    ada_time = ghost_time = part_time = 0;
    for (istats = 0; istats < count_adapt_stats; istats++) {
      ada_time += adap_stats[istats].sum_values;
    }
    for (istats = 0; istats < count_ghost_stats; istats++) {
      ghost_time += ghost_stats[istats].sum_values;
    }
    if (repartition) {
      for (istats = 0; istats < count_partition_stats; istats++) {
        part_time += partition_stats[istats].sum_values;
      }
    }
    sc_stats_set1 (&adap_stats[count_adapt_stats], ada_time, "forest balance: Total adapt time");
    sc_stats_set1 (&ghost_stats[count_ghost_stats], ghost_time, "forest balance: Total ghost time");
    if (repartition) {
      sc_stats_set1 (&partition_stats[count_partition_stats], part_time, "forest balance: Total partition time");
    }

    /* Compute and print the intermediate stats */
    T8_ASSERT (count_rounds + 1 <= num_stats_allocated);
    sc_stats_compute (forest->mpicomm, count_adapt_stats + 1, adap_stats);
    sc_stats_compute (forest->mpicomm, count_ghost_stats + 1, ghost_stats);
    if (repartition) {
      sc_stats_compute (forest->mpicomm, count_partition_stats + 1, partition_stats);
    }
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_adapt_stats + 1, adap_stats, 1, 1);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_ghost_stats + 1, ghost_stats, 1, 1);
    if (repartition) {
      sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_partition_stats + 1, partition_stats, 1, 1);
    }
    T8_FREE (adap_stats);
    T8_FREE (ghost_stats);
    if (repartition) {
      T8_FREE (partition_stats);
    }
  }
}

/* Check whether the local elements of a forest are balanced. */
int
t8_forest_is_balanced (t8_forest_t forest)
{
  t8_forest_t forest_from;
  t8_locidx_t num_trees, num_elements;
  t8_locidx_t itree, ielem;
  void *data_temp;
  int dummy_int;

  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  /* temporarily save forest_from */
  forest_from = forest->set_from;

  forest->set_from = forest;

  /* temporarily save forest t8code_data */
  data_temp = forest->t8code_data;
  forest->t8code_data = &dummy_int;

  num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over all trees */
  for (itree = 0; itree < num_trees; itree++) {
    num_elements = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Iterate over all elements of this tree */
    for (ielem = 0; ielem < num_elements; ielem++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielem);
      /* Test if this element would need to be refined in the balance step.
       * If so, the forest is not balanced locally. */
      if (t8_forest_balance_adapt (forest, forest, itree, tree_class, ielem, scheme, 0, 1,
                                   (t8_element_t **) (&element))) {
        forest->set_from = forest_from;
        forest->t8code_data = data_temp;
        return 0;
      }
    }
  }
  forest->set_from = forest_from;
  forest->t8code_data = data_temp;
  return 1;
}

T8_EXTERN_C_END ();
