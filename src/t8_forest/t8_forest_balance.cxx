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
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest.h>
#include <t8_data/t8_face_neighbor_hash_table.h>
#include <t8_data/t8_locidx_list.h>
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
t8_forest_balance_adapt (t8_forest_t forest, t8_forest_t forest_from,
                         t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         int num_elements, const t8_element_t * elements[])
{
  int                *pdone, iface, num_faces, num_half_neighbors, ineigh;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  const t8_element_t *element = elements[0];
  t8_element_t      **half_neighbors;

  /* We only need to check an element, if its level is smaller then the maximum
   * level in the forest minus 2.
   * Otherwise there cannot exist neighbors of level greater than the element's level plus one.
   * The variable maxlevel_existing is only set if we enter this function from t8_forest_balance.
   * If we enter from the check function is_balanced, then it may not be set.
   */

  if (forest_from->maxlevel_existing <= 0 ||
      ts->t8_element_level (element) <= forest_from->maxlevel_existing - 2) {

    pdone = (int *) forest->t8code_data;

    num_faces = ts->t8_element_num_faces (element);
    for (iface = 0; iface < num_faces; iface++) {
      /* Get the element class and scheme of the face neighbor */
      neigh_class = t8_forest_element_neighbor_eclass (forest_from,
                                                       ltree_id, element,
                                                       iface);
      neigh_scheme = t8_forest_get_eclass_scheme (forest_from, neigh_class);
      /* Allocate memory for the number of half face neighbors */
      num_half_neighbors = ts->t8_element_num_face_children (element, iface);
      half_neighbors = T8_ALLOC (t8_element_t *, num_half_neighbors);
      neigh_scheme->t8_element_new (num_half_neighbors, half_neighbors);
      /* Compute the half face neighbors of element at this face */
      neighbor_tree = t8_forest_element_half_face_neighbors (forest_from,
                                                             ltree_id,
                                                             element,
                                                             half_neighbors,
                                                             neigh_scheme,
                                                             iface,
                                                             num_half_neighbors,
                                                             NULL);
      if (neighbor_tree >= 0) {
        /* The face neighbors do exist, check for each one, whether it has
         * local or ghost leaf descendants in the forest.
         * If so, the element will be refined. */
        for (ineigh = 0; ineigh < num_half_neighbors; ineigh++) {
          if (t8_forest_element_has_leaf_desc (forest_from, neighbor_tree,
                                               half_neighbors[ineigh],
                                               neigh_scheme)) {
            /* This element should be refined */
            *pdone = 0;
            /* clean-up */
            neigh_scheme->t8_element_destroy (num_half_neighbors,
                                              half_neighbors);
            T8_FREE (half_neighbors);
            return 1;
          }
        }
      }
      /* clean-up */
      neigh_scheme->t8_element_destroy (num_half_neighbors, half_neighbors);
      T8_FREE (half_neighbors);
    }
  }

  return 0;
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
t8_forest_balance (t8_forest_t forest, int repartition)
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
    if (repartition) {
      partition_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats);
    }
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
  while (!done_global) {
    done = 1;

    T8_ASSERT (forest_from->maxlevel_existing >= 0);
    /* Initialize the temp forest to be adapted from forest_from */
    t8_forest_init (&forest_temp);
    /* Update the maximum occurring level */
    forest_temp->maxlevel_existing = forest_from->maxlevel_existing;
    /* Adapt the forest */
    t8_forest_set_adapt (forest_temp, forest_from, t8_forest_balance_adapt,
                         0);
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
      while (count >= num_stats - 2) {
        /* re-allocate memory for stats */
        num_stats += 2;
        adap_stats = T8_REALLOC (adap_stats, sc_statinfo_t, num_stats);
        ghost_stats = T8_REALLOC (ghost_stats, sc_statinfo_t, num_stats);
        if (repartition) {
          partition_stats =
            T8_REALLOC (partition_stats, sc_statinfo_t, num_stats);
        }
      }
      sc_stats_set1 (&adap_stats[count], forest_temp->profile->adapt_runtime,
                     "forest balance: Adapt time");
      if (!repartition) {
        sc_stats_set1 (&ghost_stats[count],
                       forest_temp->profile->ghost_runtime,
                       "forest balance: Ghost time");
      }
    }

    /* Compute the logical and of all process local done values, if this results
     * in 1 then all processes are finished */
    sc_MPI_Allreduce (&done, &done_global, 1, sc_MPI_INT, sc_MPI_LAND,
                      forest->mpicomm);

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
        sc_stats_set1 (&partition_stats[count],
                       forest_partition->profile->partition_runtime,
                       "forest balance: Partition time");
        sc_stats_set1 (&ghost_stats[count],
                       forest_partition->profile->ghost_runtime,
                       "forest balance: Ghost time");
      }

      forest_temp = forest_partition;
      forest_partition = NULL;
    }
    /* Adapt forest_temp in the next round */
    forest_from = forest_temp;
    count++;
  }

  T8_ASSERT (t8_forest_is_balanced (forest_temp));
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
      if (repartition) {
        part_time += partition_stats[i].sum_values;
      }
    }
    sc_stats_set1 (&adap_stats[count], ada_time,
                   "forest balance: Total adapt time");
    sc_stats_set1 (&ghost_stats[count], ghost_time,
                   "forest balance: Total ghost time");
    if (repartition) {
      sc_stats_set1 (&partition_stats[count], part_time,
                     "forest balance: Total partition time");
    }

    /* Compute and print the intermediate stats */
    sc_stats_compute (forest->mpicomm, count + 1, adap_stats);
    sc_stats_compute (forest->mpicomm, count + 1, ghost_stats);
    if (repartition) {
      sc_stats_compute (forest->mpicomm, count + 1, partition_stats);
    }
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count + 1,
                    adap_stats, 1, 1);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count + 1,
                    ghost_stats, 1, 1);
    if (repartition) {
      sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count + 1,
                      partition_stats, 1, 1);
    }
    T8_FREE (adap_stats);
    T8_FREE (ghost_stats);
    if (repartition) {
      T8_FREE (partition_stats);
    }
  }
  /* The forest is now balanced */
  forest->is_balanced = 1;
}

/* Adapt and balance a forest in the same step.
 * We do not create any new elemnts until the final step.
 * We perform these steps:
 *  - Create an array of refinement markers from the adapt function.
 *    Stating for each element whether it will get refined, coarsened or unchanged.
 *  - Iterate through all elements that will be unchanged or coarsened and check whether the balance
 *    condition will break. If so, change their refinement marker appropiately.
 *  - Repeat step 2 until no marker changes.
 *  - Carry out the refinement via the markers.
 */
/* TODO: We currenly recompute the leaf face neighbors in every step.
 *       we can spare this by storing the leaf face neighbors.
 *       Possibly with two arrays: t8_locidx_t *face_neighbors;
 *                                 size_t *offsets_into_face_neighbors;
 *       Or both sc_arrays.
 * TODO: Can we extend this function to also work when forest_from is not balanced?
 *       Does it suffice to also add refinement markers in the first round for all elements
 *       that need to get refined in order to restore balance?
 */
void
t8_forest_balance_and_adapt (t8_forest_t forest)
{
  t8_forest_t         forest_from = forest->set_from;
  /* The forest from must be committed and balanced */
  T8_ASSERT (t8_forest_is_committed (forest_from));
  T8_ASSERT (t8_forest_is_balanced (forest_from));
  /* This technique does not allow for recursive adaptation. */
  T8_ASSERT (forest->set_adapt_recursive == 0);

  /* We assume that we have ghost elements in forest_from. */
  T8_ASSERT (t8_forest_get_ghost_type (forest_from) == T8_GHOST_FACES);
  const t8_locidx_t   num_local_elements =
    t8_forest_get_num_element (forest_from);
  const t8_locidx_t   num_ghosts = t8_forest_get_num_ghosts (forest_from);
  sc_array_t          markers;
  t8_locidx_list_t    elements_that_do_not_refine;
  t8_locidx_list_iterator_t list_iterator;
  /* Keeping track of whether this process is finished and
   * whether all processes are finished. */
  int                 done_locally = 0, done_globally = 0;
  int                 balance_rounds = 0;
  /* The current communicator */
  const sc_MPI_Comm   comm = t8_forest_get_mpicomm (forest_from);
  /* A hash table to store the once computed face neighbors. */
  t8_face_neighbor_hash_table_t *face_neighbor_table;

  t8_global_productionf
    ("Into t8_forest_balance_and_adapt with %lli global elements.\n",
     (long long) t8_forest_get_global_num_elements (forest_from));
  t8_log_indent_push ();

  if (forest->profile != NULL) {
    /* Profiling is enabled, we measure the runtime of this function. */
    forest->profile->balance_runtime -= sc_MPI_Wtime ();
    /* We increment balance_rounds in each round */
    forest->profile->balance_rounds = 0;
  }

  /* Initialize refinement markers for all elements and ghosts */
  sc_array_init_size (&markers, sizeof (short),
                      num_local_elements + num_ghosts);

  /* Initialize a list of the element indices that do not get refined (markes[i] != 1).
   * These are the elements for which we have to check whether they need to get refined
   * (or not coarsened) in order to preserve the balance condition.
   * We need to check these elements repeatedly until nothing changes anymore. */

  /* initialize the linked list */
  t8_locidx_list_init (&elements_that_do_not_refine);

  /* We now create the refinement markers.
   * +1 means refine the element
   *  0 means do nothing
   * -1 means coarsen the family (must be set for all members of the family)
   * 
   * Thus, the value of the markers corresponds to the difference in refinement
   * levels.
   */
  t8_forest_adapt_build_marker_array (forest, &markers,
                                      &elements_that_do_not_refine);

  /* Communicate ghost values for the markers array */
  t8_forest_ghost_exchange_data (forest_from, &markers);

  /* Create a hash table to store the face neighbors of the elements that
   * do not refine. */
  face_neighbor_table = t8_face_neighbor_hash_table_new (forest_from);

  /* TODO: We used these for debugging. Reactivate if needed. Delete eventually. */
#if 0
#ifdef T8_ENABLE_DEBUG
  {
    /* Print vtk and marker array */
    t8_locidx_t         ielement;
    char                output[BUFSIZ + 1] = "";

    if (forest->dimension > 0) {
      t8_forest_write_vtk (forest_from, "t8_DEBUG_forest_From_balance_adapt");
    }

    for (ielement = 0; ielement < num_local_elements + num_ghosts; ++ielement) {
      snprintf (output + strlen (output), BUFSIZ - strlen (output),
                " | %hi |", *(short *) t8_sc_array_index_locidx (&markers,
                                                                 ielement));
    }
    t8_debugf ("[H] Markers for %i element(s): %s\n",
               num_local_elements + num_ghosts, output);
  }
#endif
#endif

  /* We now iterate through the not refined elements and check
   * whether the balance condition will be broken for them.
   * If so we change there marker to refinement or do nothing,
   * whichever is appropriate. 
   * If an element that was not tagged for refinement gets tagged
   * for refinement in the process, we remove it from the
   * list of unrefined elements. */
  /* We have to repeat this iteration, updating the ghost markers
   * everytime in between, until no process updates its markers 
   * anymore. */
  while (!done_globally) {
    int                 changed_a_marker = 0;   /* Keep track if we change a marker. */
    t8_locidx_t         current_tree = 0;       /* The local tree that the current element is in */
    t8_locidx_t         current_tree_offset = 0;        /* The element offset of the current tree */
    t8_locidx_t         next_tree_offset;       /* The element offset of the next tree */
    t8_eclass_scheme_c *eclass_scheme;  /* The element scheme of the current tree */
    t8_eclass_t         tree_class;     /* The element class of the current tree */
    t8_locidx_t         element_in_list;        /* Counts how many elements in the list we already handled */
    /* We set done_locally to 1 and if we need to change any elements
     * we will set it to 0. */
    done_locally = 1;
    /* TODO: If we know that we checked all local elements, are
     *       are done locally, but not globally, we can optimize by
     *       checking only the neighbors of the updated ghost elements. */

    balance_rounds++;

    if (num_local_elements > 0) {
      /* We only need to iterate over the elements if there are any. */
      /* Init next tree_offset as the offset of tree 1
       * Note that this call is only legal if we have at least one element. */
      next_tree_offset =
        current_tree_offset + t8_forest_get_tree_num_elements (forest_from,
                                                               0);
      /* Get the tree's class and element scheme */
      tree_class = t8_forest_get_tree_class (forest_from, current_tree);
      eclass_scheme = t8_forest_get_eclass_scheme (forest_from, tree_class);
      /* Iterate over all elements that do not yet refine */
      for (element_in_list = 0, t8_locidx_list_iterator_init
           (&elements_that_do_not_refine, &list_iterator);
           !t8_locidx_list_iterator_is_end (&list_iterator);
           t8_locidx_list_iterator_next (&list_iterator), ++element_in_list) {
        /* Get the index of the next non-refined element */
        const t8_locidx_t   element_index =
          t8_locidx_list_iterator_get_value (&list_iterator);
        const short         current_marker =
          *(short *) t8_sc_array_index_locidx (&markers, element_index);
        t8_locidx_t         element_index_in_tree;
        int                 element_level;
        int                 element_level_after_adapt;
        t8_element_t       *element;
        t8_locidx_t         iface;
        /* Sometimes we do not need to continue checking the rest of the face
         * neighbors. We keep track of whether to continue or not in this variable. */
        int                 continue_neigh_check = 1;
        t8_face_neighbor_hash_t *face_neighbor_hash_entry;
        /* Flag that we activate when this element unmarked its family
         * for coarsening, so that we do not do it more than once per element. */
        int                 unmarked_family_for_coarsen = 0;

        /* Check whether this element is still in the current tree.
         * If not, update the current tree. */
        if (next_tree_offset <= element_index) {
          /* This element is in one of the next trees */
          /* Find the tree this element is in */
          while (next_tree_offset <= element_index) {
            current_tree++;
            current_tree_offset =
              t8_forest_get_tree_element_offset (forest_from, current_tree);
            next_tree_offset =
              current_tree_offset +
              t8_forest_get_tree_num_elements (forest_from, current_tree);
          }
          T8_ASSERT (current_tree_offset <= element_index
                     && element_index < next_tree_offset);
          /* We can now update the eclass and eclass scheme */
          tree_class = t8_forest_get_tree_class (forest_from, current_tree);
          eclass_scheme =
            t8_forest_get_eclass_scheme (forest_from, tree_class);
        }
        /* Compute the index of the element in this tree */
        element_index_in_tree = element_index - current_tree_offset;
        /* Get the element */
        element =
          t8_forest_get_element_in_tree (forest_from, current_tree,
                                         element_index - current_tree_offset);
        /* Get the level of the element and compute its level after adaptation. */
        element_level = eclass_scheme->t8_element_level (element);
        element_level_after_adapt = element_level + current_marker;
        if (balance_rounds == 1) {
          /* This is the first time that we visit this element. We compute its face
           * neighbors and store them in the hash table.
           * This will also set the face_neighbor_hash_entry to point to the inserted values. */
          int                 retval =
            t8_face_neighbor_hash_table_insert_element (face_neighbor_table,
                                                        current_tree,
                                                        element_index_in_tree,
                                                        1,
                                                        &face_neighbor_hash_entry);
          /* Check that we could insert the element. */
          T8_ASSERT (retval);
        }
        else {
          /* The element must be contained in the hash table, we look for it. */
          face_neighbor_hash_entry =
            t8_face_neighbor_hash_table_lookup (face_neighbor_table,
                                                current_tree,
                                                element_index_in_tree);
          /* The element must have been found. */
          T8_ASSERT (face_neighbor_hash_entry != NULL);
        }
        for (iface = 0;
             iface < face_neighbor_hash_entry->num_faces
             && continue_neigh_check; ++iface) {
          const int           num_face_neighbors =
            face_neighbor_hash_entry->number_of_neighbors[iface];
          int                 ineigh;
          int                 neigh_level, neigh_level_after_adapt,
            level_diff;
          short               neigh_marker;
          const t8_locidx_t  *neighbor_indices =
            face_neighbor_hash_entry->face_neighbor_indices[iface];
          const t8_element_t **face_neighbors = (const t8_element_t **)
            face_neighbor_hash_entry->face_neighbors[iface];

          /* Get this neighbors class and scheme */
          const t8_eclass_t   neigh_class =
            t8_forest_element_neighbor_eclass (forest_from, current_tree,
                                               element,
                                               iface);
          t8_eclass_scheme_c *neigh_scheme =
            t8_forest_get_eclass_scheme (forest_from, neigh_class);

#ifdef T8_ENABLE_DEBUG
          /* In debugging mode check that the face neighbors we get from the 
           * offsets are actually the correct face neighbors */
          {
            t8_element_t      **face_neighbors_debug;
            t8_eclass_scheme_c *neigh_scheme_debug;
            int                 iface_debug;
            int                 num_face_neighbors_debug;
            int                *dual_faces_debug;
            t8_locidx_t        *face_neighbor_indices_debug;
            t8_forest_leaf_face_neighbors (forest_from, current_tree, element,
                                           &face_neighbors_debug, iface,
                                           &dual_faces_debug,
                                           &num_face_neighbors_debug,
                                           &face_neighbor_indices_debug,
                                           &neigh_scheme_debug, 1);
            T8_ASSERT (neigh_scheme == neigh_scheme_debug);
            for (iface_debug = 0; iface_debug < num_face_neighbors;
                 ++iface_debug) {
              T8_ASSERT (!neigh_scheme->t8_element_compare
                         (face_neighbors_debug[iface_debug],
                          face_neighbors[iface_debug]));
            }
            T8_FREE (face_neighbor_indices_debug);
            T8_FREE (dual_faces_debug);
            neigh_scheme->t8_element_destroy (num_face_neighbors_debug,
                                              face_neighbors_debug);
            T8_FREE (face_neighbors_debug);
          }
#endif

          if (face_neighbors != NULL) {
            for (ineigh = 0;
                 ineigh < num_face_neighbors && continue_neigh_check;
                 ++ineigh) {
              t8_locidx_t         neighbor_index;
              neighbor_index = neighbor_indices[ineigh];
              T8_ASSERT (face_neighbors[ineigh] != NULL);

              /* Get the level of this neighbor */
              neigh_level =
                neigh_scheme->t8_element_level (face_neighbors[ineigh]);
              /* Get the refinement marker of this neighbor */
              neigh_marker =
                *(short *) t8_sc_array_index_locidx (&markers,
                                                     neighbor_index);
              /* Compute the level that the neighbor will have after adapting the mesh */
              neigh_level_after_adapt = neigh_level + neigh_marker;
              level_diff =
                neigh_level_after_adapt - element_level_after_adapt;
              if (level_diff > 1) {
                T8_ASSERT (level_diff == 2 || level_diff == 3);
                /* After adaptation the current element's level would be smaller
                 * than the neighbors level.
                 * Thus, the balance condition would be broken.
                 * Depending on the element's marker, we do:
                 *   If the marker is 0:
                 *      - mark this element for adaptation (with 1) and remove from the list.
                 *      - Stop checking its face neighbors.
                 *   If the marker is -1:
                 *      If the level difference is 2:
                 *        - mark this element and all its siblings with 0
                 *        - Continue checking its face neighbors
                 *      If the level difference is 3:
                 *        - mark this element's siblings with 0
                 *        - mark this element with 1 and remove from the list
                 *        - Stop checking its face neighbors
                 */
                if (current_marker == -1) {
                  /* If the marker is -1, we mark all siblings with 0 */
                  /* If the level difference is 2, then nothing else is to be done. */

                  /* To identify the indices of this element's siblings,
                   * we need its sibid and the number of siblings. */
                  if (!unmarked_family_for_coarsen) {
                    /* We may have already done this for a previous neighbor of this
                     * element, and this can skip this if we did already unmark the family. */
                    const t8_locidx_t   childid =
                      eclass_scheme->t8_element_child_id (element);
                    const t8_locidx_t   num_siblings =
                      eclass_scheme->t8_element_num_siblings (element);
                    t8_locidx_t         isib;
                    /* The index of the first sibling in the family is this neighbor's index minus its child id */
                    const t8_locidx_t   first_sibling =
                      element_index - childid;
                    /* NOTE: We know that the whole family are local elements, since our 
                     *       marker is -1. Thus the current elements' whole family is local and was
                     *       marker for coarsening. */
#ifdef T8_ENABLE_DEBUG
                    /* Check whether first_sibling is indeed the first sibling */
                    const t8_locidx_t   first_sibling_in_tree =
                      first_sibling - current_tree_offset;
                    const t8_element_t *first_sib =
                      t8_forest_get_element_in_tree (forest_from,
                                                     current_tree,
                                                     first_sibling_in_tree);
                    t8_element_t       *test_sib, *test_parent;
                    T8_ASSERT (first_sibling_in_tree ==
                               element_index_in_tree - childid);

                    T8_ASSERT (eclass_scheme->t8_element_child_id (first_sib)
                               == 0);
                    /* Construct the sibling and check whether it matches first_sib */
                    eclass_scheme->t8_element_new (1, &test_sib);
                    eclass_scheme->t8_element_new (1, &test_parent);
                    eclass_scheme->t8_element_parent (element, test_parent);
                    eclass_scheme->t8_element_child (test_parent, 0,
                                                     test_sib);
                    /* Check that fist_sib and test_sib are the same element. */
                    T8_ASSERT (!eclass_scheme->t8_element_compare
                               (test_sib, first_sib));
                    eclass_scheme->t8_element_destroy (1, &test_sib);
                    eclass_scheme->t8_element_destroy (1, &test_parent);
#endif
                    /* Mark all siblings with 0 (do not coarsen or refine). */
                    for (isib = 0; isib < num_siblings; ++isib) {
                      short              *sibs_marker =
                        (short *) t8_sc_array_index_locidx (&markers,
                                                            first_sibling +
                                                            isib);
                      /* The family must have been marked for coarsening. */
                      T8_ASSERT (*sibs_marker == -1);
                      *sibs_marker = 0;
                      changed_a_marker = 1;
                    }
                    unmarked_family_for_coarsen = 1;
                  }
                }               /* End if (current_marker == -1) */
                if (level_diff == 3 || current_marker == 0) {
                  /* If the level difference is 3 or the marker is 0 (and level difference is 2),
                   * we mark this element for refinement and stop checking its neighbor. */
                  T8_ASSERT ((current_marker == -1 && level_diff == 3)
                             || (current_marker == 0 && level_diff == 2));
                  /* Mark this element for refinement */
                  *(short *) t8_sc_array_index_locidx (&markers,
                                                       element_index) = 1;
                  changed_a_marker = 1;
                  /* Remove this element from the list */
                  t8_locidx_list_iterator_remove_entry (&list_iterator);
                  /* Do not continue checking the neighbors. */
                  continue_neigh_check = 0;
                }
              }                 /* End if (level_diff > 1) */
            }                   /* End neighbor for-loop */
          }                     /* End face_neighbors != NULL */
        }                       /* End face neighbor loop */
      }                         /* End element iteration */
    }                           /* if (num_local_elements > 0) ends here */
    if (!changed_a_marker) {
      /* We did not change anything and thus, this process is done. */
      done_locally = 1;
    }
    else {
      /* We changed a marker and thus this process is not done yet. */
      done_locally = 0;
    }
    /* Update the done_globally flag */
    sc_MPI_Allreduce (&done_locally, &done_globally, 1, sc_MPI_INT,
                      sc_MPI_BOR, comm);
  }

  /* We are now done with the iterations and updated all markers. */

  /* clean up the list */
  t8_locidx_list_reset (&elements_that_do_not_refine);

  /* Destroy the face neighbor hash table. */
  t8_face_neighbor_hash_table_destroy (face_neighbor_table);

  /* We can now finally refine and coarsen the elements according to the markers.
   * In order to do so, we use the existing adapt functions and have the 
   * t8_forest_adapt_marker_array_callback function as adapt callback. */
  forest->set_adapt_fn = t8_forest_adapt_marker_array_callback;

  /* Temporarily store the current user data in order to set
   * the markers array as user data for adaptation. */
  void               *user_pointer_safe = t8_forest_get_user_data (forest);
  t8_forest_set_user_data (forest, &markers);
  /* Adapt the forest */
  t8_forest_adapt (forest);
  /* Restore the user data */
  t8_forest_set_user_data (forest, user_pointer_safe);
  /* Clean up the markers */
  sc_array_reset (&markers);
  /* The forest is now balanced */
  forest->is_balanced = 1;

  if (forest->profile != NULL) {
    /* If profiling is enabled, we measure the runtime of this function */
    forest->profile->balance_runtime += sc_MPI_Wtime ();
    /* And set the number of rounds */
    forest->profile->balance_rounds = balance_rounds;
  }

  t8_log_indent_pop ();

  t8_global_productionf
    ("Done t8_forest_balance_and_adapt with %lli global elements.\n",
     (long long) forest->global_num_elements);

  /* TODO: Write a test for this function.
   *       - First test: a mesh that does not need to be balanced,
   *                     verify that this function does the same as adapt.
   *       - second test: a mesh that does need to be balanced,
   *                     verify that this function  produces the same forest
   *                     as doing first adapt and then balance.
   */
  /* TODO: Add a new timer to the profile that times the runtime of this function. */
}

/* Check whether the local elements of a forest are balanced. */
int
t8_forest_is_balanced (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  t8_locidx_t         num_trees, num_elements;
  t8_locidx_t         itree, ielem;
  t8_eclass_scheme_c *ts;
  void               *data_temp;
  int                 dummy_int;

  T8_ASSERT (t8_forest_is_committed (forest));

#ifndef T8_ENABLE_DEBUG
  /* If the is_balanced flag is set, we assume that the forest is balanced.
   * In debugging mode, we check anyways */
  if (forest->is_balanced) {
    return 1;
  }
#endif

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
    ts =
      t8_forest_get_eclass_scheme (forest,
                                   t8_forest_get_tree_class (forest, itree));
    /* Iterate over all elements of this tree */
    for (ielem = 0; ielem < num_elements; ielem++) {
      const t8_element_t *element =
        t8_forest_get_element_in_tree (forest, itree, ielem);
      /* Test if this element would need to be refined in the balance step.
       * If so, the forest is not balanced locally. */
      if (t8_forest_balance_adapt
          (forest, forest, itree, ielem, ts, 1, &element)) {
        /* The forest is not balanced */
        forest->set_from = forest_from;
        forest->t8code_data = data_temp;
        /* Check whether the is_balanced flag was set correctly */
        T8_ASSERT (forest->is_balanced == 0);
        return 0;
      }
    }
  }
  /* The forest is balanced */
  forest->set_from = forest_from;
  forest->t8code_data = data_temp;
  /* Set the is_balanced flag */
  forest->is_balanced = 1;
  return 1;
}

T8_EXTERN_C_END ();
