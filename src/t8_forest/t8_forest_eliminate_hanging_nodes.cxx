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
#include <t8_forest/t8_forest_eliminate_hanging_nodes.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
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
t8_forest_balance_adapt (t8_forest_t forest, t8_forest_t forest_from,
                         t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         int num_elements, t8_element_t * elements[])
{
  int                *pdone, iface, num_faces, num_half_neighbors, ineigh;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  t8_element_t       *element = elements[0], **half_neighbors;

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
t8_forest_eliminate_hanging_nodes (t8_forest_t forest)
{
  t8_forest_t         forest_temp, forest_from, forest_partition;
  int                 done = 0, done_global = 0;
  int                 count = 0, num_stats, i;
  double              ada_time, ghost_time, part_time;
  sc_statinfo_t      *adap_stats, *ghost_stats, *partition_stats;

  t8_global_productionf
    ("Into t8_forest_eliminate_hanging_nodes with %lli global elements.\n",
     (long long) t8_forest_get_global_num_elements (forest->set_from));
  t8_log_indent_push ();
}

/* Check whether all hanging nodes are eliminated. */
int
t8_forest_hanging_nodes_eliminated (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  t8_locidx_t         num_trees, num_elements;
  t8_locidx_t         itree, ielem;
  t8_element_t       *element;
  t8_eclass_scheme_c *ts;
  void               *data_temp;
  int                 dummy_int;

  T8_ASSERT (t8_forest_is_committed (forest));

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
      element = t8_forest_get_element_in_tree (forest, itree, ielem);
      /* Test if this element would need to be refined in the balance step.
       * If so, the forest is not balanced locally. */
      if (t8_forest_balance_adapt
          (forest, forest, itree, ielem, ts, 1, &element)) {
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
