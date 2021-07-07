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
#include <t8_schemes/t8_new_feature/t8_subelements/t8_subelements_quad_cxx.hxx>

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

int
t8_forest_subelements_adapt (t8_forest_t forest, t8_forest_t forest_from,
                             t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                             t8_eclass_scheme_c * ts,
                             int num_elements, t8_element_t * elements[])
{
  /* this function determines whether a subelement (and which type) should be used for a specific 
   * element of the forest. The return value depends on the neighbor elements and their 
   * refinement level. */
  int                num_faces, iface, *dual_faces, num_neighbors, subelement_type = 0;
  const t8_element_t *current_element;
  t8_element_t       **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;

  current_element = t8_forest_get_element_in_tree (forest_from, ltree_id, lelement_id);

  num_faces = ts->t8_element_num_faces (current_element);

  for (iface = 0; iface < num_faces; iface++) {
    t8_forest_leaf_face_neighbors (forest_from, ltree_id, current_element, &neighbor_leafs,
                                     iface, &dual_faces, &num_neighbors,
                                     &element_indices, &neigh_scheme, 1); 
    if (num_neighbors > 1) {
      subelement_type = subelement_type * 2 + 1;
    }
    else {
      subelement_type = subelement_type *2;
    }
    if (num_neighbors > 0) {
        neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

        T8_FREE (element_indices);
        T8_FREE (neighbor_leafs);
        T8_FREE (dual_faces);
    }
  }
  if (subelement_type == 0) {
    return 0;
  }
  else if (subelement_type == 15) {
    /* if all four neihbors are refined, then use standard refinement for the given element */ 
    return 1;
  }
  else {
    /* avoid refine = 1, since this value is for the recursive refinement scheme */
    return subelement_type + 1;
  }
}

#if 0
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
#endif

void
t8_forest_subelements (t8_forest_t forest)
{
  t8_global_productionf ("Into t8_forest_subelements.\n");

  forest->set_adapt_fn = t8_forest_subelements_adapt;
  forest->set_adapt_recursive = 0;
  t8_forest_copy_trees (forest, forest->set_from, 0);
  t8_forest_adapt (forest);

  t8_global_productionf ("Done t8_forest_subelements.\n");
}

#if 0
/* Check whether all hanging nodes are eliminated. */
int
t8_forest_subelements_used (t8_forest_t forest)
{
  /* NOTE do something */
  return 0;
}
#endif

T8_EXTERN_C_END ();
