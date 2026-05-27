/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_quads_hanging_nodes.cxx
 * This is an example to demonstrate hanging node resolution for quads. 
 */

#include "t8_forest_subelement.hxx"
#include <t8.h>
#include "t8_forest_general.h"
#include "t8_forest_types.h"
#include "t8_forest_private.h"
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_scheme_implementation.hxx>
#include <t8_element/t8_element.h>
#include "t8_forest_adapt.h"

/* This is the adapt function, called for each element in a balanced forest during transition.
 * We refine an element into a suitable transition cell if it has at most one hanging face */
int
discard_subelements_callback ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                              [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class,
                              [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                              [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                              t8_element_t *elements[])
{
  const t8_subelementquad_scheme *subelem_scheme = (const t8_subelementquad_scheme *) scheme;
  if (subelem_scheme->element_is_subelement (elements[0])) {
    return -1;
  }
  return 0;
}

/* This is the adapt function, called for each element in a balanced forest during transition.
 * We refine an element into a suitable transition cell if it has at most one hanging face */
int
t8_remove_hanging_nodes_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                  [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                  const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                                  [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int subelement_type = 0;
  /* We use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
     * Every face has a flag parameter, which is set to 1, if there is a neighbour with a higher level 
     * and to 0, if the level of the neighbour is at most the level of the element.   
     *             
     *              f0                         1
     *        x - - x - - x              x - - x - - x       
     *        |           |              | \   |   / |
     *        |           |              |   \ | /   |                                                            | f3 | f2 | f1 | f0 |
     *    f3  x           | f2   -->   1 x - - x     | 0   -->   binary code (according to the face enumeration): |  1 |  0 |  0 |  1 | = 9 in base 10  
     *        |           |              |   /   \   |
     *        | elem      |              | /       \ |
     *        x - - - - - x              x - - - - - x
     *              f1                         0 
     *                      
     */
  const int num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  for (int iface = 0; iface < num_faces; iface++) {
    const t8_element_t **neighbors; /**< Neighboring elements. */
    int *dual_faces_internal;       /**< Face indices of the neighbor elements. */
    int num_neighbors;              /**< Number of neighboring elements. */
    t8_locidx_t *neighids;          /**< Neighboring elements ids. */
    t8_eclass_t neigh_class;        /**< Neighboring elements tree class. */

    t8_forest_leaf_face_neighbors (forest_from, which_tree, elements[0], &neighbors, iface, &dual_faces_internal,
                                   &num_neighbors, &neighids, &neigh_class);

    if (num_neighbors > 1) {
      subelement_type += 1 << ((num_faces - 1) - iface);
    }
    /* clean-up */
    if (num_neighbors > 0) {
      // Free allocated memory.
      T8_FREE (neighbors);
      T8_FREE (dual_faces_internal);
      T8_FREE (neighids);
    }
  }

  /* returning the right subelement types */
  if (subelement_type == 0) { /* in this case, there are no hanging nodes and we do not need to do anything */
    return 0;
  }
  else if (subelement_type == 15) { /* Normal 1:8 refinement */
    return 1;
  }
  else { /* use subelements and add 1 to every type, to avoid refine = 1 */
    return subelement_type + 1;
  }
}

/** Adapt forest according to callback. */
bool
t8_forest_has_subelements (t8_forest_t forest)
{
  if (t8_eclass_scheme_is_subelement (t8_forest_get_scheme (forest), T8_ECLASS_QUAD)) {
    return false;
  }
  const t8_subelementquad_scheme *scheme = (const t8_subelementquad_scheme *) t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {

    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielem) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      if (scheme->element_is_subelement (elem)) {
        return true;
      }
    }
  }
  return false;
}

/** Adapt forest according to callback. */
void
t8_forest_discard_subelements (t8_forest_t forest)
{
  if (!t8_forest_has_subelements (forest)) {
    return;
  }
  forest = t8_forest_new_adapt (forest, discard_subelements_callback, 0, 0, NULL);
}

/** Adapt forest according to callback. */
void
t8_forest_remove_hanging_nodes (t8_forest_t forest)
{

  t8_global_productionf ("Into t8_forest_remove_hanging_nodes.\n");
  forest = t8_forest_new_adapt (forest, t8_remove_hanging_nodes_callback, 0, 0, NULL);
  t8_global_productionf ("Done t8_forest_remove_hanging_nodes.\n");
}
