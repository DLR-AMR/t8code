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

#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements/t8_subelements_quad_cxx.hxx>

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
t8_forest_remove_hanging_faces_adapt (t8_forest_t forest,
                                      t8_forest_t forest_from,
                                      t8_locidx_t ltree_id,
                                      t8_locidx_t lelement_id,
                                      t8_eclass_scheme_c * ts,
                                      int num_elements,
                                      t8_element_t * elements[])
{
  /* this function determines whether a subelement (and which type) should be used for a specific 
   * element of the forest. The return value depends on the neighbor elements and their 
   * refinement level. */
  int                 num_faces, iface, *dual_faces, num_neighbors,
    subelement_type = 0;
  const t8_element_t *current_element;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;

  /* In order to remove al hanging nodes, the forest must be balanced */
  T8_ASSERT (t8_forest_is_balanced (forest_from));

  current_element =
    t8_forest_get_element_in_tree (forest_from, ltree_id, lelement_id);
  num_faces = ts->t8_element_num_faces (current_element);

  /* We use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
   * Every face has a flag parameter, wich is set to 1, if there is a neighbour with a higher level 
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
   * Note, that this procedure is independent of the eclass (we only show an example for the quad scheme). 
   * Each neighbour-structure will lead to a unique binary code. 
   * Within the element scheme of the given eclass, this binary code is used to construct the right subelement type,
   * in order to remove hanging nodes from the mesh. */

  for (iface = 0; iface < num_faces; iface++) {
    /* We are interested in the number of neighbours of a given element and a given face of the element. */
    t8_forest_leaf_face_neighbors (forest_from, ltree_id, current_element,
                                   &neighbor_leafs, iface, &dual_faces,
                                   &num_neighbors, &element_indices,
                                   &neigh_scheme, 1);
    /* If the number of neighbours of a face is higher than 1, then we know that there must be a hanging node. */

    /* This procedure determines the decimal value of the binary representation of the neighbour structure. 
     *     
     *    binary:   |    1    |    0    |    0    |    1    |
     *    decimal:     1*2^3  +  0*2^2  +  0*2^1  +  1*2^0  =  9
     *  
     */
    int                 coefficient;
    if (num_neighbors == 0 || num_neighbors == 1) {
      coefficient = 0;
    }
    else {
      coefficient = 1;
    }
    subelement_type += coefficient * 1 << ((num_faces - 1) - iface);

    if (num_neighbors > 0) {    /* free memory */
      neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

      T8_FREE (element_indices);
      T8_FREE (neighbor_leafs);
      T8_FREE (dual_faces);
    }
  }

  /* Print some useful output in debug mode */
  t8_debugf ("element id: %i    subelement type: %i\n", lelement_id,
             subelement_type);

  /* returning the right subelement types */
  if (subelement_type == 0) {   /* in this case, there are no hanging nodes and we do not need to do anything */
    return 0;
  }
  else {                        /* use subelements and add 1 to every type, to avoid refine = 1 */
    return subelement_type + 1;
  }
}

void
t8_forest_remove_hanging_faces (t8_forest_t forest)
{
  t8_global_productionf ("Into t8_forest_subelements.\n");

  forest->set_adapt_fn = t8_forest_remove_hanging_faces_adapt;
  forest->set_adapt_recursive = 0;
  t8_forest_copy_trees (forest, forest->set_from, 0);
  t8_forest_adapt (forest);

  t8_global_productionf ("Done t8_forest_subelements.\n");
}

#if 0
/* Check whether all hanging nodes are eliminated. */
int
t8_forest_hanging_faces_removed (t8_forest_t forest)
{
  /* TODO: implement this function in the future */
  return 0;
}
#endif

T8_EXTERN_C_END ();
