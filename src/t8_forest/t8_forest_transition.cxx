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

/* Description:
 * In this file, we define the call-back function that is used to construct transition cells.
 */

#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_quads_transition/t8_transition/t8_transition_quad_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This is the adapt function, called for each element in a balanced forest during transition.
 * We refine an element into a suitable transition cell if it has at most one hanging face */
int
t8_forest_subelement_adapt_remove_hanging_faces (t8_forest_t forest,
                            t8_forest_t forest_from,
                            t8_locidx_t ltree_id,
                            t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts,
                            const int is_family,
                            int num_elements, t8_element_t *elements[])
{
  int                 iface, num_faces, neigh_face, transition_type = 0;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  t8_element_t       *element = elements[0], **face_neighbor;

#if 0                           /* use this to refine all elements of a (uniform) mesh via transition type 16 for tests */
  return 16;
#endif

  /* hanging faces can only exist at non-maxlevel elements */
  if (forest_from->maxlevel_existing <= 0 ||
      ts->t8_element_level (element) < forest_from->maxlevel) {

    num_faces = ts->t8_element_num_faces (element);

    /* We use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
     * Every face has a flag parameter, wich is set to 1, if there is a neighbor with a higher level 
     * and to 0, if the level of the neighbor is at most the level of the element.   
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
     * Each neighbor-structure will lead to a unique binary code. 
     * Within the element scheme of the given eclass, this binary code is used to construct the right subelement type,
     * in order to remove hanging nodes from the mesh. */

    for (iface = 0; iface < num_faces; iface++) {
      /* Get the element class and scheme of the face neighbor */
      neigh_class = t8_forest_element_neighbor_eclass (forest_from,
                                                       ltree_id, element,
                                                       iface);

      neigh_scheme = t8_forest_get_eclass_scheme (forest_from, neigh_class);

      /* Allocate memory for the virtual face neighbor */
      face_neighbor = T8_ALLOC (t8_element_t *, 1);

      neigh_scheme->t8_element_new (1, face_neighbor);

      /* Compute the virtual face neighbor of element at this face */
      neighbor_tree = t8_forest_element_face_neighbor (forest_from, ltree_id,
                                                       element,
                                                       face_neighbor[0],
                                                       neigh_scheme,
                                                       iface, &neigh_face);

      if (neighbor_tree >= 0) {
        if (t8_forest_element_has_leaf_desc (forest_from, neighbor_tree,
                                             face_neighbor[0],
                                             neigh_scheme)) {
          /* Compute transition type as the decimal represenation of the binary concatenation */
          transition_type += 1 << ((num_faces - 1) - iface);
        }
      }
      /* clean-up */
      neigh_scheme->t8_element_destroy (1, face_neighbor);
      T8_FREE (face_neighbor);
    }

    /* returning the right subelement types */
    if (transition_type == 0) { /* no hanging faces in this case */
      return 0;
    }
    else if (transition_type == 15) {   /* four hanging faces in this case */
      return 1;
    }
    else {                      /* use a transition cell of subelements and add 1 to every type, to avoid refine = 1 */
      return transition_type + 1;
    }
  }
  return 0;                     /* if elem has maxlevel then keep it unchanged since there will never be hanging faces */
}

void
t8_forest_transition (t8_forest_t forest)
{
  t8_global_productionf ("Into t8_forest_transition.\n");

  /* Set ghost layers of all processes in order to find hanging faces over process-boundaries */
  if (forest->set_from->ghosts == NULL) {
    forest->set_from->ghost_type = T8_GHOST_FACES;
    t8_forest_ghost_create_topdown (forest->set_from);
  }

  forest->set_adapt_fn = t8_forest_subelement_adapt_remove_hanging_faces;
  forest->set_adapt_recursive = 0;
  t8_forest_copy_trees (forest, forest->set_from, 0);
  t8_forest_adapt (forest);

  t8_global_productionf ("Done t8_forest_transition.\n");
}

/* TODO: it would be better to have a proper test if all hanging nodes are removed */
int
t8_forest_is_conformal (t8_forest_t forest)
{
  t8_global_productionf ("Warning: This test is not implemented.\n");
  return 0;
}

T8_EXTERN_C_END ();
