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

#include "t8_eclass.h"
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_conformal_quad/t8_transition_conformal_quad_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This is the conformal transition refine function for the 2D quad scheme. 
 * It will return a values p>1 in order to exchange the current element with a transition cell of type p,
 * which is defined in the conformal_quad scheme. */
int
t8_forest_transition_conformal_quad (t8_forest_t forest,
                                     t8_forest_t forest_from,
                                     t8_locidx_t ltree_id,
                                     t8_locidx_t lelement_id,
                                     t8_eclass_scheme_c *ts,
                                     const int is_family,
                                     int num_elements,
                                     t8_element_t *elements[])
{
  int                 iface, num_faces, neigh_face, transition_type = 0;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  t8_element_t       *element = elements[0], **face_neighbor;

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
}                               /* end of t8_forest_transition_conformal_quad */

int
t8_forest_transition_conformal_hex (t8_forest_t forest,
                                     t8_forest_t forest_from,
                                     t8_locidx_t ltree_id,
                                     t8_locidx_t lelement_id,
                                     t8_eclass_scheme_c *ts,
                                     const int is_family,
                                     int num_elements,
                                     t8_element_t *elements[])
{
  int                 iface, num_faces, neigh_face, transition_type = 0;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  t8_element_t       *element = elements[0], **face_neighbor;

  /* Hanging faces can only exist at non-maxlevel elements */
  if (forest_from->maxlevel_existing <= 0 ||
      ts->t8_element_level (element) < forest_from->maxlevel) {

    num_faces = ts->t8_element_num_faces (element);

    /* TODO: Update this comment to HEX. */
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
      // t8_element_t Array mit einem Element
      face_neighbor = T8_ALLOC (t8_element_t *, 1);

      neigh_scheme->t8_element_new (1, face_neighbor);

      /* Compute the virtual face neighbor of element at this face */
      neighbor_tree = t8_forest_element_face_neighbor (forest_from, ltree_id,
                                                       element,
                                                       face_neighbor[0],
                                                       neigh_scheme,
                                                       iface, &neigh_face);

      /* TODO: Update this code block for hex / sub-pyramids. */
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
    else if (transition_type == 63) {   /* Six hanging faces in this case */
      return 1;
    }
    else {                      /* use a transition cell of subelements and add 1 to every type, to avoid refine = 1 */
      return transition_type + 1;
    }
  }
  return 0;                     /* if elem has maxlevel then keep it unchanged since there will never be hanging faces */
}

/* This is the entry function for all transition schemes, called bei forest_adapt.
 * The eclass of the current element hands off to the specific refine implementation above.
 * Other schemes, for other eclasses can easily be added. */
int
t8_forest_transition_entry (t8_forest_t forest,
                            t8_forest_t forest_from,
                            t8_locidx_t ltree_id,
                            t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts,
                            const int is_family,
                            int num_elements, t8_element_t *elements[])
{
  //T8_ASSERT (forest->set_subelements == 1);
  T8_ASSERT (forest->is_transitioned == 0
             && forest->set_from->is_transitioned == 0);

  /* TODO: there may be a better way for this than using switch statements and int functions */

  /* the current element decides over the refine function. Normally, this function returns one fixed value per 
   * element scheme, but note that it would also possible to switch the refine function within one tree. */

  switch (ts->t8_element_get_transition_refine_identifier ()) {
  case T8_TRANSITION_CONFORMAL_ZERO_REFINE_FUNCTION:
    /* if no transition scheme is implemented for the given element/tree, 
     * then return zero and keep the current element unchanged. 
     * The default_common implementation of the above function returns 0. */
    return 0;
  case T8_TRANSITION_CONFORMAL_QUAD_REFINE_FUNCTION:
    return t8_forest_transition_conformal_quad (forest, forest_from, ltree_id,
                                                lelement_id, ts, is_family,
                                                num_elements, elements);
  case T8_TRANSITION_CONFORMAL_HEX_REFINE_FUNCTION:                                           
    return t8_forest_transition_conformal_hex (forest, forest_from, ltree_id,
                                                lelement_id, ts, is_family,
                                                num_elements, elements);
  default:
    SC_ABORT
      ("The given eclass scheme must specify a valid transition refine function.");
  }
}                               /* end of t8_forest_transition_entry */

void
t8_forest_transition (t8_forest_t forest)
{
  T8_ASSERT (forest->set_subelements == 1);
  T8_ASSERT (forest->is_transitioned == 0
             && forest->set_from->is_transitioned == 0);
  /* in the following, we will call forest_adapt to with the transition 
   * refinement function in order to transition the forest. The refinement is then
   * based on forest->set_from, which must be balanced. */

  t8_global_productionf ("Into t8_forest_transition.\n");

  /* Set ghost layers of all processes */
  if (forest->set_from->ghosts == NULL) {
    forest->set_from->ghost_type = T8_GHOST_FACES;
    t8_forest_ghost_create (forest->set_from);
  }

  forest->set_adapt_fn = t8_forest_transition_entry;
  forest->set_adapt_recursive = 0;
  t8_forest_copy_trees (forest, forest->set_from, 0);
  t8_forest_adapt (forest);

  t8_global_productionf ("Done t8_forest_transition.\n");
}                               /* end of t8_forest_transition */


/* This is the entry function for all untransition a forest, called bei forest_adapt 
 * in t8_forest_untransition.
 */
int
t8_forest_untransition_entry (t8_forest_t forest,
                            t8_forest_t forest_from,
                            t8_locidx_t ltree_id,
                            t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts,
                            const int is_family,
                            int num_elements, t8_element_t *elements[])
{
  T8_ASSERT (forest->is_transitioned == 1
             && forest->set_from->is_transitioned == 1);

  //Iterate through elements array and if an element is a subelement, return -1 for 
  //coarsening this element in forest_adapt. All other elements remain unchanged. 
  t8_element_t *element = elements[0];

  if(ts->t8_element_is_subelement(element)){
    return -1;
  }
  else{
    return 0;
  }

}  


void
t8_forest_untransition (t8_forest_t forest)
{
  T8_ASSERT (forest->is_transitioned == 1
             && forest->set_from->is_transitioned == 1);
  /* In the following, we will call forest_adapt to coarsen all transition cells
   * back to their regular parent element. So, forest->set_from is a transitioned forest and 
   * forest will be non-transitioned. */

  t8_global_productionf ("Into t8_forest_untransition.\n");

  forest->set_adapt_fn = t8_forest_untransition_entry;
  forest->set_adapt_recursive = 0;
  t8_forest_copy_trees (forest, forest->set_from, 0);
  t8_forest_adapt (forest);
  forest->is_transitioned = 0;
  t8_global_productionf ("Done t8_forest_untransition.\n");
}                               

/* Test whether the forest is transitioned.
 * Note 1) We allow non-committed forests in this implementation since this check is used in forest_commit() 
 * Note 2) This test does only check whether there is at least one subelement in the forest. 
 *         It is not tested whether the forest is conformal or has any other property. */
int
t8_forest_is_transitioned (t8_forest_t forest)
{
  t8_eclass_scheme_c *tscheme;
  t8_locidx_t         tree_count, num_trees;
  t8_tree_t           current_tree;
  t8_element_array_t *telements;
  t8_locidx_t         elem_count, num_elems;

  /* iterate through the forest and check for subelements */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (tree_count = 0; tree_count < num_trees; tree_count++) {
    current_tree = t8_forest_get_tree (forest, tree_count);
    telements = &current_tree->elements;
    num_elems = (t8_locidx_t) t8_element_array_get_count (telements);
    for (elem_count = 0; elem_count < num_elems; elem_count++) {
      t8_element_t       *current_element =
        t8_element_array_index_locidx (telements, elem_count);
      tscheme = forest->scheme_cxx->eclass_schemes[current_tree->eclass];
      if (tscheme->t8_element_is_subelement (current_element)) {
        /* subelement found -> return true */
        return 1;
      }
    }
  }

  /* only return false if there is no subelement in the forest */
  return 0;
}                               /* end of t8_forest_is_transitioned */

T8_EXTERN_C_END ();
