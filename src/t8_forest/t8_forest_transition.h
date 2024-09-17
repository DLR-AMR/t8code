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

/** \file t8_forest_transition.h
 * We define the eliminate_hanging_nodes routine to transform a 2:1 balanced, nonconformal forest
 * into a conformal forest. The routine relies on a 2D quad-scheme that has been balanced, such that 
 * there is a 2:1 balance between all elements.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_balance */

#ifndef T8_FOREST_TRANSITION_H
#define T8_FOREST_TRANSITION_H

#include <t8.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();

/**  In this function, use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
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
 * Note, that this procedure is independent of the eclass (we only show an example for the quad scheme). 
 * Each neighbour-structure will lead to a unique binary code. 
 * Within the element scheme of the given eclass, this binary code is used to construct the right subelement type,
 * in order to remove hanging nodes from the mesh.

 * \param [in] forest       the forest to which the new elements belong
 * \param [in] forest_from  the forest that is adapted.
 * \param [in] which_tree   the local tree containing \a elements
 * \param [in] lelement_id  the local element id in \a forest_old in the tree of the current element
 * \param [in] ts           the eclass scheme of the tree
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements the number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 * \return transition type + 1 (to avoid return value 1)
 */

int
t8_forest_transition_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts, int num_elements, t8_element_t *elements[]);

/** This function is the point of entry for the hanging-faces-removing subelements. 
 * In this function, the corresponding callback function for the refine value in forest_adapt is set
 * and the forest is adapted.
 * \param [in,out] forest       The input forest consists of hanging nodes on faces- 
 *                              The output forest contains transition cells and no hanging nodes on faces. 
 */
void
t8_forest_transition (t8_forest_t forest);

/** This function is the point of entry in order to untransition a forest before adapting it. 
 * Every transition cell gets coarsened back to its parent element.
 * \param [in,out] forest       The input forest consists of subelements- 
 *                              The output forest does not. 
 */

void
t8_forest_untransition (t8_forest_t forest);

/** Check whether the forest is transitioned, meaning that subelements exist. 
 * \param [in] forest       The input forest that may or may not contains subelements
 * \return true if forest contains subelements or false if not. 
 *                             
 */

int
t8_forest_is_transitioned (t8_forest_t forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_TRANSITION_H! */
