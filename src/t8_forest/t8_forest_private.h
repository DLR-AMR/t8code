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

/** \file t8_forest_private.h
 * We define routines for a forest of elements that are not part
 * of the official t8_forest.h interface but used internally.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_private */

#ifndef T8_FOREST_PRIVATE_H
#define T8_FOREST_PRIVATE_H

#include <t8.h>
#include <t8_forest.h>

T8_EXTERN_C_BEGIN ();

/* TODO: document */

/* For each tree in a forest compute its first and last descendant */
void                t8_forest_compute_desc (t8_forest_t forest);

/* Create the elements on this process given a uniform partition
 * of the coarse mesh. */
void                t8_forest_populate (t8_forest_t forest);

/* return nonzero if the first tree of a forest is shared with a smaller
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 */
int                 t8_forest_first_tree_shared (t8_forest_t forest);

/* Allocate memory for trees and set their values as in from.
 * For each tree allocate enough element memory to fit the elements of from.
 * If copy_elements is true, copy the elements of from into the element memory.
 */
void                t8_forest_copy_trees (t8_forest_t forest,
                                          t8_forest_t from,
                                          int copy_elements);

/** Return an element of a tree.
 * \param [in]  tree  The tree.
 * \param [in]  elem_in_tree The index of the element within the tree.
 * \return      Returns the elemen with index \a elem_in_tree of the
 *              element array of \a tree.
 */
t8_element_t       *t8_forest_get_tree_element (t8_tree_t tree,
                                                t8_locidx_t elem_in_tree);

/** Find the owner process of a given element.
 * \param [in]    forest  The forest.
 * \param [in]    gtreeid The global id of the tree in which the element lies.
 * \param [in]    element The element to look for.
 * \param [in]    eclass  The element class of the tree \a gtreeid.
 * \return                The mpirank of the process that owns \a element.
 * \note The element must exist in the forest.
 */
int                 t8_forest_element_find_owner (t8_forest_t forest,
                                                  t8_gloidx_t gtreeid,
                                                  t8_element_t * element,
                                                  t8_eclass_t eclass);

/** Construct all face neighbors of half size of a given element.
 * \param [in]    forest The forest.
 * \param [in]    ltreeid The local tree id of the tree in which the element is.
 * \param [in]    elem    The element of which to construct the neighbors.
 * \param [in,out] neighs An array of allocated elements of the correct element class.
 *                        On output the face neighbors of \a elem across \a face of one
 *                        bigger refinement level are stored.
 * \param [in]    face    The number of the face of \a elem.
 * \param [in]    num_neighs The number of allocated element in \a neighs. Must match the
 *                        number of face neighbors of one bigger refinement level.
 * \return                The local id of the tree in which the neighbors are.
 */
t8_locidx_t         t8_forest_element_half_face_neighbors (t8_forest_t forest,
                                                           t8_locidx_t
                                                           ltreeid,
                                                           const t8_element_t
                                                           * elem,
                                                           t8_element_t *
                                                           neighs[], int face,
                                                           int num_neighs);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_PRIVATE_H! */
