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

/** \file t8_forest_element_face_neighbor.h
 * Functions to query for the same level face neighbors of a forest element or ghost.
 * \ref t8_forest_element_neighbor_eclass
 * \ref t8_forest_element_face_neighbor
 * \ref t8_forest_leaf_face_orientation
 */

#ifndef T8_FOREST_ELEMENT_FACE_NEIGHBOR_H
#define T8_FOREST_ELEMENT_FACE_NEIGHBOR_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();

/** Return the eclass of the tree in which a face neighbor of a given element or ghost
 * lies.
 * \param [in]      forest      A committed forest.
 * \param [in]      ltreeid     The local tree or ghost tree in which the element lies. 0 <= \a ltreeid < num_local_trees + num_ghost_trees
 * \param [in]      elem        An element or ghost in the tree \a ltreeid.
 * \param [in]      face        A face number of \a elem.
 * \return                      The eclass of the local tree or ghost tree that
 *                              is face neighbor of \a elem across \a face.
 *                              T8_ECLASS_INVALID if no neighbor exists.
 */
t8_eclass_t
t8_forest_element_neighbor_eclass (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *elem,
                                   const int face);

/** Construct the face neighbor of an element, possibly across tree boundaries.
 * Returns the global tree-id of the tree in which the neighbor element lies in.
 *
 * \param [in] forest       The forest.
 * \param [in] ltreeid      The local tree in which the element lies.
 * \param [in] elem         The element to be considered.
 * \param [in,out] neigh    On input an allocated element of the scheme of the
 *                          face_neighbors eclass.
 *                          On output, this element's data is filled with the
 *                          data of the face neighbor. If the neighbor does not exist
 *                          the data could be modified arbitrarily.
 * \param [in] neigh_eclass The eclass of \a neigh.
 * \param [in] face         The number of the face along which the neighbor should be constructed.
 * \param [out] neigh_face  The number of the face viewed from perspective of \a neigh.
 * \return The global tree-id of the tree in which \a neigh is in.
 *        -1 if there exists no neighbor across that face. Domain boundary.
 *        -2 if the neighbor is not in a local tree or ghost tree. Process/Ghost boundary.
 */
t8_gloidx_t
t8_forest_element_face_neighbor (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem, t8_element_t *neigh,
                                 const t8_eclass_t neigh_eclass, int face, int *neigh_face);

/** Compute the leaf face orientation at given face in a forest.
 * \param [in]    forest  The forest. Must have a valid ghost layer.
 * \param [in]    ltreeid A local tree id.
 * \param [in]    scheme      The eclass scheme of the element.
 * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
 * \param [in]    face    The index of the face across which the face neighbors
 *                        are searched.
 * \return                Face orientation encoded as integer.
 *
 * For more information about the encoding of face orientation refer to \ref t8_cmesh_get_face_neighbor.
 */
int
t8_forest_leaf_face_orientation (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_scheme_c *scheme,
                                 const t8_element_t *leaf, const int face);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_ELEMENT_FACE_NEIGHBOR_H */
