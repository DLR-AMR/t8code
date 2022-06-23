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

/** \file t8_dhex_bits.h
 */

#ifndef T8_DHEX_BITS_H
#define T8_DHEX_BITS_H

#include <t8.h>
#include <t8_schemes/t8_default/t8_default_hex/t8_dhex.h>

T8_EXTERN_C_BEGIN ();

/** Compute the parent of given hexaeder.
 * \param [in]  q Input hexaeder.
 * \param [in,out] r Existing hexaeder whose Morton index will be filled
 *                   with the Morton index of the parent of \a q.
 */
void t8_dhex_parent (const t8_dhex_t * q, t8_dhex_t * r);

/** Compute the specified sibling of given hexaeder.
 * \param [in]     q  Input hexaeder.
 * \param [in,out] r  Existing hexaeder whose Morton index will be filled
 *                    with the coordinates of specified sibling.
 * \param [in]     sibling_id The id of the sibling computed: 0,...,7.
 */
void t8_dhex_sibling (const t8_dhex_t * q, t8_dhex_t * r,
                                      int sibling_id);

/** Compare two hexaeders with respect to their Morton ordering.
 * \return < 0 if \a v1 < \a v2,
 *           0 if \a v1 == \a v2,
 *         > 0 if \a v1 > \a v2
 */
int t8_dhex_compare (const void *v1, const void *v2);

/** Test if a hexaeder has valid Morton indices
 * in the 3x3x3 box around root.
 * \param [in] q Quadrant to be tested.
 * \return True if \a q is extended.
 */
int t8_dhex_is_extended (const t8_dhex_t * q);

/** Test if a hexaeder is parent of another hexaeder.
 * \param [in] q Quadrant to be tested.
 * \param [in] r Possible child hexaeder.
 * \return true if \a q is the parent of \a r.
 */
int t8_dhex_is_parent (const t8_dhex_t * q, const t8_dhex_t * r);

/** Compute the eight children of a hexaeder, array version.
 * \param [in]     q  Input hexaeder.
 * \param [in,out] c  Pointers to the eight computed children in z-order.
 *                    q may point to the same hexaeder as c[0].
 * \note The user_data of c[i] is never modified.
 */
void t8_dhex_childrenpv (const t8_dhex_t * q, t8_dhex_t * c[]);

/** Compute the position of this child within its siblings.
 * \return Child id: 0,...,7
 */
int t8_dhex_child_id (const t8_dhex_t * q);

/** Compute the position of the ancestor of this child at
 * level \a level within its siblings.
 * \return Child id: 0,...,7
 */
int t8_dhex_ancestor_id (const t8_dhex_t * q, int level);

/** Test if eight hexaeders are siblings in Morton ordering, array version.
 * \param [in] q   Array of eight pointers to hexaeders.
 */
int t8_dhex_is_familypv (t8_dhex_t * q[]);

/** Set hexaeder Morton indices based on linear position in uniform grid.
 * \param [in,out] q         Quadrant whose Morton indices will be set.
 * \param [in]     level     Level of the grid and of the resulting hexaeder.
 * \param [in]     id        Linear index of the hexaeder on a uniform grid.
 */
void t8_dhex_set_morton (t8_dhex_t * q, int level, uint64_t id);

/** Computes the linear position of a hexaeder in a uniform grid.
 * \param [in] q         Quadrant whose linear index will be computed.
 * \param [in] level     The level of the regular grid compared to which the
 *                       linear position is to be computed.
 * \return               Linear position of this hexaeder on a grid.
 * \note                 If the hexaeder is smaller than the grid (has a higher
 *                       hexaeder->level), the result is computed from its
 *                       ancestor at the grid's level.
 *                       If the hexaeder has a smaller level than the grid (it
 *                       is bigger than a grid cell), the grid cell sharing its
 *                       lower left corner is used as reference.
 */
uint64_t t8_dhex_linear_id (const t8_dhex_t * q, int level);

/** Compute the first descendant of a given hexaeder on a given level.
 * \param [in]  q      Input hexaeder.
 * \param [out] fd     First descendant of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void t8_dhex_first_descendant (const t8_dhex_t * q, t8_dhex_t * fd, int level);

/** Computes the nearest common ancestor of two hexaeders in the same tree.
 * \param [in]     q1 First input hexaeder.
 * \param [in]     q2 Second input hexaeder.
 * \param [in,out] r Existing hexaeder whose Morton index will be filled.
 * \note \a q1, \a q2, \a r may point to the same hexaeder.
 */
void t8_dhex_nearest_common_ancestor (const t8_dhex_t * q1, const t8_dhex_t * q2, t8_dhex_t * r);

/** Compute the descendant of a hexaeder touching a given corner.
 * \param [in]     q   Input hexaeder.
 * \param [in,out] r   Existing hexaeder whose Morton index will be filled.
 * \param [in]     c   The corner of \a q that \a r touches.
 * \param [in] level   The size of \a r.
 */
void t8_dhex_corner_descendant (const t8_dhex_t * q, t8_dhex_t * r, int c, int level);

/** Compute the face neighbor of a hexaeder.
 * \param [in]     q      Input hexaeder, must be valid.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] r      Existing hexaeder whose Morton index will be filled.
 * \note \a q may point to the same hexaeder as \a r.
 */
void t8_dhex_face_neighbor (const t8_dhex_t * q, int face, t8_dhex_t * r);

/** Test if a hexaeder is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return True if \a q is inside the unit tree.
 */
int t8_dhex_is_inside_root (const t8_dhex_t * q);

/** Prints one line with x, y and level of the hexaeder.
 * \param [in] log_priority  See \ref logpriorities in sc.h.
 * \param [in] q             Puadrant to print.
 */
void t8_dhex_print (int log_priority, const t8_dhex_t * q);

/** Compute the last descendant of a hexaeder on a given level.
 * \param [in]  q      Input hexaeder.
 * \param [out] ld     Last descendant of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void t8_dhex_last_descendant (const t8_dhex_t * q, t8_dhex_t * ld, int level);

/** Test if a hexaeder is used to represent a mesh node.
 * \param [in] q        Quadrant to be tested.
 * \param [in] inside   If true, boundary nodes must be clamped inside.
 *                      If false, nodes must align with the hexaeder grid.
 * \return True if \a q is a node.
 */
int t8_dhex_is_node (const t8_dhex_t * q, int inside);

/** Test if a hexaeder has valid Morton indices and is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return True if \a q is valid.
 */
int t8_dhex_is_valid (const t8_dhex_t * q);

/** Test if eight hexaeders are siblings in Morton ordering.
 */
int t8_dhex_is_family (const t8_dhex_t * q0,
                       const t8_dhex_t * q1,
                       const t8_dhex_t * q2,
                       const t8_dhex_t * q3,
                       const t8_dhex_t * q4,
                       const t8_dhex_t * q5,
                       const t8_dhex_t * q6,
                       const t8_dhex_t * q7);

/** Compute the eight children of a hexaeder.
 * \param [in]     q  Input hexaeder.
 * \param [in,out] c0 First computed child.
 * \note  \a c0 may point to the same hexaeder as \a q.
 */
void t8_dhex_children (const t8_dhex_t * q,
                    t8_dhex_t * c0, t8_dhex_t * c1,
                    t8_dhex_t * c2, t8_dhex_t * c3,
                    t8_dhex_t * c4, t8_dhex_t * c5,
                    t8_dhex_t * c6, t8_dhex_t * c7);

T8_EXTERN_C_END ();

#endif
