/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file t8_default_tet_bits.h
 */

#ifndef T8_DEFAULT_TET_BITS_H
#define T8_DEFAULT_TET_BITS_H

#include <t8_element.h>
#include "t8_default_tet.h"

/** Test if two tetrahedra have the same coordinates, type and level.
 * \return true if \a q1 describes the same quadrant as \a q2.
 */
int                 p4est_quadrant_is_equal (const t8_tet_t *t1,
                                             const t8_tet_t *t2);

/** Compute the parent of a tetrahedron.
 * \param [in]  elem Input tetrahedron.
 * \param [in,out] parent Existing tetrahedron whose data will
 *                  be filled with the data of elem's parent.
 * \note \a elem may point to the same quadrant as \a parent.
 */
void                t8_default_tet_parent (const t8_element_t * elem,
                                           t8_element_t * parent);

/** Compute the coordinates of the four vertices of a tetrahedron.
 * \param [in] t    Input tetrahedron.
 * \param [out] coordinates An array of 4x3 t8_tcoord_t that
 * 		     will be filled with the coordinates of t's vertices.
 */
void                t8_default_tet_compute_coords (const t8_tet_t * t,
                                                   t8_tcoord_t
                                                   coordinates[4][3]);

/** Compute the childid-th child in Bey order of a tetrahedron t.
 * \param [in] t    Input tetrahedron.
 * \param [in,out] childid The id of the child, 0..7 in Bey order.
 * \param [out] child  Existing tetrahedron whose data will be filled
 * 		    with the date of t's childid-th child.
 */
void                t8_default_tet_child (const t8_element_t * elem,
                                          int childid, t8_element_t * child);

/** Compute a specific sibling of a tetrahedron.
 * \param [in]     elem  Input tetrahedron.
 * \param [in,out] sibling  Existing tetrahedron whose data will be filled
 *                    with the data of sibling no. sibling_id of elem.
 * \param [in]     sibid The id of the sibling computed, 0..7 in Bey order.
 */
void                t8_default_tet_sibling (const t8_element_t * elem,
                                            int sibid,
                                            t8_element_t * sibling);

/** Compute the face neighbor of a tetrahedron.
 * \param [in]     t      Input tetrahedron.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] n      Existing tetrahedron whose data will be filled.
 * \note \a t may point to the same quadrant as \a n.
 */
int                 t8_default_tet_face_neighbour (const t8_tet_t * t,
                                                   t8_tet_t * n, int face);

#endif /* T8_DEFAULT_TET_BITS_H */
