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

/** \file t8_dtet_bits.h
 * TODO: Document this.
 * TODO: Run make doxygen and grep for files.
 *       Also document all arguments of functions.
 * TODO: Group the dtet_is functions together.
 */

#ifndef T8_DTET_BITS_H
#define T8_DTET_BITS_H

#include <t8_element.h>
#include "t8_dtet.h"

T8_EXTERN_C_BEGIN ();

/** Compute the coordinates of a vertex of a tetrahedron.
 * \param [in] t    Input tetrahedron.
 * \param [out] coordinates An array of 2 t8_dtet_coord_t that
 * 		     will be filled with the coordinates of the vertex.
 * \param [in] vertex The number of the vertex.
 */
void                t8_dtet_compute_coords (const t8_dtet_t * t,
                                            t8_dtet_coord_t coordinates[3],
                                            const int vertex);

/** Compute the coordinates of the four vertices of a tetrahedron.
 * \param [in] t    Input tetrahedron.
 * \param [out] coordinates An array of 4x3 t8_dtet_coord_t that
 * 		     will be filled with the coordinates of t's vertices.
 */
void                t8_dtet_compute_all_coords (const t8_dtet_t * t,
                                                t8_dtet_coord_t
                                                coordinates[4][3]);

/** Test if two tetrahedra have the same coordinates, type and level.
 * \return true if \a t1 describes the same tetrahedron as \a t2.
 */
int                 t8_dtet_is_equal (const t8_dtet_t * t1,
                                      const t8_dtet_t * t2);

/** Compute the parent of a tetrahedron.
 * \param [in]  elem Input tetrahedron.
 * \param [in,out] parent Existing tetrahedron whose data will
 *                  be filled with the data of elem's parent.
 * \note \a elem may point to the same tetrahedron as \a parent.
 */
void                t8_dtet_parent (const t8_dtet_t * t, t8_dtet_t * parent);

/** Compute the childid-th child in Bey order of a tetrahedron t.
 * \param [in] t    Input tetrahedron.
 * \param [in,out] childid The id of the child, 0..7 in Bey order.
 * \param [out] child  Existing tetrahedron whose data will be filled
 * 		    with the date of t's childid-th child.
 */
void                t8_dtet_child (const t8_dtet_t * elem,
                                   int childid, t8_dtet_t * child);

/** Compute a specific sibling of a tetrahedron.
 * \param [in]     elem  Input tetrahedron.
 * \param [in,out] sibling  Existing tetrahedron whose data will be filled
 *                    with the data of sibling no. sibling_id of elem.
 * \param [in]     sibid The id of the sibling computed, 0..7 in Bey order.
 */
void                t8_dtet_sibling (const t8_dtet_t * elem,
                                     int sibid, t8_dtet_t * sibling);

/** Compute the face neighbor of a tetrahedron.
 * \param [in]     t      Input tetrahedron.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] n      Existing tetrahedron whose data will be filled.
 * \note \a t may point to the same tetrahedron as \a n.
 */
int                 t8_dtet_face_neighbour (const t8_dtet_t * t,
                                            t8_dtet_t * n, int face);

/** Test if two tetrahedra are siblings.
 * \param [in] t1 First tetrahedron to be tested.
 * \param [in] t2 Second tetrahedron to be tested.
 * \return true if \a t1 is unequal to and a sibling of \a t2.
 */
int                 t8_dtet_is_sibling (const t8_dtet_t * t1,
                                        const t8_dtet_t * t2);

/** Test if a tetrahedron is the parent of another tetrahedron.
 * \param [in] t tetrahedron to be tested.
 * \param [in] c Possible child tetrahedron.
 * \return true if \a t is the parent of \a c.
 */
int                 t8_dtet_is_parent (const t8_dtet_t * t,
                                       const t8_dtet_t * c);

/** Test if a tetrahedron is an ancestor of another tetrahedron.
 * \param [in] t tetrahedron to be tested.
 * \param [in] c Descendent tetrahedron.
 * \return true if \a t is unequal to and an ancestor of \a c.
 */
int                 t8_dtet_is_ancestor (const t8_dtet_t * t,
                                         const t8_dtet_t * c);

T8_EXTERN_C_END ();

#endif /* T8_DTET_BITS_H */
