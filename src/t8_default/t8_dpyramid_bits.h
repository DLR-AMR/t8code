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


#ifndef T8_DPYRAMID_BITS_H
#define T8_DPYRAMID_BITS_H

#include "t8_element.h"
#include "t8_dpyramid.h"

T8_EXTERN_C_BEGIN();

/** Initialize a pyramid as the pyramid with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] p  Existing pyramid whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void                t8_dpyramid_init_linear_id (t8_dpyramid_t * p, int level,
                                              uint64_t id);

/** Compute the level of a pyramid.
 * \param [in] p    Line whose pyramid is computed.
 * \return          The level of \a p.
 */
int                 t8_dpyramid_get_level (const t8_dpyramid_t * p);


/** Computes the linear position of a pyramid in an uniform grid.
 * \param [in] p  pyramid whose id will be computed.
 * \return Returns the linear position of this pyramid on a grid.
 */
uint64_t t8_dpyramid_linear_id(const t8_dpyramid_t * p, int level);

/** Compare two elements. returns negativ if p1 < p2, zero if p1 equals p2
 *  and positiv if p1 > p2.
 *  If p2 is a copy of p1 then the elements are equal.
 */
int                 t8_dpyramid_compare (const t8_dpyramid_t * p1,
                                       const t8_dpyramid_t * p2);

/** Compute the first descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the smallest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] s       Existing pyramid whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void                t8_dpyramid_first_descendant (const t8_dpyramid_t * p,
                                                t8_dpyramid_t * desc,
                                                int level);

/** Compute the last descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the largest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] s       Existing pyramid whose data will be filled with the data
 *                      of \a p's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void
t8_dpyramid_last_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                             int level);

/** Compute the coordinates of a vertex of a pyramid.
 * \param [in] p    Input pyramid.
 * \param [out] coordinates An array of 3 t8_dpyramid_coord_t that
 * 		     will be filled with the coordinates of the vertex.
 * \param [in] vertex The number of the vertex.
 */
void                t8_dpyramid_compute_coords (const t8_dpyramid_t * p,
                                             int vertex, int coords[]);


T8_EXTERN_C_END ();

#endif /* T8_DPYRAMID_BITS_H */

