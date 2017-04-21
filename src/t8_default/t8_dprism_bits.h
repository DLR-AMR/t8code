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

/** \file t8_dprism_bits.h
 */

#ifndef T8_DPRISM_BITS_H
#define T8_DPRISM_BITS_H

#include <t8_element.h>
#include "t8_dprism.h"

T8_EXTERN_C_BEGIN ();

/** Compute the level of a prism.
 * \param [in] l    Line whose prism is computed.
 * \return          The level of \a l.
 */
int                 t8_dprism_get_level (const t8_dprism_t * p);

/** Copy all values from one prism to another.
 * \param [in] l    The prism to be copied.
 * \param [in,out] dest Existing prism whose data will be filled with the data
 *                   of \a l.
 */
void                t8_dprism_copy (const t8_dprism_t * l,
                                    t8_dprism_t * dest);

/** Initialize a prism as the prism with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] l  Existing prism whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void                t8_dprism_init_linear_id (t8_dprism_t * l, int level,
                                              uint64_t id);

/** Computes the successor of a prism in a uniform grid of level \a level.
 * \param [in] l  prism whose id will be computed.
 * \param [in,out] s Existing prism whose data will be filled with the
 *                data of \a l's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void                t8_dprism_successor (const t8_dprism_t * l,
                                         t8_dprism_t * succ, int level);

/** Compute the parent of a prism.
 * \param [in]  elem Input prism.
 * \param [in,out] parent Existing prism whose data will
 *                  be filled with the data of elem's parent.
 * \note \a elem may point to the same triangle as \a parent.
 */
void                t8_dprism_parent (const t8_dprism_t * l,
                                      t8_dprism_t * parent);

/** Compute the first descendant of a prism at a given level. This is the descendant of
 * the prism in a uniform level refinement that has the smallest id.
 * \param [in] l        Prism whose descendant is computed.
 * \param [out] s       Existing prism whose data will be filled with the data
 *                      of \a l's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a l's refinement
 *                      level.
 */
void                t8_dprism_first_descendant (const t8_dprism_t * elem,
                                                t8_dprism_t * desc,
                                                int level);

/** Compute the last descendant of a prism at a given level. This is the descendant of
 * the prism in a uniform level refinement that has the largest id.
 * \param [in] l        Prism whose descendant is computed.
 * \param [out] s       Existing prism whose data will be filled with the data
 *                      of \a l's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a l's refinement
 *                      level.
 */
void                t8_dprism_last_descendant (const t8_dprism_t * l,
                                               t8_dprism_t * s, int level);

/** Compute the coordinates of a vertex of a prism.
 * \param [in] t    Input prism.
 * \param [out] coordinates An array of 2 t8_dprism_coord_t that
 * 		     will be filled with the coordinates of the vertex.
 * \param [in] vertex The number of the vertex.
 */
void                t8_dprism_vertex_coords (const t8_dprism_t * t,
                                             int vertex, int coords[]);
/** Computes the linear position of a prism in an uniform grid.
 * \param [in] line  Prism whose id will be computed.
 * \return Returns the linear position of this prism on a grid.
 */
uint64_t            t8_dprism_linear_id (const t8_dprism_t * elem, int level);

T8_EXTERN_C_END ();

#endif /* T8_DPRISM_BITS_H */
