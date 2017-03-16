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

/** \file t8_dline_bits.h
 */

#ifndef T8_DLINE_BITS_H
#define T8_DLINE_BITS_H

#include <t8_element.h>
#include "t8_dline.h"

T8_EXTERN_C_BEGIN ();

/** Compute the level of a line.
 * \param [in] l    Line whose level is computed.
 * \return          The level of \a l.
 */
int                 t8_dline_get_level (const t8_dline_t * l);

/** Copy all values from one line to another.
 * \param [in] l    The line to be copied.
 * \param [in,out] dest Existing line whose data will be filled with the data
 *                   of \a l.
 */
void                t8_dline_copy (const t8_dline_t * l, t8_dline_t * dest);

/** Compute the parent of a line.
 * \param [in]  l   The input line.
 * \param [in,out] parent Existing line whose data will be filled with the parent
 *                  data of \a l.
 */
void                t8_dline_parent (const t8_dline_t * l,
                                     t8_dline_t * parent);

/** Compute the childid-th child in Morton order of a line.
 * \param [in] l    Input Line.
 * \param [in,out] childid The id of the child, 0 or 1, in Morton order.
 * \param [out] child  Existing Line whose data will be filled
 * 		    with the date of l's childid-th child.
 */
void                t8_dline_child (const t8_dline_t * l, int childid,
                                    t8_dline_t * child);

/** Initialize a line as the line with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] l  Existing line whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void                t8_dline_init_linear_id (t8_dline_t * l, int level,
                                             uint64_t id);

/** Computes the successor of a line in a uniform grid of level \a level.
 * \param [in] l  line whose id will be computed.
 * \param [in,out] s Existing line whose data will be filled with the
 *                data of \a l's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void                t8_dline_successor (const t8_dline_t * l,
                                        t8_dline_t * succ, int level);

/** Suppose we have two trees that share a common face f.
 *  Given a Line e that is a subface of f in one of the trees
 *  and given the orientation of the tree connection, construct the face
 *  Line of the respective tree neighbor that logically coincides with e
 *  but lies in the coordinate system of the neighbor tree.
 *  \param [in] elem1     The face element.
 *  \param [in,out] elem2 On return the face elment \a elem1 with respective
 *                        to the coordinate system of the other tree.
 *  \param [in] orientation The orientation of the tree-tree connection.
 *                        0 if vertex 0 of face 0 coincides with vertex 0 of face 1.
 *                        1 if vertex 0 of face 0 coincides with vertex 1 of face 1.
 */
void                t8_dline_transform_face (const t8_dline_t * line1,
                                             t8_dline_t * line2,
                                             int orientation);

/** Compute the first descendant of a line at a given level. This is the descendant of
 * the line in a uniform level refinement that has the smallest id.
 * \param [in] l        Line whose descendant is computed.
 * \param [out] s       Existing line whose data will be filled with the data
 *                      of \a l's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a l's refinement
 *                      level.
 */
void                t8_dline_first_descendant (const t8_dline_t * l,
                                               t8_dline_t * s, int level);

T8_EXTERN_C_END ();

#endif /* T8_DLINE_BITS_H */
