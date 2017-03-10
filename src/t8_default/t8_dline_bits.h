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

T8_EXTERN_C_END ();

#endif /* T8_DLINE_BITS_H */
