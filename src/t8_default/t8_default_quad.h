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

/** \file t8_default_quad.h
 * We use a p4est_quadrant_t object as storage for the T8 quadrant.
 * To record if and if yes, how this quadrant is part of a 3D octant, we use
 * the member pad8 for the surrounding toplevel dimension (2 or 3), pad16 for
 * the direction of its normal relative to a toplevel octant (0, 1, or 2), and
 * p.user_long for the p4est_qcoord_t coordinate in the normal direction.
 */

#ifndef T8_DEFAULT_QUAD_H
#define T8_DEFAULT_QUAD_H

#include <p4est.h>
#include <t8_element.h>

/** The structure holding a quadrilateral element in the default scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */
typedef p4est_quadrant_t t8_pquad_t;

/** Return the toplevel dimension. */
#define T8_QUAD_GET_TDIM(quad) ((int) (quad)->pad8)

/** Return the direction of the third dimension.
 * This is only valid to call if the toplevel dimension is three.
 */
#define T8_QUAD_GET_TNORMAL(quad)                               \
  ( T8_ASSERT (T8_QUAD_GET_TDIM(quad) == 3),                    \
    ((int) (quad)->pad16) )

/** Return the coordinate in the third dimension.
 * This is only valid to call if the toplevel dimension is three.
 */
#define T8_QUAD_GET_TCOORD(quad)                                \
  ( T8_ASSERT (T8_QUAD_GET_TDIM(quad) == 3),                    \
    ((int) (quad)->p.user_long) )

/** Set the toplevel dimension of a quadrilateral. */
#define T8_QUAD_SET_TDIM(quad,dim)                              \
  do { T8_ASSERT ((dim) == 2 || (dim) == 3);                    \
       (quad)->pad8 = (int8_t) (dim); } while (0)

/** Set the direction of the third demension. */
#define T8_QUAD_SET_TNORMAL(quad,normal)                        \
  do { T8_ASSERT ((normal) >= 0 && (normal) < 3);               \
       (quad)->pad16 = (int16_t) (normal); } while (0)

/** Set the coordinate in the third dimension. */
#define T8_QUAD_SET_TCOORD(quad,coord)                          \
  do { (quad)->p.user_long = (long) (coord); } while (0)

/** Provide an implementation for the quadrilateral element class. */
t8_eclass_scheme_t *t8_default_scheme_new_quad (void);

#endif /* !T8_DEFAULT_QUAD_H */
