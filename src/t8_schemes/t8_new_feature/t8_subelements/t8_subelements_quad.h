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

/** \file t8_subelement_quad.h
 * We use a p4est_quadrant_t object as storage for the T8 quadrant.
 * Additionally, we store some more information to use subelements on top 
 * of some recursive quad refinement. This information is 
 * 
 *     dummy_is_subelement (is a given element a subelement?)
 *     subelement_type (what type of transition cell is used?)
 *     subelement_id (what subelement of the transition cell is the given subelement?)
 * 
 * In order to refine a quad element using subelements, it is important to know these additional information. 
 * We can use them for example to determine the coordinates of the subelement vertices, 
 * since they can differ for subelements of the same id but a different subelement type. */

#ifndef T8_DEFAULT_QUAD_H
#define T8_DEFAULT_QUAD_H

#include <p4est.h>
#include <t8_element.h>

/** The structure holding a quadrilateral element in the default scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */

/* Define the struct, that stores all information needed for the quad scheme and subelements.
 * 
 *         p4est quadrant        recursive quad       refinement, using a transition 
 *                                 refinement             cell with subelements
 *         x - - - - - x         x - - x - - x                x - - - - - x
 *         |           |         |     |     |                | \   2   / |       
 *         |           |         |     |     |                | 1 \   /   |
 *         |           |   -->   x - - x - - x       or       x - - x   3 |
 *         |           |         |     |     |                | 0 / | \   |
 *         |           |         |     |     |                | / 5 | 4 \ |
 *         x - - - - - x         x - - x - - x                x - - x - - x 
 * 
 * A p4est quadrant can be refined, using either the standard quad scheme, or a transition cell, consisting of different subelements. 
 * The quad refinement scheme is recursive, whereas a transition cell can only be used once, for example to remove hanging nodes, after the mesh has been adapted and balanced. 
 * There are different types of transition cells possible, which we will refer to as subelement_type. 
 * Each transition cell consists of different subelements. The given example consists of 6 different subelements, whose ids range from 0 to 5.
 * A dummy variable will store the information, whether a given element is a subelement or a standard quad element. */

typedef struct
{
  p4est_quadrant_t    p4q;      /* p4est quadrant */
  int                 dummy_is_subelement;      /* saves the information, whether an element is a subelement */
  int                 subelement_type;  /* saves the information, which type of transition cell a subelement is associated to */
  int                 subelement_id;    /* saves the information, what children subelement the given element is */
} t8_quad_with_subelements;

typedef t8_quad_with_subelements t8_pquad_t;

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
