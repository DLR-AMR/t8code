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

#include <p4est_bits.h>
#include "t8_dline_bits.h"
#include "t8_default_common_cxx.hxx"
#include "t8_default_quad_yx_cxx.hxx"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

int
t8_default_scheme_quad_xy_c::t8_element_tree_face (const t8_element_t * elem,
                                                   int face)
{
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  /* The quadrants are rotated by 90 degrees:
   *
   * Face numbers of the quadrants
   *    1
   *   x--x
   * 2 |  | 3
   *   x--x
   *    0
   *
   * Thus quadrant_face -> tree_face
   *          0 -> 2
   *          1 -> 3    f xor 2
   *          2 -> 0
   *          3 -> 1
   */
  return face ^ 2;
}

void
t8_default_scheme_quad_xy_c::t8_element_vertex_coords (const t8_element_t * t,
                                                       int vertex,
                                                       int coords[])
{
  const p4est_quadrant_t *q1 = (const p4est_quadrant_t *) t;
  int                 len;

  T8_ASSERT (0 <= vertex && vertex < 4);
  /* Get the length of the quadrant */
  len = P4EST_QUADRANT_LEN (q1->level);
  /* Compute the x and y coordinates of the vertex depending on the
   * vertex number */
  coords[0] = q1->y + (vertex & 1 ? 1 : 0) * len;
  coords[1] = q1->x + (vertex & 2 ? 1 : 0) * len;
}

/* Constructor */
t8_default_scheme_quad_xy_c::t8_default_scheme_quad_xy_c (void)
{
  /* empty since we use the construction of t8_default_scheme_quad_c */
}

t8_default_scheme_quad_xy_c::~t8_default_scheme_quad_xy_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
