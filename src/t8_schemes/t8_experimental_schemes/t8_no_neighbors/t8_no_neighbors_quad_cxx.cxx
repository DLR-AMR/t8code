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
#include "t8_no_neighbors_quad_cxx.hxx"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

void
t8_no_neighbors_scheme_quad_c::t8_element_transform_face (const t8_element_t *
                                                     elem1,
                                                     t8_element_t * elem2,
                                                     int orientation,
                                                     int sign,
                                                     int is_smaller_face)
{
    /* There are no face neighbors. Calling this function does not make sense. */
}

int
t8_no_neighbors_scheme_quad_c::t8_element_face_neighbor_inside (const t8_element_t
                                                           * elem,
                                                           t8_element_t *
                                                           neigh, int face,
                                                           int *neigh_face)
{
  return 0;
}

/* Constructor */
t8_no_neighbors_scheme_quad_c::t8_no_neighbors_scheme_quad_c (void)
{
    /* Base class constructor suffices. Hence, this constructor is empty */
}

t8_no_neighbors_scheme_quad_c::~t8_no_neighbors_scheme_quad_c ()
{
  /* This destructor is empty since the destructor of the
   * no_neighbors_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
