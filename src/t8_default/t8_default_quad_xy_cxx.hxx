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

#ifndef T8_DEFAULT_QUAD_YX_CXX_HXX
#define T8_DEFAULT_QUAD_YX_CXX_HXX

#include <p4est.h>
#include <t8_element_cxx.hxx>
#include "t8_default_common_cxx.hxx"
#include "t8_default_quad_cxx.hxx"

struct t8_default_scheme_quad_xy_c:public t8_default_scheme_quad_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_default_scheme_quad_xy_c ();

  ~t8_default_scheme_quad_xy_c ();

  /** Return the tree face id given a boundary face. */
  virtual int         t8_element_tree_face (const t8_element_t * elem,
                                            int face);

  /** Compute the integer coordinates of a given element vertex. */
  virtual void        t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[]);
};

#endif /* !T8_DEFAULT_QUAD_YX_CXX_HXX */
