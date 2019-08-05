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

/** \file t8_no_neighbors_quad_cxx.hxx
 *
 * This quadrant uses the same implementation as \ref t8_default_scheme_quad_c
 * except that this quadrant does not have any face neighbors.
 */

#ifndef T8_NO_NEIGHBORS_QUAD_CXX_HXX
#define T8_NO_NEIGHBORS_QUAD_CXX_HXX

#include <p4est.h>
#include <t8_element_cxx.hxx>
#include "../t8_default/t8_default_quad_cxx.hxx"
#include "../t8_default/t8_default_line_cxx.hxx"


struct t8_no_neighbors_scheme_quad_c:public t8_default_scheme_quad_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_no_neighbors_scheme_quad_c ();

  ~t8_no_neighbors_scheme_quad_c ();

   /** Transform the coordinates of a quadrilateral considered as boundary element
   *  in a tree-tree connection. */
  virtual void        t8_element_transform_face (const t8_element_t * elem1,
                                                 t8_element_t * elem2,
                                                 int orientation, int sign,
                                                 int is_smaller_face);

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise. */
  virtual int         t8_element_face_neighbor_inside (const t8_element_t *
                                                       elem,
                                                       t8_element_t * neigh,
                                                       int face,
                                                       int *neigh_face);
};

#endif /* !T8_NO_NEIGHBORS_QUAD_CXX_HXX */
