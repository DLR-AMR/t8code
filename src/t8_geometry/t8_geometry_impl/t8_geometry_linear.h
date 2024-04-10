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

/** \file t8_geometry_linear.h
 * This header provides the C interface to create a linear geometry.
 */

#ifndef T8_GEOMETRY_LINEAR_H
#define T8_GEOMETRY_LINEAR_H

#include <t8.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_with_vertices.h>

T8_EXTERN_C_BEGIN ();

/**
 * Create a new linear geometry of a given dimension.
 * The geometry is only all tree types and as many vertices as the tree type
 * has. The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
 * Sets the dimension and the name to "t8_geom_linear_{dim}"
 * \param [in] dim  0 <= \a dimension <= 3. The dimension.
 * \return          A pointer to an allocated t8_geometry_linear struct, as
 *                  if the \ref t8_geometry_linear (int dim) constructor was called.
 */
t8_geometry_c *
t8_geometry_linear_new (int dim);

/** Destroy a linear geometry that was created with \ref t8_geometry_linear_new.
 * \param [in,out] geom A linear geometry. Set to NULL on output.
 */
void
t8_geometry_linear_destroy (t8_geometry_c **geom);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_LINEAR_H */
