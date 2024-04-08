/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

/** \file t8_geometry_examples.h
 * This header provides the C interface to create a cubed_sphere geometry.
 */

#ifndef T8_GEOMETRY_EXAMPLES_H
#define T8_GEOMETRY_EXAMPLES_H

#include <t8.h>
#include <t8_geometry/t8_geometry.h>

T8_EXTERN_C_BEGIN ();

/** Destroy a geometry object.
 * \param [in,out] geom   A pointer to a geometry object. Set to NULL on output.
 */
void
t8_geometry_destroy (t8_geometry_c **geom);

/** Create a new quadrangulated_disk geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_quadrangulated_disk_new ();

/** Create a new triangulated_spherical_surface geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_triangulated_spherical_surface_new ();

/** Create a new quadrangulated_spherical_surface geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_quadrangulated_spherical_surface_new ();

/** Create a new cubed_spherical_shell geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_cubed_spherical_shell_new ();

/** Create a new spherical_shell geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_prismed_spherical_shell_new ();

/** Create a new cubed sphere geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_cubed_sphere_new ();

T8_EXTERN_C_END ();

#endif /* T8_GEOMETRY_EXAMPLE_H */
