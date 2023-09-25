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

/* In this file we collect the implementations of the geometry base class
 * and its C interface.
 * */

#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_base.h>

/** Get the dimension of a geometry.
 * \param [in]  geom  A geometry.
 * \return            The dimension of \a geom.
 */
int
t8_geom_get_dimension (const t8_geometry_c *geom)
{
  T8_ASSERT (geom != NULL);

  return geom->t8_geom_get_dimension ();
}

/** Get the name of a geometry.
 * \param [in]  geom  A geometry.
 * \return            The name of \a geom.
 */
const char *
t8_geom_get_name (const t8_geometry_c *geom)
{
  T8_ASSERT (geom != NULL);

  return geom->t8_geom_get_name ();
}
