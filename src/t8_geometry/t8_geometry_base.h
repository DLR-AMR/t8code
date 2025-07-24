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

/** \file t8_geometry_base.h
 * This header file provides a C interface to the t8_geometry_c class.
 * Use these function if you need to query a t8_geometry_c object from
 * C code.
 * We recommend to use the C++ interface from t8_geometry_base.hxx if possible.
 */

#ifndef T8_GEOMETRY_BASE_H
#define T8_GEOMETRY_BASE_H

#include <t8.h>
#include <t8_geometry/t8_geometry.h>

T8_EXTERN_C_BEGIN ();

/** Get the name of a geometry.
 * \param [in]  geom  A geometry.
 * \return            The name of \a geom.
 */
const char *
t8_geom_get_name (const t8_geometry_c *geom);

/** Get the type of a geometry.
 *
 * \param [in] geom  A geometry.
 * \return           The type of \a geom.
 */
t8_geometry_type_t
t8_geom_get_type (const t8_geometry_c *geom);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_BASE_H */
