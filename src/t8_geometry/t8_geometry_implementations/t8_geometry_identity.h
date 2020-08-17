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

/** \file t8_geometry_identity.h
 * This header provides the C interface to create an identity geometry.
 */

#ifndef T8_GEOMETRY_IDENTITY_H
#define T8_GEOMETRY_IDENTITY_H

#include <t8.h>
#include <t8_geometry/t8_geometry.h>

/* Create a new identity geometry of a given dimension.
 * \param [in] dimension  0 <= \a dimension <= 3. The dimension.
 * \return          A pointer to an allocated t8_geometry_identity struct, as
 *                  if the t8_geometry_identity (int dimension) constructor was called.
 */
t8_geometry_c      *t8_geometry_identitiy_new (int dimension);

#endif /* !T8_GEOMETRY_IDENTITY_H! */
