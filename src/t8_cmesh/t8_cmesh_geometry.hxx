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

/** \file t8_cmesh_geometry.hxx
 * Internal functions that we need for the cmesh geometry.
 * TODO: document this file
 */

#ifndef T8_CMESH_GEOMETRY_H
#define T8_CMESH_GEOMETRY_H

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_hash.hxx>

T8_EXTERN_C_BEGIN ();

/** Get the hash of the geometry stored for a tree in a cmesh.
 * \param [in] cmesh   A committed cmesh.
 * \param [in] gtreeid A global tree in \a cmesh.
 * \return             The hash of the tree's geometry. If the tree does not have a geometry, returns \ref t8_geometry_empty_hash.
 */
t8_geometry_hash_t
t8_cmesh_get_tree_geom_hash (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_GEOMETRY_H */
