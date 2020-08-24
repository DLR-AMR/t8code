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

/** \file t8_cmesh_geometry.cxx
 *
 * TODO: document this file
 */

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_geometry.h>

void
t8_cmesh_register_geometry (t8_cmesh_t cmesh, const t8_geometry_c * geometry)
{
  /* Must be called before cmesh is committed. */
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not initialized, do it now. */
    t8_geom_handler_init (&cmesh->geometry_handler);
  }

  t8_geom_handler_register_geometry (cmesh->geometry_handler, geometry);
}

void
t8_cmesh_set_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                            const char *geom_name)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  /* Add the name of the geometry as an attribute to the tree.
   * We will copy the string. */
  t8_cmesh_set_attribute_string (cmesh, gtreeid, t8_get_package_id (),
                                 T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, geom_name);
}

const t8_geometry_c *
t8_cmesh_get_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  t8_geometry_handler_t *geom_handler = cmesh->geometry_handler;
  T8_ASSERT (t8_geom_handler_is_committed (geom_handler));

  const char         *geom_name =
    t8_cmesh_get_tree_geom_name (cmesh, gtreeid);
  T8_ASSERT (geom_name != NULL);
  /* Find this geometry in the handler. */
  return t8_geom_handler_find_geometry (geom_handler, geom_name);
}

const char         *
t8_cmesh_get_tree_geom_name (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  /* Look up the name of the geometry in the attributes. */
  const char         *geom_name =
    (const char *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                           T8_CMESH_GEOMETRY_ATTRIBUTE_KEY,
                                           ltreeid);
  T8_ASSERT (geom_name != NULL);
  return geom_name;
}
