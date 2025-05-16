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
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

void
t8_cmesh_register_geometry (t8_cmesh_t cmesh, t8_geometry_c *geometry)
{
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not constructed, do it now. */
    cmesh->geometry_handler = new t8_geometry_handler ();
  }
  cmesh->geometry_handler->register_geometry (geometry);
}

void
t8_cmesh_set_tree_geometry (t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, const t8_geometry_c *geom)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  /* Add the hash of the geometry as an attribute to the tree. */
  t8_geometry_hash_t hash = geom->t8_geom_get_hash ();
  t8_cmesh_set_attribute (cmesh, gtreeid, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, &hash,
                          sizeof (t8_geometry_hash_t), 0);
}

const t8_geometry_c *
t8_cmesh_get_tree_geometry (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  t8_geometry_handler *geom_handler = cmesh->geometry_handler;

  if (geom_handler == nullptr) {
    /* If no geometry handler is present, no geometries have been registered and 
    * the tree does not have a geometry. */
    return nullptr;
  }

  if (geom_handler->get_num_geometries () == 1) {
    /* The geometry handler only has one geometry and the trees 
     * thus do not need to store their geometry's hash
     * (we assume all trees have this geometry).
     */
    return geom_handler->get_unique_geometry ();
  }
  const t8_geometry_hash_t geom_hash = t8_cmesh_get_tree_geom_hash (cmesh, gtreeid);

  /* Look up the geometry in the geometry handler's hash table and return it. */
  return geom_handler->get_geometry (geom_hash);
}

t8_geometry_hash_t
t8_cmesh_get_tree_geom_hash (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  t8_geometry_handler *geom_handler = cmesh->geometry_handler;

  if (geom_handler == nullptr) {
    /* If no geometry handler is present, no geometries have been registered and 
    * the tree does not have a geometry. */
    return t8_geometry_empty_hash;
  }

  if (geom_handler->get_num_geometries () == 1) {
    /* There is only one geometry registered in this cmesh, so we assume
     * that this geometry is used for all trees. */
    auto geom = geom_handler->get_unique_geometry ();
    T8_ASSERT (geom != NULL);  // If there is exactly one geometry, then it must be active.
#if T8_ENABLE_DEBUG
    /* In debug mode, get the tree's geometry anyways and check that it is either
     * NULL or the hash of the unique geometry. */

    t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
    /* Look up the hash of the geometry in the attributes. */
    const t8_geometry_hash_t *geom_hash = (const t8_geometry_hash_t *) t8_cmesh_get_attribute (
      cmesh, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, ltreeid);
    T8_ASSERT (geom_hash != NULL);
    T8_ASSERT (*geom_hash == geom->t8_geom_get_hash ());
#endif /* T8_ENABLE_DEBUG */
    return geom->t8_geom_get_hash ();
  }

  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  /* Look up the hash of the geometry in the attributes. */
  const t8_geometry_hash_t *geometry_hash = (const t8_geometry_hash_t *) t8_cmesh_get_attribute (
    cmesh, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, ltreeid);
  return *geometry_hash;
}
