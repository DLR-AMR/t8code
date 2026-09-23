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

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_geometry_internal.hxx>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

void
t8_cmesh_register_geometry (t8_cmesh_t cmesh, t8_geometry_c *geometry)
{
  detail::t8_geometry_handler *geometry_handler = t8_cmesh_get_geometry_handler (cmesh);
  if (geometry_handler == nullptr) {
    /* The handler was not constructed, do it now. */
    t8_cmesh_add_geometry_handler (cmesh);
  }
  geometry_handler->register_geometry (geometry);
}

void
t8_cmesh_set_tree_geometry (t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, const t8_geometry_c *geom)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  /* Add the hash of the geometry as an attribute to the tree. */
  t8_geometry_hash hash = geom->t8_geom_get_hash ();
  t8_cmesh_set_attribute (cmesh, gtreeid, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, &hash,
                          sizeof (t8_geometry_hash), 0);
}

/* Return the geometry handler of the cmesh. C version.
 * \param [in] cmesh       The cmesh to be considered. Does not need be committed.
 * \return                 The geometry handler of the cmesh.
 * \note                   The return value might be NULL if no geometry handler exists.
 * \note                   Handle with care. This function should be used by t8code devs only.
 */
t8_geometry_handler_c *
t8_cmesh_get_geometry_handler_c (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh) || t8_cmesh_is_committed (cmesh));

  return cmesh->geometry_handler;
}

/* Return the geometry handler of the cmesh.
 * \param [in] cmesh       The cmesh to be considered. Does not need be committed.
 * \return                 The geometry handler of the cmesh.
 * \note                   The return value might be NULL if no geometry handler exists.
 */
detail::t8_geometry_handler *
t8_cmesh_get_geometry_handler (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh) || t8_cmesh_is_committed (cmesh));

  return detail::t8_geom_handler_from_c (t8_cmesh_get_geometry_handler_c (cmesh));
}

/* Construct a new geometry_handler for a cmesh and add it to the cmesh.
 * \param [in] cmesh      The cmesh to be considered. Must be initialized. Does not need to be committed.
 * \return                On success, the new geometry_handler. nullptr on failure (out of memory).
 * \note                  Handle with care. This function should be used by t8code devs only.
 */
detail::t8_geometry_handler *
t8_cmesh_add_geometry_handler (t8_cmesh_t cmesh)
{
  return t8_cmesh_set_geometry_handler (cmesh, nullptr);
}

/* Set a geometry handler or construct a new geometry_handler for a cmesh and add it to the cmesh.
 * \param [in] cmesh      The cmesh to be considered. Must be initialized. Does not need to be committed.
 * \param [in] new_handler  The geometry handler to be set. If nullptr then a new handler will be allocated.
 * \return                On success, the new geometry_handler. nullptr on failure (out of memory).
 */
detail::t8_geometry_handler *
t8_cmesh_set_geometry_handler (t8_cmesh_t cmesh, detail::t8_geometry_handler *new_handler)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh) || t8_cmesh_is_committed (cmesh));

  // Check that we do not overwrite an existing handler.
  T8_ASSERT (t8_cmesh_get_geometry_handler (cmesh) == nullptr);

  // If not present allocate a new handler
  if (new_handler == nullptr) {
    new_handler = new detail::t8_geometry_handler ();
  }
  else {
    // Increase reference count of handler
    new_handler->ref ();
  }
  // Convert the handler to C pointer and add to cmesh
  t8_geometry_handler_c *new_handler_c = detail::t8_geom_handler_to_c (new_handler);
  return cmesh->geometry_handler = new_handler_c;
}

const t8_geometry_c *
t8_cmesh_get_tree_geometry (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  detail::t8_geometry_handler *geom_handler = t8_cmesh_get_geometry_handler (cmesh);

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
  const t8_geometry_hash geom_hash = t8_cmesh_get_tree_geom_hash (cmesh, gtreeid);

  /* Look up the geometry in the geometry handler's hash table and return it. */
  return geom_handler->get_geometry (geom_hash);
}

t8_geometry_hash
t8_cmesh_get_tree_geom_hash (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  detail::t8_geometry_handler *geom_handler = t8_cmesh_get_geometry_handler (cmesh);

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
    const t8_geometry_hash *geom_hash = (const t8_geometry_hash *) t8_cmesh_get_attribute (
      cmesh, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, ltreeid);
    T8_ASSERT (geom_hash != NULL);
    T8_ASSERT (*geom_hash == geom->t8_geom_get_hash ());
#endif /* T8_ENABLE_DEBUG */
    return geom->t8_geom_get_hash ();
  }

  t8_locidx_t const ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  /* Look up the hash of the geometry in the attributes. */
  const t8_geometry_hash *geometry_hash = (const t8_geometry_hash *) t8_cmesh_get_attribute (
    cmesh, t8_get_package_id (), T8_CMESH_GEOMETRY_ATTRIBUTE_KEY, ltreeid);
  return *geometry_hash;
}
