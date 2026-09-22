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
 */

#ifndef T8_CMESH_GEOMETRY_H
#define T8_CMESH_GEOMETRY_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_hash.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

T8_EXTERN_C_BEGIN ();

/** Get the hash of the geometry stored for a tree in a cmesh.
 * \param [in] cmesh   A committed cmesh.
 * \param [in] gtreeid A global tree in \a cmesh.
 * \return             The hash of the tree's geometry. If the tree does not have a geometry, returns t8_geometry_empty_hash.
 */
t8_geometry_hash
t8_cmesh_get_tree_geom_hash (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

/** Return the geometry handler of the cmesh.
 * \param [in] cmesh       The cmesh to be considered. Does not need be committed.
 * \return                 The geometry handler of the cmesh.
 * \note                   The return value might be NULL if no geometry handler exists.
 */
detail::t8_geometry_handler *
t8_cmesh_get_geometry_handler (const t8_cmesh_t cmesh);

/** Construct a new geometry_handler for a cmesh and add it to the cmesh.
 * \param [in] cmesh      The cmesh to be considered. Does not need to be committed.
 * \return                On success, the new geometry_handler. nullptr on failure (out of memory).
 * \note                  Handle with care. This function should be used by t8code devs only.
 */
detail::t8_geometry_handler *
t8_cmesh_add_geometry_handler (t8_cmesh_t cmesh);

T8_EXTERN_C_END ();

/**
 * Create and register a geometry with the coarse mesh. The coarse mesh takes the ownership of the geometry.
 * @tparam geometry_type 
 * \param [in,out] cmesh The cmesh.
 * \param [in,out] args The constructor arguments of the geometry.
 * \return         A pointer to the geometry.
 */
template <typename geometry_type, typename... _args>
inline geometry_type *
t8_cmesh_register_geometry (t8_cmesh_t cmesh, _args &&...args)
{
  detail::t8_geometry_handler *geometry_handler = t8_cmesh_get_geometry_handler (cmesh);
  if (geometry_handler == nullptr) {
    /* The handler was not constructed, do it now. */
    geometry_handler = t8_cmesh_add_geometry_handler (cmesh);
    T8_ASSERT (geometry_handler != nullptr);
  }
  return geometry_handler->register_geometry<geometry_type> (std::forward<_args> (args)...);
}

#endif /* !T8_CMESH_GEOMETRY_H */
