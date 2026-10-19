/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_geometry_handler_internal.hxx
 *  Internal cmesh geoemetry handler interface functions.
 *  This file must be installed since the function are used by public
 *  t8code functions, but they should not be used by user code.
 */
#pragma once

#include <t8.h>
#include <t8_geometry/t8_geometry_handler.hxx>

/** Namespace detail to hide implementation details from the user. */
namespace detail
{

/** Return the geometry handler of the cmesh.
 * \param [in] cmesh       The cmesh to be considered. Does not need be committed.
 * \return                 The geometry handler of the cmesh.
 * \note                   The return value might be NULL if no geometry handler exists.
 */
t8_geometry_handler *
t8_cmesh_get_geometry_handler (const t8_cmesh_t cmesh);

/** Construct a new geometry_handler for a cmesh and add it to the cmesh.
 * \param [in] cmesh      The cmesh to be considered. Must be initialized. Does not need to be committed.
 * \return                On success, the new geometry_handler. nullptr on failure (out of memory).
 * \note                  Handle with care. This function should be used by t8code devs only.
 */
t8_geometry_handler *
t8_cmesh_add_geometry_handler (t8_cmesh_t cmesh);

}  // namespace detail
