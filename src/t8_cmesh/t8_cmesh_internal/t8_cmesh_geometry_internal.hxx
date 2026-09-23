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

/** \file t8_cmesh_geometry_internal.hxx
 *  Internal cmesh geoemetry interface functions
 */
#pragma once

#include <t8.h>
#include <t8_geometry/t8_geometry_handler.hxx>

/** Set a geometry handler or construct a new geometry_handler for a cmesh and add it to the cmesh.
 * \param [in] cmesh      The cmesh to be considered. Must be initialized. Does not need to be committed.
 * \param [in] new_handler  The geometry handler to be set. If nullptr then a new handler will be allocated.
 * \return                On success, the new geometry_handler. nullptr on failure (out of memory).
 */
detail::t8_geometry_handler *
t8_cmesh_set_geometry_handler (t8_cmesh_t cmesh, detail::t8_geometry_handler *new_handler);
