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

/** \file t8_extruded.hxx
 * Define the extruded scheme interface.
 * Currently provides hexahedra, which are made of extruded quads (\ref t8_extruded_scheme_hex).
 * These are only refined in x- and y-direction and always span the whole tree height.
 */

#pragma once

#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh.h>

/** Return the extruded scheme implementation. */
const t8_scheme *
t8_scheme_new_extruded ();

/** Check whether a given eclass_scheme is one of the extruded schemes.
 * \param [in] scheme   A (pointer to a) scheme.
 * \param [in] eclass   The eclass to check.
 * \return              True if \a scheme uses an extruded scheme for the element class, false otherwise.
 */
int
t8_eclass_scheme_is_extruded (const t8_scheme *scheme, const t8_eclass_t eclass);

/** Check whether the face connections of a cmesh are compatible with the extruded hex scheme.
 * The extruded hex elements span the whole tree height, so the extrusion (z-) directions of neighboring hex trees
 * have to be parallel. That is, for each local hex tree:
 *  - A lateral face (0, ..., 3) is only connected to a lateral face of another hex tree, such that the z-axes of both
 *    trees are parallel (possibly with opposite directions).
 *  - The bottom and top face (4, 5) are connected to the bottom or top face of a hex tree, or to a quad face of a
 *    tree of a different class.
 * \param [in] cmesh    A committed cmesh.
 * \return              True if all face connections of the local hex trees are compatible, false otherwise.
 * \note This function is not collective. It only checks the local trees of \a cmesh.
 */
bool
t8_cmesh_is_extrusion_compatible (t8_cmesh_t cmesh);
