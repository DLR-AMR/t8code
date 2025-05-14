/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_cmesh.hxx
 * We define the coarse mesh of trees in this file.
 */

#ifndef T8_CMESH_HXX
#define T8_CMESH_HXX

#include <vector>
#include <unordered_map>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity_types.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

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
  if (cmesh->geometry_handler == NULL) {
    /* The handler was not constructed, do it now. */
    cmesh->geometry_handler = new t8_geometry_handler ();
  }
  return cmesh->geometry_handler->register_geometry<geometry_type> (std::forward<_args> (args)...);
}

/** Get the list of global trees and local vertex ids a global vertex is connected to.
 * Cmesh Interface function.
 *  
 * \param [in] cmesh A committed cmesh.
 * \param [in] global_vertex_id The global id of a vertex in \a cmesh.
 * \return The list of global tree ids and local vertex ids of \a global_vertex_id.
 */
const tree_vertex_list &
t8_cmesh_get_vertex_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_vertex);

#endif /* T8_CMESH_HXX */
