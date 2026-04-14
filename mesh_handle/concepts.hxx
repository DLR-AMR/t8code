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

/** \file concepts.hxx
 * Concepts related to the mesh handle.
 */

#pragma once

#include <type_traits>

namespace t8_mesh_handle
{
/** Concept to restrict types to \ref mesh class instantiations.
 * \tparam TType Type that should be checked.
 */
template <typename TType>
concept T8MeshType = requires { typename TType::mesh_tag; };

/** Concept to ensure that a type is MPI safe.
 * \tparam TType Type that should be checked to be MPI safe.
 */
template <typename TType>
concept T8MPISafeType
  = std::is_void_v<TType> || (std::is_trivially_copyable_v<TType> && std::is_standard_layout_v<TType>);
}  // namespace t8_mesh_handle
