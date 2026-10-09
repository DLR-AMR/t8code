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

/**
 * \file t8_cmesh_boundary_conditions_types.hxx
 * Implements helper types for boundary conditions.
 */

#pragma once

#include <t8_eclass/t8_eclass.h>
#include <t8_data/t8_static_vector.hxx>

/**
 * A container to store boundary conditions.
 * \tparam TType The type the boundary conditions are saved in.
 */
template <typename TType>
using t8_boundary_conditions = t8_static_vector<TType, T8_ECLASS_MAX_FACES>;
