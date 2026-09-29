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
 * \file t8_forest_boundary_conditions.hxx
 * Public interface for retrieving boundary conditions of forest elements.
 */

#pragma once

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_element/t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions_types.hxx>
#include <optional>
#include <string_view>

/**
 * Retrieves the boundary conditions of a forest element.
 *
 * \param [in] forest   The forest the element lives in.
 * \param [in] ltreeid  The local id of the forest tree.
 * \param [in] element  The element.
 * \return A container with the boundary conditions. Note, that only element faces at the boundary of a
 * tree will have boundary conditions. Internal faces will return an empty optional.
 */
t8_boundary_conditions<std::optional<std::string_view>>
t8_forest_get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/**
 * Retrieves the boundary condition of a face of a forest element.
 * Retrieving all boundary conditions at once via \ref t8_forest_get_boundary_conditions() will be faster.
 *
 * \param [in] forest   The forest the element lives in.
 * \param [in] ltreeid  The local id of the forest tree.
 * \param [in] element  The element.
 * \param [in] face     The face id of the element.
 * \return The boundary condition. It will be empty if the element is not touching the boundary of the tree.
 */
std::optional<std::string_view>
t8_forest_get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face);
