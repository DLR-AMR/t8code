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
 * \file t8_cmesh_boundary_conditions.h
 * Public interface for the definition and retrieval of boundary conditions.
 */

#pragma once

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler_types.h>
#include <t8_data/t8_static_vector.hxx>

#include <vector>
#include <string_view>

/**
 * Applies boundary conditions to the faces of a cmesh cell.
 *
 * \tparam TStringRange             An iterable container filled with string like values.
 * \param [in] cmesh                The cmesh.
 * \param [in] gtreeid              The global id of the tree the boundary conditions should be set for.
 * \param [in] boundary_conditions  The boundary conditions to set. Container must have the same length
 *                                  as the eclass of the cell has faces.
 */
template <std::ranges::input_range TStringRange>
  requires std::convertible_to<std::ranges::range_reference_t<TStringRange>, std::string_view>
void
t8_cmesh_set_boundary_conditions (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, TStringRange boundary_conditions)
{
  detail::t8_cmesh_boundary_condition_handler *handler = t8_cmesh_get_boundary_condition_handler (cmesh);
  if (handler == nullptr) {
    handler = t8_cmesh_add_boundary_condition_handler (cmesh);
  }
  handler->add_boundary_conditions (gtreeid, boundary_conditions);
}

/**
 * Retrieves the boundary conditions of a cmesh cell.
 *
 * \param [in] cmesh    The cmesh the cell lives in.
 * \param [in] ltreeid  The local cmesh id of the cell.
 * \note The cmesh local cell id is a different one as the tree id inside the forest.
 * \return A container with the boundary conditions.
 */
t8_boundary_conditions<std::string_view>
t8_cmesh_get_boundary_conditions (t8_cmesh_t cmesh, t8_locidx_t ltreeid);

/**
 * Retrieves the boundary condition of one face of a cmesh cell.
 * Retrieving all boundary conditions at once via \ref t8_cmesh_get_boundary_conditions will be faster.
 *
 * \param [in] cmesh    The cmesh the cell lives in.
 * \param [in] ltreeid  The local cmesh id of the cell.
 * \param [in] face     The face id of the cell.
 * \note The cmesh local cell id is a different one as the tree id inside the forest.
 * \return The boundary condition of the tree face.
 */
std::string_view
t8_cmesh_get_boundary_condition (t8_cmesh_t cmesh, t8_locidx_t ltreeid, int face);

/**
 * Retrieves the boundary conditions of a forest element.
 *
 * \param [in] forest   The forest the element lives in.
 * \param [in] ltreeid  The local id of the forest tree.
 * \param [in] element  The element.
 * \return A container with the boundary conditions. Note, that only elements faces at the boundary of a
 * tree will have boundary conditions. Internal faces will return an empty optional.
 */
t8_boundary_conditions<std::optional<std::string_view>>
t8_forest_get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/**
 * Retrieves the boundary condition of a face of a forest element.
 * Retrieving all boundary conditions at once via \ref t8_forest_get_boundary_conditions will be faster.
 *
 * \param [in] forest   The forest the element lives in.
 * \param [in] ltreeid  The local id of the forest tree.
 * \param [in] element  The element.
 * \param [in] face     The face id of the element.
 * \return The boundary condition. It will be empty if the element is not touching the boundary of the tree.
 */
std::optional<std::string_view>
t8_forest_get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face);
