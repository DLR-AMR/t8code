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
 * \file t8_cmesh_boundary_conditions.cxx
 * Implementation context for \ref t8_cmesh_boundary_conditions.hxx and \ref t8_cmesh_boundary_conditions_c_interface.h.
 */

#include "t8_cmesh_boundary_conditions.hxx"
#include "t8_cmesh_boundary_conditions_c_interface.h"

#include <span>

/**************************************** FUNCTION DEFINITIONS ****************************************/

t8_boundary_conditions<std::string_view>
t8_cmesh_get_boundary_conditions (t8_cmesh_t cmesh, t8_locidx_t ltreeid)
{
  const detail::t8_cmesh_boundary_condition_handler *handler = t8_cmesh_get_boundary_condition_handler (cmesh);
  SC_CHECK_ABORTF (handler != NULL, "ERROR: Trying to retrieve boundary conditions, even though none were set.\n");
  return handler->get_boundary_conditions (ltreeid);
};

std::string_view
t8_cmesh_get_boundary_condition (t8_cmesh_t cmesh, t8_locidx_t ltreeid, int face)
{
  const detail::t8_cmesh_boundary_condition_handler *handler = t8_cmesh_get_boundary_condition_handler (cmesh);
  SC_CHECK_ABORTF (handler != NULL, "ERROR: Trying to retrieve boundary conditions, even though none were set.\n");
  return handler->get_boundary_condition (ltreeid, face);
};

t8_boundary_conditions<std::optional<std::string_view>>
t8_forest_get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const detail::t8_cmesh_boundary_condition_handler *handler = t8_cmesh_get_boundary_condition_handler (cmesh);
  SC_CHECK_ABORTF (handler != NULL, "ERROR: Trying to retrieve boundary conditions, even though none were set.\n");
  return handler->get_boundary_conditions (forest, ltreeid, element);
};

std::optional<std::string_view>
t8_forest_get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const detail::t8_cmesh_boundary_condition_handler *handler = t8_cmesh_get_boundary_condition_handler (cmesh);
  SC_CHECK_ABORTF (handler != NULL, "ERROR: Trying to retrieve boundary conditions, even though none were set.\n");
  return handler->get_boundary_condition (forest, ltreeid, element, face);
};

/**************************************** C INTERFACE ****************************************/

T8_EXTERN_C_BEGIN ();

void
t8_cmesh_set_boundary_conditions (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const char *boundary_conditions[],
                                  size_t length)
{
  T8_ASSERT (boundary_conditions != nullptr);
  T8_ASSERT (length > 0);
  T8_ASSERT (length <= T8_ECLASS_MAX_FACES);
  const auto boundary_conditions_cpp = std::span<const char *const> { boundary_conditions, length };
  t8_cmesh_set_boundary_conditions (cmesh, gtreeid, boundary_conditions_cpp);
}

void
t8_cmesh_get_boundary_conditions (t8_cmesh_t cmesh, t8_locidx_t ltreeid,
                                  const char *boundary_conditions[T8_ECLASS_MAX_FACES], size_t *length)
{
  const auto boundary_conditions_cpp = t8_cmesh_get_boundary_conditions (cmesh, ltreeid);
  *length = boundary_conditions_cpp.size ();
  for (size_t i_condition = 0; i_condition < *length; ++i_condition) {
    boundary_conditions[i_condition] = boundary_conditions_cpp[i_condition].data ();
  }
};

void
t8_cmesh_get_boundary_condition (t8_cmesh_t cmesh, t8_locidx_t ltreeid, int face,
                                 [[maybe_unused]] const char **boundary_condition)
{
  *boundary_condition = t8_cmesh_get_boundary_condition (cmesh, ltreeid, face).data ();
};

void
t8_forest_get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const char *boundary_conditions[T8_ECLASS_MAX_FACES], size_t *length)
{
  const auto boundary_conditions_cpp = t8_forest_get_boundary_conditions (forest, ltreeid, element);
  *length = boundary_conditions_cpp.size ();
  for (size_t i_condition = 0; i_condition < *length; ++i_condition) {
    if (boundary_conditions_cpp[i_condition].has_value ()) {
      boundary_conditions[i_condition] = boundary_conditions_cpp[i_condition]->data ();
    }
    else {
      boundary_conditions[i_condition] = nullptr;
    }
  }
};

void
t8_forest_get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                  [[maybe_unused]] const char **boundary_condition)
{
  const auto boundary_condition_cpp = t8_forest_get_boundary_condition (forest, ltreeid, element, face);
  if (boundary_condition_cpp.has_value ()) {
    *boundary_condition = boundary_condition_cpp->data ();
  }
  else {
    *boundary_condition = nullptr;
  }
};

T8_EXTERN_C_END ();
