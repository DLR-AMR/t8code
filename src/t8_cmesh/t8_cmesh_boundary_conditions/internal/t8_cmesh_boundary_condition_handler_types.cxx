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
 * \file t8_cmesh_boundary_condition_handler_types.cxx
 * Implements functionality for working with private headers and c types.
 */

#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler_types.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>

T8_EXTERN_C_BEGIN ();

t8_cmesh_boundary_condition_handler_c *
t8_cmesh_get_boundary_condition_handler (t8_cmesh_t cmesh)
{
  return cmesh->boundary_condition_handler;
}

t8_cmesh_boundary_condition_handler_c *
t8_cmesh_add_boundary_condition_handler (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh->boundary_condition_handler == nullptr);
  cmesh->boundary_condition_handler = new detail::t8_cmesh_boundary_condition_handler (cmesh);
  return cmesh->boundary_condition_handler;
}

T8_EXTERN_C_END ();
