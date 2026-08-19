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
 * \file t8_cmesh_boundary_condition_handler_types.h
 * Implements functionality for working with private headers and c types.
 */

#pragma once

#include <t8_cmesh/t8_cmesh.h>

#ifdef __cplusplus
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler.hxx>

/** This typedef is used for the opaque pointers to the handler.
 * We need it so that we can use t8_cmesh_boundary_condition_handler_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling).
 * We have one cpp version with the correct namespace and one c version pointing to nothing.
 * TODO: Delete this when the cmesh is a proper cpp class.
 */
typedef struct detail::t8_cmesh_boundary_condition_handler t8_cmesh_boundary_condition_handler_c;

#else

/** This typedef is used for the opaque pointers to the handler.
 * We need it so that we can use t8_cmesh_boundary_condition_handler_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling).
 * We have one cpp version with the correct namespace and one c version pointing to nothing.
 * TODO: Delete this when the cmesh is a proper cpp class.
 */
typedef struct t8_cmesh_boundary_condition_handler t8_cmesh_boundary_condition_handler_c;

#endif

T8_EXTERN_C_BEGIN ();

/**
 * Returns the boundary condition handler of the cmesh.
 * This is needed because we implement templated functions which need to access the handler, but the cmesh type is not installed.
 * \param [in] cmesh  The cmesh
 * \return            A pointer to the boundary condition handler. nullptr if none was set.
 */
t8_cmesh_boundary_condition_handler_c *
t8_cmesh_get_boundary_condition_handler (t8_cmesh_t cmesh);

/**
 * Adds a boundary condition handler to a cmesh. The cmesh cannot have a handler yet.
 * \param [in,out] cmesh  The cmesh
 * \return                A pointer to the newly created handler inside the cmesh.
 */
t8_cmesh_boundary_condition_handler_c *
t8_cmesh_add_boundary_condition_handler (t8_cmesh_t cmesh);

T8_EXTERN_C_END ();
