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
 * \file t8_cmesh_boundary_conditions_c_interface.h:
 * Public interface for the definition and retrieval of boundary conditions.
 */

#pragma once

T8_EXTERN_C_BEGIN ();

/**
 * Applies boundary conditions to the faces of a cmesh cell.
 *
 * \param [in]  cmesh               The cmesh.
 * \param [in]  gtreeid             The global id of the tree the boundary conditions should be set for.
 * \param [in]  boundary_conditions The boundary conditions to set.
 * \param [in]  length              The length of \a boundary_conditions. Container must have the same length
 *                                  as the eclass of the cell has faces.
 */
void
t8_cmesh_set_boundary_conditions (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const char *boundary_conditions[],
                                  size_t length);

/**
 * Retrieves the boundary conditions of a cmesh cell.
 *
 * \param [in]  cmesh               The cmesh the cell lives in.
 * \param [in]  ltreeid             The local cmesh id of the cell.
 * \param [out] boundary_conditions The boundary conditions of the faces.
 * \param [out] length              The length of \a boundary_conditions.
 * \note The cmesh local cell id is a different one as the tree id inside the forest.
 */
void
t8_cmesh_get_boundary_conditions (t8_cmesh_t cmesh, t8_locidx_t ltreeid,
                                  const char *boundary_conditions[T8_ECLASS_MAX_FACES], size_t *length);

/**
 * Retrieves the boundary condition of one face of a cmesh cell.
 * Retrieving all boundary conditions at once via \ref t8_cmesh_get_boundary_conditions() will be faster.
 *
 * \param [in]  cmesh               The cmesh the cell lives in.
 * \param [in]  ltreeid             The local cmesh id of the cell.
 * \param [in]  face                The face id of the cell.
 * \param [out] boundary_condition  The boundary condition of the tree face.
 * \note The cmesh local cell id is a different one as the tree id inside the forest.
 */
void
t8_cmesh_get_boundary_condition (t8_cmesh_t cmesh, t8_locidx_t ltreeid, int face, const char **boundary_condition);

/**
 * Retrieves the boundary conditions of a forest element.
 *
 * \param [in]  forest              The forest the element lives in.
 * \param [in]  ltreeid             The local id of the forest tree.
 * \param [in]  element             The element.
 * \param [out] boundary_conditions The boundary conditions of the element. String will be nullptr if the elements
 *                                  face is internal; if it does not touch the trees face, since only the tree faces carry boundary
 *                                  conditions. All inner element faces have neighbors anyways.
 * \param [out] length              The length of \a boundary_conditions.
 */
void
t8_forest_get_boundary_conditions (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const char *boundary_conditions[T8_ECLASS_MAX_FACES], size_t *length);

/**
 * Retrieves the boundary condition of a face of a forest element.
 * Retrieving all boundary conditions at once via \ref t8_forest_get_boundary_conditions() will be faster.
 *
 * \param [in] forest               The forest the element lives in.
 * \param [in] ltreeid              The local id of the forest tree.
 * \param [in] element              The element.
 * \param [in] face                 The face id of the element.
 * \param [out] boundary_condition  The boundary condition of the element. Will be nullptr if the elements
 *                                  face is internal; if it does not touch the trees face, since only the tree faces carry boundary
 *                                  conditions. All inner element faces have neighbors anyways.
 */
void
t8_forest_get_boundary_condition (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                  const char **boundary_condition);

T8_EXTERN_C_END ();
