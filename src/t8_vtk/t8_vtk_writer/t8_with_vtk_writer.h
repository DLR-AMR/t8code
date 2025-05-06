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


#include <t8.h>
#include <t8_vtk.h>
#include <t8_forest/t8_forest_types.h>

#include <vtkUnstructuredGrid.h>

T8_EXTERN_C_BEGIN ();


/**
 * Translate a forest into a vtkUnstructuredGrid with respect to the given flags. 
 * This function uses the vtk library. t8code must be configured with
 * "--enable-vtk" in order to use it.
 * \param [in]  forest    The forest.
 * \param[in, out] unstructuredGrid 
 * \param [in]  write_treeid If true, the global tree id is written for each element.
 * \param [in]  write_mpirank If true, the mpirank is written for each element.
 * \param [in]  write_level If true, the refinement level is written for each element.
 * \param [in]  write_element_id If true, the global element id is written for each element.
 * \param [in]  curved_flag If true, write the elements as curved element types from vtk.
 * \param [in]  write_ghosts If true, write out ghost elements as well.
 * \param [in]  num_data  Number of user defined double valued data fields to write.
 * \param [in]  data      Array of t8_vtk_data_field_t of length \a num_data
 *                        providing the user defined per element data.
 *                        If scalar and vector fields are used, all scalar fields
 *                        must come first in the array.
 */
void
t8_forest_to_vtkUnstructuredGrid (t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                  const int write_treeid, const int write_mpirank, const int write_level,
                                  const int write_element_id, const int write_ghosts, const int curved_flag,
                                  const int num_data, t8_vtk_data_field_t *data);

T8_EXTERN_C_END ();
