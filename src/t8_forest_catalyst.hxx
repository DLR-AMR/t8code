/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** The two functions t8_build_vtk_unstructured_grid
 * and t8_write_vtk_only can be used as follows:
 * First, you specify the forest in t8code. Then you build
 * the vtkUnstructuredGrid with t8_build_vtk_unstructured_grid.
 * You can then write that vtkUnstructuredGrid with t8_write_vtk_only.
 */

#ifndef T8_FOREST_CATALYST_HXX
#define T8_FOREST_CATALYST_HXX

#include <t8_vtk.h>
#include <t8_forest.h>
#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#endif

/* function declarations */

/** Builds a vtkUnstrucuredGrid from a forest.
 * This function uses the vtk library. t8code must be configured with
 * "--with-vtk" in order to use it.
 * Currently does not support pyramid elements.
 * \param [in]  forest    The forest.
 * \param [in]  write_treeid If true, the global tree id is written for each element.
 * \param [in]  write_mpirank If true, the mpirank is written for each element .
 * \param [in]  write_level If true, the refinement level is written for each element.
 * \param [in]  write_element_id If true, the global element id is written for each element.
 * \param [in]  curved_flag If true, write the elements as curved element types from vtk.
 * \param [in]  num_data  Number of user defined double valued data fields to write.
 * \param [in]  data      Array of t8_vtk_data_field_t of length \a num_data
 *                        providing the used defined per element data.
 *                        If scalar and vector fields are used, all scalar fields
 *                        must come first in the array.
 * \return  A vtkUnstrucuredGrid.
 * \note If t8code was not configured with vtk, use \ref t8_forest_vtk_write_file.
 */
#if T8_WITH_VTK
vtkSmartPointer < vtkUnstructuredGrid >
t8_build_vtk_unstructured_grid (t8_forest_t forest, int write_treeid,
                                int write_mpirank, int write_level,
                                int write_element_id, int curved_flag,
                                int num_data, t8_vtk_data_field_t * data);

/** Writes the vtu file from vtkUnstrucuredGrid with specified fileprefix.
 * This function uses the vtk library. t8code must be configured with
 * "--with-vtk" in order to use it.
 * Currently does not support pyramid elements.
 * \param [in]  forest    The forest.
 * \param [in]  fileprefix The prefix of the output files. The meta file will be named \a fileprefix.pvtu .
 * \return  True if successful, false if not (process local).
 * \note If t8code was not configured with vtk, use \ref t8_forest_vtk_write_file
 */

int
 
 
 
 
 
 
 
 
t8_write_vtk_only (vtkUnstructuredGrid * unstructuredGrid,
                   const char *fileprefix, t8_forest_t forest);
#endif
#endif /* !T8_FOREST_VTK_H */
