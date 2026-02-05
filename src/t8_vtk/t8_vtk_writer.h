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

/**
 * \file t8_vtk_writer.h
 * This file contains wrappers for the vtk writer functions.
 */

#ifndef T8_VTK_WRITER_C_INTERFACE_H
#define T8_VTK_WRITER_C_INTERFACE_H

#include <t8.h>
#include <t8_vtk.h>
#include <t8_forest/t8_forest_types.h>

#if T8_ENABLE_VTK
#include <vtkUnstructuredGrid.h>
#endif

T8_EXTERN_C_BEGIN ();

#if T8_ENABLE_VTK
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
#endif

/** Write the forest in .pvtu file format. Writes one .vtu file per
 * process and a meta .pvtu file.
 * This function uses the vtk library. t8code must be configured with
 * "--enable-vtk" in order to use it.
 * Currently does not support pyramid elements.
 * \param [in]  forest    The forest.
 * \param [in]  fileprefix The prefix of the output files. The meta file will be named \a fileprefix.pvtu .
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
 * \return  True if successful, false if not (process local).
 * \note If t8code was not configured with vtk, use \ref t8_forest_vtk_write_file
 */
int
t8_forest_vtk_write_file_via_API (t8_forest_t forest, const char *fileprefix, const int write_treeid,
                                  const int write_mpirank, const int write_level, const int write_element_id,
                                  const int curved_flag, const int write_ghosts, const int num_data,
                                  t8_vtk_data_field_t *data);

/** Write the forest in .pvtu file format. Writes one .vtu file per
 * process and a meta .pvtu file.
 * This function writes ASCII files and can be used when
 * t8code is not configure with "--enable-vtk" and
 * \ref t8_forest_vtk_write_file_via_API is not available.
 * \param [in]  forest    The forest.
 * \param [in]  fileprefix  The prefix of the output files.
 * \param [in]  write_treeid If true, the global tree id is written for each element.
 * \param [in]  write_mpirank If true, the mpirank is written for each element.
 * \param [in]  write_level If true, the refinement level is written for each element.
 * \param [in]  write_element_id If true, the global element id is written for each element.
 * \param [in]  write_ghosts If true, each process additionally writes its ghost elements.
 *                           For ghost element the treeid is -1.
 * \param [in]  num_data  Number of user defined double valued data fields to write.
 * \param [in]  data      Array of t8_vtk_data_field_t of length \a num_data
 *                        providing the used defined per element data.
 *                        If scalar and vector fields are used, all scalar fields
 *                        must come first in the array.
 * \return  True if successful, false if not (process local).
 */
int
t8_forest_vtk_write_file (t8_forest_t forest, const char *fileprefix, const int write_treeid, const int write_mpirank,
                          const int write_level, const int write_element_id, int write_ghosts, const int num_data,
                          t8_vtk_data_field_t *data);

/**
 * Write the cmesh in .pvtu file format. Writes one .vtu file per
 * process and a meta .pvtu file.
 * This function uses the vtk library. t8code must be configured with
 * "--enable-vtk" in order to use it.
 * 
 * \param[in] cmesh The cmesh
 * \param[in] fileprefix The prefix of the output files
 * \param[in] comm The communicator to use
 * \return int 
 * \note If t8code was not configured with vtk, use \ref t8_cmesh_vtk_write_file
 */
int
t8_cmesh_vtk_write_file_via_API (t8_cmesh_t cmesh, const char *fileprefix, sc_MPI_Comm comm);

/**
 * Write the cmesh in .pvtu file format. Writes one .vtu file per
 * process and a meta .pvtu file.
 * This function writes ASCII files and can be used when
 * t8code is not configure with "--enable-vtk" and
 * \ref t8_cmesh_vtk_write_file_via_API is not available. 
 * 
 * \param[in] cmesh The cmesh
 * \param[in] fileprefix The prefix of the output files 
 * \return True (nonzero) if successful, false (zero) otherwise 
 */
int
t8_cmesh_vtk_write_file (t8_cmesh_t cmesh, const char *fileprefix);

T8_EXTERN_C_END ();
#endif /* T8_VTK_WRITER_C_INTERFACE_H */
