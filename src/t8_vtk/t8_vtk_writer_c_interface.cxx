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

#include <t8_vtk/t8_vtk_writer_c_interface.h>
#include <t8_vtk/t8_vtk_writer.hxx>
#include <string>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>

T8_EXTERN_C_BEGIN ();

int
t8_forest_vtk_write_file_via_API (t8_forest_t forest, const char *fileprefix, const int write_treeid,
                                  const int write_mpirank, const int write_level, const int write_element_id,
                                  const int curved_flag, const int write_ghosts, const int num_data,
                                  t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_with_API (forest);
}

int
t8_forest_vtk_write_file (t8_forest_t forest, const char *fileprefix, const int write_treeid, const int write_mpirank,
                          const int write_level, const int write_element_id, int write_ghosts, const int num_data,
                          t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, false,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_ASCII (forest);
}

int
t8_cmesh_vtk_write_file_via_API (t8_cmesh_t cmesh, const char *fileprefix, sc_MPI_Comm comm)
{
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), comm);
  return writer.write_with_API (cmesh);
}

int
t8_cmesh_vtk_write_file (t8_cmesh_t cmesh, const char *fileprefix)
{
  /* No mpi Communicator is needed for ASCII output*/
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), sc_MPI_COMM_NULL);
  return writer.write_ASCII (cmesh);
}

#if T8_WITH_VTK
void
t8_forest_to_vtkUnstructuredGrid (t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                  const int write_treeid, const int write_mpirank, const int write_level,
                                  const int write_element_id, const int write_ghosts, const int curved_flag,
                                  const int num_data, t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (""), num_data, data, t8_forest_get_mpicomm (forest));
  writer.grid_to_vtkUnstructuredGrid (forest, unstructuredGrid);
}
#endif

T8_EXTERN_C_END ();