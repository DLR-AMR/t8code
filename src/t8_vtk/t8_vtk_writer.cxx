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

#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <string>
#include <t8_vtk/t8_vtk_writer.hxx>
#include <t8_vtk.h>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <t8_forest/t8_forest_vtk.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif

#if T8_WITH_VTK
template <>
void
t8_grid_to_vtkUnstructuredGrid<t8_forest_t> (const t8_forest_t grid,
                                             vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                             const int write_treeid, const int write_mpirank, const int write_level,
                                             const int write_element_id, const int write_ghosts, const int curved_flag,
                                             const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
  t8_forest_to_vtkUnstructuredGrid (grid, unstructuredGrid, write_treeid, write_mpirank, write_level, write_element_id,
                                    write_ghosts, curved_flag, num_data, data);
}

template <>
void
t8_grid_to_vtkUnstructuredGrid<t8_cmesh_t> (const t8_cmesh_t grid,
                                            vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                            const int write_treeid, const int write_mpirank, const int write_level,
                                            const int write_element_id, const int write_ghosts, const int curved_flag,
                                            const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
  return;
}

#endif

template <typename grid_t>
bool
t8_write_vtk_via_API (const grid_t grid, const std::string fileprefix, const int write_treeid, const int write_mpirank,
                      const int write_level, const int write_element_id, const int curved_flag, const int write_ghosts,
                      const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  T8_ASSERT (!fileprefix.empty ());

  /* 
   * Write file: First we construct the unstructured Grid 
   * that will store the points and elements. It requires
   * information about the points(coordinates, stored in the points object)
   * and the cells(cellTypes and which points belong to this cell) 
   */
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New ();

  std::string mpifilename = fileprefix + std::string (".pvtu");

  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriterObj = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New ();
  t8_grid_to_vtkUnstructuredGrid (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                  curved_flag, write_ghosts, num_data, data, comm);
  /*
    * Get/Set whether the appended data section is base64 encoded. 
    * If encoded, reading and writing will be slower, but the file 
    * will be fully valid XML and text-only. 
    * If not encoded, the XML specification will be violated, 
    * but reading and writing will be fast. The default is to do the encoding.
    * Documentation: https://vtk.org/doc/release/5.0/html/a02260.html#z3560_2
    */
  pwriterObj->EncodeAppendedDataOff ();

  /* We set the filename of the pvtu file. The filenames of the vtu files
    * are given based on the name of the pvtu file and the process number.
    */
  pwriterObj->SetFileName (mpifilename.c_str ());

/*
   * Since we want to write multiple files, the processes 
   * have to communicate. Therefore, we define the communicator
   * vtk_comm and set it as the communicator. 
   * We have to set a controller for the pwriterObj, 
   * therefore we define the controller vtk_mpi_ctrl.
   */
#if T8_ENABLE_MPI
  vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New ();
  vtkMPICommunicatorOpaqueComm vtk_opaque_comm (&comm);
  vtk_comm->InitializeExternal (&vtk_opaque_comm);

  vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New ();
  vtk_mpi_ctrl->SetCommunicator (vtk_comm);

  pwriterObj->SetController (vtk_mpi_ctrl);
#endif
  /*
    * We set the number of pieces as the number of mpi processes,
    * since we want to write a file for each process. We also
    * need to define a Start and EndPiece for the current
    * process. Then we can set the inputData for the writer:
    * We want to write the unstructured Grid, update the writer
    * and then write.
    * 
    * Note: We could write more than one file per process here, if desired.
    */
  int mpisize;
  int mpirank;
  int mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  pwriterObj->SetNumberOfPieces (mpisize);
  pwriterObj->SetStartPiece (mpirank);
  pwriterObj->SetEndPiece (mpirank);

  /* We set the input data and write the vtu files. */
  pwriterObj->SetInputData (unstructuredGrid);
  pwriterObj->Update ();
  if (pwriterObj->Write ()) {
    return true;
  }
  else {
    t8_errorf ("Error when writing vtk file.\n");
  }

  /* Return whether writing was successful */
  return false;

#else
  t8_global_errorf ("Warning: t8code is not linked against vtk library. Vtk output will not be generated.\n");
  t8_global_productionf ("Consider calling 't8_forest_write_vtk' or 't8_forest_vtk_write_file' instead.\n");
  return false;
#endif
}