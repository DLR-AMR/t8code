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

/** \file t8_vtk_writer.cxx
 * Implements functions declared in \ref t8_vtk_writer.hxx
 *  or the C interface \ref t8_vtk_writer.h.
 */

#include <t8_vtk/t8_vtk_writer.hxx>
#include <t8_vtk/t8_vtk_writer.h>

#if T8_ENABLE_VTK
#include <vtkUnstructuredGrid.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif

#if T8_ENABLE_VTK
/**
 * Template specialization for forests.
 *
 */
template <>
void
vtk_writer<t8_forest_t>::t8_grid_tree_to_vtk_cells (
  const t8_forest_t forest, [[maybe_unused]] vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkMergePoints> points, int *cellTypes,
  const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
  const bool ghosts, const t8_locidx_t itree, const bool merge_points)
{
  /* For both ghosts and pure-local trees iterate over all elements and translate them into a vtk cell. */
  if (ghosts) {
    const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_leaf_elements (forest, itree);
    for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
      const t8_element_t *element = t8_forest_ghost_get_leaf_element (forest, itree, ielem_ghost);
      this->t8_grid_element_to_vtk_cell (forest, element, itree + num_local_trees, offset, true, *elem_id, point_id,
                                         cellTypes, points, cellArray, vtk_treeid, vtk_mpirank, vtk_level,
                                         vtk_element_id, merge_points);
      (*elem_id)++;
    }
  }
  else {
    const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    /* We iterate over all elements in the tree */
    for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      this->t8_grid_element_to_vtk_cell (forest, element, itree, offset, false, *elem_id, point_id, cellTypes, points,
                                         cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, merge_points);
      (*elem_id)++;

    } /* end of loop over elements */
  }
  return;
}

/**
 * Template specialization for cmeshes.
 *
 */
template <>
void
vtk_writer<t8_cmesh_t>::t8_grid_tree_to_vtk_cells (
  const t8_cmesh_t cmesh, [[maybe_unused]] vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkMergePoints> points, int *cellTypes,
  [[maybe_unused]] const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id,
  const t8_gloidx_t offset, const bool ghosts, const t8_locidx_t itree, const bool merge_points)
{
  /* A cmesh does not have any further elements, we can call the translator directly. */
  this->t8_grid_element_to_vtk_cell (cmesh, NULL, itree, offset, ghosts, *elem_id, point_id, cellTypes, points,
                                     cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, merge_points);
  (*elem_id)++;
  return;
}
#endif /* T8_ENABLE_VTK */

/* Defined here and not in the header, so that the header does not need to
 * include vtkMPI.h, which clashes with OpenCASCADE's Handle(Class) macro. */
template <typename grid_t>
bool
vtk_writer<grid_t>::write_vtk ([[maybe_unused]] const grid_t grid)
{
#if T8_ENABLE_VTK
  T8_ASSERT (!fileprefix.empty ());

  /*
 * Write file: First we construct the unstructured Grid that will store the points and elements. It requires
 * information about the points(coordinates, stored in the points object) and the cells(cellTypes and which points
 * belong to this cell).
 */
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New ();

  std::string mpifilename = fileprefix + std::string (".pvtu");

  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriterObj = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New ();
  t8_grid_to_vtkUnstructuredGrid (grid, unstructuredGrid);
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
#endif /* T8_ENABLE_MPI */
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

// Explicit instantiation: vtk_writer is only used for forests and cmeshes.
template bool
vtk_writer<t8_forest_t>::write_vtk (const t8_forest_t grid);
template bool
vtk_writer<t8_cmesh_t>::write_vtk (const t8_cmesh_t grid);

/**
 * Write a forest to a VTK file in ASCII format.
 */
template <>
bool
vtk_writer<t8_forest_t>::write_ASCII (const t8_forest_t forest)
{
  return t8_forest_vtk_write_ASCII (forest, this->fileprefix.c_str (), this->write_treeid, this->write_mpirank,
                                    this->write_level, this->write_element_id, this->write_ghosts, this->num_data,
                                    this->data);
}

/**
 * Write a cmesh to a VTK file in ASCII format.
 */
template <>
bool
vtk_writer<t8_cmesh_t>::write_ASCII (const t8_cmesh_t cmesh)
{
  return t8_cmesh_vtk_write_ASCII (cmesh, this->fileprefix.c_str ());
}

/* Implementation of the c-interface */
T8_EXTERN_C_BEGIN ();

int
t8_forest_vtk_write_file_via_API (const t8_forest_t forest, const char *fileprefix, const int write_treeid,
                                  const int write_mpirank, const int write_level, const int write_element_id,
                                  const int curved_flag, const int write_ghosts, const int num_data,
                                  t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_with_API (forest);
}

int
t8_forest_vtk_write_file (const t8_forest_t forest, const char *fileprefix, const int write_treeid,
                          const int write_mpirank, const int write_level, const int write_element_id, int write_ghosts,
                          const int num_data, t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, false,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_ASCII (forest);
}

int
t8_cmesh_vtk_write_file_via_API (const t8_cmesh_t cmesh, const char *fileprefix, sc_MPI_Comm comm)
{
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), comm);
  return writer.write_with_API (cmesh);
}

int
t8_cmesh_vtk_write_file (const t8_cmesh_t cmesh, const char *fileprefix)
{
  /* No mpi Communicator is needed for ASCII output*/
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), sc_MPI_COMM_NULL);
  return writer.write_ASCII (cmesh);
}

T8_EXTERN_C_END ();

#if T8_ENABLE_VTK
void
t8_forest_to_vtkUnstructuredGrid (const t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                  const int write_treeid, const int write_mpirank, const int write_level,
                                  const int write_element_id, const int write_ghosts, const int curved_flag,
                                  const int num_data, t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (""), num_data, data, t8_forest_get_mpicomm (forest));
  writer.grid_to_vtkUnstructuredGrid (forest, unstructuredGrid);
}
#endif
