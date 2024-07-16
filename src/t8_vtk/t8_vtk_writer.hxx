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

#ifndef T8_VTK_WRITER_HXX
#define T8_VTK_WRITER_HXX

#include <t8_cmesh.h>
#include <t8_vtk.h>
#include <t8.h>
#include "t8_vtk/t8_vtk_writer_impl.hxx"
#include "t8_vtk/t8_vtk_writer_helper.hxx"

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#endif
template <typename grid_t>
class vtk_writer {
 public:
  vtk_writer (const bool write_treeid, const bool write_mpirank, const bool write_level, const bool write_element_id,
              const bool write_ghosts, const bool curved_flag, std::string fileprefix, const int num_data,
              t8_vtk_data_field_t *data, sc_MPI_Comm comm)
    : write_treeid (write_treeid), write_mpirank (write_mpirank), write_level (write_level),
      write_ghosts (write_ghosts), curved_flag (curved_flag), fileprefix (fileprefix), num_data (num_data), data (data),
      comm (comm) {};

  bool
  write (const grid_t grid)
  {
    return write_vtk (grid);
  }

 private:
  /**
 * Construct an unstructuredGrid from either a forest or cmesh. The flags can be used to define what parameters we want to write. 
 * 
 * \param[in] grid a forest or a cmesh
 * \param[in, out] unstructuredGrid an unstructuredGrid that we want to fill with the data of \a grid
 */
  void
  t8_grid_to_vtkUnstructuredGrid (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid)
  {
    T8_ASSERT (grid != NULL);

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New ();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New ();
    vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid = vtkSmartPointer<t8_vtk_gloidx_array_type_t>::New ();
    vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank = vtkSmartPointer<t8_vtk_gloidx_array_type_t>::New ();
    vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level = vtkSmartPointer<t8_vtk_gloidx_array_type_t>::New ();
    vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id = vtkSmartPointer<t8_vtk_gloidx_array_type_t>::New ();

    /*
   * We need the dataArray for writing double valued user defined data in the vtu files.
   * We want to write num_data many timesteps/arrays.
   * We need num_data many vtkDoubleArrays, so we need to allocate storage.
   * Later we call the constructor with: dataArrays[idata]=vtkDoubleArray::New()
   */
    vtkDoubleArray **dataArrays;
    dataArrays = T8_ALLOC (vtkDoubleArray *, num_data);

    long int point_id = 0; /* The id of the point in the points Object. */
    const t8_gloidx_t offset
      = grid_first_local_id (grid); /* offset to take the elements of the previous processes into account*/
    t8_gloidx_t elem_id = offset;

    /* check if we have to write ghosts on this process? */
    bool do_ghosts = grid_do_ghosts (grid, write_ghosts);
    /* compute the number of cells on this process */
    t8_locidx_t num_cells = num_cells_to_write (grid, do_ghosts);

    int *cellTypes = T8_ALLOC (int, num_cells);
    T8_ASSERT (cellTypes != NULL);

    /* Iterate over all trees and translate them. */
    const t8_locidx_t num_local_trees = grid_local_num_trees (grid);
    for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
      t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, this->write_treeid, this->write_mpirank, this->write_level,
                                 this->write_element_id, this->curved_flag, this->comm, vtk_treeid, vtk_mpirank,
                                 vtk_level, vtk_element_id, cellArray, points, cellTypes, num_local_trees, &elem_id,
                                 &point_id, offset, false, itree);
    }
    if (do_ghosts) {
      /* Iterate over all ghost trees and translate them. */
      const t8_locidx_t num_ghost_trees = grid_local_num_ghost_trees (grid);
      for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
        t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, this->write_treeid, this->write_mpirank, this->write_level,
                                   this->write_element_id, this->curved_flag, this->comm, vtk_treeid, vtk_mpirank,
                                   vtk_level, vtk_element_id, cellArray, points, cellTypes, num_local_trees, &elem_id,
                                   &point_id, offset, true, itree_ghost);
      }
    }

    /* construct the unstructuredGrid */
    unstructuredGrid->SetPoints (points);
    unstructuredGrid->SetCells (cellTypes, cellArray);

    if (this->write_treeid) {
      vtk_treeid->SetName ("treeid");
      unstructuredGrid->GetCellData ()->AddArray (vtk_treeid);
    }
    if (this->write_mpirank) {
      vtk_mpirank->SetName ("mpirank");
      unstructuredGrid->GetCellData ()->AddArray (vtk_mpirank);
    }
    if (this->write_level) {
      vtk_level->SetName ("level");
      unstructuredGrid->GetCellData ()->AddArray (vtk_level);
    }
    if (this->write_element_id) {
      vtk_element_id->SetName ("element_id");
      unstructuredGrid->GetCellData ()->AddArray (vtk_element_id);
    }

    /* Write the user defined data fields. 
   * For that we iterate over the idata, set the name, the array
   * and then give this data to the unstructured Grid Object.
   * We differentiate between scalar and vector data.
   */
    const t8_locidx_t num_elements = grid_local_num_elements (grid);
    for (int idata = 0; idata < num_data; idata++) {
      dataArrays[idata] = vtkDoubleArray::New ();
      const int num_components = this->data[idata].type == T8_VTK_SCALAR ? 1 : 3;
      dataArrays[idata]->SetName (this->data[idata].description); /* Set the name of the array */
      dataArrays[idata]->SetNumberOfTuples (num_elements);        /* We want number of tuples=number of elements */
      dataArrays[idata]->SetNumberOfComponents (num_components);  /* Each tuples has 3 values */
      dataArrays[idata]->SetVoidArray (this->data[idata].data, num_elements * num_components, 1);
      unstructuredGrid->GetCellData ()->AddArray (dataArrays[idata]);
    }

    /* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */
    for (int idata = 0; idata < num_data; idata++) {
      dataArrays[idata]->Delete ();
    }

    T8_FREE (cellTypes);
    T8_FREE (dataArrays);
    return;
  }

  bool
  write_vtk (const grid_t grid)
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

  bool write_treeid;
  bool write_mpirank;
  bool write_level;
  bool write_element_id;
  bool write_ghosts;
  bool curved_flag;
  std::string fileprefix;
  int num_data;
  t8_vtk_data_field_t *data;
  sc_MPI_Comm comm;
};

/**
 * \brief 
 * 
 * \param[in] grid 
 * \param[in] fileprefix 
 * \param[in] write_treeid 
 * \param[in] write_mpirank 
 * \param[in] write_level 
 * \param[in] write_element_id 
 * \param[in] curved_flag 
 * \param[in] write_ghosts 
 * \param[in] num_data 
 * \param[in] data 
 * \return true 
 * \return false 
 */
//template <typename grid_t>
//int
//t8_write_vtk_via_API (const grid_t grid, std::string fileprefix, const int write_treeid, const int write_mpirank,
//                      const int write_level, const int write_element_id, const int curved_flag, const int write_ghosts,
//                      const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
//{
//  return t8_write_vtk (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id, curved_flag,
//                       write_ghosts, num_data, data, comm);
//}

#endif /* T8_VTK_WRITER_HXX */