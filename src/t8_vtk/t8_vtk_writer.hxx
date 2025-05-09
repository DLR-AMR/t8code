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
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include "t8_forest/t8_forest_types.h"
#include "t8_vtk/t8_vtk_writer_helper.hxx"
#include "t8_vtk/t8_vtk_write_ASCII.hxx"

#include <string>
#include <t8_vtk.h>
#include <t8_types/t8_vec.hxx>

#if T8_ENABLE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkPyramid.h>
#include <vtkWedge.h>
#include <vtkQuadraticEdge.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuadraticQuad.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkQuadraticWedge.h>
#include <vtkQuadraticPyramid.h>
#include <vtkTypeInt64Array.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif /* T8_ENABLE_MPI */
#endif /* T8_ENABLE_VTK */

/**
 * A class that controls the writing of vtk files for cmeshes or forests. 
 * 
 * \tparam grid_t can be a forest or a cmesh. 
 */
template <typename grid_t>
class vtk_writer {
 public:
  /**
   * Construct a new vtk writer object. All parameters are set to false by default. By default no data is used and
   * \a num_data is set to zero. A default \a fileprefix is NOT given. 
   * 
   * \param write_treeid True, if we want to write the tree id of every element.
   * \param write_mpirank True, if we want to write the mpirankof every element.
   * \param write_level True, if we want to write the level of every element. Uses level 0 if used for a cmesh.
   * \param write_element_id True, if we want to write the element id of every element. Ignored if used for a cmesh.
   * \param write_ghosts True, if we want to write the ghost elements, too. 
   * \param curved_flag True, if we want to use quadratic vtk cells. Uses the geometry of the grid to evaluate points between corners. 
   * \param fileprefix The prefix of the output-file.
   * \param num_data The number of data-fields to print.
   * \param data The data to use.
   * \param comm The communicator for parallel output.
   */
  vtk_writer (const bool write_treeid, const bool write_mpirank, const bool write_level, const bool write_element_id,
              const bool write_ghosts, const bool curved_flag, std::string fileprefix, const int num_data,
              t8_vtk_data_field_t *data, sc_MPI_Comm comm)
    : write_treeid (write_treeid), write_mpirank (write_mpirank), write_level (write_level),
      write_element_id (write_element_id), write_ghosts (write_ghosts), curved_flag (curved_flag),
      fileprefix (fileprefix), num_data (num_data), data (data), comm (comm)
  {
  }

  /**
   * Construct a new vtk writer object. All parameters are set to false.
   * 
   * \param[in] fileprefix 
   * \param[in] comm 
   */
  vtk_writer (std::string fileprefix, sc_MPI_Comm comm): fileprefix (fileprefix), comm (comm)
  {
  }

#if T8_ENABLE_VTK
  void
  grid_to_vtkUnstructuredGrid (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid)
  {
    this->t8_grid_to_vtkUnstructuredGrid (grid, unstructuredGrid);
  }
#endif /* T8_ENABLE_VTK */

  /**
   * A vtk-writer function that uses the vtk API.
   * 
   * \param[in] grid The forest or cmesh that is translated.
   * \return true, if writing was successful. 
   * \return false if writing was not successful. 
   */
  bool
  write_with_API (const grid_t grid)
  {
    return write_vtk (grid);
  }

  /**
   * A vtk-writer function that uses the vtk API
   * 
   * \param[in] grid The forest or cmesh that is translated
   * \return true 
   * \return false 
   */
  bool
  write_ASCII (const grid_t grid);

  /**
   * Set the write treeid flag. Set to true, if you want to write the tree id of every element.
   * 
   * \param[in] write_treeid true or false
   */
  inline void
  set_write_treeid (const bool write_treeid)
  {
    this->write_treeid = write_treeid;
  }

  /**
   * Set the write mpirank flag. Set to true, if you want to write the mpirank of every element.
   * 
   * \param[in] write_mpirank true or false
   */
  inline void
  set_write_mpirank (const bool write_mpirank)
  {
    this->write_mpirank = write_mpirank;
  }

  /**
   * Set the write level flag. Set to true, if you want to write the level of every element.
   * 
   * \param[in] write_level true or false
   */
  inline void
  set_write_level (const bool write_level)
  {
    this->write_level = write_level;
  }

  /**
   * Set the write element id flag. Set to true, if you want to write the element id of every element.
   * 
   * \param[in] write_element_id true or false
   */
  inline void
  set_write_element_id (const bool write_element_id)
  {
    this->write_element_id = write_element_id;
  }

  /**
   * Set the write ghosts flag. Set to true, if you want to write the ghost elements, too.
   * 
   * \param[in] write_ghosts true or false
   */
  inline void
  set_write_ghosts (const bool write_ghosts)
  {
    this->write_ghosts = write_ghosts;
  }

  /**
   * Set the curved flag. Set to true, if you want to use quadratic vtk cells. 
   * Uses the geometry of the grid to evaluate points between corners.
   * 
   * \param[in] curved_flag true or false
   */
  inline void
  set_curved_flag (const bool curved_flag)
  {
    this->curved_flag = curved_flag;
  }

  /**
   * Set the fileprefix for the output files.
   * \param[in] fileprefix 
   */
  inline void
  set_fileprefix (std::string fileprefix)
  {
    this->fileprefix = fileprefix;
  }

 private:
#if T8_ENABLE_VTK
  /**
 * Translate a single element from the forest into a vtkCell and fill the vtkArrays with
 * the data related to the element (not element_data).
 * 
 * \tparam grid_t 
 * \param[in] grid A forest or a cmesh.
 * \param element A pointer to an element. Only necessary if a forest is used. Will be ignored if \a grid is a cmesh. 
 * \param[in] itree The local id of the current tree.
 * \param[in] offset offset the ids by the number of elements/trees of the previous processes.
 * \param[in] is_ghost Flag to decide whether we write a ghost element or not.
 * \param elem_id The id for the element to use by vtk. 
 * \param point_id The next id to use to identify vtkpoints.
 * \param[in, out] cellTypes An int array to fill with the type of each element/tree of \a grid
 * \param[in, out] points A vtk Pointarray to fill with  points representing the points in the grid
 * \param[in, out] cellArray A vtk Cellarray to fill with the cells representing the grid.
 * \param[in, out] vtk_treeid A vtk array to fill with the tree ids of \a grid.
 * \param[in, out] vtk_mpirank A vtk array to fill with the mpirank of each element/tree of \a grid.
 * \param[in, out] vtk_level A vtk array to fill with the level of each element/tree of \a grid.
 * \param[in, out] vtk_element_id A vtk array to fill with the id of each element/tree of \a grid.
 */
  void
  t8_grid_element_to_vtk_cell (const grid_t grid, const t8_element_t *element, const t8_locidx_t itree,
                               const t8_gloidx_t offset, const int is_ghost, const int elem_id, long int *point_id,
                               int *cellTypes, vtkSmartPointer<vtkPoints> points,
                               vtkSmartPointer<vtkCellArray> cellArray,
                               vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid,
                               vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                               vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level,
                               vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id)
  {
    vtkSmartPointer<vtkCell> pvtkCell = NULL;

    /* Get the shape of the current element and the respective shape of the vtk_cell. */
    const t8_element_shape_t element_shape = grid_element_shape (grid, itree, element);
    const int num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);

    if (curved_flag == 0) {
      switch (element_shape) {
      case T8_ECLASS_VERTEX:
        pvtkCell = vtkSmartPointer<vtkVertex>::New ();
        break;
      case T8_ECLASS_LINE:
        pvtkCell = vtkSmartPointer<vtkLine>::New ();
        break;
      case T8_ECLASS_QUAD:
        pvtkCell = vtkSmartPointer<vtkQuad>::New ();
        break;
      case T8_ECLASS_TRIANGLE:
        pvtkCell = vtkSmartPointer<vtkTriangle>::New ();
        break;
      case T8_ECLASS_HEX:
        pvtkCell = vtkSmartPointer<vtkHexahedron>::New ();
        break;
      case T8_ECLASS_TET:
        pvtkCell = vtkSmartPointer<vtkTetra>::New ();
        break;
      case T8_ECLASS_PRISM:
        pvtkCell = vtkSmartPointer<vtkWedge>::New ();
        break;
      case T8_ECLASS_PYRAMID:
        pvtkCell = vtkSmartPointer<vtkPyramid>::New ();
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }
    else { /* curved_flag != 0 */
      switch (element_shape) {
      case T8_ECLASS_VERTEX:
        pvtkCell = vtkSmartPointer<vtkVertex>::New ();
        break;
      case T8_ECLASS_LINE:
        pvtkCell = vtkSmartPointer<vtkQuadraticEdge>::New ();
        break;
      case T8_ECLASS_QUAD:
        pvtkCell = vtkSmartPointer<vtkQuadraticQuad>::New ();
        break;
      case T8_ECLASS_TRIANGLE:
        pvtkCell = vtkSmartPointer<vtkQuadraticTriangle>::New ();
        break;
      case T8_ECLASS_HEX:
        pvtkCell = vtkSmartPointer<vtkQuadraticHexahedron>::New ();
        break;
      case T8_ECLASS_TET:
        pvtkCell = vtkSmartPointer<vtkQuadraticTetra>::New ();
        break;
      case T8_ECLASS_PRISM:
        pvtkCell = vtkSmartPointer<vtkQuadraticWedge>::New ();
        break;
      case T8_ECLASS_PYRAMID:
        pvtkCell = vtkSmartPointer<vtkQuadraticPyramid>::New ();
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }

    /* Compute the coordinates of the element/tree. */
    double *coordinates = T8_ALLOC (double, 3 * num_node);

    grid_element_to_coords (grid, itree, element, curved_flag, coordinates, num_node, element_shape);

    for (int ivertex = 0; ivertex < num_node; ivertex++, (*point_id)++) {
      const size_t offset_3d = 3 * ivertex;
      /* Insert the point in the points array. */
      points->InsertNextPoint (coordinates[offset_3d], coordinates[offset_3d + 1], coordinates[offset_3d + 2]);
      pvtkCell->GetPointIds ()->SetId (ivertex, *point_id);
    }
    T8_FREE (coordinates);

    /* Fill the cell array. */
    cellArray->InsertNextCell (pvtkCell);

    /* Write additional information if desired. */
    if (curved_flag == 0) {
      cellTypes[elem_id - offset] = t8_eclass_vtk_type[element_shape];
    }
    else {
      cellTypes[elem_id - offset] = t8_curved_eclass_vtk_type[element_shape];
    }
    if (write_treeid == 1) {
      const t8_gloidx_t gtree_id = tree_local_to_global_id (grid, itree);
      if (is_ghost) {
        vtk_treeid->InsertNextValue (-1);
      }
      else {
        vtk_treeid->InsertNextValue (gtree_id);
      }
    }
    if (write_mpirank == 1) {
      int mpirank;
      int mpiret = sc_MPI_Comm_rank (this->comm, &mpirank);
      SC_CHECK_MPI (mpiret);
      vtk_mpirank->InsertNextValue (mpirank);
    }
    if (write_level == 1) {
      vtk_level->InsertNextValue (grid_element_level (grid, itree, element));
    }
    if (write_element_id == 1) {
      vtk_element_id->InsertNextValue (elem_id);
    }
  }

  /**
 * Iterate over all trees (and if desired ghost trees to) and call the function that translate the tree into.
 * 
 * \tparam grid_t 
 * \param[in] grid A forest or a cmesh.
 * \param[in, out] unstructuredGrid The unstructuredGrid to fill.
 * \param[in, out] vtk_treeid A vtk array to fill with the tree ids of \a grid.
 * \param[in, out] vtk_mpirank A vtk array to fill with the mpirank of each element/tree of \a grid.
 * \param[in, out] vtk_level A vtk array to fill with the level of each element/tree of \a grid.
 * \param[in, out] vtk_element_id A vtk array to fill with the id of each element/tree of \a grid.
 * \param[in, out] cellArray A vtk Cellarray to fill with the cells representing the grid.
 * \param[in, out] points A vtk Pointarray to fill with  points representing the points in the grid.
 * \param[in, out] cellTypes An int array to fill with the type of each element/tree of \a grid.
 * \param[in] num_local_trees The number of local trees.
 * \param[in, out] elem_id The id of the current element. Will be increased after the call, depending on the number of elements processed.
 * \param[in, out] point_id The id of the points. Will be increased after the call, depending on the number of elements processed.
 * \param[in] offset Offset the ids by the number of elements/trees of the previous processes.
 * \param[in] ghosts Flag to decide whether we write a ghost element or not.
 * \param[in] itree The local id of the current tree.
 */
  void
  t8_grid_tree_to_vtk_cells (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
                             vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
                             const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id,
                             const t8_gloidx_t offset, const bool ghosts, const t8_locidx_t itree);

  /**
 * Construct an unstructuredGrid from either a forest or cmesh. The flags can be used to define what parameters we want to write. 
 * 
 * \param[in] grid A forest or a cmesh.
 * \param[in, out] unstructuredGrid An unstructuredGrid that we want to fill with the data of \a grid.
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
   * Later we call the constructor with: dataArrays[idata]=vtkDoubleArray::New().
   */
    vtkDoubleArray **dataArrays;
    dataArrays = T8_ALLOC (vtkDoubleArray *, num_data);

    long int point_id = 0; /* The id of the point in the points Object. */
    const t8_gloidx_t offset
      = grid_first_local_id (grid); /* offset to take the elements of the previous processes into account*/
    t8_gloidx_t elem_id = offset;

    /* Check if we have to write ghosts on this process. */
    bool do_ghosts = grid_do_ghosts (grid, write_ghosts);
    /* Compute the number of cells on this process. */
    t8_locidx_t num_cells = num_cells_to_write (grid, do_ghosts);

    int *cellTypes = T8_ALLOC (int, num_cells);
    T8_ASSERT (cellTypes != NULL);

    /* Allocate VTK Memory for the arrays */
    int iMaxCellSize = 20;
    cellArray->AllocateEstimate(num_cells, iMaxCellSize);
    points->Allocate(num_cells * iMaxCellSize);
    vtk_treeid->Allocate(num_cells);
    vtk_mpirank->Allocate(num_cells);
    vtk_level->Allocate(num_cells);
    vtk_element_id->Allocate(num_cells);

    /* Iterate over all trees and translate them. */
    const t8_locidx_t num_local_trees = grid_local_num_trees (grid);
    for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
      t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, cellArray,
                                 points, cellTypes, num_local_trees, &elem_id, &point_id, offset, false, itree);
    }
    if (do_ghosts) {
      /* Iterate over all ghost trees and translate them. */
      const t8_locidx_t num_ghost_trees = grid_local_num_ghost_trees (grid);
      for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
        t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id,
                                   cellArray, points, cellTypes, num_local_trees, &elem_id, &point_id, offset, true,
                                   itree_ghost);
      }
    }

    /* Construct the unstructuredGrid. */
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

    /* Write the user defined data fields. For that we iterate over the idata, set the name, the array and then give 
     * this data to the unstructured Grid Object.We differentiate between scalar and vector data.
     */
    const t8_locidx_t num_elements = grid_local_num_elements (grid);
    for (int idata = 0; idata < num_data; idata++) {
      dataArrays[idata] = vtkDoubleArray::New ();
      const int num_components = this->data[idata].type == T8_VTK_SCALAR ? 1 : 3;
      dataArrays[idata]->SetName (this->data[idata].description); /* Set the name of the array. */
      dataArrays[idata]->SetNumberOfTuples (num_elements);        /* We want number of tuples=number of elements. */
      dataArrays[idata]->SetNumberOfComponents (num_components);  /* Each tuples has 3 values. */
      dataArrays[idata]->SetVoidArray (this->data[idata].data, num_elements * num_components, 1);
      unstructuredGrid->GetCellData ()->AddArray (dataArrays[idata]);
    }

    /* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */
    for (int idata = 0; idata < num_data; idata++) {
      dataArrays[idata]->Delete ();
    }

    /* Release unused memory in the arrays */
    unstructuredGrid->Squeeze();

    T8_FREE (cellTypes);
    T8_FREE (dataArrays);
    return;
  }
#endif /* T8_ENABLE_VTK */

  /**
   * Write a vtk file given a forest or a cmesh.
   * 
   * \param[in] grid a forest or a cmesh that will be translated into a vtk-file.
   * \return true if writing was successful.
   * \return false if writing was not successful.
   */
  bool
  write_vtk ([[maybe_unused]] const grid_t grid)
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

  bool write_treeid = false;
  bool write_mpirank = false;
  bool write_level = false;
  bool write_element_id = false;
  bool write_ghosts = false;
  bool curved_flag = false;
  std::string fileprefix;
  int num_data = 0;
  t8_vtk_data_field_t *data = NULL;
  sc_MPI_Comm comm;
};

#endif /* T8_VTK_WRITER_HXX */
