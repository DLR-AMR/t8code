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
#include <vtkCellArray.h>
#include <vtkCellData.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif

#if T8_WITH_VTK
/**
 * Translate a single element from the forest into a vtkCell and fill the vtkArrays with
 * the data related to the element (not element_data).
 */
static void
t8_forest_element_to_vtk_cell (
  t8_forest_t forest, const t8_element_t *element, t8_eclass_scheme_c *scheme, const t8_locidx_t itree,
  const t8_gloidx_t offset, const int write_treeid, const int write_mpirank, const int write_level,
  const int write_element_id, const int curved_flag, const int is_ghost, const int elem_id, long int *point_id,
  int *cellTypes, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> cellArray,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id)
{
  vtkSmartPointer<vtkCell> pvtkCell = NULL;

  const t8_element_shape_t element_shape = scheme->t8_element_shape (element);
  const int num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
  /* depending on the element type we choose the correct vtk cell to insert points to */
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
  double *coordinates = T8_ALLOC (double, 3 * num_node);
  /* Compute coordinates for all vertices inside the domain. */
  t8_forest_vtk_get_element_nodes (forest, itree, element, 0, curved_flag, coordinates);
  /* For each element we iterate over all points */
  for (int ivertex = 0; ivertex < num_node; ivertex++, (*point_id)++) {
    const size_t offset_3d = 3 * ivertex;
    /* Insert point in the points array */
    points->InsertNextPoint (coordinates[offset_3d], coordinates[offset_3d + 1], coordinates[offset_3d + 2]);

    pvtkCell->GetPointIds ()->SetId (ivertex, *point_id);
  }
  T8_FREE (coordinates);
  /* We insert the next cell in the cell array */
  cellArray->InsertNextCell (pvtkCell);
  /*
   * Write current cell Type in the cell Types array at the elem_id index.
   * Depending on the values of the binary inputs write_treeid,
   * write_mpirank and write_element_id we also fill the corresponding
   * arrays with the data we want(treeid,mpirank,element_id).
   * To get the element id, we have to add the local id in the tree 
   * plus theo
   */
  if (curved_flag == 0) {
    cellTypes[elem_id - offset] = t8_eclass_vtk_type[element_shape];
  }
  else {
    cellTypes[elem_id - offset] = t8_curved_eclass_vtk_type[element_shape];
  }
  if (write_treeid == 1) {
    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest, itree);
    if (is_ghost) {
      vtk_treeid->InsertNextValue (-1);
    }
    else {
      vtk_treeid->InsertNextValue (gtree_id);
    }
  }
  if (write_mpirank == 1) {
    vtk_mpirank->InsertNextValue (forest->mpirank);
  }
  if (write_level == 1) {
    vtk_level->InsertNextValue (scheme->t8_element_level (element));
  }
  if (write_element_id == 1) {
    vtk_element_id->InsertNextValue (elem_id);
  }
}

template <typename grid_t>
void
t8_grid_tree_to_vtk_cells (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                           const int write_treeid, const int write_mpirank, const int write_level,
                           const int write_element_id, const int write_ghosts, const int curved_flag,
                           const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm);

template <>
void
t8_grid_tree_to_vtk_cells<t8_forest_t> (const t8_forest_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                        const int write_treeid, const int write_mpirank, const int write_level,
                                        const int write_element_id, const int write_ghosts, const int curved_flag,
                                        const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
  /* 
     * We get the current tree, the scheme for this tree
     * and the number of elements in this tree. We need the vertices of
     * the tree to get the coordinates of the elements later. We need
     * the number of elements in this tree to iterate over all of them.
     */
  t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest, itree));
  const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
  /* We iterate over all elements in the tree */
  for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
    const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
    T8_ASSERT (element != NULL);

    t8_forest_element_to_vtk_cell (forest, element, scheme, itree, offset, write_treeid, write_mpirank, write_level,
                                   write_element_id, curved_flag, 0, elem_id, &point_id, cellTypes, points, cellArray,
                                   vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);

    elem_id++;
  } /* end of loop over elements */
  if (ghosts) {
    /* Get number of ghost-trees */
    const t8_locidx_t num_ghost_trees = t8_forest_ghost_num_trees (forest);
    for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
      /* Get Tree scheme of the ghost tree */
      t8_eclass_scheme_c *scheme
        = t8_forest_get_eclass_scheme (forest, t8_forest_ghost_get_tree_class (forest, itree_ghost));
      const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_elements (forest, itree_ghost);
      for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
        const t8_element_t *element = t8_forest_ghost_get_element (forest, itree_ghost, ielem_ghost);
        t8_forest_element_to_vtk_cell (forest, element, scheme, itree_ghost + num_local_trees, offset, write_treeid,
                                       write_mpirank, write_level, write_element_id, curved_flag, 1, elem_id, &point_id,
                                       cellTypes, points, cellArray, vtk_treeid, vtk_mpirank, vtk_level,
                                       vtk_element_id);
        elem_id++;
      }
    }
  }

  template <typename grid_t>
    void t8_grid_to_vtkUnstructuredGrid
    < (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, const int write_treeid,
       const int write_mpirank, const int write_level, const int write_element_id, const int write_ghosts,
       const int curved_flag, const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
  {
    T8_ASSERT (grid != NULL);

    long int point_id = 0; /* The id of the point in the points Object. */

    const t8_gloidx_t offset = t8_forest_get_first_local_element_id (forest);
    t8_gloidx_t elem_id = offset;

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New ();

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New ();

    int ghosts = write_ghosts;
    if (forest->ghosts == NULL || forest->ghosts->num_ghosts_elements == 0) {
      /* Never write ghost elements if there aren't any */
      ghosts = 0;
    }
    T8_ASSERT (forest->ghosts != NULL || !ghosts);
    /* 
   * The cellTypes Array stores the element types as integers(see vtk doc).
   */
    t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
    if (ghosts) {
      num_elements += t8_forest_get_num_ghosts (forest);
    }
    int *cellTypes = T8_ALLOC (int, num_elements);

    /*
   * We need the vertex coords array to be of the 
   * correct dim. Since it is always the same
   * in one mesh, we take the dim of one element.
   * We add 1 if we look at a vertext(dim=0) because 
   * an array of size 0 is not allowed. 
   * Then we allocate memory, because we do not know
   * beforehand how many entries the array needs.
   */

    /*
   * We have to define the vtkTypeInt64Array that hold 
   * metadata if wanted. 
   */

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

    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

    /* We iterate over all local trees*/
    for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
      /* 
     * We get the current tree, the scheme for this tree
     * and the number of elements in this tree. We need the vertices of
     * the tree to get the coordinates of the elements later. We need
     * the number of elements in this tree to iterate over all of them.
     */
      t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest, itree));
      const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      /* We iterate over all elements in the tree */
      for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
        const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
        T8_ASSERT (element != NULL);

        t8_forest_element_to_vtk_cell (forest, element, scheme, itree, offset, write_treeid, write_mpirank, write_level,
                                       write_element_id, curved_flag, 0, elem_id, &point_id, cellTypes, points,
                                       cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);

        elem_id++;
      } /* end of loop over elements */
    }   /* end of loop over local trees */
    if (ghosts) {
      /* Get number of ghost-trees */
      const t8_locidx_t num_ghost_trees = t8_forest_ghost_num_trees (forest);
      for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
        /* Get Tree scheme of the ghost tree */
        t8_eclass_scheme_c *scheme
          = t8_forest_get_eclass_scheme (forest, t8_forest_ghost_get_tree_class (forest, itree_ghost));
        const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_elements (forest, itree_ghost);
        for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
          const t8_element_t *element = t8_forest_ghost_get_element (forest, itree_ghost, ielem_ghost);
          t8_forest_element_to_vtk_cell (forest, element, scheme, itree_ghost + num_local_trees, offset, write_treeid,
                                         write_mpirank, write_level, write_element_id, curved_flag, 1, elem_id,
                                         &point_id, cellTypes, points, cellArray, vtk_treeid, vtk_mpirank, vtk_level,
                                         vtk_element_id);
          elem_id++;
        }
      }
    }

    unstructuredGrid->SetPoints (points);
    unstructuredGrid->SetCells (cellTypes, cellArray);

    if (write_treeid) {
      vtk_treeid->SetName ("treeid");
      unstructuredGrid->GetCellData ()->AddArray (vtk_treeid);
    }
    if (write_mpirank) {
      vtk_mpirank->SetName ("mpirank");
      unstructuredGrid->GetCellData ()->AddArray (vtk_mpirank);
    }
    if (write_level) {
      vtk_level->SetName ("level");
      unstructuredGrid->GetCellData ()->AddArray (vtk_level);
    }
    if (write_element_id) {
      vtk_element_id->SetName ("element_id");
      unstructuredGrid->GetCellData ()->AddArray (vtk_element_id);
    }

    /* Write the user defined data fields. 
   * For that we iterate over the idata, set the name, the array
   * and then give this data to the unstructured Grid Object.
   * We differentiate between scalar and vector data.
   */
    for (int idata = 0; idata < num_data; idata++) {
      dataArrays[idata] = vtkDoubleArray::New ();
      const int num_components = data[idata].type == T8_VTK_SCALAR ? 1 : 3;
      dataArrays[idata]->SetName (data[idata].description);      /* Set the name of the array */
      dataArrays[idata]->SetNumberOfTuples (num_elements);       /* We want number of tuples=number of elements */
      dataArrays[idata]->SetNumberOfComponents (num_components); /* Each tuples has 3 values */
      dataArrays[idata]->SetVoidArray (data[idata].data, num_elements * num_components, 1);
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

#endif

  template <typename grid_t>
  bool t8_write_vtk_via_API (const grid_t grid, const std::string fileprefix, const int write_treeid,
                             const int write_mpirank, const int write_level, const int write_element_id,
                             const int curved_flag, const int write_ghosts, const int num_data,
                             t8_vtk_data_field_t *data, sc_MPI_Comm comm)
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