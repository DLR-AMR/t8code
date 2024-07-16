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

#ifndef T8_VTK_WRITER_IMPL_HXX
#define T8_VTK_WRITER_IMPL_HXX

#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include "t8_forest/t8_forest_types.h"
#include "t8_vtk/t8_vtk_writer_helper.hxx"

#include <string>
#include <t8_vtk.h>
#include <t8_element.hxx>
#include <t8_vec.h>

#if T8_WITH_VTK
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
#endif
#endif

#if T8_WITH_VTK

/**
 * Translate a single element from the forest into a vtkCell and fill the vtkArrays with
 * the data related to the element (not element_data).
 * 
 * \tparam grid_t 
 * \param[in] grid a forest or a cmesh
 * \param element 
 * \param[in] itree the local id of the current tree
 * \param[in] offset offset the ids by the number of elements/trees of the previous processes
 * \param[in] write_treeid flag to decide if we write the treeid
 * \param[in] write_mpirank flag to decide if we write the mpirank
 * \param[in] write_level flag to decide if we write the level
 * \param[in] write_element_id flag to decide if we write the element_id
 * \param[in] curved_flag flag to decide if we write the elements in a curved style.
 * \param[in] is_ghost flag to decide whether we write a ghost element or not
 * \param elem_id the id of the current element
 * \param point_id the next id to use to identiry vtkpoints
 * \param[in, out] cellTypes an int array to fill with the type of each element/tree of \a grid
 * \param[in, out] points a vtk Pointarray to fill with  points representing the points in the grid
 * \param[in, out] cellArray a vtk Cellarray to fill with the cells representing the grid.
 * \param[in, out] vtk_treeid a vtk array to fill with the tree ids of \a grid
 * \param[in, out] vtk_mpirank a vtk array to fill with the mpirank of each element/tree of \a grid
 * \param[in, out] vtk_level a vtk array to fill with the level of each element/tree of \a grid
 * \param[in, out] vtk_element_id a vtk array to fill with the id of each element/tree of \a grid
 * \param[in] comm the communicator to use
 */
template <typename grid_t>
void
t8_grid_element_to_vtk_cell (const grid_t grid, const t8_element_t *element, const t8_locidx_t itree,
                             const t8_gloidx_t offset, const int write_treeid, const int write_mpirank,
                             const int write_level, const int write_element_id, const int curved_flag,
                             const int is_ghost, const int elem_id, long int *point_id, int *cellTypes,
                             vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> cellArray,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level,
                             vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id, sc_MPI_Comm comm)
{
  vtkSmartPointer<vtkCell> pvtkCell = NULL;

  /* Get the shape of the current element and the respective shape of the vtk_cell*/
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

  /* Compute the coordinates of the element/tree */
  double *coordinates = T8_ALLOC (double, 3 * num_node);

  grid_element_to_coords (grid, itree, element, curved_flag, coordinates, num_node, element_shape);

  for (int ivertex = 0; ivertex < num_node; ivertex++, (*point_id)++) {
    const size_t offset_3d = 3 * ivertex;
    /* Insert point in the points array */
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
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
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
 * Iterate over all trees (and if desired ghost trees to) and call the function that translate the tree into
 * 
 * \tparam grid_t 
 * \param[in] grid a forest or a cmesh
 * \param[in, out] unstructuredGrid the unstructuredGrid to fill
 * \param[in] write_treeid flag to decide if we write the treeid
 * \param[in] write_mpirank flag to decide if we write the mpirank
 * \param[in] write_level flag to decide if we write the level
 * \param[in] write_element_id flag to decide if we write the element_id
 * \param[in] write_ghosts flag to decide if we write the ghosts
 * \param[in] curved_flag flag to decide if we write the elements in a curved style. 
 * \param[in] comm the communicator to use
 * \param[in, out] vtk_treeid a vtk array to fill with the tree ids of \a grid
 * \param[in, out] vtk_mpirank a vtk array to fill with the mpirank of each element/tree of \a grid
 * \param[in, out] vtk_level a vtk array to fill with the level of each element/tree of \a grid
 * \param[in, out] vtk_element_id a vtk array to fill with the id of each element/tree of \a grid
 * \param[in, out] cellArray a vtk Cellarray to fill with the cells representing the grid.
 * \param[in, out] points a vtk Pointarray to fill with  points representing the points in the grid
 * \param[in, out] cellTypes an int array to fill with the type of each element/tree of \a grid
 * \param[in] num_local_trees the number of local trees 
 * \param[in, out] elem_id the id of the current element. Will be increased after the call, depending on the number of elements processed
 * \param[in, out] point_id the id of the points. Will be increased after the call, depending on the number of elements processed
 * \param[in] offset offset the ids by the number of elements/trees of the previous processes
 * \param[in] ghosts flag to decide whether we write a ghost element or not
 * \param[in] itree the local id of the current tree
 */
template <typename grid_t>
void
t8_grid_tree_to_vtk_cells (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                           const int write_treeid, const int write_mpirank, const int write_level,
                           const int write_element_id, const int curved_flag, sc_MPI_Comm comm,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
                           vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
                           const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id,
                           const t8_gloidx_t offset, const bool ghosts, const t8_locidx_t itree);

template <>
void
t8_grid_tree_to_vtk_cells<t8_forest_t> (
  const t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, const int write_treeid,
  const int write_mpirank, const int write_level, const int write_element_id, const int curved_flag, sc_MPI_Comm comm,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
  const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
  const bool ghosts, const t8_locidx_t itree)
{
  /* for both ghosts and pure-local trees iterate over all elements and translate them into a vtk cell. */
  if (ghosts) {
    const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_elements (forest, itree);
    for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
      const t8_element_t *element = t8_forest_ghost_get_element (forest, itree, ielem_ghost);
      t8_grid_element_to_vtk_cell (forest, element, itree + num_local_trees, offset, write_treeid, write_mpirank,
                                   write_level, write_element_id, curved_flag, true, *elem_id, point_id, cellTypes,
                                   points, cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, comm);
      (*elem_id)++;
    }
  }
  else {
    const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    /* We iterate over all elements in the tree */
    for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      t8_grid_element_to_vtk_cell (forest, element, itree, offset, write_treeid, write_mpirank, write_level,
                                   write_element_id, curved_flag, true, *elem_id, point_id, cellTypes, points,
                                   cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, comm);
      (*elem_id)++;

    } /* end of loop over elements */
  }
  return;
}

template <>
void
t8_grid_tree_to_vtk_cells<t8_cmesh_t> (
  const t8_cmesh_t cmesh, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, const int write_treeid,
  const int write_mpirank, const int write_level, const int write_element_id, const int curved_flag, sc_MPI_Comm comm,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
  const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
  const bool ghosts, const t8_locidx_t itree)
{
  /* a cmesh does not have any further elements, we can call the translatore directly. */
  t8_grid_element_to_vtk_cell (cmesh, NULL, itree, offset, write_treeid, write_mpirank, write_level, write_element_id,
                               curved_flag, ghosts, *elem_id, point_id, cellTypes, points, cellArray, vtk_treeid,
                               vtk_mpirank, vtk_level, vtk_element_id, comm);
  (*elem_id)++;
  return;
}

/**
 * Construct an unstructuredGrid from either a forest or cmesh. The flags can be used to define what parameters we want to write. 
 * 
 * \tparam grid_t 
 * \param[in] grid a forest or a cmesh
 * \param[in, out] unstructuredGrid an unstructuredGrid that we want to fill with the data of \a grid
 * \param[in] write_treeid flag to decide if we write the treeid
 * \param[in] write_mpirank flag to decide if we write the mpirank
 * \param[in] write_level flag to decide if we write the level
 * \param[in] write_element_id flag to decide if we write the element_id
 * \param[in] write_ghosts flag to decide if we write the ghosts
 * \param[in] curved_flag flag to decide if we write the elements in a curved style. 
 * \param[in] num_data the number of data fields. 
 * \param[in] data the data fields. 
 * \param[in] comm the communicator to use. 
 */
template <typename grid_t>
void
t8_grid_to_vtkUnstructuredGrid (const grid_t grid, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                const int write_treeid, const int write_mpirank, const int write_level,
                                const int write_element_id, const int write_ghosts, const int curved_flag,
                                const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
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
    t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, write_treeid, write_mpirank, write_level, write_element_id,
                               curved_flag, comm, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, cellArray, points,
                               cellTypes, num_local_trees, &elem_id, &point_id, offset, false, itree);
  }
  if (do_ghosts) {
    /* Iterate over all ghost trees and translate them. */
    const t8_locidx_t num_ghost_trees = grid_local_num_ghost_trees (grid);
    for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
      t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, write_treeid, write_mpirank, write_level, write_element_id,
                                 curved_flag, comm, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, cellArray,
                                 points, cellTypes, num_local_trees, &elem_id, &point_id, offset, true, itree_ghost);
    }
  }

  /* construct the unstructuredGrid */
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
  const t8_locidx_t num_elements = grid_local_num_elements (grid);
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
int
t8_write_vtk (const grid_t grid, std::string fileprefix, const int write_treeid, const int write_mpirank,
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
  t8_grid_to_vtkUnstructuredGrid<grid_t> (grid, unstructuredGrid, write_treeid, write_mpirank, write_level,
                                          write_element_id, write_ghosts, curved_flag, num_data, data, comm);
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

#endif /* T8_VTK_WRITER_IMPL_HXX */