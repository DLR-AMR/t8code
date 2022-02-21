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

#include <t8_forest_vtk.h>
#include <t8_forest_catalyst.hxx>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#if T8_WITH_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkWedge.h>
#include <vtkQuadraticEdge.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuadraticQuad.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkQuadraticWedge.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkTypeInt64Array.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif
#include <t8.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
/*T8_EXTERN_C_BEGIN ();*/

/* TODO:  */

/* In this file, we have the two functions t8_build_vtk_unstructured_grid,
 * which builds a vtk unstructured grid from a forest, and t8_write_vtk_only,
 * which then writes the vtu/pvtu files of that unstructured grid. The main 
 * reason for this seperation of the original function t8_forest_write_vtk_via_API
 * in t8_forest_vtk.cxx is the use of t8code in combination with the Catalyst 
 * Software. For further information on Catalyst, see: 
 * https://www.paraview.org/Wiki/ParaView/Catalyst/Overview
 * 
 * 
 * 
 * 
 */
/* TODO:  
 * - t8_build_vtk_unstructured_grid return 0, return type is wrong
 * - use of T8_EXTERN_C_BEGIN?
 */

int
t8_write_vtk_only (vtkUnstructuredGrid * unstructuredGrid,
                   const char *fileprefix, t8_forest_t forest)
{
#if T8_WITH_VTK
  int                 freturn = 0;
  char                mpifilename[BUFSIZ];
  snprintf (mpifilename, BUFSIZ, "%s.pvtu", fileprefix);

  vtkSmartPointer < vtkXMLPUnstructuredGridWriter > pwriterObj =
    vtkSmartPointer < vtkXMLPUnstructuredGridWriter >::New ();
/*
 * Get/Set whether the appended data section is base64 encoded. 
 * If encoded, reading and writing will be slower, but the file 
 * will be fully valid XML and text-only. 
 * If not encoded, the XML specification will be violated, 
 * but reading and writing will be fast. The default is to do the encoding.
 * Documentation: https://vtk.org/doc/release/5.0/html/a02260.html#z3560_2
 * 
 * We set the filename of the pvtu file. The filenames of the vtu files
 * are given based on the name of the pvtu file and the process number.
 */
  pwriterObj->EncodeAppendedDataOff ();
  pwriterObj->SetFileName (mpifilename);

/*
 * Since we want to write multiple files, the processes 
 * have to communicate. Therefore, we define the communicator
 * vtk_comm and set it as the communicator. 
 * We have to set a controller for the pwriterObj, 
 * therefore we define the controller vtk_mpi_ctrl.
 */
#if T8_ENABLE_MPI
  vtkSmartPointer < vtkMPICommunicator > vtk_comm =
    vtkSmartPointer < vtkMPICommunicator >::New ();
  vtkMPICommunicatorOpaqueComm vtk_opaque_comm (&forest->mpicomm);
  vtk_comm->InitializeExternal (&vtk_opaque_comm);

  vtkSmartPointer < vtkMPIController > vtk_mpi_ctrl =
    vtkSmartPointer < vtkMPIController >::New ();
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
  pwriterObj->SetNumberOfPieces (forest->mpisize);
  pwriterObj->SetStartPiece (forest->mpirank);
  pwriterObj->SetEndPiece (forest->mpirank);

/* Write the user defined data fields. 
 * For that we iterate over the idata, set the name, the array
 * and then give this data to the unstructured Grid Object.
 * We differentiate between scalar and vector data.
 */

/* We set the input data and write the vtu files. */
  pwriterObj->SetInputData (unstructuredGrid);
  pwriterObj->Update ();
  if (pwriterObj->Write ()) {
    /* Writing failed */
    freturn = 1;
  }
  else {
    t8_errorf ("Error when writing vtk file.\n");
  }

/* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */
/* Return whether writing was successful */
  return freturn;
#else
  t8_global_errorf
    ("Warning: t8code is not linked against vtk library. Vtk output will not be generated.\n");
  t8_global_productionf
    ("Consider calling 't8_forest_write_vtk' or 't8_forest_vtk_write_file' instead.\n");
  return 0;
#endif
}

vtkSmartPointer < vtkUnstructuredGrid >
t8_build_vtk_unstructured_grid (t8_forest_t forest, int write_treeid,
                                int write_mpirank, int write_level,
                                int write_element_id, int curved_flag,
                                int num_data, t8_vtk_data_field_t * data)
{
#if T8_WITH_VTK
  /*Check assertions: forest and fileprefix are not NULL and forest is commited */
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);

  long int            point_id = 0;     /* The id of the point in the points Object. */
  t8_locidx_t         ielement; /* The iterator over elements in a tree. */
  t8_locidx_t         itree, ivertex;
  double              coordinates[3];
  double              vertex_coords[3] = { 0, 0, 0 };
  int                 elem_id = 0;
  t8_locidx_t         num_elements;
  t8_gloidx_t         gtreeid;
  t8_cmesh_t          cmesh;
  int                 num_corners;

/* Since we want to use different element types and a points Array and cellArray 
 * we have to declare these vtk objects. The cellArray stores the Elements.
 * The points and cellArray are needed to store the data we want to write in the Unstructured Grid. 
 */
  vtkNew < vtkPoints > points;
  vtkNew < vtkCellArray > cellArray;
  vtkNew < vtkHexahedron > hexa;
  vtkNew < vtkVertex > vertex;
  vtkNew < vtkLine > line;
  vtkNew < vtkQuad > quad;
  vtkNew < vtkTriangle > tri;
  vtkNew < vtkWedge > prism;
  vtkNew < vtkTetra > tet;

  vtkNew < vtkQuadraticEdge > quadraticedge;
  vtkNew < vtkQuadraticTriangle > quadratictri;
  vtkNew < vtkQuadraticQuad > quadraticquad;
  vtkNew < vtkQuadraticTetra > quadratictet;
  vtkNew < vtkQuadraticHexahedron > quadratichexa;
  vtkNew < vtkQuadraticWedge > quadraticprism;

  /* 
   * The cellTypes Array stores the element types as integers(see vtk doc).
   */
  num_elements = t8_forest_get_local_num_elements (forest);
  int                *cellTypes = T8_ALLOC (int, num_elements);

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

  t8_vtk_gloidx_array_type_t *vtk_treeid = t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_mpirank =
    t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_level = t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_element_id =
    t8_vtk_gloidx_array_type_t::New ();

/*
 * We need the dataArray for writing double valued user defined data in the vtu files.
 * We want to write num_data many timesteps/arrays.
 * We need num_data many vtkDoubleArrays, so we need to allocate storage.
 * Later we call the constructor with: dataArrays[idata]=vtkDoubleArray::New()
 */
  vtkDoubleArray    **dataArrays;
  dataArrays = T8_ALLOC (vtkDoubleArray *, num_data);

  cmesh = t8_forest_get_cmesh (forest);
/* We iterate over all local trees*/
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
/* 
 * We get the current tree, the scheme for this tree
 * and the number of elements in this tree. We need the vertices of
 * the tree to get the coordinates of the elements later. We need
 * the number of elements in this tree to iterate over all of them.
 */
    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    t8_locidx_t         elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    t8_locidx_t         offset =
      t8_forest_get_tree_element_offset (forest, itree);
    /* We iterate over all elements in the tree */
    /* Compute the global tree id */
    gtreeid = t8_forest_global_tree_id (forest, itree);
    for (ielement = 0; ielement < elems_in_tree; ielement++) {
      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      vtkSmartPointer < vtkCell > pvtkCell = NULL;
      t8_element_shape_t  element_shape = scheme->t8_element_shape (element);
      num_corners = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
      /* depending on the element type we choose the correct vtk cell to insert points to */
      if (curved_flag == 0) {
        switch (element_shape) {
        case T8_ECLASS_VERTEX:
          pvtkCell = vertex;
          break;
        case T8_ECLASS_LINE:
          pvtkCell = line;
          break;
        case T8_ECLASS_QUAD:
          pvtkCell = quad;
          break;
        case T8_ECLASS_TRIANGLE:
          pvtkCell = tri;
          break;
        case T8_ECLASS_HEX:
          pvtkCell = hexa;
          break;
        case T8_ECLASS_TET:
          pvtkCell = tet;
          break;
        case T8_ECLASS_PRISM:
          pvtkCell = prism;
          break;
        case T8_ECLASS_PYRAMID:
          SC_CHECK_ABORT (element_shape != T8_ECLASS_PYRAMID,
                          "Pyramids are not supported in vtk output");
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {
        switch (element_shape) {
        case T8_ECLASS_VERTEX:
          pvtkCell = vertex;
          break;
        case T8_ECLASS_LINE:
          pvtkCell = quadraticedge;
          break;
        case T8_ECLASS_QUAD:
          pvtkCell = quadraticquad;
          break;
        case T8_ECLASS_TRIANGLE:
          pvtkCell = quadratictri;
          break;
        case T8_ECLASS_HEX:
          pvtkCell = quadratichexa;
          break;
        case T8_ECLASS_TET:
          pvtkCell = quadratictet;
          break;
        case T8_ECLASS_PRISM:
          pvtkCell = quadraticprism;
          break;
        case T8_ECLASS_PYRAMID:
          SC_CHECK_ABORT (element_shape != T8_ECLASS_PYRAMID,
                          "Pyramids are not supported in vtk output");
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }

      /* For each element we iterate over all points */
      for (ivertex = 0; ivertex < num_corners; ivertex++, point_id++) {
        /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
        if (curved_flag) {
          t8_curved_element_get_reference_node_coords (element, element_shape,
                                                       scheme, ivertex,
                                                       vertex_coords);
        }
        else {
          scheme->t8_element_vertex_reference_coords (element,
                                                      t8_eclass_vtk_corner_number
                                                      [element_shape]
                                                      [ivertex],
                                                      vertex_coords);
        }

        /* Evalute the geometry */
        t8_geometry_evaluate (cmesh, gtreeid, vertex_coords, coordinates);

        /* Insert point in the points array */
        points->InsertNextPoint (coordinates[0], coordinates[1],
                                 coordinates[2]);
        /* Set the point ids to the vtk cell */
        pvtkCell->GetPointIds ()->SetId (ivertex, point_id);
      }                         /* end loop over all vertices of the element */

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

      /* *INDENT-OFF* */
      if(curved_flag==0){
        cellTypes[elem_id] = t8_eclass_vtk_type[element_shape];
      }
      else{
        cellTypes[elem_id] = t8_curved_eclass_vtk_type[element_shape];
      }
      if (write_treeid == 1) {
        vtk_treeid->InsertNextValue (itree);
      }
      if (write_mpirank == 1) {
        vtk_mpirank->InsertNextValue (forest->mpirank);
      }
      if (write_level == 1) {
        vtk_level->InsertNextValue (scheme->t8_element_level (element));
      }
      if (write_element_id == 1) {
        vtk_element_id->InsertNextValue (elem_id + offset +
                                         t8_forest_get_first_local_element_id
                                         (forest));
      /* *INDENT-ON* */
    }
    elem_id++;
  }                             /* end of loop over elements */
}                               /* end of loop over local trees */

  /* 
   * Write file: First we construct the unstructured Grid 
   * that will store the points and elements. It requires
   * information about the points(coordinates, stored in the points object)
   * and the cells(cellTypes and which points belong to this cell) 
   */
vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid =
  vtkSmartPointer < vtkUnstructuredGrid >::New ();
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
  if (data[idata].type == T8_VTK_SCALAR) {
    dataArrays[idata]->SetName (data[idata].description);       /* Set the name of the array */
    dataArrays[idata]->SetVoidArray (data[idata].data, num_elements, 1);        /* We write the data in the array from the input array */
    unstructuredGrid->GetCellData ()->AddArray (dataArrays[idata]);     /* We add the array to the cell data object */
  }
  else {
    dataArrays[idata]->SetName (data[idata].description);       /* Set the name of the array */
    dataArrays[idata]->SetNumberOfTuples (num_elements);        /* We want number of tuples=number of elements */
    dataArrays[idata]->SetNumberOfComponents (3);       /* Each tuples has 3 values */
    dataArrays[idata]->SetVoidArray (data[idata].data, num_elements * 3, 1);
    unstructuredGrid->GetCellData ()->SetVectors (dataArrays[idata]);
  }
}

/* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */

vtk_treeid->Delete ();
vtk_mpirank->Delete ();
vtk_level->Delete ();
vtk_element_id->Delete ();
for (int idata = 0; idata < num_data; idata++) {
  dataArrays[idata]->Delete ();
}

T8_FREE (cellTypes);
T8_FREE (dataArrays);
/* Return whether writing was successful */
return unstructuredGrid;
#else
  t8_global_errorf
    ("Warning: t8code is not linked against vtk library. Vtk output will not be generated.\n");
  t8_global_productionf
    ("Consider calling 't8_forest_write_vtk' or 't8_forest_vtk_write_file' instead.\n");
  return 0;
#endif
}
