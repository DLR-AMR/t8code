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

/** file t8_forest_vtk.h
 */

/* TODO: Document this file */

#include <t8_forest/t8_forest_vtk_helper.hxx>
#include <t8_vec.h>
#include <t8_element_cxx.hxx>
#if T8_WITH_VTK
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
#include <vtkDoubleArray.h>
#include <vtkTypeInt64Array.h>

#endif /* T8_WITH_VTK */

#if T8_WITH_VTK
/* lookup table for number of nodes for curved eclasses. */
const int           t8_curved_eclass_num_nodes[T8_ECLASS_COUNT] =
  { 1, 3, 8, 6, 20, 10, 15, 13 };

/* lookup table for vtk types of curved elements */
const int           t8_curved_eclass_vtk_type[T8_ECLASS_COUNT] =
  { 1, 21, 23, 22, 25, 24, 26, 27 };
#endif

/* 
 * depending on whether we want to write curved or non-curved elements
 * we need the right number of points, so we choose the right lookup table
 */
#if T8_WITH_VTK
static int
t8_get_number_of_vtk_nodes (t8_element_shape_t eclass, int curved_flag)
{
  /* use the lookup table of the eclasses. */
  if (curved_flag) {
    return t8_curved_eclass_num_nodes[eclass];
  }
  return t8_eclass_num_vertices[eclass];
}
#endif

/* If we want to write curved elements, we need to calculate 
 * the reference coordinates. For the vertices(end points)
 * of the elements, we can use t8_element_vertex_reference_coords 
 * to get them. But for curved elements, we also need nodes at the 
 * middle points of lines of elements. We get those coordinates by 
 * adding the vertices and multiplying by 0.5. To get the 
 * correct node, we use e.g. (vertex - 3) % 4, for each 
 * element there is a correct order, therefore we have those
 * formulas. For more information look into the vtk documentation.
 * TODO: Add Pyramids when they are merged into the dev branch.
 * */
#if T8_WITH_VTK
static void
t8_curved_element_get_reference_node_coords (const t8_element_t *elem,
                                             t8_element_shape_t eclass,
                                             t8_eclass_scheme_c *scheme,
                                             int vertex, double *coords)
{
  double              vertex_coords[3] = { 0, 0, 0 };
  int                 i;
  int                 j;

  switch (eclass) {
  case T8_ECLASS_VERTEX:
    scheme->t8_element_vertex_reference_coords (elem,
                                                t8_eclass_vtk_corner_number
                                                [eclass][vertex], coords);
    break;
  case T8_ECLASS_LINE:
    if (vertex < 2) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex - 1],
                                                  vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex - 2],
                                                  coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_QUAD:
    if (vertex < 4) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 4) % 4];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }

    break;
  case T8_ECLASS_TRIANGLE:
    if (0 <= vertex && vertex <= 2) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      i = (vertex - 3) % 3;
      j = (vertex - 2) % 3;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_HEX:
    if (vertex < 8) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (7 < vertex && vertex < 12) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 8) % 4];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 7) % 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else if (11 < vertex && vertex < 16) {
      i = t8_eclass_vtk_corner_number[eclass][((vertex - 8) % 4) + 4];
      j = t8_eclass_vtk_corner_number[eclass][((vertex - 7) % 4) + 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 16];
      j = i + 4;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }

    break;
  case T8_ECLASS_TET:
    if (vertex < 4) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (3 < vertex && vertex < 7) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 4) % 3];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 7];
      j = 3;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_PRISM:
    if (vertex < 6) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (5 < vertex && vertex < 9) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 3];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 2) % 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else if (8 < vertex && vertex < 12) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex % 3) + 3];
      j = t8_eclass_vtk_corner_number[eclass][((vertex + 1) % 3) + 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 12];
      j = t8_eclass_vtk_corner_number[eclass][(vertex % 12) + 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  default:
    scheme->t8_element_vertex_reference_coords (elem,
                                                t8_eclass_vtk_corner_number
                                                [eclass][vertex], coords);
    break;
  }
}
#endif

#if T8_WITH_VTK
void
t8_forest_to_vtkUnstructuredGrid (t8_forest_t forest,
                                  vtkSmartPointer < vtkUnstructuredGrid >
                                  vtkGrid, int write_treeid,
                                  int write_mpirank, int write_level,
                                  int write_element_id, int curved_flag,
                                  int num_data, t8_vtk_data_field_t *data)
{
  /* Since we want to use different element types and a points Array and cellArray 
   * we have to declare these vtk objects. The cellArray stores the Elements.
   * The points and cellArray are needed to store the data we want to write in the Unstructured Grid. 
   */
  long int            point_id = 0;     /* The id of the point in the points Object. */
  t8_locidx_t         elem_id = 0;
  double              coordinates[3];
  double              vertex_coords[3] = { 0, 0, 0 };
  vtkNew < vtkPoints > points;
  vtkNew < vtkCellArray > cellArray;
  vtkNew < vtkHexahedron > hexa;
  vtkNew < vtkVertex > vertex;
  vtkNew < vtkLine > line;
  vtkNew < vtkQuad > quad;
  vtkNew < vtkTriangle > tri;
  vtkNew < vtkWedge > prism;
  vtkNew < vtkTetra > tet;
  vtkNew < vtkPyramid > pyra;

  vtkNew < vtkQuadraticEdge > quadraticedge;
  vtkNew < vtkQuadraticTriangle > quadratictri;
  vtkNew < vtkQuadraticQuad > quadraticquad;
  vtkNew < vtkQuadraticTetra > quadratictet;
  vtkNew < vtkQuadraticHexahedron > quadratichexa;
  vtkNew < vtkQuadraticWedge > quadraticprism;

  /* 
   * The cellTypes Array stores the element types as integers(see vtk doc).
   */
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest);
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

  t8_cmesh_t          cmesh = t8_forest_get_cmesh (forest);

  /* We iterate over all local trees */
  const t8_locidx_t   num_local_trees =
    t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    /* 
     * We get the current tree, the scheme for this tree
     * and the number of elements in this tree. We need the vertices of
     * the tree to get the coordinates of the elements later. We need
     * the number of elements in this tree to iterate over all of them.
     */
    t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest,
                                                              t8_forest_get_tree_class
                                                              (forest,
                                                               itree));
    const t8_locidx_t   elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    const t8_locidx_t   offset =
      t8_forest_get_tree_element_offset (forest, itree);
    /* We iterate over all elements in the tree */
    /* Compute the global tree id */
    const t8_gloidx_t   gtreeid = t8_forest_global_tree_id (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      vtkSmartPointer < vtkCell > pvtkCell = NULL;
      const t8_element_shape_t element_shape =
        scheme->t8_element_shape (element);
      const int           num_node =
        t8_get_number_of_vtk_nodes (element_shape, curved_flag);
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
          pvtkCell = pyra;
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {                    /* curved_flag != 0 */
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
                          "Quadratic Pyramids are not supported in vtk output");
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      for (int ivertex = 0; ivertex < num_node; ivertex++, point_id++) {
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
        /* Evaluate the geometry */
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
        vtk_treeid->InsertNextValue (gtreeid);
      }
      if (write_mpirank == 1) {

        vtk_mpirank->InsertNextValue (t8_forest_get_mpicomm(forest));
      }
      if (write_level == 1) {
        vtk_level->InsertNextValue (scheme->t8_element_level (element));
      }
      if (write_element_id == 1) {
        vtk_element_id->InsertNextValue (elem_id + offset +
                                         t8_forest_get_first_local_element_id
                                         (forest));
      }
      /* *INDENT-ON* */
      elem_id++;
    }                           /* end of loop over elements */
  }                             /* end of loop over local trees */

  /* Fill the vtkUnstructuredGrid */
  vtkGrid->SetPoints (points);
  vtkGrid->SetCells (cellTypes, cellArray);

  /* *INDENT-OFF* */
  if (write_treeid) {
    vtk_treeid->SetName ("treeid");
    vtkGrid->GetCellData()->AddArray(vtk_treeid);
  }
  if (write_mpirank) {
    vtk_mpirank->SetName ("mpirank");
    vtkGrid->GetCellData()->AddArray(vtk_mpirank);
  }
  if (write_level) {
    vtk_level->SetName ("level");
    vtkGrid->GetCellData()->AddArray(vtk_level);
  }
  if (write_element_id) {
    vtk_element_id->SetName ("element_id");
    vtkGrid->GetCellData()->AddArray(vtk_element_id);
  }
  /* *INDENT-ON* */

  /* Write the user defined data fields. 
   * For that we iterate over the idata, set the name, the array
   * and then give this data to the unstructured Grid Object.
   * We differentiate between scalar and vector data.
   */
  for (int idata = 0; idata < num_data; idata++) {
    dataArrays[idata] = vtkDoubleArray::New ();
    if (data[idata].type == T8_VTK_SCALAR) {
      dataArrays[idata]->SetName (data[idata].description);     /* Set the name of the array */
      dataArrays[idata]->SetVoidArray (data[idata].data, num_elements, 1);      /* We write the data in the array from the input array */
      vtkGrid->GetCellData ()->AddArray (dataArrays[idata]);    /* We add the array to the cell data object */
    }
    else {
      dataArrays[idata]->SetName (data[idata].description);     /* Set the name of the array */
      dataArrays[idata]->SetNumberOfTuples (num_elements);      /* We want number of tuples=number of elements */
      dataArrays[idata]->SetNumberOfComponents (3);     /* Each tuples has 3 values */
      dataArrays[idata]->SetVoidArray (data[idata].data, num_elements * 3, 1);  /*  */
      vtkGrid->GetCellData ()->SetVectors (dataArrays[idata]);  /*  */
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
}
#endif
