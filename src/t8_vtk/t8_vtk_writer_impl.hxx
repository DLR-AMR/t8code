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
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include "t8_forest/t8_forest_types.h"

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

#define T8_FOREST_VTK_QUADRATIC_ELEMENT_MAX_CORNERS 20
/** Lookup table for number of nodes for curved eclasses. */
const int t8_curved_eclass_num_nodes[T8_ECLASS_COUNT] = { 1, 3, 8, 6, 20, 10, 15, 13 };

/** Lookup table for vtk types of curved elements */
const int t8_curved_eclass_vtk_type[T8_ECLASS_COUNT] = { 1, 21, 23, 22, 25, 24, 26, 27 };

/** Map vtk element corners to element reference coordinates. The reference
 * coordinates are defined in such a way, that the linear vtk corners are listed
 * first and then the curved coords. This way, this array can be used for linear
 * vtk elements as well as quadratic vtk elements.
 */
const double t8_forest_vtk_point_to_element_ref_coords[T8_ECLASS_COUNT][T8_FOREST_VTK_QUADRATIC_ELEMENT_MAX_CORNERS][3]
  = { { /* T8_ECLASS_VERTEX */
        { 0, 0, 0 },    { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_LINE */
        { 0, 0, 0 },    { 1, 0, 0 },    { 0.5, 0, 0 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_QUAD */
        { 0, 0, 0 },    { 1, 0, 0 },    { 1, 1, 0 },    { 0, 1, 0 },    { 0.5, 0, 0 },  { 1, 0.5, 0 },  { 0.5, 1, 0 },
        { 0, 0.5, 0 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_TRIANGLE */
        { 0, 0, 0 },    { 1, 0, 0 },    { 1, 1, 0 },    { 0.5, 0, 0 },  { 1, 0.5, 0 },  { 0.5, 0.5, 0 }, { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 },
        { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_HEX */
        { 0, 0, 0 },   { 1, 0, 0 },   { 1, 1, 0 },   { 0, 1, 0 },   { 0, 0, 1 },   { 1, 0, 1 },   { 1, 1, 1 },
        { 0, 1, 1 },   { 0.5, 0, 0 }, { 1, 0.5, 0 }, { 0.5, 1, 0 }, { 0, 0.5, 0 }, { 0.5, 0, 1 }, { 1, 0.5, 1 },
        { 0.5, 1, 1 }, { 0, 0.5, 1 }, { 0, 0, 0.5 }, { 1, 0, 0.5 }, { 1, 1, 0.5 }, { 0, 1, 0.5 } },
      { /* T8_ECLASS_TET */
        { 0, 0, 0 },     { 1, 0, 0 },       { 1, 1, 1 },     { 1, 0, 1 },    { 0.5, 0, 0 },
        { 1, 0.5, 0.5 }, { 0.5, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 1, 0, 0.5 },  { 1, 0.5, 1 },
        { -1, -1, -1 },  { -1, -1, -1 },    { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 },  { -1, -1, -1 },    { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 } },
      { /* T8_ECLASS_PRISM */
        { 0, 0, 0 },   { 1, 0, 0 },     { 1, 1, 0 },    { 0, 0, 1 },    { 1, 0, 1 },     { 1, 1, 1 },   { 0.5, 0, 0 },
        { 1, 0.5, 0 }, { 0.5, 0.5, 0 }, { 0.5, 0, 1 },  { 1, 0.5, 1 },  { 0.5, 0.5, 1 }, { 0, 0, 0.5 }, { 1, 0, 0.5 },
        { 1, 1, 0.5 }, { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 } },
      { /* T8_ECLASS_PYRAMID */
        { 0, 0, 0 },     { 1, 0, 0 },    { 1, 1, 0 },     { 0, 1, 0 },    { 1, 1, 1 },
        { 0.5, 0, 0 },   { 1, 0.5, 0 },  { 0.5, 1, 0 },   { 0, 0.5, 0 },  { 0.5, 0.5, 0.5 },
        { 1, 0.5, 0.5 }, { 1, 1, 0.5 },  { 0.5, 1, 0.5 }, { -1, -1, -1 }, { -1, -1, -1 },
        { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 },  { -1, -1, -1 }, { -1, -1, -1 } } };

static int
t8_get_number_of_vtk_nodes (const t8_element_shape_t eclass, const int curved_flag)
{
  /* use the lookup table of the eclasses. */
  if (curved_flag) {
    return t8_curved_eclass_num_nodes[eclass];
  }
  return t8_eclass_num_vertices[eclass];
}

static void
t8_forest_vtk_get_element_nodes (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int vertex,
                                 const int curved_flag, double *out_coords)
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_shape_t element_shape = scheme->t8_element_shape (element);
  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[element_shape][vertex];
  const int num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
  t8_forest_element_from_ref_coords (forest, ltreeid, element, ref_coords, num_node, out_coords);
}

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

static void
t8_cmesh_tree_to_vtk_cell (t8_cmesh_t cmesh, const t8_locidx_t itree, const t8_gloidx_t offset, const int write_treeid,
                           const int write_mpirank, const int write_level, const int write_element_id,
                           const int curved_flag, const int is_ghost, long int *point_id, int *cellTypes,
                           vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkCellArray> cellArray,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid,
                           vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank, sc_MPI_Comm comm)
{
  vtkSmartPointer<vtkCell> pvtkCell = NULL;

  const t8_element_shape_t element_shape = (t8_element_shape_t) t8_cmesh_get_tree_class (cmesh, itree);
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
  double *coordinates = T8_ALLOC (double, 3 * num_node);

  const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[element_shape][curved_flag];
  t8_geometry_evaluate (cmesh, offset + itree, ref_coords, num_node, coordinates);

  for (int ivertex = 0; ivertex < num_node; ivertex++, (*point_id)++) {
    const size_t offset_3d = 3 * ivertex;
    /* Insert point in the points array */
    points->InsertNextPoint (coordinates[offset_3d], coordinates[offset_3d + 1], coordinates[offset_3d + 2]);
    pvtkCell->GetPointIds ()->SetId (ivertex, *point_id);
  }
  T8_FREE (coordinates);
  cellArray->InsertNextCell (pvtkCell);
  if (curved_flag == 0) {
    cellTypes[itree] = t8_eclass_vtk_type[element_shape];
  }
  else {
    cellTypes[itree] = t8_curved_eclass_vtk_type[element_shape];
  }
  if (write_treeid == 1) {
    const t8_gloidx_t gtree_id = t8_cmesh_get_global_id (cmesh, itree);
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
}

template <typename grid_t>
t8_locidx_t
grid_local_num_elements (const grid_t grid);

template <>
t8_locidx_t
grid_local_num_elements<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_local_num_elements (grid);
}

template <>
t8_locidx_t
grid_local_num_elements<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_local_trees (grid);
}

template <typename grid_t>
t8_locidx_t
grid_local_num_trees (const grid_t grid);

template <>
t8_locidx_t
grid_local_num_trees<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_num_local_trees (grid);
}

template <>
t8_locidx_t
grid_local_num_trees<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_local_trees (grid);
}

template <typename grid_t>
t8_locidx_t
grid_local_num_ghost_trees (const grid_t grid);

template <>
t8_locidx_t
grid_local_num_ghost_trees<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_num_ghost_trees (grid);
}

template <>
t8_locidx_t
grid_local_num_ghost_trees<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_num_ghosts (grid);
}

template <typename grid_t>
t8_gloidx_t
grid_first_local_id (const grid_t grid);

template <>
t8_gloidx_t
grid_first_local_id<t8_forest_t> (const t8_forest_t grid)
{
  return t8_forest_get_first_local_element_id (grid);
}

template <>
t8_gloidx_t
grid_first_local_id<t8_cmesh_t> (const t8_cmesh_t grid)
{
  return t8_cmesh_get_first_treeid (grid);
}

template <typename grid_t>
bool
grid_do_ghosts (const grid_t grid, const int write_ghosts);

template <>
bool
grid_do_ghosts<t8_forest_t> (const t8_forest_t grid, const int write_ghosts)
{
  bool ghosts = write_ghosts;
  if (grid->ghosts == NULL || grid->ghosts->num_ghosts_elements == 0) {
    /* Never write ghost elements if there aren't any */
    ghosts = false;
  }
  T8_ASSERT (grid->ghosts != NULL || !ghosts);
  return ghosts;
}

template <>
bool
grid_do_ghosts<t8_cmesh_t> (const t8_cmesh_t grid, const int write_ghosts)
{
  bool ghosts = write_ghosts;
  if (t8_cmesh_get_num_ghosts (grid) == 0) {
    /* Never write ghost elements if there aren't any */
    ghosts = false;
  }
  return ghosts;
}

template <typename grid_t>
t8_locidx_t
num_cells_to_write (const grid_t grid, const int write_ghosts);

template <>
t8_locidx_t
num_cells_to_write<t8_forest_t> (const t8_forest_t grid, const int write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_forest_get_num_ghost_trees (grid) : 0);
}

template <>
t8_locidx_t
num_cells_to_write<t8_cmesh_t> (const t8_cmesh_t grid, const int write_ghosts)
{
  return grid_local_num_elements (grid) + (write_ghosts ? t8_cmesh_get_num_ghosts (grid) : 0);
}

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
  if (ghosts) {
    t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, t8_forest_ghost_get_tree_class (forest, itree));
    const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_elements (forest, itree);
    for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
      const t8_element_t *element = t8_forest_ghost_get_element (forest, itree, ielem_ghost);
      t8_forest_element_to_vtk_cell (forest, element, scheme, itree + num_local_trees, offset, write_treeid,
                                     write_mpirank, write_level, write_element_id, curved_flag, 1, *elem_id, point_id,
                                     cellTypes, points, cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);
      (*elem_id)++;
    }
  }
  else {
    t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest, itree));
    const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    /* We iterate over all elements in the tree */
    for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);

      t8_forest_element_to_vtk_cell (forest, element, scheme, itree, offset, write_treeid, write_mpirank, write_level,
                                     write_element_id, curved_flag, 0, *elem_id, point_id, cellTypes, points, cellArray,
                                     vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);
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
  t8_cmesh_tree_to_vtk_cell (cmesh, itree, offset, write_treeid, write_mpirank, write_level, write_element_id,
                             curved_flag, ghosts, point_id, cellTypes, points, cellArray, vtk_treeid, vtk_mpirank,
                             comm);
  (*elem_id)++;
  return;
}

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
  const t8_gloidx_t offset = grid_first_local_id (grid);
  t8_gloidx_t elem_id = offset;

  bool do_ghosts = grid_do_ghosts (grid, write_ghosts);
  t8_locidx_t num_cells = num_cells_to_write (grid, do_ghosts);

  int *cellTypes = T8_ALLOC (int, num_cells);
  T8_ASSERT (cellTypes != NULL);

  const t8_locidx_t num_local_trees = grid_local_num_trees (grid);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, write_treeid, write_mpirank, write_level, write_element_id,
                               write_ghosts, curved_flag, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id, cellArray,
                               points, cellTypes, num_local_trees, &elem_id, &point_id, offset, false, itree);
  }
  if (do_ghosts) {
    const t8_locidx_t num_ghost_trees = grid_local_num_ghost_trees (grid);
    for (t8_locidx_t itree_ghost = 0; itree_ghost < num_ghost_trees; itree_ghost++) {
      t8_grid_tree_to_vtk_cells (grid, unstructuredGrid, write_treeid, write_mpirank, write_level, write_element_id,
                                 write_ghosts, curved_flag, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id,
                                 cellArray, points, cellTypes, num_local_trees, &elem_id, &point_id, offset, true,
                                 itree_ghost);
    }
  }

  T8_ASSERT (cellTypes != NULL);
  /* call grid specific function to translate to vtk*/

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