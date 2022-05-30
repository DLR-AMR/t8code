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

/** \This file provides header for functions that can help with vtk-reader 
*/

#ifndef T8_CMESH_VTK_HELPER
#define T8_CMESH_VTK_HELPER

#if T8_WITH_VTK
#include<vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>

/* look_up Table to transform vtkCellType into T8_ECLASS
   T8_ECLASS_INVALID, if a the vtkCellType is not supported by t8code.
   see https://vtk.org/doc/nightly/html/vtkCellType_8h.html to check.*/
const t8_eclass_t   t8_cmesh_vtk_type_to_t8_type[82] = {
  T8_ECLASS_INVALID, T8_ECLASS_VERTEX, T8_ECLASS_INVALID, T8_ECLASS_LINE,
  T8_ECLASS_INVALID, T8_ECLASS_TRIANGLE, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_TET,
  T8_ECLASS_HEX, T8_ECLASS_HEX, T8_ECLASS_PRISM, T8_ECLASS_PYRAMID,
  T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID
};

/** iterate over a vtkDataSet via a Celliterator and construct a tree
 * for every cell. All trees are then commited in a cmesh. For each cell
 * add the CellData that lays on the cell (if existent). 
 * 
 * TODO: Use parallel 
 * 
 * \param [in]      cells             The Input cells
 * \param [in]      cell_data         Data lying on the cells
 * \param [in]      comm              The communicator to use
 * \return t8_cmesh_t                 The cmesh constructed using the \a cells.
 */
t8_cmesh_t          t8_vtk_iterate_cells (vtkSmartPointer < vtkDataSet >
                                          cells,
                                          vtkSmartPointer < vtkCellData >
                                          cell_data, sc_MPI_Comm comm);

/**
 * Read the Poly-data of a file containing vtkPolyData
 * \param [in]      filename            The file containing the Data
 * \returns         vtkSmartPointer<vtkPolyData>    A pointer to vtkPolyData
 */
vtkSmartPointer < vtkPolyData > t8_read_poly (const char *filename);

/**
 * Read the unstructured grid of a file containing vtkPolyData
 * \param [in]      filename            The file containing the Data
 * \returns         vtkSmartPointer<vtkUnstructuredGrid>    A pointer to vtkPolyData
 */
vtkSmartPointer < vtkUnstructuredGrid >
t8_read_unstructured (const char *filename);
#endif

#endif /* T8_CMESH_VTK_HELPER */
