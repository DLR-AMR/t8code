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

/** \file t8_cmesh_vtk_reader.cxx
* Implementation of a Reader for vtk/vtu files using the vtk-library.
* The functions can only be used when t8code is linked with the vtk-library.
*/

#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_vtk_helper.hxx>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkTriangleFilter.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#endif
T8_EXTERN_C_BEGIN ();

/*Construct a cmesh given a filename and a*/
t8_cmesh_t
t8_cmesh_read_from_vtk_unstructured (const char *filename,
                                     const int num_files,
                                     const int compute_face_neigh,
                                     sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  /*The Incoming data must be an unstructured Grid */
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkSmartPointer < vtkCellData > cellData;
  /*Prepare grid for translation */
  unstructuredGrid = t8_read_unstructured (filename);

  /* Get the Data of the all cells */
  cellData = unstructuredGrid->GetCellData ();

  /*Actual translation */
  return t8_vtk_iterate_cells (unstructuredGrid, cellData, comm);

#else
  /*Return empty cmesh if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return t8_cmesh_new_empty (comm, 0, 0);
#endif
}

t8_cmesh_t
t8_cmesh_read_from_vtk_poly (const char *filename, const int num_files,
                             const int compute_face_neigh, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  vtkSmartPointer < vtkPolyData > poly_data;
  vtkSmartPointer < vtkCellArray > cells;
  vtkSmartPointer < vtkCellData > cell_data;
  vtkSmartPointer < vtkPolyData > triangulated;
  vtkNew < vtkTriangleFilter > tri_filter;

  /* Prepare the poly-data for the translation from vtk to t8code.
   * We split all polygons (which are not supported by t8code) to
   * triangles, vertices and lines. */
  poly_data = t8_read_poly (filename);
  tri_filter->SetInputData (poly_data);
  /*PolyVertex to vertex */
  tri_filter->PassVertsOn ();
  /*PolyLines to lines */
  tri_filter->PassLinesOn ();
  tri_filter->Update ();
  triangulated = tri_filter->GetOutput ();

  cell_data = triangulated->GetCellData ();

  return t8_vtk_iterate_cells (triangulated, cell_data, comm);
#else
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return t8_cmesh_new_empty (comm, 0, 0);
#endif
}

T8_EXTERN_C_END ();
