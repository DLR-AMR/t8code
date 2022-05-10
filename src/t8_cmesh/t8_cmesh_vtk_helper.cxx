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

/** This file provied helper functions for vtk-reader
*/
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_cmesh/t8_cmesh_reader_helper.hxx>

#if T8_WITH_VTK
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPolyData.h>
#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>

#include <t8_cmesh/t8_cmesh_vtk_helper.hxx>

vtkSmartPointer < vtkUnstructuredGrid >
t8_read_unstructured (const char *filename)
{
  char                tmp[BUFSIZ], *extension;
  /*Get the file-extension to decide which reader to use */
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");

  /*Read the file */
  if (strcmp (extension, "vtu") == 0) {
    vtkSmartPointer < vtkXMLUnstructuredGridReader > reader =
      vtkSmartPointer < vtkXMLUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkSmartPointer < vtkUnstructuredGridReader > reader =
      vtkSmartPointer < vtkUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else {
    t8_global_errorf ("Please use .vtk or .vtu file\n");
    return NULL;
  }
}

vtkSmartPointer < vtkPolyData > t8_read_poly (const char *filename)
{
  char                tmp[BUFSIZ], *extension;
  /*Get the file-extension to decide which reader to use */
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");
  T8_ASSERT (strcmp (extension, ""));
  t8_debugf ("[D] %s\n", extension);
  /*Read the file depending on the extension */
  if (strcmp (extension, "ply") == 0) {
    vtkNew < vtkPLYReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "vtp") == 0) {
    vtkNew < vtkXMLPolyDataReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "obj") == 0) {
    vtkNew < vtkOBJReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "stl") == 0) {
    vtkNew < vtkSTLReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkNew < vtkPolyDataReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, "g") == 0) {
    vtkNew < vtkBYUReader > reader;
    reader->SetGeometryFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else {
    t8_global_errorf ("Please use .ply, .vtp, .obj, .stl, .vtk or .g file\n");
    return NULL;
  }
}

t8_cmesh_t
t8_vtk_iterate_cells (vtkSmartPointer < vtkDataSet > cells,
                      vtkSmartPointer < vtkCellData > cell_data,
                      sc_MPI_Comm comm)
{
  vtkCellIterator    *cell_it;
  vtkSmartPointer < vtkPoints > points;
  t8_gloidx_t         cell_id;
  double             *vertices;
  int                 cell_type;
  int                 num_data_arrays, num_points;
  double            **tuples;
  size_t             *data_size;
  t8_gloidx_t         tree_id = 0;
  int                 max_dim = -1;
  t8_cmesh_t          cmesh;
  int                 max_cell_points = -1;

  max_cell_points = cells->GetMaxCellSize ();
  T8_ASSERT (max_cell_points > 0);
  t8_cmesh_init (&cmesh);
  vertices = T8_ALLOC (double, 3 * max_cell_points);
  /*Get cell iterator */
  cell_it = cells->NewCellIterator ();
  /* get the number of data-arrays per cell */
  num_data_arrays = cell_data->GetNumberOfArrays ();
  T8_ASSERT (num_data_arrays >= 0);
  t8_debugf ("[D] read %i data-arrays\n", num_data_arrays);
  /*Prepare attributes */
  if (num_data_arrays > 0) {
    int                 tuple_size;
    tuples = T8_ALLOC (double *, num_data_arrays);
    data_size = T8_ALLOC (size_t, num_data_arrays);
    for (int i = 0; i < num_data_arrays; i++) {
      vtkDataArray       *data = cell_data->GetArray (i);
      tuple_size = data->GetNumberOfComponents ();
      data_size[i] = sizeof (double) * tuple_size;
      /*Allocate memory for a tuple in array i */
      tuples[i] = T8_ALLOC (double, tuple_size);
    }
  }

  /*Iterate over all cells */
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal ();
       cell_it->GoToNextCell ()) {

    /*Set the t8_eclass of the cell */
    cell_type = t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    SC_CHECK_ABORTF (t8_eclass_is_valid ((t8_eclass_t) cell_type),
                     "vtk-cell-type %i not supported by t8code\n", cell_type);
    t8_cmesh_set_tree_class (cmesh, tree_id, (t8_eclass_t) cell_type);
    /*Get the points of the cell */
    num_points = cell_it->GetNumberOfPoints ();
    T8_ASSERT (num_points > 0);
    points = cell_it->GetPoints ();
    for (int i = 0; i < num_points; i++) {
      points->GetPoint (i, &vertices[3 * i]);
    }
    /*The order of the vertices in vtk might give a tree with negative volume */
    if (t8_cmesh_tree_vertices_negative_volume
        ((t8_eclass_t) cell_type, vertices, num_points)) {
      t8_cmesh_correct_volume (vertices, (t8_eclass_t) cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_points);

    /*Get and set the data of each cell */
    for (int dtype = 0; dtype < num_data_arrays; dtype++) {
      cell_id = cell_it->GetCellId ();
      vtkDataArray       *data = cell_data->GetArray (dtype);
      data->GetTuple (cell_id, tuples[dtype]);
      t8_cmesh_set_attribute (cmesh, cell_id, t8_get_package_id (), dtype + 1,
                              tuples[dtype], data_size[dtype], 0);
    }
    /*Check geometry-dimension */
    if (max_dim < cell_it->GetCellDimension ()) {
      max_dim = cell_it->GetCellDimension ();
    }
    tree_id++;
  }
  t8_debugf ("[D] read %li trees\n", tree_id++);
  /*Set the geometry */
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (max_dim);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);

  /*Clean-up */
  cell_it->Delete ();
  if (num_data_arrays > 0) {
    T8_FREE (data_size);
    for (int i = num_data_arrays - 1; i >= 0; i--) {
      T8_FREE (tuples[i]);
    }
    T8_FREE (tuples);
  }
  T8_FREE (vertices);
  return cmesh;
}

#endif
