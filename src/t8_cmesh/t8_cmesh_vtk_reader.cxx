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
#include <t8_cmesh/t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh/t8_cmesh_reader_helper.hxx>
#include <t8_eclass.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

#if T8_WITH_VTK
#include <vtkCellIterator.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
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
#endif
T8_EXTERN_C_BEGIN ();

vtkSmartPointer < vtkUnstructuredGrid >
t8_read_unstructured (const char *filename)
{
  char               *tmp, *extension;
  /*Get the file-extension to decide which reader to use */
  tmp = T8_ALLOC (char, BUFSIZ);
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");
  T8_ASSERT (strcmp (extension, ""));

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
  T8_FREE (tmp);
}

vtkSmartPointer < vtkPolyData > t8_read_poly (const char *filename)
{
  char               *tmp, *extension;
  /*Get the file-extension to decide which reader to use */
  tmp = T8_ALLOC (char, BUFSIZ);
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");
  T8_ASSERT (strcmp (extension, ""));

  /*Read the file */
  if (strcmp (extension, ".ply") == 0) {
    vtkNew < vtkPLYReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, ".vtp") == 0) {
    vtkNew < vtkXMLPolyDataReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, ".obj") == 0) {
    vtkNew < vtkOBJReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, ".stl") == 0) {
    vtkNew < vtkSTLReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, ".vtk") == 0) {
    vtkNew < vtkPolyDataReader > reader;
    reader->SetFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else if (strcmp (extension, ".g") == 0) {
    vtkNew < vtkBYUReader > reader;
    reader->SetGeometryFileName (filename);
    reader->Update ();
    return reader->GetOutput ();
  }
  else {
    t8_global_errorf ("Please use .vtk or .vtu file\n");
  }
  T8_FREE (tmp);
}

/*Construct a cmesh given a filename and a*/
t8_cmesh_t
t8_cmesh_read_from_vtk_unstructured (const char *filename,
                                     const int num_files,
                                     const int compute_face_neigh,
                                     sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
#if T8_WITH_VTK
  /*The Incoming data must be an unstructured Grid */
  int                 cell_type;
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkCellIterator    *cell_it;
  vtkCellData        *cellData;
  int                 num_data_arrays = 0;
  size_t             *data_size;
  double            **tuples;   /*vtk stores data as doubles */
  int                 max_dim = -1;     /*max dimenstion of the cells for geometry */
  int                 num_cell_points, max_cell_points;
  t8_gloidx_t         cell_id;
  vtkSmartPointer < vtkPoints > points =
    vtkSmartPointer < vtkPoints >::New ();
  double             *vertices;

  unstructuredGrid = t8_read_unstructured (filename);

  t8_cmesh_init (&cmesh);
  /* New Iterator to iterate over all cells in the grid */
  cell_it = unstructuredGrid->NewCellIterator ();
  /* Get the Data of the all cells */
  cellData = unstructuredGrid->GetCellData ();
  /* Get the number of data per cell */
  num_data_arrays = cellData->GetNumberOfArrays ();

  /*Prepare data-structures for reading data */
  t8_debugf ("[D] found %i data arrays\n", num_data_arrays);
  //num_data_arrays = 0;
  if (num_data_arrays > 0) {
    int                 tuple_size;
    tuples = T8_ALLOC (double *, num_data_arrays);
    data_size = T8_ALLOC (size_t, num_data_arrays);
    for (int i = 0; i < num_data_arrays; i++) {
      vtkDataArray       *data = cellData->GetArray (i);
      tuple_size = data->GetNumberOfComponents ();
      data_size[i] = sizeof (double) * tuple_size;
      /*Allocate memory for a tuple in array i */
      tuples[i] = T8_ALLOC (double, tuple_size);
    }
  }
  max_cell_points = unstructuredGrid->GetMaxCellSize ();
  /*Allocate maximal possible points to avoid reallocation in the loop */
  vertices = T8_ALLOC (double, 3 * max_cell_points);
  /*Each cell in vtk will be a tree in the cmesh */
  t8_gloidx_t         tree_id = 0;
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal ();
       cell_it->GoToNextCell ()) {
    /*Set the t8_eclass of the cell */
    cell_type = t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    /* Id of the cell in vtk */
    cell_id = cell_it->GetCellId ();
    /*Update dimension of the cmesh. */
    if (max_dim < t8_eclass_to_dimension[cell_type]) {
      max_dim = t8_eclass_to_dimension[cell_type];
    }
    SC_CHECK_ABORTF (cell_type != T8_ECLASS_INVALID,
                     "vtk-cell-type %i not supported by t8code\n",
                     cell_it->GetCellType ());
    t8_cmesh_set_tree_class (cmesh, tree_id, (t8_eclass_t) cell_type);
    /*Set the vertices of the tree */
    num_cell_points = cell_it->GetNumberOfPoints ();
    /*Get the actuall points of the cell */
    points = cell_it->GetPoints ();
    /* For every corner of the cell the the points */
    for (int i = 0; i < num_cell_points; i++) {
      /*Every Point has 3 coords */
      points->GetPoint (i, &vertices[3 * i]);
    }
    /*The order of the vertices in vtk might give a tree with negative volume */
    if (t8_cmesh_tree_vertices_negative_volume
        ((t8_eclass_t) cell_type, vertices, num_cell_points)) {
      t8_cmesh_correct_volume (vertices, (t8_eclass_t) cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_cell_points);
    for (int dtype = 0; dtype < num_data_arrays; dtype++) {
      vtkDataArray       *data = cellData->GetArray (dtype);
      data->GetTuple (cell_id, (double *) tuples[dtype]);
      /*key 0 is reserved for tree-vertices */
      t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), dtype + 1,
                              tuples[dtype], data_size[dtype], 0);
    }
    tree_id++;
  }
  cell_it->Delete ();
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (max_dim);
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_commit (cmesh, comm);
  if (num_data_arrays > 0) {
    T8_FREE (data_size);
    for (int i = num_data_arrays - 1; i >= 0; i--) {
      T8_FREE (tuples[i]);
    }
    T8_FREE (tuples);
  }
#else
  /*TODO: Proper return value to prevent compiler-errors */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
#endif
  T8_FREE (vertices);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_read_from_vtk_poly (const char *filename, const int num_files,
                             const int compute_face_neigh, sc_MPI_Comm comm)
{
  vtkSmartPointer < vtkPolyData > poly_data;
  t8_cmesh_t          cmesh;
  poly_data = t8_read_poly (filename);

  t8_cmesh_init (&cmesh);

}

T8_EXTERN_C_END ();
