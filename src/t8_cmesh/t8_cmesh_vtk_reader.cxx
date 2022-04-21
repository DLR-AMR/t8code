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

#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh/t8_cmesh_reader_helper.hxx>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCellIterator.h>
#include <vtkIdList.h>
#endif
T8_EXTERN_C_BEGIN ();

#if T8_WITH_VTK
/* Given the pointIds in vtk-Order of a cell, its eclass, and a number of face
 * Compute for each point of that face the number of cells that touch this point
 * (via unstructuredGrid->GetPointCells() ). Fill the ids into a vtkIdList and use
 * IntersectWith(), to get the cell-id that touches all points of the face.
 * This cell is unique (face-neighbors are unique). 
 * 
 * This Operation is probably very costly and should be optional.*/
t8_gloidx_t
t8_cmesh_neighbour_at_face (t8_eclass_t eclass, int t8_face_num)
{

}
#endif

/*Construct a cmesh given a filename and a*/
t8_cmesh_t
t8_cmesh_read_from_vtk (const char *filename, const int num_files,
                        const int compute_face_neigh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
#if T8_WITH_VTK
  /*The Incoming data must be an unstructured Grid */
  t8_eclass_t         cell_type;
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkCellIterator    *cell_it;
  char               *tmp, *extension;
  int                 max_dim = 0;      /*max dimenstion of the cells for geometry */
  int                 num_cell_points, max_cell_points;
  t8_gloidx_t         cell_id;
  vtkSmartPointer < vtkPoints > points =
    vtkSmartPointer < vtkPoints >::New ();
  double             *vertices;
  /*Get the file-extensioni to decide which reader to use */
  tmp = T8_ALLOC (char, BUFSIZ);
  strcpy (tmp, filename);
  extension = strtok (tmp, ".");
  extension = strtok (NULL, ".");
  T8_FREE (tmp);
  T8_ASSERT (strcmp (extension, ""));

  /*Read the file */
  if (strcmp (extension, "vtu") == 0) {
    t8_debugf ("[D] use xml unstructured\n");
    vtkSmartPointer < vtkXMLUnstructuredGridReader > reader =
      vtkSmartPointer < vtkXMLUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    unstructuredGrid = reader->GetOutput ();
  }
  else if (strcmp (extension, "vtk") == 0) {
    t8_debugf ("[D] use unstructured\n");
    vtkSmartPointer < vtkUnstructuredGridReader > reader =
      vtkSmartPointer < vtkUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    unstructuredGrid = reader->GetOutput ();
  }
  else {
    t8_global_errorf ("Please use .vtk or .vtu file\n");
  }

  t8_cmesh_init (&cmesh);

  t8_debugf ("[D] num_cells: %lli\n", unstructuredGrid->GetNumberOfCells ());
  /*New Iteratore to iterate over all cells in the grid */
  cell_it = unstructuredGrid->NewCellIterator ();
  max_cell_points = unstructuredGrid->GetMaxCellSize ();
  vertices = T8_ALLOC (double, 3 * max_cell_points);
  t8_gloidx_t         tree_id = 0;
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal ();
       cell_it->GoToNextCell ()) {
    /*Set the t8_eclass of the cell */
    cell_type = t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    T8_ASSERT (cell_type != T8_ECLASS_ZERO);
    t8_debugf ("[D] cell has type %i\n", cell_type);
    t8_cmesh_set_tree_class (cmesh, tree_id, cell_type);

    /*Set the vertices of the tree */
    num_cell_points = cell_it->GetNumberOfPoints ();
    t8_debugf ("[D] cell has %i points\n", num_cell_points);
    points = cell_it->GetPoints ();
    for (int i = 0; i < num_cell_points; i++) {
      points->GetPoint (i, &vertices[3 * i]);
      t8_debugf ("[D] %i: %f %f %f\n", i, vertices[3 * i],
                 vertices[3 * i + 1], vertices[3 * i + 2]);
    }
    /*The order of the vertices in vtk might give a tree with negative volume */
    if (t8_cmesh_tree_vertices_negative_volume
        (cell_type, vertices, num_cell_points)) {
      t8_cmesh_correct_volume (vertices, cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_cell_points);
    if (compute_face_neigh) {
      unstructuredGrid->BuildLinks ();
      cell_id = cell_it->GetCellId ();
      int                 num_faces = cell_it->GetNumberOfFaces ();
      vtkSmartPointer < vtkIdList > faces =
        vtkSmartPointer < vtkIdList >::New ();
      faces = cell_it->GetFaces ();
      t8_debugf ("[D] numfaces %i\n", num_faces);
      t8_debugf ("[D] numfaces %i\n",
                 cell_it->GetFaces ()->GetNumberOfIds ());
    }

    tree_id++;
  }
  t8_cmesh_commit (cmesh, comm);
#else
  /*TODO: Proper return value to prevent compiler-errors */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
#endif
  T8_FREE (vertices);
  return cmesh;
}

T8_EXTERN_C_END ();
