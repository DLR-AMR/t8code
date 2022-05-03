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

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCellIterator.h>
#include <vtkIdList.h>
#endif
T8_EXTERN_C_BEGIN ();

#if T8_WITH_VTK
/* Given the point-Ids of the cell with id cell_id and a face-number of that 
 * cell, we compute the neighbor of the cell along the face defined by the facenumber. .*/
t8_gloidx_t
t8_cmesh_neighbour_at_face (vtkSmartPointer < vtkUnstructuredGrid >
                            unstructuredGrid,
                            vtkSmartPointer < vtkIdList > pointIds,
                            t8_eclass_t eclass, int face_num,
                            vtkIdType cell_id)
{
  /* Points defining the face */
  vtkIdType          *face_points;
  /* cell_ids, will be filled by neighbor-computation */
  vtkSmartPointer < vtkIdList > cell_ids =
    vtkSmartPointer < vtkIdList >::New ();
  /*Get the eclass of the face */
  t8_eclass_t         face_eclass =
    (t8_eclass_t) t8_eclass_face_types[eclass][face_num];
  T8_ASSERT (face_eclass >= -1);
  /*Look up how many vertices the face has */
  vtkIdType           num_points = t8_eclass_num_vertices[face_eclass];
  face_points = T8_ALLOC (vtkIdType, num_points);
  /*Iterate over all corners of the face and store the pointid in face_points */
  for (int fp = 0; fp < num_points; fp++) {
    int                 face_corner =
      t8_vtk_cell_face_to_vertex_num[eclass][face_num][fp];
    T8_ASSERT (face_corner != -1);
    face_points[fp] = pointIds->GetId (face_corner);
  }
  /* Compute all cells touching the the given points (without cell with id cell_id
   * The computed id is unique, as only the face-neighbor is able to touch all points
   * of the face.*/
  unstructuredGrid->GetCellNeighbors (cell_id, num_points, face_points,
                                      cell_ids);
  T8_ASSERT (0 <= cell_ids->GetNumberOfIds ()
             && cell_ids->GetNumberOfIds () <= 1);
  T8_FREE (face_points);
  /*There is no neighbore -> element is at the boundary of the mesh */
  if (cell_ids->GetNumberOfIds () == 0) {
    return -1;
  }
  else {
    return cell_ids->GetId (0);
  }
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
  /*Get the file-extension to decide which reader to use */
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
  /*New Iterator to iterate over all cells in the grid */
  cell_it = unstructuredGrid->NewCellIterator ();
  max_cell_points = unstructuredGrid->GetMaxCellSize ();
  /*Allocate maximal possible points to avoid reallocation in the loop */
  vertices = T8_ALLOC (double, 3 * max_cell_points);
  /*Each cell in vtk will be a tree in the cmesh */
  t8_gloidx_t         tree_id = 0;
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal ();
       cell_it->GoToNextCell ()) {
    /*Set the t8_eclass of the cell */
    cell_type = t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    SC_CHECK_ABORTF (cell_type != T8_ECLASS_INVALID,
                     "vtk-cell-type %i not supported by t8code\n",
                     cell_it->GetCellType ());
    t8_debugf ("[D] cell is a %s\n", t8_eclass_to_string[cell_type]);
    t8_cmesh_set_tree_class (cmesh, tree_id, cell_type);

    /*Set the vertices of the tree */
    num_cell_points = cell_it->GetNumberOfPoints ();
    t8_debugf ("[D] cell has %i points\n", num_cell_points);
    points = cell_it->GetPoints ();
    for (int i = 0; i < num_cell_points; i++) {
      /*Every Point has 3 coords */
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

    /* We don't always need to know the connection of the cell. Therefore,
     * this computation is optional.*/
    if (compute_face_neigh) {
      /* WIP */
      unstructuredGrid->BuildLinks ();
      cell_id = cell_it->GetCellId ();
      vtkSmartPointer < vtkIdList > pointIds =
        vtkSmartPointer < vtkIdList >::New ();
      pointIds = cell_it->GetPointIds ();

      for (int face = 0; face < t8_eclass_num_faces[cell_type]; face++) {
        t8_gloidx_t         n =
          t8_cmesh_neighbour_at_face (unstructuredGrid, pointIds, cell_type,
                                      face, cell_id);
        t8_debugf ("[D] neigh: %li\n", n);
      }
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
