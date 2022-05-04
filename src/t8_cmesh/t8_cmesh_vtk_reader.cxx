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
                            vtkIdType cell_id, vtkIdType * face_points)
{

  /* cell_ids, will be filled by neighbor-computation */
  vtkSmartPointer < vtkIdList > cell_ids =
    vtkSmartPointer < vtkIdList >::New ();
  /*Get the eclass of the face */
  t8_eclass_t         face_eclass =
    (t8_eclass_t) t8_eclass_face_types[eclass][face_num];
  int                 face_corner, vtk_corner;
  T8_ASSERT (face_eclass >= -1);
  /*Look up how many vertices the face has */
  vtkIdType           num_points = t8_eclass_num_vertices[face_eclass];
  /*Iterate over all corners of the face and store the pointid in face_points */
  for (int fp = 0; fp < num_points; fp++) {
    /*tree-corners of the face (t8code-numeration of the faces) in t8code-order */
    face_corner = t8_face_vertex_to_tree_vertex[eclass][face_num][fp];
    T8_ASSERT (face_corner != -1);
    vtk_corner = t8_eclass_vtk_corner_number[eclass][face_corner];
    T8_ASSERT (vtk_corner != -1);
    face_points[fp] = pointIds->GetId (face_corner);
  }
  /* Compute all cells touching the the given points (without cell with id cell_id
   * The computed id is unique, as only the face-neighbor is able to touch all points
   * of the face.*/
  unstructuredGrid->GetCellNeighbors (cell_id, num_points, face_points,
                                      cell_ids);
  T8_ASSERT (0 <= cell_ids->GetNumberOfIds ()
             && cell_ids->GetNumberOfIds () <= 1);
  /*There is no neighbor -> element is at the boundary of the mesh */
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
  int                 cell_type, neigh_type, face_type;
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkCellIterator    *cell_it;
  char               *tmp, *extension;
  int                 max_dim = 0;      /*max dimenstion of the cells for geometry */
  int                 num_face_points, num_cell_points, max_cell_points;
  t8_gloidx_t         cell_id, neigh_id;
  vtkSmartPointer < vtkPoints > points =
    vtkSmartPointer < vtkPoints >::New ();
  double             *vertices;
  t8_msh_file_face_t  face_a, face_b;
  /* Points defining the face */
  vtkIdType          *face_points;
  face_points = T8_ALLOC (vtkIdType, T8_ECLASS_MAX_CORNERS_2D);
  /*We only allow 2D-Faces */
  face_a.vertices = T8_ALLOC (long, T8_ECLASS_MAX_CORNERS_2D);
  face_b.vertices = T8_ALLOC (long, T8_ECLASS_MAX_CORNERS_2D);
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
    t8_cmesh_set_tree_class (cmesh, tree_id, (t8_eclass_t) cell_type);

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
        ((t8_eclass_t) cell_type, vertices, num_cell_points)) {
      t8_cmesh_correct_volume (vertices, (t8_eclass_t) cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_cell_points);

    /* We don't always need to know the connection of the cell. Therefore,
     * this computation is optional.*/
    if (compute_face_neigh) {
      /* WIP */
      vtkSmartPointer < vtkGenericCell > neigh_cell =
        vtkSmartPointer < vtkGenericCell >::New ();
      unstructuredGrid->BuildLinks ();
      cell_id = cell_it->GetCellId ();
      vtkSmartPointer < vtkIdList > pointIds =
        vtkSmartPointer < vtkIdList >::New ();
      vtkSmartPointer < vtkIdList > neighIds =
        vtkSmartPointer < vtkIdList >::New ();
      vtkSmartPointer < vtkIdList > cellpointIds =
        vtkSmartPointer < vtkIdList >::New ();
      vtkSmartPointer < vtkIdList > neigh_point =
        vtkSmartPointer < vtkIdList >::New ();
      pointIds = cell_it->GetPointIds ();
      face_a.ltree_id = tree_id;
      for (int face = 0; face < t8_eclass_num_faces[cell_type]; face++) {
        face_type = t8_eclass_face_types[cell_type][face];
        num_face_points = t8_eclass_num_vertices[face_type];
        //face_a.face_number = face;
        face_a.num_vertices = num_face_points;
        /* vtk_point_ids of the face in face_points */
        neigh_id =
          t8_cmesh_neighbour_at_face (unstructuredGrid, pointIds,
                                      (t8_eclass_t) cell_type, face, cell_id,
                                      face_points);
        if (neigh_id != -1) {
          int                 orientation, vtk_corner, face_corner,
            tree_corner;
          /* Be carefull here with parallel reader or partitioned cmesh
           * This will probably not work for partitioned cmeshes.*/
          face_b.ltree_id = neigh_id;
          face_b.num_vertices = num_face_points;
          t8_debugf ("[D] neigh: %li of current tree: %i\n", neigh_id,
                     cell_id);;
          unstructuredGrid->GetCell (neigh_id, neigh_cell);
          neigh_point = neigh_cell->GetPointIds ();
          neigh_type =
            t8_cmesh_vtk_type_to_t8_type[neigh_cell->GetCellType ()];
          t8_debugf ("[D] neigh_type is: %i %s\n", neigh_type,
                     t8_eclass_to_string[neigh_type]);
          neighIds->SetNumberOfIds (num_face_points);
          cellpointIds->SetNumberOfIds (num_face_points);
          for (int fp = 0; fp < num_face_points; fp++) {
            vtk_corner = pointIds->FindIdLocation (face_points[fp]);
            face_corner = t8_eclass_vtk_corner_number[cell_type][vtk_corner];
            cellpointIds->SetId (fp, face_corner);
            //face_a.vertices[fp] = face_points[fp];
            /*tree-corners of the face (t8code-numeration of the faces) in t8code-order */
            T8_ASSERT (face_corner != -1);
            /* vtk_point_ids in t8code-order of the face */
            face_a.vertices[fp] = face_corner;
            vtk_corner = neigh_point->FindIdLocation (face_points[fp]);
            face_corner = t8_eclass_vtk_corner_number[neigh_type][vtk_corner];
            face_b.vertices[fp] = face_corner;
            neighIds->SetId (fp, face_corner);
          }
          neighIds->Sort ();
          cellpointIds->Sort ();
          for (int fp = 0; fp < num_face_points; fp++) {
            face_b.vertices[fp] = neighIds->GetId (fp);
            face_a.vertices[fp] = cellpointIds->GetId (fp);
          }
          t8_cmesh_set_face_num (&face_b, (t8_eclass_t) neigh_type);
          t8_cmesh_set_face_num (&face_a, (t8_eclass_t) cell_type);
          for (int fp = 0; fp < num_face_points; fp++) {
            tree_corner =
              t8_face_vertex_to_tree_vertex[neigh_type][face][neighIds->GetId
                                                              (fp)];
            face_corner =
              t8_eclass_vtk_corner_number[neigh_type][tree_corner];
            T8_ASSERT (face_corner != -1);
            face_a.vertices[fp] = face_points[fp];
            face_b.vertices[fp] = neigh_point->GetId (face_corner);
            t8_debugf ("[D] face_b.v(%i): %i\n", fp, face_b.vertices[fp]);
          }
          orientation =
            t8_msh_file_face_orientation (&face_a, &face_b,
                                          (t8_eclass_t) cell_type,
                                          (t8_eclass_t) neigh_type);
          t8_debugf
            ("[D] set face of tree %i, along face: %i to neigh: %i with face %i and orientation: %i\n",
             tree_id, face_a.face_number, neigh_id, face_b.face_number,
             orientation);
          t8_cmesh_set_join (cmesh, tree_id, neigh_id, face_a.face_number,
                             face_b.face_number, orientation);
        }
        else {
          t8_debugf ("[D] tree %i has no neighbor along vtk-face: %i\n",
                     tree_id, face);
        }

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
  T8_FREE (face_a.vertices);
  T8_FREE (face_b.vertices);
  T8_FREE (vertices);
  T8_FREE (face_points);
  return cmesh;
}

T8_EXTERN_C_END ();
