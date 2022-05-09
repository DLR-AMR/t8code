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
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCellIterator.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
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

/**
 * Given a list of corner-numbers and the eclass of the element with these corners
 * compute the face-number of the face described by the corner-numbers.
 * t8_face_vertex_to_tree_vertex stores the indices in increasing order, hence the incoming
 * ids have to be in order.
 * 
 * \param[in] facePointIds  The corner-numbers of a face of an element in ascending order
 * \param[in] eclass        The class of the element with face described by \a facePointIds
 * \param[in] num_points    The number of corners of the face
 * \return                  The number of the face in t8code-order.
 */
int
t8_cmesh_set_face_num (vtkSmartPointer < vtkIdList > facePointIds,
                       t8_eclass_t eclass, int num_points)
{
  int                 face_class, face_check;
  /*iterate over all faces of the class */
  for (int face_iter = 0; face_iter < t8_eclass_num_faces[eclass];
       face_iter++) {
    face_check = 0;
    face_class = t8_eclass_face_types[eclass][face_iter];
    T8_ASSERT (face_class != -1);
    /*iterate over all tree-vertices of that class */
    /*If the number of points does not match the number of points of the current face, we can skip the check */
    if (num_points == t8_eclass_num_vertices[face_class]) {
      for (int i = 0; i < num_points; i++) {
        /*check, if all indices of that face coincide with the ids given by vertex_ids */
        if (facePointIds->GetId (i) !=
            t8_face_vertex_to_tree_vertex[eclass][face_iter][i]) {
          face_check = 1;
        }
      }
      /*If all vertices coincide, the current face is the correct face_number */
      if (face_check == 0) {
        t8_debugf ("[D] neigh_face: %i\n", face_iter);
        return face_iter;
      }
    }
  }
  /*Error, no face_number found */
  SC_ABORT ("No matching face found");
  return -1;
}

/** Given the points of a vtk-cell and a set of corner-numbers of a face in
 *  t8code order find the vtk-ids of the face and set the face-number
 * 
 * \param[in]       pointIds        The vtk-ids of the points of the Cell     
 * \param[in, out]  facePointIds    Fill with the t8code-corner-numbers of the face given by \a face_points
 * \param[in, out]  face            A face, set the the face_number
 * \param[in]       cell_type       The eclass of the vtk-cell
 * \param[in]       num_points      The number of points of the face
 * \param[in]       face_points     The vtk-ids of a face of the Cell described by pointIds
 */
void
vtk_face_to_t8_face (vtkSmartPointer < vtkIdList > pointIds,
                     vtkSmartPointer < vtkIdList > facePointIds,
                     t8_msh_file_face_t * face,
                     t8_eclass_t cell_type, int num_face_points,
                     vtkIdType * face_points)
{
  int                 vtk_corner, t8_corner;
  /* Iterate over all points on that face */
  for (int fp = 0; fp < num_face_points; fp++) {
    /* Get the vtk-corner number of the face_point */
    vtk_corner = pointIds->FindIdLocation (face_points[fp]);
    /* From vtk-corner to t8code-corner */
    t8_corner = t8_eclass_vtk_corner_number[cell_type][vtk_corner];
    /*Fill the face */
    facePointIds->SetId (fp, t8_corner);
  }
  /*t8_cmesh_set_face_num needs the ids in ascending order */
  facePointIds->Sort ();
  /*Set the face-number */
  face->face_number =
    t8_cmesh_set_face_num (facePointIds, cell_type, num_face_points);
  T8_ASSERT (face->face_number >= 0);
}

/**
 * Set the vertice-arrays of \a face_a and its neighbourface \a face_b. 
 * 
 * \param[in, out]  face_a            A face of an element, the vertices will be filled with the vtk-ids of the face , the vertices will be filled with the vtk-ids of the face
 * \param[in, out]  face_b            The face of a neighbouring element, touching \a face_a , the vertices will be filled with the vtk-ids of the face , the vertices will be filled with the vtk-ids of the face
 * \param[in]       face_points       The vtk-ids describing the face
 * \param[in]       neighFacePointIds The t8code corners of the face of the neighbor
 * \param[in]       neighPointIds     The vtk-ids describing the neighbor
 * \param[in]       num_face_points   The number of points of the face
 * \param[in]       neigh_type        The type of the neighbor
 */
void
t8_set_face_and_neigh_face (t8_msh_file_face_t * face_a,
                            t8_msh_file_face_t * face_b,
                            vtkIdType * face_points,
                            vtkSmartPointer < vtkIdList > neighFacePointIds,
                            vtkSmartPointer < vtkIdList > neighPointIds,
                            int num_face_points, t8_eclass_t neigh_type)
{
  int                 tree_corner, face_corner;
  /* Iterate over all points of the face */
  for (int fp = 0; fp < num_face_points; fp++) {
    /* get the t8code-id of the vertex */
    tree_corner =
      t8_face_vertex_to_tree_vertex[neigh_type][face_b->face_number]
      [neighFacePointIds->GetId (fp)];
    /*transform to local vtk-id */
    face_corner = t8_eclass_vtk_corner_number[neigh_type][tree_corner];
    face_a->vertices[fp] = face_points[fp];
    /*set the point-id according to the local vtk-id */
    face_b->vertices[fp] = neighPointIds->GetId (face_corner);
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
  int                 cell_type;
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkCellIterator    *cell_it;
  vtkCellData        *cellData;
  int                 num_data_arrays = 0;
  size_t             *data_size;
  double            **tuples;   /*vtk stores data as doubles */
  char               *tmp, *extension;
  int                 max_dim = -1;     /*max dimenstion of the cells for geometry */
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
  T8_ASSERT (strcmp (extension, ""));

  /*Read the file */
  if (strcmp (extension, "vtu") == 0) {
    vtkSmartPointer < vtkXMLUnstructuredGridReader > reader =
      vtkSmartPointer < vtkXMLUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    unstructuredGrid = reader->GetOutput ();
  }
  else if (strcmp (extension, "vtk") == 0) {
    vtkSmartPointer < vtkUnstructuredGridReader > reader =
      vtkSmartPointer < vtkUnstructuredGridReader >::New ();
    reader->SetFileName (filename);
    reader->Update ();
    unstructuredGrid = reader->GetOutput ();
  }
  else {
    t8_global_errorf ("Please use .vtk or .vtu file\n");
  }
  T8_FREE (tmp);
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

T8_EXTERN_C_END ();
