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

#include <t8_cmesh_vtk_writer.h>        /* To set tree vertices */

#include "t8_vtk_reader.hxx"
#include "t8_vtk_unstructured.hxx"
#include "t8_vtk_polydata.hxx"
#include "t8_vtk_types.h"
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

#if T8_WITH_VTK
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
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

/**
 * If the vertices of a tree describe a negative \param, 
 * permute the tree vertices. 
 * 
 * \param[in, out] tree_vertices The vertices of a tree
 * \param[in] eclass             The eclass of the tree.
 */
void
t8_cmesh_correct_volume (double *tree_vertices, t8_eclass_t eclass)
{
  /* The \param described is negative. We need to change vertices.
   * For tets we switch 0 and 3.
   * For prisms we switch 0 and 3, 1 and 4, 2 and 5.
   * For hexahedra we switch 0 and 4, 1 and 5, 2 and 6, 3 and 7.
   * For pyramids we switch 0 and 4 */
  double              temp;
  int                 num_switches = 0;
  int                 switch_indices[4] = { 0 };
  int                 iswitch;
  T8_ASSERT (t8_eclass_to_dimension[eclass] == 3);
  t8_debugf ("Correcting negative volume.\n");
  switch (eclass) {
  case T8_ECLASS_TET:
    /* We switch vertex 0 and vertex 3 */
    num_switches = 1;
    switch_indices[0] = 3;
    break;
  case T8_ECLASS_PRISM:
    num_switches = 3;
    switch_indices[0] = 3;
    switch_indices[1] = 4;
    switch_indices[2] = 5;
    break;
  case T8_ECLASS_HEX:
    num_switches = 4;
    switch_indices[0] = 4;
    switch_indices[1] = 5;
    switch_indices[2] = 6;
    switch_indices[3] = 7;
    break;
  case T8_ECLASS_PYRAMID:
    num_switches = 1;
    switch_indices[0] = 4;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  for (iswitch = 0; iswitch < num_switches; ++iswitch) {
    /* We switch vertex 0 + iswitch and vertex switch_indices[iswitch] */
    for (int i = 0; i < 3; i++) {
      temp = tree_vertices[3 * iswitch + i];
      tree_vertices[3 * iswitch + i] =
        tree_vertices[3 * switch_indices[iswitch] + i];
      tree_vertices[3 * switch_indices[iswitch] + i] = temp;
    }
  }
  T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume
             (eclass, tree_vertices, t8_eclass_num_vertices[eclass]));
}

#if T8_WITH_VTK

vtk_read_success_t
t8_file_to_vtkGrid (const char *filename,
                    vtkSmartPointer < vtkDataSet > vtkGrid,
                    const int partition, const int main_proc,
                    sc_MPI_Comm comm, const vtk_file_type_t vtk_file_type)
{
  vtk_read_success_t  main_proc_read_successful = read_failure;
  int                 mpirank;
  int                 mpisize;
  int                 mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (filename != NULL);
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);
  /* Read the file and set the pointer to the vtkGrid */
  if (!partition || mpirank == main_proc) {
    switch (vtk_file_type) {
    case VTK_UNSTRUCTURED_FILE:
      main_proc_read_successful = t8_read_unstructured (filename, vtkGrid);
      break;
    case VTK_POLYDATA_FILE:
      main_proc_read_successful = t8_read_poly (filename, vtkGrid);
      break;
    default:
      vtkGrid = NULL;
      t8_errorf ("Filetype not supported.\n");
      break;
    }
    if (partition) {
      /* Communicate the success/failure of the reading process. */
      sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc,
                    comm);
    }
  }
  if (partition) {
    sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
  }
  return main_proc_read_successful;
}

/**
 * Get the dimension of a vtkDataSet
 * 
 * \param[in] vtkGrid The vtkDataSet
 * \return The dimension of \a vtkGrid. 
 */
int
t8_get_dimension (vtkSmartPointer < vtkDataSet > vtkGrid)
{
  /* This array contains the type of each cell */
  vtkSmartPointer < vtkCellTypes > cell_type_of_each_cell =
    vtkSmartPointer < vtkCellTypes >::New ();
  vtkGrid->GetCellTypes (cell_type_of_each_cell);
  vtkSmartPointer < vtkUnsignedCharArray > cell_types =
    cell_type_of_each_cell->GetCellTypesArray ();
  int                 max_cell_type = -1;
  const vtkIdType     num_types = cell_types->GetNumberOfTuples ();
  /* Iterate over all types and compare the dimension. They are stored
   * in ascending order. */
  for (vtkIdType cell_types_it = 0; cell_types_it < num_types;
       cell_types_it++) {
    const vtkIdType     type = cell_types->GetValue (cell_types_it);
    if (type > max_cell_type) {
      max_cell_type = type;
    }
  }
  T8_ASSERT (0 <= max_cell_type && max_cell_type < 82);
  const int           ieclass = t8_cmesh_vtk_type_to_t8_type[max_cell_type];
  T8_ASSERT (ieclass != T8_ECLASS_INVALID);

  return t8_eclass_to_dimension[ieclass];
}

/**
 * Iterate over all cells of a vtkDataset and construct a cmesh representing
 * The vtkGrid. Each cell in the vtkDataSet becomes a tree in the cmesh. This 
 * function construct a cmesh on a single process. 
 * 
 * \param[in] vtkGrid The vtkGrid that gets tranlated
 * \param[in, out] cmesh   An empty cmesh that is filled with the data. 
 * \param[in] comm        A communicator. 
 * \return  The number of elements that have been read by the process.  
 */

t8_gloidx_t
t8_vtk_iterate_cells (vtkSmartPointer < vtkDataSet > vtkGrid,
                      t8_cmesh_t cmesh, sc_MPI_Comm comm)
{

  double             *vertices;
  double            **tuples;
  size_t             *data_size;
  t8_gloidx_t         tree_id = 0;
  int                 max_dim = -1;

  vtkCellIterator    *cell_it;
  vtkSmartPointer < vtkPoints > points;
  vtkSmartPointer < vtkCellData > cell_data = vtkGrid->GetCellData ();
  const int           max_cell_points = vtkGrid->GetMaxCellSize ();

  T8_ASSERT (max_cell_points >= 0);
  vertices = T8_ALLOC (double, 3 * max_cell_points);
  /* Get cell iterator */
  cell_it = vtkGrid->NewCellIterator ();
  /* get the number of data-arrays per cell */
  const int           num_data_arrays = cell_data->GetNumberOfArrays ();
  T8_ASSERT (num_data_arrays >= 0);

  t8_debugf ("[D] read %i data-arrays\n", num_data_arrays);
  /* Prepare attributes */
  if (num_data_arrays > 0) {
    size_t              tuple_size;
    tuples = T8_ALLOC (double *, num_data_arrays);
    data_size = T8_ALLOC (size_t, num_data_arrays);
    for (int idata = 0; idata < num_data_arrays; idata++) {
      vtkDataArray       *data = cell_data->GetArray (idata);
      tuple_size = data->GetNumberOfComponents ();
      data_size[idata] = sizeof (double) * tuple_size;
      t8_debugf ("[D] data_size[%i] = %li, tuple_size %li\n", idata,
                 data_size[idata], tuple_size);
      /* Allocate memory for a tuple in array i */
      tuples[idata] = T8_ALLOC (double, tuple_size);
    }
  }

  /* Iterate over all cells */
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal ();
       cell_it->GoToNextCell ()) {

    /* Set the t8_eclass of the cell */
    const t8_eclass_t   cell_type =
      t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    SC_CHECK_ABORTF (t8_eclass_is_valid (cell_type),
                     "vtk-cell-type %i not supported by t8code\n", cell_type);
    t8_cmesh_set_tree_class (cmesh, tree_id, cell_type);
    /* Get the points of the cell */
    const int           num_points = cell_it->GetNumberOfPoints ();
    T8_ASSERT (num_points > 0);
    points = cell_it->GetPoints ();

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      points->GetPoint (t8_element_shape_vtk_corner_number
                        (cell_type, ipoint), &vertices[3 * ipoint]);

    }
    /* The order of the vertices in vtk might give a tree with negative \param */
    if (t8_cmesh_tree_vertices_negative_volume
        (cell_type, vertices, num_points)) {
      t8_cmesh_correct_volume (vertices, cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_points);

    /* TODO: Avoid magic numbers in the attribute setting. */
    /* Get and set the data of each cell */
    for (int dtype = 0; dtype < num_data_arrays; dtype++) {
      const t8_gloidx_t   cell_id = cell_it->GetCellId ();
      vtkDataArray       *data = cell_data->GetArray (dtype);
      data->GetTuple (cell_id, tuples[dtype]);
      t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), dtype + 1,
                              tuples[dtype], data_size[dtype], 0);
    }
    /* Check geometry-dimension */
    if (max_dim < cell_it->GetCellDimension ()) {
      max_dim = cell_it->GetCellDimension ();
    }
    tree_id++;
  }
  t8_debugf ("[D] read %li trees\n", tree_id);

  /* Clean-up */
  cell_it->Delete ();
  if (num_data_arrays > 0) {
    T8_FREE (data_size);
    for (int idata = num_data_arrays - 1; idata >= 0; idata--) {
      T8_FREE (tuples[idata]);
    }
    T8_FREE (tuples);
  }
  T8_FREE (vertices);
  return tree_id;
}

t8_cmesh_t
t8_vtkGrid_to_cmesh (vtkSmartPointer < vtkDataSet > vtkGrid,
                     const int partition, const int main_proc,
                     sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  int                 mpisize;
  int                 mpirank;
  int                 mpiret;
  /* Get the size of the communicator and the rank of the process. */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Ensure that the main-proc is a valid proc. */
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);

  /* Already declared here, because we might use them during communication */
  t8_gloidx_t         num_trees;
  int                 dim;

  t8_cmesh_init (&cmesh);
  if (!partition || mpirank == main_proc) {
    num_trees = t8_vtk_iterate_cells (vtkGrid, cmesh, comm);
    dim = t8_get_dimension (vtkGrid);
    t8_cmesh_set_dimension (cmesh, dim);
    t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);
    t8_cmesh_register_geometry (cmesh, linear_geom);
  }
  if (partition) {
    t8_gloidx_t         first_tree;
    t8_gloidx_t         last_tree;
    if (mpirank == main_proc) {
      first_tree = 0;
      last_tree = num_trees - 1;
    }
    /* Communicate the dimension to all processes */
    sc_MPI_Bcast (&dim, 1, sc_MPI_INT, main_proc, comm);
    t8_debugf ("[D] dim: %i\n", dim);
    /* Communicate the number of trees to all processes. 
     * TODO: This probably crashes when a vtkGrid is distributed in many 
     * files. */
    sc_MPI_Bcast (&num_trees, 1, T8_MPI_GLOIDX, main_proc, comm);
    t8_cmesh_set_dimension (cmesh, dim);
    /* Build the partition. */
    if (mpirank < main_proc) {
      first_tree = 0;
      last_tree = -1;
      t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);
      t8_cmesh_register_geometry (cmesh, linear_geom);
    }
    else if (mpirank > main_proc) {
      first_tree = num_trees;
      last_tree = num_trees - 1;
      t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);
      t8_cmesh_register_geometry (cmesh, linear_geom);
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  return cmesh;
}

vtkSmartPointer < vtkDataSet >
t8_vtk_reader (const char *filename, const int partition,
               const int main_proc, sc_MPI_Comm comm,
               const vtk_file_type_t vtk_file_type)
{
  int                 mpisize;
  int                 mpirank;
  int                 mpiret;
  /* Get the size of the communicator and the rank of the process. */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Ensure that the main-proc is a valid proc. */
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);
  T8_ASSERT (filename != NULL);
  vtk_read_success_t  main_proc_read_successful = read_failure;

  vtkSmartPointer < vtkDataSet > vtkGrid;
  switch (vtk_file_type) {
  case VTK_UNSTRUCTURED_FILE:
    vtkGrid = vtkSmartPointer < vtkUnstructuredGrid >::New ();
    break;
  case VTK_POLYDATA_FILE:
    vtkGrid = vtkSmartPointer < vtkPolyData >::New ();
    break;
  default:
    t8_errorf ("Filetype is not supported.\n");
    break;
  }
  T8_ASSERT (partition == 0 || (main_proc >= 0 && main_proc < mpisize));
  /* Read the file and set the pointer. */
  main_proc_read_successful =
    t8_file_to_vtkGrid (filename, vtkGrid, partition, main_proc, comm,
                        vtk_file_type);

  if (!main_proc_read_successful) {
    t8_global_errorf
      ("Main process (Rank %i) did not read the file successfully.\n",
       main_proc);
    return NULL;
  }
  else {
    return vtkGrid;
  }
}

#endif /* T8_WITH_VTK */

t8_cmesh_t
t8_vtk_reader_cmesh (const char *filename, const int partition,
                     const int main_proc, sc_MPI_Comm comm,
                     const vtk_file_type_t vtk_file_type)
{
#if T8_WITH_VTK
  vtkSmartPointer < vtkDataSet > vtkGrid =
    t8_vtk_reader (filename, partition, main_proc, comm, vtk_file_type);
  if (vtkGrid != NULL) {
    t8_cmesh_t          cmesh =
      t8_vtkGrid_to_cmesh (vtkGrid, partition, main_proc, comm);
    T8_ASSERT (cmesh != NULL);
    return cmesh;
  }
  else {
    t8_global_errorf ("Error translating file %s\n", filename);
    return NULL;
  }
#else
  /* Return NULL if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
#endif
  return NULL;
}

T8_EXTERN_C_END ();
