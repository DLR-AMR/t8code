/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8_vtk/t8_with_vtk/t8_vtk_reader.hxx>
#include <t8_vtk/t8_with_vtk/t8_vtk_reader.hxx>
#include <t8_vtk/t8_vtk_types.h>
#include <t8_cmesh.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>

#include <t8_vtk/t8_with_vtk/t8_vtk_polydata.hxx>
#include <t8_vtk/t8_with_vtk/t8_vtk_parallel.hxx>
#include <t8_vtk/t8_with_vtk/t8_vtk_unstructured.hxx>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellTypes.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>

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
  double temp;
  int num_switches = 0;
  int switch_indices[4] = { 0 };
  int iswitch;
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
      tree_vertices[3 * iswitch + i] = tree_vertices[3 * switch_indices[iswitch] + i];
      tree_vertices[3 * switch_indices[iswitch] + i] = temp;
    }
  }
  T8_ASSERT (!t8_cmesh_tree_vertices_negative_volume (eclass, tree_vertices, t8_eclass_num_vertices[eclass]));
}

vtk_read_success_t
t8_file_to_vtkGrid (const char *filename, vtkSmartPointer<vtkDataSet> vtkGrid, const int partition, const int main_proc,
                    sc_MPI_Comm comm, const vtk_file_type_t vtk_file_type)
{
  vtk_read_success_t main_proc_read_successful = read_failure;
  int mpirank;
  int mpisize;
  int mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  T8_ASSERT (filename != NULL);
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);
  /* Read the file and set the pointer to the vtkGrid
   * We read the file if:
   * - We do not use a partitioned read, every process reads the vtk-file, or if
   * - We use a partitioned read for a non-parallel filetype and the main-process reaches the code-block, or if
   * - We use a parallel file-type and use a partitioned read, every proc reads its chunk of the files. 
   */
  if (!partition || mpirank == main_proc || vtk_file_type & VTK_PARALLEL_FILE) {
    switch (vtk_file_type) {
    case VTK_UNSTRUCTURED_FILE:
      main_proc_read_successful = t8_read_unstructured (filename, vtkGrid);
      break;
    case VTK_POLYDATA_FILE:
      main_proc_read_successful = t8_read_polyData (filename, vtkGrid);
      break;
    case VTK_PARALLEL_UNSTRUCTURED_FILE:
      if (!partition) {
        main_proc_read_successful = t8_read_unstructured (filename, vtkGrid);
      }
      else {
        main_proc_read_successful = t8_read_parallel_unstructured (filename, vtkGrid, comm);
        break;
      }
      break;
    case VTK_PARALLEL_POLYDATA_FILE:
      if (!partition) {
        main_proc_read_successful = t8_read_polyData (filename, vtkGrid);
      }
      else {
        main_proc_read_successful = t8_read_parallel_polyData (filename, vtkGrid, comm);
        break;
      }
      break;
    default:
      vtkGrid = NULL;
      t8_errorf ("Filetype not supported.\n");
      break;
    }
    if (partition && vtk_file_type & VTK_SERIAL_FILE) {
      /* Communicate the success/failure of the reading process. */
      sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
      return main_proc_read_successful;
    }
  }
  if (partition) {
    if (vtk_file_type & VTK_SERIAL_FILE) {
      sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
    }
    else {
      int recv_buf;
      sc_MPI_Allreduce (&main_proc_read_successful, &recv_buf, 1, sc_MPI_INT, sc_MPI_LOR, comm);
      main_proc_read_successful = (vtk_read_success_t) recv_buf;
    }
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
t8_get_dimension (vtkSmartPointer<vtkDataSet> vtkGrid)
{
  /* This array contains the type of each cell */
  vtkSmartPointer<vtkCellTypes> cell_type_of_each_cell = vtkSmartPointer<vtkCellTypes>::New ();
  vtkGrid->GetCellTypes (cell_type_of_each_cell);
  vtkSmartPointer<vtkUnsignedCharArray> cell_types = cell_type_of_each_cell->GetCellTypesArray ();
  int max_cell_type = -1;
  const vtkIdType num_types = cell_types->GetNumberOfTuples ();
  /* Iterate over all types and compare the dimension. They are stored
   * in ascending order. */
  for (vtkIdType cell_types_it = 0; cell_types_it < num_types; cell_types_it++) {
    const vtkIdType type = cell_types->GetValue (cell_types_it);
    if (type > max_cell_type) {
      max_cell_type = type;
    }
  }
  T8_ASSERT (0 <= max_cell_type && max_cell_type < 82);
  const int ieclass = t8_cmesh_vtk_type_to_t8_type[max_cell_type];
  T8_ASSERT (ieclass != T8_ECLASS_INVALID);

  return t8_eclass_to_dimension[ieclass];
}

/**
 * Iterate over all cells of a vtkDataset and construct a cmesh representing
 * the vtkGrid. Each cell in the vtkDataSet becomes a tree in the cmesh. This 
 * function constructs a cmesh on a single process. 
 * 
 * \param[in] vtkGrid       The vtkGrid that gets translated
 * \param[in, out] cmesh    An empty cmesh that is filled with the data. 
 * \param[in] first_tree    The global id of the first tree. Will be the global id of the first tree on this proc. 
 * \param[in] comm        A communicator. 
 * \return  The number of elements that have been read by the process.
 */

static void
t8_vtk_iterate_cells (vtkSmartPointer<vtkDataSet> vtkGrid, t8_cmesh_t cmesh, const t8_gloidx_t first_tree,
                      [[maybe_unused]] sc_MPI_Comm comm)
{
  double **tuples = NULL;
  size_t *data_size = NULL;
  t8_gloidx_t tree_id = first_tree;

  vtkCellIterator *cell_it;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellData> cell_data = vtkGrid->GetCellData ();
  const int max_cell_points = vtkGrid->GetMaxCellSize ();

  T8_ASSERT (max_cell_points >= 0);
  double *vertices = T8_ALLOC (double, 3 * max_cell_points);
  /* Get cell iterator */
  cell_it = vtkGrid->NewCellIterator ();
  /* get the number of data-arrays per cell */
  const int num_data_arrays = cell_data->GetNumberOfArrays ();
  T8_ASSERT (num_data_arrays >= 0);

  /* Prepare attributes */
  if (num_data_arrays > 0) {
    size_t tuple_size;
    tuples = T8_ALLOC (double *, num_data_arrays);
    data_size = T8_ALLOC (size_t, num_data_arrays);
    for (int idata = 0; idata < num_data_arrays; idata++) {
      vtkDataArray *data = cell_data->GetArray (idata);
      tuple_size = data->GetNumberOfComponents ();
      data_size[idata] = sizeof (double) * tuple_size;
      /* Allocate memory for a tuple in array i */
      tuples[idata] = T8_ALLOC (double, tuple_size);
    }
  }

  /* Iterate over all cells */
  for (cell_it->InitTraversal (); !cell_it->IsDoneWithTraversal (); cell_it->GoToNextCell ()) {

    /* Set the t8_eclass of the cell */
    const t8_eclass_t cell_type = t8_cmesh_vtk_type_to_t8_type[cell_it->GetCellType ()];
    SC_CHECK_ABORTF (t8_eclass_is_valid (cell_type), "vtk-cell-type %i not supported by t8code\n", cell_type);
    t8_cmesh_set_tree_class (cmesh, tree_id, cell_type);
    /* Get the points of the cell */
    const int num_points = cell_it->GetNumberOfPoints ();
    T8_ASSERT (num_points > 0);
    points = cell_it->GetPoints ();

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      points->GetPoint (t8_element_shape_vtk_corner_number (cell_type, ipoint), &vertices[3 * ipoint]);
    }
    /* The order of the vertices in vtk might give a tree with negative \param */
    if (t8_cmesh_tree_vertices_negative_volume (cell_type, vertices, num_points)) {
      t8_cmesh_correct_volume (vertices, cell_type);
    }
    t8_cmesh_set_tree_vertices (cmesh, tree_id, vertices, num_points);

    /* TODO: Avoid magic numbers in the attribute setting. */
    /* Get and set the data of each cell */
    for (int dtype = 0; dtype < num_data_arrays; dtype++) {
      const t8_gloidx_t cell_id = cell_it->GetCellId ();
      vtkDataArray *data = cell_data->GetArray (dtype);
      data->GetTuple (cell_id, tuples[dtype]);
      t8_cmesh_set_attribute (cmesh, tree_id, t8_get_package_id (), dtype + 1, tuples[dtype], data_size[dtype], 0);
    }
    tree_id++;
  }
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
  return;
}

/**
 * Set the partition for cmesh coming from a distributed vtkGrid (like pvtu)
 * 
 * \param[in, out] cmesh On input a cmesh, on output a cmesh with a partition according to the number of trees read on each proc
 * \param[in] mpirank The mpirank of this proc
 * \param[in] mpisize The size of the mpicommunicator
 * \param[in] num_trees The number of trees read on this proc. Must be >= 0
 * \param[in] dim     The dimension of the cells
 * \param[in] comm    The mpi-communicator to use. 
 * \return            the global id of the first tree on this proc. 
 */
static t8_gloidx_t
t8_vtk_partition (t8_cmesh_t cmesh, const int mpirank, const int mpisize, t8_gloidx_t num_trees,
                  [[maybe_unused]] int dim, sc_MPI_Comm comm)
{
  t8_gloidx_t first_tree = 0;
  t8_gloidx_t last_tree = 1;
  /* Compute the global id of the first tree on each proc. */
  t8_shmem_init (comm);
  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
  t8_shmem_array_t offsets = NULL;
  t8_shmem_array_init (&offsets, sizeof (t8_gloidx_t), mpisize + 1, comm);

  t8_shmem_array_prefix ((void *) &num_trees, offsets, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);

  first_tree = t8_shmem_array_get_gloidx (offsets, mpirank);
  /* Set the partition of the cmesh. */
  if (num_trees == 0) {
    last_tree = first_tree - 1;
  }
  else {
    last_tree = first_tree + num_trees - 1;
  }
  const int set_face_knowledge = 3; /* Expect face connection of local and ghost trees. */
  t8_cmesh_set_partition_range (cmesh, set_face_knowledge, first_tree, last_tree);
  t8_shmem_array_destroy (&offsets);
  return first_tree;
}

t8_cmesh_t
t8_vtkGrid_to_cmesh (vtkSmartPointer<vtkDataSet> vtkGrid, const int partition, const int main_proc,
                     const int distributed_grid, sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  int mpisize;
  int mpirank;
  int mpiret;
  t8_cmesh_init (&cmesh);
  /* Get the size of the communicator and the rank of the process. */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Ensure that the main-proc is a valid proc. */
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);

  /* Already declared here, because we might use them during communication */
  t8_gloidx_t num_trees = 0;
  if (!partition || mpirank == main_proc || distributed_grid) {
    num_trees = vtkGrid->GetNumberOfCells ();
  }

  /* Set the dimension on all procs (even empty procs). */
  int dim = num_trees > 0 ? t8_get_dimension (vtkGrid) : 0;
  int dim_buf = dim;
  mpiret = sc_MPI_Allreduce ((void *) &dim, &dim_buf, 1, sc_MPI_INT, sc_MPI_MAX, comm);
  SC_CHECK_MPI (mpiret);
  t8_cmesh_set_dimension (cmesh, dim_buf);

  /* Set the geometry. */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* Global-id of the first local tree */
  t8_gloidx_t first_tree = 0;

  /* Set the partition first, so we know the global id of the first tree on all procs. */
  if (partition) {
    first_tree = t8_vtk_partition (cmesh, mpirank, mpisize, num_trees, dim, comm);
  }

  /* Translation of vtkGrid to cmesh 
   * We translate the file if:
   * - We do not use a partitioned read, every process has read the vtk-file and shall now translate it, or if
   * - We use a partitioned read for a non-parallel filetype and the main-process reaches the code-block, or if
   * - We use a parallel file-type and use a partitioned read, every proc translates its chunk of the grid. 
   */
  if (!partition || mpirank == main_proc || distributed_grid) {
    t8_vtk_iterate_cells (vtkGrid, cmesh, first_tree, comm);
  }

  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  return cmesh;
}

vtkSmartPointer<vtkPointSet>
t8_vtkGrid_to_vtkPointSet (vtkSmartPointer<vtkDataSet> vtkGrid)
{
  if (vtkGrid == NULL) {
    return NULL;
  }
  /* Set points */
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New ();
  const vtkIdType num_points = vtkGrid->GetNumberOfPoints ();
  points->SetDataType (VTK_DOUBLE);
  points->SetNumberOfPoints (num_points);

  for (vtkIdType ipoint = 0; ipoint < num_points; ipoint++) {
    double vp[3];
    vtkGrid->GetPoint (ipoint, vp);
    points->SetPoint (ipoint, vp);
  }
  points->Modified ();
  T8_ASSERT (points->GetNumberOfPoints () == num_points);
  vtkSmartPointer<vtkPointSet> cloud = vtkSmartPointer<vtkPointSet>::New ();
  cloud->SetPoints (points);

  /* Map cell data to point data */
  vtkSmartPointer<vtkCellDataToPointData> c2p = vtkCellDataToPointData::New ();
  c2p->PassCellDataOff ();
  c2p->SetInputData (vtkGrid);
  c2p->Update ();
  cloud->DeepCopy (c2p->GetOutput ());
  c2p->Delete ();
  //cloud->DeepCopy (vtkPointSet::SafeDownCast (c2p->GetOutput ()));

  return cloud;
}

vtkSmartPointer<vtkDataSet>
t8_vtk_reader (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
               const vtk_file_type_t vtk_file_type)
{
  int mpisize;
  int mpirank;
  int mpiret;
  /* Get the size of the communicator and the rank of the process. */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  /* Ensure that the main-proc is a valid proc. */
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);
  T8_ASSERT (filename != NULL);
  vtk_read_success_t main_proc_read_successful = read_failure;

  vtkSmartPointer<vtkDataSet> vtkGrid = NULL;
  switch (vtk_file_type) {
  case VTK_UNSTRUCTURED_FILE:
    vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New ();
    break;
  case VTK_POLYDATA_FILE:
    vtkGrid = vtkSmartPointer<vtkPolyData>::New ();
    break;
  case VTK_PARALLEL_UNSTRUCTURED_FILE:
    vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New ();
    break;
  case VTK_PARALLEL_POLYDATA_FILE:
    vtkGrid = vtkSmartPointer<vtkPolyData>::New ();
    break;
  default:
    t8_errorf ("Filetype is not supported.\n");
    break;
  }
  T8_ASSERT (partition == 0 || (main_proc >= 0 && main_proc < mpisize));
  /* Read the file and set the pointer. */
  main_proc_read_successful = t8_file_to_vtkGrid (filename, vtkGrid, partition, main_proc, comm, vtk_file_type);

  if (!main_proc_read_successful) {
    t8_global_errorf ("Reading process (Rank %i) did not read the file successfully.\n", main_proc);
    return NULL;
  }
  else {
    return vtkGrid;
  }
}

vtkSmartPointer<vtkPointSet>
t8_vtk_reader_pointSet ([[maybe_unused]] const char *filename, [[maybe_unused]] const int partition,
                        [[maybe_unused]] const int main_proc, [[maybe_unused]] sc_MPI_Comm comm,
                        [[maybe_unused]] const vtk_file_type_t vtk_file_type)
{
  vtkSmartPointer<vtkDataSet> vtkGrid = t8_vtk_reader (filename, partition, main_proc, comm, vtk_file_type);
  return t8_vtkGrid_to_vtkPointSet (vtkGrid);
}
t8_cmesh_t
t8_vtk_reader_cmesh ([[maybe_unused]] const char *filename, [[maybe_unused]] const int partition,
                     [[maybe_unused]] const int main_proc, [[maybe_unused]] sc_MPI_Comm comm,
                     [[maybe_unused]] const vtk_file_type_t vtk_file_type)
{
  vtkSmartPointer<vtkDataSet> vtkGrid = t8_vtk_reader (filename, partition, main_proc, comm, vtk_file_type);
  if (vtkGrid != NULL) {
    const int distributed_grid = (vtk_file_type & VTK_PARALLEL_FILE) && partition;
    t8_cmesh_t cmesh = t8_vtkGrid_to_cmesh (vtkGrid, partition, main_proc, distributed_grid, comm);
    T8_ASSERT (cmesh != NULL);
    return cmesh;
  }
  else {
    t8_global_errorf ("Error translating file \n");
    return NULL;
  }
}
