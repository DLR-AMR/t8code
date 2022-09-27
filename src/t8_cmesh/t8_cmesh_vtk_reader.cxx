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
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkTriangleFilter.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#endif

T8_EXTERN_C_BEGIN ();

/* Construct a cmesh given a filename a communicator */
t8_cmesh_t
t8_cmesh_read_from_vtk_unstructured (const char *filename,
                                     const int partition, const int main_proc,
                                     sc_MPI_Comm comm)
{
  int                 mpirank;
  int                 mpisize;
  int                 mpiret;
  int                 dim;
  t8_cmesh_t          cmesh;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  int                 main_proc_read_successful = 0;
  t8_gloidx_t         num_trees;

  T8_ASSERT (partition == 0 || (main_proc >= 0 && main_proc < mpisize));

  t8_cmesh_init (&cmesh);

#if T8_WITH_VTK
  if (!partition || mpirank == main_proc) {
    /* The Incoming data must be an unstructured Grid */
    vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
    vtkSmartPointer < vtkCellData > cellData;
    /* Prepare grid for translation */
    unstructuredGrid = t8_read_unstructured (filename);
    if (unstructuredGrid == NULL) {
      t8_errorf ("Could not read file.\n");
      if (partition) {
        main_proc_read_successful = 0;
        /* Reading failed, communicate to other processes */
        sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc,
                      comm);
      }
      return NULL;
    }

    /* Get the Data of the all cells */
    cellData = unstructuredGrid->GetCellData ();
    if (cellData == NULL) {
      t8_productionf ("No cellData found.\n");
    }

    /* Actual translation */
    num_trees =
      t8_vtk_iterate_cells (unstructuredGrid, cellData, cmesh, comm);
    dim = t8_get_dimension (unstructuredGrid);
    t8_cmesh_set_dimension (cmesh, dim);
    t8_geometry_c      *linear_geom = t8_geometry_linear_new (dim);
    t8_cmesh_register_geometry (cmesh, linear_geom);
    main_proc_read_successful = 1;
  }
  if (partition) {
    t8_gloidx_t         first_tree;
    t8_gloidx_t         last_tree;

    sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
    if (!main_proc_read_successful) {
      t8_debugf ("Main process did not read the file successfully.\n");
      t8_cmesh_destroy (&cmesh);
      return NULL;
    }
    if (mpirank == main_proc) {
      first_tree = 0;
      last_tree = num_trees - 1;
    }
    sc_MPI_Bcast (&dim, 1, sc_MPI_INT, main_proc, comm);
    t8_debugf ("[D] dim: %i\n", dim);
    sc_MPI_Bcast (&num_trees, 1, T8_MPI_GLOIDX, main_proc, comm);
    t8_cmesh_set_dimension (cmesh, dim);

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
  T8_ASSERT (cmesh != NULL);
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
  }
  return cmesh;
#else
  /* Return NULL if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return NULL;
#endif
}

t8_cmesh_t
t8_cmesh_read_from_vtk_poly (const char *filename, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  vtkSmartPointer < vtkPolyData > poly_data;
  vtkSmartPointer < vtkCellArray > cells;
  vtkSmartPointer < vtkCellData > cell_data;
  vtkSmartPointer < vtkPolyData > triangulated;
  vtkNew < vtkTriangleFilter > tri_filter;
  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);

  /* Prepare the poly-data for the translation from vtk to t8code.
   * We split all polygons (which are not supported by t8code) to
   * triangles, vertices and lines. */
  poly_data = t8_read_poly (filename);
  if (poly_data == NULL) {
    t8_errorf ("Could not read file.\n");
    return NULL;
  }
  tri_filter->SetInputData (poly_data);
  /* PolyVertex to vertex */
  tri_filter->PassVertsOn ();
  /* PolyLines to lines */
  tri_filter->PassLinesOn ();
  tri_filter->Update ();
  triangulated = tri_filter->GetOutput ();
  if (triangulated == NULL) {
    t8_errorf ("Unable to triangulate mesh.\n");
  }

  cell_data = triangulated->GetCellData ();
  if (cell_data == NULL) {
    t8_productionf ("No cellData found.\n");
  }

  t8_vtk_iterate_cells (triangulated, cell_data, cmesh, comm);
  T8_ASSERT (cmesh != NULL);
  if (cmesh != NULL) {
    t8_cmesh_commit (cmesh, comm);
    return cmesh;
  }
  else {
    return NULL;
  }
#else
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return NULL;
#endif
}

T8_EXTERN_C_END ();
