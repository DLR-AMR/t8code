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

#include <t8.h>
#include <t8_cmesh.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_vtk_reader.hxx>

#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest.h>
#if T8_WITH_VTK
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#endif

typedef struct
{
  sc_array_t         *point_ids;
  double              average;
} element_data_t;

typedef struct
{
  t8_shmem_array_t    vtk_points;
  t8_shmem_array_t    point_data;
  int                 num_data;
  int                 num_points;
  int                 data_size;
  int                 max_level;
  element_data_t     *element_data;
} t8_user_data_t;

int
t8_non_empty_adapt (t8_forest_t forest, t8_forest_t forest_from,
                    t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                    t8_eclass_scheme_c *ts, const int is_family,
                    const int num_elements, t8_element_t *elements[])
{
  t8_user_data_t     *user_data;
  user_data = (t8_user_data_t *) t8_forest_get_user_data (forest);
  const int           level = ts->t8_element_level (elements[0]);
  if (level == user_data->max_level && !is_family) {
    /* It is not possible to refine this level */
    return 0;
  }
  const t8_locidx_t   offset =
    t8_forest_get_tree_element_offset (forest_from, ltree_id);
  element_data_t     *elem_data =
    &(user_data->element_data[offset + lelement_id]);
  const int           num_points = (int) elem_data->point_ids->elem_count;
  if (num_points > 1) {
    return 1;
  }
  else {
    return 0;
  }
}

double             *
t8_shmem_array_get_point (t8_shmem_array_t array, int index)
{
  T8_ASSERT (t8_shmem_array_is_initialized (array));
  T8_ASSERT (0 <= index
             && (size_t) index < t8_shmem_array_get_elem_count (array) / 3);

  return (double *) t8_shmem_array_index (array, 3 * index);
}

/* We start without any refinement, hence all indices are put into the root-element */
static void
t8_init_element_data (t8_user_data_t * user_data)
{
  user_data->element_data = T8_ALLOC (element_data_t, 1);
  user_data->element_data[0].point_ids =
    sc_array_new_count (sizeof (int), user_data->num_points);
  user_data->element_data[0].average = 0;
  for (int ipoint = 0; ipoint < user_data->num_points; ipoint++) {
    int                *point_id =
      (int *) sc_array_index_int (user_data->element_data[0].point_ids,
                                  ipoint);
    *point_id = ipoint;
    const double       *my_data =
      (double *) t8_shmem_array_index (user_data->point_data, ipoint);
    user_data->element_data[0].average += *my_data;
  }
  user_data->element_data[0].average /= user_data->num_points;
  return;
}

static void
t8_destroy_element_data (t8_forest_t forest, t8_user_data_t * user_data)
{
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest);
  t8_debugf ("[D] num_elements: %i\n", num_elements);
  for (t8_locidx_t ielem = num_elements - 1; ielem >= 0; ielem--) {
    t8_debugf ("[D] ielem: %i\n", ielem);
    sc_array_destroy_null (&(user_data->element_data[ielem].point_ids));
  }
  t8_debugf ("[D] sc_array destroyed\n");
  T8_FREE (user_data->element_data);
}

static void
t8_init_user_data (t8_user_data_t * user_data,
                   vtkSmartPointer < vtkDataSet > data, sc_MPI_Comm comm)
{
  /* Get the point-coordinates */
  int                 num_local_points = 0;
  /* Get number of local points, might be 0 */
  if (data != NULL) {
    num_local_points = (int) data->GetNumberOfPoints ();
  }
  double             *local_points = T8_ALLOC (double, 3 * num_local_points);

  /* Write point-coordinates */
  for (vtkIdType ipoint = 0; ipoint < num_local_points; ipoint++) {
    data->GetPoint (ipoint, &(local_points[3 * ipoint]));
  }

  /*Get global number of points */
  int                 num_global_points = 0;
  int                 mpiret = sc_MPI_Allreduce ((void *) (&num_local_points),
                                                 (void
                                                  *) (&num_global_points), 1,
                                                 sc_MPI_INT,
                                                 sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  /* Fill shmem with point-data. */
  t8_shmem_init (comm);
  t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);

  t8_shmem_array_init (&(user_data->vtk_points), sizeof (double),
                       3 * num_global_points, comm);

  t8_shmem_array_allgatherv ((void *) local_points, 3 * num_local_points,
                             sc_MPI_DOUBLE, user_data->vtk_points,
                             sc_MPI_DOUBLE, comm);

  if (data != NULL) {
    /* Map cell data to point data */
    vtkSmartPointer < vtkCellDataToPointData > c2p =
      vtkCellDataToPointData::New ();
    c2p->PassCellDataOff ();
    c2p->SetInputData (data);
    c2p->Update ();
  }

  int                 num_data_arrays = 0;

  vtkSmartPointer < vtkPointData > point_data = NULL;
  if (data != NULL) {
    point_data = data->GetPointData ();
    num_data_arrays = point_data->GetNumberOfArrays ();
  }
  t8_debugf ("[D] num_data_arrays: %i\n", num_data_arrays);
  if (num_data_arrays > 0) {
    // Currently only a single data array.
    int                 tuple_dim = 0;
    if (data != NULL) {
      vtkDataArray       *my_data = point_data->GetArray (0);
      tuple_dim = my_data->GetNumberOfComponents ();
      t8_debugf ("[D] data has dim: %i\n", tuple_dim);
    }
    mpiret = sc_MPI_Bcast ((void *) &tuple_dim, 1, sc_MPI_INT, 0, comm);
    t8_debugf ("[D] tuple_dim2: %i\n", tuple_dim);
    SC_CHECK_MPI (mpiret);
    double             *data_array =
      T8_ALLOC (double, tuple_dim * num_local_points);
    if (data != NULL) {
      vtkDataArray       *my_data = point_data->GetArray (0);
      const size_t        data_size = sizeof (double) * tuple_dim;
      for (vtkIdType ipoint = 0; ipoint < num_local_points; ipoint++) {
        my_data->GetTuple (ipoint, &(data_array[ipoint * tuple_dim]));
      }
    }

    t8_shmem_array_init (&(user_data->point_data), sizeof (double),
                         tuple_dim * num_global_points, comm);
    t8_shmem_array_allgatherv ((void *) data_array,
                               tuple_dim * num_local_points, sc_MPI_DOUBLE,
                               user_data->point_data, sc_MPI_DOUBLE, comm);
    t8_debugf ("[D] point_data allgatherved\n");
    user_data->data_size = tuple_dim;
    user_data->num_data = num_data_arrays;
    user_data->num_points = num_global_points;
    user_data->max_level = 6;
    T8_FREE (data_array);
  }
  T8_FREE (local_points);

  t8_init_element_data (user_data);
  t8_debugf ("[D] successfully initialized user_data\n");
}

static void
t8_user_data_destroy (t8_forest_t forest, sc_MPI_Comm comm)
{
  t8_user_data_t     *user_data =
    (t8_user_data_t *) t8_forest_get_user_data (forest);
  t8_shmem_array_destroy (&(user_data->point_data));
  t8_shmem_array_destroy (&(user_data->vtk_points));
  t8_debugf ("[D] destroy shmem\n");
  t8_destroy_element_data (forest, user_data);
  t8_debugf ("[D] Destroy user Data\n");
  return;
}

static void
t8_pipeline (t8_forest_t forest, vtkSmartPointer < vtkDataSet > data,
             sc_MPI_Comm comm)
{
  t8_user_data_t      user_data;
  t8_init_user_data (&user_data, data, comm);
  t8_forest_set_user_data (forest, (void *) &user_data);
  for (int ipoint = 0; ipoint < user_data.num_points; ipoint++) {
    const double       *my_point =
      t8_shmem_array_get_point (user_data.vtk_points, ipoint);
    t8_debugf ("[D] point: %3.3f %3.3f %3.3f\n", my_point[0], my_point[1],
               my_point[2]);
    const double       *my_data =
      (double *) t8_shmem_array_index (user_data.point_data, ipoint);
    for (int ituple = 0; ituple < user_data.data_size; ituple++) {
      t8_debugf ("%2.2f ", my_data[ituple]);
    }
    t8_debugf ("\n");
  }
  for (int ipoint = 0; ipoint < user_data.num_points; ipoint++) {
    int                *point =
      (int *) sc_array_index_int (user_data.element_data->point_ids, ipoint);
    t8_debugf ("[D] point_id: %i\n", *point);
  }
  t8_debugf ("[D] average: %f\n", user_data.element_data[0].average);
  t8_user_data_destroy (forest, comm);
  return;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level = 0;
  mpiret = sc_MPI_Init (&argc, &argv);
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_cmesh_t          cmesh =
    t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  vtkSmartPointer < vtkDataSet > vtk_grid =
    t8_vtk_reader
    ("/localdata1/knap_da/projects/t8code/t8code/test/testfiles/test_vtk_tri.vtu",
     0, 0, comm, VTK_UNSTRUCTURED_FILE);

  t8_pipeline (forest, vtk_grid, comm);
  t8_forest_unref (&forest);
  sc_finalize ();
  return 0;
}
