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
#include <example/vtk/interpolate.hxx>

#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
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
  int                 point_id;
  int                 has_been_set;
} interpolate_search_point_t;

static int
t8_search_callback (t8_forest_t forest,
                    t8_locidx_t ltreeid,
                    const t8_element_t *element,
                    const int is_leaf,
                    t8_element_array_t *leaf_elements,
                    t8_locidx_t tree_leaf_index,
                    void *query, size_t query_index)
{
  T8_ASSERT (query == NULL);
  /* Get user_data?? */
  return 1;
}

static int
t8_search_query_callback (t8_forest_t forest,
                          t8_locidx_t ltreeid,
                          const t8_element_t *element,
                          const int is_leaf,
                          t8_element_array_t *leaf_elements,
                          t8_locidx_t
                          tree_leaf_index, void *query, size_t query_index)
{
  interpolate_search_point_t *point = (interpolate_search_point_t *) query;
  const double        tolerance = 1e-8;
  MeshAdapter        *adapter =
    (MeshAdapter *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapter != NULL);
  double             *vtk_point =
    (double *) t8_shmem_array_index (adapter->GetVTKPoints (),
                                     3 * point->point_id);
  const int           point_is_inside =
    t8_forest_element_point_inside (forest,
                                    ltreeid, element, vtk_point, tolerance);
  if (point_is_inside) {
    if (is_leaf && !point->has_been_set) {
      /* The point is inside the element, and the element is a leaf. */
      t8_locidx_t         element_index =
        t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
      sc_array_t         *points_per_elem =
        adapter->get_point_id_per_element (element_index);
      int                *new_point = (int *) sc_array_push (points_per_elem);
      *new_point = point->point_id;
      point->has_been_set = 1;
    }
    return 1;
  }
  /* The point is not inside the element, deactivate the query */
  return 0;
}

MeshAdapter::MeshAdapter (vtkSmartPointer < vtkDataSet > input,
                          sc_MPI_Comm comm)
{
  levelMaximum = 6;

  /* Set MPI Communicator */
  m_iComm = comm;

  /* Create just a point set from the input and ignore the topology */
  vtkSmartPointer < vtkPointSet > pointSet =
    t8_vtkGrid_to_vtkPointSet (input);

  /* Bounding Box */
    /*-----------------------------------------------------------------------------------------------*/

  /* Reduce the global bounds. Each MPI rank just has a portion of data */
  vtkBoundingBox      globalBounds;
  const vtkBoundingBox localBounds (pointSet->GetBounds ());

  /* Collect and reduce bounds from all processes */
  //vtkMultiProcessController::GetGlobalController()->AllReduce(localBounds, globalBounds);

  /* Get and Convert Boundaries */
  double              vtkBounds[6];
  double             *t8Bounds = T8_ALLOC_ZERO (double, 24);

  globalBounds.GetBounds (vtkBounds);

  ConvertVTKBoundariesToT8Boundaries (vtkBounds, t8Bounds);

  /* Initialize */
    /*-----------------------------------------------------------------------------------------------*/

  int                 dimensionCell = GetDimensionCell (input);
  //int dimensionData = GetDimensionData(vtkBounds);
  int                 dimensionData = 2;

  t8_debugf ("( t8SeriesWriter ) >> Dimension Cell: %i\n", dimensionCell);
  t8_debugf ("( t8SeriesWriter ) >> Dimension Data: %i\n", dimensionData);

  /* Create Cmesh based on PointCloud */
  t8_eclass_t         eclass = T8_ECLASS_LINE;

  if (dimensionData == 2) {
    eclass = T8_ECLASS_QUAD;
  }

  if (dimensionData == 3) {
    eclass = T8_ECLASS_HEX;
  }

  /* Create Cmesh with correct eclass */
  t8_cmesh_t          cmesh =
    t8_cmesh_new_hypercube_pad (eclass, m_iComm, t8Bounds, 2, 3, 1);

  T8_FREE (t8Bounds);

  /* Create Initial Forest */
  t8_forest_t         forest_tmp =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0, comm);

  /* Initialize, Partition and Commit */
  t8_forest_init (&forest);
  t8_forest_set_partition (forest, forest_tmp, 1);
  t8_forest_commit (forest);

  /* Get Local Number of Elements */
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest);
  const t8_gloidx_t   num_global_elements =
    t8_forest_get_global_num_elements (forest);

  /* TODO: Write vtk-point data stuff only where forest exits. */

  /* Get MPI Rank and process if Rank 0 */
  int                 mpirank;
  int                 mpisize;
  int                 mpiret = sc_MPI_Comm_rank (m_iComm, &mpirank);

  /* Check if MPI is fine */
  SC_CHECK_MPI (mpiret);

  mpisize = sc_MPI_Comm_size (m_iComm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* Get number of local points, might be 0 */
  int                 num_local_points = 0;

  num_local_points = (int) pointSet->GetNumberOfPoints ();
  t8_debugf ("[D] num_local_points: %i\n", num_local_points);

  double             *local_points = T8_ALLOC (double, 3 * num_local_points);

  /* Write point-coordinates */
  for (vtkIdType i = 0; i < num_local_points; i++) {
    pointSet->GetPoint (i, &(local_points[3 * i]));
  }

  /*Get global number of points */
  num_global_points = 0;
  mpiret = sc_MPI_Allreduce ((void *) (&num_local_points),
                             (void *) (&num_global_points), 1,
                             sc_MPI_INT, sc_MPI_SUM, m_iComm);

  SC_CHECK_MPI (mpiret);

  /* Fill shmem with point-data. */
  t8_shmem_init (m_iComm);
  t8_shmem_set_type (m_iComm, T8_SHMEM_BEST_TYPE);

  t8_shmem_array_init (&vtk_points, sizeof (double),
                       3 * num_global_points, m_iComm);

  t8_shmem_array_allgatherv ((void *) local_points, 3 * num_local_points,
                             sc_MPI_DOUBLE, vtk_points,
                             sc_MPI_DOUBLE, m_iComm);

  num_data = 0;

  vtkSmartPointer < vtkPointData > vtk_point_data = NULL;
  if (num_local_points > 0) {
    vtk_point_data = pointSet->GetPointData ();
    num_data = vtk_point_data->GetNumberOfArrays ();
  }

  /* Currently proc 0 is the main proc. */
  mpiret = sc_MPI_Bcast ((void *) &num_data, 1, sc_MPI_INT, 0, m_iComm);

  SC_CHECK_MPI (mpiret);
  point_data = T8_ALLOC (t8_shmem_array_t, num_data);

  data_dim = T8_ALLOC_ZERO (int, num_data);
  if (num_local_points > 0) {
    for (int i = 0; i < num_data; i++) {
      data_dim[i] =
        (int) (vtk_point_data->GetArray (i)->GetNumberOfComponents ());
    }
  }
  mpiret = sc_MPI_Bcast ((void *) data_dim, num_data, sc_MPI_INT, 0, m_iComm);

  SC_CHECK_MPI (mpiret);

  for (int i = 0; i < num_data; i++) {
    /* Allocate memory depending on data dimension. */
    double             *data_array =
      T8_ALLOC (double, data_dim[i] * num_local_points);

    if (num_local_points > 0) {
      vtkDataArray       *my_data = vtk_point_data->GetArray (i);

      for (vtkIdType j = 0; j < num_local_points; j++) {
        my_data->GetTuple (j, &(data_array[j * data_dim[i]]));
      }
    }
    t8_shmem_array_init (&(point_data[i]), sizeof (double),
                         data_dim[i] * num_global_points, m_iComm);
    t8_shmem_array_allgatherv ((void *) data_array,
                               data_dim[i] * num_local_points, sc_MPI_DOUBLE,
                               point_data[i], sc_MPI_DOUBLE, m_iComm);
    T8_FREE (data_array);
  }

  t8_shmem_array_init (&point_ids, sizeof (int), num_global_points, m_iComm);
  sc_array_t         *point_queries;
  point_queries =
    sc_array_new_count (sizeof (interpolate_search_point_t),
                        num_global_points);
  for (int ipoint = 0; ipoint < num_global_points; ipoint++) {
    interpolate_search_point_t *query =
      (interpolate_search_point_t *) sc_array_index_int (point_queries,
                                                         ipoint);
    query->point_id = ipoint;
    query->has_been_set = 0;
  }

  if (t8_shmem_array_start_writing (point_ids)) {
    for (int ipoint = 0; ipoint < num_global_points; ipoint++) {
      int                *index =
        (int *) t8_shmem_array_index_for_writing (point_ids, ipoint);

      /* initialize points with -1 to ensure, that every points is associated
       * once and only once with an element. */
      *index = -1;
    }
  }
  t8_shmem_array_end_writing (point_ids);

  points_per_element = T8_ALLOC (sc_array_t *, num_elements);
  for (int ielem = 0; ielem < num_elements; ielem++) {
    points_per_element[ielem] = sc_array_new (sizeof (int));
  }

  average = T8_ALLOC (sc_array_t *, num_data);
  for (int idata = 0; idata < num_data; idata++) {
    average[idata] =
      sc_array_new_count (sizeof (element_data_t), num_elements);
  }

  element_points =
    sc_array_new_count (sizeof (element_point_t), num_elements);
  element_point_t    *first_elem =
    (element_point_t *) sc_array_index_int (element_points, 0);
  first_elem->offset = 0;
  if (num_global_elements == 1 && num_elements == 1) {
    if (t8_shmem_array_start_writing (point_ids)) {
      for (int ipoint = 0; ipoint < num_global_points; ipoint++) {
        int                *index =
          (int *) t8_shmem_array_index_for_writing (point_ids, ipoint);
        *index = ipoint;
      }
    }
    t8_shmem_array_end_writing (point_ids);
  }
  else {
    t8_forest_set_user_data (forest, (void *) this);

    t8_forest_search (forest, t8_search_callback, t8_search_query_callback,
                      point_queries);

    int                 num_set_points = 0;
    for (int ielem = 0; ielem < num_elements; ielem++) {
      for (int ipoint = 0; ipoint < points_per_element[ielem]->elem_count;
           ipoint++) {
        int                *ipoint_id =
          (int *) sc_array_index (points_per_element[ielem], ipoint);
        T8_ASSERT (*ipoint_id < 121);
      }
      num_set_points += points_per_element[ielem]->elem_count;
    }

    int                *set_point_ids = T8_ALLOC_ZERO (int, num_set_points);
    int                 dest = 0;
    for (int ielem = 0; ielem < num_elements; ielem++) {
      const int           num_points_per_elem =
        points_per_element[ielem]->elem_count;
      memcpy (&set_point_ids[dest], points_per_element[ielem]->array,
              num_points_per_elem * sizeof (int));
      dest += num_points_per_elem;
    }

    t8_shmem_array_allgatherv (set_point_ids, num_set_points, sc_MPI_INT,
                               point_ids, sc_MPI_INT, comm);

    if (mpisize > 0) {
      t8_shmem_array_t    offsets;      /* Offsets of the point-ids */
      t8_debugf ("[D] mpirank: %i, mpisize: %i\n", mpirank, mpisize);
      t8_shmem_array_init (&offsets, sizeof (int), mpisize, comm);

      t8_shmem_array_prefix (&dest, offsets, 1, sc_MPI_INT, sc_MPI_SUM, comm);

      first_elem->offset = *(int *) t8_shmem_array_index (offsets, mpirank);
      t8_shmem_array_destroy (&offsets);
    }
    first_elem->num_points = points_per_element[0]->elem_count;

    /* Update ielem_ins offset */
    for (t8_locidx_t ielem = 1; ielem < num_elements; ielem++) {
      element_point_t    *ielem_point_in =
        get_element_point (GetElementPoints (), ielem);
      element_point_t    *ielem_point_in_prev =
        get_element_point (GetElementPoints (), ielem - 1);
      ielem_point_in->offset =
        ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
      ielem_point_in->num_points = points_per_element[ielem]->elem_count;
    }

    for (int ielem = num_elements - 1; ielem >= 0; ielem--) {
      sc_array_destroy (points_per_element[ielem]);
    }
    T8_FREE (set_point_ids);
    T8_FREE (points_per_element);
  }

  this->SetElements ();
  sc_array_destroy (point_queries);
  T8_FREE (local_points);
  t8_forest_set_user_data (forest, (void *) this);
}

void
MeshAdapter::WritePVTU (const char *fileprefix, const char **data_names,
                        const t8_vtk_data_type_t * data_types)
{
  t8_debugf ("[D] WritePVTU\n");

  t8_vtk_data_field_t *vtk_data = T8_ALLOC (t8_vtk_data_field_t, num_data);
  double            **data_array = T8_ALLOC (double *, num_data);
  const t8_locidx_t   num_local_elements =
    t8_forest_get_local_num_elements (forest);

  /*Linearise each data-field */
  if (num_local_elements > 0) {
    for (int idata = 0; idata < num_data; idata++) {
      /*for(t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++){
         ("[D] elem: %i\n", ielem);
         element_point_t *elem_p = get_element_point(element_points, ielem);
         for(int ipoint = elem_p->offset; ipoint < elem_p->offset + elem_p->num_points; ipoint++){
         int *point_id = (int *) t8_shmem_array_index(point_ids, ipoint);
         printf("%i, ", *point_id);
         }
         } */
      sc_array_t         *iaverage = average[idata];

      data_array[idata] =
        T8_ALLOC_ZERO (double, num_local_elements * data_dim[idata]);
      for (t8_locidx_t ielem = 0; ielem < num_local_elements; ielem++) {

        const int           current_dim = data_dim[idata];
        element_data_t     *ielem_data =
          (element_data_t *) sc_array_index_int (iaverage, ielem);

        for (int idim = 0; idim < current_dim; idim++) {

          data_array[idata][current_dim * ielem + idim] =
            ielem_data->data[idim];
          /* printf("%f\n", ielem_data->data[idim]); */
        }
      }
      snprintf (vtk_data[idata].description, BUFSIZ, "%s", data_names[idata]);
      vtk_data[idata].data = data_array[idata];
      vtk_data[idata].type = data_types[idata];
    }
  }

  if (t8_forest_write_vtk_ext
      (forest, fileprefix, 1, 1, 1, 1, 0, 0, 0, num_data, vtk_data)) {
    t8_debugf ("[Interpolate] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Interpolate] Error writing to files %s\n", fileprefix);
  }
  for (int idata = num_data - 1; idata >= 0; idata--) {
    T8_FREE (data_array[idata]);
  }

  T8_FREE (data_array);
  T8_FREE (vtk_data);
}

//----------------------------------------------------------------------------
void
MeshAdapter::SetElements ()
{
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest);

  for (int idata = 0; idata < num_data; idata++) {
    const int           idata_dim = data_dim[idata];
    for (t8_locidx_t ielem = 0; ielem < num_elements; ielem++) {
      element_data_t     *ielem_data =
        get_element_data (average, ielem, idata);
      ielem_data->data = T8_ALLOC_ZERO (double, idata_dim);
      const int           offset =
        *(get_element_point_offset (element_points, ielem));
      const int           num_points =
        *(get_element_num_points (element_points, ielem));
      for (int ipoint = offset; ipoint < offset + num_points; ipoint++) {
        const int           ipoint_id =
          *((int *) t8_shmem_array_index (point_ids, ipoint));
        const double       *my_data =
          (double *) t8_shmem_array_index (point_data[idata],
                                           ipoint_id * idata_dim);

        for (int idim = 0; idim < idata_dim; idim++) {
          ielem_data->data[idim] += my_data[idim];
        }
      }

      for (int idim = 0; idim < idata_dim; idim++) {
        ielem_data->data[idim] /= ((num_points == 0) ? 1 : num_points);
      }
    }
  }                             /*
                                   for (int idata = 0; idata < num_data; idata++) {
                                   for (t8_locidx_t ielem = 0; ielem < num_elements; ielem++) {
                                   element_data_t     *ielem_data =
                                   get_element_data (average, ielem, idata);
                                   t8_debugf ("[D] idata: %i elem %i\n", idata, ielem);
                                   for (int idim = 0; idim < data_dim[idata]; idim++) {
                                   printf ("%f, ", ielem_data->data[idim]);
                                   }
                                   printf ("\n");

                                   }
                                   } */
}

//----------------------------------------------------------------------------
void
MeshAdapter::Adapt (t8_forest_adapt_t adaptCallback,
                    t8_forest_replace_t replaceCallback)
{
  sc_MPI_Comm         comm = this->m_iComm;
  t8_forest_ref (forest);
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, (void *) this);

  /* Set the Adapt Callback */

  t8_forest_set_adapt (forest_adapt, forest, adaptCallback, 0);
  /* Todo: Balance und Ghost hier?? */
  t8_forest_commit (forest_adapt);

  const t8_locidx_t   adapt_num_elems =
    t8_forest_get_local_num_elements (forest_adapt);
  const t8_locidx_t   old_num_elems =
    t8_forest_get_local_num_elements (forest);

  element_points_adapt =
    sc_array_new_count (sizeof (element_point_t), adapt_num_elems);

  average_adapt = T8_ALLOC (sc_array_t *, num_data);

  for (int idata = 0; idata < num_data; idata++) {
    average_adapt[idata] =
      sc_array_new_count (sizeof (element_data_t), adapt_num_elems);
  }

  points_per_element = T8_ALLOC (sc_array_t *, adapt_num_elems);
  for (int ielem = 0; ielem < adapt_num_elems; ielem++) {
    points_per_element[ielem] = sc_array_new (sizeof (int));
  }

  /* Set the Replace Callback */
  t8_forest_iterate_replace (forest_adapt, forest, replaceCallback);
  int                 local_num_points = 0;
  for (int ielem = 0; ielem < adapt_num_elems; ielem++) {
    local_num_points += points_per_element[ielem]->elem_count;
  }
  int                *set_point_ids = T8_ALLOC_ZERO (int, local_num_points);
  int                 dest = 0;
  for (int ielem = 0; ielem < adapt_num_elems; ielem++) {
    const int           num_points_per_elem =
      points_per_element[ielem]->elem_count;
    memcpy (&set_point_ids[dest], points_per_element[ielem]->array,
            num_points_per_elem * sizeof (int));
    dest += num_points_per_elem;
  }

  t8_shmem_array_allgatherv (set_point_ids, local_num_points, sc_MPI_INT,
                             point_ids, sc_MPI_INT, comm);

  element_point_t    *first_elem =
    (element_point_t *) sc_array_index_int (element_points_adapt, 0);
  first_elem->offset = 0;

  int                 mpisize;
  int                 mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  int                 mpirank;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  if (mpisize > 0) {
    t8_shmem_array_t    offsets;
    t8_shmem_array_init (&offsets, sizeof (int), mpisize, comm);

    t8_shmem_array_prefix (&dest, offsets, 1, sc_MPI_INT, sc_MPI_SUM, comm);

    first_elem->offset = *(int *) t8_shmem_array_index (offsets, mpirank);
    t8_shmem_array_destroy (&offsets);
  }
  first_elem->num_points = points_per_element[0]->elem_count;
  /* Update ielem_ins offset */
  for (t8_locidx_t ielem = 1; ielem < adapt_num_elems; ielem++) {
    element_point_t    *ielem_point_in =
      get_element_point (GetElementPointsAdapt (), ielem);
    const element_point_t *ielem_point_in_prev =
      get_element_point (GetElementPointsAdapt (), ielem - 1);
    ielem_point_in->offset =
      ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
    ielem_point_in->num_points = points_per_element[ielem]->elem_count;
  }

  T8_FREE (set_point_ids);

  /* Set forest to forest_adapt and delete all data from forest. */
  t8_forest_unref (&forest);
  forest = forest_adapt;
  forest_adapt = NULL;
  for (int idata = num_data - 1; idata >= 0; idata--) {
    sc_array_t         *idata_array = average[idata];
    for (t8_locidx_t ielem = old_num_elems - 1; ielem >= 0; ielem--) {
      element_data_t     *ielem_data =
        (element_data_t *) sc_array_index_int (idata_array, ielem);
      T8_FREE (ielem_data->data);
    }
    sc_array_destroy (idata_array);
  }

  T8_FREE (average);
  average = average_adapt;
  average_adapt = NULL;
  sc_array_destroy (element_points);
  element_points = element_points_adapt;
  element_points_adapt = NULL;
}

//----------------------------------------------------------------------------
void
MeshAdapter::ConvertVTKBoundariesToT8Boundaries (double *vtkBounds,
                                                 double *t8Bounds)
{

  t8Bounds[0] = 0;
  t8Bounds[1] = 0;
  t8Bounds[2] = 0;

  t8Bounds[3] = 10;
  t8Bounds[4] = 0;
  t8Bounds[5] = 0;

  t8Bounds[6] = 0;
  t8Bounds[7] = 10;
  t8Bounds[8] = 0;

  t8Bounds[9] = 10;
  t8Bounds[10] = 10;
  t8Bounds[11] = 0;

  // XMin, XMax, YMin, YMax, ZMin, ZMax

  t8Bounds[12] = vtkBounds[0];
  t8Bounds[13] = vtkBounds[2];
  t8Bounds[14] = vtkBounds[5];

  t8Bounds[15] = vtkBounds[1];
  t8Bounds[16] = vtkBounds[2];
  t8Bounds[17] = vtkBounds[5];

  t8Bounds[18] = vtkBounds[0];
  t8Bounds[19] = vtkBounds[3];
  t8Bounds[20] = vtkBounds[5];

  t8Bounds[21] = vtkBounds[1];
  t8Bounds[22] = vtkBounds[3];
  t8Bounds[23] = vtkBounds[5];
}

//----------------------------------------------------------------------------
int
MeshAdapter::GetDimensionData (double *vtkBounds)
{
  int                 dimensionData = -1;

  /* Get DataSet Dimension */
  auto                dimX =
    ((int) (vtkBounds[1] - vtkBounds[0]) > 0) ? 1 : 0;
  auto                dimY =
    ((int) (vtkBounds[3] - vtkBounds[2]) > 0) ? 1 : 0;
  auto                dimZ =
    ((int) (vtkBounds[5] - vtkBounds[4]) > 0) ? 1 : 0;

  dimensionData = dimX + dimY + dimZ;

  if (dimensionData < 1 || dimensionData > 3) {
    t8_global_errorf
      ("( t8SeriesWriter ) >>  Error: Dimension of Data is not 1, 2 or 3\n");
  }

  return dimensionData;
}

//----------------------------------------------------------------------------
int
MeshAdapter::GetDimensionCell (vtkSmartPointer < vtkDataSet > dataSet)
{
  int                 dimensionCell = -1;

  /* Get Cell Dimension */
  if (dataSet->GetNumberOfCells () > 0) {
    dimensionCell = dataSet->GetCell (0)->GetCellDimension ();
  }

  if (dimensionCell < 1 || dimensionCell > 3) {
    t8_global_errorf
      ("( t8SeriesWriter ) >>  Error: Dimension of Cell is not 1, 2 or 3\n");
  }

  return dimensionCell;
}

//----------------------------------------------------------------------------
MeshAdapter::~MeshAdapter ()
{
  t8_shmem_array_destroy (&vtk_points);
  t8_shmem_array_destroy (&point_ids);

  for (int idata = num_data - 1; idata >= 0; idata--) {
    t8_shmem_array_destroy (&point_data[idata]);
  }

  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest);

  for (int idata = num_data - 1; idata >= 0; idata--) {
    sc_array_t         *idata_array = average[idata];

    for (t8_locidx_t ielem = num_elements - 1; ielem >= 0; ielem--) {
      element_data_t     *ielem_data =
        (element_data_t *) sc_array_index_int (idata_array, ielem);
      T8_FREE (ielem_data->data);
    }

    sc_array_destroy (idata_array);
  }

  T8_FREE (average);
  T8_FREE (data_dim);
  T8_FREE (point_data);

  sc_array_destroy (element_points);
  t8_forest_unref (&forest);
  t8_debugf ("[D] destroyed class\n");
}

void
MeshAdapter::partition ()
{
  t8_forest_t         forest_partition;
  t8_forest_ref (forest);
  t8_forest_init (&forest_partition);
  t8_forest_set_user_data (forest_partition, this);
  t8_forest_set_partition (forest_partition, forest, 0);
  /*TODO: Ghosts?? */
  t8_forest_commit (forest_partition);

  const t8_locidx_t   num_local_elements =
    t8_forest_get_local_num_elements (forest);
  const t8_locidx_t   num_local_elements_part =
    t8_forest_get_local_num_elements (forest_partition);

  sc_array_t          points_view;
  sc_array_init_view (&points_view, element_points, 0, num_local_elements);
  sc_array_t         *data_view = T8_ALLOC (sc_array, num_data);

  sc_array_t         *partition_points
    = sc_array_new_count (sizeof (element_point_t), num_local_elements_part);
  sc_array_t        **partition_data = T8_ALLOC (sc_array_t *, num_data);

  for (int idata = 0; idata < num_data; idata++) {
    partition_data[idata] =
      sc_array_new_count (sizeof (element_data_t), num_local_elements_part);
  }

  sc_array_t          points_part_view;
  sc_array_init_view (&points_part_view, partition_points, 0,
                      num_local_elements_part);
  sc_array_t         *data_part_view = T8_ALLOC (sc_array, num_data);

  for (int idata = 0; idata < num_data; idata++) {
    sc_array_init_view (&data_view[idata], average[idata], 0,
                        num_local_elements);
    sc_array_init_view (&data_part_view[idata], partition_data[idata], 0,
                        num_local_elements_part);
  }

  t8_forest_partition_data (forest, forest_partition, &points_view,
                            &points_part_view);

  for (int idata = 0; idata < num_data; idata++) {
    t8_forest_partition_data (forest, forest_partition, &data_view[idata],
                              &data_part_view[idata]);
  }
  t8_forest_unref (&forest);

  forest = forest_partition;
  forest_partition = NULL;

  sc_array_destroy (element_points);
  element_points = partition_points;
  partition_points = NULL;

  for (int idata = num_data - 1; idata >= 0; idata--) {
    sc_array_t         *idata_array = average[idata];
    sc_array_t         *idata_partion = partition_data[idata];
    for (t8_locidx_t ielem = num_local_elements - 1; ielem >= 0; ielem--) {
      element_data_t     *ielem_data =
        (element_data_t *) sc_array_index_int (idata_array, ielem);
      T8_FREE (ielem_data->data);
    }

    sc_array_destroy (idata_array);
    idata_array = idata_partion;
    idata_partion = NULL;
  }
  T8_FREE (average);
  average = partition_data;
  partition_data = NULL;
}
