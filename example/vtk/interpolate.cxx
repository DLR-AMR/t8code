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

int
t8_non_empty_adapt (t8_forest_t forest, t8_forest_t forest_from,
                    t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                    t8_eclass_scheme_c *ts, const int is_family,
                    const int num_elements, t8_element_t *elements[])
{
  interpolate        *user_data =
    (interpolate *) t8_forest_get_user_data (forest);
  const int           level = ts->t8_element_level (elements[0]);
  if (level == user_data->max_level && !is_family) {
    /* It is not possible to refine this level */
    return 0;
  }
  const t8_locidx_t   offset =
    t8_forest_get_tree_element_offset (forest_from, ltree_id);
  element_point_t    *elem_point =
    user_data->get_element_point (user_data->element_points,
                                  offset + lelement_id);
  const int           num_points = elem_point->num_points;
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

static void
t8_interpolate_replace (t8_forest_t forest_old,
                        t8_forest_t forest_new,
                        t8_locidx_t which_tree,
                        t8_eclass_scheme_c *ts,
                        int refine,
                        int num_outgoing,
                        t8_locidx_t first_outgoing,
                        int num_incoming, t8_locidx_t first_incoming)
{
  /* Get Metadata */
  interpolate        *inter =
    (interpolate *) t8_forest_get_user_data (forest_old);
  T8_ASSERT (forest_old == inter->forest);
  T8_ASSERT (forest_new == inter->forest_adapt);

  /* ID to Data from the new forest (Data will go in there) */
  t8_locidx_t         first_incoming_data =
    first_incoming + t8_forest_get_tree_element_offset (forest_new,
                                                        which_tree);
  /* Data from the old forest (Data is going out) */
  t8_locidx_t         first_outgoing_data =
    first_outgoing + t8_forest_get_tree_element_offset (forest_old,
                                                        which_tree);

  /* Pointer to element data */
  element_point_t    *elem_point_in =
    inter->get_element_point (inter->element_points_adapt,
                              first_incoming_data);
  element_point_t    *elem_point_out =
    inter->get_element_point (inter->element_points, first_outgoing_data);

  if (refine == 0) {
    /* point_ids array stays the same */
    T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
    /* Allocate and copy point_ids array */
    memcpy (elem_point_in, elem_point_out, sizeof (element_point_t));
  }
  else if (refine == 1) {
    /* New offsets and new num_points for each ielem_in */

    /* Tempory array to hold the point-ids for ielem_in */
    sc_array_t        **index = T8_ALLOC (sc_array_t *, num_incoming);
    for (t8_locidx_t ielem = 0; ielem < num_incoming; ielem++) {
      index[ielem] = sc_array_new (sizeof (int));
      element_point_t    *ielem_point_in =
        inter->get_element_point (inter->element_points_adapt,
                                  first_incoming_data + ielem);
      ielem_point_in->num_points = 0;
      ielem_point_in->offset = elem_point_out->offset;
    }
    const int           num_outgoing_points = elem_point_out->num_points;
    const int           offset_outgoing = elem_point_out->offset;

    /* Fill the array with the point_ids of elem_out. 
     * we pop the id from the list as soon as it is insed of ielem_in */
    sc_array_t         *point_indices =
      sc_array_new_count (sizeof (int), num_outgoing_points);
    for (int ipoint = 0; ipoint < num_outgoing_points; ipoint++) {
      const int           ipoint_id = *((int *)
                                        t8_shmem_array_index
                                        (inter->point_ids,
                                         ipoint + offset_outgoing));
      int                *point_id =
        (int *) sc_array_index_int (point_indices, ipoint);
      *point_id = ipoint_id;
    }

    for (int ipoint = 0; ipoint < num_outgoing_points; ipoint++) {
      /* Ensures that no points is associated twice. */
      const int           ipoint_id = *((int *) sc_array_pop (point_indices));
      for (t8_locidx_t ielem = 0; ielem < num_incoming; ielem++) {
        t8_element_t       *elem =
          t8_forest_get_element_in_tree (forest_new, which_tree,
                                         first_incoming + ielem);
        element_point_t    *ielem_point_in =
          inter->get_element_point (inter->element_points_adapt,
                                    first_incoming_data + ielem);
        double             *vtk_point =
          (double *) t8_shmem_array_index (inter->vtk_points, 3 * ipoint_id);
        if (t8_forest_element_point_inside
            (forest_new, which_tree, elem, vtk_point, 0.001)) {
          int                *new_point_id =
            (int *) sc_array_push (index[ielem]);
          *new_point_id = ipoint_id;
          ielem_point_in->num_points++;
          break;
        }
      }
    }
    sc_array_destroy (point_indices);

    /* Update ielem_ins offset */
    for (t8_locidx_t ielem = 1; ielem < num_incoming; ielem++) {
      element_point_t    *ielem_point_in =
        inter->get_element_point (inter->element_points_adapt,
                                  first_incoming_data + ielem);
      element_point_t    *ielem_point_in_prev =
        inter->get_element_point (inter->element_points_adapt,
                                  first_incoming_data + ielem - 1);
      ielem_point_in->offset =
        ielem_point_in_prev->offset + ielem_point_in_prev->num_points;
    }
    /* Update point_ids. */
    if (t8_shmem_array_start_writing (inter->point_ids)) {
      for (t8_locidx_t ielem = 0; ielem < num_incoming; ielem++) {
        element_point_t    *ielem_point_in =
          inter->get_element_point (inter->element_points_adapt,
                                    first_incoming_data + ielem);
        //t8_debugf("[D] write from: %i to %i\n", ielem_point_in->offset, ielem_point_in->offset + ielem_point_in->num_points);
        for (int ipoint = 0; ipoint < ielem_point_in->num_points; ipoint++) {
          int                *point_index =
            (int *) t8_shmem_array_index_for_writing (inter->point_ids,
                                                      ipoint +
                                                      ielem_point_in->offset);
          *point_index = *((int *) sc_array_index_int (index[ielem], ipoint));
        }
      }
    }
    else {
      SC_ABORTF ("Can't write point_ids");
    }
    t8_shmem_array_end_writing (inter->point_ids);
    for (int ielem = 0; ielem < num_incoming; ielem++) {
      sc_array_destroy (index[ielem]);
    }
    T8_FREE (index);
  }
  else {
    /* point-ids array stays the same. Offset of element_in is set to the
     * offset of the first element_out. Sum over the num_points to get the
     * total of points in the coarsend array. */
    elem_point_in->offset = elem_point_out->offset;
    elem_point_in->num_points = 0;
    for (int ielem = 0; ielem < num_outgoing; ielem++) {
      element_point_t    *ielem_point_out =
        inter->get_element_point (inter->element_points,
                                  first_outgoing_data + ielem);
      elem_point_in->num_points += ielem_point_out->num_points;
    }
  }
}

void
interpolate::set_elements ()
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

void
interpolate::adapt ()
{

  t8_forest_ref (forest);
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, this);
  t8_forest_set_adapt (forest_adapt, forest, t8_non_empty_adapt, 0);
  /* Balance und Ghost hier?? */
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

  t8_forest_iterate_replace (forest_adapt, forest, t8_interpolate_replace);

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

void
interpolate::partition ()
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
    }

    sc_array_destroy (idata_array);
    idata_array = idata_partion;
    idata_partion = NULL;
  }
  T8_FREE (average);
  average = partition_data;
  partition_data = NULL;
}

static void
t8_pipeline (vtkSmartPointer < vtkDataSet > data, sc_MPI_Comm comm)
{
  interpolate         test = interpolate (data, comm);
  t8_debugf ("[D] Constructed class\n");

  const char          filename[11] = "class_test";
  const char          first_data[13] = "class_test_0";

  const char         *data_name[3] = {
    "average",
    "data1",
    "data2"
  };
  t8_vtk_data_type_t  types[3];

  types[0] = T8_VTK_SCALAR;
  types[1] = T8_VTK_SCALAR;
  types[2] = T8_VTK_SCALAR;

  test.interpolate_write_vtk (first_data, data_name, types);

  for (int i = 0; i < 4; i++) {
    test.adapt ();

    //test.partition ();
    test.set_elements ();
    char                fileprefix[12];
    snprintf (fileprefix, BUFSIZ, "%s_%i", filename, i + 1);
    test.interpolate_write_vtk (fileprefix, data_name, types);
  }
  t8_debugf ("[D] wrote vtk file\n");
  return;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  mpiret = sc_MPI_Init (&argc, &argv);
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  vtkSmartPointer < vtkDataSet > vtk_grid =
    t8_vtk_reader
    ("/localdata1/knap_da/projects/t8code/t8code/test/testfiles/test_vtk_tri.vtu",
     0, 0, comm, VTK_UNSTRUCTURED_FILE);
  t8_debugf ("[D] read successfull\n");
  t8_pipeline (vtk_grid, comm);
  sc_finalize ();
  return 0;
}
