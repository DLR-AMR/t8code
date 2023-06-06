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
  interpolate        *inter =
    (interpolate *) t8_forest_get_user_data (forest_old);
  T8_ASSERT (forest_old == inter->forest);
  T8_ASSERT (forest_new == inter->forest_adapt);

  t8_locidx_t         first_incoming_data =
    first_incoming + t8_forest_get_tree_element_offset (forest_new,
                                                        which_tree);
  t8_locidx_t         first_outgoing_data =
    first_outgoing + t8_forest_get_tree_element_offset (forest_old,
                                                        which_tree);

  element_data_t     *elem_data_in =
    (element_data_t *) & (inter->element_data_adapt[first_incoming_data]);
  element_data_t     *elem_data_out =
    (element_data_t *) & (inter->element_data[first_outgoing_data]);

  if (refine == 0) {
    T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
    elem_data_in->point_ids =
      sc_array_new_count (sizeof (int), elem_data_out->point_ids->elem_count);
    sc_array_copy (elem_data_in->point_ids, elem_data_out->point_ids);
    const int           num_data = inter->num_data;

    elem_data_in->average = T8_ALLOC (double *, num_data);
    for (int idata = 0; idata < num_data; idata++) {
      elem_data_in->average[idata] =
        T8_ALLOC (double, inter->data_dim[idata]);
      memcpy (elem_data_in->average[idata], elem_data_out->average[idata],
              sizeof (double) * inter->data_dim[idata]);
    }
  }
  else if (refine == 1) {
    for (int ielem = 0; ielem < num_incoming; ielem++) {
      element_data_t     *ielem_data_in =
        (element_data_t *) & (inter->element_data_adapt[first_incoming_data +
                                                        ielem]);
      t8_element_t       *elem_in =
        t8_forest_get_element_in_tree (forest_new, which_tree,
                                       first_incoming + ielem);
      ielem_data_in->average = T8_ALLOC (double *, inter->num_data);
      for (int idata = 0; idata < inter->num_data; idata++) {
        ielem_data_in->average[idata] =
          T8_ALLOC_ZERO (double, inter->data_dim[idata]);
      }
      ielem_data_in->point_ids = sc_array_new (sizeof (int));
      for (int ipoint = 0;
           (size_t) ipoint < elem_data_out->point_ids->elem_count; ipoint++) {
        int                *out_point_id =
          (int *) sc_array_index_int (elem_data_out->point_ids, ipoint);
        const double       *point =
          t8_shmem_array_get_point (inter->vtk_points, *out_point_id);
        if (t8_forest_element_point_inside
            (forest_new, which_tree, elem_in, point, 0.001) != 0) {
          int                *point_id =
            (int *) sc_array_push (ielem_data_in->point_ids);
          *point_id = *out_point_id;
          for (int idata = 0; idata < inter->num_data; idata++) {
            const double       *my_data =
              (double *) t8_shmem_array_index (inter->point_data[idata],
                                               inter->data_dim[idata] *
                                               *point_id);
            for (int idim = 0; idim < inter->data_dim[idata]; idim++) {
              ielem_data_in->average[idata][idim] += my_data[idim];
            }
          }
        }
      }
      for (int idata = 0; idata < inter->num_data; idata++) {
        for (int idim = 0; idim < inter->data_dim[idata]; idim++) {
          ielem_data_in->average[idata][idim] /=
            ielem_data_in->point_ids->elem_count;
        }
      }
    }
  }
  else {
    elem_data_in->average = T8_ALLOC (double *, inter->num_data);
    int                 num_merge_points =
      elem_data_out->point_ids->elem_count;
    size_t             *offsets = T8_ALLOC_ZERO (size_t, num_outgoing + 1);
    offsets[1] = num_merge_points;
    for (int ielem = 1; ielem < num_outgoing; ielem++) {
      element_data_t     *ielem_data_out =
        (element_data_t *) &
        (inter->element_data[first_outgoing_data + ielem]);
      num_merge_points += ielem_data_out->point_ids->elem_count;
      offsets[ielem + 1] = num_merge_points;
    }
    sc_array_init_count (elem_data_in->point_ids, sizeof (int),
                         num_merge_points);
    for (int ielem = 0; ielem < num_outgoing; ielem++) {
      element_data_t     *ielem_data_out =
        (element_data_t *) &
        (inter->element_data[first_outgoing_data + ielem]);
      sc_array_copy_into (elem_data_in->point_ids, offsets[ielem],
                          ielem_data_out->point_ids);
      for (int idata = 0; idata < inter->num_data; idata++) {
        elem_data_in->average[idata] =
          T8_ALLOC_ZERO (double, inter->data_dim[idata]);
        for (int idim = 0; idim < inter->data_dim[idata]; idim++) {
          elem_data_in->average[idata][idim] +=
            ielem_data_out->average[idata][idim] *
            ielem_data_out->point_ids->elem_count;
        }
      }
    }
    for (int idata = 0; idata < inter->num_data; idata++) {
      for (int idim = 0; idim < inter->data_dim[idata]; idim++) {
        elem_data_in->average[idata][idim] /=
          elem_data_in->point_ids->elem_count;
      }
    }
  }
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
  element_data_adapt = T8_ALLOC (element_data_t, adapt_num_elems);

  t8_forest_iterate_replace (forest_adapt, forest, t8_interpolate_replace);

  t8_forest_unref (&forest);
  forest = forest_adapt;
  forest_adapt = NULL;
  for (t8_locidx_t ielem = old_num_elems - 1; ielem >= 0; ielem--) {
    sc_array_destroy (element_data[ielem].point_ids);
    for (int idata = num_data - 1; idata >= 0; idata--) {
      T8_FREE (element_data[ielem].average[idata]);
    }
    T8_FREE (element_data[ielem].average);
  }
  T8_FREE (element_data);
  element_data = element_data_adapt;
  element_data_adapt = NULL;
}

static void
t8_pipeline (vtkSmartPointer < vtkDataSet > data, sc_MPI_Comm comm)
{
  interpolate         test = interpolate (data, comm);
  t8_debugf ("[D] Constructed class\n");

  const char          filename[12] = "class_test";

  const char         *data_name[3] = {
    "average",
    "data1",
    "data2"
  };
  t8_vtk_data_type_t  types[3];

  types[0] = T8_VTK_SCALAR;
  types[1] = T8_VTK_SCALAR;
  types[2] = T8_VTK_SCALAR;

  test.interpolate_write_vtk (filename, data_name, types);
  for (int i = 0; i < 6; i++) {
    test.adapt ();
    char                fileprefix[12];
    snprintf (fileprefix, BUFSIZ, "%s_%i", filename, i);
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
