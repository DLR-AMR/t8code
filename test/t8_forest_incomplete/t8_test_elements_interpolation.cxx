#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_iterate.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

#include <t8_forest/t8_forest_partition.h>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  double              midpoint[3];
  double              radius;
  double              *element_data;
};

int
t8_adapt_non (t8_forest_t forest,
              t8_forest_t forest_from,
              t8_locidx_t which_tree,
              t8_locidx_t lelement_id,
              t8_eclass_scheme_c *ts,
              const int is_family,
              const int num_elements, t8_element_t *elements[])
{
  return 0;
}

int
t8_adapt_refine (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts,
                 const int is_family,
                 const int num_elements, t8_element_t *elements[])
{
  double              centroid[3];

  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  double              dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->radius) {
    return 1;
  }

  return 0;
}

int
t8_adapt_remove (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts,
                 const int is_family,
                 const int num_elements, t8_element_t *elements[])
{
  double              centroid[3];

  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  double              dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->radius) {
    return -2;
  }

  return 0;
}

int
t8_adapt_coarse (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts,
                 const int is_family,
                 const int num_elements, t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

void
t8_forest_replace (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c *ts,
                   int refine,
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming,
                   t8_locidx_t first_incoming)
{

  struct t8_adapt_data *adapt_data_new =
    (struct t8_adapt_data *) t8_forest_get_user_data (forest_new);
  const struct t8_adapt_data *adapt_data_old =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest_old);

  for (t8_locidx_t t = 0; t < which_tree; t++)
  {
    first_incoming += t8_forest_get_tree_num_elements (forest_new, t);
    first_outgoing += t8_forest_get_tree_num_elements (forest_old, t);
  }

  if (refine == 0) {
    T8_ASSERT (num_incoming == num_outgoing && num_incoming == 1);
    adapt_data_new->element_data[first_incoming] = 
      adapt_data_old->element_data[first_outgoing];
  }
  else if (refine == -1) {
    T8_ASSERT (num_incoming == 1 && num_outgoing > 0);
    adapt_data_new->element_data[first_incoming] = 0;
    for (t8_locidx_t i = 0; i < num_outgoing; i++) {
          adapt_data_new->element_data[first_incoming] +=
            adapt_data_old->element_data[first_outgoing+i];
    }
  }
  else if (refine == -2) {
    /* Do nothing */
  }
  else {
    T8_ASSERT (false);
  }

  t8_debugf("[IL] old %lf, new %lf \n ",
    adapt_data_old->element_data[first_outgoing], 
    adapt_data_new->element_data[first_incoming]);

}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 int do_partition, int recursive, int do_face_ghost, void *user_data)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  t8_forest_set_ghost (forest_new, do_face_ghost, T8_GHOST_FACES);
  if (do_partition) {
    t8_forest_set_partition (forest_new, NULL, 0);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

static void
t8_output_data (t8_forest_t forest,
                struct t8_adapt_data *data,
                const char *prefix)
{
  t8_locidx_t         num_elements =
    t8_forest_get_local_num_elements (forest);
  t8_locidx_t         ielem;
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double             *element_data = T8_ALLOC (double, num_elements);
  /* The number of user defined data fields to write. */
  int                 num_data = 1;
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Element own data");
  vtk_data.data = element_data;
  /* Copy the elment's data from our data array to the output array. */
  for (ielem = 0; ielem < num_elements; ++ielem) {
    element_data[ielem] = data->element_data[ielem];
  }

  {
    /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
     * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
     * properties of the forest to write. */
    int                 write_treeid = 1;
    int                 write_mpirank = 1;
    int                 write_level = 1;
    int                 write_element_id = 1;
    int                 write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank,
                             write_level, write_element_id, write_ghosts,
                             0, 0, num_data, &vtk_data);
  }
  T8_FREE (element_data);
}

void
t8_test_linear_interpolation() {
  int                 level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest_adapt, forest;
  t8_scheme_cxx_t    *scheme;
  double             *data, *data_adapt;

  double              radius = 0.85;
  struct t8_adapt_data adapt_data = {
    {0, 0, 0},
    radius,
    NULL
  };

  scheme = t8_scheme_new_default_cxx ();
  level = 4;

  /* Construct a cmesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_PYRAMID, sc_MPI_COMM_WORLD, 0, 0, 0);

  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  /* Set user data */
  forest = t8_adapt_forest (forest, t8_adapt_remove, 0, 0, 0, &adapt_data);
  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (double, t8_forest_get_local_num_elements (forest));
  for (t8_locidx_t i = 0; i < t8_forest_get_local_num_elements (forest); i++) {
    data[i] = 1;
  }

  /* print forest */
  adapt_data = {
    {0, 0, 0},
    radius,
    data
  };
  forest = t8_adapt_forest (forest, t8_adapt_non, 0, 0, 0, &adapt_data);
  t8_output_data (forest, &adapt_data, 
    "/home/ioannis/VBshare/paraview_export/t8_test_interpolation_forest");
 

  /* Coarse forest */
  t8_forest_ref (forest);
  struct t8_adapt_data adapt_data_adapt = {
    {0, 0, 0},
    radius,
    NULL
  };
  forest_adapt = t8_adapt_forest (forest, t8_adapt_coarse, 0, 0, 0, &adapt_data_adapt);
  SC_CHECK_ABORT (t8_forest_no_overlap(forest_adapt),
                  "The forest has overlapping elements");

  /* Interpolate element data */
  data_adapt = T8_ALLOC (double, t8_forest_get_local_num_elements (forest_adapt));
  adapt_data_adapt = {
    {0, 0, 0},
    radius,
    data_adapt
  };

  t8_forest_iterate_replace (forest_adapt, forest, t8_forest_replace);

  t8_output_data (forest_adapt, &adapt_data_adapt, 
    "/home/ioannis/VBshare/paraview_export/t8_test_interpolation_forest_adapt");
    

  T8_FREE (data);
  T8_FREE (data_adapt);
  t8_forest_unref (&forest);
  t8_forest_unref (&forest_adapt);
}



int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (SC_LP_DEFAULT);

  t8_test_linear_interpolation();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
