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
  double              *element_data_adapt;
};

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
                   int num_outgoing,
                   t8_locidx_t first_outgoing,
                   int num_incoming,
                   t8_locidx_t first_incoming)
{
 /* TODO */
}

static t8_forest_t
t8_build_forest (sc_MPI_Comm comm, int level)
{
  t8_cmesh_t          cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  struct t8_adapt_data adapt_data = {
    {0.5, 0.5, 0},
    0.2
  };

  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
  t8_forest_t         forest_apbg;

  t8_forest_init (&forest_apbg);
  t8_forest_set_user_data (forest_apbg, &adapt_data);
  t8_forest_set_adapt (forest_apbg, forest, t8_adapt_refine, 0);
  t8_forest_set_partition (forest_apbg, NULL, 0);
  //t8_forest_set_balance (forest_apbg, NULL, 0);
  //t8_forest_set_ghost (forest_apbg, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_apbg);

  SC_CHECK_ABORT (t8_forest_no_overlap(forest),
            "The forest has overlapping elements");

  return forest_apbg;
}


t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 int do_partition, int recursive, int do_face_ghost, void *user_data)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  if (adapt_fn != NULL) {
    t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  }
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

void
t8_test_linear_interpolation() {
  int                 level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest_new, forest_old;
  t8_scheme_cxx_t    *scheme;
  double             *data, *data_adapt;

  struct t8_adapt_data adapt_data = {
    {0.5, 0.5, 0},
    0.2,
    data,
    data_adapt
  };

  scheme = t8_scheme_new_default_cxx ();

  /* Construct a cmesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Compute the first level, such that no process is empty */
  level = t8_forest_min_nonempty_level (cmesh, scheme);
  level = SC_MAX (level, 5);

  
  t8_cmesh_ref (cmesh);
  forest_old = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  /* Set user data*/
  forest_old = t8_adapt_forest (forest_old, NULL, 0, 0, 0, &adapt_data);
  /* partition forest */
  //forest_old = t8_adapt_forest (forest_old, NULL, 1, 0, 0, &adapt_data);
  SC_CHECK_ABORT (t8_forest_no_overlap(forest_old),
                "The forest has overlapping elements");

  /* Build initial data array and set data for the local elements. */
  data = T8_ALLOC (double, t8_forest_get_local_num_elements (forest_old));
  for (size_t i = 0; i < t8_forest_get_local_num_elements (forest_old); i++) {
    data[i] = 1;
  }
  
  /* Coarse forest */
  t8_forest_ref (forest_old);
  forest_new = t8_adapt_forest (forest_old, t8_adapt_coarse, 0, 0, 0, &adapt_data);
  SC_CHECK_ABORT (t8_forest_no_overlap(forest_new),
                  "The forest has overlapping elements");

  /* Interpolate element data */
  data_adapt = T8_ALLOC (double, t8_forest_get_local_num_elements (forest_new));
  t8_forest_iterate_replace (forest_new, forest_old, t8_forest_replace);
    

  T8_FREE (data);
  T8_FREE (data_adapt);
  t8_forest_unref (&forest_old);
  t8_forest_unref (&forest_new);
  t8_cmesh_destroy (&cmesh);
}



int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;
  t8_forest_t         forest;

  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_test_linear_interpolation();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
