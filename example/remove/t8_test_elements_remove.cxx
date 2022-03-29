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
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

struct t8_adapt_data
{
  double  midpoint[6][3];
};

int
t8_adapt_callback_init (t8_forest_t forest,
                        t8_forest_t forest_from,
                        t8_locidx_t which_tree,
                        t8_locidx_t lelement_id,
                        t8_eclass_scheme_c * ts,
                        int is_family,
                        int num_elements, 
                        t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  
  /* remove every element outside an central ball of radius 0.65 */
  //double center[3] = {0.5, 0.5, 0.5};
  //dist = t8_vec_dist(center, centroid);
  //if (dist > 0.5) {
  //  return -2;
  //}

  for (int i = 0; i < 6; i++) {
    dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
    /* remove core of every ball */
    if (dist < 0.4) {
      return -2;
    }
    /* refine shell of every ball */
    if (dist < 0.5) {
      return 1;
    }
  }

  return 0;
}

int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          int is_family,
                          int num_elements, 
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist;
  
  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  
  for (int i = 0; i < 6; i++) {
    dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
    if (dist < 0.45) {
      return -2;
    }
  }

  return 0;
}

/* Coarse if at least one element of a family is within a given radius. */
int
t8_adapt_callback_coarse (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          int is_family,
                          int num_elements, 
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist;

  T8_ASSERT (adapt_data != NULL);
  
  if (is_family) {
    /* all Balls */
    for (int i = 0; i < 6; i++) {
      /* all member of family */
      for (int j = 0; j < num_elements; j++) {
        t8_forest_element_centroid (forest_from, which_tree, elements[j], centroid);
        dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
        /* if one family member satisfies the condition, coarse */
        if (dist < 0.5) {
          return -1;
        }
      }
    }
  }
  return 0;
}

int
t8_adapt_callback_refine (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          int is_family,
                          int num_elements, 
                          t8_element_t * elements[])
{
  const struct t8_adapt_data *adapt_data = (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  

  for (int i = 0; i < 6; i++) {
    dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
    if (dist < 0.5) {
      return 1;
    }
  }

  return 0;
}


void
t8_test_emelemts_remove (int cmesh_id)
{
  int                 level, min_level, max_level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest, forest_1, forest_2;
  t8_scheme_cxx_t    *scheme;

  struct t8_adapt_data adapt_data = {{{1.0, 0.5, 0.5},
                                      {0.5, 1.0, 0.5},
                                      {0.5, 0.5, 1.0},
                                      {0.0, 0.5, 0.5},
                                      {0.5, 0.0, 0.5},
                                      {0.5, 0.5, 0.0}}};

  scheme = t8_scheme_new_default_cxx ();
  /* Construct a cmesh */
  cmesh = t8_test_create_cmesh (cmesh_id);
  /* Compute the first level, such that no process is empty */
  min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* Use one level with empty processes */
  min_level = SC_MAX (min_level, 3);
  max_level = min_level + 1;
  
  for (level = min_level; level < max_level; level++) {
    t8_cmesh_ref (cmesh);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    //t8_forest_write_vtk (forest, "t8_example_test");
    T8_ASSERT (t8_forest_is_committed (forest));
    //t8_global_productionf("### %i\n", cmesh_id);
    forest_1 = t8_forest_new_adapt (forest  , t8_adapt_callback_init  , 0, 0, &adapt_data);
    forest_1 = t8_forest_new_adapt (forest_1, t8_adapt_callback_remove, 0, 0, &adapt_data);

    t8_forest_ref (forest_1);
    forest_2 = t8_forest_new_adapt (forest_1, t8_adapt_callback_coarse, 0, 0, &adapt_data);
    forest_2 = t8_forest_new_adapt (forest_2, t8_adapt_callback_refine, 0, 0, &adapt_data);
    forest_2 = t8_forest_new_adapt (forest_2, t8_adapt_callback_remove, 0, 0, &adapt_data);

    SC_CHECK_ABORT (t8_forest_is_equal(forest_1, forest_2),
                    "The forests are not equal");

    t8_scheme_cxx_ref (scheme);
    t8_forest_unref (&forest_1);
    t8_forest_unref (&forest_2);
  }
  t8_scheme_cxx_unref (&scheme);
  t8_cmesh_destroy (&cmesh);
}

void
test_cmesh_emelemts_remove_all ()
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    /* This if statement is necessary to make the test work by avoiding specific cmeshes which do not work yet for this test.
     * When the issues are gone, remove the if statement. */
    if (cmesh_id != 6 && cmesh_id != 89 && (cmesh_id < 237 || cmesh_id > 256)) {
      /* Skip all t8_test_create_new_bigmesh_cmesh */
      if (cmesh_id < 97 || cmesh_id > 256) {
        t8_test_emelemts_remove(cmesh_id);
      }
    }
  }
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
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  test_cmesh_emelemts_remove_all ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
