/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <sc_refcount.h>
#include <sc_shmem.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <t8_schemes/t8_new_feature/t8_subelements_cxx.hxx>
/* to validate the results with the original quad scheme */
// #include <t8_schemes/t8_default_cxx.hxx> 
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>

typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Compute the distance to a sphere arount a mid_point with given radius. */
static double
t8_basic_level_set_sphere (const double x[3], double t, void *data)
{
  t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
  double             *M = sdata->mid_point;

  return t8_vec_dist (M, x) - sdata->radius;
}

static void
t8_basic_refine_test (t8_eclass_t eclass)
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  int                 initlevel = 0;                
  int                 minlevel = initlevel; 
  int                 maxlevel = initlevel + 1;

  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);

  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());
  /* to validate the results with the original quad scheme */
  // t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  
  t8_forest_set_level (forest, initlevel);
  t8_forest_commit (forest);

  /* Output to vtk */
  snprintf (filename, BUFSIZ, "e_s_forest_uniform_%s",
            t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest, filename);

    t8_example_level_set_struct_t ls_data;
    t8_basic_sphere_data_t sdata;

    sdata.mid_point[0] = 0.5;
    sdata.mid_point[1] = 0.5;
    sdata.mid_point[2] = 0;
    sdata.radius = 0.3;

    ls_data.band_width = 0.5;
    ls_data.L = t8_basic_level_set_sphere;
    ls_data.min_level = minlevel;
    ls_data.max_level = maxlevel;
    ls_data.udata = &sdata;
    t8_forest_set_user_data (forest_adapt, &ls_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);
    // t8_forest_set_balance (forest_adapt, NULL, 0);

  t8_forest_commit (forest_adapt);

  /* Output to vtk */
  snprintf (filename, BUFSIZ, "e_s_forest_adapt_%s", t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest_adapt, filename);

  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* use T8_ECLASS_QUAD (our quad_w_sub scheme) */
  t8_basic_refine_test (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
