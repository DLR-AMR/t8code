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

#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>

/* In this example, subelements are used to remove hanging nodes from a refined 2D quad scheme. 
 * At first, a cmesh is adapted using the standard 2D refinement scheme AND balance. 
 * In the following step, the new subelement functions are used to identify elements that have hanging nodes, 
 * which are then adapted once more, using transition cells with subelements in order to remove the hanging nodes. */

/* Define the data structure for the refinement criteria of this example (a circle with some user given midpoint and radius).
 * Elements whose anchor node is closer to the circle will be refined to a higher level than elements whose anchor node is farther away. */
typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Compute the distance to a sphere around a mid_point with given radius. */
static double
t8_basic_level_set_sphere (const double x[3], double t, void *data)
{
  t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
  double             *M = sdata->mid_point;

  return t8_vec_dist (M, x) - sdata->radius;
}

static void
t8_refine_with_subelements (t8_eclass_t eclass)
{
  t8_productionf ("Into the t8_refine_with_subelements function");

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  int                 initlevel = 2;    /* level of the cmesh */
  int                 minlevel = 1;     /* lowest level allowed for coarsening */
  int                 maxlevel = 3;     /* highest level allowed for refining */

  t8_forest_init (&forest);
  t8_forest_init (&forest_adapt);

  /* building the cmesh, using the initlevel */
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());
  t8_forest_set_level (forest, initlevel);

  t8_forest_commit (forest);

  snprintf (filename, BUFSIZ, "forest_uniform_%s",
            t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest, filename);

  /* adapting the cmesh, using user-data, up to a given max and min level */
  t8_example_level_set_struct_t ls_data;
  t8_basic_sphere_data_t sdata;

  sdata.mid_point[0] = 0.5;
  sdata.mid_point[1] = 0.5;
  sdata.mid_point[2] = 0;
  sdata.radius = 0.3;

  ls_data.band_width = 1;
  ls_data.L = t8_basic_level_set_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* Adapt the mesh according to the user data. At the moment, balancing the adapted mesh afterwards is important for 
   * the subelement functions to work. */
  t8_forest_set_user_data (forest_adapt, &ls_data);
  t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);
  t8_forest_set_balance (forest_adapt, NULL, 0);

  /* This function is the point of entry for the subelements.
   * Analogue to the other set-functions, we add subelements to the from_method and will therefore 
   * use the corresponding functions later during the adaption. */
  t8_forest_set_subelements (forest_adapt, NULL);

  t8_forest_commit (forest_adapt);

  snprintf (filename, BUFSIZ, "forest_adapt_no_hanging_nodes_%s",
            t8_eclass_to_string[eclass]);
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

  /* At the moment, subelements are only implemented for the quad scheme */
  t8_refine_with_subelements (T8_ECLASS_QUAD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
