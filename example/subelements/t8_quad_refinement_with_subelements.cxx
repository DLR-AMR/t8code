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

#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include <t8_schemes/t8_quads_w_rectangular_subelements/t8_subelements_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>

/* In this example, subelements are used to remove hanging nodes from a refined 2D quad scheme. 
 * At first, a cmesh is adapted using the standard 2D refinement scheme AND balance (where balance is 
 * automatically set in set_remove_hanging_faces if not done before). 
 * In the following step, the new subelement functions are used to identify elements that have hanging faces, 
 * which are then adapted once more, using transition cells with subelements in order to remove the hanging faces. 
 * The integer value "timesteps" determines the number of times, the mesh is adapted. During this process, 
 * a circle spreads within the mesh and elements near this circle will be refined to a higher level. */

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

/* Recommended settings for the refinement test with subelements: 
 *   initlevel = 1
 *   minlevel = initlevel 
 *   maxlevel = 3
 *   do_subelements = 1 */
static void
t8_refine_with_subelements (t8_eclass_t eclass)
{
  t8_productionf ("Into the t8_refine_with_subelements function");

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* refinement settings */
  int                 initlevel = 2;    /* initial uniform refinement level */
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening */
  int                 maxlevel = 6;     /* highest level allowed for refining */

  /* cmesh settings (only one of the following suggestions should be one, the others 0) */
  int                 single_tree = 0;
  int                 multiple_tree = 1, num_x_trees = 2, num_y_trees = 2;
  int                 hybrid_cmesh = 0;

  /* adaptation setting */
  int                 do_balance = 0;
  int                 do_subelements = 1;

  /* timestep settings */
  int                 timesteps = 1;    /* Number of times, the mesh is refined */
  double              delta = 0.3;      /* The value, the radius increases after each timestep */
  int                 i;

  /* initializing the forest */
  t8_forest_init (&forest);

  /* building the cmesh, using the initlevel */
  if (single_tree) {            /* single quad cmesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }

  if (multiple_tree) {          /* p4est_connectivity_new_brick (num_x_trees, num_y_trees, 0, 0) -> cmesh of (num_x_trees x num_y_trees) many quads */
    p4est_connectivity_t *brick =
      p4est_connectivity_new_brick (num_x_trees, num_y_trees, 0, 0);
    cmesh = t8_cmesh_new_from_p4est (brick, sc_MPI_COMM_WORLD, 0);
    p4est_connectivity_destroy (brick);
  }

  if (hybrid_cmesh) {           /* TODO: this does not work at the moment */
    cmesh = t8_cmesh_new_hypercube_hybrid (2, sc_MPI_COMM_WORLD, 0, 0);
  }

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_subelement_cxx ());
  t8_forest_set_level (forest, initlevel);

  t8_forest_commit (forest);

  /* print cmesh file */
  snprintf (filename, BUFSIZ, "forest_cmesh_%s",
            t8_eclass_to_string[eclass]);
  t8_cmesh_vtk_write_file (cmesh, filename, 1);

  /* print uniform mesh file */
  snprintf (filename, BUFSIZ, "forest_uniform_%s",
            t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest, filename);

  /* user-data */
  t8_example_level_set_struct_t ls_data;
  t8_basic_sphere_data_t sdata;

  /* midpoint and radius of a sphere 
   * TODO: check if, for symmetry, the midpoint should be on an elements corner */
  /* shift the midpoiunt of the circle by (shift_x,shift_y) to ensure midpoints on corners of the uniform mesh */
  int shift_x = 0; /* shift_x should be smaler than 2^minlevel / 2 such that midpoint stays in the quadrilateral tree */
  int shift_y = 0; 
  sdata.mid_point[0] = 0; //1.0 / 2.0 + shift_x * 1.0/(1 << (minlevel));
  sdata.mid_point[1] = 0; //1.0 / 2.0 + shift_y * 1.0/(1 << (minlevel)); 
  sdata.mid_point[2] = 0;
  sdata.radius = 1.2;

  /* refinement parameter */
  ls_data.band_width = 1;
  ls_data.L = t8_basic_level_set_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* TODO: the time_forest examples use several refine_rounds within each time-loop. 
   * This might resolve the problem with appearing artifacts during adaptation with multiple timesteps.
   * This does not seem to be a problem of subelements but of the timestep scheme implemented in this example. */

  /* Adapting the mesh for different timesteps */
  for (i = 0; i < timesteps; i++) {

    t8_productionf
      ("This is t8_refine_with_subelements. Into timestep %i of %i\n", i + 1,
       timesteps);

    t8_forest_init (&forest_adapt);

    /* Adapt the mesh according to the user data */
    t8_forest_set_user_data (forest_adapt, &ls_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);

    if (do_balance) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (do_subelements) {
      /* Analogue to the other set-functions, this function adds subelements to the from_method. 
       * The forest will therefore use subelements while adapting in order to remove hanging faces from the mesh. */
      t8_forest_set_remove_hanging_faces (forest_adapt, NULL);
    }

    t8_forest_commit (forest_adapt);

    /* print to vtk */
    snprintf (filename, BUFSIZ, "forest_adapt_no_hanging_nodes_timestep%i_%s",
              i, t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest_adapt, filename);

    /* Set forest to the partitioned forest, so it gets adapted
     * in the next timestep. */
    forest = forest_adapt;

    /* Increasing the radius of the sphere for the next timestep */
    sdata.radius += delta;

    t8_productionf
      ("This is t8_refine_with_subelements. Done timestep %i of %i\n", i + 1,
       timesteps);
  }                             /* end of time-loop */
  t8_forest_unref (&forest_adapt);
}                               /* end of function */

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
