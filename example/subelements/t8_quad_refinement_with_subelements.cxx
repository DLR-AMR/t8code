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
#include <time.h>

/* In this example, subelements are used to remove hanging faces from a refined 2D quad scheme. 
 * At first, a cmesh is created and a uniform version initialized. 
 * Then, it is adapted using some geometric refinement criterion and balanced (balance is automatically set in set_remove_hanging_faces if not done before). 
 * In the following step, the new subelement functions are used to identify elements that have hanging faces, 
 * which are then adapted once more, using transition cells with subelements in order to remove the hanging faces. 
 * The integer value "num_timesteps" determines the number of times, the mesh is adapted. During this process, the refinement criterion can be adjusted (for example an increasing radius) */

/* Define the data structure for the refinement criteria of this example (a circle with some user given midpoint and radius).
 * Elements whose anchor node is closer to the circle will be refined to a higher level than elements whose anchor node is farther away. */
typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Exemplary geometric refinement criteria */
static int
t8_advect_adapt_tests (t8_forest_t forest,
                       t8_forest_t forest_from,
                       t8_locidx_t which_tree,
                       t8_locidx_t lelement_id,
                       t8_eclass_scheme_c * ts,
                       int num_elements, t8_element_t * elements[])
{
  /* Get the minimum and maximum x-coordinate from the user data pointer of forest */
  t8_example_level_set_struct_t *data;
  int                 level;
  data = (t8_example_level_set_struct_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);

#if 0                           /* refine all lower right elements */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);

  if (coord[0] > coord[1]) {
    return 1;
  }
  return 0;
#endif

#if 1                           /* refinement diag */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);

  if (coord[0] == coord[1] && level < data->max_level) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* refinement every second element */
  if (lelement_id % 2 == 0) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* refinement all left elements */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);
  int                 len = ts->t8_element_root_len (elements[0]);
  if (coord[0] < len / 2) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* nested grid reifnement */
  int                 coord[3] = { };
  ts->t8_element_anchor (elements[0], coord);

  if (((coord[0] == 1 << 29 && coord[1] == 1 << 29) ||
       (coord[0] == 1 << 29 && coord[1] + (1 << (30 - level)) == 1 << 29) ||
       (coord[0] + (1 << (30 - level)) == 1 << 29 && coord[1] == 1 << 29) ||
       (coord[0] + (1 << (30 - level)) == 1 << 29
        && coord[1] + (1 << (30 - level)) == 1 << 29))
      && level < data->max_level) {
    return 1;
  }
  return 0;
#endif

#if 0                           /* random refinement */
  int                 r = rand () % 99; /* random number between 0 and 99 */

  if (level < data->max_level && r < 50) {
    return 1;
  }
  return 0;
#endif
}

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
static double
t8_refine_with_subelements (t8_eclass_t eclass, int initlevel, int adaptlevel)
{
  double adapt_time = 0;
  t8_productionf ("Into the t8_refine_with_subelements function\n");

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* refinement settings */
                      initlevel = 3;    /* initial uniform refinement level */
  int                 minlevel = initlevel;   /* lowest level allowed for coarsening */
                      adaptlevel = 3;   /* number of allowed adapt levels */
  int                 maxlevel = initlevel + adaptlevel;    /* highest level allowed for refining */

  int                 refine_recursive = 1;
  int                 do_exemplary_refinement = 0;

  /* cmesh settings (only one of the following suggestions should be one, the others 0) */
  int                 single_tree = 1;
  int                 multiple_tree = 0, num_x_trees = 2, num_y_trees = 1;
  int                 hybrid_cmesh = 0;

  /* adaptation setting */
  int                 do_balance = 0;
  int                 do_transition = 1;

  /* MPI settings */
  int                 do_ghost = 1;
  int                 ghost_version = 3;

  /* timestep settings */
  int                 do_different_refinements = 0;   /* to change the refinement during multiple num_timesteps */
  int                 num_timesteps = 1;    /* Number of times, the mesh is refined */
  double              delta = 0.2;   /* The value, the radius increases after each timestep */

  /* VTK settings */
  int                 do_vtk = 1;   /* print results */
  int                 do_vtk_cmesh = 0;

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
  if (do_vtk && do_vtk_cmesh) {
    snprintf (filename, BUFSIZ, "forest_cmesh_%s", t8_eclass_to_string[eclass]);
    t8_cmesh_vtk_write_file (cmesh, filename, 1);

    /* print uniform mesh file */
    snprintf (filename, BUFSIZ, "forest_uniform_%s",
              t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest, filename);
  }

  /* user-data */
  t8_example_level_set_struct_t ls_data;
  t8_basic_sphere_data_t sdata;

  /* Midpoint and radius of a sphere */
  /* shift the midpoiunt of the circle by (shift_x,shift_y) to ensure midpoints on corners of the uniform mesh */
  // int  shift_x = 0;      /* shift_x, shift_y should be smaler than 2^minlevel / 2 such that midpoint stays in the quadrilateral tree */
  // int  shift_y = 0;
  sdata.mid_point[0] = 0;    // 1.0 / 2.0 + shift_x * 1.0/(1 << (minlevel));
  sdata.mid_point[1] = 0;    // 1.0 / 2.0 + shift_y * 1.0/(1 << (minlevel)); 
  sdata.mid_point[2] = 0;
  sdata.radius = 0.4;

  /* refinement parameter */
  ls_data.band_width = 1;
  ls_data.L = t8_basic_level_set_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* TODO: the time_forest examples use several refine_rounds within each time-loop. 
   * This might resolve the problem with appearing artifacts during adaptation with multiple num_timesteps.
   * This does not seem to be a problem of subelements but of the timestep scheme implemented in this example. */

  /* Adapting the mesh for different num_timesteps */
  int i;
  for (i = 0; i < num_timesteps; i++) {

    /* change the refinement */
    if (i == do_different_refinements - 1) {
      do_transition = 1;
      adaptlevel = 1;
      maxlevel = initlevel + adaptlevel;
      do_exemplary_refinement = 1;
    }

    t8_productionf
      ("This is t8_refine_with_subelements. Into timestep %i of %i\n", i + 1,
       num_timesteps);

    t8_forest_init (&forest_adapt);

    /* Adapt the mesh according to the user data */
    t8_forest_set_user_data (forest_adapt, &ls_data);

    if (do_exemplary_refinement) {
      time_t              t;    /* we might use a random refinement and set the randomizer seed */
      srand ((unsigned) time (&t));
      t8_forest_set_adapt (forest_adapt, forest, t8_advect_adapt_tests,
                           refine_recursive);
    }
    else {
      t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set,
                           refine_recursive);
    }

    if (do_balance) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (do_transition) {
      /* Analogue to the other set-functions, this function adds subelements to the from_method. 
       * The forest will therefore use subelements while adapting in order to remove hanging faces from the mesh. */
      t8_forest_set_remove_hanging_faces (forest_adapt, NULL);
      ghost_version = 1;
    }

    if (do_ghost) {
      /* set ghosts after adaptation/balancing/transitioning */
      t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES, ghost_version);
    }

    adapt_time -= sc_MPI_Wtime (); 
    t8_forest_commit (forest_adapt);
    adapt_time += sc_MPI_Wtime ();
    t8_productionf ("Commit runtime: %f\n", adapt_time);

    if (do_vtk) { /* print to vtk */
      snprintf (filename, BUFSIZ, "forest_adapt_no_hanging_nodes_timestep%i_%s",
                i, t8_eclass_to_string[eclass]);
      t8_forest_write_vtk (forest_adapt, filename);
    }

    /* Set forest to the partitioned forest, so it gets adapted
     * in the next timestep. */
    forest = forest_adapt;

    /* Increasing the radius of the sphere for the next timestep */
    sdata.radius += delta;

    t8_productionf
      ("This is t8_refine_with_subelements. Done timestep %i of %i\n", i + 1,
       num_timesteps);
  }                             /* end of time-loop */
  t8_forest_unref (&forest_adapt);

  return adapt_time;
}                               /* end of function */

int
main (int argc, char **argv)
{
#if 1 /* standard: only one run */
  int                 mpiret;
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  double commit_time = t8_refine_with_subelements (T8_ECLASS_QUAD, 0, 0);
  t8_productionf ("Commit time: %f\n", commit_time);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);


#else /* Can be used for benchmarks. Executes the adaptation several times and with possibly differing init and adapt levels and returns the mean runtime */
  int init_min = 9, init_max = 9, adapt_min = 2, adapt_max = 5;
  int repeat, init_count, adapt_count, repitition = 5;
  double commit_time[(adapt_max - adapt_min + 1) * (init_max - init_min + 1) * repitition] = {0};
  int h = 0;

  for (init_count = init_min; init_count <= init_max; init_count++) {
    for (adapt_count = adapt_min; adapt_count <= adapt_max; adapt_count++) {
      for (repeat = 0; repeat < repitition; repeat++) {
        int                 mpiret;
        mpiret = sc_MPI_Init (&argc, &argv);
        SC_CHECK_MPI (mpiret);

        sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
        t8_init (SC_LP_DEFAULT);

        /* At the moment, subelements are only implemented for the quad scheme */
        t8_productionf ("init: %i, adapt: %i\n", init_count, adapt_count);
        commit_time[h] = t8_refine_with_subelements (T8_ECLASS_QUAD, init_count, adapt_count);
        h++;
        sc_finalize ();

        mpiret = sc_MPI_Finalize ();
        SC_CHECK_MPI (mpiret);
      }
    }
  }

  int i,j;
  int count = 0;
  for (i = 0; i < (adapt_max - adapt_min + 1) * (init_max - init_min + 1); i++) {
    double commit_time_sum = 0;
    for (j = 0; j < repitition; j++) {
      // t8_productionf ("commit time: %f\n", commit_time[count + j]);
      commit_time_sum += commit_time[count + j];
    }
    count += repitition;
    t8_productionf ("init: %i  time commit sum: %f  commit time mean: %f\n", i/(adapt_max - adapt_min + 1) + init_min, commit_time_sum, commit_time_sum/repitition);
  }
#endif
  
  return 0;
}

