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

#include <sc_options.h>
#include <t8_default_cxx.hxx>
#include <t8_forest.h>
#include <example/common/t8_example_common.h>

typedef struct
{
  t8_scalar_function_3d_fn u; /**< Fluid field */
  t8_scalar_function_3d_fn phi_0; /**< Initial condition for phi */
  t8_forest_t         forest; /**< The forest in use */
  sc_array            element_data; /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  sc_MPI_Comm         comm; /**< MPI communicator used */
  double              t; /**< Current simulation time */
  double              T; /**< End time */
  double              delta_t; /**< Current time step */
  int                 level; /**< Initial refinement level */
  int                 maxlevel; /**< Maximum refinement level */
} t8_advect_problem_t;

typedef struct
{
  double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
  double              delta_x; /**< Width of this element */
  double              phi; /**< Value of solution at midpoint */
} t8_advect_element_data_t;

double
t8_advect_flux_lax_friedrich (t8_advect_problem_t * problem,
                              t8_advect_element_data_t * el_data_plus,
                              t8_advect_element_data_t * el_data_minus)
{
  double              alpha = 1;        /* TODO: Choose alpha according to a reasonable criterion */
  double              x_j_half[3];
  int                 idim;
  double              u_at_x_j_half;
  double              phi_sum, phi_diff;

  /*
   *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
   *       x_j     x_j+1
   *          x_j_half
   */
  /* Compute x_j_half */
  for (idim = 0; idim < 3; idim++) {
    x_j_half[idim] =
      (el_data_plus->midpoint[idim] + el_data_minus->midpoint[idim]) / 2;
  }

  /* Compute u at the interval boundary. */
  u_at_x_j_half = problem->u (x_j_half, problem->t);

  /* Compute the sum of both phi values */
  phi_sum = el_data_minus->phi + el_data_plus->phi;
  /* Compute the difference of both */
  phi_diff = el_data_plus->phi - el_data_minus->phi;
  return .5 * (u_at_x_j_half * phi_sum - alpha * phi_diff);
}

void
t8_advect_advance_element (t8_advect_problem_t * problem,
                           t8_advect_element_data_t * elem,
                           double flux_left, double flux_right)
{
  /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
  elem->phi = (problem->delta_t / elem->delta_x) * (flux_left - flux_right)
    + elem->phi;
}

t8_advect_problem_t *
t8_advect_problem_init (t8_scalar_function_3d_fn u,
                        t8_scalar_function_3d_fn phi_0, int level,
                        int maxlevel, double T, double delta_t,
                        sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_advect_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;

  /* Construct new hypercube cmesh (unit interval) */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, comm, 0, 0);

  /* allocate problem */
  problem = T8_ALLOC (t8_advect_problem_t, 1);
  /* Fill problem parameters */
  problem->u = u;
  problem->phi_0 = phi_0;
  problem->level = level;
  problem->maxlevel = maxlevel;
  problem->t = 0;
  problem->T = T;
  problem->delta_t = delta_t;
  problem->comm = comm;

  /* Contruct uniform forest with ghosts */
  default_scheme = t8_scheme_new_default_cxx ();
  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  /* Initialize the element array with num_local_elements + num_ghosts entries. */
  sc_array_init_size (&problem->element_data,
                      sizeof (t8_advect_element_data_t),
                      t8_forest_get_num_element (problem->forest) +
                      t8_forest_get_num_ghosts (problem->forest));
  return problem;
}

void
t8_advect_compute_element_data (t8_advect_problem_t * problem,
                                t8_advect_element_data_t * elem_data,
                                t8_element_t * element, t8_locidx_t ltreeid,
                                t8_eclass_scheme_c * ts,
                                double *tree_vertices)
{
  /* Compute the midpoint coordinates of element */
  t8_forest_element_centroid (problem->forest, ltreeid, element,
                              tree_vertices, elem_data->midpoint);
  /* Compute the length of this element */
  elem_data->delta_x =
    1. / (((uint64_t) 1) << ts->t8_element_level (element));
}

void
t8_advect_problem_init_elements (t8_advect_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element;
  t8_advect_element_data_t *elem_data;
  t8_eclass_scheme_c *ts;
  double             *tree_vertices;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    ts =
      t8_forest_get_eclass_scheme (problem->forest,
                                   t8_forest_get_tree_class (problem->forest,
                                                             itree));
    num_elems_in_tree =
      t8_forest_get_tree_num_elements (problem->forest, itree);
    /* TODO: A forest get tree vertices function */
    tree_vertices =
      t8_cmesh_get_tree_vertices (t8_forest_get_cmesh (problem->forest),
                                  t8_forest_ltreeid_to_cmesh_ltreeid
                                  (problem->forest, itree));
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element =
        t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_data = (t8_advect_element_data_t *)
        t8_sc_array_index_locidx (&problem->element_data, idata);
      /* Initialize the element's midpoint and length */
      t8_advect_compute_element_data (problem, elem_data, element, itree,
                                      ts, tree_vertices);
      /* Set the initial condition */
      elem_data->phi = problem->phi_0 (elem_data->midpoint, 0);
    }
  }
}

void
t8_advect_problem_destroy (t8_advect_problem_t ** pproblem)
{
  t8_advect_problem_t *problem;
  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }
  /* Unref the forest */
  t8_forest_unref (&problem->forest);
  /* Free the element array */
  sc_array_reset (&problem->element_data);
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

void
t8_advect_solve (t8_scalar_function_3d_fn u,
                 t8_scalar_function_3d_fn phi_0, int level, int maxlevel,
                 double T, double delta_t, sc_MPI_Comm comm)
{
  t8_advect_problem_t *problem;
  problem =
    t8_advect_problem_init (u, phi_0, level, maxlevel, T, delta_t, comm);
  t8_advect_problem_init_elements (problem);
  t8_advect_problem_destroy (&problem);
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 parsed, helpme;
  double              T, delta_t;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "This program solves the 1D advection equation on "
            "the interval [0,1].\nThe user can choose the initial uniform "
            "refinement level.\n\n%s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_double (opt, 'T', "end-time", &T, 1,
                         "The duration of the simulation. Default: 1");
  sc_options_add_double (opt, 't', "delta-t", &delta_t, 0.1,
                         "The length of ont time-step. Default: 0.1");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n ", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level) {
    /* Computation */
    t8_advect_solve (constant_one, step_function, level, level + 4, T,
                     delta_t, sc_MPI_COMM_WORLD);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
