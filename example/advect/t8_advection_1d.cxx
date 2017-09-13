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
  t8_scalar_function_1d_fn u; /**< Fluid field */
  t8_forest_t         forest; /**< The forest in use */
  sc_array            element_data; /**< Array of type t8_advect_element_data_t of length
                              num_local_elements + num_ghosts */
  sc_MPI_Comm         comm; /**< MPI communicator used */
  int                 level; /**< Initial refinement level */
  int                 maxlevel; /**< Maximum refinement level */
} t8_advect_problem_t;

typedef struct
{
  double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
  double              phi; /**< Value of solution at midpoint */
} t8_advect_element_data_t;

t8_advect_problem_t *
t8_advect_problem_init (t8_scalar_function_1d_fn u, int level, int maxlevel,
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
  problem->level = level;
  problem->maxlevel = maxlevel;
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

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 parsed, helpme;
  t8_advect_problem_t *problem;

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

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n ", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level) {
    /* Computation */
    problem =
      t8_advect_problem_init (constant_one, level, level + 4,
                              sc_MPI_COMM_WORLD);
    t8_advect_problem_destroy (&problem);
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
