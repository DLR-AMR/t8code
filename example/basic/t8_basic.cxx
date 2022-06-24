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
#include <sc_refcount.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_version.h>

/** Construct the cmesh given the dimension of the examples. 
 * For 1D a line is constructed, for 2D a mesh consisting of 2 triangles
 * forming a quad and in 3D a cube made out of hexs, tets and prisms.
 * \param[in] dim           The dimension of the example. 1 <= \a dim <= 3.  
 * \return                  The cmesh that is specified by the parameters*/
static t8_cmesh_t
t8_basic_create_cmesh (const int dim)
{
  t8_cmesh_t          cmesh;
  switch (dim) {
  case 1:
    {
      cmesh =
        t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 0);
      break;
    }
  case 2:
    {
      cmesh =
        t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0,
                                0);
      break;
    }
  case 3:
    {
      cmesh = t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);
      break;
    }

  default:

    SC_ABORT_NOT_REACHED ();
  }
  return cmesh;
}

/* The adapt Callback used in this example. We compute the center of the first given element
 * If the x-coordinate of the center lays between 0.25 and 0.75 we refine. If it is smaller than 0.25
 * we coarsen the element. Otherwise, we do nothing. */
static int
t8_basic_adapt (t8_forest_t forest, t8_forest_t forest_from,
                t8_locidx_t which_tree, t8_locidx_t lelement_id,
                t8_eclass_scheme_c *ts, const int is_family,
                const int num_elements, t8_element_t *elements[])
{
  int                 level, i;
  double              coords[3] = { 0, 0, 0 };
  const int           maxlevel = 5;
  level = ts->t8_element_level (elements[0]);
  /* Check, if the element is not finer than the maximal refinement level */
  if (level >= maxlevel) {
    /* We do not refine if the maximum level is reached */
    return 0;
  }
  t8_forest_element_centroid (forest_from, which_tree, elements[0], coords);
  /* If the x-coordinate lays between 0.25 and 0.75, we refine the element. */
  if (0.25 <= coords[0] && coords[0] <= 0.75) {
    return 1;
  }
  /* If the x-coordinate is smaller than 0.25 we coarsen the element. */
  else if (coords[0] < 0.25 && is_family) {
    /* Check that every element of the family should be coarsen. */
    for (i = 1; i < num_elements; i++) {
      t8_forest_element_centroid (forest_from, which_tree, elements[i],
                                  coords);
      if (coords[0] < 0.25) {
        return -1;
      }
    }
    return 0;
  }
  /* Do not refine or coarsen the element. */
  else {
    return 0;
  }
}

/** This function creates, adaptively refines, partitions and optionally balances a
 * forest of the hypercube-mesh. For 1D a line is constructed, for 2D a mesh consisting of 2 triangles
 * forming a quad and in 3D a hybrid hypercube is created.
 * \param[in] dim         The dimension of the example
 * \param[in] do_balance  Option to balance the cmesh. If true (non-zero) a 2:1 balance will be established.  
*/
static void
t8_basic_hypercube (const int dim, const int do_balance)
{
  t8_forest_t         forest, forest_adapt;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  int                 mpirank, mpiret;
  const int           uniform_lvl = 2;

  t8_global_productionf ("Constructing %i dimensional hypercube mesh \n",
                         dim);
  /* Create the cmesh */
  cmesh = t8_basic_create_cmesh (dim);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (vtuname, BUFSIZ, "cmesh_basic_%i__dim", dim);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
    t8_global_productionf ("Output to %s\n", vtuname);
  }
  else {
    t8_errorf ("Error when trying to write cmesh to %s.\n", vtuname);
  }
  /* Initialize the forest */
  t8_forest_init (&forest);
  /* Initialize the cmesh of the forest */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  /* Set the scheme of the forest. In this case, the default schemes are used */
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());

  /* Set uniform refinement level */
  t8_forest_set_level (forest, uniform_lvl);
  /* Generate the uniform refinement. */
  t8_forest_commit (forest);
  t8_global_productionf ("Successfully committed forest.\n");

  /*  Write ouptut to a vtk file */
  snprintf (vtuname, BUFSIZ, "forest_basic_%id_uniform", dim);
  t8_forest_write_vtk (forest, vtuname);
  t8_global_productionf ("Output to %s\n", vtuname);

  /* Initialize a second forest, that will be the adaptively refined forest */
  t8_forest_init (&forest_adapt);
  /* Construct forest_adapt from forest */
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, 0);
  /* Construct balanced and partitioned forest */
  t8_forest_set_partition (forest_adapt, forest, 0);
  if (do_balance) {
    t8_forest_set_balance (forest_adapt, forest, 0);
  }
  /* Generate the adapted, partitioned (and balanced) forest. */
  t8_forest_commit (forest_adapt);
  t8_global_productionf ("Successfully committed forest.\n");

  /* Write output to a vtk file */
  snprintf (vtuname, BUFSIZ, "forest_hypercube_%id_adapt", dim);
  t8_forest_write_vtk (forest_adapt, vtuname);
  t8_global_productionf ("Output to %s\n", vtuname);

  /* Destroy the forest */
  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 dim, do_balance;
  int                 parsed, helpme;
  int                 sreturn;
  int                 print_version = 0;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* brief help message */
  sreturn = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                      "for a brief overview of all options.",
                      basename (argv[0]), basename (argv[0]));

  if (sreturn >= BUFSIZ) {
    /* Usage string was truncated. */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string to '%s'\n", usage);
  }

  /* long help message */
  sreturn =
    snprintf (help, BUFSIZ,
              "This program constructs a level 2 uniformly refined "
              "cubical mesh, which is then adaptively refined with "
              "level 5 as the maximal level.\n"
              "The user can choose the dimension of the mesh "
              "and whether it should be 2:1 balanced or not.\n\n%s\n", usage);
  if (sreturn >= BUFSIZ) {
    /* help message was truncated. */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'd', "dimension", &dim, 1,
                      "The dimension of the mesh. Choose 1 <= dimension <= 3.");
  sc_options_add_switch (opt, 'b', "balance", &do_balance,
                         "Additionally 2:1 balance the forest.");

  sc_options_add_switch (opt, 'v', "version", &print_version,
                         "Print the version number of t8code and exit.");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  /* Print usage */
  sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
  }
  else if (!print_version && parsed >= 0 && 1 <= dim && dim <= 3) {
    t8_basic_hypercube (dim, do_balance);
  }
  else if (print_version) {
    t8_global_productionf ("This is t8code version '%s'\n",
                           t8_get_package_string ());
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
