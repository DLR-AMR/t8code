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
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>

/** Construct the cmesh given the dimension of the examples. 
 * For 1D a line is constructed, for 2D a mesh consisting of 2 triangles
 * forming a quad and in 3D a cube made out of hexs, tets and prisms.
 * \param[in] dim           The dimension of the example. 1 <= \a dim <= 3. 
 * \param[in] do_partition  Option to partition the cmesh. If true (non-zero) the cmesh will be partitioned across the processes.  
 * \return                  The cmesh that is specified by the parameters*/
static t8_cmesh_t
t8_basic_create_cmesh (const int dim, const int do_partition)
{
  t8_cmesh_t          cmesh;
  switch (dim) {
  case 1:
    {
      cmesh =
        t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0,
                                do_partition, 0);
      break;
    }
  case 2:
    {
      cmesh =
        t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0,
                                do_partition, 0);
      break;
    }
  case 3:
    {
      cmesh =
        t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, do_partition, 0);
      break;
    }

  default:

    SC_ABORT_NOT_REACHED ();
  }
  return cmesh;
}

/* The adapt Callback used in this examples. Every element that has an odd
 * x-coordinate (in its reference element) wrt to its level is refined up to the maximal refinement level */
static int
t8_basic_adapt (t8_forest_t forest, t8_forest_t forest_from,
                t8_locidx_t which_tree, t8_locidx_t lelement_id,
                t8_eclass_scheme_c * ts, int num_elements,
                t8_element_t * elements[])
{
  int                 level, max_lvl, shift;
  double              coords[3] = { 0, 0, 0 };
  int                 scaled_x_coord;
  const int           rootlen = ts->t8_element_root_len (elements[0]);
  level = ts->t8_element_level (elements[0]);
  /* Check, if the element is not finer than the maximal refinement level */
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  /* Compute shift to take current level into account */
  max_lvl = ts->t8_element_maxlevel ();
  shift = max_lvl - level;
  ts->t8_element_vertex_reference_coords (elements[0], 0, coords);
  /* Scale x-Coord to integer coordinate */
  scaled_x_coord = coords[0] * rootlen;
  if ((scaled_x_coord >> shift) % 2 == 1) {
    return 1;
  }
  else {
    return 0;
  }
}

/** This function creates, adaptivly refines, partitions and optionally balances a
 * forest of the hypercube-mesh. In 3D a hybrid hypercube is created.
 * \param[in] dim         The dimension of the example
 * \param[in] do_balacne  Option to balance the mesh
*/
static void
t8_basic_hypercube (const int dim, const int do_balance)
{
  t8_forest_t         forest, forest_adapt, forest_partition;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ], cmesh_file[BUFSIZ];
  int                 mpirank, mpiret, adapt_level = 5;
  const int           uniform_lvl = 2;

  t8_global_productionf ("Contructing %i dimensional hypercube mesh \n", dim);
  /* Create and save the cmesh */
  cmesh = t8_basic_create_cmesh (dim, 0);
  snprintf (cmesh_file, BUFSIZ, "cmesh_basic_%i_dim", dim);
  t8_cmesh_save (cmesh, cmesh_file);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (vtuname, BUFSIZ, "cmesh_basic_%i__dim", dim);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
    t8_debugf ("Output to %s\n", vtuname);
  }
  else {
    t8_debugf ("Error in output\n");
  }
  /* Initialise the forest */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());

  /* Set uniform refinement level */
  t8_forest_set_level (forest, uniform_lvl);
  t8_forest_commit (forest);
  t8_debugf ("Successfully committed forest.\n");
  snprintf (vtuname, BUFSIZ, "forest_basic_%i_dim_uniform", dim);
  t8_forest_write_vtk (forest, vtuname);
  t8_debugf ("Output to %s\n", vtuname);

  t8_forest_init (&forest_adapt);
  /* Set maximum refinement level */
  t8_forest_set_user_data (forest_adapt, &adapt_level);
  /* Construct forest_adapt from forest */
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, 1);

  forest_partition = forest_adapt;
  /*Construct balanced and partitioned forest */
  t8_forest_set_partition (forest_partition, NULL, 0);
  if (do_balance) {
    t8_forest_set_balance (forest_partition, NULL, 0);
  }

  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);
  t8_forest_print_profile (forest_partition);

  snprintf (vtuname, BUFSIZ, "forest_hypercube_%i_dim_adapt", dim);
  t8_forest_write_vtk (forest_partition, vtuname);
  t8_debugf ("Output to %s\n", vtuname);
  t8_forest_unref (&forest_partition);
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
              "This program constructs an up to level 2 uniformly refined "
              "cubical mesh, which adaptivly refines up to level 5.\n"
              "The user can choose the dimension of the mesh "
              "and whether it should be balanced or not.\n\n%s\n", usage);
  if (sreturn >= BUFSIZ) {
    /* help message was truncated. */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'd', "dimension", &dim, 1,
                      "The dimension of the mesh.");
  sc_options_add_switch (opt, 'b', "balance", &do_balance,
                         "Additionally balance the forest.");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 1 <= dim && dim <= 3) {
    t8_basic_hypercube (dim, do_balance);
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
