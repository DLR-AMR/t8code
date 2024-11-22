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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_profiling.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <t8_cmesh/t8_cmesh_examples.h>

static int
t8_basic_adapt_refine_type (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                            t8_locidx_t lelement_id, t8_scheme *ts, const int is_family, const int num_elements,
                            t8_element_t *elements[])
{
  int level;
  int type;

  T8_ASSERT (!is_family || num_elements == ts->element_get_num_children (tree_class, elements[0]));

  level = ts->element_get_level (tree_class, elements[0]);
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  /* get the type of the current element */
  type = ((t8_dprism_t *) elements[0])->tri.type;
  /* refine type 0 */
  if (type == 0) {
    return 1;
  }
  return 0;
}

static int
t8_basic_adapt_refine_tet (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                           t8_locidx_t lelement_id, t8_scheme *ts, const int is_family, const int num_elements,
                           t8_element_t *elements[])
{
  int level;
  int type;

  T8_ASSERT (!is_family || num_elements == ts->element_get_num_children (tree_class, elements[0]));

  level = ts->element_get_level (tree_class, elements[0]);
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  /* get the type of the current element */
  type = ((t8_dtet_t *) elements[0])->type;
  /* refine type 0 */
  if (type == 0 || type == 2 || type == 4) {
    return 1;
  }
  return 0;
}

static void
t8_time_refine (int start_level, int end_level, int create_forest, int cube, int adapt, int do_balance,
                t8_eclass_t eclass)
{
  t8_forest_t forest, forest_adapt, forest_partition;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];
  char vtuname[BUFSIZ];

  T8_ASSERT (eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_TET);
  t8_forest_init (&forest);

  if (cube == 0) {
    t8_forest_set_cmesh (forest, t8_cmesh_new_bigmesh (eclass, 512, sc_MPI_COMM_WORLD), sc_MPI_COMM_WORLD);
  }
  else {
    t8_forest_set_cmesh (forest, t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0), sc_MPI_COMM_WORLD);
  }
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, start_level);
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  t8_forest_commit (forest);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");
  if (cube == 1) {
    snprintf (vtuname, BUFSIZ, "forest_hypercube_%s", t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest, vtuname);
    t8_debugf ("Output to %s\n", vtuname);
  }

  if (adapt == 1) {
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, &end_level);

    t8_forest_set_profiling (forest_adapt, 1);
    if (eclass == T8_ECLASS_PRISM) {
      t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt_refine_type, 1);
    }
    else {
      t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt_refine_tet, 1);
    }
    forest_partition = forest_adapt;
    /* partition the adapted forest */
    t8_forest_set_partition (forest_partition, NULL, 0);
    /* enable profiling for the partitioned forest */
    t8_forest_set_profiling (forest_partition, 1);
    /* if desired do balance */
    if (do_balance) {
      t8_forest_set_balance (forest_partition, NULL, 0);
    }
    t8_forest_commit (forest_partition);
    t8_forest_print_profile (forest_partition);
    if (cube == 1) {
      snprintf (vtuname, BUFSIZ, "forest_hypercube_adapt_%s", t8_eclass_to_string[eclass]);
      t8_forest_write_vtk (forest_partition, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
    }
    t8_forest_unref (&forest_partition);
  }
  else {
    t8_forest_print_profile (forest);
    t8_forest_unref (&forest);
  }
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_options_t *opt;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int create_forest;
  int start_level = 0, end_level = 1, cube = 0, adapt = 0, do_balance = 0, eclass_int;
  int parsed, helpme;
  int sreturnA, sreturnB;

  /* brief help message */
  sreturnA = snprintf (usage, BUFSIZ,
                       "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                       "for a brief overview of all options.",
                       basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturnB = snprintf (help, BUFSIZ,
                       "This program constructs a prism mesh of 512 prisms. "
                       "\nThe user can choose the initial refinement level and the final\n"
                       "refinement level of the mesh. If not set, the initial level is 0,\n"
                       "the final level is 1. The program has no visual output, if desired,\n"
                       "the user can switch to a hpyercube mesh.\n\n%s\n",
                       usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 's', "slevel", &start_level, 0, "initial refine level");
  sc_options_add_int (opt, 'f', "flevel", &end_level, 0,
                      "Final refine level: greater or equal to initial refine level");
  sc_options_add_int (opt, 'a', "adapt", &adapt, 0, "adapt = 1 -> adaptive refining is used");
  sc_options_add_switch (opt, 'b', "balance", &do_balance, "Establish a 2:1 balance in the forest.");
  sc_options_add_int (opt, 'c', "cube", &cube, 0, "cube = 1 -> use the hypercube mesh and visual output.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 6,
                      "This option specifies the type of elements to use.\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism");

  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (end_level < start_level) {
    t8_debugf ("Wrong usage of end and start level, end level set to start level + 1\n");
    end_level = start_level;
  }
  if (end_level == start_level) {
    adapt = 0;
    t8_debugf ("End_level = start_level, adapt is set to zero\n");
  }
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= start_level && start_level <= end_level && (eclass_int == 5 || eclass_int == 6)) {
    create_forest = 1;
    t8_time_refine (start_level, end_level, create_forest, cube, adapt, do_balance, (t8_eclass_t) eclass_int);
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
