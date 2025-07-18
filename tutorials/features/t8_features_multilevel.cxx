/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2025 the developers

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

#include <t8.h>                                       /* General t8code header, always include this. */
#include <sc_options.h>                               /* CLI parser */
#include <t8_cmesh.hxx>                               /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>               /* example cmeshes */
#include <t8_forest/t8_forest_general.h>              /* forest definition and basic interface. */
#include <t8_schemes/t8_standalone/t8_standalone.hxx> /* standalone refinement scheme. */
#include <string>                                     /* std::string */
#include <array>                                      /* std::array */

/* Refine a tree with quadrilateral elements.
 * At every second adaptcall (level is even) remove the central elements
 * of the mesh or leave them untouched. Refine the remaining elements.
 *
 *  |x|x|x|x|       |x|x|x|x|
 *  |x|x|x|x|  -->  |x| | |x|
 *  |x|x|x|x|       |x| | |x|
 *  |x|x|x|x|       |x|x|x|x|
 *
 */
static int
t8_adapt_menger_quad ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                      [[maybe_unused]] t8_locidx_t which_tree, t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                      [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                      t8_element_t *elements[])
{
  if (scheme->element_get_level (tree_class, elements[0])) {
    const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
    const int level_element = scheme->element_get_level (tree_class, elements[0]);
    const int ancestor_id = scheme->element_get_ancestor_id (tree_class, elements[0], level_element - 1);

    if (0 == level_element % 2) {
      /* ancestor_id == 0 && child_id == 3
       ancestor_id == 1 && child_id == 1
       ancestor_id == 2 && child_id == 2
       ancestor_id == 3 && child_id == 0
       + 2 because all ids are shifted in a multilevel scheme */
      if (ancestor_id + child_id == 5) {
        return 0;
      }
    }
  }
  return 1;
}

/* Refine a tree with triangular elements.
 * At every adaptcall remove the central element with child id 3
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_tri ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                         [[maybe_unused]] t8_locidx_t which_tree, t8_eclass_t tree_class,
                         [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                         [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                         t8_element_t *elements[])
{
  const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (level == 0) {
    return 1;
  }
  if (child_id == 3) {
    return 0;
  }
  return 1;
}

void
t8_multilevel_tutorial (const t8_eclass_t eclass, const int level)
{
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  t8_cmesh_t cmesh = t8_cmesh_new_from_class (eclass, comm);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone_multilevel (), 0, 0, comm);
  for (int i_level = 0; i_level < level; ++i_level) {
    t8_forest_t forest_new;
    t8_forest_init (&forest_new);
    switch (eclass) {
    case T8_ECLASS_TRIANGLE:
      t8_forest_set_adapt (forest_new, forest, t8_adapt_sierpinski_tri, 0);
      break;
    case T8_ECLASS_QUAD:
      t8_forest_set_adapt (forest_new, forest, t8_adapt_menger_quad, 0);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_commit (forest_new);
    forest = forest_new;
  }
  std::string filename = "multilevel_tutorial";
  filename += "_" + std::string (t8_eclass_to_string[eclass]);
  t8_forest_write_vtk (forest, filename.c_str ());
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  sc_options_t *opt;  /* The options we want to parse */
  char usage[BUFSIZ]; /* Usage message */
  char help[BUFSIZ];  /* Help message */
  int helpme;         /* Print help message */
  int parsed;         /* Return value of sc_options_parse */
  int sreturn;        /* Return value of sc functions */
  int mpiret;         /* Return value of MPI functions */
  t8_eclass_t eclass;
  int level;

  /* brief help message */
  snprintf (usage, BUFSIZ,
            "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options. \n",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the multilevel scheme conversion of t8code.\n"
                      "Usage: %s\n",
                      usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* We will use MPI_COMM_WORLD as a communicator. */
  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'e', "eclass", (int *) &eclass, 2, "Tree class used in the cmesh.");
  sc_options_add_int (opt, 'l', "level", &level, 5, "Max refinement level of the resulting forest.");
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed == 0) {
    /*wrong usage */
    t8_global_productionf ("\n\tERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    t8_multilevel_tutorial (eclass, level);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
