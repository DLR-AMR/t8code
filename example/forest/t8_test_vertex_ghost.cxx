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

#include <t8.h>
#include <t8_eclass.h>
#include <sc_options.h>
#include <sc_mpi.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_vtk/t8_vtk_writer.h>

/* Build a forest on a hypercube mesh
 * and refine the first tree of a process once.
 * Create ghost layer and print it.
 * partition the forest, create ghost layer and print it. */




 static void
t8_test_ghost_hypercube (t8_eclass_t eclass, int level, sc_MPI_Comm comm, int no_vtk)
{
  t8_cmesh_t cmesh;

  if (eclass < T8_ECLASS_COUNT) {
    // Build a hypercube out of the given eclass
    cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
  }
  else if (eclass == T8_ECLASS_COUNT) {
    // Build a 3D hybrid hypercube with tets, hexes and prisms
    cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  }
  t8_forest_t forest;

  /* Initialize the forest */
  t8_forest_init (&forest);

  /* Set the cmesh, scheme and level */
  t8_forest_set_cmesh (forest, cmesh, comm);
  t8_forest_set_scheme (forest, t8_scheme_new_standalone());
  t8_forest_set_level (forest, level);
  t8_forest_set_ghost (forest, 1, T8_GHOST_VERTICES);
  /* commit the forest */
  t8_forest_commit (forest);

  char fileprefix[BUFSIZ] ="ghost";

  t8_forest_vtk_write_file(forest,fileprefix,1,1,1,1,1,0,nullptr);
  t8_forest_unref(&forest);

}


int
main (int argc, char **argv)
{
  int mpiret, parsed, eclass_int, level, helpme;
  int dim, no_vtk;
  sc_options_t *opt;
  const char *prefix;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int sreturnA, sreturnB;

  sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>", basename (argv[0]));
  sreturnB = snprintf (help, BUFSIZ, "help string\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);
  /*
   * COMMAND LINE OPTION SETUP
   */
  opt = sc_options_new (argv[0]);
  /* Level -l */
  sc_options_add_int (opt, 'l', "level", &level, 0, "The refinement level of the mesh.");
  /* enable/disable vtk -o */
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "disable vtk output");
  /* Use a hypercube mesh -e */
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 3,
                      "The element type.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid\n"
                      "\t\t8 - hex/tet/prism hybrid");
  /* Print help -h */
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  /*
   * END OF COMMAND LINE OPTION SETUP
   */
  /* parse command line options */
  parsed = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  /* check for wrong usage of arguments */
  if (parsed < 0 || parsed != argc || eclass_int < T8_ECLASS_VERTEX || eclass_int > T8_ECLASS_COUNT) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (helpme) {
    /* Print help string and then exit */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    t8_global_productionf ("Testing ghost on a hypercube cmesh with %s elements\n",
                            eclass_int < T8_ECLASS_COUNT ? t8_eclass_to_string[eclass_int] : "hybrid");
    t8_test_ghost_hypercube ((t8_eclass_t) eclass_int, level, sc_MPI_COMM_WORLD, no_vtk);
  }

  /* Clean-up */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
