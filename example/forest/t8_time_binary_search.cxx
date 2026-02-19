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
#include "sc.h"
#include "t8_data/t8_containers.h"
#include "t8_forest/t8_forest_private.h"
#include "t8_schemes/t8_scheme.h"
#include "t8_schemes/t8_scheme.hxx"
#include <t8_vtk/t8_vtk_writer.h>

static void
t8_test_binary_search (t8_eclass_t eclass, int level, int num_elements)
{
  const t8_scheme *scheme = t8_scheme_new_standalone ();
  t8_element_array_t elements;
  t8_element_array_init_size (&elements,scheme, eclass, num_elements);
  //create element array, first element first descendant of root at level, then successor
  t8_element_t *first_elem = t8_element_array_index_int_mutable (&elements, 0);
  scheme->set_to_root (eclass, first_elem);
  scheme->element_get_first_descendant(eclass, first_elem, first_elem, level);
  //check that num_elements <= num_descendants
  for (int i = 1; i < num_elements; i++) {
    scheme->element_construct_successor (eclass, t8_element_array_index_int (&elements, i - 1),
                                         t8_element_array_index_int_mutable (&elements, i));
  }

  double time = -sc_MPI_Wtime ();
  for (int i = 0; i < num_elements; i++) {
    t8_element_t *element = t8_element_array_index_int_mutable (&elements, i);
    t8_linearidx_t linear_id = scheme->element_get_linear_id (eclass, element, level);
    int found = t8_forest_bin_search_lower (&elements, linear_id, level);
    if (found != i) {
      SC_ABORTF ("i: %i, found: %i\n", i, found);
    }
  }
  time += sc_MPI_Wtime ();
  t8_infof ("time of linearid binary search: %f\n", time);
  // for each element, binary search for it in 2 ways and check that the right index was calculated.
  time = -sc_MPI_Wtime ();
  for (int i = 0; i < num_elements; i++) {
    t8_element_t *element = t8_element_array_index_int_mutable (&elements, i);
    t8_debugf("################################ \n");
    t8_debugf("search for element with index %i \n",i);
    int found = t8_forest_bin_search_lower_compare (&elements, element);
    if (found != i) {
      SC_ABORTF ("i: %i, found: %i\n", i, found);
    }
  }
  time += sc_MPI_Wtime ();
  t8_infof ("time of compare  binary search: %f\n", time);
  t8_element_array_reset(&elements);
  delete scheme;
}

int
main (int argc, char **argv)
{
  int mpiret, parsed, eclass_int, helpme;
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
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 2,
                      "The element type.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid\n");
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
    t8_global_productionf ("Testing binary search with %s elements\n",
                           t8_eclass_to_string[eclass_int]);
    int start_level, end_level, num_elements;
    int dim = t8_eclass_to_dimension[eclass_int];
    switch (dim) {
      case 3:
      start_level = 5;
      num_elements = 1<<15;
      end_level = 20;
      if(eclass_int == 7){
        end_level = 19;
      }
        break;
      case 2:
      start_level = 5;
      num_elements = 1<<10;
      // start_level = 7;
      // num_elements = 1<<14;
      end_level = 29;
        break;
      default:
        SC_ABORT("no parameters given for d<=1!\n");
    }
    for(int level = start_level; level < end_level; level++){
      t8_global_productionf("level: %i:\n",level);
      t8_test_binary_search ((t8_eclass_t) eclass_int, level, num_elements);
    }
    for(int level = end_level-1; level >= start_level; level--){
      t8_global_productionf("level: %i:\n",level);
      t8_test_binary_search ((t8_eclass_t) eclass_int, level, num_elements);
    }
  }
  /* Clean-up */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
