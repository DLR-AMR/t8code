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

#include <sc_refcount.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <t8_schemes/t8_scheme.hxx>

static const t8_scheme *
t8_time_elements_scheme (int scheme_option)
{
  const t8_scheme *scheme = NULL;
  switch (scheme_option) {
  case 1:
    scheme = t8_scheme_new_standalone ();
    break;
  case 2:
    scheme = t8_scheme_new_default ();
    break;
  default:
    t8_global_errorf (" [search] Unknown scheme option %d.\n", scheme_option);
    SC_ABORT ("Not implemented.");
    break;
  }
  return scheme;
}

static int
adapt_all ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
           [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
           [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme_c *scheme,
           [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t *elements[])
{
  return 1;
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  t8_cmesh_t cmesh;
  t8_forest_t forest;
  sc_options_t *opt;
  int particle_option;
  int scheme_option;
  int eclass_option;
  int help = 0;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);
  comm = sc_MPI_COMM_WORLD;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is the search benchmark of t8code.\n");
  t8_global_productionf (
    " [search] We will search for all elements in a forest that contain randomly created particles.\n");
  t8_global_productionf (" [search] \n");

  opt = sc_options_new (argv[1]);
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_int (opt, 'p', "particle_fill", &particle_option, 2,
                      "Option to fill the particles, 1: random, 2: one per leaf, 3: only one particle");
  sc_options_add_int (opt, 's', "scheme", &scheme_option, 2,
                      "Option to choose the scheme, 1: standalone scheme, 2: default scheme");
  sc_options_add_int (opt, 'e', "elements", &eclass_option, 4,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron (default)\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");

  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_from_class ((t8_eclass) eclass_option, comm);
  /* Build a uniform forest on it. */
  const t8_scheme *scheme = t8_time_elements_scheme (scheme_option);
  forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, comm);
  int num_leaves, level;
  while (1) {
    num_leaves = t8_forest_get_tree_num_leaf_elements (forest, 0);
    forest = t8_forest_new_adapt (forest, adapt_all, 0, 0, 0);
    t8_element_array_t *element_array = t8_element_array_new_count (scheme, (t8_eclass) eclass_option, num_leaves);
    element_array = t8_forest_tree_get_leaf_elements (forest, 0);
    const t8_element_t *element = t8_element_array_index_int (element_array, 0);
    level = scheme->element_get_level ((t8_eclass) eclass_option, element);
    t8_productionf ("level: %i \n", level);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
