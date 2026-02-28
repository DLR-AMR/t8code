/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/* This application allows to measure the forest creation, adaption,
 * partitioning, and ghost layer creation runtimes for different element
 * classes. As a stress test the cubical domain is refined into a fractal
 * pattern.
 *
 * NOTE:
 *  This program is in many regards a duplication of `t8_time_fractal.cxx`. However,
 *  it is tailored towards profiling the ghost layer creation. `t8_time_fractal.cxx` is
 *  tailored towards profiling the element deletion feature of t8code. So called 'incomplete'
 *  trees, however, do not support the creation of ghost layers yet. When this feature is
 *  available, this file may be merged with `t8_time_fractal.cxx`.
 *
 * TODO:
 *  - Allow for multiple trees.
 *  - Allow periodic boundaries.
 *  - Allow for arbitrary cmeshes (nice-to-have).
 */

#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

/* Refine a tree with quadrilateral elements into a Menger carpet.
 * At every second adaptcall (level is even) remove the central elements 
 * of the mesh or leave them untouched. Refine the remaining elements.
 */
static int
t8_adapt_menger_quad (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                      t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level_element = ts->t8_element_level (elements[0]);
  const int ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element - 1);

  if (0 == level_element % 2) {
    /* ancestor_id == 0 && child_id == 3
       ancestor_id == 1 && child_id == 1
       ancestor_id == 2 && child_id == 2
       ancestor_id == 3 && child_id == 0 */
    if (ancestor_id + child_id == 3) {
      return 0;
    }
  }
  if (level_element < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with triangular elements into a Sierpinski carpet.
 * At every adaptcall remove the central element with child id 2 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_tri (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2) {
    return 0;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with hexahedral elements into a Menger sponge.
 * At every second adaptcall (level is even) remove the central elements 
 * of the mesh or leave them untouched. Refine the remaining elements.
 */
static int
t8_adapt_menger_hex (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                     t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level_element = ts->t8_element_level (elements[0]);
  const int ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element - 1);

  if (0 == level_element % 2) {
    if (ancestor_id < 4) {
      if (child_id > 3) {
        if (4 != child_id - ancestor_id) {
          return 0;
        }
      }
      else {
        if (3 == child_id + ancestor_id) {
          return 0;
        }
      }
    }
    else {
      if (child_id > 3) {
        if (11 == child_id + ancestor_id) {
          return 0;
        }
      }
      else {
        if (4 != ancestor_id - child_id) {
          return 0;
        }
      }
    }
  }
  if (level_element < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with tetrahedral elements into a Sierpinski carpet.
 * At every adaptcall remove the elements with child id 2, 3, 5 and 6 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_tet (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2 || child_id == 3 || child_id == 5 || child_id == 6) {
    return 0;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with prism elements into a Sierpinski sponge.
 * At every adaptcall remove the elements with child id 2 and 6 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_prism (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_eclass_scheme_c *ts, const int is_family, const int num_elements,
                           t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2 || child_id == 6) {
    return 0;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with pyramid elements into a Sierpinski sponge.
 * At every adaptcall remove the elements with child id 1, 3, 5, 6 and 8 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_pyramid (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                             const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 1 || child_id == 3 || child_id == 5 || child_id == 6 || child_id == 8) {
    return 0;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/** Constructing different fractals on each element type.
 * \param [in] level_initial Initial level of the uniform forest.
 * \param [in] level_end     Final level of the fractal.
 * \param [in] eclass        Element type for each tree.
 * \param [in] do_balance    Balance the forest after adaptation.
 * \param [in] do_output     1 if VTU output is needed.
 */
static void
t8_construct_fractal (int level_initial, int level_end, const t8_eclass_t eclass, const int do_balance,
                      const int do_output)
{
  // Set up userdata to adapt forest.
  int user_data[1] = { level_end };

  t8_forest_t forest;
  t8_forest_t forest_adapt;
  t8_cmesh_t cmesh;

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);

  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, level_initial);
  t8_forest_commit (forest);

  t8_forest_init (&forest_adapt);
  t8_forest_set_profiling (forest_adapt, 1);
  t8_forest_set_user_data (forest_adapt, &user_data);
  if (do_balance) {
    t8_forest_set_balance (forest_adapt, NULL, 1);
  }
  t8_forest_set_partition (forest_adapt, NULL, 0);
  t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);

  switch (eclass) {
  case T8_ECLASS_QUAD:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger_quad, 1);
    break;
  case T8_ECLASS_TRIANGLE:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tri, 1);
    break;
  case T8_ECLASS_HEX:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger_hex, 1);
    break;
  case T8_ECLASS_TET:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet, 1);
    break;
  case T8_ECLASS_PRISM:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_prism, 1);
    break;
  case T8_ECLASS_PYRAMID:
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_pyramid, 1);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  t8_forest_commit (forest_adapt);

  t8_forest_print_profile (forest_adapt);

  if (do_output) {
    char vtuname[BUFSIZ];
    snprintf (vtuname, BUFSIZ, "forest_fractal_adapt_%s", t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest_adapt, vtuname);
    t8_debugf ("Output to %s\n", vtuname);
  }

  t8_forest_unref (&forest_adapt);
}

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  // Brief help message.
  int sreturnA = snprintf (usage, BUFSIZ,
                           "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                           "for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  // Long help message.
  int sreturnB = snprintf (help, BUFSIZ, "This program constructs a fractal mesh.\n\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    // The usage string or help message was truncated
    // Note: gcc >= 7.1 prints a warning if we
    // do not check the return value of snprintf.
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  // Parameter for t8_construct_fractal and command line.
  int level_initial = -1;
  int level_end = -1;
  int do_output = -1;
  int do_balance = -1;
  int eclass_int = -1;
  int helpme;

  // initialize command line argument parser.
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'u', "uniform_level", &level_initial, 0, "Initial uniform refinement level.");
  sc_options_add_int (opt, 'f', "final_level", &level_end, 3,
                      "Final refine level, greater to initial refinement level.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 4,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron (default)\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");
  sc_options_add_int (opt, 'b', "balance", &do_balance, 0, "Balance the forest after adaption.\n");
  sc_options_add_int (opt, 'o', "output", &do_output, 0,
                      "Specify if mesh should be outputted.\n"
                      "\t\t\t\t     If yes, the forest contains only one tree.\n"
                      "\t\t\t\t\t1 - visual output\n"
                      "\t\t\t\t\t0 - no visual output (default)");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  T8_ASSERT (level_initial >= 0);
  T8_ASSERT (level_initial <= level_end);
  T8_ASSERT (eclass_int >= T8_ECLASS_QUAD);
  T8_ASSERT (eclass_int < T8_ECLASS_COUNT);

  if (helpme) {
    // Display help message and usage.
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {
    t8_construct_fractal (level_initial, level_end, (t8_eclass_t) eclass_int, do_balance, do_output);
  }
  else {
    // Wrong usage.
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
