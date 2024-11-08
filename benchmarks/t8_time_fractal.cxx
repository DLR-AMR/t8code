/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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
t8_adapt_menger_quad (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                      t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements get removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level_element = ts->t8_element_level (elements[0]);
  const int ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element - 1);

  if (0 == level_element % 2) {
    /* ancestor_id == 0 && child_id == 3
       ancestor_id == 1 && child_id == 1
       ancestor_id == 2 && child_id == 2
       ancestor_id == 3 && child_id == 0 */
    if (ancestor_id + child_id == 3) {
      return ret;
    }
  }
  if (level_element < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with triangular elements.
 * At every adaptcall remove the central element with child id 2 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_tri (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements shall be removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2) {
    return ret;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with hexahedral elements.
 * At every second adaptcall (level is even) remove the central elements 
 * of the mesh or leave them untouched. Refine the remaining elements.
 */
static int
t8_adapt_menger_hex (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                     t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements shall be removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level_element = ts->t8_element_level (elements[0]);
  const int ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element - 1);

  if (0 == level_element % 2) {
    if (ancestor_id < 4) {
      if (child_id > 3) {
        if (4 != child_id - ancestor_id) {
          return ret;
        }
      }
      else {
        if (3 == child_id + ancestor_id) {
          return ret;
        }
      }
    }
    else {
      if (child_id > 3) {
        if (11 == child_id + ancestor_id) {
          return ret;
        }
      }
      else {
        if (4 != ancestor_id - child_id) {
          return ret;
        }
      }
    }
  }
  if (level_element < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with tetrahedral elements.
 * At every adaptcall remove the elements with child id 2, 3, 5 and 6 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_tet (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements shall be removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2 || child_id == 3 || child_id == 5 || child_id == 6) {
    return ret;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with prism elements.
 * At every adaptcall remove the elements with child id 2 and 6 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_prism (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements shall be removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 2 || child_id == 6) {
    return ret;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Refine a tree with pyramid elements.
 * At every adaptcall remove the elements with child id 1, 3, 5, 6 and 8 
 * of the mesh or leave it untouched. Refine the remaining elements.
 */
static int
t8_adapt_sierpinski_pyramid (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             t8_locidx_t lelement_id, t8_scheme *ts, const int is_family, const int num_elements,
                             t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int ret = adapt_data[1]; /* -2 if elements shall be removed, 0 else */

  const int child_id = ts->t8_element_child_id (elements[0]);
  const int level = ts->t8_element_level (elements[0]);

  if (child_id == 1 || child_id == 3 || child_id == 5 || child_id == 6 || child_id == 8) {
    return ret;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

/* Coarse every family in the mesh. */
static int
t8_adapt_coarse (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                 t8_scheme *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

/** Constructing different fractals on each element type.
 * \param [in] level_initial Initial level of the uniform forest.
 * \param [in] level_end     Final level of the fractal.
 * \param [in] iterative     1 if fractal is constructed iterative
 *                           0 if recursive
 * \param [in] remove        1 if elements that will not get refined will be 
 *                           removed instead.
 * \param [in] trees         Number of trees the forest will contain.
 * \param [in] eclass        Element type for each tree.
 * \param [in] output        1 if VTU output is needed.
 * \param [in] coarse        Number of times the final fractal shall be coarsen.
 * \param [in] runs          Number of times the hole procedure will be repeated.
 */
static void
t8_construct_fractal (int level_initial, int level_end, const int iterative, const int remove, const int trees,
                      const t8_eclass_t eclass, const int output, const int coarse, const int runs)
{
  T8_ASSERT (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_HEX
             || eclass == T8_ECLASS_TET || eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_PYRAMID);

  /* Quadrilateral and hexahedron elements must have an even initial level 
   * greater 0, such that the refinement pattern can be applied. */
  T8_ASSERT ((eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) && level_initial > 1 && 0 == level_initial % 2
             && 0 == level_end % 2);

  /* Set up userdata to adapt forest. 
   * user_data[0] -> level_end
   * user_data[1] -> {0, -2} -> return for adapt callback */
  int user_data[2] = { level_end, (remove == 1) ? -2 : remove };

  /* Time measurement */
  double time_refine = 0;
  double time_coarse = 0;
  sc_statinfo_t times[3];
  sc_stats_init (&times[0], "refine");
  sc_stats_init (&times[1], "coarse");

  t8_forest_t forest;
  t8_forest_t forest_adapt;
  t8_cmesh_t cmesh;

  for (int run = 0; run < runs; run++) {
    if (output) {
      cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    }
    else {
      cmesh = t8_cmesh_new_bigmesh (eclass, trees, sc_MPI_COMM_WORLD);
    }

    t8_forest_init (&forest);
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default ());
    t8_forest_set_level (forest, level_initial);
    t8_forest_commit (forest);

    T8_ASSERT (level_initial < level_end);
    for (int level = level_initial; level <= level_end; level++) {
      /* Adapt - refine (and remove) */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_user_data (forest_adapt, &user_data);
      switch (eclass) {
      case T8_ECLASS_QUAD:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger_quad, iterative);
        break;
      case T8_ECLASS_TRIANGLE:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tri, iterative);
        break;
      case T8_ECLASS_HEX:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger_hex, iterative);
        break;
      case T8_ECLASS_TET:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet, iterative);
        break;
      case T8_ECLASS_PRISM:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_prism, iterative);
        break;
      default:
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_pyramid, iterative);
        break;
      }
      t8_forest_commit (forest_adapt);
      time_refine += t8_forest_profile_get_adapt_time (forest_adapt);
      forest = forest_adapt;

      /* In case of recursive adapt, only one iteration is required. */
      if (iterative == 0) {
        break;
      }
    }

    for (int level = 0; level < coarse; level++) {
      /* Adapt - coarse */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_adapt (forest_adapt, forest, t8_adapt_coarse, 0);
      t8_forest_commit (forest_adapt);
      time_coarse += t8_forest_profile_get_adapt_time (forest_adapt);
      forest = forest_adapt;
    }

    if (run < runs - 1) {
      t8_forest_unref (&forest);
    }
  }

  /* vtu output */
  if (output) {
    char vtuname[BUFSIZ];
    snprintf (vtuname, BUFSIZ, "forest_fractal_adapt_%s", t8_eclass_to_string[eclass]);
    t8_forest_write_vtk (forest, vtuname);
    t8_debugf ("Output to %s\n", vtuname);
  }

  sc_stats_accumulate (&times[0], time_refine);
  sc_stats_accumulate (&times[1], time_coarse);
  sc_stats_compute (sc_MPI_COMM_WORLD, 2, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 2, times, 1, 1);

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* brief help message */
  int sreturnA = snprintf (usage, BUFSIZ,
                           "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                           "for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* long help message */
  int sreturnB = snprintf (help, BUFSIZ, "This program constructs a fractal mesh.\n\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Parameter for t8_construct_fractal and command line */
  int level_initial = -1;
  int level_end = -1;
  int iterative = -1;
  int trees = 1;
  int output = -1;
  int coarse = -1;
  int remove = -1;
  int runs = -1;
  int eclass_int = -1;
  int helpme;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'u', "uniform_level", &level_initial, -1, "Initial uniform refinement level.");
  sc_options_add_int (opt, 'f', "final_level", &level_end, -1,
                      "Final refine level, greater to initial refinement level.");
  sc_options_add_int (opt, 'i', "iterative", &iterative, 0,
                      "Specify if refining is recursive or iterative.\n"
                      "\t\t\t\t\t1 - refine iterative\n"
                      "\t\t\t\t\t0 - refine recursive (default)");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 4,
                      "Specify the type of elements to use.\n"
                      "\t\t\t\t\t2 - quadrilateral\n"
                      "\t\t\t\t\t3 - triangle\n"
                      "\t\t\t\t\t4 - hexahedron (default)\n"
                      "\t\t\t\t\t5 - tetrahedron\n"
                      "\t\t\t\t\t6 - prism\n"
                      "\t\t\t\t\t7 - pyramid");
  sc_options_add_int (opt, 'd', "delete", &remove, 1,
                      "Specify if elements in fractal should be removed.\n"
                      "\t\t\t\t\t1 - delete elements after refining (default)\n"
                      "\t\t\t\t\t0 - never delete elements");
  sc_options_add_int (opt, 't', "trees", &trees, 512, "Number of trees the forest will contain. The default is 512.");
  sc_options_add_int (opt, 'o', "output", &output, 0,
                      "Specify if mesh should be outputted.\n"
                      "\t\t\t\t     If yes, the forest contains only one tree.\n"
                      "\t\t\t\t\t1 - visual output\n"
                      "\t\t\t\t\t0 - no visual output (default)");
  sc_options_add_int (opt, 'c', "coarse", &coarse, 0, "Number of times to coarse hole mesh.");
  sc_options_add_int (opt, 'r', "runs", &runs, 1,
                      "Number of times the fractal gets constructed. The default is 1.\n"
                      "\t\t\t\t     Note, the runntime summs up.");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if ((eclass_int == 2 || eclass_int == 4)
           && (level_initial < 2 || 0 != level_initial % 2 || 0 != level_end % 2)) {
    t8_global_productionf ("\n\t ERROR: Quadrilateral and hexahedron elements must have an even levels greater 0,\n"
                           "\t\tsuch that the refinement pattern can be applied.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && level_initial >= 0 && level_initial < level_end && (iterative == 0 || iterative == 1)
           && (remove == 0 || remove == 1) && (output == 0 || output == 1) && coarse >= 0 && trees > 0
           && (eclass_int > 1 || eclass_int < 8) && runs > 0) {
    t8_construct_fractal (level_initial, level_end, iterative, remove, trees, (t8_eclass_t) eclass_int, output, coarse,
                          runs);
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
