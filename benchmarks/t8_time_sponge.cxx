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
#include <t8_forest.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

static int
t8_adapt_menger (t8_forest_t forest,
                t8_forest_t forest_from,
                t8_locidx_t which_tree,
                t8_locidx_t lelement_id,
                t8_eclass_scheme_c * ts,
                const int is_family,
                const int num_elements, 
                t8_element_t * elements[])
{
  int child_id;
  int ancestor_id;
  int level_element;

  level_element = ts->t8_element_level (elements[0]);

  if (0 == level_element%2) {
    child_id = ts->t8_element_child_id (elements[0]);
    ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element-1);
    if (ancestor_id < 4) {
      if (child_id > 3) {
        if (4 != child_id - ancestor_id) {
          return -2;
        }
      }
      else {
        if (3 == child_id + ancestor_id) {
          return -2;
        }
      }
    }
    else {
      if (child_id > 3) {
        if (11 == child_id + ancestor_id) {
          return -2;
        }
      }
      else {
        if (4 != ancestor_id - child_id) {
          return -2;
        }
      }
    }
  }
  if (level_element < *(int *) t8_forest_get_user_data (forest)) {
    return 1;
  }

  return 0;
}

static int
t8_adapt_sierpinski_tet (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         const int is_family,
                         const int num_elements, 
                         t8_element_t * elements[])
{
  int child_id, level;
  child_id = ts->t8_element_child_id (elements[0]);  
  if (child_id == 2 || child_id == 3 
      || child_id == 5 || child_id == 6) {
    return -2;
  }
  level = ts->t8_element_level (elements[0]);
  if (level < *(int *) t8_forest_get_user_data (forest)) {
    return 1;
  }
  return 0;
}

static int
t8_adapt_sierpinski_tri (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c * ts,
                         const int is_family,
                         const int num_elements, 
                         t8_element_t * elements[])
{
  int child_id, level;
  child_id = ts->t8_element_child_id (elements[0]);  
  if (child_id == 2) {
    return -2;
  }
  level = ts->t8_element_level (elements[0]);
  if (level < *(int *) t8_forest_get_user_data (forest)) {
    return 1;
  }
  return 0;
}

static void
t8_construct_sponge (int end_level, t8_eclass_t eclass, int output)
{
  t8_forest_t         forest, forest_adapt, forest_partition;
  t8_cmesh_t          cmesh;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[1];
  char                vtuname[BUFSIZ];
  int                 start_level;

  if (eclass == T8_ECLASS_TRIANGLE) {
    start_level = 0;
  }
  else if (eclass == T8_ECLASS_HEX) {
    start_level = 2;
  }
  else {
    start_level = 0;
  }

  T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_TET);
  t8_forest_init (&forest);
  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  //cmesh = t8_cmesh_new_bigmesh (eclass, 1, sc_MPI_COMM_WORLD);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, start_level);

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  t8_forest_commit (forest);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");

  /* Adapt */
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &end_level);
  t8_forest_set_profiling (forest_adapt, 1);
  if (eclass == T8_ECLASS_TRIANGLE) {
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tri, 1);
  }
  else if (eclass == T8_ECLASS_HEX) {
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger, 1);
  }
  else {
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet, 1);
  }

  forest_partition = forest_adapt;
  /* partition the adapted forest */
  t8_forest_set_partition (forest_partition, NULL, 0);
  /* enable profiling for the partitioned forest */
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);
  t8_forest_print_profile (forest_partition);
  if (output == 1) {
      snprintf (vtuname, BUFSIZ, "forest_sponge_adapt_%s",
                t8_eclass_to_string[eclass]);
      t8_forest_write_vtk (forest_partition, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
  }
  t8_forest_unref (&forest_partition);
}


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 create_forest;
  int                 end_level = 4;
  int                 output = 0;
  int                 eclass_int;
  int                 parsed;
  int                 helpme;
  int                 sreturnA;
  int                 sreturnB;

  /* brief help message */
  sreturnA = snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                       "for a brief overview of all options.",
                       basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturnB = snprintf (help, BUFSIZ,
                       "This program constructs a sponge mesh."
                       "\nThe user can choose the initial refinement level and the final\n"
                       "refinement level of the mesh. If not set, the final level is 4.\n"
                       "The program has no visual output, if desired.\n\n%s\n",
                       usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf
      ("Warning: Truncated usage string and help message to '%s' and '%s'\n",
       usage, help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'f', "flevel", &end_level, 4,
                      "Final refine level: greater or equal to initial refine level");
  sc_options_add_int (opt, 'o', "output", &output, 0,
                      "output = 1 -> visual output.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 4,
                      "This option specifies"
                      " the type of elements to use.\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n\t\t5 - tetrahedron");                    

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 4 <= end_level
          && (eclass_int == 3 || eclass_int == 4 || eclass_int == 5)) {
    t8_construct_sponge (end_level, (t8_eclass_t) eclass_int, output);
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