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
#include <t8_vec.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

static int
t8_adapt_carve_and_refine (t8_forest_t forest,
                           t8_forest_t forest_from,
                           t8_locidx_t which_tree,
                           t8_locidx_t lelement_id,
                           t8_eclass_scheme_c * ts,
                           const int is_family,
                           const int num_elements,  
                           t8_element_t * elements[])
{
 
  double  centroid[3];
  double  midpoint[3] = {0.5, 0.5, 0.5};
  double  dist;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist(midpoint, centroid);
  if (dist < 0.39) {
    return -2;
  }
  else if (dist < 0.425) {
    return 1;
  }
  return 0;
}

static int
t8_adapt_coarse_all (t8_forest_t forest,
                     t8_forest_t forest_from,
                     t8_locidx_t which_tree,
                     t8_locidx_t lelement_id,
                     t8_eclass_scheme_c * ts,
                     const int is_family,
                     const int num_elements, 
                     t8_element_t * elements[])
{
  if(is_family) {
    return -1;
  }
  return 0;
}

static void
t8_carve_ball (int start_level, int end_level, int eclass_int, int output, int coarse)
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_forest_t         forest_partition;
  t8_cmesh_t          cmesh;
  double              init_time = 0;
  double              init_time_s;
  double              adapt_time = 0;
  double              partition_time = 0;
  double              adapt_coarse_time = 0;
  double              partition_coarse_time = 0;
  sc_statinfo_t       times[5];
  char                vtuname[BUFSIZ];
  int                 level;
  int                 runs = 5;
  int                 procs_sent;
  int64_t             num_elements = INT64_MAX;

  sc_stats_init (&times[0], "init");
  sc_stats_init (&times[1], "carve");
  sc_stats_init (&times[2], "carve_partition");
  sc_stats_init (&times[3], "coarse");
  sc_stats_init (&times[4], "coarse_partition");

  if (start_level >= end_level) {
    end_level = start_level + 1;
  }
  if (coarse == -1) {
    coarse = INT_MAX;
  }

  for (int run = 0; run < runs; run++) {

    t8_forest_init (&forest);
    if (eclass_int == 0){
      cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
    }
    else {
      cmesh = t8_cmesh_new_hypercube ((t8_eclass_t) eclass_int, sc_MPI_COMM_WORLD, 0, 0, 0);
    }

    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
    t8_forest_set_level (forest, start_level);
    init_time_s = sc_MPI_Wtime ();
    t8_forest_commit (forest);
    init_time += sc_MPI_Wtime () - init_time_s;

    for (level = start_level; level < end_level; level++) {
      /* Adapt - refine, remove */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_adapt (forest_adapt, forest, t8_adapt_carve_and_refine, 0);
      t8_forest_commit (forest_adapt);
      adapt_time += t8_forest_profile_get_adapt_time(forest_adapt);

      /* Partition the adapted forest */
      t8_forest_init (&forest_partition);
      t8_forest_set_profiling (forest_partition, 1);
      t8_forest_set_partition (forest_partition, forest_adapt, 0);
      t8_forest_commit (forest_partition);
      partition_time += t8_forest_profile_get_partition_time (forest_partition, &procs_sent);

      forest = forest_partition;
    } 
      
    for (level = 0; level < coarse; level++) {
      /* Adapt - coarse */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_adapt (forest_adapt, forest, t8_adapt_coarse_all, 0);
      t8_forest_commit (forest_adapt);
      adapt_coarse_time += t8_forest_profile_get_adapt_time(forest_adapt);

      if (num_elements > t8_forest_get_global_num_elements(forest_adapt)) {
        num_elements = t8_forest_get_global_num_elements(forest_adapt);
      }
      else {
        coarse = 0;
      }

      /* Partition the adapted forest */
      t8_forest_init (&forest_partition);
      t8_forest_set_profiling (forest_partition, 1);
      t8_forest_set_partition (forest_partition, forest_adapt, 0);
      t8_forest_commit (forest_partition);
      partition_coarse_time += t8_forest_profile_get_partition_time (forest_partition, &procs_sent);

      forest = forest_partition;
    }

    if (run < runs-1) {
      t8_forest_unref (&forest);
    }

  }
  
  init_time = init_time * (1/(double) runs);
  adapt_time = adapt_time * (1/(double) runs);
  partition_time = partition_time * (1/(double) runs);
  adapt_coarse_time = adapt_coarse_time * (1/(double) runs);
  partition_coarse_time = partition_coarse_time * (1/(double) runs);

  /* vtu output */
  if (output) {
      snprintf (vtuname, BUFSIZ, "/home/ioannis/VBshare/paraview_export/forest_carved_balls");
      t8_forest_write_vtk (forest, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
  }

  sc_stats_accumulate (&times[0], init_time);
  sc_stats_accumulate (&times[1], adapt_time);
  sc_stats_accumulate (&times[2], partition_time);
  sc_stats_accumulate (&times[3], adapt_coarse_time);
  sc_stats_accumulate (&times[4], partition_coarse_time);
  sc_stats_compute (sc_MPI_COMM_WORLD, 5, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 5, times, 1, 1);

  t8_forest_unref (&forest);
}


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 start_level = 0;
  int                 end_level = 0;
  int                 output = 0;
  int                 coarse = 0;
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
                       "This program constructs a hybrid hypercube \n\n%s\n",
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
  sc_options_add_int (opt, 's', "slevel", &start_level, 0,
                      "Start refine level: greater to 0");
  sc_options_add_int (opt, 'f', "flevel", &end_level, 0,
                      "Final refine level: greater to start level");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 0,
                      "This option specifies the type of elements to use.\n"
                      "\t\t0 - hybrid\n"
                      "\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n"
                      "\t\t6 - prism\n"
                      "\t\t7 - pyramid");
  sc_options_add_int (opt, 'o', "output", &output, 0,
                      "output = 1 -> visual output.");                    
  sc_options_add_int (opt, 'c', "coarse", &coarse, 0,
                      "number of times to coarse all elements, if possible"); 
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {
    t8_carve_ball (start_level, end_level, eclass_int, output, coarse);
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