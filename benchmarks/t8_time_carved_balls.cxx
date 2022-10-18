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

struct t8_adapt_data
{
  double  midpoint[6][3];
};

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
  const struct t8_adapt_data *adapt_data = 
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);
  
  double  centroid[3];
  double  dist;
  T8_ASSERT (adapt_data != NULL);
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  for (int i = 0; i < 6; i++) {
    dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
    /* remove core of every ball */
    if (dist < 0.4) {
      return -2;
    }
  }
  for (int i = 0; i < 6; i++) {
    dist = t8_vec_dist(adapt_data->midpoint[i], centroid);
    /* refine shell of every ball */
    if (dist < 0.49) {
      return 1;
    }
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
t8_carve_balls (int start_level, int end_level, int output, int coarse)
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_forest_t         forest_partition;
  t8_cmesh_t          cmesh;
  double              adapt_time = 0;
  double              partition_time = 0;
  double              adapt_coarse_time = 0;
  double              partition_coarse_time = 0;
  sc_statinfo_t       times[4];
  char                vtuname[BUFSIZ];
  int                 level;
  int                 procs_sent;

  sc_stats_init (&times[0], "carve");
  sc_stats_init (&times[1], "carve_partition");
  sc_stats_init (&times[2], "coarse");
  sc_stats_init (&times[3], "coarse_partition");

  t8_forest_init (&forest);
  cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 1, 0);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, start_level);
  t8_forest_commit (forest);

  struct t8_adapt_data adapt_data = {{{1.0, 0.5, 0.5},
                                      {0.5, 1.0, 0.5},
                                      {0.5, 0.5, 1.0},
                                      {0.0, 0.5, 0.5},
                                      {0.5, 0.0, 0.5},
                                      {0.5, 0.5, 0.0}}};

  for (level = start_level; level <= end_level; level++) {
    /* Adapt - refine, remove */
    t8_forest_init (&forest_adapt);
    t8_forest_set_profiling (forest_adapt, 1);
    t8_forest_set_user_data (forest_adapt, &adapt_data);
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

    /* Partition the adapted forest */
    t8_forest_init (&forest_partition);
    t8_forest_set_profiling (forest_partition, 1);
    t8_forest_set_partition (forest_partition, forest_adapt, 0);
    t8_forest_commit (forest_partition);
    partition_coarse_time += t8_forest_profile_get_partition_time (forest_partition, &procs_sent);

    forest = forest_partition;
  }


  /* vtu output */
  if (output) {
      snprintf (vtuname, BUFSIZ, "/home/ioannis/VBshare/paraview_export/forest_carved_balls");
      t8_forest_write_vtk (forest, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
  }

  sc_stats_accumulate (&times[0], adapt_time);
  sc_stats_accumulate (&times[1], partition_time);
  sc_stats_accumulate (&times[2], adapt_coarse_time);
  sc_stats_accumulate (&times[3], partition_coarse_time);
  sc_stats_compute (sc_MPI_COMM_WORLD, 4, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 4, times, 1, 1);

  t8_forest_unref (&forest);
}


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 create_forest;
  int                 start_level = 0;
  int                 end_level = 1;
  int                 output = 0;
  int                 coarse = 0;
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
  sc_options_add_int (opt, 'f', "flevel", &end_level, 1,
                      "Final refine level: greater or equal to start level");
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
  else {
    t8_carve_balls (start_level, end_level, output, coarse);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;

}