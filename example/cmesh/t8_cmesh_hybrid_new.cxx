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

#include <sc_options.h>
#include <sc_flops.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>
#include <sc_statistics.h>

/**
 * Adaptation criterion for a forest containing hexahedra, tetrahedra, prisms and/or pryramids. 
 * Refines every second element, except for pyramids. 
 * Refines every pyramid
 * 
 * \param [in] forest the forest
 * \param [in] forest_from 
 * \param [in] which_tree The local id of the tree
 * \param [in] lelement_id the local id of the element
 * \param [in] ts the scheme to use
 * \param [in] is_family flag, if the \a elements form a family
 * \param [in] num_elements number of elements
 * \param [in] elements A single element or a collection of \a num_elements
 * \return int 
 */
static int
t8_basic_hybrid_refine (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                        [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] t8_eclass_t tree_class,
                        [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                        [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                        t8_element_t *elements[])
{
  /*If the level is equal or higher than given by the user do not refine*/
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  else {
    const int id = scheme->element_get_child_id (tree_class, elements[0]);
    /*Refine every second element */
    switch (scheme->element_get_shape (tree_class, elements[0])) {
    case T8_ECLASS_HEX:
      return id % 2 == 0 ? 1 : 0;
    case T8_ECLASS_TET:
      return id % 2 == 0 ? 1 : 0;
    case T8_ECLASS_PRISM:
      return id % 2 == 0 ? 1 : 0;
    case T8_ECLASS_PYRAMID:
      /* Refine every element. */
      return 1;
    default:
      return 1;
    }
  }
}

/**
 * Creates a cmesh and forest using this cmesh, executes all high-level algorithms (new, adapt, ghost, partition, ...) and measures their runtime
 * 
 * \param[in] level The initial level of forest
 * \param[in] endlvl The maximal level of the forest
 * \param[in] do_vtk Flag, if vtk-output should be produced
 */
static void
t8_basic_hybrid (const int level, int endlvl, const int do_vtk)
{
  /* Create and initialize timing stats */
  sc_statinfo_t times[6];
  sc_stats_init (&times[0], "new");
  sc_stats_init (&times[1], "adapt");
  sc_stats_init (&times[2], "ghost");
  sc_stats_init (&times[3], "partition");
  sc_stats_init (&times[4], "total");
  double total_time = -sc_MPI_Wtime ();

  t8_cmesh_t cmesh;
  /* Create the cmesh */
  t8_global_productionf ("Constructing full hybrid mesh.\n");
  cmesh = t8_cmesh_new_pyramid_deformed (sc_MPI_COMM_WORLD);
  //cmesh = t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);

  /* Partition the cmesh */
  t8_cmesh_t cmesh_partition;
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_uniform (cmesh_partition, level, t8_scheme_new_default ());
  t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);
  cmesh = cmesh_partition;

  /* Save the cmesh and create vtk output if the flag is set. */
  char vtuname[BUFSIZ];
  char cmesh_file[BUFSIZ];
  snprintf (cmesh_file, BUFSIZ, "cmesh_hybrid");
  snprintf (vtuname, BUFSIZ, "cmesh_hybrid");
  t8_cmesh_save (cmesh, cmesh_file);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname) == 0) {
    t8_debugf ("Output to %s\n", vtuname);
  }
  else {
    t8_debugf ("Error in output\n");
  }

  /* Create the forest and measure the time needed for the initial refinement*/
  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_profiling (forest, 1);
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, level);
  double new_time = -sc_MPI_Wtime ();
  t8_forest_commit (forest);
  new_time += sc_MPI_Wtime ();

  if (do_vtk) {
    /* write the initially refined forest into a vtu-file*/
    snprintf (vtuname, BUFSIZ, "forest_hybrid");
    t8_forest_write_vtk (forest, vtuname);
  }

  /* Initialize the adapted forest*/
  t8_forest_t forest_adapt;
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &endlvl);
  t8_forest_set_profiling (forest_adapt, 1);

  /* Set the adaptation criterion to use. */
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_hybrid_refine, 1);

  /* Measure the time needed to adapt */
  double adapt_time = -sc_MPI_Wtime ();
  t8_forest_commit (forest_adapt);
  adapt_time += sc_MPI_Wtime ();

  /* Initialize the partitioned forest */
  t8_forest_t forest_partition;
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);

  /* We want to compute the ghosts, too */
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  /* Use the profiling algorithm to measure the time needed and also how many ghost elements are sent, 
   *  how many procs we sent to and how many rounds we need to balance the forest. */
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);

  /* Write vtk-output*/
  if (do_vtk) {
    snprintf (vtuname, BUFSIZ, "forest_hybrid_partition");
    t8_forest_write_vtk (forest_partition, vtuname);
  }
  total_time += sc_MPI_Wtime ();
  int procs_sent;         /* Procs we sent to during forest_partition */
  t8_locidx_t ghost_sent; /* Number of ghost-elements we sent */
  t8_forest_print_profile (forest_partition);
  const double partition_time = t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
  const double ghost_time = t8_forest_profile_get_ghost_time (forest_partition, &ghost_sent);

  /* Accumulate the timings over all procs and print the stats */
  sc_stats_accumulate (&times[0], new_time);
  sc_stats_accumulate (&times[1], adapt_time);
  sc_stats_accumulate (&times[2], ghost_time);
  sc_stats_accumulate (&times[3], partition_time);
  sc_stats_accumulate (&times[4], total_time);
  sc_stats_compute (sc_MPI_COMM_WORLD, 5, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 5, times, 1, 1);

  t8_forest_unref (&forest_partition);
}

int
main (int argc, char **argv)
{
  int level;
  int endlvl;
  int helpme;
  int do_vtk;
  sc_options_t *opt;
  char help[BUFSIZ];

  /* long help message */
  snprintf (help, BUFSIZ,
            "ADD EXAMPLE DESCRIPTION.\n\n"
            "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.\n",
            basename (argv[0]), basename (argv[0]));
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The refinement level of the mesh.");
  sc_options_add_int (opt, 'f', "final-level", &endlvl, 1, "The final refinement level of the mesh.");
  sc_options_add_switch (opt, 'v', "vtk", &do_vtk, "Enable vtk-output.");

  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level) {
    t8_basic_hybrid (level, endlvl, do_vtk);
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
