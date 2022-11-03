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
  int level_element;
  int level_max;
  int ancestor_id;

  child_id = ts->t8_element_child_id (elements[0]);  
  level_element = ts->t8_element_level (elements[0]);
  level_max = *(int *) t8_forest_get_user_data (forest);
  ancestor_id = ts->t8_element_ancestor_id (elements[0], level_element-1);
  
  if (0 == level_element%2) {
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
  if (level_element < level_max) {
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
  int child_id;
  int level;
  int level_max;

  child_id = ts->t8_element_child_id (elements[0]);  
  level = ts->t8_element_level (elements[0]);
  level_max = *(int *) t8_forest_get_user_data (forest);

  if (child_id == 2 || child_id == 3 
      || child_id == 5 || child_id == 6) {
    return -2;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

static int
t8_adapt_sierpinski_prism (t8_forest_t forest,
                           t8_forest_t forest_from,
                           t8_locidx_t which_tree,
                           t8_locidx_t lelement_id,
                           t8_eclass_scheme_c * ts,
                           const int is_family,
                           const int num_elements, 
                           t8_element_t * elements[])
{
  int child_id;
  int level;
  int level_max;

  child_id = ts->t8_element_child_id (elements[0]);  
  level = ts->t8_element_level (elements[0]);
  level_max = *(int *) t8_forest_get_user_data (forest);

  if (child_id == 2 || child_id == 6) {
    return -2;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}

static int
t8_adapt_sierpinski_pyramid (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t which_tree,
                             t8_locidx_t lelement_id,
                             t8_eclass_scheme_c * ts,
                             const int is_family,
                             const int num_elements, 
                             t8_element_t * elements[])
{
  int child_id;
  int level;
  int level_max;

  child_id = ts->t8_element_child_id (elements[0]);  
  level = ts->t8_element_level (elements[0]);
  level_max = *(int *) t8_forest_get_user_data (forest);
  
  if (child_id == 2 || child_id == 3 
      || child_id == 5 || child_id == 6) {
    return -2;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}


static int
t8_adapt_sierpinski_tet_full (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t which_tree,
                              t8_locidx_t lelement_id,
                              t8_eclass_scheme_c * ts,
                              const int is_family,
                              const int num_elements, 
                              t8_element_t * elements[])
{
  int child_id;
  int level;
  int level_max;

  child_id = ts->t8_element_child_id (elements[0]);  
  level = ts->t8_element_level (elements[0]);
  level_max = *(int *) t8_forest_get_user_data (forest);
  
  if (child_id == 2 || child_id == 3 
      || child_id == 5 || child_id == 6) {
    return 0;
  }
  if (level < level_max) {
    return 1;
  }
  return 0;
}


static int
t8_adapt_coarse (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c * ts,
                 const int is_family,
                 const int num_elements, 
                 t8_element_t * elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

static void
t8_construct_sponge (int end_level_it, int end_level_rec, t8_eclass_t eclass, int remove, int output, int coarse)
{
  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_forest_t         forest_partition;
  t8_cmesh_t          cmesh;
  double              adapt_time = 0;
  double              partition_time = 0;
  double              adapt_rec_time = 0;
  double              adapt_coarse_time = 0;
  double              partition_coarse_time = 0;
  sc_statinfo_t       times[5];
  char                vtuname[BUFSIZ];
  int                 level;
  int                 start_level;
  int                 procs_sent;
  int                 runs = 30;

  sc_stats_init (&times[0], "refine_it");
  sc_stats_init (&times[1], "refine_it_partition");
  sc_stats_init (&times[2], "refine_rec");
  sc_stats_init (&times[3], "coarse");
  sc_stats_init (&times[4], "coarse_partition");

  if (eclass == T8_ECLASS_HEX) {
    start_level = 2;
    if (end_level_it < start_level) {
      end_level_it = 4;
    }
    else if (0 != end_level_it%2) {
      end_level_it--;
    }
  }
  else {
    start_level = 0;
  }

  if (end_level_it > end_level_rec) {
    end_level_rec = end_level_it;
  }

  T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_TET
             || eclass == T8_ECLASS_PRISM || eclass == T8_ECLASS_PYRAMID);

  for (int run = 0; run < runs; run++) {

#if 1
    t8_forest_init (&forest);
    //cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    cmesh = t8_cmesh_new_bigmesh (eclass, 512, sc_MPI_COMM_WORLD),
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
    t8_forest_set_level (forest, end_level_it);
    t8_forest_commit (forest);
#else
    t8_forest_init (&forest);
    //cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    cmesh = t8_cmesh_new_bigmesh (eclass, 128, sc_MPI_COMM_WORLD),
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
    t8_forest_set_level (forest, start_level);
    t8_forest_commit (forest);

    for (level = start_level; level <= end_level_it; level++) {
      /* Adapt - refine, remove */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_user_data (forest_adapt, &end_level_it);
      if (remove == 1) {
        if (eclass == T8_ECLASS_HEX) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger, 0);
        }
        else if (eclass == T8_ECLASS_TET) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet, 0);
        }
        else if (eclass == T8_ECLASS_PRISM) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_prism, 0);
        }
        else {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_pyramid, 0);
        }
      }
      else {
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet_full, 0);
      }

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
#endif

    t8_global_productionf("##############################################################\n");
    t8_global_productionf("######################### END ITERATE ########################\n");
    t8_global_productionf("##############################################################\n");

    if (end_level_it < end_level_rec) {
    /* Adapt - refine, remove */
      t8_forest_init (&forest_adapt);
      t8_forest_set_profiling (forest_adapt, 1);
      t8_forest_set_user_data (forest_adapt, &end_level_rec);
      if (remove == 1) {
        if (eclass == T8_ECLASS_HEX) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_menger, 1);
        }
        else if (eclass == T8_ECLASS_TET) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet, 1);
        }
        else if (eclass == T8_ECLASS_PRISM) {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_prism, 1);
        }
        else {
          t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_pyramid, 1);
        }
      }
      else {
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_sierpinski_tet_full, 1);
      }

      t8_forest_commit (forest_adapt);
      adapt_rec_time += t8_forest_profile_get_adapt_time(forest_adapt);

      forest = forest_adapt;
    }

    t8_global_productionf("##############################################################\n");
    t8_global_productionf("######################## END RECURSIVE #######################\n");
    t8_global_productionf("##############################################################\n");

    if (coarse) {
      for (level = 0; level <= end_level_rec; level++) {
        /* Adapt - coarse */
        t8_forest_init (&forest_adapt);
        t8_forest_set_profiling (forest_adapt, 1);
        t8_forest_set_adapt (forest_adapt, forest, t8_adapt_coarse, 0);
        t8_forest_commit (forest_adapt);
        adapt_coarse_time += t8_forest_profile_get_adapt_time(forest_adapt);

        /* Partition the adapted forest */
#if 0
        t8_forest_init (&forest_partition);
        t8_forest_set_profiling (forest_partition, 1);
        t8_forest_set_partition (forest_partition, forest_adapt, 0);
        t8_forest_commit (forest_partition);
        partition_coarse_time += t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
        
        forest = forest_partition;
#else
        forest = forest_adapt;
#endif
      }
    }

    t8_global_productionf("##############################################################\n");
    t8_global_productionf("########################## END COARSE ########################\n");
    t8_global_productionf("##############################################################\n");

    if (run < runs-1) {
      t8_forest_unref (&forest);
    }

  }
  
  adapt_time = adapt_time * (1/(double) runs);
  partition_time = partition_time * (1/(double) runs);
  adapt_rec_time = adapt_rec_time * (1/(double) runs);

  adapt_coarse_time = adapt_coarse_time * (1/(double) runs);
  partition_coarse_time = partition_coarse_time * (1/(double) runs);

  /* vtu output */
  if (output) {
      snprintf (vtuname, BUFSIZ, "/home/ioannis/VBshare/paraview_export/forest_sponge_adapt_%s",
                t8_eclass_to_string[eclass]);
      t8_forest_write_vtk (forest, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
  }

  sc_stats_accumulate (&times[0], adapt_time);
  sc_stats_accumulate (&times[1], partition_time);
  sc_stats_accumulate (&times[2], adapt_rec_time);
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
  int                 end_level_it = 4;
  int                 end_level_rec = 4;
  int                 output = 0;
  int                 coarse = 0;
  int                 remove = 1;
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
                       "\nThe user can choose the the final refinement level\n"
                       "of the mesh. If not set, the final level is 4.\n"
                       "The program has a visual output, if desired.\n\n%s\n",
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
  sc_options_add_int (opt, 'i', "iteration", &end_level_it, 4,
                      "Final refine level of iteration: greater to 0");
  sc_options_add_int (opt, 'r', "recursion", &end_level_rec, 4,
                      "Final refine level of recursion: greater to 0");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 4,
                      "This option specifies the type of elements to use.\n"
                      "\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n"
                      "\t\t6 - prism\n"
                      "\t\t7 - pyramid");
  sc_options_add_int (opt, 'd', "remove", &remove, 1,
                      "remove = 0 -> do not remove elements\n"
                      "\t\t\t\t\tOnly works for tetrahedrons.");
  sc_options_add_int (opt, 'o', "output", &output, 0,
                      "output = 1 -> visual output.");                    
  sc_options_add_int (opt, 'c', "coarse", &coarse, 0,
                      "coarse = 1 -> coarse all elements f times, if possible"); 

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= end_level_it
          && (eclass_int > 3 || eclass_int < 8)) {
    t8_construct_sponge (end_level_it, end_level_rec,(t8_eclass_t) eclass_int, remove, output, coarse);
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