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
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_dtri.h>
#include <t8_schemes/t8_default/t8_dtet.h>
#include <t8_schemes/t8_default/t8_dprism.h>
#include <t8_schemes/t8_default/t8_dpyramid.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>

static int
t8_basic_hypercube_refine(t8_forest_t forest, t8_forest_t forest_from,
                          t8_locidx_t which_tree, t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts, int num_elements,
                          t8_element_t * elements[])
{
    int         type, level;
    level = ts->t8_element_level(elements[0]);
    if (level >= *(int *) t8_forest_get_user_data (forest) || level == 0) {
      return 0;
    }
    switch(ts->t8_element_shape(elements[0])){
    case T8_ECLASS_TRIANGLE:
        type = ((t8_dtri_t *)elements[0])->type;
        break;
    case T8_ECLASS_TET:
        type = ((t8_dtet_t *)elements[0])->type;
        break;
    case T8_ECLASS_PRISM:
        type = ((t8_dprism_t *)elements[0])->tri.type;
        break;
    case T8_ECLASS_PYRAMID:
        type = ((t8_dpyramid_t *)elements[0])->type;
        break;
    }
/*TODO: Implement a better refinement argument, one that effects all elements*/
    if(type == 6){
        return 1;
    }
    return 0;
}

static int
t8_basic_hypercube_coarsen(t8_forest_t forest, t8_forest_t forest_from,
                          t8_locidx_t which_tree, t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts, int num_elements,
                          t8_element_t * elements[])
{
    int         type, level;
    level = ts->t8_element_level(elements[0]);
    if(level > 1 && num_elements >1) return -1;
    return 0;
}

static void
t8_basic_hypercube (t8_eclass_t eclass, int set_level, int do_adapt,
                    int create_forest, int do_partition, int do_balance)
{
  t8_forest_t         forest, forest_adapt, forest_coarsen;
  t8_cmesh_t          cmesh, cmesh_partition;
  char                vtuname[BUFSIZ], cmesh_file[BUFSIZ];
  int                 mpirank, mpiret, endlvl = set_level+3;

  t8_global_productionf ("Contructing hypercube mesh with element class %s\n",
                         t8_eclass_to_string[eclass]);
  cmesh =
        t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, do_partition, 0);
        //t8_cmesh_new_from_class(eclass, sc_MPI_COMM_WORLD);

  snprintf (cmesh_file, BUFSIZ, "cmesh_hcube_%s",
            t8_eclass_to_string[eclass]);
  t8_cmesh_save (cmesh, cmesh_file);

  if (do_partition) {
    /* repartition the cmesh to match the desired forest partition */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_partition_uniform (cmesh_partition, set_level,
                                    t8_scheme_new_default_cxx ());
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);
    cmesh = cmesh_partition;
  }


  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  snprintf (vtuname, BUFSIZ, "cmesh_hypercube_%s",
            t8_eclass_to_string[eclass]);
  if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
    t8_debugf ("Output to %s\n", vtuname);
  }
  else {
    t8_debugf ("Error in output\n");
  }
  if (create_forest) {
    t8_forest_init (&forest);
    t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());

    t8_forest_set_level (forest, set_level);

    if (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX
        || eclass == T8_ECLASS_TRIANGLE || eclass == T8_ECLASS_TET
        || eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_PRISM
        || eclass == T8_ECLASS_PYRAMID) {
      t8_forest_commit (forest);
      t8_debugf ("Successfully committed forest.\n");
      snprintf (vtuname, BUFSIZ, "forest_hypercube_%s",
                t8_eclass_to_string[eclass]);
      t8_forest_write_vtk (forest, vtuname);
      t8_debugf ("Output to %s\n", vtuname);
      if(do_adapt){
          t8_forest_init(&forest_adapt);
          t8_forest_set_user_data(forest_adapt, &endlvl);
          t8_forest_set_profiling(forest_adapt, 1);
          t8_forest_set_adapt(forest_adapt, forest, t8_basic_hypercube_refine, 1);
          t8_forest_set_ghost_ext(forest_adapt, 1, T8_GHOST_FACES, 2);
          t8_forest_commit(forest_adapt);
          t8_debugf ("Successfully adapted forest.\n");
          snprintf (vtuname, BUFSIZ, "forest_hypercube_adapt_%s",
                    t8_eclass_to_string[eclass]);
          t8_forest_write_vtk (forest_adapt, vtuname);
          t8_debugf ("Output to %s\n", vtuname);
          /* Ensure that the correct forest is passed to unref later */

          t8_forest_init(&forest_coarsen);
          t8_forest_set_user_data(forest_coarsen, &endlvl);
          t8_forest_set_profiling(forest_coarsen, 1);
          t8_forest_set_adapt(forest_coarsen, forest_adapt, t8_basic_hypercube_coarsen, 0);
          t8_forest_commit(forest_coarsen);
          t8_debugf ("Successfully coarsened forest.\n");
          snprintf (vtuname, BUFSIZ, "forest_hypercube_coarsen_%s",
                    t8_eclass_to_string[eclass]);
          t8_forest_write_vtk (forest_coarsen, vtuname);
          t8_debugf ("Output to %s\n", vtuname);
      }
      if (do_balance) {
        t8_forest_t         forest_balance;

        t8_forest_init (&forest_balance);
        t8_forest_set_balance (forest_balance, forest, 1);
        t8_forest_commit (forest_balance);

        t8_debugf ("Successfully balanced forest.\n");
        snprintf (vtuname, BUFSIZ, "forest_hypercube_balanced_%s",
                  t8_eclass_to_string[eclass]);
        t8_forest_write_vtk (forest_balance, vtuname);
        t8_debugf ("Output to %s\n", vtuname);
        /* Ensure that the correct forest is passed to unref later */
        forest = forest_balance;
      }
    }
    if(do_adapt){
        t8_forest_unref(&forest_coarsen);
    }
    else{
        t8_forest_unref (&forest);
    }
  }
  else {
    t8_cmesh_unref (&cmesh);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level, do_partition, create_forest, do_balance, adapt;
  int                 eclass_int;
  int                 parsed, helpme;
  t8_eclass_t         eclass;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "This program constructs a uniformly refined "
            "cubical mesh.\nThe user can choose the type of mesh elements to "
            "use and the refinement level of the mesh.\n\n%s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_switch (opt, 'p', "partition", &do_partition,
                         "Enable coarse mesh partitioning.");
  sc_options_add_switch (opt, 'b', "balance", &do_balance,
                         "Additionally balance the forest.");
  sc_options_add_switch (opt, 'a', "adapt", &adapt, "Refine and coarsen the mesh "
                                                    "adaptivly.");
  sc_options_add_int (opt, 'e', "elements", &eclass_int, 0,
                      "The type of elements to use.\n"
                      "\t\t0 - vertex\n\t\t1 - line\n\t\t2 - quad\n"
                      "\t\t3 - triangle\n\t\t4 - hexahedron\n"
                      "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && 0 <= eclass_int
           && eclass_int < T8_ECLASS_COUNT) {
    create_forest = 1;
    eclass = (t8_eclass_t) eclass_int;
    t8_basic_hypercube (eclass, level, adapt, create_forest, do_partition,
                        do_balance);
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
