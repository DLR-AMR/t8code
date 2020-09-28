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
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_dtet.h>
#include <t8_schemes/t8_default/t8_dprism.h>
#include <t8_schemes/t8_default/t8_dpyramid.h>
#include <t8_forest.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include <sc_statistics.h>


static int
t8_basic_hybrid_refine(t8_forest_t forest, t8_forest_t forest_from,
                          t8_locidx_t which_tree, t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts, int num_elements,
                          t8_element_t * elements[])
{
    int         level, id;
    level = ts->t8_element_level(elements[0]);
    if (level >= *(int *) t8_forest_get_user_data (forest)) {
      return 0;
    }
    else{
        switch (ts->t8_element_shape(elements[0])) {
        case T8_ECLASS_HEX:
            id=ts->t8_element_child_id(elements[0]);
            return id%2==0?1:0;
        case T8_ECLASS_TET:
            id=ts->t8_element_child_id(elements[0]);
            return id%2==0?1:0;
        case T8_ECLASS_PRISM:
            id=ts->t8_element_child_id(elements[0]);
            return id%2==0?1:0;
        case T8_ECLASS_PYRAMID:
            return 1;
        default:
            return 1;
        }
    }
}

static int
t8_basic_cake_refine(t8_forest_t forest, t8_forest_t forest_from,
                          t8_locidx_t which_tree, t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts, int num_elements,
                          t8_element_t * elements[])
{
    int         level, type;
    level = ts->t8_element_level(elements[0]);
    if (level >= *(int *) t8_forest_get_user_data (forest)) {
      return 0;
    }
    else{
        int32_t h = T8_DPYRAMID_LEN(level);
        switch (ts->t8_element_shape(elements[0])) {
        case T8_ECLASS_TET:
            type = ((t8_dtet_t *)elements[0])->type;
            if(type == 0 || type == 2 || type == 4){
                return 1;
            }
            else{
                return 0;
            }
        case T8_ECLASS_PYRAMID:

            if(!(((t8_dpyramid_t *) elements[0])->x & h)){
                type = ((t8_dpyramid_t *) elements[0])->type;
                return type == 6 ? 1:0;
            }
            else{
                return 0;
            }
        default:
            return 1;
        }
    }
}

static void
t8_basic_hybrid(int level, int endlvl, int do_vtk, t8_eclass_t eclass,
                int num_elements, int mesh, const char* prefix)
{
    t8_forest_t forest, forest_adapt, forest_partition;
    t8_cmesh_t  cmesh, cmesh_partition;
    char        vtuname[BUFSIZ], cmesh_file[BUFSIZ];
    int         mpirank, mpiret;
    double      new_time = 0, adapt_time = 0, ghost_time = 0, partition_time = 0, total_time = 0;
    sc_statinfo_t   times[5];
    int         procs_sent;
    t8_locidx_t ghost_sent;
    sc_stats_init(&times[0], "new");
    sc_stats_init(&times[1], "adapt");
    sc_stats_init(&times[2], "ghost");
    sc_stats_init(&times[3], "partition");
    sc_stats_init(&times[4], "total");
    total_time -= sc_MPI_Wtime();
    switch (mesh) {
    case 0:
        t8_global_productionf ("Contructing cake mesh with %i pyramids.\n",
                               num_elements);
        cmesh = t8_cmesh_new_pyramid_cake(sc_MPI_COMM_WORLD, num_elements);
        break;
    case 1:
        t8_global_productionf ("Contructing full hybrid mesh.\n");
        cmesh = t8_cmesh_new_full_hybrid(sc_MPI_COMM_WORLD);
        break;
    case 2:
        t8_global_productionf ("Contructing long brick out of %i pyramids.\n",
                               num_elements);
        cmesh = t8_cmesh_new_long_brick_pyramid(sc_MPI_COMM_WORLD, num_elements);
        break;
    case 3:
        t8_global_productionf ("Contructing mesh from: %s.\n", prefix);
        cmesh = t8_cmesh_from_msh_file((char *) prefix, 1, sc_MPI_COMM_WORLD, 3, 0);
        break;
    default:
        t8_global_productionf ("Contructing a single %s.\n", t8_eclass_to_string[eclass]);
        cmesh = t8_cmesh_new_from_class(eclass, sc_MPI_COMM_WORLD);
        break;
    }

    snprintf(cmesh_file, BUFSIZ,"cmesh_hybrid");
    t8_cmesh_save(cmesh, cmesh_file);
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);

    snprintf (vtuname, BUFSIZ, "cmesh_hybrid");
    if (t8_cmesh_vtk_write_file (cmesh, vtuname, 1.0) == 0) {
      t8_debugf ("Output to %s\n", vtuname);
    }
    else {
      t8_debugf ("Error in output\n");
    }

    t8_cmesh_init(&cmesh_partition);
    t8_cmesh_set_derive(cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform(cmesh_partition, level, t8_scheme_new_default_cxx());
    t8_cmesh_commit(cmesh_partition, sc_MPI_COMM_WORLD);
    cmesh = cmesh_partition;

    t8_forest_init(&forest);
    t8_forest_set_profiling(forest, 1);
    t8_forest_set_cmesh(forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme(forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(forest, level);
    new_time -= sc_MPI_Wtime();
    t8_forest_commit(forest);
    new_time += sc_MPI_Wtime();
    if(do_vtk){
        snprintf (vtuname, BUFSIZ, "forest_hybrid");
        t8_forest_write_vtk (forest, vtuname);
        t8_debugf("[D] output to %s\n", vtuname);
    }
    t8_forest_init(&forest_adapt);
    t8_forest_set_user_data(forest_adapt, &endlvl);
    t8_forest_set_profiling(forest_adapt, 1);
    if(eclass == T8_ECLASS_PYRAMID){
        t8_debugf("Use cake-adapt\n");
         t8_forest_set_adapt(forest_adapt, forest, t8_basic_cake_refine, 1);
    }else{
        t8_forest_set_adapt(forest_adapt, forest, t8_basic_hybrid_refine, 1);
    }
    //t8_forest_set_ghost_ext(forest_adapt, 1, T8_GHOST_FACES, 2);
    adapt_time -= sc_MPI_Wtime();
    t8_forest_commit(forest_adapt);
    adapt_time += sc_MPI_Wtime();
    //t8_debugf ("Successfully adapted forest.\n");
    //snprintf (vtuname, BUFSIZ, "forest_hybrid_refine");
    //t8_forest_write_vtk (forest_adapt, vtuname);
    //t8_debugf ("Output to %s\n", vtuname);
    //t8_forest_unref(&forest_adapt);
    /* Ensure that the correct forest is passed to unref later */

    t8_forest_init(&forest_partition);
    t8_forest_set_partition(forest_partition, forest_adapt, 0);
    t8_forest_set_ghost_ext(forest_partition, 1, T8_GHOST_FACES, 2);
    t8_forest_set_profiling(forest_partition, 1);
    t8_forest_commit(forest_partition);
    if(do_vtk){
        snprintf (vtuname, BUFSIZ, "forest_hybrid_partition");
        t8_forest_write_vtk (forest_partition, vtuname);
        t8_debugf ("Output to %s\n", vtuname);
    }
    t8_forest_print_profile(forest_partition);
    partition_time += t8_forest_profile_get_partition_time(forest_partition, &procs_sent);
    ghost_time += t8_forest_profile_get_ghost_time(forest_partition, &ghost_sent);
    total_time += sc_MPI_Wtime();

    sc_stats_accumulate(&times[0], new_time);
    sc_stats_accumulate(&times[1], adapt_time);
    sc_stats_accumulate(&times[2], ghost_time);
    sc_stats_accumulate(&times[3], partition_time);
    sc_stats_accumulate(&times[4], total_time);
    sc_stats_compute(sc_MPI_COMM_WORLD, 5, times);
    sc_stats_print(t8_get_package_id(), SC_LP_ESSENTIAL, 5, times, 1,1);

    t8_forest_unref(&forest_partition);
}

int
main (int argc, char **argv)
{
  int                 mpiret, parsed;
  int                 level, endlvl, helpme, do_vtk, eclass_int, mesh, elements;
  t8_eclass_t         eclass;
  sc_options_t        *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  const char          *file;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "ADD EXAMPLE DESCRIPTION.\n\n%s\n", usage);
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  opt = sc_options_new(argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_int (opt, 'f', "final-level", &endlvl, 1,
                      "The final refinement level of the mesh.");
  sc_options_add_switch(opt, 'v', "vtk", &do_vtk, "Enable vtk-output.");
  sc_options_add_int(opt, 'e', "element", &eclass_int, 4, "Given an element-class, the programm will "
                    " construct a single element of this class. Is ignored, if the option cake is chosen."
                  "The type of elements to use.\n"
                    "\t\t4 - hexahedron\n"
                  "\t\t5 - tetrahedron\n\t\t6 - prism\n\t\t7 - pyramid");
  sc_options_add_int(opt, 'm', "mesh", &mesh, 4,"A mesh to choose from."
                     " The meshes to chose from are: \n"
                     "\t\t0 - cake of pyramids\n\t\t1 - hybrid mesh with all 3D elements"
                     "\n\t\t2 - a long row of bricks out ot pyramids"
                     "\n\t\t3 - user-specific mesh-file"
                     "\n\t\tdefault - a single cube. The element can be changed with option -e");
  sc_options_add_int(opt, 'n', "num_elements", &elements, 3, "The number of elements to use"
                     " if a cake is build. Has to be larger than 2.");
  sc_options_add_string(opt, 'f', "mshfile", &file, "NULL", "Prefix of the msh-file.");

  parsed = sc_options_parse(t8_get_package_id(), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if(parsed >= 0 && 0 <= level && 4 <= eclass_int && eclass_int < T8_ECLASS_COUNT &&
          elements > 2 && 0 <=mesh && mesh < 5){
      eclass = (t8_eclass_t) eclass_int;
      t8_basic_hybrid (level, endlvl, do_vtk, eclass, elements, mesh, file);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  sc_options_destroy(opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
