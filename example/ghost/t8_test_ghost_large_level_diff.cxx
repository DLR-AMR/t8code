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

#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk.h>
#include <t8_default/t8_dprism.h>
#include <t8_default/t8_dprism_bits.h>
#include <t8_default/t8_dtri.h>
#include <t8_default/t8_dtet.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_default_cxx.hxx>

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

static int
t8_ghost_fractal_adapt (t8_forest_t forest, t8_forest_t forest_from,
                            t8_locidx_t which_tree, t8_locidx_t lelement_id,
                            t8_eclass_scheme_c * ts, int num_elements,
                            t8_element_t * elements[])
{
  int                   level;
  int                   type, child_id;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));

  level = ts->t8_element_level (elements[0]);
  if (level >= *(int *) t8_forest_get_user_data (forest)) {
    return 0;
  }
  if(ts->eclass == T8_ECLASS_PRISM){
      type = ((t8_dprism_t *) elements[0])->tri.type;
      /* refine type 0 */
      if (type == 0) {
        child_id = (((t8_dprism_t *) elements[0])->line.level == 0) ?
                    (which_tree % T8_DPRISM_CHILDREN):
                    t8_dprism_child_id((t8_dprism_t *) elements[0]);
        if(child_id == 3 || child_id == 4){
            return 0;
        }
        else{
            return 1;
        }
      }
      return 0;
  }
  else if(ts->eclass == T8_ECLASS_TET){
      type = ((t8_dtet_t * ) elements[0])->type;
      if (type == 0 || type == 3 || type == 5){
          return 1;
      }
      return 0;
  }
  else if(ts->eclass == T8_ECLASS_HEX){
      child_id = (((p4est_quadrant_t *)elements[0])->level == 0) ? (which_tree % P4EST_CHILDREN)
                              : p4est_quadrant_child_id((p4est_quadrant_t *)elements[0]);
      if(child_id == 0 || child_id == 3 || child_id == 5 || child_id == 6){
          return 1;
      };
      return 0;
  }

}

static void
t8_ghost_large_level_diff(const char* path, int dim, int level, int refine,
                          int no_vtk, sc_MPI_Comm comm){
    t8_forest_t     forest, forest_adapt, forest_partition;
    t8_cmesh_t      cmesh, cmesh_partition;
    sc_flopinfo_t       fi, snapshot;
    sc_statinfo_t       stats[1];

    //T8_ASSERT(path == NULL);
    //get cmesh
    //cmesh = t8_cmesh_from_msh_file((char *) path, 1, comm, dim, 0);
    cmesh = t8_cmesh_new_hypercube(T8_ECLASS_PRISM, comm, 0, 0, 0);
    t8_cmesh_init(&cmesh_partition);
    t8_cmesh_set_derive(cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform(cmesh_partition, level);
    t8_cmesh_commit(cmesh_partition, comm);
    if(!no_vtk){
        t8_cmesh_vtk_write_file(cmesh_partition, "partitioned_cmesh", 1.0);
    }

    //New
    t8_forest_init(&forest);
    t8_forest_set_cmesh(forest, cmesh_partition, comm);
    t8_forest_set_scheme(forest, t8_scheme_new_default_cxx() );
    t8_forest_set_level(forest, level);
    sc_flops_start (&fi);
    sc_flops_snap (&fi, &snapshot);
    t8_forest_commit (forest);
    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[0], snapshot.iwtime, "New");
    if(!no_vtk){
        t8_forest_write_vtk(forest, "Forest_new");
    }
    t8_global_productionf("Successfully commited forest.\n");

    //Adapt
    t8_forest_init(&forest_adapt);
    refine += level;
    t8_forest_set_user_data(forest_adapt, &refine);
    t8_forest_set_profiling(forest_adapt, 1);
    t8_forest_set_adapt(forest_adapt, forest, t8_ghost_fractal_adapt, 1);
    t8_forest_commit(forest_adapt);
    if(!no_vtk){
        t8_forest_write_vtk(forest_adapt, "Forest_adapt");
    }
    t8_global_productionf("Successfully refined forest adaptivly.\n");

    //Partition
    t8_forest_init(&forest_partition);
    t8_forest_set_partition(forest_partition, forest_adapt, 0);
    t8_forest_set_ghost_ext(forest_partition, 1, T8_GHOST_FACES, 3);
    t8_forest_set_profiling(forest_partition, 1);
    t8_forest_commit(forest_partition);
    if(!no_vtk){
        t8_forest_write_vtk(forest_partition, "Forest_partition");
    }
    t8_global_productionf("Successfully partitioned forest.\n");
    t8_forest_ghost_print(forest_partition);

    t8_forest_print_profile (forest_partition);
    t8_forest_unref (&forest_partition);

    sc_stats_compute (comm, 1, stats);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);


}

int
main (int argc, char *argv[])
{
  int                 mpiret, parsed, dim, level, refine, mpisize, helpme, no_vtk;
  sc_options_t       *opt;
  const char         *path;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  t8_cmesh_t          cmesh;

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>",
            basename (argv[0]));
  snprintf (help, BUFSIZ,
            "This program reads a .msh file "
            "created by the GMSH program and constructs a "
            "t8code coarse mesh from them.\n\n%s\n\nExample: %s -f A1\nTo open the file A1.ms."
            "\n\nThe default dimension of the mesh to read is 2. Since the "
            ".msh format stores elements of all (lower) dimensions "
            "the user must provide the argument for a different dimension by hand, if "
            "desired.\n", usage, basename (argv[0]));

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  opt = sc_options_new (argv[0]);
  sc_options_add_switch(opt, 'h', "help", &helpme, "Display a short help message:");
  sc_options_add_string (opt, 'p', "path", &path, "", "The path to the mesh-file. "
                         "The file must end in .msh and be created with gmsh.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "The dimension of the mesh.");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The intial refinement level of the mesh.");
  sc_options_add_int (opt, 'r', "refine", &refine, 0, "The number of levels that the forest "
                       "is refined from the initial level.");
  sc_options_add_switch(opt, 'o', "no-vtk", &no_vtk, "Suppress vtk output.");
  parsed =sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if(helpme){
      /* display help message and usage */
      t8_global_productionf ("%s\n", help);
      sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && 0 <= refine) {
    t8_ghost_large_level_diff(path, dim, level, refine, no_vtk,sc_MPI_COMM_WORLD);
    return 1;
  }
  else {
    /*wrong usage*/
      t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
      sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    /*
    cmesh = t8_read_msh_file_build_cmesh (path, partition, dim, master);
    t8_cmesh_destroy (&cmesh);
    sc_options_print_summary (t8_get_package_id (), SC_LP_PRODUCTION, opt);
    */
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
