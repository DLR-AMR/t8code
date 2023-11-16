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

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h> 
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>
#include <t8_cmesh_vtk_writer.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>


int
t8_adapt_callback (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts,
                         const int is_family,
                         const int num_elements, t8_element_t *elements[])
{
    return 1;
}


void
t8_time_new_uniform (t8_cmesh_t cmesh, t8_scheme_cxx_t *scheme, int level, sc_MPI_Comm comm)
{
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[1];
  t8_forest_t         forest;

  /* Start timer */
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  /* commit (= partition) the new cmesh */
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
  /* measure passed time */
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New uniform");
  /* print stats */
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);

  /* cleanup */
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level;
  sc_MPI_Comm         comm;
  sc_options_t        *opt;
  t8_cmesh_t          cmesh;
  int         eclass;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_STATISTICS);
  comm = sc_MPI_COMM_WORLD;

  /** OPTIONS */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'l', "level", &level, 1,
                      "Refinement level.");
  sc_options_add_int (opt, 'e', "eclass", &eclass, T8_ECLASS_HEX,
                      "Element Class");
  /* parse command line options */
  int first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                                 opt, argc, argv);
  /* check for wrong usage of arguments */
  if (first_argc < 0 || first_argc != argc || level < 0) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_destroy (opt);


  cmesh = t8_cmesh_new_from_class ((t8_eclass_t) eclass, comm);
//  cmesh = t8_cmesh_new_multi_class_3D (comm);
  t8_cmesh_ref(cmesh);

  t8_scheme_cxx_t    *scheme_standalone = t8_scheme_new_standalone_cxx ();
  t8_scheme_cxx_t    *scheme_default = t8_scheme_new_default_cxx ();
  t8_time_new_uniform(cmesh, scheme_standalone, level, comm);
  t8_time_new_uniform(cmesh, scheme_default, level, comm);


  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
