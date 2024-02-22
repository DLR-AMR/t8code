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

#include <sc_refcount.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

void
t8_timings_new (int level, int eclass_int)
{
  t8_forest_t forest;
  t8_eclass_t eclass;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  T8_ASSERT (level >= 0);

  eclass = eclass_int;

  t8_global_productionf ("=P= Starting forest_new\n");

  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0), sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, level);

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_forest_commit (forest);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");

  t8_global_productionf ("=P= Done forest_new.\n");

  t8_forest_unref (&forest);

  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

int
main (int argc, char **argv)
{
  int mpiret, mpisize;
  int start_level, eclass_int;
  int first_argc;
  sc_options_t *opt;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 's', "slevel", &start_level, 0, "initial refine level");
  sc_options_add_int (opt, 'e', "eclass", &eclass_int, (int) T8_ECLASS_QUAD, "eclass");

  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);

  if (first_argc < 0 || first_argc != argc || 0 > eclass_int || (int) T8_ECLASS_COUNT <= eclass_int) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }

  t8_timings_new (start_level, eclass_int);

  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
