/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
#include <t8_default.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

/* This function refines every element */
static int
t8_basic_adapt_refine (t8_forest_t forest, t8_topidx_t which_tree,
                t8_eclass_scheme_t * ts,
                int num_elements, t8_element_t * elements[])
{
#if 0
  int                 level;
#endif
  T8_ASSERT (num_elements == 1 || num_elements ==
             t8_eclass_num_children[ts->eclass]);
#if 0
  level = t8_element_level (ts, elements[0]);
  /* coarsen */
  if (num_elements > 1) {
    if (level > 0)
      return -1;
    return 0;
  }
#endif
  return 1;
}

/* This function coarsens each element */
static int
t8_basic_adapt_coarsen (t8_forest_t forest, t8_topidx_t which_tree,
                t8_eclass_scheme_t * ts,
                int num_elements, t8_element_t * elements[])
{
  if (num_elements > 1) {
    return -1;
  }
  return 0;
}

static void
t8_timings_adapt (int start_l, int end_l, int runs, int dim)
{
  t8_forest_t        *forests;
  int                 li, num_levels, cur_for, run;
  t8_eclass_t         eclass;

  num_levels = end_l - start_l + 1;
  T8_ASSERT (num_levels > 0);
  T8_ASSERT (runs > 0);
  if (num_levels > 0) {
    forests = T8_ALLOC (t8_forest_t, num_levels * runs);
  }
  t8_forest_init (&forests[0]);

  eclass = dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET;
/*
  t8_forest_set_cmesh (forests[0],
                       t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0));
*/
  t8_forest_set_cmesh (forests[0],
                       t8_cmesh_new_bigmesh (eclass, 512, sc_MPI_COMM_WORLD, 0));
  t8_forest_set_scheme (forests[0], t8_scheme_new_default ());
  t8_forest_set_level (forests[0], start_l);
  t8_forest_commit (forests[0]);

  for (run = 0, cur_for = 0;run < runs;run++) {
    for (li = 1; li < num_levels; li++, cur_for++) {
      t8_forest_init (&forests[cur_for]);
      t8_forest_set_adapt_temp (forests[cur_for], forests[cur_for - 1], t8_basic_adapt_refine,
                                NULL, 0);
      t8_forest_commit (forests[cur_for]);
    }
    for (li = 1; li < num_levels; li++, cur_for++) {
      t8_forest_init (&forests[cur_for]);
      t8_forest_set_adapt_temp (forests[cur_for], forests[cur_for - 1], t8_basic_adapt_coarsen,
                                NULL, 0);
      t8_forest_commit (forests[cur_for]);
    }
  }

  t8_forest_unref (&forests[num_levels - 1]);
  T8_FREE (forests);
}

int
main (int argc, char **argv)
{
  int                 mpiret, mpisize;
  int                 start_level, end_level, dim;
  int                 first_argc;
  int		      repeat;
  sc_options_t       *opt;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[1];

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);
 
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 's', "slevel", &start_level, 0,
                      "initial refine level");
  sc_options_add_int (opt, 'e', "elevel", &end_level, 0,
                      "Final refine level: greater or equal to initial refine level");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "dimension: 2 or 3");
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc
      || 2 > dim || dim > 3 || end_level < start_level) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_timings_adapt (start_level, end_level, 20, dim);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "Adapt");

  t8_global_productionf ("Timings for NP=%i\n", mpisize);
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);

  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
