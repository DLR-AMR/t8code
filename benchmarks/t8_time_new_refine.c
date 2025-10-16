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
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_scheme.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8_forest/t8_forest_types.h> /* TODO: This file should not be included from an application */
/* This function refines every element */
static int
t8_basic_adapt_refine (__attribute__ ((unused)) t8_forest_t forest, __attribute__ ((unused)) t8_forest_t forest_from,
                       __attribute__ ((unused)) t8_locidx_t which_tree, t8_eclass_t tree_class,
                       __attribute__ ((unused)) t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                       __attribute__ ((unused)) const int is_family, const int num_elements, t8_element_t *elements[])
{
  const int level = t8_element_get_level (scheme, tree_class, elements[0]);
  /* coarsen */
  if (num_elements > 1) {
    if (level > 0)
      return -1;
    return 0;
  }
  return 1;
}

/* This function coarsens each element */
static int
t8_basic_adapt_coarsen (__attribute__ ((unused)) t8_forest_t forest, __attribute__ ((unused)) t8_forest_t forest_from,
                        __attribute__ ((unused)) t8_locidx_t which_tree,
                        __attribute__ ((unused)) t8_eclass_t tree_class,
                        __attribute__ ((unused)) t8_locidx_t lelement_id,
                        __attribute__ ((unused)) const t8_scheme_c *scheme, const int is_family,
                        __attribute__ ((unused)) int num_elements, __attribute__ ((unused)) t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

static void
t8_timings_adapt (int start_l, int end_l, int runs, int dim)
{
  t8_forest_t *forests;
  int li, num_levels, cur_for, run;
  t8_eclass_t eclass;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  num_levels = end_l - start_l + 1;
  T8_ASSERT (num_levels > 0);
  T8_ASSERT (runs > 0);

  forests = T8_ALLOC (t8_forest_t, 2 * num_levels * runs);

  t8_forest_init (&forests[0]);

  eclass = dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET;
  /*
  t8_forest_set_cmesh (forests[0],
                       t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0), 0);
*/
  t8_forest_set_cmesh (forests[0], t8_cmesh_new_bigmesh (eclass, 512, sc_MPI_COMM_WORLD), sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forests[0], t8_scheme_new_default ());
  t8_forest_set_level (forests[0], start_l);
  t8_forest_commit (forests[0]);

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  for (run = 0, cur_for = 1; run < runs; run++) {
    for (li = 1; li < num_levels; li++, cur_for++) {
      t8_forest_init (&forests[cur_for]);
      t8_forest_set_adapt (forests[cur_for], forests[cur_for - 1], t8_basic_adapt_refine, 0);
      t8_forest_commit (forests[cur_for]);
    }
    for (li = 1; li < num_levels; li++, cur_for++) {
      t8_forest_init (&forests[cur_for]);
      t8_forest_set_adapt (forests[cur_for], forests[cur_for - 1], t8_basic_adapt_coarsen, 0);
      t8_forest_commit (forests[cur_for]);
    }
  }

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "Adapt");

  t8_forest_unref (&forests[cur_for - 1]);
  T8_FREE (forests);

  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

void
t8_timings_new (int level, int dim)
{
  t8_forest_t forest;
  t8_eclass_t eclass;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  T8_ASSERT (level >= 0);
  T8_ASSERT (dim == 2 || dim == 3);

  eclass = dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET;

  t8_global_productionf ("=P= Starting forest_new with %.0f elements.\n", 512 * pow (2, dim * level));

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0, 0), sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  t8_forest_set_level (forest, level);
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
  int start_level, end_level, dim;
  int first_argc;
  int use_refine, use_new;
  sc_options_t *opt;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 's', "slevel", &start_level, 0, "initial refine level");
  sc_options_add_int (opt, 'e', "elevel", &end_level, 0,
                      "Final refine level: greater or equal to initial refine level");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "dimension: 2 or 3");
  sc_options_add_switch (opt, 'r', "refine", &use_refine,
                         "Time refining from start_level to end_level - this is default. ");
  sc_options_add_switch (opt, 'n', "new", &use_new,
                         "Time new with start_level - If this is given -r has to be switched on by hand if desired.");

  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);

  if (end_level < start_level) {
    end_level = start_level;
  }

  if (first_argc < 0 || first_argc != argc || 2 > dim || dim > 3 || (use_refine && end_level < start_level)) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }

  if (use_refine || !use_new) {
    t8_timings_adapt (start_level, end_level, 1, dim);
  }
  if (use_new) {
    t8_timings_new (start_level, dim);
  }

  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
