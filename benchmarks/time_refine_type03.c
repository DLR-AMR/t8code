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
#include <t8_schemes/t8_default/t8_dtri.h>
#include <t8_schemes/t8_default/t8_dtet.h>
#include <t8_schemes/t8_default.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_general.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8_forest/t8_forest_types.h> /* TODO: This file should not be included from an application */

int max_ref_level = 0;

/* This function refines every element */
static int
t8_basic_adapt_refine_type (t8_forest_t forest, t8_locidx_t which_tree, t8_eclass_t tree_class, const t8_scheme_c *ts,
                            const int is_family, int num_elements, t8_element_t *elements[])
{
  int level;
  int type;
  int dim;

  T8_ASSERT (!is_family || num_elements == t8_eclass_num_children[tree_class]);

  dim = t8_eclass_to_dimension[tree_class];
  level = t8_element_get_level (forest, tree_class, elements[0]);
  if (level >= max_ref_level) {
    return 0;
  }
  /* get the type of the current element */
  type = dim == 2 ? ((t8_dtri_t *) elements[0])->type : ((t8_dtet_t *) elements[0])->type;
  /* refine type 0 and 3 */
  if (type == 0 || type == 3) {
    return 1;
  }
  return 0;
}

static void
t8_timings_adapt_type (int start_l, int dim)
{
  t8_forest_t forests[2];
  t8_eclass_t eclass;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];
  long long num_el;

  t8_forest_init (&forests[0]);

  eclass = dim == 2 ? T8_ECLASS_TRIANGLE : T8_ECLASS_TET;

  t8_forest_set_cmesh (forests[0], t8_cmesh_new_bigmesh (eclass, 512, sc_MPI_COMM_WORLD, 0), sc_MPI_COMM_WORLD);
  t8_forest_set_scheme (forests[0], t8_scheme_new_default ());
  t8_forest_set_level (forests[0], start_l);
  t8_forest_commit (forests[0]);

  t8_forest_init (&forests[1]);
  t8_forest_set_adapt (forests[1], forests[0], t8_basic_adapt_refine_type, 1);

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_forest_commit (forests[1]);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");

  num_el = (long long) forests[1]->local_num_elements;

  t8_forest_unref (&forests[1]);

  t8_debugf ("=P= I have %lli elements\n", num_el);
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
}

int
main (int argc, char **argv)
{
  int mpiret, mpisize;
  int start_level, end_level, dim;
  int first_argc;
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
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  if (first_argc < 0 || first_argc != argc || 2 > dim || dim > 3 || end_level < start_level) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }

  max_ref_level = end_level;
  t8_timings_adapt_type (start_level, dim);

  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
