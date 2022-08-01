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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <sc_functions.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* *INDENT-OFF* */
class linear_id:public testing::TestWithParam <t8_eclass > {
protected:
  void SetUp () override {
    eclass = GetParam ();
    scheme = t8_scheme_new_default_cxx ();
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test);
    ts->t8_element_set_linear_id (element, 0, 0);
  }
  void TearDown ()override {
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &test);
    t8_scheme_cxx_unref (&scheme);
  }
  t8_element_t *element, *child, *test;
  t8_scheme_cxx * scheme; t8_eclass_scheme_c *ts;
  t8_eclass_t eclass; sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

/* *INDENT-ON* */
static int
t8_test_init_linear_id_refine_everything (t8_forest_t forest,
                                          t8_forest_t forest_from,
                                          t8_locidx_t which_tree,
                                          t8_locidx_t lelement_id,
                                          t8_eclass_scheme_c *ts,
                                          const int is_family,
                                          const int num_elements,
                                          t8_element_t *elements[])
{
  return 1;
}

/* *INDENT-OFF* */
TEST_P (linear_id, uniform_forest) {
  t8_forest_t forest, forest_adapt;
  t8_cmesh_t cmesh;
  t8_locidx_t first_id, last_id, j, id;
  t8_locidx_t first_tid, last_tid, tree_id; t8_element_t *element;
  int i;
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 5;
#else
  const int maxlvl = 6;
#endif
  cmesh = t8_cmesh_new_from_class (ts->eclass, comm);
  t8_cmesh_ref (cmesh);
  forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  t8_scheme_cxx_ref (scheme);
  for (i = 0; i < maxlvl; i++) {
    /*Get the id of the first and last tree on this process */
    first_tid = t8_forest_get_first_local_tree_id (forest);
    last_tid = first_tid + t8_forest_get_num_local_trees (forest);
    /*Iterate over trees */
    for (tree_id = first_id; tree_id < last_tid; tree_id++) {
      /*Get id of the first and last element on this tree */
      first_id = t8_forest_get_first_local_element_id (forest);
      last_id = first_id + t8_forest_get_local_num_elements (forest);
      /*Iterate over elements */
      for (j = first_id; j < last_id; j++) {
        /*Get the j-th element and check the computed linear id */
        element = t8_forest_get_element (forest, j, &tree_id);
        /*Get the ID of the element at current level */
        id = ts->t8_element_get_linear_id (element, i);
        EXPECT_EQ (id, j);
      }
    }
    t8_forest_init (&forest_adapt);
    t8_forest_set_level (forest_adapt, i + 1);
    t8_forest_set_adapt (forest_adapt, forest,
                          t8_test_init_linear_id_refine_everything,
                          0); t8_forest_commit (forest_adapt);
    forest = forest_adapt;
  }
  t8_cmesh_unref (&cmesh); 
  t8_forest_unref (&forest_adapt);
}

/* *INDENT-ON* */

#if 0
/*Check, if all descendants of an element at level maxlvl have the same id on
 * the level of the input element as the input element*/
static void
t8_id_at_other_lvl_check (t8_element_t *element,
                          t8_element * child,
                          t8_eclass_scheme_c *ts, int maxlvl)
{
  int                 level = ts->t8_element_level (element);
  t8_linearidx_t      current_id =
    ts->t8_element_get_linear_id (element, level);
  t8_linearidx_t      num_descendants =
    ts->t8_element_count_leafs (element, level);
  t8_linearidx_t      id_at_lvl =
    ts->t8_element_get_linear_id (element, maxlvl);
  t8_linearidx_t      i;
  t8_linearidx_t      id;
  for (i = 0; i < num_descendants; i++) {
    ts->t8_element_set_linear_id (child, maxlvl, id_at_lvl + i);
    id = ts->t8_element_get_linear_id (child, level);
    SC_CHECK_ABORT (id == current_id, "Wrong id");
  }
}

static void
t8_check_linear_id (const int maxlvl)
{
  t8_element_t       *element, *child, *test;
  t8_scheme_cxx_t    *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  int                 level;
  t8_linearidx_t      num_desc;
  t8_linearidx_t      id;
  t8_eclass_t         eclass;
  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_COUNT; eclassi++) {
    t8_productionf ("Testing linear_id for eclass %s\n",
                    t8_eclass_to_string[eclassi]);
    eclass = (t8_eclass_t) eclassi;
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
    t8_check_uniform_forest (ts, scheme, sc_MPI_COMM_WORLD, maxlvl);
    /* Get element and initialize it */
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test);
    ts->t8_element_set_linear_id (element, 0, 0);
    for (level = 0; level < maxlvl; level++) {
      num_desc = ts->t8_element_count_leafs_from_root (level);
      for (id = 0; id < num_desc; id++) {
        ts->t8_element_set_linear_id (child, level, id);
        t8_id_at_other_lvl_check (child, test, ts, maxlvl);
      }
    }
    /* Destroy element */
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &test);
  }

  /* Destroy scheme */
  t8_scheme_cxx_unref (&scheme);
}
#endif

INSTANTIATE_TEST_SUITE_P (t8_test_init_linear_id, linear_id,
                          testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));
#if 0
int
main (int argc, char **argv)
{
  int                 mpiret;
#ifdef T8_ENABLE_DEBUG
  const int           maxlvl = 3;
#else
  const int           maxlvl = 4;
#endif
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);
  t8_check_linear_id (maxlvl);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
#endif
