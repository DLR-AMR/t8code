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

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* A search function that matches all elements.
 * This function assumes that the forest user pointer is an sc_array
 * with one int for each local leaf.
 * If this function is called for a leaf, it sets the corresponding entry to 1.
 */
static int
t8_test_search_all_fn (t8_forest_t forest,
                       t8_locidx_t ltreeid,
                       const t8_element_t *element,
                       const int is_leaf,
                       t8_element_array_t *leaf_elements,
                       t8_locidx_t tree_leaf_index, void *query,
                       size_t query_index)
{
  SC_CHECK_ABORT (query == NULL,
                  "Search callback must not be called with query argument.");

  sc_array_t         *matched_leafs =
    (sc_array_t *) t8_forest_get_user_data (forest);
  if (is_leaf) {
    t8_locidx_t         tree_offset;
    t8_locidx_t         test_ltreeid;
    t8_element_t       *test_element;
    t8_eclass_t         tree_class =
      t8_forest_get_tree_class (forest, ltreeid);
    t8_eclass_scheme_c *ts;
    ts = t8_forest_get_eclass_scheme (forest, tree_class);

    tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    /* Set the corresponding entry to 1 */
    *(int *) t8_sc_array_index_locidx (matched_leafs,
                                       tree_offset + tree_leaf_index) = 1;
    /* Test whether tree_leaf_index is actually the index of the element */
    test_element =
      t8_forest_get_element (forest, tree_offset + tree_leaf_index,
                             &test_ltreeid);
    SC_CHECK_ABORT (ts->t8_element_compare (element, test_element) == 0,
                    "Element and index passed to search callback do not match.");
    SC_CHECK_ABORT (ltreeid == test_ltreeid, "Tree missmatch in search.");
  }
  return 1;
}

static int
t8_test_search_query_all_fn (t8_forest_t forest,
                             t8_locidx_t ltreeid,
                             const t8_element_t *element,
                             const int is_leaf,
                             t8_element_array_t *leaf_elements,
                             t8_locidx_t tree_leaf_index, void *query,
                             size_t query_index)
{
  /* The query callback is allways called with a query */
  SC_CHECK_ABORT (query != NULL,
                  "query callback must be called with query argument.");
  /* The query is an int with value 42 (see below) */
  SC_CHECK_ABORT (*(int *) query == 42,
                  "Wrong query argument passed to query callback.");
  /* The query index gives the position of the query in the queries array
   * of the calling search forest_search. Since there is only one query in the
   * array in this test, the index must always be 0. */
  SC_CHECK_ABORT (query_index == 0,
                  "Wrong query index passed to query callback.");
  if (is_leaf) {
    /* Test whether tree_leaf_index is actually the index of the element */
    t8_locidx_t         tree_offset;
    t8_locidx_t         test_ltreeid;
    t8_element_t       *test_element;
    t8_eclass_t         tree_class =
      t8_forest_get_tree_class (forest, ltreeid);
    t8_eclass_scheme_c *ts;
    ts = t8_forest_get_eclass_scheme (forest, tree_class);

    tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    test_element =
      t8_forest_get_element (forest, tree_offset + tree_leaf_index,
                             &test_ltreeid);
    SC_CHECK_ABORT (ts->t8_element_compare (element, test_element) == 0,
                    "Element and index passed to search callback do not match.");
    SC_CHECK_ABORT (ltreeid == test_ltreeid, "Tree missmatch in search.");
  }

  return 1;
}

static void
t8_test_search_one_query_matches_all (sc_MPI_Comm comm, t8_eclass_t eclass,
                                      int level)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  const int           query = 42;
  t8_locidx_t         ielement;
  t8_locidx_t         num_elements;
  sc_array_t          queries;
  sc_array_t          matched_leafs;

  default_scheme = t8_scheme_new_default_cxx ();
  /* Construct a cube coarse mesh */
  cmesh = t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
  /* Build a uniform forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);

  /* set up a single query containing our query */
  sc_array_init_size (&queries, sizeof (int), 1);
  *(int *) sc_array_index (&queries, 0) = query;

  num_elements = t8_forest_get_local_num_elements (forest);
  /* set up an array in which we flag whether an element was matched in the
   * search */
  sc_array_init_size (&matched_leafs, sizeof (int), num_elements);
  /* write 0 in every entry */
  for (ielement = 0; ielement < num_elements; ++ielement) {
    *(int *) t8_sc_array_index_locidx (&matched_leafs, ielement) = 0;
  }

  /* Set the array as user data so that we can access it in the search callback */
  t8_forest_set_user_data (forest, &matched_leafs);
  /* Call search. This search matches all elements. After this call we expect
   * all entries in the matched_leafs array to be set to 1. */
  t8_forest_search (forest, t8_test_search_all_fn,
                    t8_test_search_query_all_fn, &queries);

  /* Check whether matched_leafs entries are all 1 */
  for (ielement = 0; ielement < num_elements; ++ielement) {
    SC_CHECK_ABORTF (*(int *)
                     t8_sc_array_index_locidx (&matched_leafs, ielement) == 1,
                     "Search did not match all leafs. First missmatch at leaf %i.",
                     ielement);
  }

  t8_forest_unref (&forest);
  sc_array_reset (&matched_leafs);
  sc_array_reset (&queries);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;
  int                 ieclass;
  int                 ilevel;
  const int           maxlevel = 6;     /* the maximum refinement level to which we test */

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  for (ieclass = T8_ECLASS_VERTEX; ieclass < T8_ECLASS_COUNT; ieclass++) {
    for (ilevel = 0; ilevel <= maxlevel; ++ilevel) {
      t8_global_productionf
        ("Testing search that matches all with eclass %s, level %i\n",
         t8_eclass_to_string[ieclass], ilevel);
      t8_test_search_one_query_matches_all (mpic, (t8_eclass_t) ieclass,
                                            ilevel);
    }
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
