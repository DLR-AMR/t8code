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

#include <sc_options.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_adapt.h>

/* Build a cmesh according to which a later forest shall be refined. */
t8_forest_t
t8_adapt_forest_init_adapt_geometry (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh = t8_cmesh_new_from_class (T8_ECLASS_TET, comm);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  const int           level = 0;
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_SELF);
  T8_ASSERT (!t8_cmesh_is_partitioned (cmesh));
  return forest;
}

t8_forest_t
t8_adapt_cmesh_init_forest (sc_MPI_Comm comm, int level)
{
  t8_cmesh_t          cmesh =
    t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  return forest;
}

/*Idee: cmesh als user_data mitgeben. */
typedef struct
{
  t8_forest_t         forest_to_adapt_from;
} t8_adapt_cmesh_user_data_t;
typedef struct
{
  t8_locidx_t         tree_id;  /* Tree id of the element */
  t8_locidx_t         element_id;       /* Id of the current element. */
} t8_adapt_cmesh_search_query_t;

static int
t8_adapt_cmesh_search_element_callback (t8_forest_t forest,
                                        t8_locidx_t ltreeid,
                                        const t8_element_t *
                                        element,
                                        const int is_leaf,
                                        t8_element_array_t *
                                        leaf_elements,
                                        t8_locidx_t
                                        tree_leaf_index, void *query,
                                        size_t query_index)
{
  T8_ASSERT (query == NULL);
  /* We want to continue the search as long as there are queries left for the element. */
  return 1;
}

static int
t8_adapt_cmesh_search_query_callback (t8_forest_t forest,
                                      t8_locidx_t ltreeid,
                                      const t8_element_t *
                                      element,
                                      const int is_leaf,
                                      t8_element_array_t *
                                      leaf_elements,
                                      t8_locidx_t
                                      tree_leaf_index, void *query,
                                      size_t query_index)
{
  /*Identifiziere Element in forest, die einen Mittelpunkt eines
     Elements in  cmesh_to_adapt_from enthalten. 
     Falls ja, verfeinere und suche weiter (bis maxlevel)
     Falls nein -> fertig */
  t8_adapt_cmesh_user_data_t *user_data =
    (t8_adapt_cmesh_user_data_t *) t8_forest_get_user_data (forest);
  const t8_forest_t    forest_to_adapt_from = user_data->forest_to_adapt_from;
  t8_adapt_cmesh_search_query_t *search_query =
    (t8_adapt_cmesh_search_query_t *) query;
  const t8_locidx_t   forest_to_adapt_from_tree_id = search_query->tree_id;
  const t8_locidx_t   forest_to_adapt_from_element_id = search_query->element_id;

  const double      tolerance = 0.2;

  /* TODO: Get tree id and element id from query */


  /* TODO: Compute midpoint of query. Return true only if midpoint is in current element. */
  /* Compute midpoint of tree */
  double              midpoint[3];
  t8_element_t *      forest_to_adapt_from_element = t8_forest_get_element_in_tree(forest_to_adapt_from, 
                                                    forest_to_adapt_from_tree_id,
                                                    forest_to_adapt_from_element_id);

  t8_forest_element_centroid (forest_to_adapt_from, 
                              forest_to_adapt_from_tree_id,
                              forest_to_adapt_from_element,
                              midpoint);

  /* Check if midpoint is inside element */
  return t8_forest_element_point_inside(forest, ltreeid, element, midpoint, tolerance);

}

static void
t8_adapt_cmesh_search (t8_forest_t forest, t8_forest_t forest_to_adapt_from,
                       sc_array_t * markers)
{
  sc_array_t          search_queries;
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest_to_adapt_from);
  const t8_locidx_t   num_trees =
    t8_forest_get_num_local_trees (forest_to_adapt_from);
  /* Initialize the queries array */
  sc_array_init_count (&search_queries,
                       sizeof (t8_adapt_cmesh_search_query_t), num_elements);
  /* Fill queries with correct ids. */

#if 0
  /* TODO Fill tree ids and element ids */
  for (t8_locidx_t itree; itree < num_trees; ++itree) {
    *(t8_sc_array_index_locidx (&search_queries, itree)) = itree;
  }
#endif
  /* Set the cmesh as forest user data */
  t8_adapt_cmesh_user_data_t search_user_data;
  search_user_data.forest_to_adapt_from = forest_to_adapt_from;
  t8_forest_set_user_data (forest, &search_user_data);
  /* Fill marker array.
   * elements that should be refined are set to 1. 0 for no refinemnet. -1 for coarsening. */
  t8_forest_search (forest, t8_adapt_cmesh_search_element_callback,
                    t8_adapt_cmesh_search_query_callback, &search_queries);
}

static              t8_forest_t
t8_adapt_cmesh_adapt_forest (t8_forest_t forest, t8_forest_t forest_to_adapt_from)
{
  /* TODO ... */

  const t8_locidx_t   num_local_elements =
    t8_forest_get_local_num_elements (forest);
  /* Create marker array  to mark elements for refinement */
  sc_array_t          marker_array;

  sc_array_init_count (&marker_array, sizeof (short), num_local_elements);

  t8_adapt_cmesh_search (forest, forest_to_adapt_from, &marker_array);

  /* Adapt the forest according to the markers */
  t8_forest_t         forest_adapt =
    t8_forest_new_adapt (forest, t8_forest_adapt_marker_array_callback,
                         0, 0, &marker_array);
}



int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  int                 helpme;
  int                 parsed;
  int                 level;
  const sc_MPI_Comm   comm = sc_MPI_COMM_WORLD;

  /* long help message */
  snprintf (help, BUFSIZ,
            "This program adapts a forest according to a provided forest mesh.\n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
#ifdef T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");

  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The minimum refinement level of the forest.");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 /* TODO: Check correct arguments */ ) {
    t8_forest_t         forest_to_adapt_from =
      t8_adapt_forest_init_adapt_geometry (comm);
    t8_forest_t         forest = t8_adapt_cmesh_init_forest (comm, level);

    /* Identifiziere Element in forest, die einen Mittelpunkt eines
       Elements in  cmesh_to_adapt_from enthalten. */
    /*[D] Wäre nicht besser: Identifiziere Elemente in forest, 
       deren Mittelpunkt in einem Element des cmesh_to_adapt_from liegen? Im Grobgitter 
       beschreibt, der Mittelpunkt das Element sehr schlecht. */
    /* Adaptiere forest, so dass alle identifizierten Elemente verfeienrt werden. */
    t8_forest_unref (&forest_to_adapt_from);
    t8_forest_unref (&forest);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
