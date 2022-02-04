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
#include <t8_cmesh_readmshfile.h>
#include <sc_statistics.h>

/* Build a cmesh according to which a later forest shall be refined. */
t8_forest_t
t8_adapt_forest_init_adapt_geometry (sc_MPI_Comm comm, const char *meshfile,
                                     const int dim, const int level)
{
  t8_cmesh_t          cmesh;
  if (meshfile != NULL) {
    cmesh = t8_cmesh_from_msh_file (meshfile, 0, sc_MPI_COMM_SELF, dim, 0);
  }
  else {
    //cmesh = t8_cmesh_new_from_class (T8_ECLASS_PRISM, sc_MPI_COMM_SELF);
    cmesh = t8_cmesh_new_line_zigzag (sc_MPI_COMM_SELF);
  }
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_SELF);
  T8_ASSERT (!t8_cmesh_is_partitioned (cmesh));
  return forest;
}

t8_forest_t
t8_adapt_cmesh_init_forest (sc_MPI_Comm comm, const int level,
                            const double scale[3],
                            const double displacement[3],
                            const char *meshfile)
{

  t8_cmesh_t          cmesh;
  if (meshfile != NULL) {
    cmesh = t8_cmesh_from_msh_file (meshfile, 0, comm, 3, 0);
    t8_cmesh_scale (cmesh, scale);
    t8_cmesh_translate (cmesh, displacement);
  }
  else {
    cmesh =
      t8_cmesh_new_hypercube_ext (T8_ECLASS_PRISM, comm, 0, 0, 0, scale,
                                  displacement);
  }

  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();

  t8_cmesh_t          cmesh_partition;
  t8_cmesh_init (&cmesh_partition);
  t8_scheme_cxx_ref (scheme);
  t8_cmesh_set_partition_uniform (cmesh_partition, level, scheme);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_commit (cmesh_partition, comm);

  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh_partition, scheme, level, 1, comm);

  return forest;
}

/*Idee: cmesh als user_data mitgeben. */
typedef struct
{
  t8_forest_t         forest_to_adapt_from;
  sc_array_t         *refinement_markers;
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
  sc_array_t         *refinement_markers = user_data->refinement_markers;
  const t8_forest_t   forest_to_adapt_from = user_data->forest_to_adapt_from;
  const t8_adapt_cmesh_search_query_t *search_query =
    (const t8_adapt_cmesh_search_query_t *) query;
  const t8_locidx_t   forest_to_adapt_from_tree_id = search_query->tree_id;
  const t8_locidx_t   forest_to_adapt_from_element_id =
    search_query->element_id;
  int                 query_is_in_element = 0;
  const double        tolerance = 1e-5;

  /* TODO: Get tree id and element id from query */

  /* TODO: Compute midpoint of query. Return true only if midpoint is in current element. */
  /* Compute midpoint of tree */
  double              midpoint[3];
  const t8_element_t *forest_to_adapt_from_element =
    t8_forest_get_element_in_tree (forest_to_adapt_from,
                                   forest_to_adapt_from_tree_id,
                                   forest_to_adapt_from_element_id);

  t8_forest_element_centroid (forest_to_adapt_from,
                              forest_to_adapt_from_tree_id,
                              forest_to_adapt_from_element, midpoint);

  /* Check if midpoint is inside element */

  const t8_eclass_t   adapt_from_eclass =
    t8_forest_get_tree_class (forest_to_adapt_from,
                              forest_to_adapt_from_tree_id);
  t8_eclass_scheme_c *scheme_to_adapt_from =
    t8_forest_get_eclass_scheme (forest_to_adapt_from, adapt_from_eclass);
  const t8_element_shape_t shape_to_adapt_from =
    scheme_to_adapt_from->t8_element_shape (element);

  if (shape_to_adapt_from == T8_ECLASS_LINE) {
    /* The element class of the forest to adapt from is line and the element's class is a hex
     * we check whether the line cuts our hex element. */
    /* Note that this only works, if the hex is axis aligned. */
    query_is_in_element =
      /* t8_forest_line_cuts_aligned_hex (forest_to_adapt_from,
         forest_to_adapt_from_tree_id,
         forest_to_adapt_from_element, forest,
         ltreeid, element); */
      t8_forest_line_cuts (forest_to_adapt_from, forest_to_adapt_from_tree_id,
                           forest_to_adapt_from_element, forest,
                           ltreeid, element);
  }
  else {
    query_is_in_element =
      t8_forest_element_point_inside (forest, ltreeid, element, midpoint,
                                      tolerance);
  }

  if (!query_is_in_element) {
    /* The query is not detected in the element. 
     * remove query from search */
    return 0;
  }

  if (is_leaf) {
    /* This element is a leaf in the searched forest.
     * Hence, we mark it for refinement. */
    t8_locidx_t         element_index =
      t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
    *(short *) t8_sc_array_index_locidx (refinement_markers,
                                         element_index) = 1;

  }
  /* Keep this query in the search */
  return 1;
}

/* Set new size of markers array and set all markers to 0. */
static void
t8_adapt_template_update_markers (t8_forest_t forest, sc_array_t * markers)
{
  const t8_locidx_t   num_local_elements =
    t8_forest_get_local_num_elements (forest);

  T8_ASSERT (markers->elem_size == sizeof (short));

  sc_array_resize (markers, num_local_elements);
  for (t8_locidx_t ielement = 0; ielement < num_local_elements; ++ielement) {
    *(short *) t8_sc_array_index_locidx (markers, ielement) = 0;
  }
}

/*Compute the number of elements according to their shape*/
void
t8_adapt_cmesh_element_count (t8_forest_t forest,
                              t8_gloidx_t * element_of_class)
{
  const t8_locidx_t   num_local_trees =
    t8_forest_get_num_local_trees (forest);
  /*Iterate over all local trees */
  for (t8_locidx_t itree_id = 0; itree_id < num_local_trees; ++itree_id) {
    const t8_eclass_t   tree_class =
      t8_forest_get_tree_class (forest, itree_id);
    if (tree_class != T8_ECLASS_PYRAMID) {
      /*if the tree is not a pyramid, all elements have the same shape */
      element_of_class[tree_class] +=
        t8_forest_get_tree_num_elements (forest, itree_id);
    }
    else {
      t8_eclass_scheme_c *pyra_scheme =
        t8_forest_get_eclass_scheme (forest, tree_class);
      t8_locidx_t         num_elems_in_tree =
        t8_forest_get_tree_num_elements (forest, itree_id);
      /* Iterate over all elements and increase the counter according to the
       * shape of the element*/
      for (t8_locidx_t ielem = 0; ielem < num_elems_in_tree; ++ielem) {
        t8_element_t       *element =
          t8_forest_get_element_in_tree (forest, itree_id, ielem);
        t8_eclass_t         element_shape =
          pyra_scheme->t8_element_shape (element);
        element_of_class[element_shape]++;
      }
    }
  }
}

static              t8_forest_t
t8_adapt_cmesh_adapt_forest (t8_forest_t forest,
                             t8_forest_t forest_to_adapt_from,
                             int num_refinement_steps, int balance)
{
  sc_array_t          search_queries;
  sc_statinfo_t       total_times[2];
  double              non_search_time_total = 0, search_time_total = 0;
  t8_gloidx_t         element_of_class[T8_ECLASS_COUNT] =
    { 0, 0, 0, 0, 0, 0, 0, 0 };
  t8_gloidx_t         recv_buff[T8_ECLASS_COUNT] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  int                 mpiret, rank;

  sc_stats_init (&total_times[0], "non-search-total");
  sc_stats_init (&total_times[1], "search-total");
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (forest_to_adapt_from);
  const t8_locidx_t   num_trees =
    t8_forest_get_num_local_trees (forest_to_adapt_from);
  /* Initialize the queries array */
  sc_array_init_count (&search_queries,
                       sizeof (t8_adapt_cmesh_search_query_t), num_elements);
  /* Fill queries with correct ids. */

  for (t8_locidx_t itree = 0, iquery = 0; itree < num_trees; ++itree) {
    const t8_locidx_t   num_elements_in_tree =
      t8_forest_get_tree_num_elements (forest_to_adapt_from, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree;
         ++ielement, ++iquery) {
      t8_adapt_cmesh_search_query_t *query = (t8_adapt_cmesh_search_query_t *)
        t8_sc_array_index_locidx (&search_queries, iquery);
      query->tree_id = itree;
      query->element_id = ielement;
    }
  }

  /* Set the cmesh as forest user data */
  t8_adapt_cmesh_user_data_t search_user_data;
  search_user_data.forest_to_adapt_from = forest_to_adapt_from;

  sc_array_t          markers;
  sc_array_init (&markers, sizeof (short));

  for (int refinement_step = 0; refinement_step < num_refinement_steps;
       ++refinement_step) {
    sc_statinfo_t       times[2];
    double              non_search_time, search_time;

    sc_stats_init (&times[0], "non-search");
    sc_stats_init (&times[1], "search");

    t8_adapt_template_update_markers (forest, &markers);
    search_user_data.refinement_markers = &markers;
    t8_forest_set_user_data (forest, &search_user_data);
    /* Fill marker array.
     * elements that should be refined are set to 1. 0 for no refinemnet. -1 for coarsening. */
    search_time = -sc_MPI_Wtime ();
    t8_forest_search (forest, t8_adapt_cmesh_search_element_callback,
                      t8_adapt_cmesh_search_query_callback, &search_queries);
    search_time += sc_MPI_Wtime ();

    /* Adapt the forest according to the markers */
    t8_forest_t         forest_adapt;
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, &markers);
    t8_forest_set_adapt (forest_adapt, forest,
                         t8_forest_adapt_marker_array_callback, 0);

    t8_forest_set_partition (forest_adapt, NULL, 0);
    t8_forest_set_profiling (forest_adapt, 1);
    if (balance) {
      t8_forest_set_balance_ext (forest_adapt, forest, 0, 1);
    }
    non_search_time = -sc_MPI_Wtime ();
    t8_forest_commit (forest_adapt);
    non_search_time += sc_MPI_Wtime ();
    forest = forest_adapt;
    t8_forest_print_profile (forest_adapt);

    search_time_total += search_time;
    non_search_time_total += non_search_time;

    sc_stats_accumulate (&times[0], non_search_time);
    sc_stats_accumulate (&times[1], search_time);
    sc_stats_compute (sc_MPI_COMM_WORLD, 2, times);
    sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 2, times, 1, 1);
  }

#if 0
  if (balance) {
    t8_forest_t         forest_balance;
    t8_forest_init (&forest_balance);
    t8_forest_set_profiling (forest_balance, 1);
    t8_forest_set_balance_ext (forest_balance, forest, 0, 1);
    t8_forest_commit (forest_balance);
    forest = forest_balance;
    t8_forest_print_profile (forest_balance);
  }
#endif

  t8_adapt_cmesh_element_count (forest, element_of_class);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Reduce ((void *) element_of_class, (void *) recv_buff,
                          (int) T8_ECLASS_COUNT, T8_MPI_GLOIDX, sc_MPI_SUM, 0,
                          sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (mpiret);

  if (rank == 0) {
    for (int i = 0; i < T8_ECLASS_COUNT; ++i) {
      t8_global_essentialf ("%s %li\n", t8_eclass_to_string[i],
                            ((t8_gloidx_t *) recv_buff)[i]);
    }
  }

  t8_global_productionf ("\n\tSummarize timings.\n\n");
  sc_stats_accumulate (&total_times[0], non_search_time_total);
  sc_stats_accumulate (&total_times[1], search_time_total);
  sc_stats_compute (sc_MPI_COMM_WORLD, 2, total_times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 2, total_times, 1,
                  1);

  sc_array_reset (&markers);
  sc_array_reset (&search_queries);

  return forest;
}

static void
t8_adapt_cmesh_write_vtk (t8_forest_t forest,
                          t8_forest_t forest_to_adapt_from,
                          const char *vtu_prefix_path, sc_MPI_Comm comm)
{
  /* Write forest to adapt from to vtk.
   * Since this forest is not partitioned (using MPI_COMM_SELF),
   * only one rank should write the files. */
  int                 mpirank, mpiret;
  char                forest_output[BUFSIZ];
  int                 retval;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /*
   * Forest to adapt from
   */
  if (mpirank == 0) {
    retval =
      snprintf (forest_output, BUFSIZ - 1, "%sforest_to_adapt_from",
                vtu_prefix_path);
    if (retval >= BUFSIZ - 1) {
      t8_errorf ("Cannot write vtk output. File path too long.\n");
    }
    else {
#if T8_WITH_VTK
      /* Use VTK library for output if possible */
      t8_forest_write_vtk_via_API (forest_to_adapt_from, forest_output, 1, 1,
                                   1, 1, 0, 0, NULL);
#else
      /* Use standart ascii output if not linked against vtk. */
      t8_forest_write_vtk (forest_to_adapt_from, forest_output);
#endif
    }
  }

  /*
   * Forest
   */

  retval =
    snprintf (forest_output, BUFSIZ - 1, "%sforest_adapt", vtu_prefix_path);
  if (retval >= BUFSIZ - 1) {
    t8_errorf ("Cannot write vtk output. File path too long.\n");
  }
  else {
#if T8_WITH_VTK
    /* Use VTK library for output if possible */
    t8_forest_write_vtk_via_API (forest, forest_output, 1, 1, 1, 1, 0, 0,
                                 NULL);
#else
    /* Use standart ascii output if not linked against vtk. */
    t8_forest_write_vtk (forest, forest_output);
#endif
  }

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
  int                 reflevel;
  int                 template_level;
  int                 dim;
  double              scale[3];
  double              displacement[3];
  const char         *mshfile = NULL;
  const char         *mshfile_forest = NULL;
  const char         *vtu_prefix_path = NULL;
  const sc_MPI_Comm   comm = sc_MPI_COMM_WORLD;
  int                 no_vtk = 0;
  int                 balance = 0;

  /* long help message */
  snprintf (help, BUFSIZ,
            "This program adapts a forest according to a provided forest mesh.\n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (comm, 1, 1, NULL, SC_LP_STATISTICS);
#ifdef T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_STATISTICS);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");

  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The minimum refinement level of the forest that will be refined.");

  sc_options_add_int (opt, 'L', "template_level", &template_level, 0,
                      "The uniform refinement level of the forest that provides the refinement template.");

  sc_options_add_int (opt, 'r', "rlevel", &reflevel, 0,
                      "The maximum refinement level of the forest that will be refined.");

  sc_options_add_switch (opt, 'b', "balance", &balance,
                         "Balance the forest.");

  sc_options_add_string (opt, 'f', "mshfile-template", &mshfile, NULL,
                         "If specified, the forest to adapt from is constructed from a .msh file with "
                         "the given prefix.\n\t\t\t\t     The files must end in .msh "
                         "and be in ASCII format version 2. -d must be specified.");

  sc_options_add_int (opt, 'd', "dim", &dim, -1,
                      "In combination with -f: The dimension of the coarse mesh to read. 1 <= d <= 3.");

  sc_options_add_string (opt, 'm', "mshfile-forest", &mshfile_forest, NULL,
                         "If specified, the forest that is adapted is constructed from a .msh file with "
                         "the given prefix.\n\t\t\t\t     The files must end in .msh "
                         "and be in ASCII format version 2. The dimension is expected to be 3.");

  sc_options_add_double (opt, '\0', "S0", &scale[0], 1,
                         "Scaling in x axis of the cube forest mesh.");
  sc_options_add_double (opt, '\0', "S1", &scale[1], 1,
                         "Scaling in y axis of the cube forest mesh.");
  sc_options_add_double (opt, '\0', "S2", &scale[2], 1,
                         "Scaling in z axis of the cube forest mesh.");
  sc_options_add_double (opt, '\0', "D0", &displacement[0], 0,
                         "Displacement in x axis of the cube forest mesh.");
  sc_options_add_double (opt, '\0', "D1", &displacement[1], 0,
                         "Displacement in y axis of the cube forest mesh.");
  sc_options_add_double (opt, '\0', "D2", &displacement[2], 0,
                         "Displacement in z axis of the cube forest mesh.");

  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Suppress vtk output. "
                         "Overwrites any -O setting.");
  sc_options_add_string (opt, 'O', "output-prefix", &vtu_prefix_path, NULL,
                         "Prefix of vtu output files. Example: \"/home/vtu/prefix_\" will result in the file name \"/home/vtu/prefix_forest_adapt\"\n"
                         "Any folders must already exist.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && level >= 0 && template_level >= 0 && reflevel >= 0 && ((dim >= 1 && dim <= 3) || mshfile == NULL)     /* If dim is provided, then mshfile must be provided as well. */
           &&scale[0] != 0 && scale[1] != 0 && scale[2] != 0) {
    t8_forest_t         forest_to_adapt_from =
      t8_adapt_forest_init_adapt_geometry (comm, mshfile, dim,
                                           template_level);
    t8_forest_t         forest =
      t8_adapt_cmesh_init_forest (comm, level, scale, displacement,
                                  mshfile_forest);

    forest =
      t8_adapt_cmesh_adapt_forest (forest, forest_to_adapt_from,
                                   reflevel - level, balance);

    if (!no_vtk) {
      t8_adapt_cmesh_write_vtk (forest, forest_to_adapt_from, vtu_prefix_path,
                                comm);
    }

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
