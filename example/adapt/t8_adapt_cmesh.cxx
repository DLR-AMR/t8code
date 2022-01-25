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
#include <t8_schemes/t8_default_cxx.hxx>

/* Build a cmesh according to which a later forest shall be refined. */
t8_cmesh_t
t8_adapt_cmesh_init_adapt_geometry (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh = t8_cmesh_new_from_class (T8_ECLASS_TET, comm);
  T8_ASSERT (!t8_cmesh_is_partitioned (cmesh));
  return cmesh;
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
  t8_cmesh_t          cmesh_to_adapt_from;
} t8_adapt_cmesh_user_data_t;

static int
t8_adapt_cmesh_search_callback (t8_forest_t forest,
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
  t8_locidx_t         current_tree_id, num_trees;
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
            "This program adapts a forest according to a provided coarse mesh.\n");
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
    t8_cmesh_t          cmesh_to_adapt_from =
      t8_adapt_cmesh_init_adapt_geometry (comm);
    t8_forest_t         forest = t8_adapt_cmesh_init_forest (comm, level);

    /* Identifiziere Element in forest, die einen Mittelpunkt eines
       Elements in  cmesh_to_adapt_from enthalten. */
    /*[D] Wäre nicht besser: Identifiziere Elemente in forest, 
       deren Mittelpunkt in einem Element des cmesh_to_adapt_from liegen? Im Grobgitter 
       beschreibt, der Mittelpunkt das Element sehr schlecht. */
    /* Adaptiere forest, so dass alle identifizierten Elemente verfeienrt werden. */
    t8_cmesh_destroy (&cmesh_to_adapt_from);
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
