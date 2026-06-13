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
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#include <sc_options.h>
#include <string>

/** Removes all elements of a local tree if they belong to the corresponding
 *  global trees which is given by the user_data. */
static int
t8_adapt_remove (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                 [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const t8_locidx_t lelement_id,
                 [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                 [[maybe_unused]] const int num_elements, [[maybe_unused]] t8_element_t *elements[])
{
  const t8_gloidx_t *tree_id = (const t8_gloidx_t *) t8_forest_get_user_data (forest);
  const t8_gloidx_t global_tree_id = t8_forest_global_tree_id (forest_from, which_tree);
  if (global_tree_id == *tree_id) {
    return -2;
  }
  return 0;
}

void
t8_strip_of_quads (t8_gloidx_t num_trees, t8_gloidx_t empty_tree, const char **vtuname)
{

  const double boundary_coords[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 };
  const int use_axis_alined = 1;

  t8_cmesh_t cmesh
    = t8_cmesh_new_hypercube_pad (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, boundary_coords, num_trees, 1, 0, use_axis_alined);

  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, sc_MPI_COMM_WORLD);

  t8_forest_write_vtk (forest, *vtuname);
  t8_debugf ("Output to %s\n", *vtuname);

  t8_forest_ref (forest);
  t8_forest_t forest_adapt;
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_adapt_remove, 0);
  t8_forest_set_partition (forest_adapt, NULL, 0);
  t8_forest_set_user_data (forest_adapt, &empty_tree);
  t8_forest_commit (forest_adapt);

  std::string vtuname_adapt (*vtuname);
  vtuname_adapt = vtuname_adapt + "_adapted";
  t8_forest_write_vtk (forest_adapt, vtuname_adapt.c_str ());
  t8_debugf ("Output to %s\n", vtuname_adapt.c_str ());

  t8_productionf ("The initial uniform forest:\n"
                  "\tfirst_local_tree: %" T8_GLOIDX_FORMAT "\n"
                  "\tlast_local_tree:  %" T8_GLOIDX_FORMAT "\n"
                  "\tlocal_num_trees:  %" T8_LOCIDX_FORMAT "\n"
                  "\tglobal_num_trees: %" T8_GLOIDX_FORMAT "\n",
                  forest->first_local_tree, forest->last_local_tree, t8_forest_get_num_local_trees (forest),
                  t8_forest_get_num_global_trees (forest));

  t8_productionf ("The adapted forest with one empty tree:\n"
                  "\tfirst_local_tree: %" T8_GLOIDX_FORMAT "\n"
                  "\tlast_local_tree:  %" T8_GLOIDX_FORMAT "\n"
                  "\tlocal_num_trees:  %" T8_LOCIDX_FORMAT "\n"
                  "\tglobal_num_trees: %" T8_GLOIDX_FORMAT "\n",
                  forest_adapt->first_local_tree, forest_adapt->last_local_tree,
                  t8_forest_get_num_local_trees (forest_adapt), t8_forest_get_num_global_trees (forest_adapt));

  t8_forest_unref (&forest_adapt);
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* brief help message */
  int sreturnA = snprintf (usage, BUFSIZ,
                           "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                           "for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* long help message */
  int sreturnB = snprintf (help, BUFSIZ,
                           "We create a forest with a strip of quad trees.\n"
                           "One tree of this strip does not contain any elements.\n\n%s\n",
                           usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Parameter for t8_construct_fractal and command line */
  t8_locidx_t num_trees;
  t8_locidx_t empty_tree;
  const char *vtuname[BUFSIZ];
  int helpme;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 't', "number of trees", &num_trees, 3, "The number of trees in row. Default is 3.");
  sc_options_add_int (opt, 'e', "empty tree", &empty_tree, 1,
                      "The global index of the tree to be empty. Default is 1.");
  sc_options_add_string (opt, 'p', "output path", vtuname, "t8_strip_of_quads", "Path of outputfiles.\n");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 < num_trees && 0 <= empty_tree && empty_tree < num_trees) {
    t8_strip_of_quads ((t8_gloidx_t) num_trees, (t8_gloidx_t) empty_tree, vtuname);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
