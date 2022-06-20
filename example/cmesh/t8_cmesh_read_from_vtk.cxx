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

#include <t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh.h>
#include <sc_options.h>
#include <t8.h>
#include <t8_vtk.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest.h>

void
t8_cmesh_construct (const char *prefix, sc_MPI_Comm comm, int num_cell_values)
{
  //t8_cmesh_t          cmesh =
  //  t8_cmesh_read_from_vtk_unstructured (prefix, 1, 0, comm);
  t8_cmesh_t          cmesh =
    t8_cmesh_read_from_vtk_poly (prefix, 1, 0, comm);
  t8_forest_t         forest;
  int                 num_trees;
  t8_vtk_data_field_t *vtk_data;
  double            **cell_values;
  double             *tree_data;
  //t8_cmesh_vtk_write_file (cmesh, "test_cmesh", 1.0);
  /* Initialize the forest */
  t8_forest_init (&forest);
  /* Initialize the cmesh of the forest */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  /* Set the scheme of the forest. In this case, the default schemes are used */
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  num_trees = t8_cmesh_get_num_trees (cmesh);
  t8_forest_commit (forest);

  if (num_cell_values > 0) {
    vtk_data = T8_ALLOC (t8_vtk_data_field_t, num_cell_values);
    cell_values = T8_ALLOC (double *, num_cell_values);
    for (int i = 0; i < num_cell_values; i++) {
      cell_values[i] = T8_ALLOC (double, num_trees);
      vtk_data[i].data = cell_values[i];
      /*TODO: Arbitrary type of data */
      vtk_data[i].type = T8_VTK_SCALAR;
      snprintf (vtk_data[i].description, BUFSIZ, "cell_data_%i", i);
    }

    for (int i = 0; i < num_trees; i++) {
      for (int j = 1; j <= num_cell_values; j++) {
        tree_data =
          (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), j,
                                             i);
        cell_values[j - 1][i] = tree_data[0];
      }
    }
  }
  else {
    vtk_data = NULL;
  }

  t8_forest_write_vtk_ext (forest, "forest", 1, 1, 1, 1, 0, 0, 0,
                           num_cell_values, vtk_data);

  if (num_cell_values > 0) {
    for (int i = num_cell_values - 1; i >= 0; i--) {
      T8_FREE (cell_values[i]);
    }
    T8_FREE (cell_values);
    T8_FREE (vtk_data);
  }
  t8_forest_unref (&forest);
  return;
}

int
main (int argc, char **argv)
{
  int                 mpiret, helpme = 0, parsed, num_keys;
  const char         *vtk_file;
  sc_options_t       *opt;
  char                usage[BUFSIZ], help[BUFSIZ];
  int                 sreturn;

  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS> <ARGUMENTS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));
  sreturn = snprintf (help, BUFSIZ,
                      "This program reads a .vtk-file and constructs a mesh representing the given Data."
                      "Arguments can be passed via:\n%s\n\n", usage);
  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_string (opt, 'f', "vtk-file", &vtk_file, "",
                         "The prefix of the .vtk file.");
  sc_options_add_int (opt, 'c', "num_cell_values", &num_keys, 0,
                      "Number of values per cell");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed < 0 || (strcmp (vtk_file, "") == 0)) {
    fprintf (stderr, "%s", help);
    return 1;
  }
  else {
    t8_cmesh_construct (vtk_file, sc_MPI_COMM_WORLD, num_keys);

  }
  sc_options_destroy (opt);
  sc_finalize ();

  return 0;
}
