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

/**
 * Construct a cmesh read from a VTK-file type supported by our vtk-reader.
 * Given the number of cell-values in the file we read the data again from the forest
 * and write it into the file forest.(p)vtu . 
 * 
 * \param[in] prefix  The prefix of the file to read the mesh from
 * \param[in] comm    The communicator used in this example
 * \param[in] values_per_cell   The number of values per cell in the mesh.
 */
void
t8_forest_construct_from_vtk (const char *prefix, sc_MPI_Comm comm,
                              const int values_per_cell)
{
  /* Read a .ply-data file and construct a cmesh representing the mesh. If 
   * there is any cell-data, it will be read too. Triangle-strips and polygons will
   * be broken down to multiple triangles. */
  t8_cmesh_t          cmesh = t8_cmesh_read_from_vtk_poly (prefix, comm);
  t8_forest_t         forest;
  t8_vtk_data_field_t *vtk_data;
  double            **cell_values;
  double             *tree_data;
  /* Initialize the forest */
  t8_forest_init (&forest);
  /* Initialize the cmesh of the forest */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  /* Set the scheme of the forest. In this case, the default schemes are used */
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  const t8_locidx_t   num_trees = t8_cmesh_get_num_trees (cmesh);
  t8_forest_commit (forest);

  /* Read the cell-data if there is any */
  if (values_per_cell > 0) {
    vtk_data = T8_ALLOC (t8_vtk_data_field_t, values_per_cell);
    cell_values = T8_ALLOC (double *, values_per_cell);
    for (int ivalues = 0; ivalues < values_per_cell; ivalues++) {
      cell_values[ivalues] = T8_ALLOC (double, num_trees);
      vtk_data[ivalues].data = cell_values[ivalues];
      /*TODO: Arbitrary type of data */
      vtk_data[ivalues].type = T8_VTK_SCALAR;
      snprintf (vtk_data[ivalues].description, BUFSIZ, "cell_data_%i",
                ivalues);
    }

    for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
      for (int ivalues = 1; ivalues <= values_per_cell; ivalues++) {
        tree_data =
          (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                             ivalues, itree);
        cell_values[ivalues - 1][itree] = tree_data[0];
      }
    }
  }
  else {
    vtk_data = NULL;
  }
  /* Write the forest */
  t8_forest_write_vtk_ext (forest, "forest", 1, 1, 1, 1, 0, 0, 0,
                           values_per_cell, vtk_data);

  /* Free the cell-data */
  if (values_per_cell > 0) {
    for (int ivalues = values_per_cell - 1; ivalues >= 0; ivalues--) {
      T8_FREE (cell_values[ivalues]);
    }
    T8_FREE (cell_values);
    T8_FREE (vtk_data);
  }
  t8_forest_unref (&forest);
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
                      "Number of values per cell the vtk-file stores per cell.");
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed < 0 || (strcmp (vtk_file, "") == 0) || num_keys < 0) {
    fprintf (stderr, "%s", help);
    return 1;
  }
  else {
    t8_forest_construct_from_vtk (vtk_file, sc_MPI_COMM_WORLD, num_keys);

  }
  sc_options_destroy (opt);
  sc_finalize ();

  return 0;
}
