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

#include <t8_interactive_visualization/t8_vis.hxx>

#include <t8.h>
#include <sc_options.h>

void
t8_forest_use_vis_handler (const char *vtk_file, sc_MPI_Comm comm,
                           const int num_keys, const int partition)
{
#if T8_WITH_VTK

  vtkSmartPointer < vtkUnstructuredGrid > grid =
    vtkSmartPointer < vtkUnstructuredGrid >::New ();
  t8_interactive_vis_vtk *test =
    new t8_interactive_vis_vtk (grid, comm, vtk_file);
  test->t8_interactive_vis_source_to_forest ();
  test->t8_interactive_vis_write ();
  t8_debugf ("[D] have read: %lli cells\n",
             test->t8_interactive_vis_get_num_cells ());
#else
  t8_global_errorf
    ("Warning: t8code is not linked againt vtk library. Link t8code with the vtk library to enable this example.\n");
#endif
  test->~t8_interactive_vis_vtk ();
}

int
main (int argc, char **argv)
{
  int                 mpiret, helpme = 0, parsed;
  sc_options_t       *opt;
  char                usage[BUFSIZ], help[BUFSIZ];
  const char         *vtk_file;
  int                 sreturn;
  int                 partition;
  int                 num_keys = 0;

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
  sc_options_add_bool (opt, 'p', "partition", &partition, 0,
                       "If set, partition the cmesh uniformly.\n");
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
    t8_forest_use_vis_handler (vtk_file, sc_MPI_COMM_WORLD, num_keys,
                               partition);

  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;

}
