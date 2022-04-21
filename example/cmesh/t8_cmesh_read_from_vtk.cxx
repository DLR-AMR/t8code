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

#include <t8_cmesh/t8_cmesh_vtk_reader.hxx>
#include <sc_options.h>
#include <t8.h>

void
t8_cmesh_construct (const char *prefix, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh = t8_cmesh_read_from_vtk (prefix, 1, 1, comm);
  t8_cmesh_destroy (&cmesh);
  return;
}

int
main (int argc, char **argv)
{
  int                 mpiret, helpme = 0, parsed;
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
    t8_cmesh_construct (vtk_file, sc_MPI_COMM_WORLD);

  }
  sc_options_destroy (opt);
  sc_finalize ();

  return 0;
}
