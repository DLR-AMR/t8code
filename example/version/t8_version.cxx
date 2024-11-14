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

#include <sc_options.h>
#include <t8.h>
#include <t8_version.h>

int
main (int argc, char **argv)
{
  sc_options_t *opt;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int helpme;
  int verbose;
  int sreturn;

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
  t8_init (SC_LP_ERROR);

  /* brief help message */
  sreturn = snprintf (usage, BUFSIZ,
                      "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                      "for a brief overview of all options.",
                      basename (argv[0]), basename (argv[0]));

  if (sreturn >= BUFSIZ) {
    /* Usage string was truncated. */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string to '%s'\n", usage);
  }

  /* long help message */
  sreturn = snprintf (help, BUFSIZ, "This program prints the version number of t8code.\n\n%s\n", usage);
  if (sreturn >= BUFSIZ) {
    /* help message was truncated. */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_switch (opt, 'v', "verbose", &verbose,
                         "Print more information. In particular major, minor and patch version.");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (helpme) {
    /* display help message and usage */
    t8_global_errorf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed < 0) {
    /* wrong usage */
    t8_global_errorf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (!verbose) {
    t8_global_errorf ("%s\n", t8_get_version_number ());
  }
  else {
    /* Verbose print */
    t8_global_errorf ("t8code package name:\t\t%s\n", t8_get_package_string ());
    t8_global_errorf ("t8code version:\t\t\t%s\n", t8_get_version_number ());
    t8_global_errorf ("t8code major version number:\t%i\n", t8_get_version_major ());
    t8_global_errorf ("t8code minor version number:\t%i\n", t8_get_version_minor ());
    t8_global_errorf ("t8code patch version number:\t%i\n", t8_get_version_patch ());
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
