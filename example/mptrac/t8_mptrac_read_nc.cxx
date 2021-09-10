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

/* See also: https://github.com/holke/t8code/wiki/Step-0---Hello-World
 *
 * In this example we initialize t8code and print a small welcome message.
 * This is the t8code equivalent of HelloWorld. */

#include <t8.h>
#include <sc_options.h>
#include "libtrac.h"

/** Split a string that contains command line parameters for MPTRAC into individual
 *  tokens.
 *  Example Input: "INIT_T0 0 INIT_T1 1"
 *          Output: "INIT_T0" "0" "INIT_T1" "1" as output and 4 as num_output.
 * \param [in] input_string The string of command line arguments.
 */
void
t8_mptrac_split_input_string (const char *input_string, char ***output,
                              int *num_output)
{
  char               *next_token;
  char                copy_of_input[BUFSIZ];
  int                 num_tokens_read = 0;
  int                 buffer_size = 10;

  /* Basic assertions for wrong usage */
  T8_ASSERT (output != NULL);
  T8_ASSERT (num_output != NULL);

  /* Check if input is NULL */
  if (input_string == NULL) {
    /* No input, hence no output. */
    *num_output = 0;
    *output = NULL;
    return;
  }

  /* Copy the input string since using strtok will modify it. */
  strncpy (copy_of_input, input_string, BUFSIZ - 1);
  /* If input string is too long, then no terminating \0 will be written, so
   * we do it by hand. */
  copy_of_input[BUFSIZ - 1] = '\0';
  if (strlen (copy_of_input) != strlen (input_string)) {
    SC_ABORTF ("Error: Input string was truncated to %s. Aborting.\n",
               copy_of_input);
  }

  /* Split the input string */
  t8_debugf ("Splitting string \"%s\" into tokens:\n", input_string);
  /* Allocate buffer_size many tokens */
  *output = T8_ALLOC (char *, buffer_size);
  next_token = strtok (copy_of_input, " ");
  (*output)[0] = next_token;
  while (next_token != NULL) {
    num_tokens_read++;
    /* Check if we need to allocate more strings */
    if (num_tokens_read >= buffer_size) {
      buffer_size *= 2;
      *output = T8_REALLOC (*output, char *, buffer_size);
    }
    t8_debugf ("%s\n", next_token);
    next_token = strtok (NULL, " ");
    (*output)[num_tokens_read] = next_token;
  }
  *num_output = num_tokens_read;
}

void
t8_mptrac_read_nc (const char *filename, char *mptrac_input)
{
  ctl_t               mptrac_control;
  met_t              *mptrac_meteo;
  mptrac_meteo = T8_ALLOC (met_t, 1);
  int                 num_arguments;
  char              **output;

  /* Split command line argument string to be passed to mptrac routines. */
  t8_mptrac_split_input_string (mptrac_input, &output, &num_arguments);

  read_ctl (filename, num_arguments, output, &mptrac_control);

  T8_FREE (mptrac_meteo);
  T8_FREE (output);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  const char         *netcdf_filename = NULL;
  char               *mptrac_input = NULL;
  int                 parsed, helpme;

  /* help message, prints when called with '-h' option */
  snprintf (help, BUFSIZ,
            "This program uses MPTRAC to read data from a netcdf file.\n");

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  /* Add command line arguments */
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_string (opt, 'n', "netcdffile", &netcdf_filename, NULL,
                         "The netcdf-file that should be read.");
  sc_options_add_string (opt, 'm', "mptrac-input",
                         (const char **) &mptrac_input, NULL,
                         "String of command line arguments passed onto mptrac. Example \"INIT_T0 0 INIT_T1 1\".");

  /* Parse the command line arguments from the input */
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (netcdf_filename != NULL) {
    /* Read the netcdf file */
    t8_global_productionf ("Reading nc file %s.\n", netcdf_filename);
    t8_mptrac_read_nc (netcdf_filename, mptrac_input);
  }
  else {
    /* Error when parsing the arguments */
    /* wrong usage */
    t8_global_essentialf ("\n\tERROR: Wrong usage.\n\n");
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
