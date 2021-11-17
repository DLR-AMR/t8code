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
#include "t8_mptrac_interpolate.h"

t8_mptrac_context_t *
t8_mptrac_context_new (const int chunk_mode, const char *filename,
                       const char *mptrac_input, int dimension,
                       int uniform_level, sc_MPI_Comm comm)
{
  t8_mptrac_context_t *context =
    (t8_mptrac_context_t *) T8_ALLOC (t8_mptrac_context_t, 1);

  T8_ASSERT (dimension == 2 || dimension == 3);

  if (sc_shmem_get_type (comm) == SC_SHMEM_NOT_SET) {
    /* Set the shmem type to the best availble. */
    //t8_shmem_set_type (comm, T8_SHMEM_BEST_TYPE);
    t8_shmem_set_type (comm, SC_SHMEM_WINDOW);
  }
  context->mptrac_meteo1 =
    (met_t *) sc_shmem_malloc (t8_get_package_id (), sizeof (met_t), 1, comm);
  context->mptrac_meteo2 =
    (met_t *) sc_shmem_malloc (t8_get_package_id (), sizeof (met_t), 1, comm);
  context->mptrac_control =
    (ctl_t *) sc_shmem_malloc (t8_get_package_id (), sizeof (ctl_t), 1, comm);
  context->mptrac_input = mptrac_input;

  /* Set filename */
  context->filename = filename;
  /* Set chunk_mode to 1 or 0 */
  context->chunk_mode = chunk_mode ? 1 : 0;
  /* Set dimension */
  context->dimension = dimension;
  /* Set default values */
  context->data = NULL;
  context->data_array = NULL;
  context->forest = NULL;
  context->level = uniform_level;
  context->missing_value = DBL_MAX;

  return context;
}

void
t8_mptrac_context_destroy (t8_mptrac_context_t ** pcontext, sc_MPI_Comm comm)
{
  T8_ASSERT (pcontext != NULL);
  T8_ASSERT (*pcontext != NULL);

  t8_mptrac_context_t *context = *pcontext;

  sc_shmem_free (t8_get_package_id (), context->mptrac_meteo1, comm);
  sc_shmem_free (t8_get_package_id (), context->mptrac_meteo2, comm);
  sc_shmem_free (t8_get_package_id (), context->mptrac_control, comm);
  if (!context->chunk_mode && context->data_array != NULL) {
    T8_FREE (context->data_array);
  }
  if (context->chunk_mode && context->data != NULL) {
    t8_latlon_chunk_destroy (&context->data);
  }
  t8_forest_unref (&context->forest);
  T8_FREE (context);
  *pcontext = NULL;
}

/** Split a string that contains command line parameters for MPTRAC into individual
 *  tokens.
 *  Example Input: "INIT_T0 0 INIT_T1 1"
 *          Output: "INIT_T0" "0" "INIT_T1" "1" as output and 4 as num_output.
 * \param [in] input_string The string of command line arguments.
 */
static void
t8_mptrac_split_input_string (const char *input_string, char ***poutput,
                              int *num_output)
{
  char               *next_token;
  char                copy_of_input[BUFSIZ];
  int                 num_tokens_read = 0;
  int                 buffer_size = 10;

  /* Basic assertions for wrong usage */
  T8_ASSERT (poutput != NULL);
  T8_ASSERT (num_output != NULL);

  /* Check if input is NULL */
  if (input_string == NULL) {
    /* No input, hence no output. */
    *num_output = 0;
    *poutput = NULL;
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
  *poutput = T8_ALLOC (char *, buffer_size);
  char              **output = *poutput;
  next_token = strtok (copy_of_input, " ");
  while (next_token != NULL) {
    /* Check if we need to allocate more strings */
    if (num_tokens_read >= buffer_size) {
      buffer_size *= 2;
      output = T8_REALLOC (output, char *, buffer_size);
    }
    /* Copy current string */
    output[num_tokens_read] = T8_ALLOC (char, BUFSIZ);
    strcpy (output[num_tokens_read], next_token);
    num_tokens_read++;
    t8_debugf ("%s\n", next_token);
    /* Read next string */
    next_token = strtok (NULL, " ");
  }
  *num_output = num_tokens_read;
}

void
t8_mptrac_read_nc (t8_mptrac_context_t * mptrac_context,
                   int read_ctl_parameters, double seconds, sc_MPI_Comm comm)
{
  int                 num_arguments;
  char              **output;

  /* Set the missing value. 1e30 seems to be one of the largest values that Paraview can
   * still except as input. */
  mptrac_context->missing_value = -1e30;
  mptrac_context->data = NULL;

  /* Now read the NC files. Since we are using shared arrays, only the 
   * processes with write permission (thus, one per shared memory domain)
   * reads the file. */
  const int           write_permission_to_meteo1 =
    sc_shmem_write_start (mptrac_context->mptrac_meteo1, comm);
  const int           write_permission_to_meteo2 =
    sc_shmem_write_start (mptrac_context->mptrac_meteo2, comm);
  const int           write_permission_to_ctl =
    sc_shmem_write_start (mptrac_context->mptrac_control, comm);

  if (!write_permission_to_meteo1) {
    /* This process does not have shared memory write permission.
     * We immediately return. */
    return;
  }
  SC_CHECK_ABORT (write_permission_to_meteo2 && write_permission_to_ctl,
                  "Shared memory error. Process does not have write access to "
                  "both meteo entries.\n");

  t8_debugf ("I have %swrite permission to meteo1\n",
             write_permission_to_meteo1 ? "" : "no");
  t8_debugf ("I have %swrite permission to meteo2\n",
             write_permission_to_meteo2 ? "" : "no");
  t8_debugf ("I have %swrite permission to control\n",
             write_permission_to_ctl ? "" : "no");

  if (read_ctl_parameters) {
    /* Split command line argument string to be passed to mptrac routines. */
    t8_mptrac_split_input_string (mptrac_context->mptrac_input, &output,
                                  &num_arguments);

    T8_ASSERT (num_arguments > 0);
    read_ctl ("-", num_arguments, output, mptrac_context->mptrac_control);
    /* We need to set the start time by hand. */
    mptrac_context->mptrac_control->t_start = seconds;
    /* Clean up split string. */
    for (int i = 0; i < num_arguments; ++i) {
      T8_FREE (output[i]);
    }
    T8_FREE (output);
  }

  /* Since we are using MPI shared memory, only one process per
   * shared memory region reads the file. */
  get_met (mptrac_context->mptrac_control, seconds,
           &mptrac_context->mptrac_meteo1, &mptrac_context->mptrac_meteo2);

  /* End writing to shared memory */
  sc_shmem_write_end (mptrac_context->mptrac_meteo1, comm);
  sc_shmem_write_end (mptrac_context->mptrac_meteo2, comm);
  sc_shmem_write_end (mptrac_context->mptrac_control, comm);
}

/* Interpolate between val1 and val2 at 0 <= interpol <= 1 */
static void
t8_mptrac_interpol_helper (const double interpol, const double val1,
                           const double val2, double *output)
{
  *output = (1 - interpol) * val1 + interpol * val2;
}

/* Convert 3D coordinates in [0,1]^3 to lat,lon,pressure coordinates. */
void
t8_mptrac_coords_to_lonlatpressure (const t8_mptrac_context_t * context,
                                    const double point[3], double *lon,
                                    double *lat, double *pressure)
{
  /* Interpolate lon coordinate */
  const int           max_lon_idx = context->mptrac_meteo1->nx;
  T8_ASSERT (max_lon_idx >= 1);
  t8_mptrac_interpol_helper (point[0], context->mptrac_meteo1->lon[0],
                             context->mptrac_meteo1->lon[max_lon_idx - 1],
                             lon);
  /* Interpolate lat coordinate */
  const int           max_lat_idx = context->mptrac_meteo1->ny;
  T8_ASSERT (max_lat_idx >= 1);
  t8_mptrac_interpol_helper (point[1], context->mptrac_meteo1->lat[0],
                             context->mptrac_meteo1->lat[max_lat_idx - 1],
                             lat);
  /* Interpolate pressure coordinate */
  const int           max_p_idx = context->mptrac_meteo1->np;
  T8_ASSERT (max_p_idx >= 1);
  t8_mptrac_interpol_helper (point[2], context->mptrac_meteo1->p[0],
                             context->mptrac_meteo1->p[max_p_idx - 1],
                             pressure);
}
