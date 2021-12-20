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

  t8_shmem_init (comm);

  t8_shmem_array_init (&context->mptrac_control, sizeof (ctl_t), 1, comm);
  t8_shmem_array_init (&context->mptrac_meteo, sizeof (met_t), 2, comm);

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

  t8_shmem_array_destroy (&context->mptrac_control);
  t8_shmem_array_destroy (&context->mptrac_meteo);

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
  int                 write_permission_to_ctl =
    t8_shmem_array_start_writing (mptrac_context->mptrac_control);
  if (t8_shmem_array_start_writing (mptrac_context->mptrac_meteo)) {
    /* This process has write permission and continues. */

    /* Double check that one process has write permission to all
     * fields. */
    SC_CHECK_ABORT (write_permission_to_ctl,
                    "Shared memory error. Process does not have write access to "
                    "both control and meteo entries.\n");

    t8_debugf ("I have write permission to meteo and control.\n");

    ctl_t              *control = (ctl_t *)
      t8_shmem_array_index_for_writing (mptrac_context->mptrac_control, 0);
    met_t              *meteo1 = (met_t *)
      t8_shmem_array_index_for_writing (mptrac_context->mptrac_meteo, 0);
    met_t              *meteo2 = (met_t *)
      t8_shmem_array_index_for_writing (mptrac_context->mptrac_meteo, 1);
    if (read_ctl_parameters) {
      /* Split command line argument string to be passed to mptrac routines. */
      t8_mptrac_split_input_string (mptrac_context->mptrac_input, &output,
                                    &num_arguments);

      T8_ASSERT (num_arguments > 0);
      read_ctl ("-", num_arguments, output, control);
      /* We need to set the start time by hand. */
      control->t_start = seconds;
      /* Clean up split string. */
      for (int i = 0; i < num_arguments; ++i) {
        T8_FREE (output[i]);
      }
      T8_FREE (output);
    }

    /* Since we are using MPI shared memory, only one process per
     * shared memory region reads the file. */
    get_met (control, seconds, &meteo1, &meteo2);
  }

  /* End writing to shared memory */
  t8_shmem_array_end_writing (mptrac_context->mptrac_control);
  t8_shmem_array_end_writing (mptrac_context->mptrac_meteo);
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
  T8_ASSERT (context != NULL);
  const met_t        *meteo1 =
    (const met_t *) t8_shmem_array_index (context->mptrac_meteo, 0);
  /* Interpolate lon coordinate */
  const int           max_lon_idx = meteo1->nx;
  T8_ASSERT (max_lon_idx >= 1);
  t8_mptrac_interpol_helper (point[0], meteo1->lon[0],
                             meteo1->lon[max_lon_idx - 1], lon);
  /* Interpolate lat coordinate */
  const int           max_lat_idx = meteo1->ny;
  T8_ASSERT (max_lat_idx >= 1);
  /* Note that we switch max and min for latitude.
   * We do this since meteo1->lat[0] = 90, meteo1->lat[max_lat_idx - 1] = -90 in the standard case.
   * But, we want our interpolation to go in positive y axis, not negative.
   * Thus, we need to switch them around.
   */
  t8_mptrac_interpol_helper (point[1],
                             meteo1->lat[max_lat_idx - 1],
                             meteo1->lat[0], lat);

  t8_debugf ("Interpolate (%f,%f) lon range (%f,%f) lat range (%f,%f)\n",
             point[0], point[1], meteo1->lon[0], meteo1->lon[max_lon_idx - 1],
             meteo1->lat[0], meteo1->lat[max_lat_idx - 1]);
  /* Interpolate pressure coordinate */
  const int           max_p_idx = meteo1->np;
  T8_ASSERT (max_p_idx >= 1);
  /* Compute upper and lower bound in km */
  const double        pressure_min_in_km = Z (meteo1->p[0]);
  const double        pressure_max_in_km = Z (meteo1->p[max_p_idx - 1]);
  double              pressure_in_km;
  /* Interpolate pressure in km */
  t8_mptrac_interpol_helper (point[2], pressure_min_in_km,
                             pressure_max_in_km, &pressure_in_km);
  t8_debugf ("Interpolate Z %f pressure range (%f,%f) in km (%f,%f)\n",
             point[2], meteo1->p[0], meteo1->p[max_p_idx - 1],
             pressure_min_in_km, pressure_max_in_km);
  /* Convert back to hpa */
  *pressure = P (pressure_in_km);

  t8_debugf ("Interpolated to (%f, %f, %f)\n", *lon, *lat, *pressure);
}

/* Return coordinate in [0,1] that will map to pressure (in hPa) under t8_mptrac_coords_to_lonlatpressure */
void
t8_mptrac_pressure_to_coord (const t8_mptrac_context_t * context,
                             const double pressure, double *z)
{
  T8_ASSERT (context != NULL);
  const met_t        *meteo1 =
    (const met_t *) t8_shmem_array_index (context->mptrac_meteo, 0);

  /* Convert pressure to km */
  const int           max_p_idx = meteo1->np;
  T8_ASSERT (max_p_idx >= 1);
  const double        pressure_in_km = Z (pressure);
  const double        pressure_min_in_km = Z (meteo1->p[0]);
  const double        pressure_max_in_km = Z (meteo1->p[max_p_idx - 1]);

  /* Map to [0,1] using x in [a,b] maps to (x-a)/(b-a) */
  *z =
    (pressure_in_km - pressure_min_in_km) / (pressure_max_in_km -
                                             pressure_min_in_km);
}
