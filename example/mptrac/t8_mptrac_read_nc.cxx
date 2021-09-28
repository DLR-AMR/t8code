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
#include <t8_messy/t8_latlon_data.h>
#include <t8_messy/t8_latlon_refine.h>
#include <t8_messy/t8_messy_coupler.h>
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include "libtrac.h"

typedef struct
{
  const char         *filename;
  const char         *mptrac_input;
  ctl_t              *mptrac_control;
  met_t              *mptrac_meteo1;
  met_t              *mptrac_meteo2;
  t8_forest_t         forest;
  int                 level;    /*< The maximum refinement level in \a forest. */
  double              missing_value;    /*< The constant we use to denote an invalid data value. */
  int                 chunk_mode;       /*< True if data is stored as latlon_data_chunk_t. Otherwise data is stored in Z-order double array. */
  t8_latlon_data_chunk_t *data; /*< The actual data that we store here, if in chunk mode. */
  double             *data_array;       /*< The actual data that we store here, if not in chunk mode. */
} t8_mptrac_context_t;

typedef struct
{
  int                 z_layer;  /*< The z-layer to consider. */
  double              threshold_coarsen;        /*< Elements with u values lower will get coarsened. */
  double              threshold_refine; /*< Elements with u values larger will get refined. */
  const t8_mptrac_context_t *context;   /*< t8_mptrac context. */
} t8_mptrac_adapt_context_t;

/* Refine mesh when z layer value is large.
 * coarsen if small. */
int
t8_mptrac_adapt_callback_2d (t8_forest_t forest,
                             t8_forest_t forest_from,
                             int which_tree,
                             int lelement_id,
                             t8_eclass_scheme_c * ts,
                             int num_elements, t8_element_t * elements[])
{
  const t8_mptrac_adapt_context_t *adapt_ctx =
    (const t8_mptrac_adapt_context_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_ctx != NULL);
  T8_ASSERT (adapt_ctx->context != NULL);
  T8_ASSERT (adapt_ctx->z_layer >= 0);
  T8_ASSERT (adapt_ctx->z_layer <= adapt_ctx->context->data->z_length);

  /* Get the absolute value of the current element at the requested z layer */
  const int           num_tracers = adapt_ctx->context->data->num_tracers;
  const int           entries_per_element =
    num_tracers * adapt_ctx->context->data->z_length;
  const double        value =
    fabs (adapt_ctx->context->data->data[lelement_id * entries_per_element +
                                         adapt_ctx->z_layer * num_tracers]);

  if (value < adapt_ctx->threshold_coarsen) {
    /* If we are a family and lower than the coarsen threshhold, we coarsen (return -1).
     * Otherwise, do nothing (return 0). */
    return num_elements > 1 ? -1 : 0;
  }
  else if (value > adapt_ctx->threshold_refine) {
    /* We are larger than the refine threshold. This element gets refined. */
    return 1;
  }
  /* Keep the element. */
  return 0;
}

#if 0
/* TODO: This is a modified copy of t8_messy_interpolate_callback.
 * make it publicly available?
 * Interpolate, but only copy value/use value of first family member.
 */
static void
t8_mptrac_onlycopy_interpolate_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, t8_eclass_scheme_c * ts, int num_outgoing,     /* previously number of cells, only interesting when 4 */
                                         t8_locidx_t first_outgoing,    /* index  of first cell in forest_old */
                                         int num_incoming,      /* number of cells to be.., should be 1 */
                                         t8_locidx_t first_incoming)
{                               /* index of new cell in forest_new */

  t8_messy_interpolation_data_t *data =
    (t8_messy_interpolation_data_t *) t8_forest_get_user_data (forest_new);

  int                 index_incoming, index_outgoing;
  int                 element_data_length = data->element_length;

  index_incoming = first_incoming * element_data_length;
  index_outgoing = first_outgoing * element_data_length;

  if (num_outgoing == num_incoming) {
    // no refining was done, just copy data over
    memcpy (data->adapt + index_incoming,
            data->data + index_outgoing,
            element_data_length * sizeof (double));
    return;
  }

  if (num_outgoing < num_incoming) {
    // cell is being split, so copy data to all new cells
    for (int i = 0; i < num_incoming; ++i) {
      // no refining was done, just copy data over
      memcpy (data->adapt + index_incoming + i * element_data_length,
              data->data + index_outgoing,
              element_data_length * sizeof (double));
    }
    return;
  }
}
#endif

/** Split a string that contains command line parameters for MPTRAC into individual
 *  tokens.
 *  Example Input: "INIT_T0 0 INIT_T1 1"
 *          Output: "INIT_T0" "0" "INIT_T1" "1" as output and 4 as num_output.
 * \param [in] input_string The string of command line arguments.
 */
void
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
                   int read_ctl_parameters, double seconds)
{
  int                 num_arguments;
  char              **output;

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
  get_met (mptrac_context->mptrac_control, seconds,
           &mptrac_context->mptrac_meteo1, &mptrac_context->mptrac_meteo2);

  /* Set the missing value. 1e30 seems to be one of the largest values that Paraview can
   * still except as input. */
  mptrac_context->missing_value = -1e30;
  mptrac_context->data = NULL;

}

/* Build a 2D quad forest matching the x and y extend of mptrac meteo data. */
void
t8_mptrac_build_2d_forest (t8_mptrac_context_t * mptrac_context)
{
  mptrac_context->forest =
    t8_latlon_refine (mptrac_context->mptrac_meteo1->nx,
                      mptrac_context->mptrac_meteo1->ny, T8_LATLON_REFINE, 1,
                      &mptrac_context->level);

  T8_ASSERT (mptrac_context->forest != NULL);
  /* Write to vtk */
  t8_forest_write_vtk (mptrac_context->forest, "test_mptrac_forest_2d");
}

/* Build a 3D cube forest and store it at an mptrac_context */
void
t8_mptrac_build_3d_forest (t8_mptrac_context_t * mptrac_context, int level,
                           sc_MPI_Comm comm)
{
  /* TODO: Add periodicity here */
  t8_cmesh_t          cmesh =
    t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  t8_scheme_cxx_t    *scheme = t8_scheme_new_default_cxx ();
  mptrac_context->forest =
    t8_forest_new_uniform (cmesh, scheme, level, 1, comm);
  /* Write to vtk */
  t8_forest_write_vtk (mptrac_context->forest, "test_mptrac_forest_3d");
}

/* Extract 2D slice from specific level into latlon_data chunk.
 * Currenty, we copy only the zonal wind */
void
t8_mptrac_build_latlon_data_for_u_original_coords (t8_mptrac_context_t *
                                                   context, double time)
{
  T8_ASSERT (context != NULL);
  T8_ASSERT (context->forest != NULL);
  const int           num_tracers = 1;
  const int           nx = context->mptrac_meteo1->nx;
  const int           ny = context->mptrac_meteo1->ny;
  const int           np = context->mptrac_meteo1->np;
  int                *shape = T8_ALLOC (int, 4);
  shape[0] = context->mptrac_meteo1->nx;
  shape[1] = context->mptrac_meteo1->ny;
  shape[2] = 1;                 /* Should this be z? */
  shape[3] = 1;                 /* Uncertain about meaning of [3], probably time. */

  T8_ASSERT (context->chunk_mode);
  T8_ASSERT (context->data_array == NULL);

  if (context->data != NULL) {
    t8_latlon_chunk_destroy (&context->data);
  }
  t8_latlon_data_chunk_t *chunk = t8_latlon_new_chunk ("MPTRAC-TESTCHUNK",
                                                       0, 0, nx, ny, np,
                                                       shape, num_tracers,
                                                       'X',
                                                       'Y', 'Z',
                                                       context->level,
                                                       context->missing_value,
                                                       T8_LATLON_DATA_MESSY);

/* TODO: VTK not shown properly, "the data array may be too short" */

  for (int itracer = 0; itracer < num_tracers; ++itracer) {
    snprintf (chunk->tracer_names + itracer * BUFSIZ, BUFSIZ,
              "mptrac_test1_%i", itracer);
  }

  /* Copy the data over to the chunk */
  for (int ix = 0; ix < nx; ++ix) {
    const double        lon = context->mptrac_meteo1->lon[ix];
    for (int iy = 0; iy < ny; ++iy) {
      const double        lat = context->mptrac_meteo1->lat[iy];
      for (int iz = 0; iz < np; ++iz) {
        const double        pressure = context->mptrac_meteo1->p[iz];

        double              value;
        /* Compute interpolation of zonal wind. */
        int                 ci[3];
        double              cw[3];
        intpol_met_time_3d (context->mptrac_meteo1, context->mptrac_meteo1->u,
                            context->mptrac_meteo2, context->mptrac_meteo2->u,
                            time, pressure, lon, lat, &value, ci, cw, 1);

        /* Copy zonal wind */
        chunk->data[(iy * nx + ix) * np + iz] = value;
#if 0                           /* Old manual interpolation routine. */
        (1 - time) * context->mptrac_meteo1->u[ix][iy][iz]
          + time * context->mptrac_meteo2->u[ix][iy][iz];
#endif
        if (ix == 10 && iy == 10 && iz == 10) {
          t8_global_productionf ("%f %f %f %f\n", time,
                                 context->mptrac_meteo1->u[ix][iy][iz],
                                 context->mptrac_meteo2->u[ix][iy][iz],
                                 chunk->data[(iy * nx + ix) * np + iz]);
        }
      }
    }
  }
  /* Convert to morton order */
  t8_latlon_data_apply_morton_order (context->forest, chunk);
  context->data = chunk;
}

/* Interpolate betwenn val1 and val2 at 0 <= interpol <= 1 */
void
t8_mptrac_interpol_helper (const double interpol, const double val1,
                           const double val2, double *output)
{
  *output = (1 - interpol) * val1 + interpol * val2;
}

/* Convert 3D coordinates in [0,1]^3 to lat,lon,pressure coordinates. */
void
t8_mptrac_coords_to_latlonpressure (const t8_mptrac_context_t * context,
                                    const double point[3], double *lat,
                                    double *lon, double *pressure)
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

/* Interpolate data form mptrac to a 3D forest.
 * Currenty, we copy only the zonal wind */
void
t8_mptrac_build_latlon_data_for_u_3D (t8_mptrac_context_t * context,
                                      double time)
{
  T8_ASSERT (context != NULL);
  T8_ASSERT (context->forest != NULL);
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (context->forest);
  const t8_locidx_t   num_trees =
    t8_forest_get_num_local_trees (context->forest);

  T8_ASSERT (!context->chunk_mode);
  T8_ASSERT (context->data == NULL);
  T8_ASSERT (context->data_array != NULL);

  t8_locidx_t         data_index = 0;
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    const double       *tree_vertices =
      t8_forest_get_tree_vertices (context->forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements; ielement++) {
      double              midpoint[3];
      const t8_element_t *element =
        t8_forest_get_element_in_tree (context->forest, itree, ielement);
      t8_forest_element_centroid (context->forest, itree, element,
                                  tree_vertices, midpoint);
      double              lat, lon, pressure;
      t8_mptrac_coords_to_latlonpressure (context, midpoint, &lat, &lon,
                                          &pressure);
      double              interpol_value;
      /* Compute interpolation of zonal wind. */
      int                 ci[3];
      double              cw[3];
      intpol_met_time_3d (context->mptrac_meteo1, context->mptrac_meteo1->u,
                          context->mptrac_meteo2, context->mptrac_meteo2->u,
                          time, pressure, lon, lat, &interpol_value, ci, cw,
                          1);

      context->data_array[data_index++] = interpol_value;
    }
  }
}

void
t8_mptrac_context_write_vtk (const t8_mptrac_context_t * context,
                             const char *vtk_filename)
{
  T8_ASSERT (context != NULL);
  T8_ASSERT (context->data != NULL);
  /* temporarily create messy data struct and use its vtk output. */
  t8_messy_data       data;
  data.adapt_data = NULL;
  data.chunk = context->data;
  data.coarsen = NULL;
  data.interpolation = NULL;
  data.forest = context->forest;
  data.errors = NULL;
  data.errors_global = NULL;
  data.errors_adapt = NULL;
  data.missing_value = context->missing_value;
  data.max_local_error = 0;
  data.max_global_error = 0;
  data.counter = 0;
  data.num_elements = t8_forest_get_local_num_elements (data.forest);
  t8_messy_write_forest (context->forest, vtk_filename, &data);
}

/* Refine the forest stored in context, but keep it alive */
t8_forest_t
t8_mptrac_refine_forest_2d (const t8_mptrac_context_t * context, int z_level,
                            const double threshold_coarsen,
                            const double threshold_refine)
{
  t8_mptrac_adapt_context_t adapt_context;

  adapt_context.context = context;
  adapt_context.threshold_coarsen = threshold_coarsen;
  adapt_context.threshold_refine = threshold_refine;
  adapt_context.z_layer = z_level;
  t8_forest_ref (context->forest);
  t8_forest_t         forest_adapt =
    t8_forest_new_adapt (context->forest, t8_mptrac_adapt_callback_2d, 0, 0,
                         &adapt_context);
  t8_forest_t         forest_partition;
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_commit (forest_partition);
  return forest_partition;
}

t8_mptrac_context_t *
t8_mptrac_context_new (const int chunk_mode, const char *filename,
                       const char *mptrac_input)
{
  t8_mptrac_context_t *context =
    (t8_mptrac_context_t *) T8_ALLOC (t8_mptrac_context_t, 1);

  context->mptrac_meteo1 = T8_ALLOC (met_t, 1);
  context->mptrac_meteo2 = T8_ALLOC (met_t, 1);
  context->mptrac_control = T8_ALLOC (ctl_t, 1);
  context->mptrac_input = mptrac_input;

  /* Set filename */
  context->filename = filename;
  /* Set chunk_mode to 1 or 0 */
  context->chunk_mode = chunk_mode ? 1 : 0;
  /* Set default values */
  context->data = NULL;
  context->data_array = NULL;
  context->forest = NULL;
  context->level = -1;
  context->missing_value = DBL_MAX;

  return context;
}

void
t8_mptrac_context_destroy (t8_mptrac_context_t ** pcontext)
{
  T8_ASSERT (pcontext != NULL);
  T8_ASSERT (*pcontext != NULL);

  t8_mptrac_context_t *context = *pcontext;

  T8_FREE (context->mptrac_meteo1);
  T8_FREE (context->mptrac_meteo2);
  T8_FREE (context->mptrac_control);
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

void
t8_mptrac_compute_example (const char *filename, const char *mptrac_input,
                           const double simulation_hours, const int dimension)
{
  t8_mptrac_context_t *context;
  double              hours;
  double              physical_time;
  const int           chunk_mode = 1;
  int                 start_six_hours = 0;

  T8_ASSERT (dimension == 2 || dimension == 3);

  /* build context */
  context = t8_mptrac_context_new (chunk_mode, filename, mptrac_input);
  /* Compute start time */
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time);
  /* Read NC files to context */
  t8_mptrac_read_nc (context, 1, physical_time);
  /* Build the forest */
  if (dimension == 2) {
    t8_mptrac_build_2d_forest (context);
  }

  const double        deltat = 1;
  T8_ASSERT (deltat < 6);
  double              time_since_last_six_hours = 0;
  int                 itime = 0;
  int                 six_hours_passed = 0;
  for (hours = 0; hours < simulation_hours; hours += deltat) {
    time2jsec (2011, 06, 05, start_six_hours + hours, 00, 00,
               00, &physical_time);
    if (time_since_last_six_hours > 6) {
      time_since_last_six_hours -= 6;
      six_hours_passed++;
      t8_mptrac_read_nc (context, 0, physical_time);
    }
    t8_global_productionf ("Interpolating time %f\n", hours);
    char                vtk_filename[BUFSIZ];
    t8_mptrac_build_latlon_data_for_u_original_coords (context,
                                                       physical_time);
    snprintf (vtk_filename, BUFSIZ, "MPTRAC_test_%04i", itime);
    t8_mptrac_context_write_vtk (context, vtk_filename);
    t8_global_productionf ("Wrote file %s\n", vtk_filename);
    t8_forest_t         forest_adapt =
      t8_mptrac_refine_forest_2d (context, 20, 5, 15);
    snprintf (vtk_filename, BUFSIZ, "MPTRAC_test_adapt_%04i", itime);
    t8_forest_write_vtk_via_API (forest_adapt, vtk_filename, 1, 1, 1, 1, 0,
                                 NULL);
    t8_global_productionf ("Wrote file %s\n", vtk_filename);
    t8_forest_unref (&forest_adapt);
    itime++;
    time_since_last_six_hours += deltat;
  }

  /* clean-up */
  t8_mptrac_context_destroy (&context);
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
  double              simulation_hours;

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
  sc_options_add_double (opt, 'e', "simulation_hours", &simulation_hours,
                         6, "Simulation time in hours.");

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
    t8_mptrac_compute_example (netcdf_filename, mptrac_input,
                               simulation_hours, 2);
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
