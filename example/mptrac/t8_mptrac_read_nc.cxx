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
#include <sc_options.h>
#include <t8_messy/t8_latlon_data.h>
#include <t8_messy/t8_latlon_refine.h>
#include <t8_messy/t8_messy_coupler.h>
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>
#include "t8_mptrac_interpolate.h"
#include "thirdparty/mptrac/libtrac.h"

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
                             int itree,
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

  if (value > adapt_ctx->threshold_refine) {
    /* We are larger than the refine threshold. This element gets refined. */
    return 1;
  }
  if (num_elements > 1 && value < adapt_ctx->threshold_coarsen) {
    /* If we are a family and lower than the coarsen threshhold, we coarsen (return -1).
     * Otherwise, do nothing (return 0). */
    return -1;
  }
  /* Keep the element. */
  return 0;
}

/* Refine mesh when u,v,w absolute value is large.
 * coarsen if small. */
int
t8_mptrac_adapt_callback_3d (t8_forest_t forest,
                             t8_forest_t forest_from,
                             int itree,
                             int lelement_id,
                             t8_eclass_scheme_c * ts,
                             int num_elements, t8_element_t * elements[])
{
  const t8_mptrac_adapt_context_t *adapt_ctx =
    (const t8_mptrac_adapt_context_t *) t8_forest_get_user_data (forest);
  T8_ASSERT (adapt_ctx != NULL);
  T8_ASSERT (adapt_ctx->context != NULL);

  /* Get the absolute value of the current element at the requested z layer */
  const int           entries_per_element =
    adapt_ctx->context->data_per_element;
  const int           element_level = ts->t8_element_level (elements[0]);
  const t8_locidx_t   element_offset =
    t8_forest_get_tree_element_offset (forest_from, itree);
  const t8_locidx_t   local_element_id = element_offset + lelement_id;  /* The index of the element in the forest. */
  /* Get a pointer to the u,v,w values */
  const double       *values =
    &adapt_ctx->context->data_array[local_element_id * entries_per_element];

  const double        norm = t8_vec_norm (values);

  /* How many levels we refine/coarsen beyond the initial uniform level. */
  T8_ASSERT (adapt_ctx->context->level >= 0);
  const int           number_of_additional_levels = 1;
  const int           maxlevel =
    adapt_ctx->context->level + number_of_additional_levels;
  const int           minlevel =
    adapt_ctx->context->level - number_of_additional_levels;

  if (norm > adapt_ctx->threshold_refine && element_level < maxlevel) {
    /* We are larger than the refine threshold. This element gets refined. */
    return 1;
  }
  if (num_elements > 1 && norm < adapt_ctx->threshold_coarsen
      && element_level > minlevel) {
    /* If we are a family and lower than the coarsen threshhold, we coarsen (return -1).
     * Otherwise, do nothing (return 0). */
    return -1;
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
t8_mptrac_onlycopy_interpolate_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t itree, t8_eclass_scheme_c * ts, int num_outgoing,  /* previously number of cells, only interesting when 4 */
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

/* Interpolate data from mptrac to a 3D forest.
 * Currenty, we copy only the zonal wind */
void
t8_mptrac_build_latlon_data_for_uvw_3D (t8_mptrac_context_t * context,
                                        double time)
{
  T8_ASSERT (context != NULL);
  T8_ASSERT (context->forest != NULL);
  const t8_locidx_t   num_elements =
    t8_forest_get_local_num_elements (context->forest);
  const t8_locidx_t   num_trees =
    t8_forest_get_num_local_trees (context->forest);
  context->data_per_element = 3;        /* We store u,v,w per element */

  T8_ASSERT (!context->chunk_mode);
  T8_ASSERT (context->data == NULL);

  /* Assure that we have num_elements many entries in our data array */
  /* TODO: When we want ghosts, we need to add num ghosts here as well. */
  context->data_array =
    (double *) T8_REALLOC (context->data_array, double,
                           context->data_per_element * num_elements);

  t8_locidx_t         data_index = 0;
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    for (t8_locidx_t ielement = 0; ielement < num_elements; ielement++) {
      double              midpoint[3];
      const t8_element_t *element =
        t8_forest_get_element_in_tree (context->forest, itree, ielement);
      t8_forest_element_centroid (context->forest, itree, element, midpoint);
      double              lat, lon, pressure;
      t8_mptrac_coords_to_lonlatpressure (context, midpoint, &lon, &lat,
                                          &pressure);
      double              interpol_value;
      int                 ci[3];
      double              cw[3];
      /* Compute interpolation of u */
      intpol_met_time_3d (context->mptrac_meteo1, context->mptrac_meteo1->u,
                          context->mptrac_meteo2, context->mptrac_meteo2->u,
                          time, pressure, lon, lat, &interpol_value, ci, cw,
                          1);
      context->data_array[data_index++] = interpol_value;
      /* Compute interpolation of v */
      intpol_met_time_3d (context->mptrac_meteo1, context->mptrac_meteo1->v, context->mptrac_meteo2, context->mptrac_meteo2->v, time, pressure, lon, lat, &interpol_value, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
      context->data_array[data_index++] = interpol_value;
      /* Compute interpolation of w */
      intpol_met_time_3d (context->mptrac_meteo1, context->mptrac_meteo1->w, context->mptrac_meteo2, context->mptrac_meteo2->w, time, pressure, lon, lat, &interpol_value, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
      context->data_array[data_index++] = interpol_value;

    }
  }
}

void
t8_mptrac_context_write_vtk (const t8_mptrac_context_t * context,
                             const char *vtk_filename)
{
  T8_ASSERT (context != NULL);

  if (context->chunk_mode) {
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
  else {
    /* Not in chunk mode, write data by hand. */
    T8_ASSERT (context->data_array != NULL);
    int                 num_data = 1;
    t8_vtk_data_field_t *vtk_data =
      (t8_vtk_data_field_t *) T8_ALLOC (t8_vtk_data_field_t, num_data);
    vtk_data->data = context->data_array;
    if (context->data_per_element == 1) {
      vtk_data->type = T8_VTK_SCALAR;
      strncpy (vtk_data->description, "mptext_test_3d_u", BUFSIZ);
    }
    else {
      T8_ASSERT (context->data_per_element == 3);
      vtk_data->type = T8_VTK_VECTOR;
      strncpy (vtk_data->description, "mptext_test_3d_uvw", BUFSIZ);
    }
    t8_forest_write_vtk_via_API (context->forest, vtk_filename, 1, 1, 1, 1,
                                 num_data, vtk_data);
    T8_FREE (vtk_data);
  }
}

/* Refine the forest stored in context, but keep it alive.
 * If the dimension is 2 than the z_level input specifies the z_level
 * to consider for refinement.
 * If the dimension is 3 than z_level input is ignored, since all levels
 * are present in the mesh. */
t8_forest_t
t8_mptrac_refine_forest (const t8_mptrac_context_t * context, int z_level,
                         const double threshold_coarsen,
                         const double threshold_refine)
{
  t8_mptrac_adapt_context_t adapt_context;

  adapt_context.context = context;
  adapt_context.threshold_coarsen = threshold_coarsen;
  adapt_context.threshold_refine = threshold_refine;
  adapt_context.z_layer = z_level;
  t8_forest_ref (context->forest);
  /* Adapt the forest. Use 2d or 3d callback depending on dimension. */
  t8_forest_t         forest_adapt = NULL;
  if (context->dimension == 2) {
    forest_adapt =
      t8_forest_new_adapt (context->forest, t8_mptrac_adapt_callback_2d, 0, 0,
                           &adapt_context);
  }
  else {
    forest_adapt =
      t8_forest_new_adapt (context->forest, t8_mptrac_adapt_callback_3d, 0, 0,
                           &adapt_context);
  }
  t8_forest_t         forest_partition;
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  t8_forest_set_balance (forest_partition, forest_adapt, 0);
  t8_forest_commit (forest_partition);
  return forest_partition;
}

void
t8_mptrac_compute_example (const char *filename, const char *mptrac_input,
                           const double simulation_hours, const int dimension,
                           const int level_3d, sc_MPI_Comm comm)
{
  t8_mptrac_context_t *context;
  double              hours;
  double              physical_time;
  const int           chunk_mode = dimension == 2 ? 1 : 0;
  int                 start_six_hours = 0;

  T8_ASSERT (dimension == 2 || dimension == 3);

  /* build context */
  context =
    t8_mptrac_context_new (chunk_mode, filename, mptrac_input, dimension,
                           level_3d, comm);
  /* Compute start time */
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time);

  /* Read NC files to context */
  t8_mptrac_read_nc (context, 1, physical_time, comm);

  /* Build the forest */
  if (dimension == 2) {
    t8_mptrac_build_2d_forest (context);
  }
  else {
    t8_mptrac_build_3d_forest (context, level_3d, sc_MPI_COMM_WORLD);
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
      /* TODO: Can we call read_nc in every time step, but it wont do anything if no file update? */
      time_since_last_six_hours -= 6;
      six_hours_passed++;
      t8_mptrac_read_nc (context, 0, physical_time, comm);
    }
    t8_global_productionf ("Interpolating time %f\n", hours);
    char                vtk_filename[BUFSIZ];
    if (dimension == 2) {
      t8_mptrac_build_latlon_data_for_u_original_coords (context,
                                                         physical_time);
    }
    else {
      t8_mptrac_build_latlon_data_for_uvw_3D (context, physical_time);
    }
    snprintf (vtk_filename, BUFSIZ, "MPTRAC_test_%04i", itime);
    t8_mptrac_context_write_vtk (context, vtk_filename);
    t8_global_productionf ("Wrote file %s\n", vtk_filename);
    t8_forest_t         forest_adapt =
      t8_mptrac_refine_forest (context, 50, 15, 15);
    snprintf (vtk_filename, BUFSIZ, "MPTRAC_test_adapt_%04i", itime);
    t8_forest_write_vtk_via_API (forest_adapt, vtk_filename, 1, 1, 1, 1, 0,
                                 NULL);
    t8_global_productionf ("Wrote file %s\n", vtk_filename);
    if (dimension == 2) {
      /* In 2D mode, we need to go back to the uniform forest, because
       * we require the uniform grid for the interpolation. */
      t8_forest_unref (&forest_adapt);
    }
    else {
      /* In 3D mode, we can continue with the adapted forest, since we interpolate
       * according to physical coordinates (and not the logical grid). */
      t8_forest_unref (&context->forest);
      context->forest = forest_adapt;
    }
    itime++;
    time_since_last_six_hours += deltat;
  }

  /* clean-up */
  t8_mptrac_context_destroy (&context, comm);
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
  int                 dimension;
  int                 level_3d;

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
  sc_options_add_int (opt, 'd', "dimension", &dimension, 2,
                      "The dimension of forest to build. Either 2 or 3 (Default = 2).");
  sc_options_add_int (opt, 'l', "level", &level_3d, 3,
                      "If dimension = 3, the initial level of the 3D forest (Default = 3).");

  /* Parse the command line arguments from the input */
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (netcdf_filename != NULL && (dimension == 2 || dimension == 3)
           && simulation_hours >= 0 && level_3d >= 0) {
    /* Read the netcdf file */
    t8_global_productionf ("Reading nc file %s.\n", netcdf_filename);
    t8_mptrac_compute_example (netcdf_filename, mptrac_input,
                               simulation_hours, dimension, level_3d,
                               sc_MPI_COMM_WORLD);
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
