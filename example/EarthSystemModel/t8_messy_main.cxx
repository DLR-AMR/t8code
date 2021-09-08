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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#ifdef T8_ENABLE_DEBUG
/* In debugging mode, we check whether we use the correct scheme.
 * To do so, we have to include the quad default scheme directly. */
#include <t8_schemes/t8_default/t8_default_quad_cxx.hxx>
#endif
#include <p4est_bits.h>
#include <t8_messy/t8_messy_helper.h>
#include <t8_messy/t8_messy_coupler.h>
#include <t8_messy/t8_latlon_refine.h>
#include <t8_messy/t8_latlon_data.h>

/* TODO: Rename those functions to t8_messy_NAME */

int custom_coarsening(t8_messy_custom_func_t* arguments) {
  t8_debugf("custom coarsening\n");
  return -1;
}

double custom_interpolation(t8_messy_custom_func_t* arguments) {
  t8_debugf("custom interpolating\n");
  return 100.0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 x_length, y_length;
  int                 parsed, helpme;
  int                 mode_int;
  int                 partition;
  

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "Given input dimension x and y, we construct a\n"
            "forest on the unit square that is the coarsest forest such\n"
            "that an x times y grid fits in the lower left corner.\n%s\n",
            usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'x', "x-length", &x_length, 32,
                      "The type of elements to use.\n");
  sc_options_add_int (opt, 'y', "y-length", &y_length, 32,
                      "The type of elements to use.\n");
  sc_options_add_switch (opt, 'p', "partition", &partition,
                         "Repartition the forest after each level of refinement/coarsening.\n");
  sc_options_add_int (opt, 'm', "modus", &mode_int, 0,
                      "The adaptation modus to use\n"
                      "\t\t0 - refine modus: We start with level 0 and refine until the final forest ist constructed.\n"
                      "\t\t1 - coarsen modus: We start with the final uniform level and coarsen elements until the final forest ist constructed.\n");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 < x_length && 0 < y_length && mode_int >= 0
           && mode_int <= 1) {


    /* number of datapoints per grid cell */
    int num_dims = 2, x, y;


    /* allocate data array */
    double ****data = T8_ALLOC(double***, x_length);
    for(x=0; x<x_length; ++x) {
      data[x] = T8_ALLOC(double**, y_length);
      for(y=0; y<y_length; ++y) {
        data[x][y] = T8_ALLOC(double*, 1);
        data[x][y][0] = T8_ALLOC(double, 1);
      }
    }

    int* shape = T8_ALLOC(int, 4);
    shape[0] = x_length;
    shape[1] = y_length;
    shape[2] = 1;

    t8_messy_coarsen_t *coarsen = t8_messy_new_coarsen_config("mean_higher", "vm1", 0, 0.8, NULL);
    t8_messy_interpolate_t *interpolation = t8_messy_new_interpolate_config("mean", NULL);

    /* initialize forest and data chunk */
    t8_messy_data_t* messy = t8_messy_initialize("test", "XYZ", shape, 0, 0, num_dims, -1.0, 0.0, coarsen, interpolation);

    /* set data for every dimension */
    //char name[BUFSIZ];

    t8_messy_gaussian(data, x_length, y_length);
    //sprintf(name, "gaussian");
    // t8_messy_add_dimension(messy, "gaussian", data);
    
    // t8_messy_sine_2d(data, x_length, y_length);
    //sprintf(name, );
    // t8_messy_add_dimension(messy, "sine_2d", data);

    /* bring input data into SFC format */
    t8_messy_apply_sfc(messy);

    t8_messy_coarsen_t *coarsen_config = T8_ALLOC(t8_messy_coarsen_t, 1);
    t8_messy_interpolate_t *interpolate_config = T8_ALLOC(t8_messy_interpolate_t, 1);

    // coarsen_config->method = T8_MESSY_COARSEN_FUNCTION;
    // coarsen_config->func = custom_coarsening;
    coarsen_config->method = T8_MESSY_COARSEN_THRESHOLD_MIN_HIGHER;
    coarsen_config->z_layer = 0;
    coarsen_config->tracer = "gaussian";
    coarsen_config->threshold = 0.8;

    //interpolate_config->method = T8_MESSY_INTERPOLATE_MEAN;
    interpolate_config->method = T8_MESSY_INTERPOLATE_FUNCTION;
    interpolate_config->func = custom_interpolation;

    messy->coarsen = coarsen_config;
    messy->interpolation = interpolate_config;

    /* coarsen data */
    t8_messy_coarsen(messy);

    t8_forest_unref (&(messy->forest));
    /* t8_forest_unref (&(messy->forest_adapt)); */
    t8_latlon_chunk_destroy(&(messy->chunk));


  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}





