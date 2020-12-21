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
#include "t8_messy_coupler.h"
#include "t8_latlon_refine.h"
#include "t8_latlon_data.h"

/* generate a random floating point number from min to max */
double randfrom(double min, double max) {
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

/**
 * Function filling data array with random values
 */
void generate_data(double ****data, int x_length, int y_length, double value) {
  int x, y;
  for(y=0; y<y_length; ++y) {
    for(x=0; x<x_length; ++x) {
      data[x][y][0][0] = randfrom(y * x_length, (y + 1) * x_length);
    }
  }
}

void sine_2d(double ****data, int x_length, int y_length) {
  double T_x = x_length / 5.0;
  double T_y = y_length / 5.0;
  int x, y;
  for(y=0; y<y_length; ++y) {
    for(x=0; x<x_length; ++x) {
      data[x][y][0][0] = sin( (2 * M_PI * x) / T_x + (2 * M_PI* y) / T_y );
    }
  }
}

void gaussian(double ****data, int x_length, int y_length) {
  double x0 = x_length * 1.0 / 2.0;
  double y0 = y_length * 1.0 / 2.0;
  double xd, yd, A, ox, oy;
  int x, y;
  
  ox = x0;
  oy = y0;
  A = 100;
  for(y=0; y<y_length; ++y) {
    for(x=0; x<x_length; ++x) {
      xd = x - x0;
      yd = y - y0;
      t8_debugf("%f\n", A * exp(- (((xd*xd) / (2.0 * (ox * ox)) + ((yd*yd)/ (2.0 * (oy * oy) ))))));
      data[x][y][0][0] = A * exp(- (((xd*xd) / (2.0 * (ox * ox)) + ((yd*yd)/ (2.0 * (oy * oy) )))));
    }
  }
}

/**
 * Simple coarsening test.
 * If the AVG of the first dimension of all four cells is even we coarsen.
 */
int
t8_messy_coarsen_callback (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          int num_elements, t8_element_t * elements[])
{


  t8_latlon_data_chunk_t *chunk = (t8_latlon_data_chunk_t*) t8_forest_get_user_data(forest);
  
  /* since we don't want to refine, 
     we can stop if we only have one element */
  if (num_elements == 1) {
    return 0;
  }

  double avg = 0.0;
  for (int i=0; i < num_elements; ++i) {
    int offset = (lelement_id + i) * chunk->dimension;
    avg += chunk->data[offset];
  }
  avg /= num_elements;

  bool isEven = (((int)avg) % 2) == 0;

  t8_debugf ("lelement_id %d avg %.4f, is even? %s \n", lelement_id, avg, (isEven ? "yes" : "no"));

  return isEven ? -1 : 0;
}

static void
t8_messy_replace_callback (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c * ts,
                   int num_outgoing, /* previously number of cells, only interesting when 4 */
                   t8_locidx_t first_outgoing, /* index  of first cell in forest_old */
                   int num_incoming, /* number of cells to be.., should be 1 */
                   t8_locidx_t first_incoming) /* index of new cell in forest_new */
{

  t8_latlon_data_chunk_t *data_chunk = (t8_latlon_data_chunk_t *) t8_forest_get_user_data(forest_new);
  int dimensions, z_length, element_data_length;
  dimensions = data_chunk->dimension;
  z_length = data_chunk->z_length;
  element_data_length = dimensions * z_length;

  int index_incoming = first_incoming * element_data_length;
  int index_outgoing = first_outgoing * element_data_length;

  t8_debugf("num_out %i, num_in %i\n", num_outgoing, num_incoming);
  if(num_outgoing > num_incoming) {
    t8_debugf("interpolating\n");
  /* when the number of previous elements (num_outgoing) is larger than the number of created cell from it (num_incoming)
   * we interpolate,
   */

    int d, z, e, offset;
    double value;
    for(z = 0; z < z_length; ++z) {
      for(d = 0; d < dimensions; ++d) {
        offset = z * d + d;
        value = 0.0;
        for(e = 0; e < num_outgoing; ++e) {
          value += data_chunk->data[index_outgoing + e * element_data_length + offset];
        }
        value /= num_outgoing;
        data_chunk->data_adapt[index_incoming + offset] = value;
      }
    }
  } else {
    t8_debugf("not interpolating\n");
    /* else just copy data over to new array */
    memcpy (data_chunk->data_adapt + index_incoming,
            data_chunk->data       + index_outgoing,
              element_data_length * sizeof (double));
  }

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
    int num_dims = 1, x, y, z;


    /* allocate data array */
    double ****data = T8_ALLOC(double***, x_length);
    for(x=0; x<x_length; ++x) {
      data[x] = T8_ALLOC(double**, y_length);
      for(y=0; y<y_length; ++y) {
        data[x][y] = T8_ALLOC(double*, 1);
        data[x][y][0] = T8_ALLOC(double, 1);
      }
    }

    /* initialize forest and data chunk */
    t8_messy_data* messy = t8_messy_initialize("test", "XYZ", 0, 0, x_length, y_length, 1, num_dims);

    /* set data for every dimension */
    for (int dim=0; dim<num_dims; ++dim) {
      /* generate dummy data */
      // generate_data(data, x_length, y_length, dim * 1.0);
      // sine_2d(data, x_length, y_length);
      gaussian(data, x_length, y_length);
      t8_messy_set_dimension(messy, data, dim);
    }

    /* bring input data into SFC format */
    t8_messy_apply_sfc(messy);

    /* coarsen data */
    t8_messy_coarsen(messy, t8_messy_coarsen_callback, t8_messy_replace_callback);

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





