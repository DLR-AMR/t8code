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

#include <t8.h>
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_forest/t8_forest_iterate.h>
#include "t8_latlon_refine.h"
#include "t8_latlon_data.h"
#include "t8_messy_coupler.h"


t8_messy_data* t8_messy_initialize(
  const char* description,
  const char* axis,
  int x_start, 
  int y_start, 
  int x_length, 
  int y_length,
  int z_length, 
  int dimensions) {

  t8_global_productionf("Initializing MESSy coupler\n");

  // create forest for smallest mesh which completely contains given messy mesh
  t8_forest_t forest = t8_latlon_refine(x_length, y_length, T8_LATLON_COARSEN, 0);
  t8_latlon_adapt_data_t *adapt_data =
    (t8_latlon_adapt_data_t *) t8_forest_get_user_data (forest);
  
  // determine axes
  int x, y, z;

  const char *c = strchr(axis, 'X');
  x = c - axis;

  const char *d = strchr(axis, 'Y');
  y = d - axis;

  const char *e = strchr(axis, 'Z');
  z = e - axis;

  /* create data chunk */
  t8_latlon_data_chunk_t *chunk = t8_latlon_new_chunk(
    description,
    x_start, y_start,
    x_length, y_length, z_length,
    dimensions, x, y, z, adapt_data->max_level,
    T8_LATLON_DATA_MESSY);

  t8_messy_data* messy_data = T8_ALLOC(t8_messy_data, 1);
  messy_data->chunk = chunk;
  messy_data->forest = forest;

  t8_global_productionf("MESSy coupler initialized\n");

  return messy_data;
}


void t8_messy_set_dimension(t8_messy_data *messy_data, double ****data, int dimension) {
  t8_latlon_data_chunk_t *chunk = messy_data->chunk;

  /* TODO: add safe guards */

  int axis = chunk->axis;
  int x, y;
  double value;
  for(x=0; x < chunk->x_length; ++x) {
    for(y=0; y < chunk->y_length; ++y) {
      value = t8_latlon_get_dimension_value(axis, data, x, y, 0, 0);
      t8_latlon_set_dimension_value(axis, chunk->in, x, y, 0, dimension, value);
    }
  }
}

void t8_messy_apply_sfc(t8_messy_data *messy_data) {
  t8_latlon_data_apply_morton_order(messy_data->chunk);
}

void t8_messy_coarsen(t8_messy_data *messy_data, t8_forest_adapt_t coarsen_callback, t8_forest_replace_t interpolate_callback) {
  char vtu_prefix[BUFSIZ];

  t8_global_productionf("MESSy coarsen grid \n");

  t8_latlon_data_chunk_t *chunk = messy_data->chunk;

  t8_forest_t forest = messy_data->forest;
  t8_forest_t forest_adapt = messy_data->forest_adapt;

  t8_forest_ref(forest);
  t8_forest_init(&forest_adapt);

  #ifdef T8_ENABLE_DEBUG
    /* In debugging mode write the forest */
    snprintf (vtu_prefix, BUFSIZ, "t8_messy_grid");
    t8_messy_write_forest(forest, vtu_prefix, chunk);
  #endif

  t8_forest_set_user_data(forest_adapt, chunk);
  t8_forest_set_adapt(forest_adapt, forest, coarsen_callback, 0);

  t8_forest_set_partition (forest_adapt, NULL, 0);

  t8_forest_commit(forest_adapt);
  t8_global_productionf("MESSy coarsen done\n");

  t8_global_productionf("MESSy interpolate grid \n");

  int num_element = t8_forest_get_num_element(forest_adapt);

  chunk->data_ids_adapt = T8_ALLOC(t8_linearidx_t, num_element);
  chunk->data_adapt = T8_ALLOC(double, num_element * chunk->z_length * chunk->dimension);

  t8_forest_iterate_replace(forest_adapt, forest, interpolate_callback);

  T8_FREE(chunk->data_ids);
  T8_FREE(chunk->data);

  chunk->data_ids = chunk->data_ids_adapt;
  chunk->data = chunk->data_adapt;

  chunk->data_ids_adapt = NULL;
  chunk->data_adapt = NULL;

  #ifdef T8_ENABLE_DEBUG
    /* In debugging mode write the forest */
    snprintf (vtu_prefix, BUFSIZ, "t8_messy_grid_interpolated");
    t8_messy_write_forest(forest_adapt, vtu_prefix, chunk);
  #endif

  t8_global_productionf("MESSy grid interpolated \n");

}

void t8_messy_write_forest(t8_forest_t forest, const char* prefix, t8_latlon_data_chunk_t *chunk) {

  int num_elements = t8_forest_get_num_element(forest);
  int num_data = chunk->dimension * chunk->z_length;
  t8_vtk_data_field_t vtk_data[num_data];
  double *dim_data_array[num_data];

  int z, d, e, offset;
  for(z = 0; z < chunk->z_length; ++z) {
    for(d = 0; d < chunk->dimension; ++d) {
      offset = z * d + d;
      dim_data_array[offset] = T8_ALLOC_ZERO (double, num_elements);
      for(e = 0; e < num_elements; ++e) {
        dim_data_array[offset][e] = chunk->data[e * num_data + z * chunk->dimension + d];
      }
      snprintf (vtk_data[offset].description, BUFSIZ, "Z %d Dim %d", z, d);
      vtk_data[offset].type = T8_VTK_SCALAR;
      vtk_data[offset].data = (double*) dim_data_array[offset];
    }
  }

  t8_forest_vtk_write_file (forest, prefix, 1, 1, 1, 1, 0, num_data, vtk_data);

  for(offset = 0; offset < num_data; ++offset) {
    T8_FREE(dim_data_array[offset]);
  }
}