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

inline 
double get_mean(int first, int num_elements, int element_length, int dimension, double* data) {
  double mean_value = 0.0;
  int offset, i;
  for (i=0; i < num_elements; ++i) {
    offset = first + i * element_length + dimension;
    mean_value += data[offset];
  }
  mean_value /= num_elements;
  return mean_value;
}

inline 
double get_max(int first, int num_elements, int element_length, int dimension, double* data) {
  double max_value = LONG_MIN;
  int offset, i;
  for (i=0; i < num_elements; ++i) {
    offset = first + i * element_length + dimension;
    max_value = fmax(max_value, data[offset]);
  }
  return max_value;
}

inline 
double get_min(int first, int num_elements, int element_length, int dimension, double* data) {
  double min_value = LONG_MAX;
  int offset, i;
  for (i=0; i < num_elements; ++i) {
    offset = first + i * element_length + dimension;
    min_value = fmin(min_value, data[offset]);
  }
  return min_value;
}

inline
void get_values(int first, int num_elements, int element_length, int dimension, double* values, double* data) {
  int offset, i;
  for (i=0; i < num_elements; ++i) {
    offset = first + i * element_length + dimension;
    values[i] = data[offset];
  }
}

t8_messy_custom_func_t* t8_messy_new_custom_func(int num_elements) {
  t8_messy_custom_func_t* func_data = T8_ALLOC(t8_messy_custom_func_t, 1);
  func_data->num_elements = num_elements;
  func_data->x_coords = T8_ALLOC_ZERO(int, num_elements);
  func_data->y_coords = T8_ALLOC_ZERO(int, num_elements);
  func_data->latitudes = T8_ALLOC_ZERO(double, num_elements);
  func_data->longitudes = T8_ALLOC_ZERO(double, num_elements);
  func_data->values = T8_ALLOC_ZERO(double, num_elements);
  func_data->dimension = T8_ALLOC(char, BUFSIZ);

  return func_data;
}

void t8_messy_destroy_custom_func(t8_messy_custom_func_t* custom) {
  T8_FREE(custom->x_coords);
  T8_FREE(custom->y_coords);
  T8_FREE(custom->latitudes);
  T8_FREE(custom->longitudes);
  T8_FREE(custom->values);
  T8_FREE(custom->dimension);
  T8_FREE(custom);
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
  int ret = 0;

  /* since we don't want to refine, 
     we can stop if we only have one element */
  if (num_elements == 1) {
    return ret;
  }

  t8_messy_data_t *messy_data = (t8_messy_data_t*) t8_forest_get_user_data(forest);
  t8_messy_coarsen_t *coarsen = messy_data->coarsen;
  t8_latlon_data_chunk_t *data_chunk = messy_data->chunk;

  t8_messy_custom_func_t* func_data = NULL;

  /* calculate how many values one element has */
  int element_length = data_chunk->z_length * data_chunk->dimensions;
  /* calculate offset for z_layer */
  int z_offset = coarsen->z_layer * data_chunk->dimensions;
  /* calculate start index for first element */
  int start = lelement_id * data_chunk->z_length * data_chunk->dimensions + z_offset;
  /* get dimension index */
  int dimension = t8_latlon_get_dimension_idx(data_chunk, coarsen->dimension, false);

  

  double value;
  switch(coarsen->method) {
    case T8_MESSY_COARSEN_FUNCTION:
      /* TODO: implement  custom coarsening functions */
      func_data = t8_messy_new_custom_func(num_elements);
      func_data->z_layer = coarsen->z_layer;
      strcpy(func_data->dimension, coarsen->dimension);

      /* set lat, lon, x, y,  */

      get_values(start, num_elements, element_length, dimension, func_data->values, data_chunk->data);
      ret = (coarsen->func)(func_data);
      t8_messy_destroy_custom_func(func_data);
      func_data = NULL;
      break;
    case T8_MESSY_COARSEN_AREA_INSIDE:
      /* TODO: area coarsening */
      break;
    case T8_MESSY_COARSEN_AREA_OUTSIDE:
      /* TODO: area coarsening */
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MIN_LOWER:
      value = get_min(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value < coarsen->threshold ? -1 : 0;
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MIN_HIGHER:
      value = get_min(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value > coarsen->threshold ? -1 : 0;
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MAX_LOWER:
      value = get_max(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value < coarsen->threshold ? -1 : 0;
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MAX_HIGHER:
      value = get_max(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value > coarsen->threshold ? -1 : 0;
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MEAN_LOWER:
      value = get_mean(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value < coarsen->threshold ? -1 : 0;
      break;
    case T8_MESSY_COARSEN_THRESHOLD_MEAN_HIGHER:
      value = get_mean(start, num_elements, element_length, dimension, data_chunk->data);
      ret = value > coarsen->threshold ? -1 : 0;
      break;
  }


  return ret;
}

static void
t8_messy_interpolate_callback (t8_forest_t forest_old,
                   t8_forest_t forest_new,
                   t8_locidx_t which_tree,
                   t8_eclass_scheme_c * ts,
                   int num_outgoing, /* previously number of cells, only interesting when 4 */
                   t8_locidx_t first_outgoing, /* index  of first cell in forest_old */
                   int num_incoming, /* number of cells to be.., should be 1 */
                   t8_locidx_t first_incoming) /* index of new cell in forest_new */
{

  t8_messy_data_t *messy_data = (t8_messy_data_t *) t8_forest_get_user_data(forest_new);
  t8_messy_interpolate_t *interpolation = messy_data->interpolation;
  t8_latlon_data_chunk_t *data_chunk = messy_data->chunk;

  int num_dimensions, z_length, element_data_length;
  int index_incoming, index_outgoing;

  num_dimensions = data_chunk->dimensions;
  z_length = data_chunk->z_length;
  element_data_length = num_dimensions * z_length;

  index_incoming = first_incoming * element_data_length;
  index_outgoing = first_outgoing * element_data_length;

  // t8_debugf("num_out %i, num_in %i\n", num_outgoing, num_incoming);

  if(num_outgoing > num_incoming) {
    // t8_debugf("interpolating\n");
    /* when the number of previous elements (num_outgoing) is larger than the number of created cell from it (num_incoming)
     * we interpolate,
     */
    t8_messy_custom_func_t* func_data = NULL;

    int d, z, z_offset, start;
    double value;
    
    switch(interpolation->method) {
      case T8_MESSY_INTERPOLATE_FUNCTION:
        /* TODO: implement custom interpolation */
        func_data = t8_messy_new_custom_func(num_outgoing);
        
        /* loop over elements and fill lat, lon, x, y in func data */

        for(z = 0; z < z_length; ++z) {
          func_data->z_layer = z;
          /* calculate offset for z_layer */
          z_offset = z * num_dimensions;
          /* calculate start index for first element */
          start = index_outgoing + z_offset;
          for(d = 0; d < num_dimensions; ++d) {
            strcpy(func_data->dimension, data_chunk->dimension_names + d * BUFSIZ);
            get_values(start, num_outgoing, element_data_length, d, func_data->values, data_chunk->data);
            value = (double)(interpolation->func)(func_data);
            data_chunk->data_adapt[index_incoming + z_offset + d] = value;
          }
        }
        t8_messy_destroy_custom_func(func_data);
        func_data = NULL;
        break;

      case T8_MESSY_INTERPOLATE_MIN:
        for(z = 0; z < z_length; ++z) {
          /* calculate offset for z_layer */
          z_offset = z * num_dimensions;
          /* calculate start index for first element */
          start = index_outgoing + z_offset;
          for(d = 0; d < num_dimensions; ++d) {
            value = get_min(start, num_outgoing, element_data_length, d, data_chunk->data);
            data_chunk->data_adapt[index_incoming + z_offset + d] = value;
          }
        }
        break;

      case T8_MESSY_INTERPOLATE_MAX:
        for(z = 0; z < z_length; ++z) {
          /* calculate offset for z_layer */
          z_offset = z * num_dimensions;
          /* calculate start index for first element */
          start = index_outgoing + z_offset;
          for(d = 0; d < num_dimensions; ++d) {
            value = get_max(start, num_outgoing, element_data_length, d, data_chunk->data);
            data_chunk->data_adapt[index_incoming + z_offset + d] = value;
          }
        }
        break;

      case T8_MESSY_INTERPOLATE_MEAN:
        for(z = 0; z < z_length; ++z) {
          /* calculate offset for z_layer */
          z_offset = z * num_dimensions;
          /* calculate start index for first element */
          start = index_outgoing + z_offset;
          for(d = 0; d < num_dimensions; ++d) {
            value = get_mean(start, num_outgoing, element_data_length, d, data_chunk->data);
            data_chunk->data_adapt[index_incoming + z_offset + d] = value;
          }
        }
        break;
    }

  } else {
    t8_debugf("not interpolating\n");
    /* else just copy data over to new array */
    memcpy (data_chunk->data_adapt + index_incoming,
            data_chunk->data       + index_outgoing,
              element_data_length * sizeof (double));
  }

}

t8_messy_coarsen_t* t8_messy_new_coarsen_config(
  const char* method,
  char* dimension,
  int z_layer,
  double threshold,
  int (*func)(t8_messy_custom_func_t *)
 ) {
  t8_messy_coarsen_t* config = T8_ALLOC(t8_messy_coarsen_t, 1);
  
  config->dimension = dimension;
  config->z_layer = z_layer;

  config->threshold = threshold;
  
  if(strcmp(method, "mean_lower") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MEAN_LOWER;
  } else if(strcmp(method, "mean_higher") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MEAN_HIGHER;
  } else if(strcmp(method, "min_lower") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MIN_LOWER;
  } else if(strcmp(method, "min_higher") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MIN_HIGHER;
  } else if(strcmp(method, "max_lower") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MAX_LOWER;
  } else if(strcmp(method, "max_higher") == 0) {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MAX_HIGHER;
  } else if(strcmp(method, "custom") == 0) {
    config->method = T8_MESSY_COARSEN_FUNCTION;
  } else {
    config->method = T8_MESSY_COARSEN_THRESHOLD_MEAN_LOWER;
  }
  
  config->func = func;

  return config;
}


t8_messy_interpolate_t* t8_messy_new_interpolate_config(
  const char* method,
  double (*func)(t8_messy_custom_func_t *)
) {
  t8_messy_interpolate_t *config = T8_ALLOC(t8_messy_interpolate_t, 1);

  if(strcmp(method, "mean") == 0) {
    config->method = T8_MESSY_INTERPOLATE_MEAN;
  } else if(strcmp(method, "min") == 0) {
    config->method = T8_MESSY_INTERPOLATE_MIN;
  } else if(strcmp(method, "max") == 0) {
    config->method = T8_MESSY_INTERPOLATE_MAX;
  } else if(strcmp(method, "custom") == 0) {
    config->method = T8_MESSY_INTERPOLATE_FUNCTION;
  } else {
    config->method = T8_MESSY_INTERPOLATE_MEAN;
  }

  config->func = func;

  return config;
}


t8_messy_data_t* t8_messy_initialize(
  const char* description,
  const char* axis,
  int* shape,
  int x_start, 
  int y_start,
  int dimensions,
  t8_messy_coarsen_t *coarsen,
  t8_messy_interpolate_t *interpolation
  ) {

  #ifdef T8_ENABLE_DEBUG
    t8_global_productionf("Initializing MESSy coupler\n");
  #endif

  /* determine axes */
  int x, y, z;
  int x_length, y_length, z_length;
  int x_axis, y_axis, z_axis;

  const char *c = strchr(axis, 'X');
  x = c - axis;

  const char *d = strchr(axis, 'Y');
  y = d - axis;

  const char *e = strchr(axis, 'Z');
  z = e - axis;

  /* TODO safeguard if one axis is not found */

  #ifdef T8_ENABLE_DEBUG
    t8_global_productionf("x: %d, y: %d, z: %d\n", x, y, z);
  #endif

  /* assign correct axis axis */
  x_axis =  x >= 0 && x <= 3 ? x : fmax(y, z) + 1;
  y_axis =  y >= 0 && y <= 3 ? y : fmax(x, z) + 1;
  z_axis =  z >= 0 && z <= 3 ? z : fmax(x, y) + 1;

  /* assign correct axis length */
  x_length = x >= 0 && x <= 3 ? shape[x] : 1;
  y_length = y >= 0 && y <= 3 ? shape[y] : 1;
  z_length = z >= 0 && z <= 3 ? shape[z] : 1;
  
  #ifdef T8_ENABLE_DEBUG
    t8_global_productionf("xaxis: %d, yaxis: %d, zaxis: %d\n", x_axis, y_axis, z_axis);
    t8_global_productionf("x_length: %d, y_length: %d, z_length: %d\n", x_length, y_length, z_length);
  #endif

  /* create forest for smallest mesh which completely contains given messy mesh */
  t8_forest_t forest = t8_latlon_refine(fmax(x_length, y_length), fmax(x_length, y_length), T8_LATLON_COARSEN, 0);
  t8_latlon_adapt_data_t *adapt_data =
    (t8_latlon_adapt_data_t *) t8_forest_get_user_data (forest);
  
  x_length = adapt_data->x_length;
  y_length = adapt_data->y_length;

  int* lshape = T8_ALLOC_ZERO(int, 3);
  memcpy(lshape, shape, sizeof(int) * 3);

  lshape[x_axis] = x_length;
  lshape[y_axis] = y_length;

  /* create data chunk */
  t8_latlon_data_chunk_t *chunk = t8_latlon_new_chunk(
    description,
    x_start, y_start,
    x_length, y_length, z_length,
    lshape,
    dimensions, x_axis, y_axis, z_axis, adapt_data->max_level,
    T8_LATLON_DATA_MESSY);

  t8_messy_data_t* messy_data = T8_ALLOC(t8_messy_data_t, 1);
  messy_data->chunk = chunk;
  messy_data->forest = forest;
  messy_data->coarsen = coarsen;
  messy_data->interpolation = interpolation;

  #ifdef T8_ENABLE_DEBUG
    t8_global_productionf("MESSy coupler initialized\n");
  #endif

  return messy_data;
}


void t8_messy_add_dimension(t8_messy_data_t *messy_data, char* dimension_name, double ****data) {
  t8_latlon_set_dimension(messy_data->chunk, dimension_name, data);
}

void t8_messy_set_dimension_values(t8_messy_data_t *messy_data, char* dimension_name, double *data) {
  t8_latlon_data_chunk_t *data_chunk = messy_data->chunk;
  int dimension_index = t8_latlon_get_dimension_idx(data_chunk, dimension_name, true);
  
  t8_debugf("set values for dimension %s at index %d\n", dimension_name, dimension_index);

  int size = data_chunk->x_length * data_chunk->y_length * data_chunk->z_length;
  int len = data_chunk->shape[0] * data_chunk->shape[1];
  int i, l, x, y, z, cell_index, data_index;
  int *idx = T8_ALLOC_ZERO(int, 3);

  for(i=0; i<size; ++i) {
    idx[0] = i / len;
    l      = i % len;
    idx[1] = l / data_chunk->shape[0];
    idx[2] = l % data_chunk->shape[0];
    
    x = idx[2 - data_chunk->x_axis];
    y = idx[2 - data_chunk->y_axis];
    z = idx[2 - data_chunk->z_axis];

    /* TODO!!!!!!! you know what */
    if(x < 64 && y < 32) {
      cell_index = y * data_chunk->x_length + x;
      data_index = ((y * data_chunk->z_length * data_chunk->x_length + x * data_chunk->z_length + z) * data_chunk->dimensions) + dimension_index;

      if(data_index >= 38912) {
        t8_debugf("(%d)[%d, %d, %d][%d] cell index: %d linear id: %d data index: %d \n", i, x, y, z, dimension_index, cell_index, t8_latlon_to_linear_id (x, y, data_chunk->level), data_index);
        
      }

      /* set data and morton id */
      //data_chunk->data[data_index] = 1; //data[i];
      //t8_debugf("(%d)[%d, %d, %d][%d] cell index: %d linear id: %d data index: %d \n", i, x, y, z, dimension_index, cell_index, t8_latlon_to_linear_id (x, y, data_chunk->level), data_index);
      memcpy((data_chunk->data) + data_index, data + i, sizeof(double));
    }
    
    //data_chunk->data_ids[cell_index] = t8_latlon_to_linear_id (x, y, data_chunk->level);
  }
}

void t8_messy_apply_sfc(t8_messy_data_t *messy_data) {
  t8_latlon_data_chunk_t *data_chunk = messy_data->chunk;
  t8_latlon_data_apply_morton_order(data_chunk);
}

void t8_messy_coarsen(t8_messy_data_t *messy_data) {
  char vtu_prefix[BUFSIZ];

  t8_global_productionf("MESSy coarsen grid \n");

  /* check if coarsening and interpolation configuration is set */
  T8_ASSERT(messy_data->coarsen != NULL);
  T8_ASSERT(messy_data->interpolation != NULL);

  t8_debugf("check if custom coarsen function is supplied \n");
  if(messy_data->coarsen->method == T8_MESSY_COARSEN_FUNCTION) {
    /* if coarsening method is custom function, check one is given */
    T8_ASSERT(messy_data->coarsen->func != NULL);
  }

  t8_debugf("check if custom interpolation function is supplied \n");
  if(messy_data->interpolation->method == T8_MESSY_INTERPOLATE_FUNCTION) {
    /* if interpolation method is custom function, check one is given */
    T8_ASSERT(messy_data->interpolation->func != NULL);
  }

  t8_debugf("extract pieces \n");
  t8_latlon_data_chunk_t *data_chunk = messy_data->chunk;

  t8_forest_t forest = messy_data->forest;
  t8_forest_t* forest_adapt = T8_ALLOC(t8_forest_t, 1);

  //messy_data->forest_adapt = forest_adapt;

  t8_debugf("ref forest \n");
  t8_forest_ref(forest);
  t8_debugf("init forest \n");
  t8_forest_init(forest_adapt);

  #ifdef T8_ENABLE_DEBUG
    /* In debugging mode write the forest */
    t8_debugf("write forest \n");
    snprintf (vtu_prefix, BUFSIZ, "t8_messy_grid");
    t8_messy_write_forest(forest, vtu_prefix, data_chunk);
  #endif

  t8_debugf("set user data \n");
  t8_forest_set_user_data(*forest_adapt, messy_data);
  t8_debugf("set adapt callback \n");
  t8_forest_set_adapt(*forest_adapt, forest, t8_messy_coarsen_callback, 0);

  t8_forest_set_partition (*forest_adapt, NULL, 0);

  t8_debugf("set commit forest \n");
  t8_forest_commit(*forest_adapt);
  t8_global_productionf("MESSy coarsen done\n");

  t8_global_productionf("MESSy interpolate grid \n");

  int num_element = t8_forest_get_num_element(*forest_adapt);

  data_chunk->data_ids_adapt = T8_ALLOC(t8_linearidx_t, num_element);
  data_chunk->data_adapt = T8_ALLOC(double, num_element * data_chunk->z_length * data_chunk->dimensions);

  t8_forest_iterate_replace(*forest_adapt, forest, t8_messy_interpolate_callback);

  T8_FREE(data_chunk->data_ids);
  T8_FREE(data_chunk->data);

  data_chunk->data_ids = data_chunk->data_ids_adapt;
  data_chunk->data = data_chunk->data_adapt;

  data_chunk->data_ids_adapt = NULL;
  data_chunk->data_adapt = NULL;

  #ifdef T8_ENABLE_DEBUG
    /* In debugging mode write the forest */
    snprintf (vtu_prefix, BUFSIZ, "t8_messy_grid_interpolated");
    t8_messy_write_forest(*forest_adapt, vtu_prefix, data_chunk);
  #endif

  t8_forest_unref (&forest);
  t8_global_productionf("MESSy grid interpolated \n");

}

void t8_messy_write_forest(t8_forest_t forest, const char* prefix, t8_latlon_data_chunk_t *data_chunk) {

  int num_elements = t8_forest_get_num_element(forest);
  int num_data = data_chunk->dimensions * data_chunk->z_length;

  t8_debugf("num elements %d, num data %d \n", num_elements, num_data);
  /* TODO: Do not use static array with variable as length */
  t8_vtk_data_field_t vtk_data[num_data];
  double *dim_data_array[num_data];
  
  int z, d, e, offset;
  for(z = 0; z < data_chunk->z_length; ++z) { 
    for(d = 0; d < data_chunk->dimensions; ++d) {
      offset = z * data_chunk->dimensions + d;
      
      dim_data_array[offset] = T8_ALLOC_ZERO (double, num_elements);
      for(e = 0; e < num_elements; ++e) {
        t8_debugf("el %d, index %d, value %.4f \n", e, e * num_data + z * data_chunk->dimensions + d, data_chunk->data[e * num_data + z * data_chunk->dimensions + d]);
        dim_data_array[offset][e] = data_chunk->data[e * num_data + z * data_chunk->dimensions + d];
      }
      snprintf (vtk_data[offset].description, BUFSIZ, "z%d_%s", z, data_chunk->dimension_names + d * BUFSIZ);
      vtk_data[offset].type = T8_VTK_SCALAR;
      vtk_data[offset].data = (double*) dim_data_array[offset];
    }
  }

  t8_forest_vtk_write_file (forest, prefix, 1, 1, 1, 1, 0, num_data, vtk_data);

  for(offset = 0; offset < num_data; ++offset) {
    T8_FREE(dim_data_array[offset]);
  }
}