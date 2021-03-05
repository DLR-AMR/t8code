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
#ifndef T8_MESSY_COUPLER_H
#define T8_MESSY_COUPLER_H

#include <t8.h>
#include <t8_forest.h>
#include "t8_latlon_data.h"

typedef struct {
  int     num_elements;
  int     z_layer;
  int*    x_coords;
  int*    y_coords;
  double* latitudes;
  double* longitudes;
  double* values;
  char*   tracer;
} t8_messy_custom_func_t;

typedef enum {
  T8_MESSY_COARSEN_THRESHOLD_MIN_LOWER,
  T8_MESSY_COARSEN_THRESHOLD_MIN_HIGHER,
  T8_MESSY_COARSEN_THRESHOLD_MAX_LOWER,
  T8_MESSY_COARSEN_THRESHOLD_MAX_HIGHER,
  T8_MESSY_COARSEN_THRESHOLD_MEAN_LOWER,
  T8_MESSY_COARSEN_THRESHOLD_MEAN_HIGHER,
  T8_MESSY_COARSEN_AREA_INSIDE,
  T8_MESSY_COARSEN_AREA_OUTSIDE,
  /* custom function determine weather to coarsen or not
   * @param t8_messy_custom_func_t* 
   * @return -1 (coarsen) or 0 (not coarsen) 
   */
  T8_MESSY_COARSEN_FUNCTION 
} T8_MESSY_COARSEN_METHOD;

typedef struct {
  T8_MESSY_COARSEN_METHOD method; /* method used for coarsen */
  char* tracer;  /* tracer by which to coarsen */
  int z_layer;
  double threshold; /* threshold for threshold coarsening */
  double *points;   /* points array for area coarsening */
  int (*func)(t8_messy_custom_func_t *) = NULL;
} t8_messy_coarsen_t;

typedef enum {
  /* custom function calculating interpolatet value
   * @param  t8_messy_custom_func_t*
   * @return (double) interpolated value
   */
  T8_MESSY_INTERPOLATE_FUNCTION,
  T8_MESSY_INTERPOLATE_MIN,
  T8_MESSY_INTERPOLATE_MAX,
  T8_MESSY_INTERPOLATE_MEAN
} T8_MESSY_INTERPOLATE_METHOD;

typedef struct {
  T8_MESSY_INTERPOLATE_METHOD method; /* method used for interpolation */
  double (*func)(t8_messy_custom_func_t *) = NULL;
} t8_messy_interpolate_t;

/* MESSy coupling object */
typedef struct t8_messy_data {
  t8_latlon_data_chunk_t *chunk;
  t8_messy_coarsen_t *coarsen;
  t8_messy_interpolate_t *interpolation;
  t8_forest_t forest;
  double* errors;
  double* errors_adapt;
  int counter;
} t8_messy_data_t;


T8_EXTERN_C_BEGIN ();


t8_messy_coarsen_t* t8_messy_new_coarsen_config(  
  const char* method,
  char* tracer,
  int z_layer,
  double threshold,
  int (*func)(t8_messy_custom_func_t *)
);

t8_messy_interpolate_t* t8_messy_new_interpolate_config(
  const char* method,
  double (*func)(t8_messy_custom_func_t *)
);

t8_messy_custom_func_t* t8_messy_new_custom_func(int num_elements);

void t8_messy_destroy_custom_func(t8_messy_custom_func_t* custom);

/* Initialize forest for messy reprensentation */
t8_messy_data_t* t8_messy_initialize(
  const char* description,
  const char* axis,
  int* shape, 
  int x_start, 
  int y_start, 
  int num_tracers,
  t8_messy_coarsen_t *coarsen,
  t8_messy_interpolate_t *interpolation);

void t8_messy_reset(t8_messy_data_t* messy_data);

/* set tracer values */
void t8_messy_set_tracer_values(t8_messy_data_t *messy_data, char* tracer_name, double *data);


/* Bring input data into SFC format */
void t8_messy_apply_sfc(t8_messy_data_t *messy_data);

/* coarsen grid with given callback */
void t8_messy_coarsen(t8_messy_data_t *messy_data);

void t8_messy_write_forest(t8_forest_t forest, const char* prefix, t8_messy_data_t* messy_data);

void t8_messy_destroy(t8_messy_data_t* messy_data);

T8_EXTERN_C_END ();

#endif /* !T8_MESSY_COUPLER_H */