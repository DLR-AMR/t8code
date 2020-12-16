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

/* MESSy coupling object */
typedef struct {
  t8_latlon_data_chunk_t *chunk;
  t8_forest_t forest;
  t8_forest_t forest_adapt;
} t8_messy_data;

/* Initialize forest for messy reprensentation */
t8_messy_data* t8_messy_initialize(
  const char* description,
  const char* axis,
  int x_start, 
  int y_start, 
  int x_length, 
  int y_length,
  int z_length,
  int dimensions);

T8_EXTERN_C_BEGIN ();

/* Add channel object data to forest */
void t8_messy_set_dimension(t8_messy_data *messy_data, double ****data, int dimension);


/* Bring input data into SFC format */
void t8_messy_apply_sfc(t8_messy_data *messy_data);

/* coarsen grid with given callback */
void t8_messy_coarsen(t8_messy_data *messy_data, t8_forest_adapt_t coarsen_callback, t8_forest_replace_t interpolate_callback);

T8_EXTERN_C_END ();

#endif /* !T8_MESSY_COUPLER_H */