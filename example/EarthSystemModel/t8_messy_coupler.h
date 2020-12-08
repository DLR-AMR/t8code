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


/* MESSy representation */
typedef struct
{
  const char         *name;       /* The name of this dataset. */
  t8_locidx_t         x_start;    /* Starting x coordinate. */
  t8_locidx_t         y_start;    /* Starting y coordinate. */
  t8_locidx_t         x_length;   /* Number of subgrid cells in x dimension. */
  t8_locidx_t         y_length;   /* Number of subgrid cells in y dimension. */
  int                 dimension;  /* Dimensionality of the data (1, 2, 3). */
  int                 axis;       /* Flag combining x, y and z axis */
  int                 x_axis;     /* X axis index in data vector */
  int                 y_axis;     /* Y axis index in data vector */
  int                 z_axis;     /* Z axis index in data vector */
} t8_messy_repr;


/* MESSy coupling object */
typedef struct {
  t8_messy_repr* repr
} t8_messy_data;


T8_EXTERN_C_BEGIN ();

/* Initialize forest for messy reprensentation */
t8_messy_data* t8_messy_initialize(
  const char* name,
  const char* axis,
  t8_locidx_t x_start, 
  t8_locidx_t y_start, 
  t8_locidx_t x_length, 
  t8_locidx_t y_length, 
  int dimension);

/* Add channel object data to forest */
void t8_messy_add_object();

T8_EXTERN_C_END ();

#endif /* !T8_MESSY_COUPLER_H */