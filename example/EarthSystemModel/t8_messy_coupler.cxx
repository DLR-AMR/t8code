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
#include "t8_latlon_refine.h"
#include "t8_messy_coupler.h"

t8_messy_data* t8_messy_initialize(
  char* name,
  char* axis,
  t8_locidx_t x_start, 
  t8_locidx_t y_start, 
  t8_locidx_t x_length, 
  t8_locidx_t y_length, 
  int dimension) {

  t8_messy_data messy_data;
  t8_messy_repr repr;

  repr.name = name;
  repr.x_length = x_length;
  repr.y_length = y_length;
  repr.x_start = x_start;
  repr.y_start = y_start;
  repr.dimension = dimension;

  // determine axes
  char *c;
  int x, y, z;

  c = strchr(name, 'x');
  x = (int)(c - name);

  c = strchr(name, 'y');
  y = (int)(c - name);

  c = strchr(name, 'z');
  z = (int)(c - name);

  repr.x_axis = x;
  repr.y_axis = y;
  repr.z_axis = z;

  messy_data.repr = *repr;

  return *messy_data;
}
