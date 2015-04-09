/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <t8_type.h>

/* *INDENT-OFF* */
const int t8_type_to_dimension[T8_TYPE_LAST] =
  { 0, 1, 2, 2, 3, 3, 3, 3 };

const int t8_type_boundary_count[T8_TYPE_LAST][T8_TYPE_LAST] =
  {{ 0,  0, 0, 0, 0, 0, 0, 0 },
   { 2,  0, 0, 0, 0, 0, 0, 0 },
   { 4,  4, 0, 0, 0, 0, 0, 0 },
   { 3,  3, 0, 0, 0, 0, 0, 0 },
   { 8, 12, 6, 0, 0, 0, 0, 0 },
   { 4,  6, 0, 4, 0, 0, 0, 0 },
   { 6,  9, 3, 2, 0, 0, 0, 0 },
   { 5,  8, 1, 4, 0, 0, 0, 0 }};
/* *INDENT-ON* */

int
t8_type_count_boundary (t8_type_t thetype, int min_dim, int *per_type)
{
  int                 t;
  int                 sum;

  sum = 0;
  for (t = 0; t < T8_TYPE_LAST; ++t) {
    if (t8_type_to_dimension[t] >= min_dim) {
      sum += (per_type[t] = t8_type_boundary_count[thetype][t]);
    }
    else {
      per_type[t] = 0;
    }
  }

  return sum;
}
