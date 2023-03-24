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

#include <t8_vec.h>

double
t8_vec_norm (const double vec[3])
{
  double              norm = 0;

  for (int i = 0; i < 3; i++) {
    norm += vec[i] * vec[i];
  }
  return sqrt (norm);
}

double
t8_vec_dist (const double vec_x[3], const double vec_y[3])
{
  double              dist = 0;

  for (int i = 0; i < 3; i++) {
    dist += SC_SQR (vec_x[i] - vec_y[i]);
  }
  return sqrt (dist);
}

void
t8_vec_ax (double vec_x[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_x[i] *= alpha;
  }
}

void
t8_vec_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] = vec_x[i] * alpha;
  }
}

/* y = ax + b */
void
t8_vec_axb (const double vec_x[3], double vec_y[3], const double alpha,
            const double b)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] = alpha * vec_x[i] + b;
  }
}

/* y = y + alpha * x */
void
t8_vec_axpy (const double vec_x[3], double vec_y[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] += alpha * vec_x[i];
  }
}

/* z = y + alpha * x */
void
t8_vec_axpyz (const double vec_x[3], const double vec_y[3], double vec_z[3],
              const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_z[i] = vec_y[i] + alpha * vec_x[i];
  }
}

double
t8_vec_dot (const double vec_x[3], const double vec_y[3])
{
  double              dot = 0;

  for (int i = 0; i < 3; i++) {
    dot += vec_x[i] * vec_y[i];
  }
  return dot;
}

void
t8_vec_cross (const double vec_x[3], const double vec_y[3], double cross[3])
{
  for (int i = 0; i < 3; i++) {
    cross[i] =
      vec_x[(i + 1) % 3] * vec_y[(i + 2) % 3] -
      vec_x[(i + 2) % 3] * vec_y[(i + 1) % 3];
  }
}
