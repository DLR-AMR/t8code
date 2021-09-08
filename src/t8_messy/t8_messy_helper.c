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
#include <t8_messy/t8_messy_helper.h>

/* generate a random floating point number from min to max */
static double
t8_messy_randfrom (double min, double max)
{
  double              range = (max - min);
  double              div = RAND_MAX / range;
  return min + (rand () / div);
}

/**
 * Function filling data array with random values
 */
void
t8_messy_generate_data (double ****data, int x_length, int y_length,
                        double value)
{
  int                 x, y;
  for (y = 0; y < y_length; ++y) {
    for (x = 0; x < x_length; ++x) {
      data[x][y][0][0] = t8_messy_randfrom (y * x_length, (y + 1) * x_length);
    }
  }
}

void
t8_messy_sine_2d (double ****data, int x_length, int y_length)
{
  double              T_x = x_length / 5.0;
  double              T_y = y_length / 5.0;
  int                 x, y;
  for (y = 0; y < y_length; ++y) {
    for (x = 0; x < x_length; ++x) {
      data[x][y][0][0] = sin ((2 * M_PI * x) / T_x + (2 * M_PI * y) / T_y);
    }
  }
}

void
t8_messy_gaussian (double ****data, int x_length, int y_length)
{
  double              x0 = x_length * 1.0 / 2.0;
  double              y0 = y_length * 1.0 / 2.0;
  double              xd, yd, A, ox, oy;
  int                 x, y;

  ox = x0;
  oy = y0;
  A = 1;
  for (y = 0; y < y_length; ++y) {
    for (x = 0; x < x_length; ++x) {
      xd = x - x0;
      yd = y - y0;
      data[x][y][0][0] =
        A *
        exp (-
             (((xd * xd) / (2.0 * (ox * ox)) +
               ((yd * yd) / (2.0 * (oy * oy))))));
    }
  }
}
