/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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

#include <t8_vec.hxx>
#include <t8_vec.h>

double
t8_vec_norm_c_interface (const double vec[3])
{
  return t8_vec_norm (vec);
}

void
t8_vec_normalize_c_interface (double vec[3])
{
  t8_vec_normalize (vec);
}

void
t8_vec_copy_c_interface (const double vec_in[3], double vec_out[3])
{
  t8_vec_copy (vec_in, vec_out);
}

double
t8_vec_dist_c_interface (const double vec_x[3], const double vec_y[3])
{
  return t8_vec_dist (vec_x, vec_y);
}

void
t8_vec_ax_c_interface (double vec_x[3], const double alpha)
{
  t8_vec_ax (vec_x, alpha);
}

void
t8_vec_axy_c_interface (const double vec_x[3], double vec_y[3], const double alpha)
{
  t8_vec_axy (vec_x, vec_y, alpha);
}

void
t8_vec_axb_c_interface (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  t8_vec_axb (vec_x, vec_y, alpha, b);
}

void
t8_vec_axpy_c_interface (const double vec_x[3], double vec_y[3], const double alpha)
{
  t8_vec_axpy (vec_x, vec_y, alpha);
}

void
t8_vec_axpyz_c_interface (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  t8_vec_axpyz (vec_x, vec_y, vec_z, alpha);
}

double
t8_vec_dot_c_interface (const double vec_x[3], const double vec_y[3])
{
  return t8_vec_dot (vec_x, vec_y);
}

void
t8_vec_cross_c_interface (const double vec_x[3], const double vec_y[3], double cross[3])
{
  t8_vec_cross (vec_x, vec_y, cross);
}

void
t8_vec_diff_c_interface (const double vec_x[3], const double vec_y[3], double diff[3])
{
  t8_vec_diff (vec_x, vec_y, diff);
}

int
t8_vec_eq_c_interface (const double vec_x[3], const double vec_y[3], const double tol)
{
  return t8_vec_eq (vec_x, vec_y, tol);
}

void
t8_vec_rescale_c_interface (double vec[3], const double new_length)
{
  t8_vec_rescale (vec, new_length);
}

void
t8_vec_tri_normal_c_interface (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  t8_vec_tri_normal (p1, p2, p3, normal);
}

void
t8_vec_orthogonal_tripod_c_interface (const double v1[3], double v2[3], double v3[3])
{
  t8_vec_orthogonal_tripod (v1, v2, v3);
}

void
t8_vec_swap_c_interface (double p1[3], double p2[3])
{
  t8_vec_swap (p1, p2);
}
