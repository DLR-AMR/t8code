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

#include <t8_types/t8_vec.hxx>
#include <t8_types/t8_vec.h>
#include <memory>

double
t8_norm (const double vec[3])
{
  const t8_3D_vec *vec_array = reinterpret_cast<const t8_3D_vec *> (vec);
  return t8_norm (*vec_array);
}

void
t8_normalize (double vec[3])
{
  t8_3D_vec *vec_array = reinterpret_cast<t8_3D_vec *> (vec);
  t8_normalize (*vec_array);
}

void
t8_copy (const double vec_in[3], double vec_out[3])
{
  const t8_3D_vec *vec_array_in = reinterpret_cast<const t8_3D_vec *> (vec_in);
  t8_3D_vec *vec_array_out = reinterpret_cast<t8_3D_vec *> (vec_out);
  t8_copy (*vec_array_in, *vec_array_out);
}

double
t8_dist (const double vec_x[3], const double vec_y[3])
{
  const t8_3D_point *vec_array_x = reinterpret_cast<const t8_3D_point *> (vec_x);
  const t8_3D_point *vec_array_y = reinterpret_cast<const t8_3D_point *> (vec_y);
  return t8_dist (*vec_array_x, *vec_array_y);
}

void
t8_ax (double vec_x[3], const double alpha)
{
  t8_3D_vec *vec_array_x = reinterpret_cast<t8_3D_vec *> (vec_x);
  t8_ax (*vec_array_x, alpha);
}

void
t8_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  t8_3D_vec *vec_array_y = reinterpret_cast<t8_3D_vec *> (vec_y);
  t8_axy (*vec_array_x, *vec_array_y, alpha);
}

void
t8_axb_c_interface (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  t8_3D_vec *vec_array_y = reinterpret_cast<t8_3D_vec *> (vec_y);
  t8_axb (*vec_array_x, *vec_array_y, alpha, b);
}

void
t8_axpy_c_interface (const double vec_x[3], double vec_y[3], const double alpha)
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  t8_3D_vec *vec_array_y = reinterpret_cast<t8_3D_vec *> (vec_y);
  t8_axpy (*vec_array_x, *vec_array_y, alpha);
}

void
t8_axpyz_c_interface (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  const t8_3D_vec *vec_array_y = reinterpret_cast<const t8_3D_vec *> (vec_y);
  t8_3D_vec *vec_array_z = reinterpret_cast<t8_3D_vec *> (vec_z);
  t8_axpyz (*vec_array_x, *vec_array_y, *vec_array_z, alpha);
}

double
t8_dot_c_interface (const double vec_x[3], const double vec_y[3])
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  const t8_3D_vec *vec_array_y = reinterpret_cast<const t8_3D_vec *> (vec_y);
  return t8_dot (*vec_array_x, *vec_array_y);
}

void
t8_cross_3D_c_interface (const double vec_x[3], const double vec_y[3], double cross[3])
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  const t8_3D_vec *vec_array_y = reinterpret_cast<const t8_3D_vec *> (vec_y);
  t8_3D_vec *cross_array = reinterpret_cast<t8_3D_vec *> (cross);
  t8_cross_3D (*vec_array_x, *vec_array_y, *cross_array);
}

double
t8_cross_2D_c_interface (const double vec_x[2], const double vec_y[2])
{
  const t8_vec<2> *vec_array_x = reinterpret_cast<const t8_vec<2> *> (vec_x);
  const t8_vec<2> *vec_array_y = reinterpret_cast<const t8_vec<2> *> (vec_y);
  return t8_cross_2D (*vec_array_x, *vec_array_y);
}

void
t8_diff_c_interface (const double vec_x[3], const double vec_y[3], double diff[3])
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  const t8_3D_vec *vec_array_y = reinterpret_cast<const t8_3D_vec *> (vec_y);
  t8_3D_vec *diff_array = reinterpret_cast<t8_3D_vec *> (diff);
  t8_diff (*vec_array_x, *vec_array_y, *diff_array);
}

int
t8_eq_c_interface (const double vec_x[3], const double vec_y[3], const double tol)
{
  const t8_3D_vec *vec_array_x = reinterpret_cast<const t8_3D_vec *> (vec_x);
  const t8_3D_vec *vec_array_y = reinterpret_cast<const t8_3D_vec *> (vec_y);
  return t8_eq (*vec_array_x, *vec_array_y, tol);
}

void
t8_rescale_c_interface (double vec[3], const double new_length)
{
  t8_3D_vec *vec_array = reinterpret_cast<t8_3D_vec *> (vec);
  t8_rescale (*vec_array, new_length);
}

void
t8_normal_of_tri_c_interface (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  const t8_3D_vec *p1_array = reinterpret_cast<const t8_3D_vec *> (p1);
  const t8_3D_vec *p2_array = reinterpret_cast<const t8_3D_vec *> (p2);
  const t8_3D_vec *p3_array = reinterpret_cast<const t8_3D_vec *> (p3);
  t8_3D_vec *normal_array = reinterpret_cast<t8_3D_vec *> (normal);
  t8_normal_of_tri (*p1_array, *p2_array, *p3_array, *normal_array);
}

void
t8_orthogonal_tripod_c_interface (const double v1[3], double v2[3], double v3[3])
{
  const t8_3D_vec *v1_array = reinterpret_cast<const t8_3D_vec *> (v1);
  t8_3D_vec *v2_array = reinterpret_cast<t8_3D_vec *> (v2);
  t8_3D_vec *v3_array = reinterpret_cast<t8_3D_vec *> (v3);
  t8_orthogonal_tripod (*v1_array, *v2_array, *v3_array);
}

void
t8_swap_c_interface (double p1[3], double p2[3])
{
  t8_3D_vec *p1_array = reinterpret_cast<t8_3D_vec *> (p1);
  t8_3D_vec *p2_array = reinterpret_cast<t8_3D_vec *> (p2);
  std::swap (*p1_array, *p2_array);
}
