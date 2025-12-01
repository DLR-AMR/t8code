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

T8_EXTERN_C_BEGIN ();

double
t8_norm (const double vec[3])
{
  const t8_3D_vec_view vec_view = make_t8_3D_vec_view (vec);
  return t8_norm (vec_view);
}

void
t8_normalize (double vec[3])
{
  t8_3D_vec_view vec_view = make_t8_3D_vec_view (vec);
  t8_normalize (vec_view);
}

void
t8_copy (const double dimensional_in[3], double dimensional_out[3])
{
  const t8_3D_vec_view vec_view_in = make_t8_3D_vec_view (dimensional_in);
  t8_3D_vec_view vec_view_out = make_t8_3D_vec_view (dimensional_out);
  t8_copy (vec_view_in, vec_view_out);
}

double
t8_dist (const double vec_x[3], const double vec_y[3])
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_dist (vec_x_view, vec_y_view);
}

void
t8_ax (double vec_x[3], const double alpha)
{
  t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  t8_ax (vec_x_view, alpha);
}

void
t8_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axy (vec_x_view, vec_y_view, alpha);
}

void
t8_axb (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axb (vec_x_view, vec_y_view, alpha, b);
}

void
t8_axpy (const double vec_x[3], double vec_y[3], const double alpha)
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axpy (vec_x_view, vec_y_view, alpha);
}

void
t8_axpyz (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  const t8_3D_vec_view vec_z_view = make_t8_3D_vec_view (vec_z);
  t8_axpyz (vec_x_view, vec_y_view, vec_z_view, alpha);
}

double
t8_dot (const double vec_x[3], const double vec_y[3])
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_dot (vec_x_view, vec_y_view);
}

void
t8_cross_3D (const double vec_x[3], const double vec_y[3], double cross[3])
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_3D_vec_view cross_view = make_t8_3D_vec_view (cross);
  t8_cross_3D (vec_x_view, vec_y_view, cross_view);
}

double
t8_cross_2D (const double vec_x[2], const double vec_y[2])
{
  const t8_2D_vec_view vec_x_view = make_t8_2D_vec_view (vec_x);
  const t8_2D_vec_view vec_y_view = make_t8_2D_vec_view (vec_y);
  return t8_cross_2D (vec_x_view, vec_y_view);
}

void
t8_diff (const double vec_x[3], const double vec_y[3], double diff[3])
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_3D_vec_view diff_view = make_t8_3D_vec_view (diff);
  t8_diff (vec_x_view, vec_y_view, diff_view);
}

int
t8_eq (const double vec_x[3], const double vec_y[3], const double tol)
{
  const t8_3D_vec_view vec_x_view = make_t8_3D_vec_view (vec_x);
  const t8_3D_vec_view vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_eq (vec_x_view, vec_y_view, tol);
}

void
t8_rescale (double vec[3], const double new_length)
{
  t8_3D_vec_view vec_view = make_t8_3D_vec_view (vec);
  t8_rescale (vec_view, new_length);
}

void
t8_normal_of_tri (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  const t8_3D_vec_view p1_view = make_t8_3D_vec_view (p1);
  const t8_3D_vec_view p2_view = make_t8_3D_vec_view (p2);
  const t8_3D_vec_view p3_view = make_t8_3D_vec_view (p3);
  const t8_3D_vec_view normal_view = make_t8_3D_vec_view (normal);
  t8_normal_of_tri (p1_view, p2_view, p3_view, normal_view);
}

void
t8_orthogonal_tripod (const double v1[3], double v2[3], double v3[3])
{
  const t8_3D_vec_view v1_view = make_t8_3D_vec_view (v1);
  t8_3D_vec_view v2_view = make_t8_3D_vec_view (v2);
  t8_3D_vec_view v3_view = make_t8_3D_vec_view (v3);
  t8_orthogonal_tripod (v1_view, v2_view, v3_view);
}

void
t8_swap (double p1[3], double p2[3])
{
  t8_3D_vec_view p1_view = make_t8_3D_vec_view (p1);
  t8_3D_vec_view p2_view = make_t8_3D_vec_view (p2);
  std::swap (p1_view, p2_view);
}

T8_EXTERN_C_END ();
