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

double
t8_norm_c_interface (const double vec[3])
{
  t8_3D_vec vec_array ({ vec[0], vec[1], vec[2] });
  return t8_norm (vec_array);
}

void
t8_normalize_c_interface (double vec[3])
{
  t8_3D_vec vec_array ({ vec[0], vec[1], vec[2] });
  t8_normalize (vec_array);
}

void
t8_copy_c_interface (const double vec_in[3], double vec_out[3])
{
  const t8_3D_vec vec_array_in ({ vec_in[0], vec_in[1], vec_in[2] });
  t8_3D_vec vec_array_out ({ vec_out[0], vec_out[1], vec_out[2] });
  t8_copy (vec_array_in, vec_array_out);
}

double
t8_dist_c_interface (const double vec_x[3], const double vec_y[3])
{
  t8_3D_point vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_point vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  return t8_dist (vec_array_x, vec_array_y);
}

void
t8_ax_c_interface (double vec_x[3], const double alpha)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_ax (vec_array_x, alpha);
}

void
t8_axy_c_interface (const double vec_x[3], double vec_y[3], const double alpha)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_axy (vec_array_x, vec_array_y, alpha);
}

void
t8_axb_c_interface (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_axb (vec_array_x, vec_array_y, alpha, b);
}

void
t8_axpy_c_interface (const double vec_x[3], double vec_y[3], const double alpha)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_axpy (vec_array_x, vec_array_y, alpha);
}

void
t8_axpyz_c_interface (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_3D_vec vec_array_z ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_axpyz (vec_array_x, vec_array_y, vec_array_z, alpha);
}

double
t8_dot_c_interface (const double vec_x[3], const double vec_y[3])
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  return t8_dot (vec_array_x, vec_array_y);
}

void
t8_cross_3D_c_interface (const double vec_x[3], const double vec_y[3], double cross[3])
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_3D_vec cross_array ({ cross[0], cross[1], cross[2] });
  t8_cross_3D (vec_array_x, vec_array_y, cross_array);
}

double
t8_cross_2D_c_interface (const double vec_x[2], const double vec_y[2])
{
  t8_vec<2> vec_array_x ({ vec_x[0], vec_x[1] });
  t8_vec<2> vec_array_y ({ vec_y[0], vec_y[1] });
  return t8_cross_2D (vec_array_x, vec_array_y);
}

void
t8_diff_c_interface (const double vec_x[3], const double vec_y[3], double diff[3])
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  t8_3D_vec diff_array ({ diff[0], diff[1], diff[2] });
  t8_diff (vec_array_x, vec_array_y, diff_array);
}

int
t8_eq_c_interface (const double vec_x[3], const double vec_y[3], const double tol)
{
  t8_3D_vec vec_array_x ({ vec_x[0], vec_x[1], vec_x[2] });
  t8_3D_vec vec_array_y ({ vec_y[0], vec_y[1], vec_y[2] });
  return t8_eq (vec_array_x, vec_array_y, tol);
}

void
t8_rescale_c_interface (double vec[3], const double new_length)
{
  t8_3D_vec vec_array ({ vec[0], vec[1], vec[2] });
  t8_rescale (vec_array, new_length);
}

void
t8_normal_of_tri_c_interface (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  t8_3D_vec p1_array ({ p1[0], p1[1], p1[2] });
  t8_3D_vec p2_array ({ p2[0], p2[1], p2[2] });
  t8_3D_vec p3_array ({ p3[0], p3[1], p3[2] });
  t8_3D_vec normal_array ({ normal[0], normal[1], normal[2] });
  t8_normal_of_tri (p1_array, p2_array, p3_array, normal_array);
}

void
t8_orthogonal_tripod_c_interface (const double v1[3], double v2[3], double v3[3])
{
  t8_3D_vec v1_array ({ v1[0], v1[1], v1[2] });
  t8_3D_vec v2_array ({ v2[0], v2[1], v2[2] });
  t8_3D_vec v3_array ({ v3[0], v3[1], v3[2] });
  t8_orthogonal_tripod (v1_array, v1_array, v1_array);
}

void
t8_swap_c_interface (double p1[3], double p2[3])
{
  t8_3D_vec p1_array ({ p1[0], p1[1], p1[2] });
  t8_3D_vec p2_array ({ p2[0], p2[1], p2[2] });
  std::swap (p1_array, p2_array);
}
