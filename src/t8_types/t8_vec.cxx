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
  const auto vec_view = make_t8_3D_vec_view (vec);
  return t8_norm (vec_view);
}

void
t8_normalize (double vec[3])
{
  auto vec_view = make_t8_3D_vec_view (vec);
  t8_normalize (vec_view);
}

void
t8_copy (const double dimensional_in[3], double dimensional_out[3])
{
  const auto vec_view_in = make_t8_3D_vec_view (dimensional_in);
  auto vec_view_out = make_t8_3D_vec_view (dimensional_out);
  t8_copy (vec_view_in, vec_view_out);
}

double
t8_dist (const double vec_x[3], const double vec_y[3])
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_dist (vec_x_view, vec_y_view);
}

void
t8_ax (double vec_x[3], const double alpha)
{
  auto vec_x_view = make_t8_3D_vec_view (vec_x);
  t8_ax (vec_x_view, alpha);
}

void
t8_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  auto vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axy (vec_x_view, vec_y_view, alpha);
}

void
t8_axb (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  auto vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axb (vec_x_view, vec_y_view, alpha, b);
}

void
t8_axpy (const double vec_x[3], double vec_y[3], const double alpha)
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  auto vec_y_view = make_t8_3D_vec_view (vec_y);
  t8_axpy (vec_x_view, vec_y_view, alpha);
}

void
t8_axpyz (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  const auto vec_z_view = make_t8_3D_vec_view (vec_z);
  t8_axpyz (vec_x_view, vec_y_view, vec_z_view, alpha);
}

double
t8_dot (const double vec_x[3], const double vec_y[3])
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_dot (vec_x_view, vec_y_view);
}

void
t8_cross_3D (const double vec_x[3], const double vec_y[3], double cross[3])
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  auto cross_view = make_t8_3D_vec_view (cross);
  t8_cross_3D (vec_x_view, vec_y_view, cross_view);
}

double
t8_cross_2D (const double vec_x[2], const double vec_y[2])
{
  const auto vec_x_view = make_t8_2D_vec_view (vec_x);
  const auto vec_y_view = make_t8_2D_vec_view (vec_y);
  return t8_cross_2D (vec_x_view, vec_y_view);
}

void
t8_diff (const double vec_x[3], const double vec_y[3], double diff[3])
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  auto diff_view = make_t8_3D_vec_view (diff);
  t8_diff (vec_x_view, vec_y_view, diff_view);
}

int
t8_eq (const double vec_x[3], const double vec_y[3], const double tol)
{
  const auto vec_x_view = make_t8_3D_vec_view (vec_x);
  const auto vec_y_view = make_t8_3D_vec_view (vec_y);
  return t8_eq (vec_x_view, vec_y_view, tol);
}

void
t8_rescale (double vec[3], const double new_length)
{
  auto vec_view = make_t8_3D_vec_view (vec);
  t8_rescale (vec_view, new_length);
}

void
t8_normal_of_tri (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  const auto p1_view = make_t8_3D_vec_view (p1);
  const auto p2_view = make_t8_3D_vec_view (p2);
  const auto p3_view = make_t8_3D_vec_view (p3);
  const auto normal_view = make_t8_3D_vec_view (normal);
  t8_normal_of_tri (p1_view, p2_view, p3_view, normal_view);
}

void
t8_orthogonal_tripod (const double v1[3], double v2[3], double v3[3])
{
  const auto v1_view = make_t8_3D_vec_view (v1);
  auto v2_view = make_t8_3D_vec_view (v2);
  auto v3_view = make_t8_3D_vec_view (v3);
  t8_orthogonal_tripod (v1_view, v2_view, v3_view);
}

void
t8_swap (double p1[3], double p2[3])
{
  auto p1_view = make_t8_3D_vec_view (p1);
  auto p2_view = make_t8_3D_vec_view (p2);
  std::swap (p1_view, p2_view);
}

/**
 * Test whether four given points in 3D are coplanar up to a given tolerance.
 * \param [in]  p_0         First point to check.
 * \param [in]  p_1         Second point to check.
 * \param [in]  p_2         Third point to check.
 * \param [in]  p_3         Fourth point to check.
 * \param [in]  tolerance   The tolerance
 * \return true if points are coplanar.
 */
int
t8_four_points_coplanar (const double p_0[3], const double p_1[3], const double p_2[3], const double p_3[3],
                         const double tolerance)
{
  /* Let p0, p1, p2, p3 be the four points.
   * The four points are coplanar if the normal vectors to the triangles
   * p0, p1, p2 and p0, p2, p3 are pointing in the same direction.
   *
   * We build the vectors A = p1 - p0, B = p2 - p0 and C = p3 - p0.
   * The normal vectors to the triangles are n1 = A x B and n2 = A x C.
   * These are pointing in the same direction if their cross product is 0.
   * Hence we check if || n1 x n2 || < tolerance. */

  /* A = p1 - p0 */
  double A[3];
  t8_axpyz (p_0, p_1, A, -1);

  /* B = p2 - p0 */
  double B[3];
  t8_axpyz (p_0, p_2, B, -1);

  /* C = p3 - p0 */
  double C[3];
  t8_axpyz (p_0, p_3, C, -1);

  /* n1 = A x B */
  double A_cross_B[3];
  t8_cross_3D (A, B, A_cross_B);

  /* n2 = A x C */
  double A_cross_C[3];
  t8_cross_3D (A, C, A_cross_C);

  /* n1 x n2 */
  double n1_cross_n2[3];
  t8_cross_3D (A_cross_B, A_cross_C, n1_cross_n2);

  /* || n1 x n2 || */
  const double norm = t8_norm (n1_cross_n2);
  return norm < tolerance;
}

T8_EXTERN_C_END ();
