/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/** \file t8_example_common_functions.cxx Provide real valued functions
  * that are used in more than one example. */

#include <example/common/t8_example_common.h>
#include <t8_vec.h>

T8_EXTERN_C_BEGIN ();

double
t8_constant_one (const double x[3], double t)
{
  return 1;
}

double
t8_constant_zero (const double x[3], double t)
{
  return 0;
}

double
t8_project_x (const double x[3], double t)
{
  return x[0];
}

double
t8_exp_distribution (const double x[3], double t)
{
  double              dummy, X;

  /* Get fractional part of t. t is thus periodically
   * mapped to the unit interval */
  t = modf (t, &dummy);
  X = x[0] - .5;
  return exp (-4 * X * X);
}

/* This function is =1 if 0.25 <= x <= 0.75 and 0 else */
double
t8_step_function (const double x[3], double t)
{
  return 0.25 <= x[0] && x[0] <= 0.75;
}

/* This function is =1 if 0.25 <= x <= 0.75,
 * it is 0 outside of 0.25-eps and 0.75+eps,
 * it interpolates linearly in between. */
double
t8_almost_step_function (const double x[3], double t)
{
  double              eps = 0.1;

  /* interpolate in [0.25-eps,0.25+eps] */
  if (0.25 - eps < x[0] && x[0] < 0.25) {
    return (x[0] - 0.25 + eps) / (eps);
  }
  /* interpolate in [0.75-eps,0.75+eps] */
  else if (0.75 < x[0] && x[0] < 0.75 + eps) {
    return 1 - (x[0] - 0.75) / (eps);
  }
  /* 1 inside [0.25,0.75], 0 outside */
  return 0.25 <= x[0] && x[0] <= 0.75;
}

double
t8_sinx (const double x[3], double t)
{
  return sin (2 * M_PI * x[0]) + 1;
}

double
t8_sinx_cosy (const double x[3], double t)
{
  return sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]);
}

double
t8_sinx_cosy_z (const double x[3], double t)
{
  return 10 * sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]) * x[3];
}

double
t8_sint (const double x[3], double t)
{
  return sin (2 * M_PI * t);
}

double
t8_sphere_75_radius (const double x[3], double t)
{
  return t8_vec_norm (x) - 0.75;
}

double
t8_sphere_05_midpoint_375_radius (const double x[3], double t)
{
  double              M[3] = { 0.5, 0.5, 0.5 };

  /* Compute M - x */
  t8_vec_axpy (x, M, -1);

  /* return |M-x| - 0.375 */
  return t8_vec_norm (M) - 0.375;
}

double
t8_sphere_05_0z_midpoint_375_radius (const double x[3], double t)
{
  double              M[3] = { 0.5, 0.5, 0 };

  /* Compute M - x */
  t8_vec_axpy (x, M, -1);

  /* return |M-x| - 0.375 */
  return t8_vec_norm (M) - 0.375;
}

void
t8_constant_one_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = x_out[1] = x_out[2] = 1;
}

void
t8_constant_one_x_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = 1;
  x_out[1] = x_out[2] = 0;
}

void
t8_constant_one_xy_vec (const double x[3], double t, double x_out[3])
{
  x_out[0] = 1;
  x_out[1] = 1;
  x_out[2] = 0;
}

void
t8_rotation_2d (const double x_in[3], double t, double x_out[3])
{
  double              x = x_in[0], y = x_in[1];

  x -= 0.5;
  y -= 0.5;

  x_out[0] = y;
  x_out[1] = -x;
  x_out[2] = 0;
}

void
t8_compressible (const double x_in[3], double t, double x_out[3])
{
  x_out[0] = (1. / 2 - x_in[0]);
  x_out[1] = 0;
  x_out[2] = 0;
}

/* The following functions model a solution to the stokes equation on
 * a spherical shell. See
 * Analytical solution for viscous incompressible Stokes flow in a
 * spherical shell
 * by Cedric Thieulot
 */

static double
t8_stokes_sphere_alpha_beta (double R_1, double R_2, double gamma, int m,
                             double *alpha, double *beta)
{
  /* We define two constants alpha and beta */
  *alpha =
    gamma * (m + 1) * (pow (R_1, -3) - pow (R_2, -3)) / (pow (R_1, -m - 4) -
                                                         pow (R_2, -m - 4));
  *beta =
    -3 * gamma * (pow (R_1, m + 1) - pow (R_2, m + 1)) / (pow (R_1, m + 4) -
                                                          pow (R_2, m + 4));

}

/* A component of the flow that depends on the inner radius R_2, the outer radius R_1,
 * a constant gamma, and a control parameter m with m != -1, m != -4 */
static double
t8_stokes_sphere_g_component (double radius, double alpha, double beta,
                              double gamma, int m)
{
  T8_ASSERT (m != -1 && m != -4);

  return -2 / (radius * radius) * (-alpha / (m + 1) * pow (radius, -m - 1) +
                                   beta / 3 * pow (radius, 3) + gamma);
}

static double
t8_stokes_sphere_f_component (double radius, double alpha, double beta, int m)
{
  return alpha * pow (radius, -m - 3) + beta * radius;
}

void
t8_stokes_flow_sphere_shell (const double x[3], double t, double x_out[])
{
  double              radius;
  double              theta, phi;
  double              alpha, beta;
  double              vel_r;
  double              vel_theta;
  double              vel_phi;
  const double        r_1 = .5, r_2 = 1, gamma = 1, m = 3;

  /* Compute spherical coordinates */
  radius = t8_vec_norm (x);
  theta = acos (x[2] / radius);
  phi = atan2 (x[1], x[0]);

  if (radius < r_1) {
    /* If there are points in the geometry that lie in the inside radius,
     * set the flow to zero. */
    x_out[0] = x_out[1] = x_out[2] = 0;
    return;
  }

  /* Compute alpha and beta */
  t8_stokes_sphere_alpha_beta (r_1, r_2, gamma, m, &alpha, &beta);
  /* Compute radial velocity and theta velocity */
  vel_r =
    t8_stokes_sphere_g_component (radius, alpha, beta, gamma,
                                  m) * cos (theta);
  vel_theta =
    t8_stokes_sphere_f_component (radius, alpha, beta, m) * sin (theta);
  /* Set phi velocity */
  vel_phi = 0;

  /* Compute euclidean coordinates */
  x_out[0] = vel_r * sin (vel_theta) * cos (vel_phi);
  x_out[1] = vel_r * sin (vel_theta) * sin (vel_phi);
  x_out[2] = vel_r * cos (vel_theta);
}

T8_EXTERN_C_END ();
