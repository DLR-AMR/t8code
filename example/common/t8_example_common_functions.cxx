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

#include <example/common/t8_example_common.hxx>
#include <t8_types/t8_vec.hxx>

double
t8_levelset_sphere (const t8_3D_point &x, [[maybe_unused]] const double t, void *data)
{
  t8_levelset_sphere_data_t *ls_data = (t8_levelset_sphere_data_t *) data;

  T8_ASSERT (ls_data->radius > 0);
  return t8_dist (x, ls_data->M) - ls_data->radius;
}

double
t8_scalar3d_constant_one ([[maybe_unused]] const double x[3], [[maybe_unused]] double t)
{
  return 1;
}

double
t8_scalar3d_constant_zero ([[maybe_unused]] const double x[3], [[maybe_unused]] double t)
{
  return 0;
}

double
t8_scalar3d_project_x (const double x[3], [[maybe_unused]] double t)
{
  return x[0];
}

double
t8_scalar3d_exp_distribution (const double x[3], double t)
{
  double dummy, X;

  /* Get fractional part of t. t is thus periodically
   * mapped to the unit interval */
  t = modf (t, &dummy);
  X = x[0] - .5;
  return exp (-4 * X * X);
}

/* This function is =1 if 0.25 <= x <= 0.75 and 0 else */
double
t8_scalar3d_step_function (const double x[3], [[maybe_unused]] double t)
{
  return 0.25 <= x[0] && x[0] <= 0.75;
}

/* This function is =1 if 0.25 <= x <= 0.75,
 * it is 0 outside of 0.25-eps and 0.75+eps,
 * it interpolates linearly in between. */
double
t8_scalar3d_almost_step_function (const double x[3], [[maybe_unused]] double t)
{
  double eps = 0.1;

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
t8_scalar3d_sinx (const double x[3], [[maybe_unused]] double t)
{
  return sin (2 * M_PI * x[0]) + 1;
}

double
t8_scalar3d_sinx_cosy (const double x[3], [[maybe_unused]] double t)
{
  return sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]);
}

double
t8_scalar3d_sinx_cosy_z (const double x[3], [[maybe_unused]] double t)
{
  return 10 * sin (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]) * x[3];
}

double
t8_scalar3d_sint ([[maybe_unused]] const double x[3], double t)
{
  return sin (2 * M_PI * t);
}

/* general level set function for a sphere with given midpoint and radius. */
static double
t8_scalar3d_sphere (const t8_3D_vec &x, t8_3D_vec &M, const double radius)
{

  /* Compute M - x */
  t8_axpy (x, M, -1);

  /* return |M-x| - radius */

  return t8_norm (M) - radius;
}

double
t8_scalar3d_sphere_75_radius (const t8_3D_vec x, [[maybe_unused]] const double t)
{
  t8_3D_vec M ({ 0, 0, 0 });
  return t8_scalar3d_sphere (x, M, 0.75);
}

double
t8_scalar3d_sphere_05_midpoint_375_radius (const t8_3D_vec x, [[maybe_unused]] double t)
{
  t8_3D_vec M ({ 0.5, 0.5, 0.5 });

  return t8_scalar3d_sphere (x, M, 0.375);
}

double
t8_scalar3d_sphere_03_midpoint_25_radius (const t8_3D_vec &x, [[maybe_unused]] const double t)
{
  t8_3D_vec M ({ 0.3, 0.3, 0.3 });

  return t8_scalar3d_sphere (x, M, 0.25);
}

double
t8_scalar3d_sphere_05_0z_midpoint_375_radius (const t8_3D_vec &x, [[maybe_unused]] const double t)
{
  t8_3D_vec M ({ 0.5, 0.5, 0 });

  return t8_scalar3d_sphere (x, M, 0.375);
}

void
t8_flow_constant_one_vec ([[maybe_unused]] const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  x_out[0] = x_out[1] = x_out[2] = 1;
}

void
t8_flow_constant_one_x_vec ([[maybe_unused]] const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  x_out[0] = 1;
  x_out[1] = x_out[2] = 0;
}

void
t8_flow_constant_one_xy_vec ([[maybe_unused]] const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  x_out[0] = 1;
  x_out[1] = 0.8;
  x_out[2] = 0;
}

void
t8_flow_constant_one_xyz_vec ([[maybe_unused]] const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  x_out[0] = 1;
  x_out[1] = 0.8;
  x_out[2] = 0.9;
}

void
t8_flow_rotation_2d (const t8_3D_point &x_in, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  double x = x_in[0], y = x_in[1];

  x -= 0.5;
  y -= 0.5;

  x_out[0] = y;
  x_out[1] = -x;
  x_out[2] = 0;

  t8_ax (x_out, 2 * M_PI);
}

void
t8_flow_compressible (const t8_3D_point &x_in, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  x_out[0] = (1. / 2 - x_in[0]);
  x_out[1] = 0;
  x_out[2] = 0;
}

/* The following function is a incompressible flow on the unit cube.
 * It is constructed from any function f with f(0) = f(1) = 0.
 */

static double
t8_incomp_cube_f_sin (double x)
{
  return sin (M_PI * x);
}

static double
t8_incomp_cube_df_sin (double x)
{
  return M_PI * cos (M_PI * x);
}

void
t8_flow_incomp_cube_flow (const t8_3D_point &x, const double t, t8_3D_vec &x_out)
{
  double (*f) (double) = t8_incomp_cube_f_sin;
  double (*df) (double) = t8_incomp_cube_df_sin;

  x_out[0] = f (x[0]) * (df (x[1]) - df (x[2]));
  x_out[1] = -1. * f (x[1]) * df (x[0]);
  x_out[2] = f (x[2]) * df (x[0]);

  t8_ax (x_out, 1. / 2);
  /* We reverse the flow at time 0.5 */
  if (t > 0.5) {
    t8_ax (x_out, -1);
  }
}

/* Convert the first two entries of a vector into 2D polar
 * coordinates. x = r cos(phi)
 *              y = r sin(phi)
 *
 * On output: polar[0] = r, polar[1] = phi
 */
static void
t8_flow_2d_polar_coords (const t8_2D_vec &x, t8_2D_vec &polar)
{
  polar[0] = sqrt (SC_SQR (x[0]) + SC_SQR (x[1]));
  polar[1] = atan2 (x[1], x[0]);
}

/* Convert a 2D vector from polar coordinates to cartesian
 * coordinates.
 * On input: polar[0] = r, polar[1] = phi
 *
 * On output: cart[0] = r cos(phi)
 *            cart[1] = r sin(phi)
 *
 */
static void
t8_flow_2d_cart_coords (const t8_2D_vec &polar_values, const t8_2D_vec &polar_coords, t8_2D_vec &cart)
{
  cart[0] = cos (polar_coords[1]) * polar_values[0] - sin (polar_coords[1]) * polar_values[1];
  cart[1] = sin (polar_coords[1]) * polar_values[0] + cos (polar_coords[1]) * polar_values[1];
}

/* 2d flow around a circle with radius R = 1 and
 * constant inflow with x-speed U = 1. */
void
t8_flow_around_circle (const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  t8_2D_vec polar;
  t8_2D_vec polar_speed;
  const t8_2D_vec x_2D = t8_2D_vec ({ x[0], x[1] });
  t8_2D_vec x_out_2D ({ x_out[0], x_out[1] });
  const double R = 0.15;

  t8_axb (x_2D, x_out_2D, 1, -0.5);
  /* Set the z-coordinate to zero */
  if (t8_norm (x_out) < R) {
    /* Set the velocity inside the circle to 0 */
    x_out[0] = x_out[1] = 0;
    return;
  }
  /* Convert x,y coordinates to polar r,phi coordinates */
  t8_flow_2d_polar_coords (x_out_2D, polar);
  /* Compute v_r (r,phi) = U (1-R^2/r^2)cos(phi) */
  polar_speed[0] = (1 - SC_SQR (R) / SC_SQR (polar[0])) * cos (polar[1]);
  /* Compute v_phi(r,phi) = -U (1+ R^2/r^2) sin (phi) */
  polar_speed[1] = -(1 + SC_SQR (R) / SC_SQR (polar[0])) * sin (polar[1]);
  t8_flow_2d_cart_coords (polar_speed, polar, x_out_2D);
  x_out[0] = x_out_2D[0];
  x_out[1] = x_out_2D[1];
  x_out[2] = 0;
}

/* The following functions model a solution to the stokes equation on
 * a spherical shell. See
 * Analytical solution for viscous incompressible Stokes flow in a
 * spherical shell
 * by Cedric Thieulot
 */

static void
t8_flow_stokes_sphere_alpha_beta (double R_1, double R_2, double gamma, int m, double *alpha, double *beta)
{
  /* We define two constants alpha and beta */
  *alpha = gamma * (m + 1) * (pow (R_1, -3) - pow (R_2, -3)) / (pow (R_1, -m - 4) - pow (R_2, -m - 4));
  *beta = -3 * gamma * (pow (R_1, m + 1) - pow (R_2, m + 1)) / (pow (R_1, m + 4) - pow (R_2, m + 4));
}

/* A component of the flow that depends on the inner radius R_2, the outer radius R_1,
 * a constant gamma, and a control parameter m with m != -1, m != -4 */
static double
t8_flow_stokes_sphere_g_component (double radius, double alpha, double beta, double gamma, int m)
{
  T8_ASSERT (m != -1 && m != -4);

  return -2 / (radius * radius) * (-alpha / (m + 1) * pow (radius, -m - 1) + beta / 3 * pow (radius, 3) + gamma);
}

static double
t8_flow_stokes_sphere_f_component (double radius, double alpha, double beta, int m)
{
  return alpha * pow (radius, -m - 3) + beta * radius;
}

void
t8_flow_stokes_flow_sphere_shell (const t8_3D_point &x_in, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  double radius;
  double theta, phi;
  double alpha, beta;
  double vel_r;
  double vel_theta;
  double vel_phi;
  const double r_1 = .5, r_2 = 1, gamma = 1, m = 3;
  /* translate unit cube to cube centered around origin */
  t8_3D_vec x ({ (x_in[0] - 0.5 * 2), (x_in[1] - 0.5 * 2), (x_in[2] - 0.5 * 2) });

  /* Compute spherical coordinates */
  radius = t8_norm (x);
  theta = acos (x[2] / radius);
  /* Phi component, not used */
  phi = atan2 (x[1], x[0]);

  if (radius < r_1) {
    /* If there are points in the geometry that lie in the inside radius,
     * set the flow to zero. */
    x_out[0] = x_out[1] = x_out[2] = 0;
    return;
  }

  /* Compute alpha and beta */
  t8_flow_stokes_sphere_alpha_beta (r_1, r_2, gamma, m, &alpha, &beta);
  /* Compute radial velocity and theta velocity */
  vel_r = t8_flow_stokes_sphere_g_component (radius, alpha, beta, gamma, m) * cos (theta);
  vel_theta = t8_flow_stokes_sphere_f_component (radius, alpha, beta, m) * sin (theta);
  /* Set phi velocity */
  vel_phi = 0;

  /* Compute euclidean coordinates */
  x_out[0] = vel_r * sin (theta) * cos (phi) + vel_theta * cos (theta) * cos (phi) - vel_phi * sin (phi);
  x_out[1] = vel_r * sin (theta) * sin (phi) + vel_theta * cos (theta) * sin (phi) + vel_phi * cos (phi);
  x_out[2] = vel_r * cos (theta) - vel_theta * cos (theta);
}

void
t8_flow_around_circle_with_angular_velocity (const t8_3D_point &x, [[maybe_unused]] const double t, t8_3D_vec &x_out)
{
  const double radius = 0.5;
  const double omega = 1.5 * M_PI;
  // convert to polar coordinates
  const double r = sqrt (x[0] * x[0] + x[1] * x[1]);
  const double theta = atan2 (x[1], x[0]);

  //calculate flow
  const double u_r = (1 - (radius * radius) / (r * r)) * cos (theta);
  const double u_theta = -(1 + (radius * radius) / (r * r)) * sin (theta) - omega / (2 * M_PI * r);

  // convert back to Cartesian
  x_out[0] = cos (theta) * u_r - sin (theta) * u_theta;
  x_out[1] = sin (theta) * u_r + cos (theta) * u_theta;
  x_out[2] = 0;
}
