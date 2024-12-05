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

/** file t8_example_common.hxx
 * This header declares datatypes and functions that are used by multiple
 * examples of t8code. This includes adaptation function to adapt at a zero
 * level-set of a level-set function and various 2 and 3 dimensional vector
 * fields that can serve as fluid velocities.
 */

#ifndef T8_EXAMPLE_COMMON_H
#define T8_EXAMPLE_COMMON_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_vec.hxx>

/** A levelset function in 3+1 space dimensions. */
typedef double (*t8_example_level_set_fn) (const t8_3D_point &, double, void *);

/** Struct to handle refinement around a level-set function. */
typedef struct
{
  t8_example_level_set_fn L; /**< The level set function. */
  void *udata;               /**< Data pointer that is passed to L */
  double band_width;         /**< Width of max_level elements around the zero-level set */
  double t;                  /**< Time value passed to levelset function */
  int min_level;             /**< The minimal refinement level. Elements with this level will not be coarsened. */
  int max_level;             /**< The maximum refinement level. Elements with this level will not be refined. */
} t8_example_level_set_struct_t;

/** Function pointer for real valued functions from d+1 space dimensions
 * functions f: R^d x R -> R */
typedef double (*t8_scalar_function_1d_fn) (double x, double t);
typedef double (*t8_scalar_function_2d_fn) (const t8_point<2> &x, double t);
typedef double (*t8_scalar_function_3d_fn) (const t8_3D_point &x, double t);
/** Function pointer for a vector valued function f: R^3 x R -> R */
typedef void (*t8_flow_function_3d_fn) (const t8_3D_point &x_in, double t, t8_3D_vec &x_out);

/* function declarations */

/** Query whether a given element is within a prescribed distance to the zero level-set
 *  of a level-set function.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     A local tree in \a forest.
 * \param [in]      element     An element of tree \a ltreeid in \a forest.
 * \param [in]      levelset    The level-set function.
 * \param [in]      band_width  Check whether the element is within a band of
 *                              \a band_width many elements of its size.
 * \param [in]      t           Time value passed to \a levelset.
 * \param [in]      udata       User data passed to \a levelset.
 * \return                      True if the absolute value of \a levelset at \a element's midpoint
 *                              is smaller than \a band_width * \a element's diameter.
 *                              False otherwise.
 *                              If \a band_width = 0 then the return value is true if and only if
 *                              the zero level-set passes through \a element.
 */
int
t8_common_within_levelset (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                           t8_example_level_set_fn levelset, double band_width, double t, void *udata);

/** Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 */
int
t8_common_adapt_balance (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                         const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                         const int is_family, const int num_elements, t8_element_t *elements[]);

/** Adapt a forest along a given level-set function.
 * The user data of forest must be a pointer to a \a t8_example_level_set_struct_t.
 * An element in the forest is refined, if it is in a band of \a band_with many
 * \a max_level elements around the zero level-set Gamma = { x | L(x) = 0}
 */
/* TODO: Currently the band_width control is not working yet.
 *        if band_with = 0, then all elements that are touched by the zero LS are refined. */
int
t8_common_adapt_level_set (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                           const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                           const int is_family, const int num_elements, t8_element_t *elements[]);

/** Real valued functions defined in t8_example_common_functions.h */

typedef struct
{
  t8_3D_point M;
  /**< midpoint */
  double radius;
  /**< radius */
} t8_levelset_sphere_data_t;

/** Distance to a sphere with given midpoint and radius.
  * data is interpreted as t8_levelset_sphere_data_t.
  * \return     dist (x,data->M) - data->radius
  */
double
t8_levelset_sphere (const t8_3D_point &x, double t, void *data);

/** Returns always 1.
 * \return 1
 */
double
t8_scalar3d_constant_one (const double x[3], double t);

/** Returns always 0.
 * \return 0
 */
double
t8_scalar3d_constant_zero (const double x[3], double t);

/** Return the x-coordinate of the input.
 * \return x[0]
 */
double
t8_scalar3d_project_x (const double x[3], double t);

/** This function is =1 if the 0.25 <= x <= 0.75 and 0 else. */
double
t8_scalar3d_step_function (const double x[3], double t);

/** This function is =1 if 0.25 <= x <= 0.75,
 * it is 0 outside of 0.25-eps and 0.75+eps,
 * it interpolates linearly in between. eps = 0.1 **/
double
t8_scalar3d_almost_step_function (const double x[3], double t);

/** A 1-d Bell-curve centered around 0.5 */
double
t8_scalar3d_exp_distribution (const double x[3], double t);

/** Sinus of 2pi x_0
 * \return sin (2pi x[0])
 */
double
t8_scalar3d_sinx (const double x[3], double t);

/** Sinus of x times cosinus of y
 * \return sin (2pi x[0]) * cos (2pi x[1])
 */
double
t8_scalar3d_sinx_cosy (const double x[3], double t);

/** Sinus of 10 * x times cosinus of y times z
 * \return 10 * sin (2pi x[0]) * cos (2pi x[1]) * x[3]
 */
double
t8_scalar3d_sinx_cosy_z (const double x[3], double t);

/** Sinus of t
 * \return sin (2pi t)
 */
double
t8_scalar3d_sint (const double x[3], double t);

/** Level-set function of a sphere around origin with radius 0.75
 * \return |x| - 0.75
 */
double
t8_scalar3d_sphere_75_radius (const t8_3D_vec x, double t);

/** Level-set function of a sphere around M = (0.5,0.5,0.5) with radius 0.375
 * \return |x - M| - 0.375
 */
double
t8_scalar3d_sphere_05_midpoint_375_radius (const t8_3D_vec x, double t);

/** Level-set function of a sphere around M = (0.3,0.3,0.3) with radius 0.25
 * \return |x - M| - 0.25
 */
double
t8_scalar3d_sphere_03_midpoint_25_radius (const t8_3D_vec x, double t);

/** Level-set function of a sphere around M = (0.5,0.5,0) with radius 0.375
 * \return |x - M| - 0.375
 */
double
t8_scalar3d_sphere_05_0z_midpoint_375_radius (const t8_3D_vec x, double t);
/** Flow functions */

/** Returns always 1 in each coordinate.
 */
void
t8_flow_constant_one_vec (const t8_3D_point &x, double t, t8_3D_vec &x_out);

/** Sets the first coordinate to 1, all other to 0. */
void
t8_flow_constant_one_x_vec (const t8_3D_point &x, double t, t8_3D_vec &x_ou);

/** Sets the first and second coordinate to 1, the third to 0. */
void
t8_flow_constant_one_xy_vec (const t8_3D_point &x, double t, t8_3D_vec &x_out);

/** Sets all coordinates to a nonzero constant. */
void
t8_flow_constant_one_xyz_vec (const t8_3D_point &x, double t, t8_3D_vec &x_out);

/** Transform the unit square to [-0.5,0.5]^2 and computes
 * x = 2pi*y, y = -2pi*x
 */
void
t8_flow_rotation_2d (const t8_3D_point &x, double t, t8_3D_vec &x_out);

void
t8_flow_compressible (const t8_3D_point &x_in, double t, t8_3D_vec &x_out);
/** Incompressible flow in unit cube */
void
t8_flow_incomp_cube_flow (const t8_3D_point &x, double t, t8_3D_vec &x_out);

/** 2d flow around a circle with radius R = 1 and
 * constant inflow with x-speed U = 1. 
 * See https://doi.org/10.13140/RG.2.2.34714.11203 */
void
t8_flow_around_circle (const t8_3D_point &x, double t, t8_3D_vec &x_out);

void
t8_flow_stokes_flow_sphere_shell (const t8_3D_point &x_in, double t, t8_3D_vec &x_out);

void
t8_flow_around_circle_with_angular_velocity (const t8_3D_point &x, double t, t8_3D_vec &x_out);

#endif /* !T8_EXAMPLE_COMMON_H */
