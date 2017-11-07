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

/** file t8_vtk.h
 * This header file collects macros that are needed for
 * the forest and cmesh vtk routines.
 * \see t8_forest_vtk.h \see t8_cmesh_vtk.h
 */

#ifndef T8_EXAMPLE_COMMON_H
#define T8_EXAMPLE_COMMON_H

#include <t8.h>
#include <t8_forest.h>

/** A levelset function in 3 space dimensions. */
typedef double      (*t8_example_level_set_fn) (double, double, double,
                                                void *);

/** Struct to handle refinement around a level-set function. */
typedef struct
{
  t8_example_level_set_fn L;  /**< The level set function. */
  void               *udata; /**< Data pointer that is passed to L */
  double              band_width; /**< Width of max_level elements around the zero-level set */
  int                 min_level; /**< The minimal refinement level. Elements with this level will not be coarsened. */
  int                 max_level; /**< The maximum refinement level. Elements with this level will not be refined. */
} t8_example_level_set_struct_t;

/** Function pointer for real valued functions from d+1 space dimensions
 * functions f: R^d x R -> R */
typedef double      (*t8_scalar_function_1d_fn) (double x, double t);
typedef double      (*t8_scalar_function_2d_fn) (const double x[2], double t);
typedef double      (*t8_scalar_function_3d_fn) (const double x[3], double t);
/** Function pointer for a vector valued function
 *  f: R^3 x R -> R */
typedef void        (*t8_flow_function_3d_fn) (const double x_in[3], double t,
                                               double x_out[3]);

T8_EXTERN_C_BEGIN ();

/* function declarations */

/** Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 */
int                 t8_common_adapt_balance (t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             t8_locidx_t which_tree,
                                             t8_locidx_t lelement_id,
                                             t8_eclass_scheme_c * ts,
                                             int num_elements,
                                             t8_element_t * elements[]);

/** Adapt a forest along a given level-set function.
 * The user data of forest must be a pointer to a \a t8_example_level_set_struct_t.
 * An element in the forest is refined, if it is in a band of \a band_with many
 * \a max_level elements around the zero level-set Gamma = { x | L(x) = 0}
 */
/* TODO: Currently the band_width control is not working yet.
 *        if band_with = 0, then all elements that are touched by the zero LS are refined. */
int                 t8_common_adapt_level_set (t8_forest_t forest,
                                               t8_forest_t forest_from,
                                               t8_locidx_t which_tree,
                                               t8_locidx_t lelement_id,
                                               t8_eclass_scheme_c * ts,
                                               int num_elements,
                                               t8_element_t * elements[]);

/** Compute the coordinates of the midpoint of an element.
 * \param [in]  forest  The forest in which the element is in (must be committed).
 * \param [in]  which_tree The local tree id of tree in which the element is in.
 * \param [in]  ts      The eclass scheme associated to the element.
 * \param [in]  element The element.
 * \param [in,out] elem_midpoint_f An array of 3 doubles. On output the coordinates
 *              of the midpoint of \a element are stored.
 * \note \a forest must be committed before calling this function.
 */
void                t8_common_midpoint (t8_forest_t forest,
                                        t8_locidx_t which_tree,
                                        t8_eclass_scheme_c * ts,
                                        t8_element_t * element,
                                        double elem_midpoint_f[3]);

/** Real valued functions defined in t8_example_common_functions.h */

/** Returns always 1.
 * \return 1
 */
double              t8_constant_one (const double x[3], double t);

/** Returns always 0.
 * \return 0
 */
double              t8_constant_zero (const double x[3], double t);

/** Return the x-coordinate of the input.
 * \return x[0]
 */
double              t8_project_x (const double x[3], double t);

/** This function is =1 if the 0.25 <= x <= 0.75 and 0 else. */
double              t8_step_function (const double x[3], double t);

/** This function is =1 if 0.25 <= x <= 0.75,
 * it is 0 outside of 0.25-eps and 0.75+eps,
 * it interpolates linearly in between. eps = 0.1 **/
double              t8_almost_step_function (const double x[3], double t);

/** A 1-d Bell-curve centered around 0.5 */
double              t8_exp_distribution (const double x[3], double t);

/** Sinus of 2pi x_0
 * \return sin (2pi x[0])
 */
double              t8_sinx (const double x[3], double t);

/** Sinus of x times cosinus of y
 * \return sin (2pi x[0]) * cos (2pi x[1])
 */
double              t8_sinx_cosy (const double x[3], double t);

/** Sinus of 10 * x times cosinus of y times z
 * \return 10 * sin (2pi x[0]) * cos (2pi x[1]) * x[3]
 */
double              t8_sinx_cosy_z (const double x[3], double t);

/** Sinus of t
 * \return sin (2pi t)
 */
double              t8_sint (const double x[3], double t);

/** Level-set function of a sphere around origin with radius 0.75
 * \return |x| - 0.75
 */
double              t8_sphere_75_radius (const double x[3], double t);

/** Level-set function of a sphere around M = (0.5,0.5,0.5) with radius 0.375
 * \return |x - M| - 0.375
 */
double              t8_sphere_05_midpoint_375_radius (const double x[3],
                                                      double t);

/** Level-set function of a sphere around M = (0.5,0.5,0) with radius 0.375
 * \return |x - M| - 0.375
 */
double              t8_sphere_05_0z_midpoint_375_radius (const double x[3],
                                                         double t);
/** Flow functions */

/** Returns always 1 in each coordinate.
 */
void                t8_constant_one_vec (const double x[3], double t,
                                         double x_out[3]);

/** Sets the first coordinate to 1, all other to 0. */
void                t8_constant_one_x_vec (const double x[3], double t,
                                           double x_out[3]);

/** Sets the first and second coordinate to 1, the third to 0. */
void                t8_constant_one_xy_vec (const double x[3], double t,
                                            double x_out[3]);
/** Transform the unit square to [-0.5,0.5]^2 and computes
 * x = y, y = -x
 */
void                t8_rotation_2d (const double x[3], double t,
                                    double x_out[3]);

void                t8_compressible (const double x_in[3], double t,
                                     double x_out[3]);
/** Incompressible flow in unit cube */
void                t8_incomp_cube_flow (const double x[3], double t,
                                         double x_out[3]);

void                t8_stokes_flow_sphere_shell (const double x[3], double t,
                                                 double x_out[3]);

T8_EXTERN_C_END ();

#endif /* !T8_EXAMPLE_COMMON_H */
