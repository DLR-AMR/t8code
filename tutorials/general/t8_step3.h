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

/** file t8_step3.h
 * This is the header file to the step3 example of t8code. It collects
 * functions of t8_step3 that we reuse in other examples.
 * In this example we discuss how to adapt a forest.
 * The main program is t8_step3_main.
 * See \ref t8_step3_adapt_forest.cxx for more details.
 */

#ifndef T8_STEP3_H
#define T8_STEP3_H

#include <t8.h>                          /* General t8code header, always include this. */
#include <t8_forest/t8_forest_general.h> /* forest definition and basic interface. */

T8_EXTERN_C_BEGIN ();

/** This is the main program of this example. It creates a coarse mesh and a forest,
 *  adapts the forest and writes some output.
 */
int
t8_step3_main (int argc, char **argv);

/* Functions used for other examples are below. */

/** Print the local and global number of elements of a forest. */
void
t8_step3_print_forest_information (t8_forest_t forest);

/* This is our own defined data that we will pass on to the
 * adaptation callback. */
struct t8_step3_adapt_data
{
  double midpoint[3];               /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
};

/** Adapt a forest according to our t8_step3_adapt_callback function.
 *  Thus, the input forest will get refined inside a sphere 
 *  of radius 0.2 around (0.5, 0.5, 0.5) and coarsened outside of radius 0.4.
 * \param [in] forest   A committed forest.
 * \return              A new forest that arises from the input \a forest via adaptation.
 */
t8_forest_t
t8_step3_adapt_forest (t8_forest_t forest);

/* The adaptation callback function. This will refine elements inside of a given sphere
 * and coarsen the elements outside of a given sphere.
 * The necessary input data is of type t8_step3_adapt_data and should be passed to forest
 * via t8_forest_set_user_data before calling this function.
 * See t8_step3.cxx for more details.
 */
int
t8_step3_adapt_callback (t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                         t8_locidx_t lelement_id, const t8_scheme_c *scheme, const int is_family,
                         const int num_elements, t8_element_t *elements[], [[maybe_unused]] void *user_data,
                         [[maybe_unused]] void *t8code_data);

T8_EXTERN_C_END ();

#endif /* !T8_STEP3_H */
