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

/** \file t8_containers.h
 * We define the t8_element_array that stores elements of a given
 * eclass scheme.
 */

#ifndef T8_CONTAINERS_H
#define T8_CONTAINERS_H

#include <t8.h>
#include <t8_element.h>
#include <sc_containers.h>

/** The t8_element_array_t is an array to store t8_element_t * of a given
 * eclass_scheme implementation. It is a wrapper around \ref sc_array_t.
 * Each time, a new element is created by the functions for \ref t8_element_array_t,
 * the eclass function either \ref t8_element_new or \ref t8_element_init is called
 * for the element.
 * Thus, each element in a \ref t8_element_array_t is automatically initialized properly.
 */
typedef struct
{
  t8_eclass_scheme_c *scheme; /**< An eclass scheme of which elements should be stored */
  sc_array_t          array;  /**< The array in which the elements are stored */
} t8_element_array_t;

T8_EXTERN_C_BEGIN ();

/** Creates a new array structure with 0 elements.
 * \param [in] scheme   The eclass scheme of which elements should be stored.
 * \return              Return an allocated array of zero length.
 */
t8_element_array_t *t8_element_array_new (t8_eclass_scheme_c * scheme);

/** Creates a new array structure with a given length (number of elements)
 * and calls \ref t8_element_new for those elements.
 * \param [in] scheme       The eclass scheme of which elements should be stored.
 * \param [in] num_elements Initial number of array elements.
 * \return                  Return an allocated array
 *                          with allocated and initialized elements for which \ref
 *                          t8_element_new was called.
 */
t8_element_array_t *t8_element_array_new_count (t8_eclass_scheme_c * scheme,
                                                size_t num_elements);

/** Initializes an already allocated (or static) array structure.
 * \param [in,out]  element_array  Array structure to be initialized.
 * \param [in]      scheme         The eclass scheme of which elements should be stored.
 */
void                t8_element_array_init (t8_element_array_t * element_array,
                                           t8_eclass_scheme_c * scheme);

/** Initializes an already allocated (or static) array structure
 * and allocates a given number of elements with \ref t8_element_new.
 * \param [in,out]  element_array Array structure to be initialized.
 * \param [in] scheme         The eclass scheme of which elements should be stored.
 * \param [in] num_elements   Number of initial array elements.
 */
void                t8_element_array_init_size (t8_element_array_t *
                                                element_array,
                                                t8_eclass_scheme_c * scheme,
                                                size_t num_elements);

/** Enlarge an array by one element.
 * \param [in, ou] element_array Array structure to be modified.
 * \return Returns a pointer to a newly added element for which \ref t8_element_init
 *         was called.
 */
t8_element_t       *t8_element_array_push (t8_element_array_t *
                                           element_array);

/** Return the number of elements stored in a t8_element_array_t.
 * \param [in]  element_array  Array structure.
 * \return                     The number of elements stored in \a element_array.
 */
size_t              t8_element_array_get_count (t8_element_array_t *
                                                element_array);

T8_EXTERN_C_END ();

#endif /* !T8_CONTAINERS_HXX */
