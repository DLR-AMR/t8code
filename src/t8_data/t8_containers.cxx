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

/** \file t8_containers.cxx
 *
 * TODO: document this file
 */

#include <t8_element_cxx.hxx>
#include <t8_data/t8_containers.h>

T8_EXTERN_C_BEGIN ();

t8_element_array_t *
t8_element_array_new (t8_eclass_scheme_c * scheme)
{
  t8_element_array_t *new_array;

  /* allocate memory */
  new_array = T8_ALLOC (t8_element_array_t, 1);
  /* initialize array */
  t8_element_array_init (new_array, scheme);

  return new_array;
}

t8_element_array_t *
t8_element_array_new_count (t8_eclass_scheme_c * scheme, size_t num_elements)
{
  t8_element_array_t *new_array;

  /* allocate memory */
  new_array = T8_ALLOC (t8_element_array_t, 1);
  /* initialize array */
  t8_element_array_init_size (new_array, scheme, num_elements);

  return new_array;
}

void
t8_element_array_init (t8_element_array_t * element_array,
                       t8_eclass_scheme_c * scheme)
{
  size_t              elem_size;

  T8_ASSERT (element_array != NULL);

  /* set the scheme */
  element_array->scheme = scheme;
  /* get the size of an element and initialize the array member */
  elem_size = scheme->t8_element_size ();
  sc_array_init (&element_array->array, elem_size);
}

void
t8_element_array_init_size (t8_element_array_t * element_array,
                            t8_eclass_scheme_c * scheme, size_t num_elements)
{
  sc_array_t         *data;
  T8_ASSERT (element_array != NULL);

  t8_element_array_init (element_array, scheme);

  /* allocate the elements by calling t8_element_new */
  data = &element_array->array;
  scheme->t8_element_new (num_elements, (t8_element_t **) & data->array);
  /* Set the elem_count and byte_alloc fields of the sc_array by hand */
  data->elem_count = num_elements;
  data->byte_alloc = num_elements * scheme->t8_element_size ();
}

t8_element_t       *
t8_element_array_push (t8_element_array_t * element_array)
{
  t8_element_t       *new_element;

  new_element = (t8_element_t *) sc_array_push (&element_array->array);
  element_array->scheme->t8_element_init (1, new_element, 0);

  return new_element;
}

size_t
t8_element_array_get_count (t8_element_array_t * element_array)
{
  T8_ASSERT (element_array != NULL);

  return element_array->array.elem_count;
}

T8_EXTERN_C_END ();
