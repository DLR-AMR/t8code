/*  This file is part of t8code.
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

#include <t8_element.hxx>
#include <sc_containers.h>
#include <t8_data/t8_containers.h>

T8_EXTERN_C_BEGIN ();

#ifdef T8_ENABLE_DEBUG
/* Query whether an element array is initialized properly. */
static int
t8_element_array_is_valid (const t8_element_array_t *element_array)
{
  int is_valid;

  /* Check that all pointers are not NULL */
  is_valid = element_array != NULL && element_array->scheme != NULL;

  if (!is_valid) {
    return 0;
  }

  /* Check that the element size of the scheme matches the size of data elements
   * stored in the array. */
  is_valid = is_valid && element_array->scheme->t8_element_size () == element_array->array.elem_size;

  return is_valid;
}
#endif

t8_element_array_t *
t8_element_array_new (t8_eclass_scheme_c *scheme)
{
  t8_element_array_t *new_array;

  /* allocate memory */
  new_array = T8_ALLOC (t8_element_array_t, 1);
  /* initialize array */
  t8_element_array_init (new_array, scheme);
  T8_ASSERT (t8_element_array_is_valid (new_array));

  return new_array;
}

t8_element_array_t *
t8_element_array_new_count (t8_eclass_scheme_c *scheme, size_t num_elements)
{
  t8_element_array_t *new_array;

  /* allocate memory */
  new_array = T8_ALLOC (t8_element_array_t, 1);
  /* initialize array */
  t8_element_array_init_size (new_array, scheme, num_elements);
  T8_ASSERT (t8_element_array_is_valid (new_array));

  return new_array;
}

void
t8_element_array_init (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme)
{
  size_t elem_size;

  T8_ASSERT (element_array != NULL);

  /* set the scheme */
  element_array->scheme = scheme;
  /* get the size of an element and initialize the array member */
  elem_size = scheme->t8_element_size ();
  sc_array_init (&element_array->array, elem_size);
  T8_ASSERT (t8_element_array_is_valid (element_array));
}

void
t8_element_array_init_size (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme, size_t num_elements)
{
  t8_element_t *first_element;
  T8_ASSERT (element_array != NULL);

  element_array->scheme = scheme;
  /* allocate the elements */
  sc_array_init_size (&element_array->array, scheme->t8_element_size (), num_elements);

  if (num_elements > 0) {
    /* Call t8_element_init for the elements */
    first_element = (t8_element_t *) sc_array_index (&element_array->array, 0);
    scheme->t8_element_init (num_elements, first_element);
  }
  T8_ASSERT (t8_element_array_is_valid (element_array));
}

void
t8_element_array_init_view (t8_element_array_t *view, const t8_element_array_t *array, size_t offset, size_t length)
{
  T8_ASSERT (t8_element_array_is_valid (array));

  /* Initialize the element array.
   * Unfortunately, we have to cast away the constness to pass to sc_array_init_view.
   */
  sc_array_init_view (&view->array, &((t8_element_array_t *) array)->array, offset, length);
  /* Set the scheme */
  view->scheme = array->scheme;
  T8_ASSERT (t8_element_array_is_valid (view));
}

void
t8_element_array_init_data (t8_element_array_t *view, t8_element_t *base, t8_eclass_scheme_c *scheme, size_t elem_count)
{
  /* Initialize the element array */
  sc_array_init_data (&view->array, (void *) base, scheme->t8_element_size (), elem_count);
  /* set the scheme */
  view->scheme = scheme;
  T8_ASSERT (t8_element_array_is_valid (view));
}

void
t8_element_array_init_copy (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme, t8_element_t *data,
                            size_t num_elements)
{
  sc_array_t *array;
  T8_ASSERT (element_array != NULL);

  t8_element_array_init (element_array, scheme);

  array = &element_array->array;
#ifdef T8_ENABLE_DEBUG
  /* Check if the elements in data are valid for scheme */
  {
    size_t ielem;
    const t8_element_t *element;
    size_t size;

    size = scheme->t8_element_size ();
    for (ielem = 0; ielem < num_elements; ielem++) {
      /* data is of incomplete type, we thus have to manually set the address
       * of the ielem-th t8_element */
      element = (const t8_element_t *) (((char *) data) + ielem * size);
      T8_ASSERT (scheme->t8_element_is_valid (element));
    }
  }
#endif
  /* Allocate enough memory for the new elements */
  sc_array_init_size (array, scheme->t8_element_size (), num_elements);
  /* Copy the elements in data */
  memcpy (array->array, data, num_elements * array->elem_size);
}

void
t8_element_array_resize (t8_element_array_t *element_array, size_t new_count)
{
  size_t old_count;
  T8_ASSERT (t8_element_array_is_valid (element_array));
  /* Store the old number of elements */
  old_count = t8_element_array_get_count (element_array);
  if (old_count < new_count) {
    /* if the new_count is larger than the previous count, we need to
    * call t8_element_init on the newly allocated elements. */
    sc_array_resize (&element_array->array, new_count);
    t8_element_t *first_new_elem;
    /* Get the first newly allocated element */
    first_new_elem = t8_element_array_index_locidx_mutable (element_array, old_count);
    /* Call t8_element_init on all new elements */
    element_array->scheme->t8_element_init (new_count - old_count, first_new_elem);
  }
  else if (old_count > new_count) {
    t8_element_t *first_old_elem;
    /* Get the first element to deinit */
    first_old_elem = t8_element_array_index_locidx_mutable (element_array, new_count);
    element_array->scheme->t8_element_deinit (old_count - new_count, first_old_elem);
    sc_array_resize (&element_array->array, new_count);
  }
  else {
    T8_ASSERT (new_count == element_array->array.elem_count);
    /* Free the allocated, but unused memory. */
    sc_array_resize (&element_array->array, new_count);
  }
}

void
t8_element_array_copy (t8_element_array_t *dest, const t8_element_array_t *src)
{
  T8_ASSERT (t8_element_array_is_valid (dest));
  T8_ASSERT (t8_element_array_is_valid (src));
  T8_ASSERT (dest->scheme == src->scheme);
  sc_array_copy (&dest->array, (sc_array_t *) &src->array); /* need to convert src->array to non-const */
}

t8_element_t *
t8_element_array_push (t8_element_array_t *element_array)
{
  t8_element_t *new_element;
  T8_ASSERT (t8_element_array_is_valid (element_array));
  new_element = (t8_element_t *) sc_array_push (&element_array->array);
  element_array->scheme->t8_element_init (1, new_element);
  return new_element;
}

t8_element_t *
t8_element_array_push_count (t8_element_array_t *element_array, size_t count)
{
  t8_element_t *new_elements;
  T8_ASSERT (t8_element_array_is_valid (element_array));
  /* grow the array */
  new_elements = (t8_element_t *) sc_array_push_count (&element_array->array, count);
  /* initialize the elements */
  element_array->scheme->t8_element_init (count, new_elements);
  return new_elements;
}

const t8_element_t *
t8_element_array_index_locidx (const t8_element_array_t *element_array, t8_locidx_t index)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  return (const t8_element_t *) t8_sc_array_index_locidx (&element_array->array, index);
}

const t8_element_t *
t8_element_array_index_int (const t8_element_array_t *element_array, int index)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  return (const t8_element_t *) sc_array_index_int ((sc_array_t *) &element_array->array,
                                                    index); /* Need to convert element_array->array to non-const */
}

t8_element_t *
t8_element_array_index_locidx_mutable (t8_element_array_t *element_array, t8_locidx_t index)
{
  return (t8_element_t *) t8_element_array_index_locidx (element_array, index);
}

t8_element_t *
t8_element_array_index_int_mutable (t8_element_array_t *element_array, int index)
{
  return (t8_element_t *) t8_element_array_index_int (element_array, index);
}

const t8_eclass_scheme_c *
t8_element_array_get_scheme (const t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  return element_array->scheme;
}

size_t
t8_element_array_get_count (const t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  return element_array->array.elem_count;
}

size_t
t8_element_array_get_size (const t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  return element_array->scheme->t8_element_size ();
}

const t8_element_t *
t8_element_array_get_data (const t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));

  if (element_array->array.elem_count > 0) {
    return (t8_element_t *) t8_element_array_index_locidx (element_array, 0);
  }
  else {
    return NULL;
  }
}

t8_element_t *
t8_element_array_get_data_mutable (t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));

  return (t8_element_t *) t8_element_array_get_data (element_array);
}

const sc_array_t *
t8_element_array_get_array (const t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));

  return &element_array->array;
}

sc_array_t *
t8_element_array_get_array_mutable (t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));

  return &element_array->array;
}

void
t8_element_array_reset (t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  size_t count = t8_element_array_get_count (element_array);
  if (count > 0) {
    t8_element_t *first_elem = t8_element_array_index_locidx_mutable (element_array, 0);
    element_array->scheme->t8_element_deinit (count, first_elem);
  }
  sc_array_reset (&element_array->array);
}

void
t8_element_array_truncate (t8_element_array_t *element_array)
{
  T8_ASSERT (t8_element_array_is_valid (element_array));
  size_t count = t8_element_array_get_count (element_array);
  if (count > 0) {
    t8_element_t *first_elem = t8_element_array_index_locidx_mutable (element_array, 0);
    element_array->scheme->t8_element_deinit (count, first_elem);
  }
  sc_array_truncate (&element_array->array);
}

T8_EXTERN_C_END ();
