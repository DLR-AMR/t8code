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
  sc_array_t array;           /**< The array in which the elements are stored */
} t8_element_array_t;

T8_EXTERN_C_BEGIN ();

/** Creates a new array structure with 0 elements.
 * \param [in] scheme   The eclass scheme of which elements should be stored.
 * \return              Return an allocated array of zero length.
 */
t8_element_array_t *
t8_element_array_new (t8_eclass_scheme_c *scheme);

/** Creates a new array structure with a given length (number of elements)
 * and calls \ref t8_element_new for those elements.
 * \param [in] scheme       The eclass scheme of which elements should be stored.
 * \param [in] num_elements Initial number of array elements.
 * \return                  Return an allocated array
 *                          with allocated and initialized elements for which \ref
 *                          t8_element_new was called.
 */
t8_element_array_t *
t8_element_array_new_count (t8_eclass_scheme_c *scheme, size_t num_elements);

/** Initializes an already allocated (or static) array structure.
 * \param [in,out]  element_array  Array structure to be initialized.
 * \param [in]      scheme         The eclass scheme of which elements should be stored.
 */
void
t8_element_array_init (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme);

/** Initializes an already allocated (or static) array structure
 * and allocates a given number of elements and initializes them with \ref t8_element_init.
 * \param [in,out]  element_array Array structure to be initialized.
 * \param [in] scheme         The eclass scheme of which elements should be stored.
 * \param [in] num_elements   Number of initial array elements.
 */
void
t8_element_array_init_size (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme, size_t num_elements);

/** Initializes an already allocated (or static) view from existing t8_element_array.
 * The array view returned does not require t8_element_array_reset (doesn't hurt though).
 * \param [in,out] view  Array structure to be initialized.
 * \param [in] array     The array must not be resized while view is alive.
 * \param [in] offset    The offset of the viewed section in element units.
 *                       This offset cannot be changed until the view is reset.
 * \param [in] length    The length of the view in element units.
 *                       The view cannot be resized to exceed this length.
 *                       It is not necessary to call sc_array_reset later.
 */
void
t8_element_array_init_view (t8_element_array_t *view, t8_element_array_t *array, size_t offset, size_t length);

/** Initializes an already allocated (or static) view from given plain C data
 * (array of t8_element_t).
 * The array view returned does not require t8_element_array_reset (doesn't hurt though).
 * \param [in,out] view     Array structure to be initialized.
 * \param [in] base         The data must not be moved while view is alive.
 *                          Must be an array of t8_element_t corresponding to \a scheme.
 * \param [in] scheme       The eclass scheme of the elements stored in \a base.
 * \param [in] elem_count   The length of the view in element units.
 *                          The view cannot be resized to exceed this length.
 *                          It is not necessary to call t8_element_array_reset later.
 */
void
t8_element_array_init_data (t8_element_array_t *view, t8_element_t *base, t8_eclass_scheme_c *scheme,
                            size_t elem_count);

/** Initializes an already allocated (or static) array structure
 * and copy an existing array of t8_element_t into it.
 * \param [in,out]  element_array Array structure to be initialized.
 * \param [in] scheme         The eclass scheme of which elements should be stored.
 * \param [in] data           An array of t8_element_t which will be copied into
 *                            \a element_array. The elements in \a data must belong to
 *                            \a scheme and must be properly initialized with either
 *                            \ref t8_element_new or \ref t8_element_init.
 * \param [in] num_elements   Number of elements in \a data to be copied.
 */
void
t8_element_array_init_copy (t8_element_array_t *element_array, t8_eclass_scheme_c *scheme, t8_element_t *data,
                            size_t num_elements);

/** Change the number of elements stored in an element array.
 * \param [in,out] element_array  The element array to be modified.
 * \param [in] new_count    The new element count of the array.
 *                          If it is zero the effect equals \ref t8_element_array_reset.
 * \note If \a new_count is larger than the number of current elements on \a element_array,
 * then \ref t8_element_init is called for the new elements.
 */
void
t8_element_array_resize (t8_element_array_t *element_array, size_t new_count);

/** Copy the contents of an array into another.
 * Both arrays must have the same eclass_scheme.
 * \param [in] dest Array will be resized and get new data.
 * \param [in] src  Array used as source of new data, will not be changed.
 */
void
t8_element_array_copy (t8_element_array_t *dest, const t8_element_array_t *src);

/** Enlarge an array by one element.
 * \param [in, ou] element_array Array structure to be modified.
 * \return Returns a pointer to a newly added element for which \ref t8_element_init
 *         was called.
 */
t8_element_t *
t8_element_array_push (t8_element_array_t *element_array);

/** Enlarge an array by a number of elements.
 * \param [in, ou] element_array Array structure to be modified.
 * \param [in]     count        The number of elements to add.
 * \return Returns a pointer to the newly added elements for which \ref t8_element_init
 *                was called.
 */
t8_element_t *
t8_element_array_push_count (t8_element_array_t *element_array, size_t count);

/** Return a given element in an array. Const version.
 * \param [in]  element_array Array of elements.
 * \param [in]  index The index of an element within the array.
 * \return            A pointer to the element stored at position \a index in
 *                    \a element_array.
 */
const t8_element_t *
t8_element_array_index_locidx (const t8_element_array_t *element_array, t8_locidx_t index);

/** Return a given element in an array. Const version.
 * \param [in]  element_array Array of elements.
 * \param [in]  index The index of an element within the array.
 * \return            A pointer to the element stored at position \a index in
 *                    \a element_array.
 */
const t8_element_t *
t8_element_array_index_int (const t8_element_array_t *element_array, int index);

/** Return a given element in an array. Mutable version.
 * \param [in]  element_array Array of elements.
 * \param [in]  index The index of an element within the array.
 * \return            A pointer to the element stored at position \a index in
 *                    \a element_array.
 */
t8_element_t *
t8_element_array_index_locidx_mutable (t8_element_array_t *element_array, t8_locidx_t index);

/** Return a given element in an array. Mutable version.
 * \param [in]  element_array Array of elements.
 * \param [in]  index The index of an element within the array.
 * \return            A pointer to the element stored at position \a index in
 *                    \a element_array.
 */
t8_element_t *
t8_element_array_index_int_mutable (t8_element_array_t *element_array, int index);

/** Return the eclass scheme associated to a t8_element_array.
 * \param [in]  element_array Array of elements.
 * \return      The eclass scheme stored at \a element_array.
 */
const t8_eclass_scheme_c *
t8_element_array_get_scheme (const t8_element_array_t *element_array);

/** Return the number of elements stored in a t8_element_array_t.
 * \param [in]  element_array  Array structure.
 * \return                     The number of elements stored in \a element_array.
 */
size_t
t8_element_array_get_count (const t8_element_array_t *element_array);
/** Return the data size of elements stored in a t8_element_array_t.
 * \param [in]  element_array  Array structure.
 * \return                     The size (in bytes) of a single element in \a element_array.
 */
size_t
t8_element_array_get_size (const t8_element_array_t *element_array);

/** Return a const pointer to the real data array stored in a t8_element_array.
 * \param [in]  element_array  Array structure.
 * \return                     A pointer to the stored data. If the number of stored
 *                             elements is 0, then NULL is returned.
 */
const t8_element_t *
t8_element_array_get_data (const t8_element_array_t *element_array);

/** Return a pointer to the real data array stored in a t8_element_array.
 * \param [in]  element_array  Array structure.
 * \return                     A pointer to the stored data. If the number of stored
 *                             elements is 0, then NULL is returned.
 */
t8_element_t *
t8_element_array_get_data_mutable (t8_element_array_t *element_array);

/** Return a const pointer to the sc_array stored in a t8_element_array.
 * \param [in]  element_array  Array structure.
 * \return                     A const pointer to the sc_array storing the data.
 * \note The data cannot be modified.
 */
const sc_array_t *
t8_element_array_get_array (const t8_element_array_t *element_array);

/** Return a mutable pointer to the sc_array stored in a t8_element_array.
 * \param [in]  element_array  Array structure.
 * \return                     A pointer to the sc_array storing the data.
 * \note The data can be modified.
 */
sc_array_t *
t8_element_array_get_array_mutable (t8_element_array_t *element_array);

/** Sets the array count to zero and frees all elements.
 * \param [in,out]  element_array  Array structure to be reset.
 * \note Calling t8_element_array_init, then any array operations,
 *       then t8_element_array_reset is memory neutral.
 */
void
t8_element_array_reset (t8_element_array_t *element_array);

/** Sets the array count to zero, but does not free elements.
 * \param [in,out]  element_array  Element array structure to be truncated.
 * \note This is intended to allow an t8_element_array to be used as a reusable
 * buffer, where the "high water mark" of the buffer is preserved, so that
 * O(log (max n)) reallocs occur over the life of the buffer.
 */
void
t8_element_array_truncate (t8_element_array_t *element_array);

T8_EXTERN_C_END ();

#endif /* !T8_CONTAINERS_HXX */
