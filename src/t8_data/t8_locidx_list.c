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

/** \file t8_locidx_list.c
 *
 * In this file we implement the routines to handle
 * t8_locidx_list_t access.
 */

#include <t8_data/t8_locidx_list.h>

/* Allocate a new t8_locidx_list object and return a poitner to it.
 * \return A pointer to a newly allocated t8_locidx_list_t
 */
t8_locidx_list_t   *
t8_locidx_list_new ()
{
  t8_locidx_list_t   *list;

  /* Allocate memory for our list */
  list = T8_ALLOC (t8_locidx_list_t, 1);

  /* initialize the new list */
  t8_locidx_list_init (list);

  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));
  /* return the initialized list */
  return list;
}

/* Initialize an allocated t8_locidx_list
 * \param [in, out] Pointer to an allocated t8_locidx_list_t
 * \note If this list is already filled the all its elements are free'd
 */
void
t8_locidx_list_init (t8_locidx_list_t * list)
{
  /* Allocate a new mempool for the list links */
  sc_mempool_t       *link_mempool = sc_mempool_new (sizeof (sc_link_t));

  /* Check that list is not NULL */
  T8_ASSERT (list != NULL);

  /* Initialize the data list */
  sc_list_init (&list->list, link_mempool);

  /* Initialize the mempool to store the t8_locidx_t */
  sc_mempool_init (&list->indices, sizeof (t8_locidx_t));

  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));
}

/* Query whether a t8_locidx_list_t is initialized.
 * \param [in] list A pointer to a list.
 * \return True (non-zero) if and only if the list is properly initialized.
 * \note Calling this with \a list = NULL is valid and will return 0.
 */
int
t8_locidx_list_is_initialized (const t8_locidx_list_t * list)
{
  if (list == NULL) {
    return 0;
  }

  /* A list is initialized if 
   * - its mempool has the proper size associated to it. 
   * - its list stores sc_link_t objects
   */
  if (list->indices.elem_size != sizeof (t8_locidx_t)
      || list->list.allocator == NULL
      || list->list.allocator->elem_size != sizeof (sc_link_t)) {
    return 0;
  }

  /* Safety check that the number of elements in the list matches
   * the number of elements in the mempool. */
  T8_ASSERT (list->list.elem_count == list->indices.elem_count);
  return 1;
}

/* Remove the first element in a list and return its entry.
 * \param [in,out] list An initliazed list with at least one entry.
 *                      On return its first entry will have been removed.
 * \return              The first entry of \a list.
 * \note It is illegal to call this function if \a list does not have any elements.
 */
t8_locidx_t
t8_locidx_list_pop (t8_locidx_list_t * list)
{
  /* List must be initialized and non-empty */
  T8_ASSERT (t8_locidx_list_is_initialized (list));
  T8_ASSERT (t8_locidx_list_count (list) > 0);
  /* Pop the last item from the internal list */
  t8_locidx_t        *ppopped_index =
    (t8_locidx_t *) sc_list_pop (&list->list);
  t8_locidx_t         popped_index = *ppopped_index;

  /* Free the memory associated to the popped index */
  sc_mempool_free (&list->indices, ppopped_index);

  /* return the index */
  return popped_index;
}

/* Remove an item after a given list position.
 * \param [in,out] list An initliazed nonempty list.
 * \param [in]     link The predecessor of the element to be removed.
 *                      If \a pred == NULL, the first element is removed,
 *                      which is equivalent to calling \ref t8_locidx_list_pop.
 * \return              The value stored at \a link.
 */
t8_locidx_t
t8_locidx_list_remove (t8_locidx_list_t * list, sc_link_t * pred)
{
  T8_ASSERT (t8_locidx_list_is_initialized (list));

  /* remove the data and store a pointer to it */
  t8_locidx_t        *pdata =
    (t8_locidx_t *) sc_list_remove (&list->list, pred);
  /* Store the value that the data points to */
  t8_locidx_t         data = *pdata;
  /* Free the data pointer */
  sc_mempool_free (&list->indices, pdata);

  return data;
}

/* Returns the number of entries in a list.
 * \param [in] list An initialized list.
 * \return          The number of entries in \a list.
 */
size_t
t8_locidx_list_count (const t8_locidx_list_t * list)
{
  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));

  return list->list.elem_count;
}

/* Append a new entry to the end of a list.
 * \param [in,out] list An initialized list.
 * \param [in]     entry The new entry.
 */
void
t8_locidx_list_append (t8_locidx_list_t * list, const t8_locidx_t entry)
{
  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));

  /* Allocated an element in the mempool. */
  t8_locidx_t        *allocated_entry =
    (t8_locidx_t *) sc_mempool_alloc (&list->indices);
  /* Add the data to the allocated element. */
  *allocated_entry = entry;
  /* Append this element to the list */
  sc_list_append (&list->list, (void *) allocated_entry);
}

/* Free all elements of a list.
 * \param [in, out] list An initialized list.
 * After calling this function \a list will not be initialized and
 * all its elements will have been free'd and its count will be zero.
 */
void
t8_locidx_list_reset (t8_locidx_list_t * list)
{
  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));

  sc_list_reset (&list->list);
  sc_mempool_reset (&list->indices);
  /* Destroy the link mempool */
  sc_mempool_destroy (list->list.allocator);
}

/** Free all elements of a list and free the list pointer.
 * \param [in, out] list An initialized list.
 * After calling this function \a list will be the NULL pointer
 * and all the elements will be freed.
 */
void
t8_locidx_list_destroy (t8_locidx_list_t ** plist)
{
  T8_ASSERT (plist != NULL);
  t8_locidx_list_t   *list = *plist;
  /* Assert that this list is initialized */
  T8_ASSERT (t8_locidx_list_is_initialized (list));

  /* Reset the list */
  t8_locidx_list_reset (list);

  /* Free the pointer */
  T8_FREE (list);
  /* Set the pointer to NULL */
  *plist = NULL;
}

#ifdef T8_ENABLE_DEBUG
/** Convert a list to a string for (debug) output.
 * \param [in]      list An initialized list.
 * \return               An allocated string that shows the contents of the list.
 * \note The user is responsible for deallocating the memory of the string with
 * T8_FREE.
 * \note The string will look like this: "NULL -> Entry0 -> Entry1 ... -> EntryN -> NULL"
 * If it is truncated then the list was too large to fit into the string. */
char               *
t8_locidx_list_to_string (t8_locidx_list_t * list)
{
  /* We will build a string in the following format:
   * "NULL -> Entry0 -> Entry1 ... -> EntryN -> NULL"
   */

  /* Calculate the length of the string */
  const size_t        num_entries = t8_locidx_list_count (list);
  const size_t        string_len = 12   /* "NULL ->" + " NULL" */
    + num_entries * 20;         /* Reserve 20 chars per entry */
  char               *string = T8_ALLOC (char, string_len);

  size_t              bytes_written = 0;
  snprintf (string, string_len + 1, "NULL ->");

  /* Iterate through the list and write the contents to the string */
  t8_locidx_list_iterator_t it;
  for (t8_locidx_list_iterator_init (list, &it);
       t8_locidx_list_iterator_is_valid (&it, list);
       t8_locidx_list_iterator_next (&it)) {
    t8_locidx_t         entry = t8_locidx_list_iterator_get_value (&it);
    bytes_written = strlen (string);
    snprintf (string + bytes_written, string_len - bytes_written, " %i ->",
              entry);
  }
  bytes_written = strlen (string);
  snprintf (string + bytes_written, string_len - bytes_written, " NULL");
  bytes_written = strlen (string);
  /* Set terminating zero */
  string[bytes_written] = '\0';

  return string;
}
#endif /* T8_ENABLE_DEBUG */

/* Initialize an allocated iterator for a list.
 *  The iterator will start at the first item in the list.
 * \param [in]  An initialized list.
 */
void
t8_locidx_list_iterator_init (t8_locidx_list_t * list,
                              t8_locidx_list_iterator_t * it)
{
  T8_ASSERT (t8_locidx_list_is_initialized (list));
  T8_ASSERT (it != NULL);

  /* Set the list pointer */
  it->list = list;
  /* Set the current entry to the first list entry */
  it->current = list->list.first;
  /* Set the prev pointer to NULL */
  it->prev = NULL;
}

/* Check whether an iterator ist valid and associated to a given list.
 * \param [in] iterator The iterator to check.
 * \param [in] list     The list to which to \a iterator is expected to be associated to.
 * \return              True (non-zero) if \a iterator is valid and points to an element in \a list.
 *                      False if either \a iterator points to the end of \a list (all items have been
 *                      iterated through.) or \a iterator is not initialized.
 */
int
t8_locidx_list_iterator_is_valid (const t8_locidx_list_iterator_t * it,
                                  const t8_locidx_list_t * list)
{
  /* Can't be valid if we are NULL */
  if (it == NULL || list == NULL) {
    return 0;
  }
  /* Cannot be valid if the list is not initialized */
  if (!t8_locidx_list_is_initialized (list)) {
    return 0;
  }

  /* The list should be the list of the iterator */
  if (it->list != list) {
    return 0;
  }

  /* If we are at the end of the list, return 0 */
  if (it->current == NULL) {
    return 0;
  }

  if (it->prev != NULL) {
    /* The next pointer of prev must be current */
    if (it->prev->next != it->current) {
      return 0;
    }
  }

  return 1;
}

/* Let an iterator point to the next entry of its list.
 * \param [in, out] it  The iterator.
 */
void
t8_locidx_list_iterator_next (t8_locidx_list_iterator_t * it)
{
  T8_ASSERT (it != NULL);
  T8_ASSERT (t8_locidx_list_iterator_is_valid (it, it->list));

  /* Store the current pointer */
  sc_link_t          *temp = it->current;

  /* Advance the current pointer */
  it->current = it->current->next;
  /* Reset the prev pointer */
  it->prev = temp;
}

/* Return the value of the item that an iterator currently points to.
 * \param [in] it   The iterator.
 * \return          The value of the t8_locidx_t in the iterator's associated
 *                  list that it currently points to.
 */
t8_locidx_t
t8_locidx_list_iterator_get_value (const t8_locidx_list_iterator_t * it)
{
  T8_ASSERT (it != NULL);
  T8_ASSERT (t8_locidx_list_iterator_is_valid (it, it->list));
  T8_ASSERT (it->current->data != NULL);

  return *(t8_locidx_t *) it->current->data;
}

/* Remove the entry that an iterator points to from the list.
 * \param [in] it   The iterator.
 */
void
t8_locidx_list_iterator_remove_entry (t8_locidx_list_iterator_t * it)
{
  T8_ASSERT (it != NULL);
  T8_ASSERT (t8_locidx_list_iterator_is_valid (it, it->list));

  /* Store the link after current */
  sc_link_t          *next = it->current->next;

  /* Remove the current link from the list */
  t8_locidx_list_remove (it->list, it->prev);

  /* Set the new current pointer */
  it->current = next;
}
