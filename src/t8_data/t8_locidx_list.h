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

/** \file t8_locidx_list.h
 * We define a linked list of t8_locidx_t as a wrapper around sc_list_t.
 */

#ifndef T8_LOCIDX_LIST_H
#define T8_LOCIDX_LIST_H

#include <t8.h>

typedef struct t8_shmem_array *t8_shmem_array_t;

T8_EXTERN_C_BEGIN ();

/** Defines a linked list of t8_locidx_t entries.
 * The entries are allocated and stored in a mempool.
 */
typedef struct
{
  sc_list_t           list;
                  /**< The linked list */
  sc_mempool_t        indices;
                        /**< The mempool that stores the indices */
} t8_locidx_list_t;

/** Allocate a new t8_locidx_list object and return a poitner to it.
 * \return A pointer to a newly allocated t8_locidx_list_t
 */
t8_locidx_list_t   *t8_locidx_list_new ();

/** Initialize an allocated t8_locidx_list
 * \param [in, out] Pointer to an allocated t8_locidx_list_t
 * \note If this list is already filled the all its elements are free'd
 */
void                t8_locidx_list_init (t8_locidx_list_t * list);

/** Query whether a t8_locidx_list_t is initialized.
 * \param [in] list A pointer to a list.
 * \return True (non-zero) if and only if the list is properly initialized.
 * \note Calling this with \a list = NULL is valid and will return 0.
 */
int                 t8_locidx_list_is_initialized (const t8_locidx_list_t *
                                                   list);

/** Remove the first element in a list and return its entry.
 * \param [in,out] list An initliazed list with at least one entry.
 *                      On return its first entry will have been removed.
 * \return              The first entry of \a list.
 * \note It is illegal to call this function if \a list does not have any elements.
 */
t8_locidx_t         t8_locidx_list_pop (t8_locidx_list_t * list);

/** Remove an item after a given list position.
 * \param [in,out] list An initliazed nonempty list.
 * \param [in]     link The predecessor of the element to be removed.
 *                      If \a pred == NULL, the first element is removed,
 *                      which is equivalent to calling \ref t8_locidx_list_pop.
 * \return              The value stored at \a link.
 */
t8_locidx_t         t8_locidx_list_remove (t8_locidx_list_t * list,
                                           sc_link_t * pred);

/** Returns the number of entries in a list.
 * \param [in] list An initialized list.
 * \return          The number of entries in \a list.
 */
size_t              t8_locidx_list_count (const t8_locidx_list_t * list);

/** Append a new entry to the end of a list.
 * \param [in,out] list An initialized list.
 * \param [in]     entry The new entry.
 */
void                t8_locidx_list_append (t8_locidx_list_t * list,
                                           const t8_locidx_t entry);

/** Free all elements of a list.
 * \param [in, out] list An initialized list.
 * After calling this function \a list will not be initialized and
 * all its elements will have been free'd and its count will be zero.
 * \note Use this to free a list that was initialized with \ref t8_locidx_list_init
 */
void                t8_locidx_list_reset (t8_locidx_list_t * list);

/** Free all elements of a list and free the list pointer.
 * \param [in, out] list An initialized list.
 * After calling this function \a list will be the NULL pointer
 * and all the elements will be freed.
 * \note Use this to free a list that was allocated with \ref t8_locidx_list_new
 */
void                t8_locidx_list_destroy (t8_locidx_list_t ** plist);

#ifdef T8_ENABLE_DEBUG
/** Convert a list to a string for (debug) output.
 * \param [in]      list An initialized list.
 * \return               An allocated string that shows the contents of the list.
 * \note The user is responsible for deallocating the memory of the string with
 * T8_FREE.
 * \note The string will look like this: "NULL -> Entry0 -> Entry1 ... -> EntryN -> NULL"
 * If it is truncated then the list was too large to fit into the string. */
char               *t8_locidx_list_to_string (t8_locidx_list_t * list);
#endif

/** The t8_locidx_list_iterator_t can be used to iterate through the 
 *  entries of a t8_locidx_list_t.
 *  The workflow should be:
 *  \ref t8_locidx_list_iterator_init
 *  while not \ref t8_locidx_list_iterator_is_end repeat
 *  \ref t8_locidx_list_iterator_get_value
 *  maybe \ref t8_locidx_list_iterator_remove_entry
 *  \ref t8_locidx_list_iterator_next
 */
typedef struct
{
  t8_locidx_list_t   *list;
                          /**< The linked list to which this iterator belongs to */
  sc_link_t          *current;
                      /**< The currently active list entry. */
  sc_link_t          *prev;
                   /**< The previously active link entry. */
} t8_locidx_list_iterator_t;

/** Initialize an allocated iterator for a list.
 *  The iterator will start at the first item in the list.
 * \param [in]  An initialized list.
 */
void                t8_locidx_list_iterator_init (t8_locidx_list_t * list,
                                                  t8_locidx_list_iterator_t *
                                                  it);

/** Check whether an iterator is valid and associated to a given list.
 * \param [in] iterator The iterator to check.
 * \param [in] list     The list to which to \a iterator is expected to be associated to.
 * \return              True (non-zero) if \a iterator is valid and points to an element in \a list.
 *                      False if either \a iterator points to the end of \a list (all items have been
 *                      iterated through.) or \a iterator is not initialized.
 */
int                 t8_locidx_list_iterator_is_valid (const
                                                      t8_locidx_list_iterator_t
                                                      * it,
                                                      const t8_locidx_list_t *
                                                      list);

/** Let an iterator point to the next entry of its list.
 * \param [in, out] it  The iterator.
 * \note If the iterator already is at the end (\ref t8_locidx_list_iterator_is_end) nothing happens.
 */
void                t8_locidx_list_iterator_next (t8_locidx_list_iterator_t *
                                                  it);

/* Check whether a valid iterator is at the end of its list.
 * \param [in] iterator The iterator to check.
 * \return              True (non-zero) if \a iterator points to the end of its list.
 *                      That is, all elements have beed iterated through.
 * \note \a iterator must be valid before calling this function.
 */
int                 t8_locidx_list_iterator_is_end (const
                                                    t8_locidx_list_iterator_t
                                                    * it);

/** Return the value of the item that an iterator currently points to.
 * \param [in] it   The iterator.
 * \return          The value of the t8_locidx_t in the iterator's associated
 *                  list that it currently points to.
 */
t8_locidx_t         t8_locidx_list_iterator_get_value (const
                                                       t8_locidx_list_iterator_t
                                                       * it);

/* *INDENT-OFF* */
/** Remove the entry that an iterator points to from the list.
 * \param [in] it   The iterator.
 */
void                t8_locidx_list_iterator_remove_entry (t8_locidx_list_iterator_t 
                                                          * it);
/* *INDENT-ON* */

T8_EXTERN_C_END ();

#endif /* !T8_LOCIDX_LIST_H */
