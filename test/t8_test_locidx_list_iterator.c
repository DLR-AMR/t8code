/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include <t8_data/t8_locidx_list.h>

/* Builds a list and fills it with the values 
 * 0, 1, 2, ..., num_items-1
 */
static t8_locidx_list_t *
test_build_list (const t8_locidx_t num_items)
{
  t8_locidx_list_t   *list;
  t8_locidx_t         i;

  /* Initialize the list */
  list = t8_locidx_list_new ();
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (list),
                  "List is not initialized.");

  /* Insert the numbers 0,1,2,...num_items - 1 into the list.
   * Check the count each time. */
  for (i = 0; i < num_items; ++i) {
    t8_locidx_list_append (list, i);
    SC_CHECK_ABORT (t8_locidx_list_count (list) == i + 1,
                    "Wrong number of items in list.");
  }
  return list;
}

/* Allocate a list and check whether it is initialized */
static void
test_locidx_list_iterator_init ()
{
  /* Build a list with 10 entries */
  t8_locidx_list_t   *list = test_build_list (10);
  t8_locidx_list_iterator_t it;

  /* Initialize an iterator */
  t8_locidx_list_iterator_init (list, &it);
  /* Check whether the iterator is valid */
  SC_CHECK_ABORT (t8_locidx_list_iterator_is_valid (&it, list),
                  "Iterator is not valid.");

  t8_locidx_list_destroy (&list);
}

/* Allocate a list and check whether we can iterate over the
 * elements and access them.
 * We iterate multiple times, reinitializing the iterator each time. */
static void
test_locidx_list_iterator_iterate ()
{
  /* Build a list with 10 entries */
  const t8_locidx_t   num_values = 10;
  t8_locidx_list_t   *list = test_build_list (num_values);
  t8_locidx_list_iterator_t it;
  t8_locidx_t         i, returned_value;
  const int           num_its = 10;
  int                 count_it;

  /* Iterate over the values and check whether we can get them back */
  for (count_it = 0; count_it < num_its; ++count_it) {
    /* Initialize the iterator */
    t8_locidx_list_iterator_init (list, &it);
    /* Check whether the iterator is valid */
    SC_CHECK_ABORT (t8_locidx_list_iterator_is_valid (&it, list),
                    "Iterator is not valid.");
    for (i = 0; t8_locidx_list_iterator_is_valid (&it, list);
         t8_locidx_list_iterator_next (&it), ++i) {
      returned_value = t8_locidx_list_iterator_get_value (&it);
      SC_CHECK_ABORTF (returned_value == i, "Got wrong value. Expected %i "
                       "got %i.", i, returned_value);
    }
    SC_CHECK_ABORT (!t8_locidx_list_iterator_is_valid (&it, list),
                    "iterator should not be valid, but is.");
    SC_CHECK_ABORT (i == num_values,
                    "Number of items that we iterated through"
                    " does not match.");
  }

  t8_locidx_list_destroy (&list);
}

/* Allocate a list and iterator and check whether we can 
 * remove elements from the list via the iterator interface.
 * We iterate multiple times, reinitializing the iterator each time. */
static void
test_locidx_list_iterator_remove ()
{
  /* Build a list with 10 entries */
  const t8_locidx_t   num_values = 10;
  t8_locidx_list_t   *list = test_build_list (num_values);
  t8_locidx_list_iterator_t it;

  /* Initialize the iterator */
  t8_locidx_list_iterator_init (list, &it);
  /* Check whether the iterator is valid */
  SC_CHECK_ABORT (t8_locidx_list_iterator_is_valid (&it, list),
                  "Iterator is not valid.");

  /* Remove the first item in the list */
  t8_locidx_list_iterator_remove_entry (&it);
  /* Check whether the iterator is valid */
  SC_CHECK_ABORT (t8_locidx_list_iterator_is_valid (&it, list),
                  "Iterator is not valid.");
  /* Check that the list contains the correct number of entries */
  SC_CHECK_ABORT (t8_locidx_list_count (list) == num_values - 1,
                  "list has invalid number of entries.");
  /* Advance the iterator */
  t8_locidx_list_iterator_next (&it);
  /* Revome the current item.
   * The list should be 1, 2, 3, 4, ... with the
   * iterator pointing to 2
   */
  /* Remove the current item in the list */
  t8_locidx_list_iterator_remove_entry (&it);
  /* Check whether the iterator is valid */
  SC_CHECK_ABORT (t8_locidx_list_iterator_is_valid (&it, list),
                  "Iterator is not valid.");
  /* Check that the list contains the correct number of entries */
  SC_CHECK_ABORT (t8_locidx_list_count (list) == num_values - 2,
                  "list has invalid number of entries.");

  /* clean-up */
  t8_locidx_list_destroy (&list);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  test_locidx_list_iterator_init ();
  test_locidx_list_iterator_iterate ();
  test_locidx_list_iterator_remove ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
