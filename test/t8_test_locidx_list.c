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

/* Allocate a list and check whether it is initialized */
static void
test_locidx_list_new ()
{
  t8_locidx_list_t   *list;

  /* Allocate the list and check whether it is initialized */
  list = t8_locidx_list_new ();
  SC_CHECK_ABORT (list != NULL, "List new returned NULL.");
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (list),
                  "List is not initialized.");
  /* Check that the list is empty */
  SC_CHECK_ABORT (t8_locidx_list_count (list) == 0, "List is not empty.");
  /* Destroy the list */
  t8_locidx_list_destroy (&list);
  /* Check that the pointer was set to NULL */
  SC_CHECK_ABORT (list == NULL, "After destruction list is not NULL.");
}

/* Initialize a list and check whether it is initialized */
static void
test_locidx_list_init ()
{
  t8_locidx_list_t    list;

  /* Initialize the list and check */
  t8_locidx_list_init (&list);
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (&list),
                  "List is not initialized.");
  /* Check that the list is empty */
  SC_CHECK_ABORT (t8_locidx_list_count (&list) == 0, "List is not empty.");
  /* Reset the list */
  t8_locidx_list_reset (&list);
}

/* Build a list, insert integers, pop integers, and count
 * the number of elements. Checking each time whether the
 * results are correct. */
static void
test_locidx_list_append_pop_count ()
{
  const t8_locidx_t   max = 1000;
  t8_locidx_t         i;
  t8_locidx_list_t   *list;

  /* Initialize an empty list */
  list = t8_locidx_list_new ();
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (list),
                  "List is not initialized.");

  /* Check that the list is empty */
  SC_CHECK_ABORT (t8_locidx_list_count (list) == 0, "List is not empty.");

  /* Insert the numbers 0,1,2,...max - 1 into the list.
   * Check the count each time. */
  for (i = 0; i < max; ++i) {
    t8_locidx_list_append (list, i);
    SC_CHECK_ABORT (t8_locidx_list_count (list) == i + 1,
                    "Wrong number of items in list.");
  }
  /* Pop the entries from the list, check whether they are correct
   * and whether the count is correct. */
  for (i = 0; i < max; ++i) {
    t8_locidx_t         popped = t8_locidx_list_pop (list);
    SC_CHECK_ABORTF (popped == i,
                     "list_pop returned wrong item. Expected %i got %i.", i,
                     popped);
    SC_CHECK_ABORT (t8_locidx_list_count (list) == max - i - 1,
                    "Wrong number of items in list.");
  }
  /* Clean up memory */
  t8_locidx_list_destroy (&list);
  /* Check that the pointer was set to NULL */
  SC_CHECK_ABORT (list == NULL, "After destruction list is not NULL.");
}

/* We initialize a list, reset it and initialize it again */
void
test_locidx_list_reinit ()
{
  t8_locidx_list_t    list;

  /* Initialize the list */
  t8_locidx_list_init (&list);
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (&list),
                  "List is not initialized.");
  /* Add an item to the list */
  t8_locidx_list_append (&list, 42);
  /* Check count */
  SC_CHECK_ABORT (t8_locidx_list_count (&list) == 1,
                  "Wrong number of items in list.");
  /* Reset the list */
  t8_locidx_list_reset (&list);
  /* Re init the list */
  t8_locidx_list_init (&list);
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (&list),
                  "List is not initialized.");
  /* Check that it is empty */
  SC_CHECK_ABORT (t8_locidx_list_count (&list) == 0, "List is not empty.");
  t8_locidx_list_reset (&list);
}

/* We initialize a list and fill it with entries.
 * We then remove several entries. */
void
test_locidx_list_remove ()
{
  t8_locidx_list_t    list;
  const t8_locidx_t   max = 10;
  t8_locidx_t         i, removed_item;

  /* Initialize the list */
  t8_locidx_list_init (&list);
  SC_CHECK_ABORT (t8_locidx_list_is_initialized (&list),
                  "List is not initialized.");

  /* Insert the numbers 0,1,2,...max - 1 into the list.
   * Check the count each time. */
  for (i = 0; i < max; ++i) {
    t8_locidx_list_append (&list, i);
    SC_CHECK_ABORT (t8_locidx_list_count (&list) == i + 1,
                    "Wrong number of items in list.");
  }

  /* The list now contains 0, 1, 2, 3, 4, ...
   * We want to remove the second entry (1) from it.
   * In order to do so, we must get a pointer to the first entry.
   */
  removed_item = t8_locidx_list_remove (&list, list.list.first);
  SC_CHECK_ABORT (removed_item == 1, "Removed wrong item from the list.");
  /* Check that one item less is in the list */
  SC_CHECK_ABORT (t8_locidx_list_count (&list) == max - 1,
                  "Wrong number of items in list.");

  /* We now remove the first item of the list until the list is empty */
  for (i = 0; i < max - 1; ++i) {
    const t8_locidx_t   expect = (i == 0 ? 0 : i + 1);  /* The item that we expect to get */
    removed_item = t8_locidx_list_remove (&list, NULL);
    /* Check that we remove the correct item */
    SC_CHECK_ABORTF (removed_item == expect,
                     "Removed wrong item from the list. Expected %i got %i\n",
                     expect, removed_item);
    /* Check that the list contains the correct number of elements */
    SC_CHECK_ABORT (t8_locidx_list_count (&list) == max - 2 - i,
                    "Wrong number of items in list.");
  }
  t8_locidx_list_reset (&list);
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

  test_locidx_list_new ();
  test_locidx_list_init ();
  test_locidx_list_append_pop_count ();
  test_locidx_list_reinit ();
  test_locidx_list_remove ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
