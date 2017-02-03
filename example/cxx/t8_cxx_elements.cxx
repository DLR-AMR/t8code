/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* In this program we test the C++ implementation of the element
 * interface versus the C implementation.
 * We create an array of NUM_ELEMENTS many quad-elements, assign a
 * unique Morton id to each one and refine each one.
 */

#include <t8.h>
#include <t8_element.h>
#include <t8_default.h>
#include <t8_element_cxx.hxx>
#include <t8_default_cxx.hxx>
#include <sc_statistics.h>
#include <sc_flops.h>

#define NUM_ELEMENTS 1e7
#define LEVEL         12        /* Must satisfy 4^LEVEL > NUM_ELEMENTS */

#define t8_cxx_print_size(_type) t8_debugf("sizeof (" #_type ") = %lu\n", \
  sizeof (_type))

/* The cxx version of the element refine test program */
void
t8_cxx_elements_test_scheme_cxx ()
{
  t8_scheme_cxx_t    *cxx_default_scheme = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *quad_scheme;
  t8_element_t      **elements = T8_ALLOC_ZERO (t8_element_t *, NUM_ELEMENTS);
  t8_element_t       *children[4];
  t8_locidx_t         ielement;

  t8_debugf ("Starting the C++ version.\n");

  /* Create the scheme class that stores the quadrant lookup functions */
  quad_scheme = cxx_default_scheme->eclass_schemes[T8_ECLASS_QUAD];
  /* Allocate memory to store the elements and 4 children */
  quad_scheme->t8_element_new (NUM_ELEMENTS, elements);
  quad_scheme->t8_element_new (4, children);

  /* Initialize each element and refine it.
   * The children are overwritten in each step */
  for (ielement = 0; ielement < NUM_ELEMENTS; ielement++) {
    quad_scheme->t8_element_set_linear_id (elements[ielement],
                                           LEVEL, ielement);
    quad_scheme->t8_element_children (elements[ielement], 4, children);
  }

  /* Clean-up */
  quad_scheme->t8_element_destroy (NUM_ELEMENTS, elements);
  quad_scheme->t8_element_destroy (4, children);
  T8_FREE (elements);
  t8_scheme_cxx_unref (&cxx_default_scheme);
}

void
t8_cxx_elements_test_scheme_c ()
{
  t8_scheme_t        *c_default_scheme = t8_scheme_new_default ();
  t8_eclass_scheme_t *quad_scheme;
  t8_element_t      **elements = T8_ALLOC_ZERO (t8_element_t *, NUM_ELEMENTS);
  t8_element_t       *children[4];
  t8_locidx_t         ielement;

  t8_debugf ("Starting the C version.\n");

  /* Create the scheme class that stores the quadrant lookup functions */
  quad_scheme = c_default_scheme->eclass_schemes[T8_ECLASS_QUAD];
  /* Allocate memory to store the elements and 4 children */
  t8_element_new (quad_scheme, NUM_ELEMENTS, elements);
  t8_element_new (quad_scheme, 4, children);

  /* Initialize each element and refine it.
   * The children are overwritten in each step */
  for (ielement = 0; ielement < NUM_ELEMENTS; ielement++) {
    t8_element_set_linear_id (quad_scheme, elements[ielement], LEVEL,
                              ielement);
    t8_element_children (quad_scheme, elements[ielement], 4, children);
  }

  /* Clean-up */
  t8_element_destroy (quad_scheme, NUM_ELEMENTS, elements);
  T8_FREE (elements);
  t8_scheme_unref (&c_default_scheme);
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[2];

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  /* Start timer */
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_cxx_elements_test_scheme_c ();

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "C Version");

  sc_flops_snap (&fi, &snapshot);

  t8_cxx_elements_test_scheme_cxx ();

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[1], snapshot.iwtime, "C++ Version");

  sc_stats_compute (sc_MPI_COMM_WORLD, 2, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_PRODUCTION, 2, stats, 1, 1);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
