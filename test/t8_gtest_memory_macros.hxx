/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_gtest_memory_macros.hxx
 * Provide memory macros.
 */

#ifndef T8_GTEST_MEMORY_MACROS_HXX
#define T8_GTEST_MEMORY_MACROS_HXX

#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <iostream>
#include <t8_schemes/t8_scheme.hxx>

/**
 * Register a package id for the t8code testsuite.
 * Used for attributes and memory diagnosis.
 */
void
t8_testsuite_register_package_id ();

/**
 * Returns the attribute package id of the t8code testsuite.
 * \return The package id.
 */
int
t8_testsuite_get_package_id ();

/** Allocate a \a t-array with \a n elements. */
#define T8_TESTSUITE_ALLOC(t, n) (t *) sc_malloc (t8_testsuite_get_package_id (), (n) * sizeof (t))

/** Allocate a \a t-array with \a n elements and init to zero. */
#define T8_TESTSUITE_ALLOC_ZERO(t, n) (t *) sc_calloc (t8_testsuite_get_package_id (), (size_t) (n), sizeof (t))

/** Deallocate a \a t-array. */
#define T8_TESTSUITE_FREE(p) sc_free (t8_testsuite_get_package_id (), (p))

/** Reallocate the \a t-array \a p with \a n elements. */
#define T8_TESTSUITE_REALLOC(p, t, n) (t *) sc_realloc (t8_testsuite_get_package_id (), (p), (n) * sizeof (t))

#endif /* T8_GTEST_MEMORY_MACROS_HXX */
