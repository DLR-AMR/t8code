/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include "t8_default_common.h"

void
t8_default_scheme_mempool_destroy (t8_eclass_scheme_t * ts)
{
  T8_ASSERT (ts->ts_context != NULL);
  sc_mempool_destroy ((sc_mempool_t *) ts->ts_context);
}

void
t8_default_mempool_alloc (void *ts_context, int length, t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    elem[i] = (t8_element_t *) sc_mempool_alloc ((sc_mempool_t *) ts_context);
  }
}

void
t8_default_mempool_free (void *ts_context, int length, t8_element_t ** elem)
{
  int                 i;

  T8_ASSERT (ts_context != NULL);
  T8_ASSERT (0 <= length);
  T8_ASSERT (elem != NULL);

  for (i = 0; i < length; ++i) {
    sc_mempool_free ((sc_mempool_t *) ts_context, elem[i]);
  }
}
