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

/** \file t8_refcount.c
 * Implements reference counting functions declared in \ref t8_refcount.h.
 */

#include <t8_refcount.h>

void
t8_refcount_init (t8_refcount_t *rc)
{
  sc_refcount_init (rc, t8_get_package_id ());
}

t8_refcount_t *
t8_refcount_new (void)
{
  t8_refcount_t *rc;

  rc = T8_ALLOC (t8_refcount_t, 1);
  t8_refcount_init (rc);

  return rc;
}

void
t8_refcount_destroy (t8_refcount_t *rc)
{
  T8_ASSERT (!sc_refcount_is_active (rc));
  T8_FREE (rc);
}
