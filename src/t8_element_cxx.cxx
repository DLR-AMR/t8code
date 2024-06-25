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

#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This belongs here since it uses c++ function,
 * see t8_element.c/.h */
void
t8_scheme_cxx_destroy (t8_scheme_cxx_t *s)
{
  int t;

  T8_ASSERT (s != NULL);
  T8_ASSERT (s->rc.refcount == 0);

  for (t = 0; t < T8_ECLASS_COUNT; ++t) {
    if (s->eclass_schemes[t] != NULL) {
      delete s->eclass_schemes[t];
    }
  }
  T8_FREE (s);
}

/* Default implementation for the element size */
size_t
t8_eclass_scheme::t8_element_size () const
{
  return element_size;
}

T8_EXTERN_C_END ();
