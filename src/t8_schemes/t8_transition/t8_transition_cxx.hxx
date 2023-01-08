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

/** \file t8_transition_cxx.hxx
 * This file is the point of entry for our transition implementation.
 * This scheme points to a consistent implementation of all element classes.
 */

#ifndef T8_TRANSITION_CXX_HXX
#define T8_TRANSITION_CXX_HXX

#include <t8_element_cxx.hxx>

T8_EXTERN_C_BEGIN ();

/** Return the subelement element implementation of t8code. */
t8_scheme_cxx_t    *t8_scheme_new_transition_cxx (void);

int                 t8_eclass_scheme_is_transition (t8_eclass_scheme_c *ts);

T8_EXTERN_C_END ();

#endif /* !T8_TRANSITION_CXX_HXX */
