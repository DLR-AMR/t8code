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

/** \file t8_different_num_child_cxx.hxx
 *
 * This file is the point of entry for our different_num_child element implementation.
 *
 * In this scheme only the line elements are implemented and each element has
 * exactly one child. Thus, there is only a single element on each refinement level
 * and this element is the whole line.
 * This scheme exists mainly to demonstrate the possibilities of the scheme
 * interface.
 */

#ifndef T8_DIFFERENT_NUM_CHILD_CXX_HXX
#define T8_DIFFERENT_NUM_CHILD_CXX_HXX

#include <t8_element_cxx.hxx>

T8_EXTERN_C_BEGIN ();

/** Return the different_num_child element implementation of t8code. */
t8_scheme_cxx_t    *t8_scheme_new_different_num_child_cxx (void);

T8_EXTERN_C_END ();

#endif /* !T8_DIFFERENT_NUM_CHILD_H */
