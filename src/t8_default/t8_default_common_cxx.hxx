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

/** \file t8_default_common_cxx.hxx
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_DEFAULT_COMMON_CXX_HXX
#define T8_DEFAULT_COMMON_CXX_HXX

#include <t8_element_cxx.hxx>

/* Macro to check whether a pointer (VAR) to a base class, comes from an
 * implementation of a child class (TYPE). */
#define T8_COMMON_IS_TYPE(VAR, TYPE) \
  ((dynamic_cast<TYPE> (VAR)) != NULL)

class               t8_default_scheme_common_c:public t8_eclass_scheme_c
{
public:
  /** Destructor for all default schemes */
  virtual ~ t8_default_scheme_common_c ();

  /** Compute the number of corners of a given element. */
  virtual int         t8_element_num_corners (const t8_element_t * elem);

  /** Allocate space for a bunch of elements. */
  virtual void        t8_element_new (int length, t8_element_t ** elem);

  /** Deallocate space for a bunch of elements. */
  virtual void        t8_element_destroy (int length, t8_element_t ** elem);
};

#endif /* !T8_DEFAULT_COMMON_CXX_HXX */
