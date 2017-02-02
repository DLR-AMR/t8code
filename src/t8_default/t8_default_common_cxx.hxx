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

/** \file t8_default_common.h
 * We provide some functions that are useful across element classes.
 */

#ifndef T8_DEFAULT_COMMON_CXX_HXX
#define T8_DEFAULT_COMMON_CXX_HXX

#include <t8_element_cxx.hxx>

class t8_default_scheme_common_c:public t8_eclass_scheme_c
{
public:
  /** Destructor for all default schemes */
  virtual ~t8_default_scheme_common_c ();
};

#endif /* !T8_DEFAULT_COMMON_CXX_HXX */
