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

#ifndef T8_DEFAULT_COMMON_CXX_H
#define T8_DEFAULT_COMMON_CXX_H

#include <t8_element_cxx.h>

/** This class independent function assumes an sc_mempool_t as context.
 * It is suitable as the ts_destroy callback in \ref t8_eclass_scheme_t.
 * We assume that the mempool has been created with the correct element size.
 * \param [in,out] ts           The element class scheme context is destroyed.
 */
void                t8_default_scheme_mempool_destroy_cxx (t8_eclass_scheme_c
                                                           * ts);

#endif /* !T8_DEFAULT_COMMON_H */
