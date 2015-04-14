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

/** \file t8_default_hex.h
 */

#ifndef T8_DEFAULT_HEX_H
#define T8_DEFAULT_HEX_H

#include <p8est.h>
#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

/** The structure holding a hexahedral element in the default scheme.
 * We make this definition public for interoperability of element classes.
 * We might want to put this into a private, scheme-specific header file.
 */
typedef p8est_quadrant_t t8_phex_t;

/** Provide an implementation for the hexahedral element class. */
t8_eclass_scheme_t *t8_default_scheme_new_hex (void);

T8_EXTERN_C_END ();

#endif /* !T8_DEFAULT_HEX_H */
