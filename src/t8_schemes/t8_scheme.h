/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_scheme.h
 * This file defines the C interface of the t8_scheme class. For more information
 * refer to the C++ interface in \ref t8_scheme.hxx.
 */

#ifndef T8_SCHEME_H
#define T8_SCHEME_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** The scheme holds implementations for one or more element classes.
 *  Opaque pointer for C interface.
 *  Detailed documentation at \ref t8_scheme.
 */
typedef struct t8_scheme t8_scheme_c;

/** Increase the reference counter of a scheme.
 * \param [in,out] scheme       On input, this scheme must be alive, that is,
 *                              exist with positive reference count.
 */
void
t8_scheme_ref (t8_scheme_c *scheme);

/** Decrease the reference counter of a scheme.
 * If the counter reaches zero, this scheme is destroyed.
 * \param [in,out] pscheme      On input, the scheme pointed to must exist
 *                              with positive reference count.  If the
 *                              reference count reaches zero, the scheme is
 *                              destroyed and this pointer set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the scheme is not modified in other ways.
 */
void
t8_scheme_unref (t8_scheme_c **pscheme);

T8_EXTERN_C_END ();

#endif /* !T8_SCHEME_H */
