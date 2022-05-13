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

/** \file t8_default_tri.h
 * The default implementation for triangles.
 */

#ifndef T8_DEFAULT_TRI_H
#define T8_DEFAULT_TRI_H

#include <t8_element.h>

T8_EXTERN_C_BEGIN ();

/** Provide an implementation for the triangle element class.
 * It is written as a self-contained library in the t8_dtri_* files.
 */
t8_eclass_scheme_t *t8_default_scheme_new_tri (void);

T8_EXTERN_C_END ();

#endif /* !T8_DEFAULT_TET_H */
