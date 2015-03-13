/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include "t8_default_quad.h"
#include "t8_default_hex.h"
#include "t8_default_tet.h"

t8_scheme_t        *
t8_scheme_new_default (void)
{
  t8_scheme_t        *s;

  s = T8_ALLOC_ZERO (t8_scheme_t, 1);
  s->type_schemes[T8_TYPE_QUAD] = t8_type_scheme_new_quad ();
  s->type_schemes[T8_TYPE_HEX] = t8_type_scheme_new_hex ();
  s->type_schemes[T8_TYPE_TET] = t8_type_scheme_new_tet ();

  return s;
}
