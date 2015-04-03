/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/** \file t8_default_tet.h
 */

#ifndef T8_DEFAULT_TET_H
#define T8_DEFAULT_TET_H

#include <t8_element.h>

typedef struct t8_tet
{
  int8_t              eclass;
          /**< We store the element class for compatibility with the pyramid. */
  int8_t              level;
  /* add coordinates etc. here */
}
t8_tet_t;

t8_eclass_scheme_t *t8_default_scheme_new_tet (void);

#endif /* !T8_DEFAULT_TET_H */
