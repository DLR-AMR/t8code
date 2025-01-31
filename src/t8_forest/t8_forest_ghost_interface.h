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

#ifndef T8_FOREST_GHOST_INTERFACE_H
#define T8_FOREST_GHOST_INTERFACE_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** This typedef holds virtual functions for a particular ghost interface.
 * We need it so that we can use t8_ghost_interface_c pointers in .c files
 * without them seeing the actual C++ code (and then not compiling)
 */
typedef struct t8_forest_ghost_interface t8_forest_ghost_interface_c;

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GHOST_INTERFACE_H */
