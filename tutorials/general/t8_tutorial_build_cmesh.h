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

/** file t8_tutorial_build_cmesh.h
 */

#ifndef T8_TUTORIAL_BUILD_CMESH_H
#define T8_TUTORIAL_BUILD_CMESH_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** This is the main program of this example.
 */
int
t8_tutorial_build_cmesh_main (int argc, char **argv);

T8_EXTERN_C_END ();

#endif /* !T8_TUTORIAL_BUILD_CMESH_H */
