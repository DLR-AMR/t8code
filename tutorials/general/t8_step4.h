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

/** file t8_step4.h
 * This is the header file to the step4 example of t8code. It collects
 * functions of t8_step4 that we reuse in other examples.
 * In this example we discuss the forest creation process in more detail including
 * partitioning and balancing a forest and creating a ghost layer.
 * The main program is t8_step4_main.
 * See \ref t8_step4_partition_balance_ghost.cxx for more details.
 */

#ifndef T8_STEP4_H
#define T8_STEP4_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** This is the main program of this example.
 */
int
t8_step4_main (int argc, char **argv);

T8_EXTERN_C_END ();

#endif /* !T8_STEP4_H */
