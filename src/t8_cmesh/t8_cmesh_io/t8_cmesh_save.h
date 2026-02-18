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

/** \file t8_cmesh_save.h
 * We define routines to save and load a cmesh to/from the file system.
 */

#ifndef T8_CMESH_SAVE_H
#define T8_CMESH_SAVE_H

#include <t8.h>

/** Increment this constant each time the file format changes.
 *  We can only read files that were written in the same format. */
#define T8_CMESH_FORMAT 0x0002

/** This enumeration contains all modes in which we can open a saved cmesh.
 * The cmesh can be loaded with more processes than it was saved and the
 * mode controls, which of the processes open files and distribute the data.
 */
typedef enum t8_load_mode {
  /** First mode. */
  T8_LOAD_FIRST = 0,
  /** In simple mode, the first n processes load the file */
  T8_LOAD_SIMPLE = T8_LOAD_FIRST,
  /** In BGQ mode, the file is loaded on n nodes and from one process of each node.
    * This needs MPI Version 3.1 or higher. */
  T8_LOAD_BGQ,
  /** Every n-th process loads a file. Handle with care, we introduce
   * it, since on Juqueen MPI-3 was not available.
   * The parameter n has to be passed as an extra parameter.
   * \see t8_cmesh_load_and_distribute */
  T8_LOAD_STRIDE,
  /** Number of modes in which we can open a saved cmesh. */
  T8_LOAD_COUNT
} t8_load_mode_t;

T8_EXTERN_C_BEGIN ();

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_SAVE_H */
