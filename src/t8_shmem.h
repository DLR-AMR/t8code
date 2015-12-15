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

/** \file t8_shmem.h
 * We define basic shared memory routines.
 */

#ifndef T8_SHMEM_H
#define T8_SHMEM_H

#include <t8.h>
#include <sc_shmem.h>

T8_EXTERN_C_BEGIN ();

/** Defines the shared memory type that is best suited for t8code and the
 * current machine.
 * \see sc_shmem.h
 */
#if defined(__bgq__)
  #define T8_SHMEM_BEST_TYPE SC_SHMEM_BGQ
#elif defined(SC_ENABLE_MPIWINSHARED)
  #define T8_SHMEM_BEST_TYPE SC_SHMEM_WINDOW
#else
  #define T8_SHMEM_BEST_TYPE SC_SHMEM_BASIC
#endif

T8_EXTERN_C_END ();

#endif /* !T8_SHMEM_H */
