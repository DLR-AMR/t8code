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

typedef struct t8_shmem_array *t8_shmem_array_t;

/** Defines the shared memory type that is best suited for t8code and the
 * current machine.
 * \see sc_shmem.h
 */
/* TODO: Change it as soon as we do not always need basic */
#if 0
#if defined(__bgq__)
#define T8_SHMEM_BEST_TYPE SC_SHMEM_BGQ
#elif defined(SC_ENABLE_MPIWINSHARED)
#define T8_SHMEM_BEST_TYPE SC_SHMEM_WINDOW
#else
#define T8_SHMEM_BEST_TYPE SC_SHMEM_BASIC
#endif
#endif
/* For testing we only use basic shmem type */
#define T8_SHMEM_BEST_TYPE SC_SHMEM_BASIC

T8_EXTERN_C_BEGIN ();

/** Initialize and allocate a shared memory array structure.
 * \param [in,out]      parray On input this pointer must be non-NULL.
 *                             On return this pointer is set to the new t8_shmem_array.
 * \param [in]          elem_size The size in bytes of an array element.
 * \param [in]          elem_count The total number of elements to allocate.
 * \param [in]          comm      The MPI communicator to be associated with the shmem_array.
 */
void                t8_shmem_array_init (t8_shmem_array_t * parray,
                                         size_t elem_size,
                                         size_t elem_count, sc_MPI_Comm comm);

/** Get the element size of a t8_shmem_array
 * \param [in]          array The array.
 * \return              The element size of \a array's elements.
 */
size_t              t8_shmem_array_get_elem_size (t8_shmem_array_t array);

/** Get the number of elements of a t8_shmem_array
 * \param [in]          array The array.
 * \return              The number of elements in \a array.
 */
size_t              t8_shmem_array_get_elem_count (t8_shmem_array_t array);

/** Return a pointer to the data of a shared memory array interpreted as
 * an t8_gloidx_t array.
 * \param [in]          array   The t8_shmem_array
 * \return              The data of \a array as t8_gloidx_t pointer.
 */
t8_gloidx_t        *t8_shmem_array_get_gloidx_array (t8_shmem_array_t array);

/** Free all memory associated with a t8_shmem_array.
 * \param [in,out]      parray  On input a pointer to a valid t8_shmem_array.
 *                      This array is freed and \a parray is set to NULL on return.
 */
void                t8_shmem_array_destroy (t8_shmem_array_t * parray);

T8_EXTERN_C_END ();

#endif /* !T8_SHMEM_H */
