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
#if defined(__bgq__)
#define T8_SHMEM_BEST_TYPE SC_SHMEM_BGQ
#elif defined(SC_ENABLE_MPIWINSHARED)
#define T8_SHMEM_BEST_TYPE SC_SHMEM_WINDOW
#else
#define T8_SHMEM_BEST_TYPE SC_SHMEM_BASIC
#endif

T8_EXTERN_C_BEGIN ();

/** Initialize shared memory usage for a communicator.
 *  This sets up the intra- and internode communicators.
 * \param [in]          comm    The MPI Communicator.
 * \note This function needs to be called to enable shared memory usage for a communicator.
 * \note Calling this function multiple times with the same communicator is safe and does
 *  not change the behaviour.
 */
void
t8_shmem_init (sc_MPI_Comm comm);

#if T8_ENABLE_DEBUG
/* If you need this function outside of debugging mode, feel free
 * to remove the macro protection. */
/** Check whether a shared memory array was properly initialized.
 * \param [in]          array   A shared memory array.
 * \return non-zero if \a array is initialized correctly.
 */
int
t8_shmem_array_is_initialized (const t8_shmem_array_t array);
#endif

/** Finalize shared memory usage for a communicator.
 *  This destroys the intra- and internode communicators.
 * \param [in]          comm    The MPI Communicator.
 * \note Call this function if you initialized the communicator for shared memory usage
 *  and you are sure that you will not use it for shared memory again.
 * \note Calling this function multiple times with the same communicator is safe and does
 *  not change the behaviour.
 */
void
t8_shmem_finalize (sc_MPI_Comm comm);

/** Set a shared memory type of a communicator.
 * If the type was already set it is overwritten.
 * \param [in,out]      comm    The MPI Communicator
 * \param [in]          type    A shared memory type.
 */
void
t8_shmem_set_type (sc_MPI_Comm comm, sc_shmem_type_t type);

/** Initialize and allocate a shared memory array structure.
 * \param [in,out]      parray On input this pointer must be non-NULL.
 *                             On return this pointer is set to the new t8_shmem_array.
 * \param [in]          elem_size The size in bytes of an array element.
 * \param [in]          elem_count The total number of elements to allocate.
 * \param [in]          comm      The MPI communicator to be associated with the shmem_array.
 *                                If not set, the shared memory type will be set to T8_SHMEM_BEST_TYPE.
 */
void
t8_shmem_array_init (t8_shmem_array_t *parray, size_t elem_size, size_t elem_count, sc_MPI_Comm comm);

/** Enable writing mode for a shmem array. Only some processes may be allowed
 *  to write into the array, which is indicated by the return value being non-zero.
 * \param [in,out]      array Initialized array. Writing will be enabled on certain processes.
 * \return                    True if the calling process can write into the array.
 * \note This function is MPI collective.
 */
int
t8_shmem_array_start_writing (t8_shmem_array_t array);

/** Disable writing mode for a shmem array.
 * \param [in,out]      array Initialized with writing mode enabled.
 * \see t8_shmem_array_start_writing.
 * \note This function is MPI collective.
 */
void
t8_shmem_array_end_writing (t8_shmem_array_t array);

/** Set an entry of a t8_shmem array that is used to store t8_gloidx_t.
 * The array must have writing mode enabled \ref t8_shmem_array_start_writing.
 * \param [in,out]      array   The array to be modified.
 * \param [in]          index   The array entry to be modified.
 * \param [in]          value   The new value to be set.
 */
void
t8_shmem_array_set_gloidx (t8_shmem_array_t array, int index, t8_gloidx_t value);

/** Copy the contents of one t8_shmem array into another.
 * \param [in,out]      dest    The array in which \a source should be copied.
 * \param [in]          source  The array to copy.
 * \note \a dest must be initialized and match in element size and element count to \a source.
 * \note \a dest must have writing mode disabled.
 */
void
t8_shmem_array_copy (t8_shmem_array_t dest, t8_shmem_array_t source);

/** Fill a t8_shmem array with an allgather.
 *
 * \param[in] sendbuf         the source from this process
 * \param[in] sendcount       the number of items to allgather
 * \param[in] sendtype        the type of items to allgather
 * \param[in,out] recvbuf     the destination shmem array
 * \param[in] recvcount       the number of items to allgather
 * \param[in] recvtype        the type of items to allgather
 * \note Writing mode must be disabled for \a recvarray.
 */
void
t8_shmem_array_allgather (const void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, t8_shmem_array_t recvarray,
                          int recvcount, sc_MPI_Datatype recvtype);

/**
 * Fill a t8_shmem array with an Allgatherv
 * Computes the recvcount-array and displacement-array for each rank of a node using the
 * sendcount.
 * The total number of items of each node is then used to compute the
 * recvcount-array and displacement-array between nodes. 
 * Use t8_shmem_array_allgather if the sendcount is equal on all procs for better scaling. 
 * 
 * \param[in] sendbuf         the source from this process
 * \param[in] sendcount       the number of items to gather on this proc
 * \param[in] sendtype        the type of items to gather
 * \param[in, out] recvarray  array of type recvtype where the data gets written to
 * \param[in] recvtype        the type of items to receive
 * \param[in] comm            the mpi communicator
 * 
 */
void
t8_shmem_array_allgatherv (void *sendbuf, const int sendcount, sc_MPI_Datatype sendtype, t8_shmem_array_t recvarray,
                           sc_MPI_Datatype recvtype, sc_MPI_Comm comm);

/**
 * Fill a t8_shmem array with an Allgather of the prefix operation over all 
 * processes. 
 * 
 * The receive array will be
 * (0, send0, send0 op send1, send0 op send1 op send2, ...)
 * 
 * \note the first entry of \a recvarray will be set to 0 using memset. 
 * The entry can be changed after calling t8_shmem_array_prefix 
 * 
 * \param[in] sendbuf           The source from this process
 * \param[in, out] recvarray    The destination shmem array
 * \param[in] count             The number of items to gather
 * \param[in] type              The type of items to gather
 * \param[in] op                The operation to be used
 * \param[in] comm              The MPI communicator
 */
void
t8_shmem_array_prefix (const void *sendbuf, t8_shmem_array_t recvarray, const int count, sc_MPI_Datatype type,
                       sc_MPI_Op op, sc_MPI_Comm comm);

/** Return the MPI communicator associated with a shmem array.
 * \param [in]          array The shmem_array to be queried.
 * \return              The MPI communicator stored at \a array.
 */
sc_MPI_Comm
t8_shmem_array_get_comm (t8_shmem_array_t array);

/** Get the element size of a t8_shmem_array
 * \param [in]          array The array.
 * \return              The element size of \a array's elements.
 */
size_t
t8_shmem_array_get_elem_size (t8_shmem_array_t array);

/** Get the number of elements of a t8_shmem_array
 * \param [in]          array The array.
 * \return              The number of elements in \a array.
 */
size_t
t8_shmem_array_get_elem_count (t8_shmem_array_t array);

/** Return a read-only pointer to the data of a shared memory array interpreted as
 * an t8_gloidx_t array.
 * \param [in]          array   The t8_shmem_array
 * \return              The data of \a array as t8_gloidx_t pointer.
 * \note Writing mode must be disabled for \a array.
 */
const t8_gloidx_t *
t8_shmem_array_get_gloidx_array (t8_shmem_array_t array);

/** Return a pointer to the data of a shared memory array interpreted as
 * an t8_gloidx_t array. The array must have writing enabled \ref t8_shmem_array_start_writing
 * and you should not write into the memory after \ref t8_shmem_array_end_writing was called.
 * \param [in]          array   The t8_shmem_array
 * \return              The data of \a array as t8_gloidx_t pointer.
 */
t8_gloidx_t *
t8_shmem_array_get_gloidx_array_for_writing (t8_shmem_array_t array);

/** Return an entry of a shared memory array that stores t8_gloidx_t.
 * \param [in]          array   The t8_shmem_array
 * \param [in]          index   The index of the entry to be queried.
 * \return              The \a index-th entry of \a array as t8_gloidx_t.
 * \note Writing mode must be disabled for \a array.
 */
t8_gloidx_t
t8_shmem_array_get_gloidx (t8_shmem_array_t array, int index);

/** Return a pointer to the data array of a t8_shmem_array.
 * \param [in]          array The t8_shmem_array.
 * \return                    A pointer to the data array of \a array.
 * \note Writing mode must be disabled for \a array.
 */
const void *
t8_shmem_array_get_array (t8_shmem_array_t array);

/** Return a read-only pointer to an element in a t8_shmem_array.
 * \param [in]          array The t8_shmem_array.
 * \param [in]          index The index of an element.
 * \return              A pointer to the element at \a index in \a array.
 * \note You should not modify the value.
 * \note Writing mode must be disabled for \a array.
 */
const void *
t8_shmem_array_index (t8_shmem_array_t array, size_t index);

/** Return a pointer to an element in a t8_shmem_array in writing mode.
 * \param [in]          array The t8_shmem_array.
 * \param [in]          index The index of an element.
 * \return              A pointer to the element at \a index in \a array.
 * \note You can modify the value before the next call to \ref t8_shmem_array_end_writing.
 * \note Writing mode must be enabled for \a array.
 */
void *
t8_shmem_array_index_for_writing (t8_shmem_array_t array, size_t index);

/* TODO: implement and comment */
/* returns true if arrays are equal 
 * \note Writing mode must be disabled for \a array_a and \a array_b.
 */
int
t8_shmem_array_is_equal (t8_shmem_array_t array_a, t8_shmem_array_t array_b);

/** Free all memory associated with a t8_shmem_array.
 * \param [in,out]      parray  On input a pointer to a valid t8_shmem_array.
 *                      This array is freed and \a parray is set to NULL on return.
 */
void
t8_shmem_array_destroy (t8_shmem_array_t *parray);

T8_EXTERN_C_END ();

#endif /* !T8_SHMEM_H */
