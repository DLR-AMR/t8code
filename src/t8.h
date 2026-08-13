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

/** \file t8.h
 * This is the administrative header file for t8code.
 * It includes standard C headers via subpackages.
 * It also provides application-level convenience code such as logging functions.
 */

#ifndef T8_H
#define T8_H

/* include config headers */
#ifndef T8_CMAKE_BUILD
#include <t8_config.h>
#endif
#include <sc_config.h>
#if (defined(T8_ENABLE_MPI) && !defined(SC_ENABLE_MPI)) || (!defined(T8_ENABLE_MPI) && defined(SC_ENABLE_MPI))
#error "MPI configured differently in t8code and libsc"
#endif
#if (defined(T8_ENABLE_MPIIO) && !defined(SC_ENABLE_MPIIO)) || (!defined(T8_ENABLE_MPIIO) && defined(SC_ENABLE_MPIIO))
#error "MPI I/O configured differently in t8code and libsc"
#endif

/* indirectly also include sc.h */
#include <sc_containers.h>

/** This macro opens the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8_EXTERN_C_BEGIN() SC_EXTERN_C_BEGIN

/** This macro closes the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8_EXTERN_C_END() SC_EXTERN_C_END

/* call this after including all headers */
T8_EXTERN_C_BEGIN ();

/** Portable way to use the const keyword determined by configure. */
#define t8_const _sc_const

/** Portable way to use the const keyword determined by configure. */
#define t8_restrict _sc_restrict

/** Widely used assertion. Only active in debug-mode. */
/* Note that we use the same definition as in sc.h. We are deliberately not using
 * #define T8_ASSERT SC_ASSERT
 * since then the assertion would not trigger if sc is not configured in debugging mode.
 * However, we want it to trigger any time t8code is in debugging mode, independent of sc.
 */
#if T8_ENABLE_DEBUG
#define T8_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define T8_ASSERT(c) SC_NOOP ()
#endif

/**Extended T8_ASSERT assertion with custom error message. Only active in debug-mode. */
#if T8_ENABLE_DEBUG
#define T8_ASSERTF(c, msg) SC_CHECK_ABORT ((c), "Assertion '" #c "': " msg)
#else
#define T8_ASSERTF(c, msg) SC_NOOP ()
#endif

/** Allocate a \a t-array with \a n elements. */
#define T8_ALLOC(t, n) (t *) sc_malloc (t8_get_package_id (), (n) * sizeof (t))

/** Allocate a \a t-array with \a n elements and init to zero. */
#define T8_ALLOC_ZERO(t, n) (t *) sc_calloc (t8_get_package_id (), (size_t) (n), sizeof (t))

/** Deallocate a \a t-array. */
#define T8_FREE(p) sc_free (t8_get_package_id (), (p))

/** Reallocate the \a t-array \a p with \a n elements. */
#define T8_REALLOC(p, t, n) (t *) sc_realloc (t8_get_package_id (), (p), (n) * sizeof (t))

/** A type for processor-local indexing. */
typedef int32_t t8_locidx_t;
/** The MPI Datatype of t8_locidx_t */
#define T8_MPI_LOCIDX sc_MPI_INT
/** Macro to get the absolute value of a t8_locidx_t */
#define T8_LOCIDX_ABS(x) ((t8_locidx_t) labs ((long) (x)))
/** Maximum possible value of a t8_locidx_t */
#define T8_LOCIDX_MAX INT32_MAX
/** Comparison function for t8_locidx_t */
#define t8_compare_locidx(v, w) sc_int32_compare (v, w)

/** A type for global indexing that holds really big numbers. */
typedef int64_t t8_gloidx_t;
/** The MPI Datatype of t8_gloidx_t */
#define T8_MPI_GLOIDX sc_MPI_LONG_LONG_INT
/** Macro to get the absolute value of a t8_gloidx_t */
#define T8_GLOIDX_ABS(x) ((t8_gloidx_t) llabs ((long long) (x)))
/** Maximum possible value of a t8_gloidx_t*/
#define T8_GLOIDX_MAX INT64_MAX
/** Comparison function for t8_gloidx_t */
#define t8_compare_gloidx(v, w) sc_int64_compare (v, w)

/** A type for storing SFC indices */
typedef uint64_t t8_linearidx_t;
/** The MPI datatype of t8_linearidx_t */
#define T8_MPI_LINEARIDX sc_MPI_UNSIGNED_LONG_LONG

/** The padding size is the size of a void pointer*/
#define T8_PADDING_SIZE (sizeof (void *))
/** Compute the number of bytes that have to be added to a given byte_count
 * such that it is a multiple of the padding size */
#define T8_ADD_PADDING(_x) ((T8_PADDING_SIZE - ((_x) % T8_PADDING_SIZE)) % T8_PADDING_SIZE)

/** Define machine precision for computations */
#define T8_PRECISION_EPS SC_EPS
/** Define square root of machine precision for computations */
#define T8_PRECISION_SQRT_EPS sqrt (T8_PRECISION_EPS)

/** Access multidimensional data on one-dimensional C arrays. */
/** Access onedimensional data on one-dimensional C arrays. */
#define T8_1D_TO_1D(nx, i) (i)
/** Access twodimensional data on one-dimensional C arrays. */
#define T8_2D_TO_1D(nx, ny, i, j) ((i) * (ny) + (j))
/** Access threedimensional data on one-dimensional C arrays. */
#define T8_3D_TO_1D(nx, ny, nz, i, j, k) (((i) * (ny) + (j)) * (nz) + (k))
/** Access fourdimensional data on one-dimensional C arrays. */
#define T8_4D_TO_1D(nx, ny, nz, nl, i, j, k, l) ((((i) * (ny) + (j)) * (nz) + (k)) * (nl) + (l))

/** Communication tags used internal to t8code. */
typedef enum {
  T8_MPI_TAG_FIRST = SC_TAG_FIRST,
  T8_MPI_PARTITION_CMESH = SC_TAG_LAST, /**< Used for coarse mesh partitioning */
  T8_MPI_PARTITION_FOREST,              /**< Used for forest partitioning */
  T8_MPI_GHOST_FOREST,                  /**< Used for for ghost layer creation */
  T8_MPI_GHOST_EXC_FOREST,              /**< Used for ghost data exchange */
  T8_MPI_CMESH_UNIFORM_BOUNDS_START,    /**< Used for cmesh uniform bounds computation. */
  T8_MPI_CMESH_UNIFORM_BOUNDS_END,      /**< Used for cmesh uniform bounds computation. */
  T8_MPI_TEST_ELEMENT_PACK_TAG,         /**< Used for testing mpi pack and unpack functionality */
  T8_MPI_TAG_LAST
} t8_MPI_tag_t;

/** Query the package identity as registered in libsc.
 * \return          This is -1 before \ref t8_init has been called
 *                  and a proper package identifier afterwards.
 */
int
t8_get_package_id (void);

/** Logging function parameterized by local/global category and priority.
 * \param [in] category     Either SC_LC_NORMAL for outputting on every rank
 *                          or SC_LC_GLOBAL for outputting on the root rank.
 * \param [in] priority     Please see sc.h for legal log priorities.
 * \param [in] fmt          Printf-style format string.
 * \param [in] ap           Argument list; see stdarg.h.
 */
void
t8_logv (int category, int priority, const char *fmt, va_list ap);

/** Logging function parameterized by local/global category and priority.
 * \param [in] category     Either SC_LC_NORMAL for outputting on every rank
 *                          or SC_LC_GLOBAL for outputting on the root rank.
 * \param [in] priority     Please see sc.h for legal log priorities.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_logf (int category, int priority, const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 3, 4)))
#endif
  ;

/** Add one space to the start of t8's default log format. */
void
t8_log_indent_push (void);

/** Remove one space from the start of a t8's default log format. */
void
t8_log_indent_pop (void);

/** Log a message on the root rank with priority SC_LP_ERROR.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_global_errorf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_ESSENTIAL.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_global_essentialf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_PRODUCTION.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_global_productionf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_INFO.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_global_infof (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_INFO.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_infof (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_PRODUCTION.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_productionf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_DEBUG.
 * \param [in] fmt          Printf-style format string.
 * \note This function does not print anything unless t8code was compiled
 * in debug mode (--enable-debug, T8_ENABLE_DEBUG was defined).
 */
void
t8_debugf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_ERROR.
 * \param [in] fmt          Printf-style format string.
 */
void
t8_errorf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/**
 * Set a custom logging function to be used by t8code.
 * When setting a custom logging function, the t8code internal logging function will be ignored.
 * \param [in] log_fcn      A function pointer to a logging function
 */
void
t8_set_external_log_fcn (void (*log_fcn) (int category, int priority, const char *msg));

/** Register t8code with libsc and print version and variable information.
 * \param [in] log_threshold Declared in sc.h.  SC_LP_DEFAULT is fine.
 *                           You can also choose from log levels SC_LP_*.
 */
void
t8_init (int log_threshold);

/** Return a pointer to an array element indexed by a t8_locidx_t.
 * \param [in] array The array of elements.
 * \param [in] index needs to be in [0]..[elem_count-1].
 * \return           A void * pointing to entry \a index in \a array.
 */
void *
t8_sc_array_index_locidx (const sc_array_t *array, const t8_locidx_t index);

/** Return values for subroutines to indicate if they fail or success. */
#define T8_SUBROUTINE_SUCCESS 1 /* true */
#define T8_SUBROUTINE_FAILURE 0 /* false */

/* call this at the end of a header file to match T8_EXTERN_C_BEGIN (). */
T8_EXTERN_C_END ();

#endif /* !T8_H */
