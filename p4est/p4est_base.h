/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file p4est_base.h
 *
 * General support types and functions
 */


/* 
 * Note, this file was sligthly adapted in order to integrate it into the
 * t8code build system. Other than that no changes were made. The t8code authors.
 */

#ifndef P4EST_BASE_H
#define P4EST_BASE_H

/* Set configs as given by the t8code build system. */
#ifndef T8_CMAKE_BUILD
#include <t8_config.h>
#endif

#if (defined(T8_ENABLE_MPI))
# define P4EST_ENABLE_MPI T8_ENABLE_MPI
# endif

#if (defined(T8_ENABLE_MPIIO))
# define P4EST_ENABLE_MPIIO T8_ENABLE_MPIIO
# endif

/* Hard coding current p4est version as of 2024-04-23. */
# define P4EST_PACKAGE_VERSION "2.8.6"

/* Emulating p4est package id with the t8code package id. */
# define p4est_package_id t8_get_package_id ()

#include <sc_config.h>
#if \
  (defined (P4EST_ENABLE_MPI) && !defined (SC_ENABLE_MPI)) || \
  (!defined (P4EST_ENABLE_MPI) && defined (SC_ENABLE_MPI))
#error "MPI configured differently in p4est and libsc"
#endif
#if \
  (defined (P4EST_ENABLE_MPIIO) && !defined (SC_ENABLE_MPIIO)) || \
  (!defined (P4EST_ENABLE_MPIIO) && defined (SC_ENABLE_MPIIO))
#error "MPI I/O configured differently in p4est and libsc"
#endif

/* indirectly also include sc.h */
#include <sc_containers.h>
#define _p4est_const _sc_const

/*--------------------------------------------------------------------*/
/*------------------------ QUERY API CHANGES -------------------------*/
/*---- definitions to allow user code to query the p4est library -----*/

/** We do no longer dereference unneeded pointers in p4est_transfer_. */
#define P4EST_COMM_TRANSFER_NULL

/** The \ref p4est_connectivity_new_disk function now accepts a bool arg.
 * The same holds for \ref p4est_wrap_new_disk. */
#define P4EST_CONN_DISK_PERIODIC

/** The \ref p4est_connectivity_reorder_newid function exists. */
#define P4EST_CONN_REORDER_NEWID

/** The \ref p4est_search_local function replaces \ref p4est_search.
 * The latter function is still available with updated internal semantics.
 * Furthermore, we have added \ref p4est_search_partition to search
 * the parallel partition without accessing any local elements,
 * and \ref p4est_search_all that combines the two.
 */
#define P4EST_SEARCH_LOCAL

/** We expose the \ref p4est_vtk_write_cell_datav function. */
#define P4EST_VTK_CELL_DATAV

/*--------------------------------------------------------------------*/

SC_EXTERN_C_BEGIN;

/** Typedef for quadrant coordinates. */
typedef int32_t     p4est_qcoord_t;
#define p4est_qcoord_compare sc_int32_compare
#define P4EST_QCOORD_BITS 32
#define P4EST_MPI_QCOORD sc_MPI_INT
#define P4EST_VTK_QCOORD "Int32"
#define P4EST_F90_QCOORD INTEGER(KIND=C_INT32_T)
#define P4EST_QCOORD_MIN INT32_MIN
#define P4EST_QCOORD_MAX INT32_MAX
#define P4EST_QCOORD_1   ((p4est_qcoord_t) 1)
#define P4EST_QCOORD_ABS(x) ((p4est_qcoord_t) labs ((long) (x)))

/** Typedef for counting topological entities (trees, tree vertices). */
typedef int32_t     p4est_topidx_t;
#define p4est_topidx_compare sc_int32_compare
#define P4EST_TOPIDX_BITS 32
#define P4EST_MPI_TOPIDX sc_MPI_INT
#define P4EST_VTK_TOPIDX "Int32"
#define P4EST_F90_TOPIDX INTEGER(KIND=C_INT32_T)
#define P4EST_TOPIDX_MIN INT32_MIN
#define P4EST_TOPIDX_MAX INT32_MAX
#define P4EST_TOPIDX_FITS_32 1
#define P4EST_TOPIDX_1   ((p4est_topidx_t) 1)
#define P4EST_TOPIDX_ABS(x) ((p4est_topidx_t) labs ((long) (x)))

/** Typedef for processor-local indexing of quadrants and nodes. */
typedef int32_t     p4est_locidx_t;
#define p4est_locidx_compare sc_int32_compare
#define P4EST_LOCIDX_BITS 32
#define P4EST_MPI_LOCIDX sc_MPI_INT
#define P4EST_VTK_LOCIDX "Int32"
#define P4EST_F90_LOCIDX INTEGER(KIND=C_INT32_T)
#define P4EST_LOCIDX_MIN INT32_MIN
#define P4EST_LOCIDX_MAX INT32_MAX
#define P4EST_LOCIDX_1   ((p4est_locidx_t) 1)
#define P4EST_LOCIDX_ABS(x) ((p4est_locidx_t) labs ((long) (x)))

/** Typedef for globally unique indexing of quadrants. */
typedef int64_t     p4est_gloidx_t;
#define p4est_gloidx_compare sc_int64_compare
#define P4EST_GLOIDX_BITS 64
#define P4EST_MPI_GLOIDX sc_MPI_LONG_LONG_INT
#define P4EST_VTK_GLOIDX "Int64"
#define P4EST_F90_GLOIDX INTEGER(KIND=C_INT64_T)
#define P4EST_GLOIDX_MIN INT64_MIN
#define P4EST_GLOIDX_MAX INT64_MAX
#define P4EST_GLOIDX_1   ((p4est_gloidx_t) 1)
#define P4EST_GLOIDX_ABS(x) ((p4est_gloidx_t) llabs ((long long) (x)))

/** Tags for MPI messages */
typedef enum p4est_comm_tag
{
  P4EST_COMM_TAG_FIRST = SC_TAG_FIRST,
  P4EST_COMM_COUNT_PERTREE = SC_TAG_LAST,
  P4EST_COMM_BALANCE_FIRST_COUNT,
  P4EST_COMM_BALANCE_FIRST_LOAD,
  P4EST_COMM_BALANCE_SECOND_COUNT,
  P4EST_COMM_BALANCE_SECOND_LOAD,
  P4EST_COMM_PARTITION_GIVEN,
  P4EST_COMM_PARTITION_WEIGHTED_LOW,
  P4EST_COMM_PARTITION_WEIGHTED_HIGH,
  P4EST_COMM_PARTITION_CORRECTION,
  P4EST_COMM_GHOST_COUNT,
  P4EST_COMM_GHOST_LOAD,
  P4EST_COMM_GHOST_EXCHANGE,
  P4EST_COMM_GHOST_EXPAND_COUNT,
  P4EST_COMM_GHOST_EXPAND_LOAD,
  P4EST_COMM_GHOST_SUPPORT_COUNT,
  P4EST_COMM_GHOST_SUPPORT_LOAD,
  P4EST_COMM_GHOST_CHECKSUM,
  P4EST_COMM_NODES_QUERY,
  P4EST_COMM_NODES_REPLY,
  P4EST_COMM_SAVE,
  P4EST_COMM_LNODES_TEST,
  P4EST_COMM_LNODES_PASS,
  P4EST_COMM_LNODES_OWNED,
  P4EST_COMM_LNODES_ALL,
  P4EST_COMM_TAG_LAST
}
p4est_comm_tag_t;

/* some error checking possibly specific to p4est */
#ifdef P4EST_ENABLE_DEBUG
#define P4EST_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#define P4EST_EXECUTE_ASSERT_FALSE(expression)                          \
  do { int _p4est_i = (int) (expression);                               \
       SC_CHECK_ABORT (!_p4est_i, "Expected false: '" #expression "'"); \
  } while (0)
#define P4EST_EXECUTE_ASSERT_TRUE(expression)                           \
  do { int _p4est_i = (int) (expression);                               \
       SC_CHECK_ABORT (_p4est_i, "Expected true: '" #expression "'");   \
  } while (0)
#define P4EST_EXECUTE_ASSERT_INT(expression,ival)                       \
  do { int _p4est_i = (int) (expression);                               \
       SC_CHECK_ABORT ((ival) == _p4est_i,                              \
                       "Expected '" #ival "': '" #expression "'");      \
  } while (0)
#define P4EST_EXECUTE_ASSERT_TOPIDX(expression,tval)                    \
  do { p4est_topidx_t _p4est_t = (p4est_topidx_t) (expression);         \
       SC_CHECK_ABORT ((tval) == _p4est_t,                              \
                       "Expected '" #tval "': '" #expression "'");      \
  } while (0)
#define P4EST_DEBUG_EXECUTE(expression)                 \
  do { (void) (expression); } while (0)
#else
#define P4EST_ASSERT(c) SC_NOOP ()
#define P4EST_EXECUTE_ASSERT_FALSE(expression)          \
  do { (void) (expression); } while (0)
#define P4EST_EXECUTE_ASSERT_TRUE(expression)           \
  do { (void) (expression); } while (0)
#define P4EST_EXECUTE_ASSERT_INT(expression,ival)       \
  do { (void) (expression); } while (0)
#define P4EST_EXECUTE_ASSERT_TOPIDX(expression,tval)    \
  do { (void) (expression); } while (0)
#define P4EST_DEBUG_EXECUTE(expression) SC_NOOP ()
#endif

/* macros for memory allocation, will abort if out of memory */
/** allocate a \a t-array with \a n elements */
#define P4EST_ALLOC(t,n)          (t *) sc_malloc (p4est_package_id,    \
                                                   (n) * sizeof(t))
/** allocate a \a t-array with \a n elements and zero */
#define P4EST_ALLOC_ZERO(t,n)     (t *) sc_calloc (p4est_package_id,    \
                                                   (size_t) (n), sizeof(t))
/** reallocate the \a t-array \a p with \a n elements */
#define P4EST_REALLOC(p,t,n)      (t *) sc_realloc (p4est_package_id,   \
                                                    (p), (n) * sizeof(t))
/** duplicate a string */
#define P4EST_STRDUP(s)                 sc_strdup (p4est_package_id, (s))
/** free an allocated array */
#define P4EST_FREE(p)                   sc_free (p4est_package_id, (p))

/* log helper macros */
#define P4EST_GLOBAL_LOG(p,s)                           \
  SC_GEN_LOG (p4est_package_id, SC_LC_GLOBAL, (p), (s))
#define P4EST_LOG(p,s)                                  \
  SC_GEN_LOG (p4est_package_id, SC_LC_NORMAL, (p), (s))
void                P4EST_GLOBAL_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
void                P4EST_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
#ifndef __cplusplus
#define P4EST_GLOBAL_LOGF(p,f,...)                                      \
  SC_GEN_LOGF (p4est_package_id, SC_LC_GLOBAL, (p), (f), __VA_ARGS__)
#define P4EST_LOGF(p,f,...)                                             \
  SC_GEN_LOGF (p4est_package_id, SC_LC_NORMAL, (p), (f), __VA_ARGS__)
#endif

/* convenience global log macros will only print if identifier <= 0 */
#define P4EST_GLOBAL_TRACE(s) P4EST_GLOBAL_LOG (SC_LP_TRACE, (s))
#define P4EST_GLOBAL_LDEBUG(s) P4EST_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define P4EST_GLOBAL_VERBOSE(s) P4EST_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define P4EST_GLOBAL_INFO(s) P4EST_GLOBAL_LOG (SC_LP_INFO, (s))
#define P4EST_GLOBAL_STATISTICS(s) P4EST_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define P4EST_GLOBAL_PRODUCTION(s) P4EST_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_GLOBAL_ESSENTIAL(s) P4EST_GLOBAL_LOG (SC_LP_ESSENTIAL, (s))
#define P4EST_GLOBAL_LERROR(s) P4EST_GLOBAL_LOG (SC_LP_ERROR, (s))
void                P4EST_GLOBAL_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define P4EST_GLOBAL_TRACEF(f,...)                      \
  P4EST_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_LDEBUGF(f,...)                     \
  P4EST_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_GLOBAL_VERBOSEF(f,...)                    \
  P4EST_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_INFOF(f,...)                       \
  P4EST_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_GLOBAL_STATISTICSF(f,...)                         \
  P4EST_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_GLOBAL_PRODUCTIONF(f,...)                         \
  P4EST_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define P4EST_GLOBAL_ESSENTIALF(f,...)                          \
  P4EST_GLOBAL_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define P4EST_GLOBAL_LERRORF(f,...)                     \
  P4EST_GLOBAL_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define P4EST_GLOBAL_NOTICE     P4EST_GLOBAL_STATISTICS
#define P4EST_GLOBAL_NOTICEF    P4EST_GLOBAL_STATISTICSF

/* convenience log macros that are active on every processor */
#define P4EST_TRACE(s) P4EST_LOG (SC_LP_TRACE, (s))
#define P4EST_LDEBUG(s) P4EST_LOG (SC_LP_DEBUG, (s))
#define P4EST_VERBOSE(s) P4EST_LOG (SC_LP_VERBOSE, (s))
#define P4EST_INFO(s) P4EST_LOG (SC_LP_INFO, (s))
#define P4EST_STATISTICS(s) P4EST_LOG (SC_LP_STATISTICS, (s))
#define P4EST_PRODUCTION(s) P4EST_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_ESSENTIAL(s) P4EST_LOG (SC_LP_ESSENTIAL, (s))
#define P4EST_LERROR(s) P4EST_LOG (SC_LP_ERROR, (s))
void                P4EST_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define P4EST_TRACEF(f,...)                     \
  P4EST_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_LDEBUGF(f,...)                    \
  P4EST_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_VERBOSEF(f,...)                   \
  P4EST_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_INFOF(f,...)                      \
  P4EST_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_STATISTICSF(f,...)                        \
  P4EST_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_PRODUCTIONF(f,...)                        \
  P4EST_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define P4EST_ESSENTIALF(f,...)                         \
  P4EST_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define P4EST_LERRORF(f,...)                    \
  P4EST_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define P4EST_NOTICE            P4EST_STATISTICS
#define P4EST_NOTICEF           P4EST_STATISTICSF

/* extern declarations */
/** the libsc package id for p4est (set in p4est_init()) */
extern int          p4est_package_id;

static inline void
p4est_log_indent_push (void)
{
  sc_log_indent_push_count (p4est_package_id, 1);
}

static inline void
p4est_log_indent_pop (void)
{
  sc_log_indent_pop_count (p4est_package_id, 1);
}

/** Registers p4est with the SC Library and sets the logging behavior.
 * This function is optional.
 * This function must only be called before additional threads are created.
 * If this function is not called or called with log_handler == NULL,
 * the default SC log handler will be used.
 * If this function is not called or called with log_threshold == SC_LP_DEFAULT,
 * the default SC log threshold will be used.
 * The default SC log settings can be changed with sc_set_log_defaults ().
 */
void                p4est_init (sc_log_handler_t log_handler,
                                int log_threshold);

/** Return whether p4est has been initialized or not.
 * Keep in mind that \ref p4est_init is an optional function
 * but it helps with proper parallel logging.
 *
 * Currently there is no inverse to \ref p4est_init, and no way to deinit it.
 * This is ok since initialization generally does no harm.
 * Just do not call libsc's finalize function while p4est is still in use.
 *
 * \return          True if p4est has been initialized with a call to
 *                  \ref p4est_init and false otherwise.
 */
int                 p4est_is_initialized (void);

/** Check for a sufficiently recent zlib installation.
 * \return          True if zlib is detected in both sc and p4est.
 */
int                 p4est_have_zlib (void);

/** Query the package identity as registered in libsc.
 * \return          This is -1 before \ref p4est_init has been called
 *                  and a proper package identifier (>= 0) afterwards.
 */
int                 p4est_get_package_id (void);

/** Compute hash value for two p4est_topidx_t integers.
 * \param [in] tt     Array of (at least) two values.
 * \return            An unsigned hash value.
 */
/*@unused@*/
static inline unsigned
p4est_topidx_hash2 (const p4est_topidx_t * tt)
{
  uint32_t            a, b, c;

#if (P4EST_TOPIDX_FITS_32)
  a = (uint32_t) tt[0];
  b = (uint32_t) tt[1];
  c = 0;
#else
  a = (uint32_t) (tt[0] && 0xFFFFFFFF);
  b = (uint32_t) (tt[0] >> 32);
  c = (uint32_t) (tt[1] && 0xFFFFFFFF);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (tt[1] >> 32);
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

/** Compute hash value for three p4est_topidx_t integers.
 * \param [in] tt     Array of (at least) three values.
 * \return            An unsigned hash value.
 */
/*@unused@*/
static inline unsigned
p4est_topidx_hash3 (const p4est_topidx_t * tt)
{
  uint32_t            a, b, c;

#if (P4EST_TOPIDX_FITS_32)
  a = (uint32_t) tt[0];
  b = (uint32_t) tt[1];
  c = (uint32_t) tt[2];
#else
  a = (uint32_t) (tt[0] && 0xFFFFFFFF);
  b = (uint32_t) (tt[0] >> 32);
  c = (uint32_t) (tt[1] && 0xFFFFFFFF);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (tt[1] >> 32);
  b += (uint32_t) (tt[2] && 0xFFFFFFFF);
  c += (uint32_t) (tt[2] >> 32);
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

/** Compute hash value for four p4est_topidx_t integers.
 * \param [in] tt     Array of (at least) four values.
 * \return            An unsigned hash value.
 */
/*@unused@*/
static inline unsigned
p4est_topidx_hash4 (const p4est_topidx_t * tt)
{
  uint32_t            a, b, c;

#if (P4EST_TOPIDX_FITS_32)
  a = (uint32_t) tt[0];
  b = (uint32_t) tt[1];
  c = (uint32_t) tt[2];
  sc_hash_mix (a, b, c);
  a += (uint32_t) tt[3];
#else
  a = (uint32_t) (tt[0] && 0xFFFFFFFF);
  b = (uint32_t) (tt[0] >> 32);
  c = (uint32_t) (tt[1] && 0xFFFFFFFF);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (tt[1] >> 32);
  b += (uint32_t) (tt[2] && 0xFFFFFFFF);
  c += (uint32_t) (tt[2] >> 32);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (tt[3] && 0xFFFFFFFF);
  b += (uint32_t) (tt[3] >> 32);
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

/*@unused@*/
static inline int
p4est_topidx_is_sorted (p4est_topidx_t * t, int length)
{
  int                 i;

  for (i = 1; i < length; ++i) {
    if (t[i - 1] > t[i]) {
      return 0;
    }
  }
  return 1;
}

/*@unused@*/
static inline void
p4est_topidx_bsort (p4est_topidx_t * t, int length)
{
  int                 i, j;
  p4est_topidx_t      tswap;

  /* go through all elements except the last */
  for (i = length - 1; i > 0; --i) {
    /* bubble up the first element until before position i */
    for (j = 0; j < i; ++j) {
      if (t[j] > t[j + 1]) {
        tswap = t[j + 1];
        t[j + 1] = t[j];
        t[j] = tswap;
      }
    }
  }
  P4EST_ASSERT (p4est_topidx_is_sorted (t, length));
}

/*@unused@*/
static inline       uint64_t
p4est_partition_cut_uint64 (uint64_t global_num, int p, int num_procs)
{
  uint64_t            result;

  /* In theory, a double * double product should never overflow
     due to the 15-bit exponent used internally on x87 and above.
     Also in theory, 80-bit floats should be used internally,
     and multiply/divide associativity goes left-to-right.
     Still checking for funny stuff just to be sure. */

  P4EST_ASSERT (0 <= p && p <= num_procs);

  if (p == num_procs) {
    /* prevent roundoff error and division by zero */
    return global_num;
  }

  result = (uint64_t)
    (((long double) global_num * (double) p) / (double) num_procs);

  P4EST_ASSERT (result <= global_num);

  return result;
}

/*@unused@*/
static inline       p4est_gloidx_t
p4est_partition_cut_gloidx (p4est_gloidx_t global_num, int p, int num_procs)
{
  p4est_gloidx_t      result;

  /* In theory, a double * double product should never overflow
     due to the 15-bit exponent used internally on x87 and above.
     Also in theory, 80-bit floats should be used internally,
     and multiply/divide associativity goes left-to-right.
     Still checking for funny stuff just to be sure. */

  P4EST_ASSERT (global_num >= 0);
  P4EST_ASSERT (0 <= p && p <= num_procs);

  if (p == num_procs) {
    /* prevent roundoff error and division by zero */
    return global_num;
  }

  result = (p4est_gloidx_t)
    (((long double) global_num * (double) p) / (double) num_procs);

  P4EST_ASSERT (0 <= result && result <= global_num);

  return result;
}

/** Return the full version of p4est.
 *
 * \return          Return the version of p4est using the format
 *                  `VERSION_MAJOR.VERSION_MINOR.VERSION_POINT`,
 *                  where `VERSION_POINT` can contain dots and
 *                  characters, e.g. to indicate the additional
 *                  number of commits and a git commit hash.
 */
const char         *p4est_version (void);

/** Return the major version of p4est.
 *
 * \return          Return the major version of p4est.
 */
int                 p4est_version_major (void);

/** Return the minor version of p4est.
 *
 * \return          Return the minor version of p4est.
 */
int                 p4est_version_minor (void);

SC_EXTERN_C_END;

#endif /* !P4EST_BASE_H */
