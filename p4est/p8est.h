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

/** \file p8est.h
 *
 * The top-level 3D p8est interface.
 *
 * \ingroup p8est
 */

/** \defgroup p8est p8est
 *
 * The 3D version of the p4est library.
 */

#ifndef P8EST_H
#define P8EST_H

/* p8est_connectivity.h includes p4est_base.h sc_containers.h */
#include <p8est_connectivity.h>

SC_EXTERN_C_BEGIN;

/** The finest level of the octree for representing nodes */
#define P8EST_OLD_MAXLEVEL 19   /* old means prior to mid-2020 */
#define P8EST_MAXLEVEL 30

/** The finest level of the octree for representing octants */
#define P8EST_OLD_QMAXLEVEL 18  /* old means prior to mid-2020 */
#define P8EST_QMAXLEVEL 29

/** The length of a side of the root quadrant */
#define P8EST_ROOT_LEN ((p4est_qcoord_t) 1 << P8EST_MAXLEVEL)

/** The length of a quadrant of level l */
#define P8EST_QUADRANT_LEN(l) ((p4est_qcoord_t) 1 << (P8EST_MAXLEVEL - (l)))

/** Create a mask of 1-bits from the left and maxlevel-level zero bits. */
#define P8EST_QUADRANT_MASK(l) (~(P8EST_QUADRANT_LEN (l) - 1))

/** The offset of the highest (farthest from the origin) quadrant at level l
 */
#define P8EST_LAST_OFFSET(l) (P8EST_ROOT_LEN - P8EST_QUADRANT_LEN (l))

/** The 3D quadrant (i.e., octant) datatype */
typedef struct p8est_quadrant
{
  /*@{ */
  p4est_qcoord_t      x, y, z;  /**< coordinates */
  /*@} */
  int8_t              level,    /**< level of refinement */
                      pad8;     /**< padding */
  int16_t             pad16;    /**< padding */
  union p8est_quadrant_data
  {
    void               *user_data;      /**< never changed by p4est */
    long                user_long;      /**< never changed by p4est */
    int                 user_int;       /**< never changed by p4est */
    p4est_topidx_t      which_tree;     /**< the tree containing the quadrant
                                             (used in auxiliary octants such
                                             as the ghost octants in
                                             p4est_ghost_t) */
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy1; /**< of ghost octants, store the tree and owner rank */
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy2; /**< of transformed octants, store the original tree and the
                 target tree */
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      local_num;
    }
    piggy3; /**< of ghost octants, store the tree and index in the owner's
                 numbering */
  }
  p; /**< a union of additional data attached to a quadrant */
}
p8est_quadrant_t;

/** The p8est tree datatype */
typedef struct p8est_tree
{
  sc_array_t          quadrants;             /**< locally stored quadrants */
  p8est_quadrant_t    first_desc,            /**< first local descendant */
                      last_desc;             /**< last local descendant */
  p4est_locidx_t      quadrants_offset;      /**< cumulative sum over earlier
                                                  trees on this processor
                                                  (locals only) */
  p4est_locidx_t      quadrants_per_level[P8EST_MAXLEVEL + 1];
                                             /**< locals only */
  int8_t              maxlevel;              /**< highest local quadrant level */
}
p8est_tree_t;

/** Data pertaining to selecting, inspecting, and profiling algorithms.
 * A pointer to this structure is hooked into the p8est main structure.
 * Declared in p8est_extended.h.  Used to profile important algorithms.
 */
typedef struct p8est_inspect p8est_inspect_t;

/** The p8est forest datatype */
typedef struct p8est
{
  sc_MPI_Comm         mpicomm;          /**< MPI communicator */
  int                 mpisize,          /**< number of MPI processes */
                      mpirank;          /**< this process's MPI rank */
  int                 mpicomm_owned;    /**< flag if communicator is owned */
  size_t              data_size;        /**< size of per-quadrant p.user_data
                     (see p8est_quadrant_t::p8est_quadrant_data::user_data) */
  void               *user_pointer;     /**< convenience pointer for users,
                                             never touched by p4est */

  long                revision;         /**< Gets bumped on mesh change */
  p4est_topidx_t      first_local_tree; /**< 0-based index of first local
                                             tree, must be -1 for an empty
                                             processor */
  p4est_topidx_t      last_local_tree;  /**< 0-based index of last local
                                             tree, must be -2 for an empty
                                             processor */
  p4est_locidx_t      local_num_quadrants;   /**< number of quadrants on all
                                                  trees on this processor */
  p4est_gloidx_t      global_num_quadrants;  /**< number of quadrants on all
                                                  trees on all processors */
  p4est_gloidx_t     *global_first_quadrant; /**< first global quadrant index
                                                  for each process and 1
                                                  beyond */
  p8est_quadrant_t   *global_first_position; /**< first smallest possible quad
                                                  for each process and 1
                                                  beyond */
  p8est_connectivity_t *connectivity; /**< connectivity structure, not owned */
  sc_array_t         *trees;          /**< array of all trees */

  sc_mempool_t       *user_data_pool; /**< memory allocator for user data */
  /*   WARNING: This is NULL if data size
     equals zero. */
  sc_mempool_t       *quadrant_pool;  /**< memory allocator for temporary
                                           quadrants */
  p8est_inspect_t    *inspect;        /**< algorithmic switches */
}
p8est_t;

/** Calculate local memory usage of a forest structure.
 * Not collective.  The memory used on the current rank is returned.
 * The connectivity structure is not counted since it is not owned;
 * use p8est_connectivity_memory_usage (p8est->connectivity).
 * \param [in] p8est    Valid forest structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_memory_used (p8est_t * p8est);

/** Return the revision counter of the forest.
 * Not collective, even though the revision value is the same on all ranks.
 * A newly created forest starts with a revision counter of zero.
 * Every refine, coarsen, partition, and balance that actually changes the mesh
 * increases the counter by one.  Operations with no effect keep the old value.
 * \param [in] p8est    The forest must be valid.
 * \return              Non-negative number.
 */
long                p8est_revision (p8est_t * p8est);

/** Callback function prototype to initialize the quadrant's user data.
 * \param [in] p8est         the forest
 * \param [in] which_tree    the tree containing \a quadrant
 * \param [in,out] quadrant  the quadrant to be initialized: if data_size > 0,
 *                           the data to be initialized is at
 *                           \a quadrant->p.user_data; otherwise, the
 *                           non-pointer user data (such as
 *                           \a quadrant->p.user_int) can be initialized
 */
typedef void        (*p8est_init_t) (p8est_t * p8est,
                                     p4est_topidx_t which_tree,
                                     p8est_quadrant_t * quadrant);

/** Callback function prototype to decide for refinement.
 * \param [in] p8est       the forest
 * \param [in] which_tree  the tree containing \a quadrant
 * \param [in] quadrant    the quadrant that may be refined
 * \return nonzero if the quadrant shall be refined.
 */
typedef int         (*p8est_refine_t) (p8est_t * p8est,
                                       p4est_topidx_t which_tree,
                                       p8est_quadrant_t * quadrant);

/** Callback function prototype to decide for coarsening.
 * \param [in] p8est       the forest
 * \param [in] which_tree  the tree containing \a quadrant
 * \param [in] quadrants   Pointers to 8 siblings in Morton ordering.
 * \return nonzero if the quadrants shall be replaced with their parent.
 */
typedef int         (*p8est_coarsen_t) (p8est_t * p8est,
                                        p4est_topidx_t which_tree,
                                        p8est_quadrant_t * quadrants[]);

/** Callback function prototype to calculate weights for partitioning.
 * \param [in] p8est       the forest
 * \param [in] which_tree  the tree containing \a quadrant
 * \return a 32bit integer >= 0 as the quadrant weight.
 * \note    Global sum of weights must fit into a 64bit integer.
 */
typedef int         (*p8est_weight_t) (p8est_t * p8est,
                                       p4est_topidx_t which_tree,
                                       p8est_quadrant_t * quadrant);

extern void        *P8EST_DATA_UNINITIALIZED;

/** set statically allocated quadrant to defined values */
#define P8EST_QUADRANT_INIT(q) \
  ((void) memset ((q), -1, sizeof (p8est_quadrant_t)))

/** Transform a quadrant coordinate into the space spanned by tree vertices.
 * \param [in] connectivity     Connectivity must provide the vertices.
 * \param [in] treeid           Identify the tree that contains x, y, z.
 * \param [in] x, y, z          Quadrant coordinates relative to treeid.
 * \param [out] vxyz            Transformed coordinates in vertex space.
 */
void                p8est_qcoord_to_vertex (p8est_connectivity_t *
                                            connectivity,
                                            p4est_topidx_t treeid,
                                            p4est_qcoord_t x,
                                            p4est_qcoord_t y,
                                            p4est_qcoord_t z, double vxyz[3]);

/** Create a new forest with an initial coarse mesh.
 * The new forest consists of equi-partitioned root quadrants.
 * When there are more processors than trees, some processors are empty.
 *
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note the p8est
 *                           does not take ownership of the memory.
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p8est
 *                           before init_fn is called the first time.
 *
 * \return This returns a valid forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this forest.
 */
p8est_t            *p8est_new (sc_MPI_Comm mpicomm,
                               p8est_connectivity_t * connectivity,
                               size_t data_size,
                               p8est_init_t init_fn, void *user_pointer);

/** Destroy a p8est.
 *
 * \note The connectivity structure is not destroyed with the p8est.
 */
void                p8est_destroy (p8est_t * p8est);

/** Make a deep copy of a p8est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 * If old and new data sizes are 0, the user_data field is copied regardless.
 * The inspect member of the copy is set to NULL.
 * The revision counter of the copy is set to zero.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \return  Returns a valid p8est that does not depend on the input,
 *                         except for borrowing the same connectivity.
 *                         Its revision counter is 0.
 */
p8est_t            *p8est_copy (p8est_t * input, int copy_data);

/** Reset user pointer and element data.
 * When the data size is changed the quadrant data is freed and allocated.
 * The initialization callback is invoked on each quadrant.
 * Old user_data content is disregarded.
 *
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 *                           May be NULL.
 * \param [in] user_pointer  Assign to the user_pointer member of the p8est
 *                           before init_fn is called the first time.
 */
void                p8est_reset_data (p8est_t * p8est, size_t data_size,
                                      p8est_init_t init_fn,
                                      void *user_pointer);

/** Refine a forest.
 * \param [in,out] p8est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] refine_fn Callback function that must return true if a quadrant
 *                       shall be refined.  If refine_recursive is true,
 *                       refine_fn is called for every existing and newly
 *                       created quadrant.  Otherwise, it is called for every
 *                       existing quadrant.  It is possible that a refinement
 *                       request made by the callback is ignored.  To catch
 *                       this case, you can examine whether init_fn gets
 *                       called, or use p8est_refine_ext in p8est_extended.h
 *                       and examine whether replace_fn gets called.
 * \param [in] init_fn   Callback function to initialize the user_data of newly
 *                       created quadrants, which is already allocated.  This
 *                       function pointer may be NULL.
 */
void                p8est_refine (p8est_t * p8est,
                                  int refine_recursive,
                                  p8est_refine_t refine_fn,
                                  p8est_init_t init_fn);

/** Coarsen a forest.
 * \param [in,out] p8est  The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p8est_coarsen (p8est_t * p8est,
                                   int coarsen_recursive,
                                   p8est_coarsen_t coarsen_fn,
                                   p8est_init_t init_fn);

/** 2:1 balance the size differences of neighboring elements in a forest.
 * \param [in,out] p8est  The p8est to be worked on.
 * \param [in] btype      Balance type (face, edge, or corner/full).  Examples:
 *                        Finite volume or discontinuous Galerkin methods only
 *                        require face balance.  Continuous finite element
 *                        methods usually require edge balance.  Corner balance
 *                        is almost never required mathematically; it just
 *                        produces a smoother mesh grading.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p8est_balance (p8est_t * p8est,
                                   p8est_connect_type_t btype,
                                   p8est_init_t init_fn);

/** Equally partition the forest.
 * The partition can be by element count or by a user-defined weight.
 *
 * The forest will be partitioned between processors such that they
 * have an approximately equal number of quadrants (or sum of weights).
 *
 * On one process, the function noops and does not call the weight callback.
 * Otherwise, the weight callback is called once per quadrant in order.
 *
 * \param [in,out] p8est      The forest that will be partitioned.
 * \param [in]     allow_for_coarsening Slightly modify partition such that
 *                            quadrant families are not split between ranks.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 *                            When running with mpisize == 1, never called.
 *                            Otherwise, called in order for all quadrants
 *                            if not NULL. A weighting function with constant
 *                            weight 1 on each quadrant is equivalent
 *                            to weight_fn == NULL but other constant weightings
 *                            may result in different uniform partitionings.
 */
void                p8est_partition (p8est_t * p8est,
                                     int allow_for_coarsening,
                                     p8est_weight_t weight_fn);

/** Compute the checksum for a forest.
 * Based on quadrant arrays only. It is independent of partition and mpisize.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p8est_checksum (p8est_t * p8est);

/** Compute a partition-dependent checksum for a forest.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p8est_checksum_partition (p8est_t * p8est);

/** Save the complete connectivity/p8est data to disk.
 *
 * This is a collective operation that all MPI processes need to call.  All
 * processes write into the same file, so the filename given needs to be
 * identical over all parallel invocations.
 *
 * By default, we write the current processor count and partition into the file
 * header.  This makes the file depend on mpisize.  For changing this see
 * p8est_save_ext() in p8est_extended.h.
 *
 * The revision counter is not saved to the file, since that would make files
 * different that come from different revisions but store the same mesh.
 *
 * \param [in] filename    Name of the file to write.
 * \param [in] p8est       Valid forest structure.
 * \param [in] save_data   If true, the element data is saved.
 *                         Otherwise, a data size of 0 is saved.
 * \note            Aborts on file errors.
 * \note            If p4est is not configured to use MPI-IO, some processes
 *                  return from this function before the file is complete, in
 *                  which case immediate read-access to the file may require a
 *                  call to sc_MPI_Barrier.
 */
void                p8est_save (const char *filename, p8est_t * p8est,
                                int save_data);

/** Load the complete connectivity/p8est structure from disk.
 *
 * This is a collective operation that all MPI processes need to call.  All
 * processes read from the same file, so the filename given needs to be
 * identical over all parallel invocations.
 *
 * By default, a file can only be loaded with the same number of processors
 * that it was stored with.  The defaults can be changed with p8est_load_ext()
 * in p8est_extended.h.
 *
 * The revision counter of the loaded p4est is set to zero.
 *
 * \param [in] filename         Name of the file to read.
 * \param [in] mpicomm          A valid MPI communicator.
 * \param [in] data_size        Size of data for each quadrant which can be
 *                              zero.  Then user_data_pool is set to NULL.
 *                              If data_size is zero, load_data is ignored.
 * \param [in] load_data        If true, the element data is loaded.  This is
 *                              only permitted if the saved data size matches.
 *                              If false, the stored data size is ignored.
 * \param [in] user_pointer     Assign to the user_pointer member of the p8est
 *                              before init_fn is called the first time.
 * \param [out] connectivity    Connectivity must be destroyed separately.
 * \return          Returns a valid forest structure. A pointer to a valid
 *                  connectivity structure is returned through the last
 *                  argument.
 * \note            Aborts on file errors or invalid file contents.
 */
p8est_t            *p8est_load (const char *filename, sc_MPI_Comm mpicomm,
                                size_t data_size, int load_data,
                                void *user_pointer,
                                p8est_connectivity_t ** connectivity);

/** Return a pointer to an array element indexed by a p4est_topidx_t.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */
/*@unused@*/
static inline p8est_tree_t *
p8est_tree_array_index (sc_array_t * array, p4est_topidx_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_tree_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p8est_tree_t *) (array->array +
                           sizeof (p8est_tree_t) * (size_t) it);
}

/** Return a pointer to a quadrant array element indexed by a size_t. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_quadrant_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_quadrant_t *) (array->array + sizeof (p8est_quadrant_t) * it);
}

/** Push the copy of a fully initialized quadrant onto a quadrant array.
 * \param [in,out] array        Valid array of quadrants is pushed to.
 * \param [in] qsrc     Pointer to a quadrant with fully initialized memory.
 *                      This means all bits in the quadrant mush have been
 *                      written to at least once, even the compiler padding.
 *                      This serves to make the function clean for valgrind.
 * \return              Newly allocated quadrant with contents of \a qsrc.
 */
static inline p8est_quadrant_t *
p8est_quadrant_array_push_copy (sc_array_t * array,
                                const p8est_quadrant_t *qsrc)
{
  p8est_quadrant_t *q;

  P4EST_ASSERT (array->elem_size == sizeof (p8est_quadrant_t));

  q = (p8est_quadrant_t *) sc_array_push (array);
  *q = *qsrc;
  return q;
}

/** Call sc_array_push for a quadrant array and fully initialize memory.
 * \param [in,out] array        Valid array of quadrants is pushed to.
 * \return              Pushed array element with fully initialized memory.
 *                      In this case, we're writing to all bits of it.
 *                      This serves to make the quadrant clean for valgrind.
 */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_array_push (sc_array_t * array)
{
  p8est_quadrant_t *q;

  P4EST_ASSERT (array->elem_size == sizeof (p8est_quadrant_t));

  q = (p8est_quadrant_t *) sc_array_push (array);
  P8EST_QUADRANT_INIT(q);
  return q;
}

/** Call sc_mempool_alloc for a mempool creating quadrants. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_mempool_alloc (sc_mempool_t * mempool)
{
  P4EST_ASSERT (mempool->elem_size == sizeof (p8est_quadrant_t));

  return (p8est_quadrant_t *) sc_mempool_alloc (mempool);
}

/** Call sc_list pop for a quadrant array. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_list_pop (sc_list_t * list)
{
  return (p8est_quadrant_t *) sc_list_pop (list);
}

SC_EXTERN_C_END;

#endif /* !P8EST_H */
