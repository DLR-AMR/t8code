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

/** \file p8est_ghost.h
 *
 * passing quadrants and data to neighboring processes
 *
 * \ingroup p8est
 */

#ifndef P8EST_GHOST_H
#define P8EST_GHOST_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** quadrants that neighbor the local domain */
typedef struct
{
  int                 mpisize;
  p4est_topidx_t      num_trees;
  p8est_connect_type_t btype; /**< which neighbors are in the ghost layer */

  /** An array of quadrants which make up the ghost layer around \a
   * forest.  Their piggy3 data member is filled with their owner's tree
   * and local number (cumulative over trees).  Quadrants are ordered in \c
   * p8est_quadrant_compare_piggy order.  These are quadrants inside the
   * neighboring tree, i.e., \c p8est_quadrant_is_inside() is true for the
   * quadrant and the neighboring tree.
   */
  sc_array_t          ghosts; /**< array of p8est_quadrant_t type */
  p4est_locidx_t     *tree_offsets;     /**< num_trees + 1 ghost indices */
  p4est_locidx_t     *proc_offsets;     /**< mpisize + 1 ghost indices */

  /** An array of local quadrants that touch the parallel boundary from the
   * inside, i.e., that are ghosts in the perspective of at least one other
   * processor.  The storage convention is the same as for \c ghosts above.
   */
  sc_array_t          mirrors; /**< array of p8est_quadrant_t type */
  p4est_locidx_t     *mirror_tree_offsets;      /**< num_trees + 1 mirror indices */
  p4est_locidx_t     *mirror_proc_mirrors;      /**< indices into mirrors grouped by
                                                   outside processor rank and
                                                   ascending within each rank */
  p4est_locidx_t     *mirror_proc_offsets;      /**< mpisize + 1 indices into 
                                                   mirror_proc_mirrors */
  p4est_locidx_t     *mirror_proc_fronts;       /**< like mirror_proc_mirrors,
                                                   but limited to the
                                                   outermost octants.  This is
                                                   NULL until
                                                   p8est_ghost_expand is
                                                   called */
  p4est_locidx_t     *mirror_proc_front_offsets;        /**< NULL until
                                                           p8est_ghost_expand is
                                                           called */
}
p8est_ghost_t;

/** Examine if a ghost structure is valid as described above.
 * Test if within a ghost-structure the arrays ghosts and mirrors are in
 * p8est_quadrant_compare_piggy order.
 * Test if local_num in piggy3 data member of the quadrants in ghosts and
 * mirrors are in ascending order (ascending within each rank for ghost).
 *
 * Test if the p4est_locidx_t arrays are in ascending order
 * (for mirror_proc_mirrors ascending within each rank)
 * \param [in] p8est    the forest.
 * \param [in] ghost    Ghost layer structure.
 * \return true if \a ghost is valid
 */
int                 p8est_ghost_is_valid (p8est_t * p8est,
                                          p8est_ghost_t * ghost);

/** Calculate the memory usage of the ghost layer.
 * \param [in] ghost    Ghost layer structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_ghost_memory_used (p8est_ghost_t * ghost);

/** Gets the processor id of a quadrant's owner.
 * The quadrant can lie outside of a tree across faces (and only faces).
 *
 * \param [in] p8est  The forest in which to search for a quadrant.
 * \param [in] treeid The tree to which the quadrant belongs.
 * \param [in] face   Supply a face direction if known, or -1 otherwise.
 * \param [in] q      The quadrant that is being searched for.
 *
 * \return Processor id of the owner
 *                or -1 if the quadrant lies outside of the mesh.
 *
 * \warning Does not work for tree edge or corner neighbors.
 */
int                 p8est_quadrant_find_owner (p8est_t * p8est,
                                               p4est_topidx_t treeid,
                                               int face,
                                               const p8est_quadrant_t * q);

/** Builds the ghost layer.
 *
 * This will gather the quadrants from each neighboring proc to build one layer
 * of face, edge and corner based ghost elements around the ones they own.
 *
 * \param [in] p8est            The forest for which the ghost layer will be
 *                              generated.
 * \param [in] btype            Which ghosts to include (across face, edge,
 *                              or corner/full).
 * \return                      A fully initialized ghost layer.
 */
p8est_ghost_t      *p8est_ghost_new (p8est_t * p8est,
                                     p8est_connect_type_t btype);

/** Generate an empty ghost layer.
 * This ghost layer pretends that there are no parallel neighbor elements.
 * It is useful if general algorithms should be run with local data only.
 * \param [in] p8est    Valid forest.
 * \param [in] ctype    Ghosts to include (none, across face, face/corner).
 *                      This variable must be valid but has no effect.
 * \return              Valid ghost layer of zero ghost elements.
 */
p8est_ghost_t      *p8est_ghost_new_local (p8est_t * p8est,
                                           p8est_connect_type_t ctype);

/** Frees all memory used for the ghost layer. */
void                p8est_ghost_destroy (p8est_ghost_t * ghost);

/** Conduct binary search for exact match on a range of the ghost layer.
 * \param [in] ghost            The ghost layer.
 * \param [in] which_proc       The owner of the searched quadrant.  Can be -1.
 * \param [in] which_tree       The tree of the searched quadrant.  Can be -1.
 * \param [in] q                Valid quadrant is searched in the ghost layer.
 * \return                      Offset in the ghost layer, or -1 if not found.
 */
ssize_t             p8est_ghost_bsearch (p8est_ghost_t * ghost,
                                         int which_proc,
                                         p4est_topidx_t which_tree,
                                         const p8est_quadrant_t * q);

/** Conduct binary search for ancestor on range of the ghost layer.
 * \param [in] ghost            The ghost layer.
 * \param [in] which_proc       The owner of the searched quadrant.  Can be -1.
 * \param [in] which_tree       The tree of the searched quadrant.  Can be -1.
 * \param [in] q                Valid quadrant's ancestor is searched.
 * \return                      Offset in the ghost layer, or -1 if not found.
 */
ssize_t             p8est_ghost_tree_contains (p8est_ghost_t * ghost,
                                               int which_proc,
                                               p4est_topidx_t which_tree,
                                               const p8est_quadrant_t * q);

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree boundaries it checks if the quadrant exists
 * across any face, but not across edges or corners.
 *
 * \param [in]  p8est        The forest in which to search for \a q.
 * \param [in]  ghost        The ghost layer in which to search for \a q.
 * \param [in]  treeid       The tree to which \a q belongs.
 * \param [in]  q            The quadrant that is being searched for.
 * \param [in,out] face      On input, face id across which \a q was created.
 *                           On output, the neighbor's face number augmented
 *                           by orientation, so face is in 0..23.
 * \param [in,out] hang      If not NULL, signals that q is bigger than
 *                           the quadrant it came from.  The child id
 *                           of that originating quadrant is passed into hang.
 *                           On output, hang holds the hanging face number
 *                           of \a q that is in contact with its originator.
 * \param [out] owner_rank   Filled with the rank of the owner if it is found
 *                           and undefined otherwise.
 *
 * \return      Returns the local number of \a q if the quadrant exists
 *              in the local forest or in the ghost_layer.  Otherwise,
 *              returns -2 for a domain boundary and -1 if not found.
 */
p4est_locidx_t      p8est_face_quadrant_exists (p8est_t * p8est,
                                                p8est_ghost_t * ghost,
                                                p4est_topidx_t treeid,
                                                const p8est_quadrant_t * q,
                                                int *face, int *hang,
                                                int *owner_rank);

/** Checks if quadrant exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners it checks if the quadrant exists
 * in any of the corner neighbors, thus it can execute multiple queries.
 *
 * \param [in]  p8est        The forest in which to search for \a q
 * \param [in]  ghost        The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree to which \a q belongs (can be extended).
 * \param [in]  q            The quadrant that is being searched for.
 * \param [in,out] exists_arr Must exist and be of of elem_size = sizeof (int)
 *                           for inter-tree corner cases.  Is resized by this
 *                           function to one entry for each corner search
 *                           and set to true/false depending on its existence
 *                           in the local forest or ghost_layer.
 * \param [in,out] rproc_arr If not NULL is filled with one rank per query.
 * \param [in,out] rquad_arr If not NULL is filled with one quadrant per query.
 *                           Its piggy3 member is defined as well.
 *
 * \return true if the quadrant exists in the local forest or in the
 *                  ghost_layer, and false if doesn't exist in either.
 */
int                 p8est_quadrant_exists (p8est_t * p8est,
                                           p8est_ghost_t * ghost,
                                           p4est_topidx_t treeid,
                                           const p8est_quadrant_t * q,
                                           sc_array_t * exists_arr,
                                           sc_array_t * rproc_arr,
                                           sc_array_t * rquad_arr);

/** Check a forest to see if it is balanced.
 *
 * This function builds the ghost layer and discards it when done.
 *
 * \param [in] p8est    The p8est to be tested.
 * \param [in] btype    Balance type (face, edge, corner or default, full).
 * \return Returns true if balanced, false otherwise.
 */
int                 p8est_is_balanced (p8est_t * p8est,
                                       p8est_connect_type_t btype);

/** Compute the parallel checksum of a ghost layer.
 * \param [in] p8est   The MPI information of this p8est will be used.
 * \param [in] ghost   A ghost layer obtained from the p8est.
 * \return             Parallel checksum on rank 0, 0 otherwise.
 */
unsigned            p8est_ghost_checksum (p8est_t * p8est,
                                          p8est_ghost_t * ghost);

/** Transfer data for local quadrants that are ghosts to other processors.
 * Send the data stored in the quadrant's user_data.  This is either the
 * pointer variable itself if \c p8est->data_size is 0, or the content of
 * the referenced memory field if p8est->data_size is positive.
 * \param [in] p8est            The forest used for reference.
 * \param [in] ghost            The ghost layer used for reference.
 * \param [in,out] ghost_data   Pre-allocated contiguous data for all ghost
 *                              quadrants in sequence.  If p8est->data_size is
 *                              0, must at least hold sizeof (void *) bytes for
 *                              each, otherwise p8est->data_size each.
 */
void                p8est_ghost_exchange_data (p8est_t * p8est,
                                               p8est_ghost_t * ghost,
                                               void *ghost_data);

/** Transient storage for asynchronous ghost exchange. */
typedef struct p8est_ghost_exchange
{
  int                 is_custom;        /**< False for p8est_ghost_exchange_data */
  int                 is_levels;        /**< Are we restricted to levels or not */
  p8est_t            *p4est;
  p8est_ghost_t      *ghost;
  int                 minlevel, maxlevel;       /**< Meaningful with is_levels */
  size_t              data_size;
  void               *ghost_data;
  int                *qactive, *qbuffer;
  sc_array_t          requests, sbuffers;
  sc_array_t          rrequests, rbuffers;
}
p8est_ghost_exchange_t;

/** Begin an asynchronous ghost data exchange by posting messages.
 * The arguments are identical to p8est_ghost_exchange_data.
 * The return type is always non-NULL and must be passed to
 * p8est_ghost_exchange_data_end to complete the exchange.
 * The ghost data must not be accessed before completion.
 * \param [in,out]  ghost_data  Must stay alive into the completion call.
 * \return          Transient storage for messages in progress.
 */
p8est_ghost_exchange_t *p8est_ghost_exchange_data_begin
  (p8est_t * p8est, p8est_ghost_t * ghost, void *ghost_data);

/** Complete an asynchronous ghost data exchange.
 * This function waits for all pending MPI communications.
 * \param [in,out]  exc Created ONLY by p8est_ghost_exchange_data_begin.
 *                      It is deallocated before this function returns.
 */
void                p8est_ghost_exchange_data_end
  (p8est_ghost_exchange_t * exc);

/** Transfer data for local quadrants that are ghosts to other processors.
 * The data size is the same for all quadrants and can be chosen arbitrarily.
 * \param [in] p8est            The forest used for reference.
 * \param [in] ghost            The ghost layer used for reference.
 * \param [in] data_size        The data size to transfer per quadrant.
 * \param [in] mirror_data      One data pointer per mirror quadrant. 
 * \param [in,out] ghost_data   Pre-allocated contiguous data for all ghosts
 *                              in sequence, which must hold at least \c
 *                              data_size for each ghost.
 */
void                p8est_ghost_exchange_custom (p8est_t * p8est,
                                                 p8est_ghost_t * ghost,
                                                 size_t data_size,
                                                 void **mirror_data,
                                                 void *ghost_data);

/** Begin an asynchronous ghost data exchange by posting messages.
 * The arguments are identical to p8est_ghost_exchange_custom.
 * The return type is always non-NULL and must be passed to
 * p8est_ghost_exchange_custom_end to complete the exchange.
 * The ghost data must not be accessed before completion.
 * The mirror data can be safely discarded right after this function returns
 * since it is copied into internal send buffers.
 * \param [in]      mirror_data Not required to stay alive any longer.
 * \param [in,out]  ghost_data  Must stay alive into the completion call.
 * \return          Transient storage for messages in progress.
 */
p8est_ghost_exchange_t *p8est_ghost_exchange_custom_begin
  (p8est_t * p8est, p8est_ghost_t * ghost,
   size_t data_size, void **mirror_data, void *ghost_data);

/** Complete an asynchronous ghost data exchange.
 * This function waits for all pending MPI communications.
 * \param [in,out]  Data created ONLY by p8est_ghost_exchange_custom_begin.
 *                  It is deallocated before this function returns.
 */
void                p8est_ghost_exchange_custom_end
  (p8est_ghost_exchange_t * exc);

/** Transfer data for local quadrants that are ghosts to other processors.
 * The data size is the same for all quadrants and can be chosen arbitrarily.
 * This function restricts the transfer to a range of refinement levels.
 * The memory for quadrants outside the level range is not dereferenced.
 * \param [in] p8est            The forest used for reference.
 * \param [in] ghost            The ghost layer used for reference.
 * \param [in] minlevel         Level of the largest quads to be exchanged.
 *                              Use <= 0 for no restriction.
 * \param [in] maxlevel         Level of the smallest quads to be exchanged.
 *                              Use >= P8EST_QMAXLEVEL for no restriction.
 * \param [in] data_size        The data size to transfer per quadrant.
 * \param [in] mirror_data      One data pointer per mirror quadrant as input. 
 * \param [in,out] ghost_data   Pre-allocated contiguous data for all ghosts
 *                              in sequence, which must hold at least \c
 *                              data_size for each ghost.
 */
void                p8est_ghost_exchange_custom_levels (p8est_t * p8est,
                                                        p8est_ghost_t * ghost,
                                                        int minlevel,
                                                        int maxlevel,
                                                        size_t data_size,
                                                        void **mirror_data,
                                                        void *ghost_data);

/** Begin an asynchronous ghost data exchange by posting messages.
 * The arguments are identical to p8est_ghost_exchange_custom_levels.
 * The return type is always non-NULL and must be passed to
 * p8est_ghost_exchange_custom_levels_end to complete the exchange.
 * The ghost data must not be accessed before completion.
 * The mirror data can be safely discarded right after this function returns
 * since it is copied into internal send buffers.
 * \param [in]      mirror_data Not required to stay alive any longer.
 * \param [in,out]  ghost_data  Must stay alive into the completion call.
 * \return          Transient storage for messages in progress.
 */
p8est_ghost_exchange_t *p8est_ghost_exchange_custom_levels_begin
  (p8est_t * p8est, p8est_ghost_t * ghost, int minlevel, int maxlevel,
   size_t data_size, void **mirror_data, void *ghost_data);

/** Complete an asynchronous ghost data exchange.
 * This function waits for all pending MPI communications.
 * \param [in,out]  exc created ONLY by p8est_ghost_exchange_custom_levels_begin.
 *                  It is deallocated before this function returns.
 */
void                p8est_ghost_exchange_custom_levels_end
  (p8est_ghost_exchange_t * exc);

/** Expand the size of the ghost layer and mirrors by one additional layer of
 * adjacency.
 * \param [in] p8est            The forest from which the ghost layer was
 *                              generated.
 * \param [in,out] ghost        The ghost layer to be expanded.
 */
void                p8est_ghost_expand (p8est_t * p8est,
                                        p8est_ghost_t * ghost);

SC_EXTERN_C_END;

#endif /* !P8EST_GHOST_H */
