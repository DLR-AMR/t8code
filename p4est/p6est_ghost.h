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

#ifndef P6EST_GHOST_H
#define P6EST_GHOST_H

/** \file p6est_ghost.h
 *
 * passing columns of layers and data to neighboring processes
 *
 * \ingroup p6est
 */

#include <p6est.h>
#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

/** columns of layers that neighbor the local domain */
typedef struct p6est_ghost
{
  int                 mpisize;
  p4est_topidx_t      num_trees;
  p4est_connect_type_t btype; /**< which neighboring columns are in the ghost layer */

  p4est_ghost_t      *column_ghost; /**< describes the ghost columns */
  sc_array_t         *column_layer_offsets; /**< array of p4est_locidx_t type:
                                              the offset of each ghost columns
                                              within the \a ghosts array of
                                              column-layers */

  /** An array of column-layers which make up the ghost layer around \a
   * p6est.  Their piggy3 data member is filled with their owner's tree
   * and local number (cumulative over trees).  Quadrants are ordered in \c
   * p4est_quadrant_compare_piggy order.  These are quadrants inside the
   * neighboring tree, i.e., \c p4est_quadrant_is_inside is true for the
   * quadrant and the neighboring tree.
   */
  sc_array_t          ghosts; /**< array of p2est_quadrant_t type */
  p4est_locidx_t     *tree_offsets;     /**< num_trees + 1 ghost indices */
  p4est_locidx_t     *proc_offsets;     /**< mpisize + 1 ghost indices */

  /** An array of local quadrants that touch the parallel boundary from the
   * inside, i.e., that are ghosts in the perspective of at least one other
   * processor.  The storage convention is the same as for \c ghosts above.
   */
  sc_array_t          mirrors; /**< array of p4est_quadrant_t type */
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
                                                   p4est_ghost_expand is
                                                   called */
  p4est_locidx_t     *mirror_proc_front_offsets;        /**< NULL until
                                                           p4est_ghost_expand is
                                                           called */

}
p6est_ghost_t;

/** Calculate the memory usage of the ghost layer.
 * \param [in] ghost    Ghost layer structure.
 * \return              Memory used in bytes.
 */
size_t              p6est_ghost_memory_used (p6est_ghost_t * ghost);

/** Builds the ghost layer.
 *
 * This will gather the quadrants from each neighboring proc to build
 * one layer of face and corner based ghost elements around the ones they own.
 *
 * \param [in] p4est            The forest for which the ghost layer will be
 *                              generated.
 * \param [in] btype            Which ghosts to include (across face, corner
 *                              or default, full).
 * \return                      A fully initialized ghost layer.
 */
p6est_ghost_t      *p6est_ghost_new (p6est_t * p4est,
                                     p4est_connect_type_t btype);

/** Expand the size of the ghost layer and mirrors by one additional layer of
 * adjacency.
 * \param [in] p6est            The forest from which the ghost layer was
 *                              generated.
 * \param [in,out] ghost        The ghost layer to be expanded.
 */
void                p6est_ghost_expand (p6est_t * p6est,
                                        p6est_ghost_t * ghost);

/** Frees all memory used for the ghost layer. */
void                p6est_ghost_destroy (p6est_ghost_t * ghost);

/** Conduct binary search for exact match on a range of the ghost layer.
 * \param [in] ghost            The ghost layer.
 * \param [in] which_proc       The owner of the searched quadrant.  Can be -1.
 * \param [in] which_tree       The tree of the searched quadrant.  Can be -1.
 * \param [in] q                Valid quadrant is searched in the ghost layer.
 * \return                      Offset in the ghost layer, or -1 if not found.
 */
ssize_t             p6est_ghost_bsearch (p6est_ghost_t * ghost,
                                         int which_proc,
                                         p4est_topidx_t which_tree,
                                         const p4est_quadrant_t * column,
                                         const p2est_quadrant_t * layer);

/** Conduct binary search for ancestor on range of the ghost layer.
 * \param [in] ghost            The ghost layer.
 * \param [in] which_proc       The owner of the searched quadrant.  Can be -1.
 * \param [in] which_tree       The tree of the searched quadrant.  Can be -1.
 * \param [in] q                Valid quadrant's ancestor is searched.
 * \return                      Offset in the ghost layer, or -1 if not found.
 */
ssize_t             p6est_ghost_contains (p6est_ghost_t * ghost,
                                          int which_proc,
                                          p4est_topidx_t which_tree,
                                          const p4est_quadrant_t * column,
                                          const p2est_quadrant_t * layer);

/** Checks if layer exists in the local forest or the ghost layer.
 *
 * For quadrants across tree corners it checks if the quadrant exists
 * in any of the corner neighbors, thus it can execute multiple queries.
 *
 * \param [in]  p4est        The forest in which to search for \a q
 * \param [in]  ghost        The ghost layer in which to search for \a q
 * \param [in]  treeid       The tree to which \a q belongs (can be extended).
 * \param [in]  column       The column that is being searched for.
 * \param [in]  layer        The layer that is being searched for.
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
int                 p6est_layer_exists (p6est_t * p6est,
                                        p6est_ghost_t * ghost,
                                        p4est_topidx_t treeid,
                                        const p4est_quadrant_t * column,
                                        const p2est_quadrant_t * layer,
                                        sc_array_t * exists_arr,
                                        sc_array_t * rproc_arr,
                                        sc_array_t * rquad_arr);

/** Check a forest to see if it is balanced.
 *
 * This function builds the ghost layer and discards it when done.
 *
 * \param [in] p4est    The p4est to be tested.
 * \param [in] btype    Balance type (face, corner or default, full).
 * \return Returns true if balanced, false otherwise.
 */
int                 p6est_is_balanced (p6est_t * p6est,
                                       p8est_connect_type_t btype);

/** Compute the parallel checksum of a ghost layer.
 * \param [in] p4est   The MPI information of this p4est will be used.
 * \param [in] ghost   A ghost layer obtained from the p4est.
 * \return             Parallel checksum on rank 0, 0 otherwise.
 */
unsigned            p6est_ghost_checksum (p6est_t * p6est,
                                          p6est_ghost_t * ghost);

SC_EXTERN_C_END;

#endif /* P6EST_GHOST_H */
