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

#ifndef P4EST_LNODES_H
#define P4EST_LNODES_H

#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

typedef int8_t      p4est_lnodes_code_t;

/** Store a parallel numbering of Lobatto points of a given degree > 0.
 *
 * Each element has degree+1 nodes per face
 * and vnodes = (degree+1)^2 nodes per volume.
 * num_local_elements is the number of local quadrants in the p4est.
 * element_nodes is of dimension vnodes * num_local_elements and lists the
 * nodes of each element in lexicographic yx-order (x varies fastest); so for
 * degree == 2, this is the layout of nodes:
 *
 *         f_3
 *  c_2           c_3
 *      6---7---8
 *      |       |
 * f_0  3   4   5  f_1
 *      |       |
 *      0---1---2
 *  c_0           c_1
 *         f_2
 *
 * element_nodes indexes into the set of local nodes, layed out as follows:
 * local nodes = [<-----owned_count----->|<-----nonlocal_nodes----->]
 *             = [<----------------num_local_nodes----------------->]
 * nonlocal_nodes contains the globally unique numbers for independent nodes
 * that are owned by other processes; for local nodes, the globally unique
 * numbers are given by i + global_offset, where i is the local number.
 * Hanging nodes are always local and don't have a global number.
 * They index the geometrically corresponding independent nodes of a neighbor.
 *
 * Whether nodes are hanging or not is decided based on the element faces.
 * This information is encoded in face_code with one int8_t per element.
 * If no faces are hanging, the value is zero, otherwise the face_code is
 * interpreted by p4est_lnodes_decode.
 *
 * Independent nodes can be shared by multiple MPI ranks.
 * The owner rank of a node is the one from the lowest numbered element
 * on the lowest numbered octree *touching* the node.
 *
 * What is meant by *touching*?
 * A quadrant is said to touch all faces/corners that are incident on it,
 * and by extension all nodes that are contained in those faces/corners.
 *
 *         X +-----------+
 *         o |           |
 *         o |           |
 * +-----+ o |     p     |
 * |  q  | o |           |
 * |     | o |           |
 * +-----+ O +-----------+
 *
 * In this example degree = 6.  There are 5 nodes that live on the face
 * between q and p, and one at each corner of that face.  The face is incident
 * on q, so q owns the nodes on the face (provided q is from a lower tree or
 * has a lower index than p).  The lower corner is incident on q, so q owns it
 * as well.  The upper corner is not incident on q, so q cannot own it.
 *
 * global_owned_count contains the number of independent nodes owned by each
 * process.
 *
 * The sharers array contains items of type p4est_lnodes_rank_t
 * that hold the ranks that own or share independent local nodes.
 * If there are no shared nodes on this processor, it is empty.
 * Otherwise, it is sorted by rank and the current process is included.
 *
 * degree < 0 indicates that the lnodes data structure is being used to number
 * the quadrant boundary object (faces and corners) rather than the $C^0$
 * Lobatto nodes:
 *
 * if degree == -1, then one node is assigned per face, and no nodes are
 * assigned per volume or per corner: this numbering can be used for low-order
 * Raviart-Thomas elements.  In this case, vnodes == 4, and the nodes are
 * listed in face-order:
 *
 *         f_3
 *  c_2           c_3
 *      +---3---+
 *      |       |
 * f_0  0       1  f_1
 *      |       |
 *      +---2---+
 *  c_0           c_1
 *         f_2
 *
 * if degree == -2, then one node is assigned per face and per corner and no
 * nodes are assigned per volume.  In this case, vnodes == 8, and the
 * nodes are listed in face-order, followed by corner-order:
 *
 *         f_3
 *  c_2           c_3
 *      6---3---7
 *      |       |
 * f_0  0       1  f_1
 *      |       |
 *      4---2---5
 *  c_0           c_1
 *         f_2
 *
 */
typedef struct p4est_lnodes
{
  sc_MPI_Comm         mpicomm;
  p4est_locidx_t      num_local_nodes;
  p4est_locidx_t      owned_count;
  p4est_gloidx_t      global_offset;
  p4est_gloidx_t     *nonlocal_nodes;
  sc_array_t         *sharers;
  p4est_locidx_t     *global_owned_count;

  int                 degree, vnodes;
  p4est_locidx_t      num_local_elements;
  p4est_lnodes_code_t *face_code;
  p4est_locidx_t     *element_nodes;
}
p4est_lnodes_t;

/** The structure stored in the sharers array.
 *
 * shared_nodes is a sorted array of p4est_locidx_t
 * that indexes into local nodes.  The shared_nodes array has a
 * contiguous (or empty) section of nodes owned by the current rank.
 * shared_mine_offset and shared_mine_count identify this section
 * by indexing the shared_nodes array, not the local nodes array.
 * owned_offset and owned_count define the section of local nodes
 * that is owned by the listed rank (the section may be empty).
 * For the current process these coincide with those in p4est_lnodes_t.
 */
typedef struct p4est_lnodes_rank
{
  int                 rank;
  sc_array_t          shared_nodes;
  p4est_locidx_t      shared_mine_offset, shared_mine_count;
  p4est_locidx_t      owned_offset, owned_count;
}
p4est_lnodes_rank_t;

/** Decode the face_code into hanging face information.
 *
 * This is mostly for demonstration purposes.  Applications probably will
 * integrate it into their own loop over the face for performance reasons.
 *
 * \param[in] face_code as in the p4est_lnodes_t structure.
 * \param[out] hanging face: if there are hanging faces,
 *             hanging_face = -1 if the face is not hanging,
 *                          = 0 if the face is the first half,
 *                          = 1 if the face is the second half.
 *             note: not touched if there are no hanging faces.
 * \return              true if any face is hanging, false otherwise.
 */
/*@unused@*/
static inline int
p4est_lnodes_decode (p4est_lnodes_code_t face_code, int hanging_face[4])
{
  P4EST_ASSERT (face_code >= 0);

  if (face_code) {
    int                 i;
    int8_t              c = face_code & 0x03;
    int                 f;
    int8_t              work = face_code >> 2;

    memset (hanging_face, -1, 4 * sizeof (int));

    for (i = 0; i < 2; ++i) {
      f = p4est_corner_faces[c][i];
      hanging_face[f] = (work & 0x01) ? p4est_corner_face_corners[c][f] : -1;
      work >>= 1;
    }

    return 1;
  }
  else {
    return 0;
  }
}

p4est_lnodes_t     *p4est_lnodes_new (p4est_t * p4est,
                                      p4est_ghost_t * ghost_layer,
                                      int degree);

void                p4est_lnodes_destroy (p4est_lnodes_t * lnodes);

/** Expand the ghost layer to include the support of all nodes supported on
 * the local partition.
 *
 * \param [in]     p4est        The forest from which the ghost layer was
 *                              generated.
 * \param [in]     lnodes       The nodes to support.
 * \param [in,out] ghost        The ghost layer to be expanded.
 */
void                p4est_ghost_support_lnodes (p4est_t * p4est,
                                                p4est_lnodes_t * lnodes,
                                                p4est_ghost_t * ghost);

/** Expand the ghost layer as in p4est_ghost_expand(), but use node support to
 * define adjacency instead of geometric adjacency.
 *
 * \param [in]     p4est        The forest from which the ghost layer was
 *                              generated.
 * \param [in]     lnodes       The nodes to support.
 * \param [in,out] ghost        The ghost layer to be expanded.
 */
void                p4est_ghost_expand_by_lnodes (p4est_t * p4est,
                                                  p4est_lnodes_t * lnodes,
                                                  p4est_ghost_t * ghost);

/** Partition using weights based on the number of nodes assigned to each
 * element in lnodes
 *
 * \param[in,out] p4est                    the forest to be repartitioned
 * \param[in]     ghost                    the ghost layer
 * \param[in]     degree                   the degree that would be passed to p4est_lnodes_new()
 * \param[in]     partition_for_coarsening whether the partition should allow
 *                                         coarsening (i.e. group siblings who
 *                                         might merge)
 */
void                p4est_partition_lnodes (p4est_t * p4est,
                                            p4est_ghost_t * ghost, int degree,
                                            int partition_for_coarsening);

/** Partition using weights that are broken down by where they reside: in
 * volumes, on faces, or on corners.
 */
void                p4est_partition_lnodes_detailed (p4est_t * p4est,
                                                     p4est_ghost_t * ghost,
                                                     int nodes_per_volume,
                                                     int nodes_per_face,
                                                     int nodes_per_corner,
                                                     int
                                                     partition_for_coarsening);

/** p4est_lnodes_buffer_t handles the communication of data associated with
 * nodes.
 *
 * \a send_buffers is an array of arrays: one buffer for each process to which
 * the current process sends node-data.  It should not be altered between
 * a shared_*_begin and a shared_*_end call.
 *
 * \a recv_buffers is an array of arrays that is used in lnodes_share_all_*.
 * \a recv_buffers[j] corresponds with lnodes->sharers[j]: it is the same
 * length as \a lnodes->sharers[j]->shared_nodes.  At the completion of
 * lnodes_share_all or lnodes_share_all_end, recv_buffers[j] contains the
 * node-data from the process lnodes->sharers[j]->rank
 * (unless j is the current rank, in which case recv_buffers[j] is empty).
 */
typedef struct p4est_lnodes_buffer
{
  sc_array_t         *requests; /* sc_MPI_Request */
  sc_array_t         *send_buffers;
  sc_array_t         *recv_buffers;
}
p4est_lnodes_buffer_t;

/** p4est_lnodes_share_owned_begin
 *
 * \a node_data is a user-defined array of arbitrary type, where each entry
 * is associated with the \a lnodes local nodes entry of matching index.
 * For every local nodes entry that is owned by a process
 * other than the current one, the value in the \a node_data array of the
 * owning process is written directly into the \a node_data array of the current
 * process.  Values of \a node_data are not guaranteed to be sent or received
 * until the \a buffer created by p4est_lnodes_share_owned_begin is passed to
 * p4est_lnodes_share_owned_end.
 *
 * To be memory neutral, the \a buffer created by
 * p4est_lnodes_share_owned_begin must be destroying with
 * p4est_lnodes_buffer_destroy (it is not destroyed by
 * p4est_lnodes_share_owned_end).
 */
p4est_lnodes_buffer_t *p4est_lnodes_share_owned_begin (sc_array_t * node_data,
                                                       p4est_lnodes_t *
                                                       lnodes);

void                p4est_lnodes_share_owned_end (p4est_lnodes_buffer_t *
                                                  buffer);

/** Equivalent to calling p4est_lnodes_share_owned_end directly after
 * p4est_lnodes_share_owned_begin.  Use if there is no local work that can be
 * done to mask the communication cost.
 */
void                p4est_lnodes_share_owned (sc_array_t * node_data,
                                              p4est_lnodes_t * lnodes);

/** p4est_lnodes_share_all_begin
 *
 * \a node_data is a user_defined array of arbitrary type, where each entry
 * is associated with the lnodes local nodes entry of matching index.
 * For every process that shares an entry with the current one, the value in
 * the \a node_data array of that process is written into a
 * \a buffer->recv_buffers entry as described above.  The user can then perform
 * some arbitrary work that requires the data from all processes that share a
 * node (such as reduce, max, min, etc.).  When the work concludes, the
 * \a buffer should be destroyed with p4est_lnodes_buffer_destroy.
 *
 * Values of \a node_data are not guaranteed to be sent, and
 * \a buffer->recv_buffer entries are not guaranteed to be received until
 * the \a buffer created by p4est_lnodes_share_all_begin is passed to
 * p4est_lnodes_share_all_end.
 */
p4est_lnodes_buffer_t *p4est_lnodes_share_all_begin (sc_array_t * node_data,
                                                     p4est_lnodes_t * lnodes);

void                p4est_lnodes_share_all_end (p4est_lnodes_buffer_t *
                                                buffer);

/** Equivalent to calling p4est_lnodes_share_all_end directly after
 * p4est_lnodes_share_all_begin.  Use if there is no local work that can be
 * done to mask the communication cost.
 * \return          A fully initialized buffer that contains the received data.
 *                  After processing this data, the buffer must be freed with
 *                  p4est_lnodes_buffer_destroy.
 */
p4est_lnodes_buffer_t *p4est_lnodes_share_all (sc_array_t * node_data,
                                               p4est_lnodes_t * lnodes);

void                p4est_lnodes_buffer_destroy (p4est_lnodes_buffer_t *
                                                 buffer);

/** Return a pointer to a lnodes_rank array element indexed by a int.
 */
/*@unused@*/
static inline p4est_lnodes_rank_t *
p4est_lnodes_rank_array_index_int (sc_array_t * array, int it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_lnodes_rank_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p4est_lnodes_rank_t *)
    (array->array + sizeof (p4est_lnodes_rank_t) * (size_t) it);
}

/** Return a pointer to a lnodes_rank array element indexed by a size_t.
 */
/*@unused@*/
static inline p4est_lnodes_rank_t *
p4est_lnodes_rank_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_lnodes_rank_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p4est_lnodes_rank_t *)
    (array->array + sizeof (p4est_lnodes_rank_t) * it);
}

/** Compute the global number of a local node number */
/*@unused@*/
static inline       p4est_gloidx_t
p4est_lnodes_global_index (p4est_lnodes_t * lnodes, p4est_locidx_t lidx)
{
  p4est_locidx_t      owned = lnodes->owned_count;
  P4EST_ASSERT (lidx >= 0 && lidx < lnodes->num_local_nodes);

  return (lidx < owned) ? lnodes->global_offset + lidx :
    lnodes->nonlocal_nodes[lidx - owned];
}

SC_EXTERN_C_END;

#endif /* !P4EST_LNODES */
