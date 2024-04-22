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

#ifndef P6EST_LNODES_H
#define P6EST_LNODES_H

#include <p6est_ghost.h>
#include <p4est_lnodes.h>
#include <p8est_lnodes.h>

SC_EXTERN_C_BEGIN;

/* A p6est_lnodes_t is exactly the same as a p8est_lnodes_t, with the only
 * difference being that the face_codes are interpreted differently to account
 * for the types of hanging faces that occur in a p6est.  Please see the
 * documentation for p8est_lnodes */

/* The only other difference is in the numbering of nodes, edges, and faces.
 *
 * Columns of nodes are numbered contiguously: this still generates a
 * partition-unique numbering.
 *
 * Although we call a p2est_quadrant_t coordinate layer->z, the orientation of
 * a layer from lnodes perspective is that the vertical axis is the X axis of
 * the 3D element, the x axis of the columns is the Y axis of the 3D element,
 * and the y axis of the columns is the Z axis of the 3D element
 */

typedef p8est_lnodes_t p6est_lnodes_t;
typedef p8est_lnodes_code_t p6est_lnodes_code_t;
typedef p8est_lnodes_rank_t p6est_lnodes_rank_t;
typedef p8est_lnodes_buffer_t p6est_lnodes_buffer_t;

/** Decode the face_code into hanging face information.
 *
 * \param[in] face_code as in the p6est_lnodes_t structure.
 * \param[out] hanging_face: if there are hanging faces or edges,
 *             hanging_face = -1 if the face is not hanging,
 *                          = the corner of the full face that it touches:
 *                            e.g. if face = i and hanging_face[i] =
 *                            j, then the interpolation operator corresponding
 *                            to corner j should be used for that face.
 *             note: not touched if there are no hanging faces or edges.
 * \param[out] hanging_edge: if there are hanging faces or edges,
 *             hanging_edge = -1 if the edge is not hanging,
 *                          =  0 if the edge is the first half of a full edge,
 *                               but neither of the two faces touching the
 *                               edge is hanging,
 *                          =  1 if the edge is the second half of a full edge,
 *                               but neither of the two faces touching the
 *                               edge is hanging,
 *                          =  2 if the edge is the first half of a full edge
 *                               and is on the boundary of a full face,
 *                          =  3 if the edge is the second half of a full edge
 *                               and is on the boundary of a full face,
 *                          =  4 if the edge is in the middle of a full face.
 *                               See the diagram below for clarification.
 *             note: not touched if there are no hanging faces or edges.
 * \return             true if any face or edge is hanging, false otherwise.
 *
 * o...............o  o...............o  +---2---+.......o  o.......+---3---+
 * :               :  :               :  |       |       :  :       |       |
 * :               :  :               :  3   2   4       :  :       4   3   3
 * :               :  :               :  |       |       :  :       |       |
 * +---4---+       :  :       +---4---+  +---4---+       :  :       +---4---+
 * |       |       :  :       |       |  :               :  :               :
 * 2   0   4       :  :       4   1   2  :               :  :               :
 * |       |       :  :       |       |  :               :  :               :
 * +---2---+.......o  o.......+---3---+  o...............o  o...............o
 *
 * o...............o  +-----(-1)------+  +---2---+.......o  o.......+---3---+
 * :               :  |               |  |       |       :  :       |       |
 * :               :  3       5       3  |       |       :  :       |       |
 * :               :  |               |  |       |       :  :       |       |
 * +-------4-------+  +-------4-------+ -1   6   4       :  :       4   7  -1
 * |               |  :               :  |       |       :  :       |       |
 * 2       4       2  :               :  |       |       :  :       |       |
 * |               |  :               :  |       |       :  :       |       |
 * +-----(-1)------+  o...............o  +---2---+.......o  o.......+---3---+
 *
 *                    o                  +-------+
 *                    :                  |\       \
 *                    :                  1 \       \
 *                    :                  |  +-------+
 *                    +-------+          +  |       |
 *                    |\       \         :\ |       |
 *                    0 \       \        : \|       |
 *                    |  +-------+       :  +-------+
 *                    +  |       |       o
 *                     \ |       |
 *                      \|       |
 *                       +-------+
 */
/*@unused@*/
static inline int
p6est_lnodes_decode (p6est_lnodes_code_t face_code, int hanging_face[6],
                     int hanging_edge[12])
{
  P4EST_ASSERT (face_code >= 0);

  if (face_code) {
    /* we pack the p4est_lnodes_code_t at the bottom, followed by a bit
     * indicating whether this layer is a first or second sibling, followed by
     * four bits indicating which of the four side faces are layerwise
     * nonconforming, followed by four bits indicating which of the four side
     * edges are layerwise nonconforming */
    p4est_lnodes_code_t fc4 = face_code & 0x000f;
    int16_t             h = (face_code & 0x0010) >> 4;
    int16_t             work = face_code >> 5;
    int                 hf;
    int                 f, e, w;

    memset (hanging_face, -1, 6 * sizeof (int));
    memset (hanging_edge, -1, 12 * sizeof (int));

    /* the first two faces are the top and bottom faces, which we know are not
     * hanging */
    p4est_lnodes_decode (fc4, hanging_face + 2);
    for (f = 0; f < 4; f++) {
      hf = hanging_face[f + 2];
      w = work & 0x0001;
      if (hf >= 0) {
        hanging_edge[p8est_face_edges[f + 2][2]] = 2 + hf;
        hanging_edge[p8est_face_edges[f + 2][3]] = 2 + hf;
        hanging_edge[p8est_face_edges[f + 2][1 ^ hf]] = 4;
        if (w) {
          hanging_edge[p8est_face_edges[f + 2][3 ^ h]] = 4;
          hanging_edge[p8est_face_edges[f + 2][1 ^ hf]] = 4;
          hanging_edge[p8est_face_edges[f + 2][hf]] = 2 + h;
          hanging_face[f + 2] = (hf << 1) | h;
        }
        else {
          hanging_face[f + 2] = 4 + hf;
        }
      }
      else if (w) {
        hanging_edge[p8est_face_edges[f + 2][3 ^ h]] = 4;
        hanging_edge[p8est_face_edges[f + 2][0]] =
          SC_MAX (hanging_edge[p8est_face_edges[f + 2][0]], 2 + h);
        hanging_edge[p8est_face_edges[f + 2][1]] =
          SC_MAX (hanging_edge[p8est_face_edges[f + 2][1]], 2 + h);
        hanging_face[f + 2] = 6 + h;
      }
      work >>= 1;
    }
    for (e = 0; e < 4; e++) {
      if (work & 0x0001) {
        if (hanging_edge[e] < 0) {
          hanging_edge[e] = h;
        }
#ifdef P4EST_ENABLE_DEBUG
        else {
          P4EST_ASSERT (hanging_edge[e] == 2 + h || hanging_edge[e] == 4);
        }
#endif
      }
      work >>= 1;
    }
    return 1;
  }
  else {
    return 0;
  }
}

p6est_lnodes_t     *p6est_lnodes_new (p6est_t * p6est,
                                      p6est_ghost_t * ghost_layer,
                                      int degree);

static inline void
p6est_lnodes_destroy (p6est_lnodes_t * lnodes)
{
  p8est_lnodes_destroy (lnodes);
}

/*@unused@*/
static inline p6est_lnodes_buffer_t *
p6est_lnodes_share_owned_begin (sc_array_t * node_data,
                                p6est_lnodes_t * lnodes)
{
  return p8est_lnodes_share_owned_begin (node_data, lnodes);
}

/*@unused@*/
static inline void
p6est_lnodes_share_owned_end (p6est_lnodes_buffer_t * buffer)
{
  p8est_lnodes_share_owned_end (buffer);
}

/*@unused@*/
static inline void
p6est_lnodes_share_owned (sc_array_t * node_data, p6est_lnodes_t * lnodes)
{
  p8est_lnodes_share_owned (node_data, lnodes);
}

/*@unused@*/
static inline p6est_lnodes_buffer_t *
p6est_lnodes_share_all_begin (sc_array_t * node_data, p6est_lnodes_t * lnodes)
{
  return p8est_lnodes_share_all_begin (node_data, lnodes);
}

/*@unused@*/
static inline void
p6est_lnodes_share_all_end (p6est_lnodes_buffer_t * buffer)
{
  p8est_lnodes_share_all_end (buffer);
}

/*@unused@*/
static inline p6est_lnodes_buffer_t *
p6est_lnodes_share_all (sc_array_t * node_data, p6est_lnodes_t * lnodes)
{
  return p8est_lnodes_share_all (node_data, lnodes);
}

/*@unused@*/
static inline void
p6est_lnodes_buffer_destroy (p6est_lnodes_buffer_t * buffer)
{
  p8est_lnodes_buffer_destroy (buffer);
}

/*@unused@*/
static inline p6est_lnodes_rank_t *
p6est_lnodes_rank_array_index_int (sc_array_t * array, int it)
{
  return p8est_lnodes_rank_array_index_int (array, it);
}

/*@unused@*/
static inline p6est_lnodes_rank_t *
p6est_lnodes_rank_array_index (sc_array_t * array, size_t it)
{
  return p8est_lnodes_rank_array_index (array, it);
}

/*@unused@*/
static inline       p4est_gloidx_t
p6est_lnodes_global_index (p6est_lnodes_t * lnodes, p4est_locidx_t lidx)
{
  return p8est_lnodes_global_index (lnodes, lidx);
}

/** For each owned node, get the global 2D number for the node-column
 * containing it.
 *
 * \param[in] p6est  The forest
 * \param[in] lnodes The nodes
 * \return an array of size \a lnodes->owned_count, giving the unique global
 * number of the node-column containing each node.  Should be free'd with
 * P4EST_FREE().
 */
p4est_gloidx_t     *p6est_lnodes_get_column_labels (p6est_t * p6est,
                                                    p8est_lnodes_t * lnodes);

SC_EXTERN_C_END;

#endif /* !P6EST_LNODES */
