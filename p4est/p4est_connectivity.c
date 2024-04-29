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

#ifndef P4_TO_P8
#include <p4est_connectivity.h>
#include <p4est.h>
#endif
#ifdef P4EST_WITH_METIS
#include <metis.h>
#endif

#ifndef P4_TO_P8

/* *INDENT-OFF* */
const int           p4est_face_corners[4][2] =
{{ 0, 2 },
 { 1, 3 },
 { 0, 1 },
 { 2, 3 }};
const int           p4est_face_dual[4] = { 1, 0, 3, 2 };

const int           p4est_corner_faces[4][2] =
{{ 0, 2 },
 { 1, 2 },
 { 0, 3 },
 { 1, 3 }};
const int           p4est_corner_face_corners[4][4] =
{{  0, -1,  0, -1 },
 { -1,  0,  1, -1 },
 {  1, -1, -1,  0 },
 { -1,  1, -1,  1 }};

const int           p4est_child_corner_faces[4][4] =
{{ -1,  2,  0, -1 },
 {  2, -1, -1,  1 },
 {  0, -1, -1,  3 },
 { -1,  1,  3, -1 }};
/* *INDENT-ON* */

#endif /* !P4_TO_P8 */

int
p4est_connectivity_face_neighbor_face_corner (int fc, int f, int nf, int o)
{
  int                 nfc;
#ifdef P4_TO_P8
  int                 pref, pset;
#endif

  /* sanity checks */
  P4EST_ASSERT (0 <= fc && fc < P4EST_HALF);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= o && o < P4EST_HALF);

#ifndef P4_TO_P8
  nfc = fc ^ o;
#else
  pref = p8est_face_permutation_refs[f][nf];
  pset = p8est_face_permutation_sets[pref][o];
  nfc = p8est_face_permutations[pset][fc];
#endif
  P4EST_ASSERT (0 <= nfc && nfc < P4EST_HALF);

  return nfc;
}

int
p4est_connectivity_face_neighbor_corner (int c, int f, int nf, int o)
{
  int                 fc, nfc;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= o && o < P4EST_HALF);

  fc = p4est_corner_face_corners[c][f];
  P4EST_ASSERT (0 <= fc && fc < P4EST_HALF);

  nfc = p4est_connectivity_face_neighbor_face_corner (fc, f, nf, o);

  P4EST_ASSERT (0 <= nfc && nfc < P4EST_HALF);

  return p4est_face_corners[nf][nfc];
}

size_t
p4est_connectivity_memory_used (p4est_connectivity_t * conn)
{
  return sizeof (p4est_connectivity_t) +
    (conn->num_vertices > 0 ?
     (conn->num_vertices * 3 * sizeof (double) +
      conn->num_trees * P4EST_CHILDREN * sizeof (p4est_topidx_t)) : 0) +
    conn->num_trees * P4EST_FACES * (sizeof (p4est_topidx_t) +
                                     sizeof (int8_t)) +
#ifdef P4_TO_P8
    conn->num_trees * P8EST_EDGES * sizeof (p4est_topidx_t) +
    (conn->num_edges + 1) * sizeof (p4est_topidx_t) +
    conn->ett_offset[conn->num_edges] * (sizeof (p4est_topidx_t) +
                                         sizeof (int8_t)) +
#endif
    conn->num_trees * P4EST_CHILDREN * sizeof (p4est_topidx_t) +
    (conn->num_corners + 1) * sizeof (p4est_topidx_t) +
    conn->ctt_offset[conn->num_corners] * (sizeof (p4est_topidx_t) +
                                           sizeof (int8_t));
}

p4est_connectivity_t *
p4est_connectivity_new_copy (p4est_topidx_t num_vertices,
                             p4est_topidx_t num_trees,
#ifdef P4_TO_P8
                             p4est_topidx_t num_edges,
#endif
                             p4est_topidx_t num_corners,
                             const double *vertices,
                             const p4est_topidx_t * ttv,
                             const p4est_topidx_t * ttt, const int8_t * ttf,
#ifdef P4_TO_P8
                             const p4est_topidx_t * tte,
                             const p4est_topidx_t * eoff,
                             const p4est_topidx_t * ett, const int8_t * ete,
#endif
                             const p4est_topidx_t * ttc,
                             const p4est_topidx_t * coff,
                             const p4est_topidx_t * ctt, const int8_t * ctc)
{
#ifdef P4_TO_P8
  p4est_topidx_t      num_ett;
#endif
  p4est_topidx_t      num_ctt;
  p4est_connectivity_t *conn;

  P4EST_ASSERT (num_vertices >= 0);
  P4EST_ASSERT (num_trees >= 0);

#ifdef P4_TO_P8
  P4EST_ASSERT (num_edges >= 0);
  P4EST_ASSERT (eoff != NULL);
  num_ett = eoff[num_edges];
#endif
  P4EST_ASSERT (num_corners >= 0);
  P4EST_ASSERT (coff != NULL);
  num_ctt = coff[num_corners];
  conn = p4est_connectivity_new (num_vertices, num_trees,
#ifdef P4_TO_P8
                                 num_edges, num_ett,
#endif
                                 num_corners, num_ctt);

  if (num_vertices > 0) {
    P4EST_ASSERT (vertices != NULL && ttv != NULL);
    memcpy (conn->vertices, vertices, sizeof (double) * 3 * num_vertices);
    memcpy (conn->tree_to_vertex, ttv,
            sizeof (p4est_topidx_t) * P4EST_CHILDREN * num_trees);
  }
  else {
    conn->vertices = NULL;
    conn->tree_to_vertex = NULL;
  }
  memcpy (conn->tree_to_tree, ttt,
          sizeof (p4est_topidx_t) * P4EST_FACES * num_trees);
  memcpy (conn->tree_to_face, ttf, sizeof (int8_t) * P4EST_FACES * num_trees);

#ifdef P4_TO_P8
  if (num_edges > 0) {
    memcpy (conn->tree_to_edge, tte,
            sizeof (p4est_topidx_t) * P8EST_EDGES * num_trees);
    memcpy (conn->edge_to_tree, ett, sizeof (p4est_topidx_t) * num_ett);
    memcpy (conn->edge_to_edge, ete, sizeof (int8_t) * num_ett);
  }
  memcpy (conn->ett_offset, eoff, sizeof (p4est_topidx_t) * (num_edges + 1));
#endif

  if (num_corners > 0) {
    memcpy (conn->tree_to_corner, ttc,
            sizeof (p4est_topidx_t) * P4EST_CHILDREN * num_trees);
    memcpy (conn->corner_to_tree, ctt, sizeof (p4est_topidx_t) * num_ctt);
    memcpy (conn->corner_to_corner, ctc, sizeof (int8_t) * num_ctt);
  }
  memcpy (conn->ctt_offset, coff,
          sizeof (p4est_topidx_t) * (num_corners + 1));

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_new (p4est_topidx_t num_vertices, p4est_topidx_t num_trees,
#ifdef P4_TO_P8
                        p4est_topidx_t num_edges, p4est_topidx_t num_ett,
#endif
                        p4est_topidx_t num_corners, p4est_topidx_t num_ctt)
{
  p4est_connectivity_t *conn;

  P4EST_ASSERT (num_vertices >= 0);
  P4EST_ASSERT (num_trees >= 0);
#ifdef P4_TO_P8
  P4EST_ASSERT (num_edges >= 0);
  P4EST_ASSERT (num_ett >= 0);
#endif
  P4EST_ASSERT (num_corners >= 0);
  P4EST_ASSERT (num_ctt >= 0);

  conn = P4EST_ALLOC_ZERO (p4est_connectivity_t, 1);

  conn->num_vertices = num_vertices;
  conn->num_trees = num_trees;
  if (num_vertices > 0) {
    conn->vertices = P4EST_ALLOC (double, 3 * num_vertices);
    conn->tree_to_vertex =
      P4EST_ALLOC (p4est_topidx_t, P4EST_CHILDREN * num_trees);
  }
  else {
    conn->vertices = NULL;
    conn->tree_to_vertex = NULL;
  }
  conn->tree_to_tree = P4EST_ALLOC (p4est_topidx_t, P4EST_FACES * num_trees);
  conn->tree_to_face = P4EST_ALLOC (int8_t, P4EST_FACES * num_trees);

#ifdef P4_TO_P8
  conn->num_edges = num_edges;
  if (num_edges > 0) {
    conn->tree_to_edge =
      P4EST_ALLOC (p4est_topidx_t, P8EST_EDGES * num_trees);
    conn->edge_to_tree = P4EST_ALLOC (p4est_topidx_t, num_ett);
    conn->edge_to_edge = P4EST_ALLOC (int8_t, num_ett);
  }
  else {
    conn->tree_to_edge = NULL;
    conn->edge_to_tree = NULL;
    conn->edge_to_edge = NULL;
  }
  conn->ett_offset = P4EST_ALLOC (p4est_topidx_t, num_edges + 1);
  conn->ett_offset[num_edges] = num_ett;
#endif

  conn->num_corners = num_corners;
  if (num_corners > 0) {
    conn->tree_to_corner =
      P4EST_ALLOC (p4est_topidx_t, P4EST_CHILDREN * num_trees);
    conn->corner_to_tree = P4EST_ALLOC (p4est_topidx_t, num_ctt);
    conn->corner_to_corner = P4EST_ALLOC (int8_t, num_ctt);
  }
  else {
    conn->tree_to_corner = NULL;
    conn->corner_to_tree = NULL;
    conn->corner_to_corner = NULL;
  }
  conn->ctt_offset = P4EST_ALLOC (p4est_topidx_t, num_corners + 1);
  conn->ctt_offset[num_corners] = num_ctt;

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_bcast (p4est_connectivity_t * conn_in, int root,
                          sc_MPI_Comm mpicomm)
{
  int                 mpirank, mpiret;
  p4est_connectivity_t *conn;
  struct
  {
    p4est_topidx_t      num_vertices, num_trees, num_corners, num_ctt;
    size_t              tree_attr_bytes;
#ifdef P4_TO_P8
    p4est_topidx_t      num_edges, num_ett;
#endif
  }
  conn_dimensions;

  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  /* fill dims_buffer on root process */
  if (mpirank == root) {
    P4EST_ASSERT (conn_in != NULL);
    memset (&conn_dimensions, -1, sizeof (conn_dimensions));
    conn = conn_in;
    conn_dimensions.num_corners = conn->num_corners;
    conn_dimensions.num_trees = conn->num_trees;
    conn_dimensions.num_vertices = conn->num_vertices;
    conn_dimensions.tree_attr_bytes = conn->tree_attr_bytes;
    conn_dimensions.num_ctt = conn->ctt_offset[conn->num_corners];
#ifdef P4_TO_P8
    conn_dimensions.num_edges = conn->num_edges;
    conn_dimensions.num_ett = conn->ett_offset[conn->num_edges];
#endif
  }
  else {
    P4EST_ASSERT (conn_in == NULL);
    conn = NULL;                /* suppress 'maybe used ininitialized' warning */
  }
  /* broadcast the dimensions to all processes */
  mpiret = sc_MPI_Bcast (&conn_dimensions, sizeof (conn_dimensions),
                         sc_MPI_BYTE, root, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* allocate memory for new connectivity */
  if (mpirank != root) {
    P4EST_ASSERT (conn == NULL);
    conn = p4est_connectivity_new (conn_dimensions.num_vertices,
                                   conn_dimensions.num_trees,
#ifdef P4_TO_P8
                                   conn_dimensions.num_edges,
                                   conn_dimensions.num_ett,
#endif
                                   conn_dimensions.num_corners,
                                   conn_dimensions.num_ctt);
    p4est_connectivity_set_attr (conn, conn_dimensions.tree_attr_bytes);
  }

  /* Broadcast the arrays if not NULL.  If a pointer is NULL on one process
   * then it is NULL on every process, therefore the if-constructions work */
  if (conn->num_vertices > 0) {
    P4EST_ASSERT (conn->vertices != NULL);
    P4EST_ASSERT (conn->tree_to_vertex != NULL);
    mpiret = sc_MPI_Bcast (conn->vertices, 3 * conn_dimensions.num_vertices,
                           sc_MPI_DOUBLE, root, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (conn->tree_to_vertex,
                           P4EST_CHILDREN * conn_dimensions.num_trees,
                           P4EST_MPI_TOPIDX, root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  mpiret =
    sc_MPI_Bcast (conn->tree_to_tree, P4EST_FACES * conn_dimensions.num_trees,
                  P4EST_MPI_TOPIDX, root, mpicomm);
  SC_CHECK_MPI (mpiret);
  mpiret =
    sc_MPI_Bcast (conn->tree_to_face, P4EST_FACES * conn_dimensions.num_trees,
                  sc_MPI_BYTE, root, mpicomm);
  SC_CHECK_MPI (mpiret);

  if (conn->num_corners > 0) {
    P4EST_ASSERT (conn->tree_to_corner != NULL);
    P4EST_ASSERT (conn->corner_to_tree != NULL);
    P4EST_ASSERT (conn->corner_to_corner != NULL);
    mpiret = sc_MPI_Bcast (conn->tree_to_corner,
                           P4EST_CHILDREN * conn_dimensions.num_trees,
                           P4EST_MPI_TOPIDX, root, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (conn->corner_to_tree, conn_dimensions.num_ctt,
                           P4EST_MPI_TOPIDX, root, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (conn->corner_to_corner, conn_dimensions.num_ctt,
                           sc_MPI_BYTE, root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  mpiret = sc_MPI_Bcast (conn->ctt_offset, conn_dimensions.num_corners,
                         P4EST_MPI_TOPIDX, root, mpicomm);
  P4EST_ASSERT (conn->ctt_offset[conn->num_corners] ==
                conn_dimensions.num_ctt);
  SC_CHECK_MPI (mpiret);
#ifdef P4_TO_P8
  if (conn->num_edges > 0) {
    P4EST_ASSERT (conn->tree_to_edge != NULL);
    P4EST_ASSERT (conn->edge_to_tree != NULL);
    P4EST_ASSERT (conn->edge_to_edge != NULL);
    mpiret = sc_MPI_Bcast (conn->tree_to_edge,
                           P8EST_EDGES * conn_dimensions.num_trees,
                           P4EST_MPI_TOPIDX, root, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (conn->edge_to_tree, conn_dimensions.num_ett,
                           P4EST_MPI_TOPIDX, root, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (conn->edge_to_edge, conn_dimensions.num_ett,
                           sc_MPI_BYTE, root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Bcast (conn->ett_offset, conn_dimensions.num_edges,
                         P4EST_MPI_TOPIDX, root, mpicomm);
  P4EST_ASSERT (conn->ett_offset[conn->num_edges] == conn_dimensions.num_ett);
  SC_CHECK_MPI (mpiret);
#endif

  if (conn->tree_attr_bytes != 0) {
    mpiret = sc_MPI_Bcast (conn->tree_to_attr,
                           conn->tree_attr_bytes * conn->num_trees,
                           sc_MPI_BYTE, root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  return conn;
}

void
p4est_connectivity_destroy (p4est_connectivity_t * conn)
{
  P4EST_FREE (conn->vertices);
  P4EST_FREE (conn->tree_to_vertex);

  P4EST_FREE (conn->tree_to_tree);
  P4EST_FREE (conn->tree_to_face);

#ifdef P4_TO_P8
  P4EST_FREE (conn->tree_to_edge);
  P4EST_FREE (conn->ett_offset);
  P4EST_FREE (conn->edge_to_tree);
  P4EST_FREE (conn->edge_to_edge);
#endif

  P4EST_FREE (conn->tree_to_corner);
  P4EST_FREE (conn->ctt_offset);
  P4EST_FREE (conn->corner_to_tree);
  P4EST_FREE (conn->corner_to_corner);

  p4est_connectivity_set_attr (conn, 0);

  P4EST_FREE (conn);
}

void
p4est_connectivity_set_attr (p4est_connectivity_t * conn,
                             size_t bytes_per_tree)
{
  if (bytes_per_tree > 0) {
    P4EST_ASSERT (conn->tree_to_attr == NULL);
    conn->tree_to_attr = P4EST_ALLOC (char, bytes_per_tree * conn->num_trees);
  }
  else {
    P4EST_FREE (conn->tree_to_attr);
    conn->tree_to_attr = NULL;
  }
  conn->tree_attr_bytes = bytes_per_tree;
}

int
p4est_connectivity_is_valid (p4est_connectivity_t * conn)
{
  int                 nvert;
  int                 face, rface, nface, orientation;
  int                 errcode, errcount;
#ifdef P4_TO_P8
  int                 edge, nedge;
  int                 flip, nflip, nflip1, nflip2;
  p4est_topidx_t      aedge, edge_begin, edge_end;
  p4est_topidx_t      nett;
#endif
  int                 corner, ncorner;
  int                 good, cfound;
  p4est_topidx_t      vertex, tree, ntree;
  p4est_topidx_t      acorner, corner_begin, corner_end;
  p4est_topidx_t      nctt;
  const p4est_topidx_t num_vertices = conn->num_vertices;
  const p4est_topidx_t num_trees = conn->num_trees;
  const p4est_topidx_t *ttv = conn->tree_to_vertex;
  const p4est_topidx_t *ttt = conn->tree_to_tree;
  const int8_t       *ttf = conn->tree_to_face;
#ifdef P4_TO_P8
  const p4est_topidx_t num_edges = conn->num_edges;
  const p4est_topidx_t *tte = conn->tree_to_edge;
  const p4est_topidx_t *eoff = conn->ett_offset;
  const p4est_topidx_t *ett = conn->edge_to_tree;
  const int8_t       *ete = conn->edge_to_edge;
  const p4est_topidx_t num_ett = eoff[num_edges];
  p8est_edge_info_t   ei;
  sc_array_t         *eta = &ei.edge_transforms;
#endif
  const p4est_topidx_t num_corners = conn->num_corners;
  const p4est_topidx_t *ttc = conn->tree_to_corner;
  const p4est_topidx_t *coff = conn->ctt_offset;
  const p4est_topidx_t *ctt = conn->corner_to_tree;
  const int8_t       *ctc = conn->corner_to_corner;
  const p4est_topidx_t num_ctt = coff[num_corners];
  p4est_corner_info_t ci;
  sc_array_t         *cta = &ci.corner_transforms;

  good = 0;
#ifdef P4_TO_P8
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  if (num_vertices == 0 && (conn->vertices != NULL || ttv != NULL)) {
    P4EST_NOTICE ("Zero vertices still with arrays\n");
    goto failure;
  }
  if (num_vertices > 0 && (conn->vertices == NULL || ttv == NULL)) {
    P4EST_NOTICE ("Nonzero vertices missing arrays\n");
    goto failure;
  }

#ifdef P4_TO_P8
  for (nett = 0; nett < num_ett; ++nett) {
    if (ett[nett] < 0 || ett[nett] >= num_trees) {
      P4EST_NOTICEF ("Edge to tree %lld out of range\n", (long long) nett);
      goto failure;
    }
    if (ete[nett] < 0 || ete[nett] >= 24) {
      P4EST_NOTICEF ("Edge to edge %lld out of range\n", (long long) nett);
      goto failure;
    }
  }
#endif

  for (nctt = 0; nctt < num_ctt; ++nctt) {
    if (ctt[nctt] < 0 || ctt[nctt] >= num_trees) {
      P4EST_NOTICEF ("Corner to tree %lld out of range\n", (long long) nctt);
      goto failure;
    }
    if (ctc[nctt] < 0 || ctc[nctt] >= P4EST_CHILDREN) {
      P4EST_NOTICEF ("Corner to corner %lld out of range\n",
                     (long long) nctt);
      goto failure;
    }
  }

  if (num_vertices > 0) {
    for (tree = 0; tree < num_trees; ++tree) {
      for (nvert = 0; nvert < P4EST_CHILDREN; ++nvert) {
        vertex = ttv[tree * P4EST_CHILDREN + nvert];
        if (vertex < 0 || vertex >= num_vertices) {
          P4EST_NOTICEF ("Tree to vertex out of range %lld %d",
                         (long long) tree, nvert);
          goto failure;
        }
      }
    }
  }

  if ((conn->tree_to_attr != NULL) != (conn->tree_attr_bytes > 0)) {
    P4EST_NOTICEF ("Tree attribute properties inconsistent %lld",
                   (long long) conn->tree_attr_bytes);
    goto failure;
  }

  for (tree = 0; tree < num_trees; ++tree) {
    for (face = 0; face < P4EST_FACES; ++face) {
      ntree = ttt[tree * P4EST_FACES + face];
      if (ntree < 0 || ntree >= num_trees) {
        P4EST_NOTICEF ("Tree to tree out of range %lld %d\n",
                       (long long) tree, face);
        goto failure;
      }
      rface = (int) ttf[tree * P4EST_FACES + face];
      if (rface < 0 || rface >= P4EST_FACES * P4EST_HALF) {
        P4EST_NOTICEF ("Tree to face out of range %lld %d\n",
                       (long long) tree, face);
        goto failure;
      }
      nface = rface % P4EST_FACES;      /* clamp to a real face index */
      orientation = rface / P4EST_FACES;        /* 0..P4EST_HALF-1 */
      if (ntree == tree) {
        /* no neighbor across this face or self-periodic */
        if (nface == face && orientation != 0) {
          P4EST_NOTICEF ("Face invalid in %lld %d\n", (long long) tree, face);
          goto failure;
        }
      }
      if (ntree != tree || nface != face) {
        /* check reciprocity */
        if (ttt[ntree * P4EST_FACES + nface] != tree) {
          P4EST_NOTICEF ("Tree to tree reciprocity in %lld %d\n",
                         (long long) tree, face);
          goto failure;
        }
        if ((int) ttf[ntree * P4EST_FACES + nface] !=
            face + P4EST_FACES * orientation) {
          P4EST_NOTICEF ("Tree to face reciprocity in %lld %d\n",
                         (long long) tree, face);
          goto failure;
        }
      }
    }

#ifdef P4_TO_P8
    for (aedge = 0; aedge < num_edges; ++aedge) {
      if (eoff[aedge + 1] < eoff[aedge]) {
        P4EST_NOTICEF ("Edge offset backwards %lld\n", (long long) aedge);
        goto failure;
      }
    }
    if (num_edges > 0) {
      for (edge = 0; edge < P8EST_EDGES; ++edge) {
        p8est_find_edge_transform (conn, tree, edge, &ei);
        aedge = tte[tree * P8EST_EDGES + edge];
        if (aedge < -1 || aedge >= num_edges) {
          P4EST_NOTICEF ("Tree to edge out of range %lld %d\n",
                         (long long) tree, edge);
          goto failure;
        }
        if (aedge == -1) {
          continue;
        }
        errcode = errcount = 0;
        flip = nflip1 = nflip2 = -1;
        edge_begin = eoff[aedge];
        edge_end = eoff[aedge + 1];
        if (edge_begin < 0 || edge_begin >= num_ett ||
            edge_end < 0 || edge_end > num_ett) {
          P4EST_NOTICEF ("Invalid edge range %lld %d\n",
                         (long long) tree, edge);
          goto failure;
        }
        for (nett = edge_begin; nett < edge_end; ++nett) {
          ntree = ett[nett];
          nedge = (int) ete[nett] % P8EST_EDGES;
          if (tte[ntree * P8EST_EDGES + nedge] != aedge) {
            P4EST_NOTICEF ("Edge to edge reciprocity in %lld %d %lld\n",
                           (long long) tree, edge, (long long) nett);
            goto failure;
          }
          nflip = (int) ete[nett] / P8EST_EDGES;
          if (ntree == tree && nedge == edge) {
            if (flip != -1 && nflip == flip) {
              errcode = 1;
              break;
            }
            flip = nflip;
            continue;
          }
          ++errcount;
        }
        if (errcode > 0) {
          P4EST_NOTICEF ("Shared edge %lld %d %lld inconsistency %d\n",
                         (long long) tree, edge, (long long) nett, errcode);
          goto failure;
        }
        if (flip == -1 ||
            !(nflip1 == -1 || nflip1 == flip) ||
            !(nflip2 == -1 || nflip2 == flip)) {
          P4EST_NOTICEF ("Shared edge %lld %d inconsistent flip\n",
                         (long long) tree, edge);
          goto failure;
        }
      }
    }
#endif

    for (acorner = 0; acorner < num_corners; ++acorner) {
      if (coff[acorner + 1] < coff[acorner]) {
        P4EST_NOTICEF ("Corner offset backwards %lld\n", (long long) acorner);
        goto failure;
      }
    }
    if (num_corners > 0) {
      for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
        p4est_find_corner_transform (conn, tree, corner, &ci);
        acorner = ttc[tree * P4EST_CHILDREN + corner];
        if (acorner < -1 || acorner >= num_corners) {
          P4EST_NOTICEF ("Tree to corner out of range %lld %d\n",
                         (long long) tree, corner);
          goto failure;
        }
        if (acorner == -1) {
          continue;
        }
        errcode = errcount = 0;
        cfound = 0;
        corner_begin = coff[acorner];
        corner_end = coff[acorner + 1];
        if (corner_begin < 0 || corner_begin >= num_ctt ||
            corner_end < 0 || corner_end > num_ctt) {
          P4EST_NOTICEF ("Invalid corner range %lld %d\n",
                         (long long) tree, corner);
          goto failure;
        }
        for (nctt = corner_begin; nctt < corner_end; ++nctt) {
          ntree = ctt[nctt];
          ncorner = (int) ctc[nctt];
          if (ttc[ntree * P4EST_CHILDREN + ncorner] != acorner) {
            P4EST_NOTICEF ("Corner to corner reciprocity in %lld %d %lld\n",
                           (long long) tree, corner, (long long) nctt);
            goto failure;
          }
          if (ntree == tree && ncorner == corner) {
            if (cfound) {
              errcode = 1;
              break;
            }
            cfound = 1;
            continue;
          }
          ++errcount;
        }
        if (errcode > 0) {
          P4EST_NOTICEF ("Shared corner %lld %d %lld inconsistency %d\n",
                         (long long) tree, corner, (long long) nctt, errcode);
          goto failure;
        }
        if (!cfound) {
          P4EST_NOTICEF ("Shared corner %lld %d inconsistent count B\n",
                         (long long) tree, corner);
          goto failure;
        }
      }
    }
  }
  good = 1;

failure:
#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);

  return good;
}

int
p4est_connectivity_is_equal (p4est_connectivity_t * conn1,
                             p4est_connectivity_t * conn2)
{
  size_t              topsize, int8size;
  size_t              tcount;
  p4est_topidx_t      num_vertices;
#ifdef P4_TO_P8
  p4est_topidx_t      num_edges, num_ett;
#endif
  p4est_topidx_t      num_corners, num_ctt;

  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);

  if (conn1->num_vertices != conn2->num_vertices ||
      conn1->num_trees != conn2->num_trees ||
#ifdef P4_TO_P8
      conn1->num_edges != conn2->num_edges ||
#endif
      conn1->num_corners != conn2->num_corners) {
    return 0;
  }

  num_vertices = conn1->num_vertices;
  if (num_vertices > 0) {
    P4EST_ASSERT (conn1->vertices != NULL && conn2->vertices != NULL);
    if (memcmp (conn1->vertices, conn2->vertices,
                sizeof (double) * 3 * num_vertices)) {
      return 0;
    }

    P4EST_ASSERT (conn1->tree_to_vertex != NULL &&
                  conn2->tree_to_vertex != NULL);
    tcount = (size_t) (P4EST_CHILDREN * conn1->num_trees);
    if (memcmp (conn1->tree_to_vertex, conn2->tree_to_vertex,
                tcount * topsize)) {
      return 0;
    }
  }

#ifdef P4_TO_P8
  tcount = (size_t) (P8EST_EDGES * conn1->num_trees);
  if (conn1->num_edges > 0 &&
      memcmp (conn1->tree_to_edge, conn2->tree_to_edge, tcount * topsize)) {
    return 0;
  }
#endif

  tcount = (size_t) (P4EST_CHILDREN * conn1->num_trees);
  if (conn1->num_corners > 0 &&
      memcmp (conn1->tree_to_corner, conn2->tree_to_corner,
              tcount * topsize)) {
    return 0;
  }

  tcount = (size_t) (P4EST_FACES * conn1->num_trees);
  if (memcmp (conn1->tree_to_tree, conn2->tree_to_tree, tcount * topsize) ||
      memcmp (conn1->tree_to_face, conn2->tree_to_face, tcount * int8size)) {
    return 0;
  }

  if ((conn1->tree_to_attr == NULL) != (conn2->tree_to_attr == NULL) ||
      conn1->tree_attr_bytes != conn2->tree_attr_bytes) {
    return 0;
  }
  tcount = (size_t) conn1->num_trees;
  if (conn1->tree_to_attr != NULL &&
      memcmp (conn1->tree_to_attr, conn2->tree_to_attr,
              tcount * conn1->tree_attr_bytes)) {
    return 0;
  }

#ifdef P4_TO_P8
  num_edges = conn1->num_edges;
  num_ett = conn1->ett_offset[num_edges];
  /* when there are no edges, the latter two ranges are zero */
  if (memcmp (conn1->ett_offset, conn2->ett_offset,
              topsize * (num_edges + 1)) ||
      memcmp (conn1->edge_to_tree, conn2->edge_to_tree,
              topsize * num_ett) ||
      memcmp (conn1->edge_to_edge, conn2->edge_to_edge, int8size * num_ett)) {
    return 0;
  }
#endif

  num_corners = conn1->num_corners;
  num_ctt = conn1->ctt_offset[num_corners];
  /* when there are no corners, the latter two ranges are zero */
  if (memcmp (conn1->ctt_offset, conn2->ctt_offset,
              topsize * (num_corners + 1)) ||
      memcmp (conn1->corner_to_tree, conn2->corner_to_tree,
              topsize * num_ctt) ||
      memcmp (conn1->corner_to_corner, conn2->corner_to_corner,
              int8size * num_ctt)) {
    return 0;
  }

  return 1;
}

int
p4est_connectivity_sink (p4est_connectivity_t * conn, sc_io_sink_t * sink)
{
  int                 retval;
  int                 has_tree_attr;
  char                magic8[8 + 1];
  char                pkgversion24[24 + 1];
  size_t              tree_attr_bytes;
  size_t              u64z, topsize, int8size;
  size_t              tcount;
  uint64_t            array10[10];
  p4est_topidx_t      num_vertices, num_trees;
  p4est_topidx_t      num_edges, num_ett, num_corners, num_ctt;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  retval = 0;
  num_vertices = conn->num_vertices;
  num_trees = conn->num_trees;
#ifdef P4_TO_P8
  num_edges = conn->num_edges;
  num_ett = conn->ett_offset[num_edges];
#else
  num_edges = num_ett = 0;
#endif
  num_corners = conn->num_corners;
  num_ctt = conn->ctt_offset[num_corners];
  has_tree_attr = ((tree_attr_bytes = conn->tree_attr_bytes) > 0);

  strncpy (magic8, P4EST_STRING, 8);
  magic8[8] = '\0';
  retval = retval || sc_io_sink_write (sink, magic8, 8);

  strncpy (pkgversion24, P4EST_PACKAGE_VERSION, 24);
  pkgversion24[24] = '\0';
  retval = retval || sc_io_sink_write (sink, pkgversion24, 24);

  u64z = sizeof (uint64_t);
  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);
  array10[0] = P4EST_ONDISK_FORMAT;
  array10[1] = (uint64_t) topsize;
  array10[2] = (uint64_t) num_vertices;
  array10[3] = (uint64_t) num_trees;
  array10[4] = (uint64_t) num_edges;
  array10[5] = (uint64_t) num_ett;
  array10[6] = (uint64_t) num_corners;
  array10[7] = (uint64_t) num_ctt;
  array10[8] = (uint64_t) conn->tree_attr_bytes;
  array10[9] = (uint64_t) 0;
  retval = retval || sc_io_sink_write (sink, array10, 10 * u64z);

  if (num_vertices > 0) {
    tcount = (size_t) (3 * num_vertices);
    retval = retval ||
      sc_io_sink_write (sink, conn->vertices, tcount * sizeof (double));
  }

#ifdef P4_TO_P8
  if (num_edges > 0) {
    tcount = (size_t) (P8EST_EDGES * num_trees);
    retval = retval ||
      sc_io_sink_write (sink, conn->tree_to_edge, tcount * topsize);
  }
#endif

  tcount = (size_t) (P4EST_CHILDREN * num_trees);
  if (num_vertices > 0) {
    retval = retval ||
      sc_io_sink_write (sink, conn->tree_to_vertex, tcount * topsize);
  }
  if (num_corners > 0) {
    retval = retval ||
      sc_io_sink_write (sink, conn->tree_to_corner, tcount * topsize);
  }

  tcount = (size_t) (P4EST_FACES * num_trees);
  retval = retval ||
    sc_io_sink_write (sink, conn->tree_to_tree, tcount * topsize) ||
    sc_io_sink_write (sink, conn->tree_to_face, tcount * int8size);

  if (has_tree_attr) {
    tcount = (size_t) num_trees;
    retval = retval ||
      sc_io_sink_write (sink, conn->tree_to_attr, tcount * tree_attr_bytes);
  }

#ifdef P4_TO_P8
  retval = retval || sc_io_sink_write (sink, conn->ett_offset,
                                       topsize * (num_edges + 1));
  if (num_edges > 0) {
    retval = retval ||
      sc_io_sink_write (sink, conn->edge_to_tree, topsize * num_ett) ||
      sc_io_sink_write (sink, conn->edge_to_edge, int8size * num_ett);
  }
#endif

  retval = retval || sc_io_sink_write (sink, conn->ctt_offset,
                                       topsize * (num_corners + 1));
  if (num_corners > 0) {
    retval = retval ||
      sc_io_sink_write (sink, conn->corner_to_tree, topsize * num_ctt) ||
      sc_io_sink_write (sink, conn->corner_to_corner, int8size * num_ctt);
  }

  return retval;
}

sc_array_t         *
p4est_connectivity_deflate (p4est_connectivity_t * conn,
                            p4est_connectivity_encode_t code)
{
  int                 retval;
  sc_array_t         *buffer;
  sc_io_sink_t       *sink;

  buffer = sc_array_new (sizeof (char));

  /* This sink writes to a memory buffer so no file errors are caught. */
  sink = sc_io_sink_new (SC_IO_TYPE_BUFFER, SC_IO_MODE_WRITE,
                         SC_IO_ENCODE_NONE, buffer);
  SC_CHECK_ABORT (sink != NULL, "sink open from buffer");

  retval = p4est_connectivity_sink (conn, sink);
  SC_CHECK_ABORT (retval == 0, "sink connectivity");

  retval = sc_io_sink_destroy (sink);
  SC_CHECK_ABORT (retval == 0, "destroy sink");

  return buffer;
}

int
p4est_connectivity_save (const char *filename, p4est_connectivity_t * conn)
{
  int                 retval;
  sc_io_sink_t       *sink;

  sink = sc_io_sink_new (SC_IO_TYPE_FILENAME, SC_IO_MODE_WRITE,
                         SC_IO_ENCODE_NONE, filename);
  if (sink == NULL) {
    return -1;
  }

  /* Close file even on earlier write error */
  retval = p4est_connectivity_sink (conn, sink);
  retval = sc_io_sink_destroy (sink) || retval;

  return retval;
}

p4est_connectivity_t *
p4est_connectivity_source (sc_io_source_t * source)
{
  int                 retval;
  int                 has_tree_attr;
  char                magic8[9];
  char                pkgversion24[25];
  size_t              tree_attr_bytes;
  size_t              u64z, topsize, int8size;
  size_t              tcount;
  uint64_t            array10[10];
  p4est_topidx_t      num_vertices, num_trees;
  p4est_topidx_t      num_edges, num_ett, num_corners, num_ctt;
  p4est_connectivity_t *conn = NULL;

  retval = sc_io_source_read (source, magic8, 8, NULL);
  magic8[8] = '\0';
  if (retval || strncmp (magic8, P4EST_STRING, 8)) {
    /* "invalid magic" */
    return NULL;
  }
  retval = sc_io_source_read (source, pkgversion24, 24, NULL);
  pkgversion24[24] = '\0';
  if (retval) {
    /* "read package version" */
    return NULL;
  }

  u64z = sizeof (uint64_t);
  topsize = sizeof (p4est_topidx_t);
  int8size = sizeof (int8_t);
  retval = sc_io_source_read (source, array10, 10 * u64z, NULL);
  if (retval) {
    /*"read header" */
    return NULL;
  }
  if (array10[0] != P4EST_ONDISK_FORMAT) {
    /* "on-disk format mismatch" */
    return NULL;
  }
  if (array10[1] != (uint64_t) topsize) {
    /* "p4est_topidx_t size mismatch" */
    return NULL;
  }
  num_vertices = (p4est_topidx_t) array10[2];
  num_trees = (p4est_topidx_t) array10[3];
  num_edges = (p4est_topidx_t) array10[4];
  num_ett = (p4est_topidx_t) array10[5];
  num_corners = (p4est_topidx_t) array10[6];
  num_ctt = (p4est_topidx_t) array10[7];
  has_tree_attr = ((tree_attr_bytes = (size_t) array10[8]) > 0);
  if (num_vertices < 0) {
    /* "negative num_vertices" */
    return NULL;
  }
  if (num_trees < 0) {
    /* "negative num_trees" */
    return NULL;
  }
#ifdef P4_TO_P8
  if (num_edges < 0) {
    /* "negative num_edges" */
    return NULL;
  }
  if (num_ett < 0) {
    /* "negative num_ett" */
    return NULL;
  }
#else
  if (num_edges != 0) {
    /* "num_edges must be zero in 2D" */
    return NULL;
  }
  if (num_ett != 0) {
    /* "num_ett must be zero in 2D" */
    return NULL;
  }
#endif
  if (num_corners < 0) {
    /* "negative num_corners" */
    return NULL;
  }
  if (num_ctt < 0) {
    /* "negative num_ctt" */
    return NULL;
  }

  conn = p4est_connectivity_new (num_vertices, num_trees,
#ifdef P4_TO_P8
                                 num_edges, num_ett,
#endif
                                 num_corners, num_ctt);
  p4est_connectivity_set_attr (conn, tree_attr_bytes);

  if (num_vertices > 0) {
    tcount = (size_t) (3 * num_vertices);
    retval = sc_io_source_read (source, conn->vertices,
                                tcount * sizeof (double), NULL);
    if (retval) {
      /* "read vertices" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }

#ifdef P4_TO_P8
  if (num_edges > 0) {
    tcount = (size_t) (P8EST_EDGES * num_trees);
    retval = sc_io_source_read (source, conn->tree_to_edge, topsize * tcount,
                                NULL);
    if (retval) {
      /* "read tte" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }
#endif

  tcount = (size_t) (P4EST_CHILDREN * num_trees);
  if (num_vertices > 0) {
    retval = sc_io_source_read (source, conn->tree_to_vertex,
                                topsize * tcount, NULL);
    if (retval) {
      /* "read ttv" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }
  if (num_corners > 0) {
    retval = sc_io_source_read (source, conn->tree_to_corner,
                                topsize * tcount, NULL);
    if (retval) {
      /* "read ttc" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }
  tcount = (size_t) (P4EST_FACES * num_trees);
  retval = sc_io_source_read (source, conn->tree_to_tree, topsize * tcount,
                              NULL);
  if (retval) {
    /* "read ttt" */
    p4est_connectivity_destroy (conn);
    return NULL;
  }
  retval = sc_io_source_read (source, conn->tree_to_face, int8size * tcount,
                              NULL);
  if (retval) {
    /* "read ttf" */
    p4est_connectivity_destroy (conn);
    return NULL;
  }
  if (has_tree_attr) {
    tcount = (size_t) num_trees;
    retval = sc_io_source_read (source, conn->tree_to_attr,
                                tcount * tree_attr_bytes, NULL);
    if (retval) {
      /* "write tree_to_attr" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }

#ifdef P4_TO_P8
  retval = sc_io_source_read (source, conn->ett_offset,
                              topsize * (num_edges + 1), NULL);
  if (retval || num_ett != conn->ett_offset[num_edges]) {
    /* "read ett_offset" */
    p4est_connectivity_destroy (conn);
    return NULL;
  }
  if (num_edges > 0) {
    retval = sc_io_source_read (source, conn->edge_to_tree, topsize * num_ett,
                                NULL);
    if (retval) {
      /* "read ett" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
    retval = sc_io_source_read (source, conn->edge_to_edge,
                                int8size * num_ett, NULL);
    if (retval) {
      /* "read ete" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }
#endif

  retval = sc_io_source_read (source, conn->ctt_offset,
                              topsize * (num_corners + 1), NULL);
  if (retval || num_ctt != conn->ctt_offset[num_corners]) {
    /* "read ctt_offset" */
    p4est_connectivity_destroy (conn);
    return NULL;
  }
  if (num_corners > 0) {
    retval = sc_io_source_read (source, conn->corner_to_tree,
                                topsize * num_ctt, NULL);
    if (retval) {
      /* "read ctt" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
    retval = sc_io_source_read (source, conn->corner_to_corner,
                                int8size * num_ctt, NULL);
    if (retval) {
      /* "read ctc" */
      p4est_connectivity_destroy (conn);
      return NULL;
    }
  }

  if (!p4est_connectivity_is_valid (conn)) {
    /* "invalid connectivity" */
    p4est_connectivity_destroy (conn);
    return NULL;
  }

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_inflate (sc_array_t * buffer)
{
  int                 retval;
  p4est_connectivity_t *conn;
  sc_io_source_t     *source;

  /* This source reads from a memory buffer so no file errors are caught. */
  source = sc_io_source_new (SC_IO_TYPE_BUFFER, SC_IO_ENCODE_NONE, buffer);
  SC_CHECK_ABORT (source != NULL, "source open from buffer");

  conn = p4est_connectivity_source (source);

  retval = sc_io_source_destroy (source);
  SC_CHECK_ABORT (retval == 0, "destroy source");

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_load (const char *filename, size_t *bytes)
{
  int                 retval;
  size_t              bytes_in;
  sc_io_source_t     *source;
  p4est_connectivity_t *conn;

  source = sc_io_source_new (SC_IO_TYPE_FILENAME,
                             SC_IO_ENCODE_NONE, filename);
  if (source == NULL) {
    return NULL;
  }

  /* Get byte length and close file even on earlier read error */
  conn = p4est_connectivity_source (source);
  retval = sc_io_source_complete (source, &bytes_in, NULL) || conn == NULL;
  retval = sc_io_source_destroy (source) || retval;
  if (retval) {
    if (conn != NULL) {
      p4est_connectivity_destroy (conn);
    }
    return NULL;
  }

  if (bytes != NULL) {
    *bytes = bytes_in;
  }
  return conn;
}

#ifndef P4_TO_P8

p4est_connectivity_t *
p4est_connectivity_new_unitsquare (void)
{
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[4 * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };
  const p4est_topidx_t tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[1 * 4] = {
    0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 4] = {
    0, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

#endif /* !P4_TO_P8 */

p4est_connectivity_t *
p4est_connectivity_new_periodic (void)
{
  const p4est_topidx_t num_vertices = P4EST_CHILDREN;
  const p4est_topidx_t num_trees = 1;
#ifdef P4_TO_P8
  const p4est_topidx_t num_edges = 3;
#endif /* P4_TO_P8 */
  const p4est_topidx_t num_corners = 1;
  const double        vertices[P4EST_CHILDREN * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
#ifdef P4_TO_P8
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
#endif /* P4_TO_P8 */
  };
  const p4est_topidx_t tree_to_vertex[1 * P4EST_CHILDREN] = {
    0, 1, 2, 3,
#ifdef P4_TO_P8
    4, 5, 6, 7,
#endif /* P4_TO_P8 */
  };
  const p4est_topidx_t tree_to_tree[1 * P4EST_FACES] = {
    0, 0, 0, 0,
#ifdef P4_TO_P8
    0, 0,
#endif /* P4_TO_P8 */
  };
  const int8_t        tree_to_face[1 * P4EST_FACES] = {
    1, 0, 3, 2,
#ifdef P4_TO_P8
    5, 4,
#endif /* P4_TO_P8 */
  };
#ifdef P4_TO_P8
  const p4est_topidx_t tree_to_edge[1 * P8EST_EDGES] = {
    0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
  };
  const p4est_topidx_t ett_offset[3 + 1] = {
    0, 4, 8, 12,
  };
  const p4est_topidx_t edge_to_tree[P8EST_EDGES] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int8_t        edge_to_edge[P8EST_EDGES] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
  };
#endif /* P4_TO_P8 */
  const p4est_topidx_t tree_to_corner[1 * P4EST_CHILDREN] = {
    0, 0, 0, 0,
#ifdef P4_TO_P8
    0, 0, 0, 0,
#endif /* P4_TO_P8 */
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, P4EST_CHILDREN,
  };
  const p4est_topidx_t corner_to_tree[P4EST_CHILDREN] = {
    0, 0, 0, 0,
#ifdef P4_TO_P8
    0, 0, 0, 0,
#endif /* P4_TO_P8 */
  };
  const int8_t        corner_to_corner[P4EST_CHILDREN] = {
    0, 1, 2, 3,
#ifdef P4_TO_P8
    4, 5, 6, 7,
#endif /* P4_TO_P8 */
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
#ifdef P4_TO_P8
                                      num_edges,
#endif /* P4_TO_P8 */
                                      num_corners, vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
#ifdef P4_TO_P8
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
#endif /* P4_TO_P8 */
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

#ifndef P4_TO_P8

p4est_connectivity_t *
p4est_connectivity_new_rotwrap (void)
{
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_corners = 1;
  const double        vertices[4 * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };
  const p4est_topidx_t tree_to_vertex[1 * 4] = {
    0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[1 * 4] = {
    0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 4] = {
    1, 0, 7, 6,
  };
  const p4est_topidx_t tree_to_corner[1 * 4] = {
    0, 0, 0, 0,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 4,
  };
  const p4est_topidx_t corner_to_tree[4] = {
    0, 0, 0, 0,
  };
  const int8_t        corner_to_corner[4] = {
    0, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_circle (void)
{
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[12 * 3] = {
    /* inner hexagon */
    0.0, 1.0, 0.0,
    0.866025404, 0.5, 0.0,
    0.866025404, -0.5, 0.0,
    0, -1.0, 0.0,
    -0.866025404, -0.5, 0.0,
    -0.866025404, 0.5, 0.0,
    /* outer hexagon */
    0.0, 2.0, 0.0,
    1.73205081, 1.0, 0.0,
    1.73205081, -1.0, 0.0,
    0, -2.0, 0.0,
    -1.73205081, -1.0, 0.0,
    -1.73205081, 1.0, 0.0,
  };
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
    7, 6, 1, 0,
    11, 5, 6, 0,
    5, 11, 4, 10,
    9, 3, 10, 4,
    2, 3, 8, 9,
    8, 7, 2, 1,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
    5, 1, 0, 0,
    1, 1, 2, 0,
    2, 2, 1, 3,
    3, 3, 4, 2,
    5, 3, 4, 4,
    4, 0, 5, 5,
  };
  const int8_t        tree_to_face[6 * 4] = {
    1, 3, 2, 3,
    0, 1, 6, 1,
    0, 1, 6, 7,
    0, 1, 5, 7,
    4, 6, 2, 3,
    4, 0, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_drop (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 10;
  const p4est_topidx_t num_trees = 5;
  const p4est_topidx_t num_ctt = 1;
  const double        vertices[10 * 3] = {
    0, 0, 0,
    1, 0, 0,
    3, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    1, 2, 0,
    2, 2, 0,
    0, 3, 0,
    3, 3, 0,
  };
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 1, 3, 4,
    1, 2, 4, 5,
    5, 2, 7, 9,
    6, 7, 8, 9,
    3, 4, 8, 6,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    0, 1, 0, 4,
    0, 2, 1, 1,
    2, 2, 1, 3,
    4, 2, 3, 3,
    4, 4, 0, 3,
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 0, 2, 2,
    1, 2, 2, 3,
    0, 1, 1, 1,
    3, 3, 2, 3,
    0, 1, 3, 0,
  };

  const p4est_topidx_t tree_to_corner[5 * 4] = {
    -1, -1, -1,  0,
    -1, -1,  0, -1,
    -1, -1, -1, -1,
    -1, -1, -1, -1,
    -1,  0, -1, -1,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 3
  };
  const p4est_topidx_t corner_to_tree[3] = {
    0, 1, 4,
  };
  const int8_t        corner_to_corner[3] = {
    3, 2, 1,
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_ctt,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_corner (void)
{
  const p4est_topidx_t num_vertices = 7;
  const p4est_topidx_t num_trees = 3;
  const p4est_topidx_t num_corners = 1;
  const double        vertices[7 * 3] = {
    -1, -1, 0,
    0, -1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1,
    0, 1, 1,
    -1, 0, 0,
  };
  const p4est_topidx_t tree_to_vertex[3 * 4] = {
    0, 1, 2, 3, 0, 2, 6, 5, 2, 3, 5, 4,
  };
  const p4est_topidx_t tree_to_tree[3 * 4] = {
    1, 0, 0, 2, 1, 2, 0, 1, 1, 2, 0, 2,
  };
  const int8_t        tree_to_face[3 * 4] = {
    2, 1, 2, 2, 0, 0, 0, 3, 1, 1, 3, 3,
  };
  const p4est_topidx_t tree_to_corner[3 * 4] = {
    -1, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, -1,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 3,
  };
  const p4est_topidx_t corner_to_tree[3] = {
    0, 1, 2,
  };
  const int8_t        corner_to_corner[3] = {
    2, 1, 0,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_pillow (void)
{
  const p4est_topidx_t num_vertices = 4;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[4 * 3] = {
    0, 0, 0,
    0, 1, 0,
    1, 0, 0,
    1, 1, 0,
  };
  const p4est_topidx_t tree_to_vertex[2 * 4] = {
    0, 1, 2, 3, 0, 1, 2, 3,
  };
  const p4est_topidx_t tree_to_tree[2 * 4] = {
    1, 1, 1, 1, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[2 * 4] = {
    0, 1, 2, 3, 0, 1, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_moebius (void)
{
  const p4est_topidx_t num_vertices = 10;
  const p4est_topidx_t num_trees = 5;
  const p4est_topidx_t num_ctt = 0;
  const double        halfsqrt3 = .5 * sqrt (3.);
  const double        vertices[10 * 3] = {
    0, 0, 0,
    0, 1, 0,
    1, 0, 0,
    1, 1, 0,
    1.5, 0, halfsqrt3,
    1.5, 1, halfsqrt3,
    .5, .5, 1.5,
    .5, .5, 2,
    -.5, 0, halfsqrt3,
    -.5, 1, halfsqrt3,
  };
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    0, 2, 1, 3,
    3, 5, 2, 4,
    4, 6, 5, 7,
    6, 7, 9, 8,
    9, 8, 1, 0,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    4, 1, 0, 0,
    0, 2, 1, 1,
    1, 3, 2, 2,
    3, 3, 2, 4,
    4, 4, 3, 0,
  };
  const int8_t        tree_to_face[5 * 4] = {
    7, 4, 2, 3,
    5, 4, 2, 3,
    5, 2, 2, 3,
    0, 1, 1, 2,
    0, 1, 3, 4,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_star (void)
{
  const p4est_topidx_t num_vertices = 13;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_corners = 1;
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
    0, 1, 3, 2,
    0, 3, 5, 4,
    5, 6, 0, 7,
    8, 7, 9, 0,
    9, 0, 10, 11,
    12, 1, 11, 0,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
    1, 0, 5, 0,
    2, 1, 0, 1,
    1, 2, 2, 3,
    3, 2, 3, 4,
    4, 5, 3, 4,
    5, 0, 5, 4,
  };
  const int8_t        tree_to_face[6 * 4] = {
    2, 1, 5, 3,
    4, 1, 0, 3,
    4, 1, 2, 5,
    0, 7, 2, 2,
    0, 7, 3, 3,
    0, 6, 2, 5,
  };
  const p4est_topidx_t tree_to_corner[6 * 4] = {
    0, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1,
    -1, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1, 0,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 6,
  };
  const p4est_topidx_t corner_to_tree[6] = {
    0, 1, 2, 3, 4, 5,
  };
  const int8_t        corner_to_corner[6] = {
    0, 0, 2, 3, 1, 3,
  };
  const double        pi = 4.0 * atan (1.0);
  const double        r1 = 1.;
  const double        r2 = 1.5;
  double              vertices[13 * 3];
  int                 i;

  vertices[0 * 3 + 0] = 0;
  vertices[0 * 3 + 1] = 0;
  vertices[0 * 3 + 2] = 0;
  for (i = 0; i < 6; ++i) {
    vertices[(2 * i + 1) * 3 + 0] = r1 * cos (i * pi / 3);
    vertices[(2 * i + 1) * 3 + 1] = r1 * sin (i * pi / 3);
    vertices[(2 * i + 1) * 3 + 2] = 0;
    vertices[(2 * i + 2) * 3 + 0] = r2 * cos ((i + .5) * pi / 3);
    vertices[(2 * i + 2) * 3 + 1] = r2 * sin ((i + .5) * pi / 3);
    vertices[(2 * i + 2) * 3 + 2] = 0;
  }

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_cubed (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[8 * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
    0, 2, 1, 3,
    2, 6, 3, 7,
    0, 4, 2, 6,
    4, 5, 6, 7,
    0, 1, 4, 5,
    1, 3, 5, 7,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
    4, 1, 2, 5,
    0, 3, 2, 5,
    0, 3, 4, 1,
    2, 5, 4, 1,
    2, 5, 0, 3,
    4, 1, 0, 3,
  };
  const int8_t        tree_to_face[6 * 4] = {
    2, 0, 0, 2,
    1, 3, 3, 1,
    2, 0, 0, 2,
    1, 3, 3, 1,
    2, 0, 0, 2,
    1, 3, 3, 1,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_icosahedron (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 22;
  const p4est_topidx_t num_trees    = 10;
  const p4est_topidx_t num_corners  =  2;
  const double         vertices[22 * 3] = {
    0.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 00 */
    1.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 01 */
    2.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 02 */
    3.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 03 */
    4.0 +   cos(M_PI/3),    sin(M_PI/3),  0.0, /* vertex 04 */
    0.0,  0.0,  0.0,                           /* vertex 05 */
    1.0,  0.0,  0.0,                           /* vertex 06 */
    2.0,  0.0,  0.0,                           /* vertex 07 */
    3.0,  0.0,  0.0,                           /* vertex 08 */
    4.0,  0.0,  0.0,                           /* vertex 09 */
    5.0,  0.0,  0.0,                           /* vertex 10 */
    0.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 11 */
    1.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 12 */
    2.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 13 */
    3.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 14 */
    4.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 15 */
    5.0 +   cos(M_PI/3), -  sin(M_PI/3),  0.0, /* vertex 16 */
    0.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 17 */
    1.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 18 */
    2.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 19 */
    3.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 20 */
    4.0 + 2*cos(M_PI/3), -2*sin(M_PI/3),  0.0, /* vertex 21 */
  };
  const p4est_topidx_t tree_to_vertex[10 * 4] = {
    5,  11,  0,  6, /* tree 0 */
    11, 17,  6, 12, /* tree 1 */
    6,  12,  1,  7, /* tree 2 */
    12, 18,  7, 13, /* tree 3 */
    7,  13,  2,  8, /* tree 4 */
    13, 19,  8, 14, /* tree 5 */
    8,  14,  3,  9, /* tree 6 */
    14, 20,  9, 15, /* tree 7 */
    9,  15,  4, 10, /* tree 8 */
    15, 21, 10, 16, /* tree 9 */
  };
  const p4est_topidx_t tree_to_tree[10 * 4] = {
    8,1,9,2, /* tree 0 */
    0,3,9,2, /* tree 1 */
    0,3,1,4, /* tree 2 */
    2,5,1,4, /* tree 3 */
    2,5,3,6, /* tree 4 */
    4,7,3,6, /* tree 5 */
    4,7,5,8, /* tree 6 */
    6,9,5,8, /* tree 7 */
    6,9,7,0, /* tree 8 */
    8,1,7,0, /* tree 9 */
  };
  const int8_t        tree_to_face[10 * 4] = {
    7,0,3,4, /* tree 0 */
    1,6,5,2, /* tree 1 */
    7,0,3,4, /* tree 2 */
    1,6,5,2, /* tree 3 */
    7,0,3,4, /* tree 4 */
    1,6,5,2, /* tree 5 */
    7,0,3,4, /* tree 6 */
    1,6,5,2, /* tree 7 */
    7,0,3,4, /* tree 8 */
    1,6,5,2, /* tree 9 */
  };
  const p4est_topidx_t tree_to_corner[10 * 4] = {
    -1,  -1,  0,  -1, /* tree 0 */
    -1,   1, -1,  -1, /* tree 1 */
    -1,  -1,  0,  -1, /* tree 2 */
    -1,   1, -1,  -1, /* tree 3 */
    -1,  -1,  0,  -1, /* tree 4 */
    -1,   1, -1,  -1, /* tree 5 */
    -1,  -1,  0,  -1, /* tree 6 */
    -1,   1, -1,  -1, /* tree 7 */
    -1,  -1,  0,  -1, /* tree 8 */
    -1,   1, -1,  -1, /* tree 9 */
  };
  const p4est_topidx_t ctt_offset[2+1] = {
    0,5,10,
  };
  /* for each corner, report the tree numbers it is attached to */
  const p4est_topidx_t corner_to_tree[10] = {
    0,2,4,6,8, /* corner 0 */
    1,3,5,7,9, /* corner 1 */
  };

  /* a given corner belong to multiple trees;
   for each tree, we report the index identifying the vertex location
   in the tree_to_vertex.
  e.g. here :
  - corner 0 is vertex 0
  - corner 1 is vertex 17
  For each entry in corner_to_tree, we report the location of the vertex in
  tree_to_vertex
  */
  const int8_t corner_to_corner[10] = {
    2, 2, 2, 2, 2,/* corner 0 (i.e vertex  0) */
    1, 1, 1, 1, 1,/* corner 1 (i.e vertex 17) */
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_shell2d (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees    = 8;
  const p4est_topidx_t num_ctt      = 0;
  const double         vertices[6 * 3] = {
    -1,  1,  0,
     0,  1,  0,
     1,  1,  0,
    -1,  2,  0,
     0,  2,  0,
     1,  2,  0,
  };
  const p4est_topidx_t tree_to_vertex[8 * 4] = {
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
    0, 1, 3, 4,
    1, 2, 4, 5,
  };
  const p4est_topidx_t tree_to_tree[8 * 4] = {
    7, 1, 0, 0,
    0, 2, 1, 1,
    1, 3, 2, 2,
    2, 4, 3, 3,
    3, 5, 4, 4,
    4, 6, 5, 5,
    5, 7, 6, 6,
    6, 0, 7, 7,
  };
  const int8_t        tree_to_face[8 * 4] = {
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
    1, 0, 2, 3,
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_disk2d (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 6;
  const p4est_topidx_t num_trees    = 4+1;
  const p4est_topidx_t num_corners  = 0;
  const p4est_topidx_t num_ctt      = 0;
  const double         vertices[6 * 3] = {
    -1,  -1,  0,
     1,  -1,  0,
    -1,   1,  0,
     1,   1,  0,
    -1,   2,  0,
     1,   2,  0,
  };
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    2, 3, 4, 5, /* tree 0 */
    2, 3, 4, 5, /* tree 1 */
    2, 3, 4, 5, /* tree 2 */
    2, 3, 4, 5, /* tree 3 */
    0, 1, 2, 3, /* tree 4  - center */
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    3, 1, 4,  0,  /* tree 0 */
    0, 2, 4,  1,  /* tree 1 */
    1, 3, 4,  2,  /* tree 2 */
    2, 0, 4,  3,  /* tree 3 */
    2, 0, 1,  3,  /* tree 4 - center */
  };
  const int8_t        tree_to_face[5 * 4] = {
    1, 0, 5, 3, /* tree 0 */
    1, 0, 6, 3, /* tree 1 */
    1, 0, 0, 3, /* tree 2 */
    1, 0, 3, 3, /* tree 3 */
    2, 6, 6, 2, /* tree 4 - center */
  };

/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_bowtie (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 7;
  const p4est_topidx_t num_trees    = 2;
  const p4est_topidx_t num_corners  = 1;
  const double         vertices[7 * 3] = {
    -0.7071,  0.7071, 0,
     0.7071,  0.7071, 0,
    -1.4142,  0,      0,
     0,       0,      0,
     1.4142,  0,      0,
    -0.7071, -0.7071, 0,
     0.7071, -0.7071, 0,
  };
  const p4est_topidx_t tree_to_vertex[2 * 4] = {
    2, 5, 0, 3, /* tree 0 */
    6, 4, 3, 1, /* tree 1 */
  };
  const p4est_topidx_t tree_to_tree[2 * 4] = {
    0, 0, 0, 0,  /* tree 0 */
    1, 1, 1, 1,  /* tree 1 */
  };
  const int8_t        tree_to_face[5 * 4] = {
    0, 1, 2, 3, /* tree 0 */
    0, 1, 2, 3, /* tree 1 */
  };
  const p4est_topidx_t tree_to_corner[2 * 4] = {
    -1, -1, -1, 0,/* tree 0 */
    -1, -1, 0, -1,/* tree 1 */
  };
  const p4est_topidx_t corner_to_tree[2] = {
    0, 1,
  };
  const p4est_topidx_t     ctt_offset[1 + 1] = {
    0, 2,
  };
  const int8_t corner_to_corner[2] = {
    3, 2,
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p4est_connectivity_new_disk_nonperiodic (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 5;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[8 * 3] = {
    -1, -1, 0,
    1, -1, 0,
    -1, 1, 0,
    1, 1, 0,
    -3, -3, 0,
    3, -3, 0,
    -3, 3, 0,
    3, 3, 0,
  };
  const p4est_topidx_t tree_to_vertex[5 * 4] = {
    4, 5, 0, 1,
    4, 0, 6, 2,
    0, 1, 2, 3,
    1, 5, 3, 7,
    2, 3, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[5 * 4] = {
    1, 3, 0, 2,
    1, 2, 0, 4,
    1, 3, 0, 4,
    2, 3, 0, 4,
    1, 3, 2, 4,
  };
  const int8_t        tree_to_face[5 * 4] = {
    2, 6, 2, 2,
    0, 0, 0, 4,
    1, 0, 3, 2,
    1, 1, 5, 1,
    7, 3, 3, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p4est_connectivity_new_disk (int periodic_a, int periodic_b)
{
  const int8_t        in_ctc[8] = { 0, 0, 1, 1, 2, 2, 3, 3 };
  const p4est_topidx_t in_ctt[8] = { 0, 1, 0, 3, 1, 4, 3, 4 };
  int                 i, j;
  int8_t             *ctc;
  p4est_topidx_t      nc;
  p4est_topidx_t     *ttc, *ctt;
  p4est_connectivity_t *conn = p4est_connectivity_new_disk_nonperiodic ();

  /* non-periodic boundary works as before */
  P4EST_ASSERT (conn->num_corners == 0);
  P4EST_ASSERT (conn->tree_to_corner == NULL);
  P4EST_ASSERT (conn->corner_to_tree == NULL);
  P4EST_ASSERT (conn->corner_to_corner == NULL);
  if (!periodic_a && !periodic_b) {
    return conn;
  }

  /* allocate arrays of proper size */
  P4EST_FREE (conn->ctt_offset);
  ttc = conn->tree_to_corner = P4EST_ALLOC (p4est_topidx_t, 5 * 4);
  ctt = conn->corner_to_tree = P4EST_ALLOC (p4est_topidx_t, 8);
  ctc = conn->corner_to_corner = P4EST_ALLOC (int8_t, 8);
  nc = conn->num_corners = (periodic_a ^ periodic_b) ? 2 : 1;
  conn->ctt_offset = P4EST_ALLOC (p4est_topidx_t, nc + 1);

  /* fill new arrays with proper data */
  conn->ctt_offset[0] = 0;
  if (nc == 1) {
    conn->ctt_offset[1] = 8;
  }
  else {
    P4EST_ASSERT (nc == 2);
    conn->ctt_offset[1] = 4;
    conn->ctt_offset[2] = 8;
  }
  for (i = 0; i < 8; ++i) {
    /* we have either 1 or 2 connecting corners */
    conn->corner_to_corner[0] = i < 4 || nc == 1 ? 0 : 1;
  }
  if (periodic_a) {
    /* tree 1, face 0 meets tree 3, face 1 */
    conn->tree_to_tree[4 * 1 + 0] = 3;
    conn->tree_to_face[4 * 1 + 0] = 1;
    conn->tree_to_tree[4 * 3 + 1] = 1;
    conn->tree_to_face[4 * 3 + 1] = 0;
  }
  if (periodic_b) {
    /* tree 0, face 2 meets tree 4, face 3 */
    conn->tree_to_tree[4 * 0 + 2] = 4;
    conn->tree_to_face[4 * 0 + 2] = 3;
    conn->tree_to_tree[4 * 4 + 3] = 0;
    conn->tree_to_face[4 * 4 + 3] = 2;
  }
  /* assign corner trees */
  memset (ttc, -1, 5 * 4 * sizeof (p4est_topidx_t));
  ttc[4 * 0 + 0] = ttc[4 * 1 + 0] = 0;
  ttc[4 * 0 + 1] = ttc[4 * 3 + 1] = !periodic_a;
  ttc[4 * 1 + 2] = ttc[4 * 4 + 2] = !periodic_b;
  ttc[4 * 3 + 3] = ttc[4 * 4 + 3] = !periodic_a || !periodic_b;
  /* assign corner trees and corners */
  for (i = 0; i < 8; ++i) {
    j = i < 2 || i >= 6 ? i : !periodic_a ? ((i - 2) ^ 2) + 2 : i;
    ctt[i] = in_ctt[j];
    ctc[i] = in_ctc[j];
  }

  /* return modified connectivity */
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  return conn;
}

#endif /* !P4_TO_P8 */

p4est_connectivity_t *
p4est_connectivity_new_twotrees (int l_face, int r_face, int orientation)
{
  int                 i;
  const p4est_topidx_t num_vertices = (P4EST_DIM - 1) * 6;      /* 6 or 12 */
  const p4est_topidx_t num_trees = 2;

  /* no tree connection via edges and corners */
#ifdef P4_TO_P8
  int                 op;
  const p4est_topidx_t num_edges = 0;
  const p4est_topidx_t num_ett = 0;
#endif /* P4_TO_P8 */
  const p4est_topidx_t num_corners = 0;
  const p4est_topidx_t num_ctt = 0;

/* *INDENT-OFF* */
  const double        vertices[(P4EST_DIM - 1) * 6 * 3] = {
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0
#ifdef P4_TO_P8
           ,
    0, 0, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1
#endif /* P4_TO_P8 */
  };

  /* define mapping from tree to vertex for each face */
  const int           leftTree[P4EST_FACES][P4EST_CHILDREN] =
#ifndef P4_TO_P8
    {{ 1, 0, 4, 3 },
     { 0, 1, 3, 4 },
     { 1, 4, 0, 3 },
     { 0, 3, 1, 4}};
#else /* !P4_TO_P8 */
    {{  1,  0,  7,  6,  4,  3, 10,  9 },
     {  0,  1,  3,  4,  6,  7,  9, 10 },
     {  1,  4,  0,  3,  7, 10,  6,  9 },
     {  0,  6,  1,  7,  3,  9,  4, 10 },
     {  1,  7,  4, 10,  0,  6,  3,  9 },
     {  0,  3,  6,  9,  1,  4,  7, 10 }};
#endif /* !P4_TO_P8 */

  const int           rightTree[P4EST_FACES][P4EST_CHILDREN] =
#ifndef P4_TO_P8
    {{ 1, 2, 4, 5 },
     { 2, 1, 5, 4 },
     { 1, 4, 2, 5 },
     { 2, 5, 1, 4 }};
#else /* !P4_TO_P8 */
    {{  1,  2,  4,  5,  7,  8, 10, 11 },
     {  2,  1,  8,  7,  5,  4, 11, 10 },
     {  1,  7,  2,  8,  4, 10,  5, 11 },
     {  2,  5,  1,  4,  8, 11,  7, 10 },
     {  1,  4,  7, 10,  2,  5,  8, 11 },
     {  2,  8,  5, 11,  1,  7,  4, 10 }};
#endif /* !P4_TO_P8 */

  /* define rotations for right tree in order to set the orientation */
  const int flip[(P4EST_DIM - 1) * 6] =
#ifndef P4_TO_P8
    { -1,  4,  5, -1,  1,  2 };
#else /* !P4_TO_P8 */
    { -1, 10, 11, -1,  7,  8, -1,  4,  5, -1,  1,  2 };
  const int rotateClockWise[(P4EST_DIM - 1) * 6] =
    { -1,  7,  8, -1,  1,  2, -1, 10, 11, -1,  4,  5 };
  const int rotateCounterClockWise[(P4EST_DIM - 1) * 6] =
    { -1,  4,  5, -1, 10, 11, -1,  1,  2, -1,  7,  8 };
#endif /* !P4_TO_P8 */
/* *INDENT-ON* */

  /* initialize values in tree_to_vertex */
  p4est_topidx_t      tree_to_vertex[P4EST_CHILDREN * 2] = {
    -1, -1, -1, -1, -1, -1, -1, -1
#ifdef P4_TO_P8
      ,
    -1, -1, -1, -1, -1, -1, -1, -1
#endif /* P4_TO_P8 */
  };

/* *INDENT-OFF* */
  /* create tree_to_tree and tree_to_face */
  p4est_topidx_t tree_to_tree[2 * P4EST_FACES] =
#ifndef P4_TO_P8
    {0, 0, 0, 0,
     1, 1, 1, 1};
#else /* !P4_TO_P8 */
    {0, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1};
#endif /* !P4_TO_P8 */
  int8_t tree_to_face[2 * P4EST_FACES] =
#ifndef P4_TO_P8
    {0, 1, 2, 3,
     0, 1, 2, 3,};
#else /* !P4_TO_P8 */
    {0, 1, 2, 3, 4, 5,
     0, 1, 2, 3, 4, 5};
#endif /* !P4_TO_P8 */
/* *INDENT-ON* */

  P4EST_ASSERT (0 <= l_face && l_face < P4EST_FACES);
  P4EST_ASSERT (0 <= r_face && r_face < P4EST_FACES);
  P4EST_ASSERT (0 <= orientation && orientation < P4EST_HALF);

  /* populate according to specified faces */
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    tree_to_vertex[i] = leftTree[l_face][i];
    tree_to_vertex[P4EST_CHILDREN + i] = rightTree[r_face][i];
  }

  /* rotate trees such that the corners fall to the respective places
     as specified */
#ifndef P4_TO_P8
  if (orientation == 1) {
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      tree_to_vertex[P4EST_CHILDREN + i] =
        flip[tree_to_vertex[P4EST_CHILDREN + i]];
    }
  }
#else /* P4_TO_P8 */
  op = -1;
  if (orientation == 3) {
    op = 2;
  }
  else if (1 <= orientation && orientation <= 2) {
    if (l_face <= r_face) {
      op = p8est_face_permutation_refs[l_face][r_face];
    }
    else {
      op = p8est_face_permutation_refs[r_face][l_face];
    }
  }
  switch (op) {
  case 0:                      /* clockwise rotation */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      tree_to_vertex[P4EST_CHILDREN + i] =
        rotateClockWise[tree_to_vertex[P4EST_CHILDREN + i]];
    }
    break;
  case 1:                      /* counterclockwise rotation */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      tree_to_vertex[P4EST_CHILDREN + i] =
        rotateCounterClockWise[tree_to_vertex[P4EST_CHILDREN + i]];
    }
    break;
  case 2:                      /* flip */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      tree_to_vertex[P4EST_CHILDREN + i] =
        flip[tree_to_vertex[P4EST_CHILDREN + i]];
    }
    break;
  default:
    /* we do nothing */
    break;
  }
#endif /* P4_TO_P8 */

  /* set values where trees are connected */
  tree_to_tree[l_face] = 1;
  tree_to_tree[P4EST_FACES + r_face] = 0;

  tree_to_face[l_face] = (int8_t) (P4EST_FACES * orientation + r_face);
  tree_to_face[P4EST_FACES + r_face] =
    (int8_t) (P4EST_FACES * orientation + l_face);

  /* create connectivity structure */
  return p4est_connectivity_new_copy (num_vertices, num_trees,
#ifdef P4_TO_P8
                                      num_edges,
#endif /* P4_TO_P8 */
                                      num_corners, vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
#ifdef P4_TO_P8
                                      NULL, &num_ett, NULL, NULL,
#endif /* P4_TO_P8 */
                                      NULL, &num_ctt, NULL, NULL);
}

static inline void
brick_linear_to_xyz (p4est_topidx_t ti, const int logx[P4EST_DIM],
                     const int rankx[P4EST_DIM], p4est_topidx_t tx[P4EST_DIM])
{
  int                 i, j, k;
  int                 lastlog = 0;

  for (i = 0; i < P4EST_DIM; i++) {
    tx[i] = 0;
  }

  for (i = 0; i < P4EST_DIM - 1; i++) {
    p4est_topidx_t      tempx[3] = { 0, 0, 0 };
    int                 logi = logx[rankx[i]] - lastlog;
    int                 idx[3] = { -1, -1, -1 };
    int                 c = 0;

    for (k = 0; k < P4EST_DIM - i; k++) {
      int                 d = rankx[i + k];

      idx[d] = 0;
    }
    for (k = 0; k < P4EST_DIM; k++) {
      if (idx[k] == 0) {
        idx[k] = c++;
      }
    }

    for (j = 0; j < logi; j++) {
      int                 base = (P4EST_DIM - i) * j;
      int                 shift = (P4EST_DIM - i - 1) * j;

      for (k = 0; k < P4EST_DIM; k++) {
        int                 id = idx[k];

        if (id >= 0) {
          tempx[k] |= (ti & (1 << (base + id))) >> (shift + id);
        }
      }
    }
    for (k = 0; k < P4EST_DIM; k++) {
      tx[k] += (tempx[k] << lastlog);
    }
    lastlog += logi;
    ti >>= (P4EST_DIM - i) * logi;
  }
  tx[rankx[P4EST_DIM - 1]] += (ti << lastlog);
}

static inline       p4est_topidx_t
brick_xyz_to_linear (const p4est_topidx_t tx[P4EST_DIM],
                     const int logx[P4EST_DIM], const int rankx[P4EST_DIM])
{
  int                 i, j, k;
  int                 lastlog = logx[rankx[P4EST_DIM - 2]];
  p4est_topidx_t      ti = tx[rankx[P4EST_DIM - 1]] >> lastlog;

  for (i = P4EST_DIM - 2; i >= 0; i--) {
    p4est_topidx_t      tempx[3] = { 0, 0, 0 };
    int                 logi =
      (i == 0) ? lastlog : lastlog - logx[rankx[i - 1]];
    int                 idx[3] = { -1, -1, -1 };
    int                 c = 0;

    for (k = 0; k < P4EST_DIM - i; k++) {
      int                 d = rankx[i + k];

      idx[d] = 0;
    }
    for (k = 0; k < P4EST_DIM; k++) {
      if (idx[k] == 0) {
        idx[k] = c++;
      }
    }

    ti <<= (P4EST_DIM - i) * logi;
    lastlog -= logi;
    for (k = 0; k < P4EST_DIM; k++) {
      tempx[k] = tx[k] >> lastlog;
    }
    for (j = 0; j < logi; j++) {
      int                 shift = (P4EST_DIM - i - 1) * j;

      for (k = 0; k < P4EST_DIM; k++) {
        int                 id = idx[k];

        if (id >= 0) {
          ti |= (tempx[k] & (1 << j)) << (shift + id);
        }
      }
    }
  }

  return ti;
}

p4est_connectivity_t *
#ifndef P4_TO_P8
p4est_connectivity_new_brick (int mi, int ni, int periodic_a, int periodic_b)
#else
p8est_connectivity_new_brick (int mi, int ni, int pi, int periodic_a,
                              int periodic_b, int periodic_c)
#endif
{
#ifndef P4_TO_P8
  const p4est_topidx_t m = (p4est_topidx_t) mi;
  const p4est_topidx_t n = (p4est_topidx_t) ni;
  const p4est_topidx_t mc = periodic_a ? m : (m - 1);
  const p4est_topidx_t nc = periodic_b ? n : (n - 1);
  const p4est_topidx_t num_trees = m * n;
  const p4est_topidx_t num_corners = mc * nc;
  const p4est_topidx_t num_ctt = P4EST_CHILDREN * num_corners;
  const p4est_topidx_t num_vertices = (m + 1) * (n + 1);
  const int           periodic[P4EST_DIM] = { periodic_a, periodic_b };
  const p4est_topidx_t max[P4EST_DIM] = { m - 1, n - 1 };
#else
  const p4est_topidx_t m = (p4est_topidx_t) mi;
  const p4est_topidx_t n = (p4est_topidx_t) ni;
  const p4est_topidx_t p = (p4est_topidx_t) pi;
  const p4est_topidx_t mc = periodic_a ? m : (m - 1);
  const p4est_topidx_t nc = periodic_b ? n : (n - 1);
  const p4est_topidx_t pc = periodic_c ? p : (p - 1);
  const p4est_topidx_t num_trees = m * n * p;
  const p4est_topidx_t num_corners = mc * nc * pc;
  const p4est_topidx_t num_ctt = P4EST_CHILDREN * num_corners;
  const p4est_topidx_t num_edges = m * nc * pc + mc * n * pc + mc * nc * p;
  const p4est_topidx_t num_ett = 4 * num_edges;
  const p4est_topidx_t num_vertices = (m + 1) * (n + 1) * (p + 1);
  const int           periodic[P4EST_DIM] = { periodic_a, periodic_b,
    periodic_c
  };
  const p4est_topidx_t max[P4EST_DIM] = { m - 1, n - 1, p - 1 };
#endif
  double             *vertices;
  p4est_topidx_t     *tree_to_vertex;
  p4est_topidx_t     *tree_to_tree;
  int8_t             *tree_to_face;
  p4est_topidx_t     *tree_to_corner;
  p4est_topidx_t     *ctt_offset;
  p4est_topidx_t     *corner_to_tree;
  int8_t             *corner_to_corner;
  p4est_topidx_t      n_iter;
  int                 logx[P4EST_DIM];
  int                 rankx[P4EST_DIM];
  int                 i, j, l;
  p4est_topidx_t      ti, tj, tk;
  p4est_topidx_t      tx, ty;
  p4est_topidx_t      tf[P4EST_FACES], tc[P4EST_CHILDREN];
  p4est_topidx_t      coord[P4EST_DIM], coord2[P4EST_DIM], ttemp;
  p4est_topidx_t     *linear_to_tree;
  p4est_topidx_t     *tree_to_corner2;
  p4est_topidx_t      vcount = 0, vicount = 0;
  int                 c[P4EST_DIM];
  p4est_connectivity_t *conn;
#ifdef P4_TO_P8
  p4est_topidx_t      tl;
  p4est_topidx_t      tz;
  p4est_topidx_t      te[P8EST_EDGES];
  p4est_topidx_t     *tree_to_edge;
  p4est_topidx_t     *ett_offset;
  p4est_topidx_t     *edge_to_tree;
  int8_t             *edge_to_edge;
  p4est_topidx_t     *tree_to_edge2;
  int                 dir1, dir2;
#endif

#ifndef P4_TO_P8
  P4EST_ASSERT (m > 0 && n > 0);
#else
  P4EST_ASSERT (m > 0 && n > 0 && p > 0);
#endif

  conn = p4est_connectivity_new (num_vertices, num_trees,
#ifdef P4_TO_P8
                                 num_edges, num_ett,
#endif
                                 num_corners, num_ctt);

  vertices = conn->vertices;
  tree_to_vertex = conn->tree_to_vertex;
  tree_to_tree = conn->tree_to_tree;
  tree_to_face = conn->tree_to_face;
#ifdef P4_TO_P8
  tree_to_edge = conn->tree_to_edge;
  ett_offset = conn->ett_offset;
  edge_to_tree = conn->edge_to_tree;
  edge_to_edge = conn->edge_to_edge;
#endif
  tree_to_corner = conn->tree_to_corner;
  ctt_offset = conn->ctt_offset;
  corner_to_tree = conn->corner_to_tree;
  corner_to_corner = conn->corner_to_corner;

#ifdef P4_TO_P8
  for (ti = 0; ti < num_edges + 1; ti++) {
    ett_offset[ti] = 4 * ti;
  }
#endif

  for (ti = 0; ti < num_corners + 1; ti++) {
    ctt_offset[ti] = P4EST_CHILDREN * ti;
  }

  for (ti = 0; ti < P4EST_CHILDREN * num_trees; ti++) {
    tree_to_vertex[ti] = -1;
  }

  logx[0] = SC_LOG2_32 (m - 1) + 1;
  logx[1] = SC_LOG2_32 (n - 1) + 1;
  n_iter = (1 << logx[0]) * (1 << logx[1]);
  if (logx[0] <= logx[1]) {
    rankx[0] = 0;
    rankx[1] = 1;
  }
  else {
    rankx[0] = 1;
    rankx[1] = 0;
  }
#ifdef P4_TO_P8
  logx[2] = SC_LOG2_32 (p - 1) + 1;
  n_iter *= (1 << logx[2]);
  if (logx[2] < logx[rankx[0]]) {
    rankx[2] = rankx[1];
    rankx[1] = rankx[0];
    rankx[0] = 2;
  }
  else if (logx[rankx[1]] <= logx[2]) {
    rankx[2] = 2;
  }
  else {
    rankx[2] = rankx[1];
    rankx[1] = 2;
  }
#endif

  linear_to_tree = P4EST_ALLOC (p4est_topidx_t, n_iter);
  tree_to_corner2 = P4EST_ALLOC (p4est_topidx_t, num_trees);
#ifdef P4_TO_P8
  tree_to_edge2 = P4EST_ALLOC (p4est_topidx_t, 3 * num_trees);
#endif

  tj = 0;
  tk = 0;
#ifdef P4_TO_P8
  tl = 0;
#endif
  for (ti = 0; ti < n_iter; ti++) {
    brick_linear_to_xyz (ti, logx, rankx, coord);
    tx = coord[0];
    ty = coord[1];
#ifdef P4_TO_P8
    tz = coord[2];
#endif
    if (tx < m && ty < n &&
#ifdef P4_TO_P8
        tz < p &&
#endif
        1) {
      linear_to_tree[ti] = tj;
      if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b) &&
#ifdef P4_TO_P8
          (tz < p - 1 || periodic_c) &&
#endif
          1) {
        tree_to_corner2[tj] = tk++;
#ifdef P4_TO_P8
        tree_to_edge2[3 * tj] = tl++;
        tree_to_edge2[3 * tj + 1] = tl++;
        tree_to_edge2[3 * tj + 2] = tl++;
#endif
      }
      else {
        tree_to_corner2[tj] = -1;
#ifdef P4_TO_P8
        if ((ty < n - 1 || periodic_b) && (tz < p - 1 || periodic_c)) {
          tree_to_edge2[3 * tj] = tl++;
        }
        else {
          tree_to_edge2[3 * tj] = -1;
        }
        if ((tx < m - 1 || periodic_a) && (tz < p - 1 || periodic_c)) {
          tree_to_edge2[3 * tj + 1] = tl++;
        }
        else {
          tree_to_edge2[3 * tj + 1] = -1;
        }
        if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b)) {
          tree_to_edge2[3 * tj + 2] = tl++;
        }
        else {
          tree_to_edge2[3 * tj + 2] = -1;
        }
#endif
      }
      tj++;
    }
    else {
      linear_to_tree[ti] = -1;
    }
  }
  P4EST_ASSERT (tj == num_trees);
  P4EST_ASSERT (tk == num_corners);
#ifdef P4_TO_P8
  P4EST_ASSERT (tl == num_edges);
#endif

  for (ti = 0; ti < n_iter; ti++) {
    brick_linear_to_xyz (ti, logx, rankx, coord);
    tx = coord[0];
    ty = coord[1];
#ifdef P4_TO_P8
    tz = coord[2];
#endif
    if (tx < m && ty < n &&
#ifdef P4_TO_P8
        tz < p &&
#endif
        1) {
      tj = linear_to_tree[ti];
      P4EST_ASSERT (tj >= 0);
      for (i = 0; i < P4EST_DIM; i++) {
        for (j = 0; j < 2; j++) {
          l = 2 * i + j;
          coord2[0] = ((tx + ((i == 0) ? (2 * j - 1) : 0)) + m) % m;
          coord2[1] = ((ty + ((i == 1) ? (2 * j - 1) : 0)) + n) % n;
#ifdef P4_TO_P8
          coord2[2] = ((tz + ((i == 2) ? (2 * j - 1) : 0)) + p) % p;
#endif
          tf[l] = brick_xyz_to_linear (coord2, logx, rankx);
          P4EST_ASSERT (tf[l] < n_iter);
          tf[l] = linear_to_tree[tf[l]];
          P4EST_ASSERT (tf[l] >= 0);
        }
#ifdef P4_TO_P8
        for (j = 0; j < 4; j++) {
          l = 4 * i + j;
          coord2[0] = ((tx + ((i == 0) ? 0 : (2 * (j & 1) - 1))) + m) % m;
          coord2[1] = ((ty + ((i == 1) ? 0 :
                              (2 * ((i == 0) ? (j & 1) : (j / 2)) - 1))) +
                       n) % n;
          coord2[2] = ((tz + ((i == 2) ? 0 : (2 * (j / 2) - 1))) + p) % p;
          te[l] = brick_xyz_to_linear (coord2, logx, rankx);
          P4EST_ASSERT (te[l] < n_iter);
          te[l] = linear_to_tree[te[l]];
          P4EST_ASSERT (te[l] >= 0);
        }
#endif
      }
      for (i = 0; i < P4EST_CHILDREN; i++) {
        coord2[0] = ((tx + (((i & 1) == 0) ? -1 : 1)) + m) % m;
        coord2[1] = ((ty + ((((i >> 1) & 1) == 0) ? -1 : 1)) + n) % n;
#ifdef P4_TO_P8
        coord2[2] = ((tz + (((i >> 2) == 0) ? -1 : 1)) + p) % p;
#endif
        tc[i] = brick_xyz_to_linear (coord2, logx, rankx);
        P4EST_ASSERT (tc[i] < n_iter);
        tc[i] = linear_to_tree[tc[i]];
        P4EST_ASSERT (tc[i] >= 0);
      }
      for (i = 0; i < P4EST_DIM; i++) {
        for (j = 0; j < 2; j++) {
          l = i * 2 + j;
          if (!periodic[i] &&
              ((coord[i] == 0 && j == 0) || (coord[i] == max[i] && j == 1))) {
            tree_to_tree[tj * P4EST_FACES + l] = tj;
            tree_to_face[tj * P4EST_FACES + l] = (int8_t) l;
          }
          else {
            tree_to_tree[tj * P4EST_FACES + l] = tf[l];
            tree_to_face[tj * P4EST_FACES + l] = (int8_t) (i * 2 + (j ^ 1));
          }
        }
#ifdef P4_TO_P8
        if (tree_to_edge != NULL) {
          /** dir1, dir2 should be in correct z order */
          dir1 = (i == 0) ? 1 : 0;
          dir2 = (i == 2) ? 1 : 2;
          for (j = 0; j < 4; j++) {
            l = i * 4 + j;
            if ((!periodic[dir1] &&
                 ((coord[dir1] == 0 && (j & 1) == 0) ||
                  (coord[dir1] == max[dir1] && (j & 1) == 1))) ||
                (!periodic[dir2] &&
                 ((coord[dir2] == 0 && (j / 2) == 0) ||
                  (coord[dir2] == max[dir2] && (j / 2) == 1)))) {
              tree_to_edge[tj * P8EST_EDGES + l] = -1;
            }
            else {
              switch (j) {
              case 0:
                ttemp = tree_to_edge2[te[l] * 3 + i];
                break;
              case 1:
                ttemp = tree_to_edge2[tf[dir2 * 2] * 3 + i];
                break;
              case 2:
                ttemp = tree_to_edge2[tf[dir1 * 2] * 3 + i];
                break;
              case 3:
                ttemp = tree_to_edge2[tj * 3 + i];
                break;
              default:
                SC_ABORT_NOT_REACHED ();
              }
              P4EST_ASSERT (ttemp >= 0);
              tree_to_edge[tj * P8EST_EDGES + l] = ttemp;
              edge_to_tree[4 * ttemp + (3 - j)] = tj;
              edge_to_edge[4 * ttemp + (3 - j)] = (int8_t) l;
            }
          }
        }
#endif
      }
      for (i = 0; i < P4EST_CHILDREN; i++) {
        if (tree_to_corner != NULL) {
          c[0] = i & 1;
          c[1] = (i >> 1) & 1;
#ifdef P4_TO_P8
          c[2] = i >> 2;
#endif
          if ((!periodic[0] &&
               ((coord[0] == 0 && c[0] == 0) ||
                (coord[0] == max[0] && c[0] == 1))) ||
              (!periodic[1] &&
               ((coord[1] == 0 && c[1] == 0) ||
                (coord[1] == max[1] && c[1] == 1))) ||
#ifdef P4_TO_P8
              (!periodic[2] &&
               ((coord[2] == 0 && c[2] == 0) ||
                (coord[2] == max[2] && c[2] == 1))) ||
#endif
              0) {
            tree_to_corner[tj * P4EST_CHILDREN + i] = -1;
          }
          else {
            switch (i) {
#ifndef P4_TO_P8
            case 0:
              ttemp = tc[0];
              break;
            case 1:
              ttemp = tf[2];
              break;
            case 2:
              ttemp = tf[0];
              break;
            case 3:
              ttemp = tj;
              break;
#else
            case 0:
              ttemp = tc[0];
              break;
            case 1:
              ttemp = te[0];
              break;
            case 2:
              ttemp = te[4];
              break;
            case 3:
              ttemp = tf[4];
              break;
            case 4:
              ttemp = te[8];
              break;
            case 5:
              ttemp = tf[2];
              break;
            case 6:
              ttemp = tf[0];
              break;
            case 7:
              ttemp = tj;
              break;
#endif
            default:
              SC_ABORT_NOT_REACHED ();
            }
            ttemp = tree_to_corner2[ttemp];
            P4EST_ASSERT (ttemp >= 0);
            tree_to_corner[tj * P4EST_CHILDREN + i] = ttemp;
            corner_to_tree[ttemp * P4EST_CHILDREN +
                           (P4EST_CHILDREN - 1 - i)] = tj;
            corner_to_corner[ttemp * P4EST_CHILDREN +
                             (P4EST_CHILDREN - 1 - i)] = (int8_t) i;
          }
        }
#ifdef P4_TO_P8
        if (tz > 0 && (i >> 2) == 0) {
          tree_to_vertex[tj * P4EST_CHILDREN + i] =
            tree_to_vertex[tf[4] * P4EST_CHILDREN + i + 4];
        }
        else
#endif
        if (ty > 0 && ((i >> 1) & 1) == 0) {
          tree_to_vertex[tj * P4EST_CHILDREN + i] =
            tree_to_vertex[tf[2] * P4EST_CHILDREN + i + 2];
        }
        else if (tx > 0 && (i & 1) == 0) {
          tree_to_vertex[tj * P4EST_CHILDREN + i] =
            tree_to_vertex[tf[0] * P4EST_CHILDREN + i + 1];
        }
        else {
          tree_to_vertex[tj * P4EST_CHILDREN + i] = vcount++;
          vertices[vicount++] = (double) (tx + (i & 1));
          vertices[vicount++] = (double) (ty + ((i >> 1) & 1));
#ifndef P4_TO_P8
          vertices[vicount++] = 0.;
#else
          vertices[vicount++] = (double) (tz + (i >> 2));
#endif
        }
      }
    }
  }

  P4EST_ASSERT (vcount == num_vertices);

  P4EST_FREE (linear_to_tree);
  P4EST_FREE (tree_to_corner2);
#ifdef P4_TO_P8
  P4EST_FREE (tree_to_edge2);
#endif

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  return conn;
}

p4est_connectivity_t *
p4est_connectivity_new_byname (const char *name)
{
#ifndef P4_TO_P8
  if (!strcmp (name, "brick23")) {
    return p4est_connectivity_new_brick (2, 3, 0, 0);
  }
  else if (!strcmp (name, "corner")) {
    return p4est_connectivity_new_corner ();
  }
  else if (!strcmp (name, "cubed")) {
    return p4est_connectivity_new_cubed ();
  }
  else if (!strcmp (name, "disk")) {
    return p4est_connectivity_new_disk (0, 0);
  }
  else if (!strcmp (name, "icosahedron")) {
    return p4est_connectivity_new_icosahedron ();
  }
  else if (!strcmp (name, "moebius")) {
    return p4est_connectivity_new_moebius ();
  }
  else if (!strcmp (name, "periodic")) {
    return p4est_connectivity_new_periodic ();
  }
  else if (!strcmp (name, "pillow")) {
    return p4est_connectivity_new_pillow ();
  }
  else if (!strcmp (name, "rotwrap")) {
    return p4est_connectivity_new_rotwrap ();
  }
  else if (!strcmp (name, "star")) {
    return p4est_connectivity_new_star ();
  }
  else if (!strcmp (name, "shell2d")) {
    return p4est_connectivity_new_shell2d ();
  }
  else if (!strcmp (name, "disk2d")) {
    return p4est_connectivity_new_disk2d ();
  }
  else if (!strcmp (name, "bowtie")) {
    return p4est_connectivity_new_bowtie ();
  }
  else if (!strcmp (name, "unit")) {
    return p4est_connectivity_new_unitsquare ();
  }
#else
  if (!strcmp (name, "brick235")) {
    return p8est_connectivity_new_brick (2, 3, 5, 0, 0, 0);
  }
  else if (!strcmp (name, "periodic")) {
    return p8est_connectivity_new_periodic ();
  }
  else if (!strcmp (name, "rotcubes")) {
    return p8est_connectivity_new_rotcubes ();
  }
  else if (!strcmp (name, "rotwrap")) {
    return p8est_connectivity_new_rotwrap ();
  }
  else if (!strcmp (name, "shell")) {
    return p8est_connectivity_new_shell ();
  }
  else if (!strcmp (name, "sphere")) {
    return p8est_connectivity_new_sphere ();
  }
  else if (!strcmp (name, "twocubes")) {
    return p8est_connectivity_new_twocubes ();
  }
  else if (!strcmp (name, "twowrap")) {
    return p8est_connectivity_new_twowrap ();
  }
  else if (!strcmp (name, "unit")) {
    return p8est_connectivity_new_unitcube ();
  }
#endif
  return NULL;
}

typedef struct
{
  p4est_topidx_t      key[P4EST_HALF];
  p4est_topidx_t      trees[2];
  int8_t              faces[2];
}
p4est_conn_face_info_t;

static unsigned
p4est_conn_face_hash (const void *v, const void *u)
{
  const p4est_conn_face_info_t *fi = (p4est_conn_face_info_t *) v;

#ifdef P4_TO_P8
  return p4est_topidx_hash4 (fi->key);
#else
  return p4est_topidx_hash2 (fi->key);
#endif
}

static int
p4est_conn_face_equal (const void *v1, const void *v2, const void *u)
{
  const p4est_conn_face_info_t *fi1 = (p4est_conn_face_info_t *) v1;
  const p4est_conn_face_info_t *fi2 = (p4est_conn_face_info_t *) v2;

  return !memcmp (fi1->key, fi2->key, P4EST_HALF * sizeof (p4est_topidx_t));
}

static void
p4est_conn_face_key (p4est_topidx_t * key, p4est_topidx_t * ttv, int face)
{
  int                 fc;

  P4EST_ASSERT (0 <= face && face < P4EST_FACES);

  for (fc = 0; fc < P4EST_HALF; ++fc) {
    key[fc] = ttv[p4est_face_corners[face][fc]];
  }
  p4est_topidx_bsort (key, P4EST_HALF);
}

#ifdef P4_TO_P8

typedef struct
{
  p4est_topidx_t      key[2];
  sc_array_t          trees;
  sc_array_t          edges;
  p4est_topidx_t      edgeid;
}
p8est_conn_edge_info_t;

static unsigned
p8est_conn_edge_hash (const void *v, const void *u)
{
  const p8est_conn_edge_info_t *ei = (p8est_conn_edge_info_t *) v;

  return p4est_topidx_hash2 (ei->key);
}

static int
p8est_conn_edge_equal (const void *v1, const void *v2, const void *u)
{
  const p8est_conn_edge_info_t *ei1 = (p8est_conn_edge_info_t *) v1;
  const p8est_conn_edge_info_t *ei2 = (p8est_conn_edge_info_t *) v2;

  return !memcmp (ei1->key, ei2->key, 2 * sizeof (p4est_topidx_t));
}

static void
p8est_conn_edge_key (p4est_topidx_t * key, p4est_topidx_t * ttv, int edge)
{
  int                 ec;

  P4EST_ASSERT (0 <= edge && edge < P8EST_EDGES);

  for (ec = 0; ec < 2; ++ec) {
    key[ec] = ttv[p8est_edge_corners[edge][ec]];
  }
  p4est_topidx_bsort (key, 2);
}

#endif /* P4_TO_P8 */

static void
p4est_expand_face_transform_internal (int iface, int target_face,
                                      int orientation, int ftransform[])
{
#ifdef P4_TO_P8
  int                 reverse;
#ifdef P4EST_ENABLE_DEBUG
  int                 i;
  int                *my_axis = &ftransform[0];
  int                *target_axis = &ftransform[3];
#endif
#endif

  P4EST_ASSERT (0 <= iface && iface < P4EST_FACES);
  P4EST_ASSERT (0 <= target_face && target_face < P4EST_FACES);
  P4EST_ASSERT (0 <= orientation && orientation < P4EST_HALF);

#ifdef P4_TO_P8
  /* the code that was here before is now in test/test_face_transform3.c */
  ftransform[0] = iface < 2 ? 1 : 0;
  ftransform[1] = iface < 4 ? 2 : 1;
  ftransform[2] = iface / 2;
  reverse =
    p8est_face_permutation_refs[0][iface] ^
    p8est_face_permutation_refs[0][target_face] ^
    (orientation == 0 || orientation == 3);
  ftransform[3 + reverse] = target_face < 2 ? 1 : 0;
  ftransform[3 + !reverse] = target_face < 4 ? 2 : 1;
  ftransform[5] = target_face / 2;
  reverse = (p8est_face_permutation_refs[iface][target_face] == 1);
  ftransform[6 + reverse] = (orientation & 1);
  ftransform[6 + !reverse] = (orientation >> 1);
  ftransform[8] = 2 * (iface & 1) + (target_face & 1);

#ifdef P4EST_ENABLE_DEBUG
  for (i = 0; i < 3; ++i) {
    P4EST_ASSERT (0 <= my_axis[i] && my_axis[i] < 3);
    P4EST_ASSERT (0 <= target_axis[i] && target_axis[i] < 3);
  }
  P4EST_ASSERT (my_axis[0] != my_axis[1] &&
                my_axis[0] != my_axis[2] && my_axis[1] != my_axis[2]);
  P4EST_ASSERT (target_axis[0] != target_axis[1] &&
                target_axis[0] != target_axis[2] &&
                target_axis[1] != target_axis[2]);
#endif
#else
  ftransform[2] = iface / 2;
  ftransform[1] = 0;
  ftransform[0] = 1 - ftransform[2];
  ftransform[5] = target_face / 2;
  ftransform[4] = 0;
  ftransform[3] = 1 - ftransform[5];
  ftransform[6] = orientation;
  ftransform[7] = 0;
  ftransform[8] = 2 * (iface & 1) + (target_face & 1);
#endif
}

void
p4est_expand_face_transform (int iface, int nface, int ftransform[])
{
  const int           target_face = nface % P4EST_FACES;
  const int           orientation = nface / P4EST_FACES;

  p4est_expand_face_transform_internal (iface, target_face, orientation,
                                        ftransform);
}

p4est_topidx_t
p4est_find_face_transform (p4est_connectivity_t * connectivity,
                           p4est_topidx_t itree, int iface, int ftransform[])
{
  int                 target_code, target_face, orientation;
  p4est_topidx_t      target_tree;

  P4EST_ASSERT (itree >= 0 && itree < connectivity->num_trees);
  P4EST_ASSERT (iface >= 0 && iface < P4EST_FACES);

  target_tree = connectivity->tree_to_tree[P4EST_FACES * itree + iface];
  target_code = (int) connectivity->tree_to_face[P4EST_FACES * itree + iface];
  target_face = target_code % P4EST_FACES;
  orientation = target_code / P4EST_FACES;

  if (target_tree == itree && target_face == iface) {
    P4EST_ASSERT (orientation == 0);
    return -1;
  }

  p4est_expand_face_transform_internal (iface, target_face, orientation,
                                        ftransform);

  return target_tree;
}

static int
p4est_find_corner_transform_internal (p4est_connectivity_t * conn,
                                      p4est_topidx_t itree, int icorner,
                                      p4est_corner_info_t * ci,
                                      p4est_topidx_t * ctt, int8_t * ctc,
                                      p4est_topidx_t corner_trees)
{
  int                 i, j;
  int                 iface, nface;
  int                 orient, fcorner;
  int                 ncorner, ncode;
  int                 fc, nc;
  p4est_topidx_t      ctree, nctree;
#ifdef P4_TO_P8
  int                 iedge, iwhich, k;
  int                 pref, pset;
  size_t              jz;
  p4est_topidx_t      aedge;
  p8est_edge_info_t   ei;
  sc_array_t         *eta;
  p8est_edge_transform_t *et;
#endif
  sc_array_t         *cta = &ci->corner_transforms;
  p4est_corner_transform_t *ct;
  int                 ndistinct = 1;
  sc_array_t          distinct;
  p4est_topidx_t      ntree;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
  P4EST_ASSERT (cta->elem_size == sizeof (p4est_corner_transform_t));

  sc_array_init_size (&distinct, sizeof (p4est_corner_transform_t), 1);
  ct = (p4est_corner_transform_t *) sc_array_index (&distinct, 0);
  ct->ntree = itree;
  ct->ncorner = icorner;

  /* find the face neighbors */
  for (i = 0; i < P4EST_DIM; ++i) {
    iface = p4est_corner_faces[icorner][i];
    ntree = conn->tree_to_tree[P4EST_FACES * itree + iface];
    ncode = (int) conn->tree_to_face[P4EST_FACES * itree + iface];
    if (ntree != itree || ncode != iface) {     /* not domain boundary */
      nface = ncode % P4EST_FACES;
      orient = ncode / P4EST_FACES;
      fcorner = p4est_corner_face_corners[icorner][iface];
      P4EST_ASSERT (fcorner >= 0);
#ifdef P4_TO_P8
      pref = p8est_face_permutation_refs[iface][nface];
      pset = p8est_face_permutation_sets[pref][orient];
      fc = p8est_face_permutations[pset][fcorner];
#else
      fc = fcorner ^ orient;
#endif
      nc = p4est_face_corners[nface][fc];
      for (j = 0; j < ndistinct; j++) {
        ct = (p4est_corner_transform_t *) sc_array_index_int (&distinct, j);
        if (ntree == ct->ntree && nc == (int) ct->ncorner) {
          break;
        }
      }
      if (j == ndistinct) {
        ct = (p4est_corner_transform_t *) sc_array_push (&distinct);
        ct->ntree = ntree;
        ct->ncorner = nc;
        ndistinct++;
      }
    }
  }

#ifdef P4_TO_P8
  /* find the three edge transforms */
  if (conn->num_edges != 0) {
    for (i = 0; i < 3; ++i) {
      iedge = p8est_corner_edges[icorner][i];
      aedge = conn->tree_to_edge[P8EST_EDGES * itree + iedge];
      if (aedge == -1) {
        continue;
      }
      iwhich = (p8est_edge_corners[iedge][1] == icorner);
      P4EST_ASSERT (p8est_edge_corners[iedge][iwhich] == icorner);

      eta = &ei.edge_transforms;
      sc_array_init (eta, sizeof (p8est_edge_transform_t));
      p8est_find_edge_transform (conn, itree, iedge, &ei);
      for (jz = 0; jz < eta->elem_count; ++jz) {
        et = p8est_edge_array_index (eta, jz);
        ntree = et->ntree;
        nc = p8est_edge_corners[et->nedge][et->nflip ^ iwhich];
        for (k = 0; k < ndistinct; k++) {
          ct = (p4est_corner_transform_t *) sc_array_index_int (&distinct, k);
          if (ntree == ct->ntree && nc == (int) ct->ncorner) {
            break;
          }
        }
        if (k == ndistinct) {
          ct = (p4est_corner_transform_t *) sc_array_push (&distinct);
          ct->ntree = ntree;
          ct->ncorner = nc;
          ndistinct++;
        }
      }
      sc_array_reset (eta);
    }
  }
#endif

  /* collect all corners that are not from face or edge neighbors */
  for (ctree = 0; ctree < corner_trees; ++ctree) {
    nctree = ctt[ctree];
    P4EST_ASSERT (0 <= nctree && nctree < conn->num_trees);
    ncorner = (int) ctc[ctree];
    P4EST_ASSERT (ncorner >= 0 && ncorner < P4EST_CHILDREN);

    /* compare against corners found via self, faces [and edges (3D)] */
    for (j = 0; j < ndistinct; j++) {
      ct = (p4est_corner_transform_t *) sc_array_index_int (&distinct, j);
      if (nctree == ct->ntree && ncorner == (int) ct->ncorner) {
        break;
      }
    }
    if (j < ndistinct) {
      continue;
    }

    /* else we have a true all-diagonal corner with ntree */
    ct = (p4est_corner_transform_t *) sc_array_push (cta);
    ct->ntree = nctree;
    ct->ncorner = (int8_t) ncorner;
  }

  sc_array_reset (&distinct);

  return ndistinct;
}

void
p4est_find_corner_transform (p4est_connectivity_t * conn,
                             p4est_topidx_t itree, int icorner,
                             p4est_corner_info_t * ci)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 ignored;
#endif
  p4est_topidx_t      corner_trees, acorner, cttac;
  sc_array_t         *cta = &ci->corner_transforms;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= icorner && icorner < P4EST_CHILDREN);
  P4EST_ASSERT (cta->elem_size == sizeof (p4est_corner_transform_t));

  /* check if this corner exists at all */
  ci->icorner = (int8_t) icorner;
  sc_array_resize (cta, 0);
  if (conn->num_corners == 0) {
    return;
  }
  acorner = conn->tree_to_corner[P4EST_CHILDREN * itree + icorner];
  if (acorner == -1) {
    return;
  }
  P4EST_ASSERT (0 <= acorner && acorner < conn->num_corners);

  /* retrieve connectivity information for this corner */
  cttac = conn->ctt_offset[acorner];
  corner_trees = conn->ctt_offset[acorner + 1] - cttac;
  P4EST_ASSERT (0 <= cttac && 1 <= corner_trees);

  /* loop through all corner neighbors and find corner connections */
#ifdef P4EST_ENABLE_DEBUG
  ignored =
#else
  (void)
#endif
    p4est_find_corner_transform_internal (conn, itree, icorner, ci,
                                          conn->corner_to_tree +
                                          cttac,
                                          conn->corner_to_corner +
                                          cttac, corner_trees);
  P4EST_ASSERT (corner_trees == (p4est_topidx_t) (cta->elem_count + ignored));
}

void
p4est_connectivity_complete (p4est_connectivity_t * conn)
{
  int                 face, corner, r;
  int                 primary, secondary, j;
  size_t              pz;
  p4est_topidx_t     *pt, treeid, nodeid, tt;
  p4est_topidx_t     *ttv, *whichttv[2];
  p4est_conn_face_info_t fikey, *fi;
  sc_hash_array_t    *face_ha;
#ifdef P4_TO_P8
  int                 edge;
  int8_t             *et;
  size_t              ez, egz;
  p4est_topidx_t     *ept, real_edges, enode[2];
  p4est_topidx_t      ett_count, ett_offset, ett_edge;
  p8est_conn_edge_info_t eikey, *ei;
  p8est_edge_info_t   einfo;
  sc_hash_array_t    *edge_ha;
  sc_array_t          edge_array, edge_to_pz;
  sc_array_t         *eta = &einfo.edge_transforms;
#endif
  int8_t             *ct;
  size_t              zcount;
  p4est_topidx_t      real_corners;
  p4est_topidx_t      ctt_count, ctt_offset, ctt_corner;
  p4est_corner_info_t cinfo;
  sc_array_t         *node_trees, *nt;
  sc_array_t         *node_corners, *nc;
  sc_array_t         *cta = &cinfo.corner_transforms;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  /* prepare data structures and remove previous connectivity information */
  face_ha = sc_hash_array_new (sizeof (p4est_conn_face_info_t),
                               p4est_conn_face_hash, p4est_conn_face_equal,
                               NULL);
#ifdef P4_TO_P8
  edge_ha = sc_hash_array_new (sizeof (p8est_conn_edge_info_t),
                               p8est_conn_edge_hash, p8est_conn_edge_equal,
                               NULL);
  P4EST_FREE (conn->tree_to_edge);
  P4EST_FREE (conn->ett_offset);
  P4EST_FREE (conn->edge_to_tree);
  P4EST_FREE (conn->edge_to_edge);
  real_edges = P8EST_EDGES * conn->num_trees;
  conn->tree_to_edge = P4EST_ALLOC (p4est_topidx_t, real_edges);
  memset (conn->tree_to_edge, -1, real_edges * sizeof (p4est_topidx_t));
  real_edges = 0;
  ett_count = 0;
  sc_array_init (&edge_to_pz, sizeof (p4est_topidx_t));
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  P4EST_FREE (conn->tree_to_corner);
  P4EST_FREE (conn->ctt_offset);
  P4EST_FREE (conn->corner_to_tree);
  P4EST_FREE (conn->corner_to_corner);
  real_corners = P4EST_CHILDREN * conn->num_trees;
  conn->tree_to_corner = P4EST_ALLOC (p4est_topidx_t, real_corners);
  memset (conn->tree_to_corner, -1, real_corners * sizeof (p4est_topidx_t));
  real_corners = 0;
  ctt_count = 0;
  node_trees = P4EST_ALLOC (sc_array_t, conn->num_vertices);
  node_corners = P4EST_ALLOC (sc_array_t, conn->num_vertices);
  for (nodeid = 0; nodeid < conn->num_vertices; ++nodeid) {
    sc_array_init (node_trees + nodeid, sizeof (p4est_topidx_t));
    sc_array_init (node_corners + nodeid, sizeof (int8_t));
  }
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  /* hash all faces and edges to identify connections, map corners */
  ttv = conn->tree_to_vertex;
  for (treeid = 0; treeid < conn->num_trees; ++treeid) {
    for (face = 0; face < P4EST_FACES; ++face) {
      p4est_conn_face_key (fikey.key, ttv, face);
      fi = (p4est_conn_face_info_t *)
        sc_hash_array_insert_unique (face_ha, &fikey, &pz);
      if (fi != NULL) {
        /* added fi to hash array as the first of two faces */
        P4EST_ASSERT (sc_array_position (&face_ha->a, fi) == pz);
        memcpy (fi->key, fikey.key, P4EST_HALF * sizeof (p4est_topidx_t));
        fi->trees[0] = treeid;
        fi->faces[0] = (int8_t) face;
        fi->trees[1] = -1;
        fi->faces[1] = -1;
      }
      else {
        /* found existing entry from the first face */
        fi = (p4est_conn_face_info_t *) sc_array_index (&face_ha->a, pz);
        P4EST_ASSERT (p4est_conn_face_equal (fi->key, fikey.key, NULL));
        P4EST_ASSERT (fi->trees[0] >= 0 && fi->faces[0] >= 0);
        P4EST_ASSERT (fi->trees[1] == -1 && fi->faces[1] == -1);
        fi->trees[1] = treeid;
        fi->faces[1] = (int8_t) face;

        /* find primary face and orientation to store it */
        primary = (fi->faces[0] <= fi->faces[1] ? 0 : 1);
        secondary = 1 - primary;
        whichttv[0] = conn->tree_to_vertex + P4EST_CHILDREN * fi->trees[0];
        whichttv[1] = ttv;
        nodeid = whichttv[primary][p4est_face_corners[fi->faces[primary]][0]];
        for (r = 0; r < P4EST_HALF; ++r) {
          corner = p4est_face_corners[fi->faces[secondary]][r];
          if (nodeid == whichttv[secondary][corner]) {
            break;
          }
        }
        P4EST_ASSERT (r < P4EST_HALF);
        for (j = 0; j < 2; ++j) {
          tt = P4EST_FACES * fi->trees[j] + fi->faces[j];
          conn->tree_to_tree[tt] = fi->trees[1 - j];
          conn->tree_to_face[tt] =
            (int8_t) (P4EST_FACES * r + fi->faces[1 - j]);
        }
      }
    }
#ifdef P4_TO_P8
    for (edge = 0; edge < P8EST_EDGES; ++edge) {
      p8est_conn_edge_key (eikey.key, ttv, edge);
      ei = (p8est_conn_edge_info_t *)
        sc_hash_array_insert_unique (edge_ha, &eikey, &pz);
      P4EST_ASSERT (conn->tree_to_edge[P8EST_EDGES * treeid + edge] == -1);
      if (ei != NULL) {
        /* added ei to hash array as the first of an edge group */
        P4EST_ASSERT (sc_array_position (&edge_ha->a, ei) == pz);
        memcpy (ei->key, eikey.key, 2 * sizeof (p4est_topidx_t));
        ei->edgeid = -1;
        sc_array_init (&ei->trees, sizeof (p4est_topidx_t));
        sc_array_init (&ei->edges, sizeof (int8_t));
      }
      else {
        /* this edge has been found before */
        ei = (p8est_conn_edge_info_t *) sc_array_index (&edge_ha->a, pz);
        P4EST_ASSERT (p8est_conn_edge_equal (ei->key, eikey.key, NULL));
        P4EST_ASSERT (ei->trees.elem_count == ei->edges.elem_count);
        if (ei->trees.elem_count == 1) {
          /* store number of this real edge and fill previous tree */
          P4EST_ASSERT (ei->edgeid == -1);
          ei->edgeid = real_edges++;
          ept = (p4est_topidx_t *) sc_array_push (&edge_to_pz);
          *ept = (p4est_topidx_t) pz;
          pt = (p4est_topidx_t *) sc_array_index (&ei->trees, 0);
          P4EST_ASSERT (0 <= *pt && *pt < conn->num_trees);
          et = (int8_t *) sc_array_index (&ei->edges, 0);
          P4EST_ASSERT (0 <= *et && *et < P8EST_EDGES);
          P4EST_ASSERT (conn->tree_to_edge[P8EST_EDGES * *pt + *et] == -1);
          conn->tree_to_edge[P8EST_EDGES * *pt + *et] = ei->edgeid;
          ett_count += 2;
        }
        else {
          /* check that edge is initialized */
          P4EST_ASSERT (ei->trees.elem_count > 1);
          P4EST_ASSERT (0 <= ei->edgeid && ei->edgeid < real_edges);
          ++ett_count;
        }
        conn->tree_to_edge[P8EST_EDGES * treeid + edge] = ei->edgeid;
      }

      /* store edge information */
      pt = (p4est_topidx_t *) sc_array_push (&ei->trees);
      *pt = treeid;
      et = (int8_t *) sc_array_push (&ei->edges);
      *et = (int8_t) edge;
    }
#endif
    for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
      nodeid = ttv[corner];
      P4EST_ASSERT (0 <= nodeid && nodeid < conn->num_vertices);
      nt = node_trees + nodeid;
      nc = node_corners + nodeid;
      zcount = nt->elem_count;
      P4EST_ASSERT (zcount == nc->elem_count);
      if (zcount == 1) {
        ctt_count += 2;
      }
      else if (zcount > 1) {
        ++ctt_count;
      }
      P4EST_ASSERT (conn->tree_to_corner[P4EST_CHILDREN * treeid + corner] ==
                    -1);
      conn->tree_to_corner[P4EST_CHILDREN * treeid + corner] = nodeid;
      pt = (p4est_topidx_t *) sc_array_push (nt);
      *pt = treeid;
      ct = (int8_t *) sc_array_push (nc);
      *ct = (int8_t) corner;
    }
    ttv += P4EST_CHILDREN;
  }

  /* single faces are on the boundary and need not be changed */
  sc_hash_array_destroy (face_ha);

  /* complete edge identification */
#ifdef P4_TO_P8
  P4EST_ASSERT (edge_to_pz.elem_count == (size_t) real_edges);
  conn->num_edges = real_edges;
  conn->ett_offset = P4EST_ALLOC (p4est_topidx_t, conn->num_edges + 1);
  conn->edge_to_tree = P4EST_ALLOC (p4est_topidx_t, ett_count);
  conn->edge_to_edge = P4EST_ALLOC (int8_t, ett_count);
  sc_hash_array_rip (edge_ha, &edge_array);
  ett_edge = 0;
  ett_offset = 0;
  real_edges = 0;

  /* loop through all connected edges */
  for (ez = 0; ez < edge_to_pz.elem_count; ++ez) {
    ept = (p4est_topidx_t *) sc_array_index (&edge_to_pz, ez);
    ei = (p8est_conn_edge_info_t *) sc_array_index (&edge_array, *ept);
    P4EST_ASSERT (ei->trees.elem_count > 1);
    P4EST_ASSERT (ei->trees.elem_count == ei->edges.elem_count);
    P4EST_ASSERT (0 <= ei->edgeid && ei->edgeid < conn->num_edges);
    P4EST_ASSERT ((size_t) ei->edgeid == ez);

    /* set up edge connection information */
    for (egz = 0; egz < ei->trees.elem_count; ++egz) {
      pt = (p4est_topidx_t *) sc_array_index (&ei->trees, egz);
      et = (int8_t *) sc_array_index (&ei->edges, egz);
      P4EST_ASSERT (0 <= *pt && *pt < conn->num_trees);
      P4EST_ASSERT (0 <= *et && *et < P8EST_EDGES);
      P4EST_ASSERT (conn->tree_to_edge[P8EST_EDGES * *pt + *et]
                    == ei->edgeid);
      if (real_edges > 0) {
        conn->tree_to_edge[P8EST_EDGES * *pt + *et] -= real_edges;
      }
      for (j = 0; j < 2; ++j) {
        enode[j] = conn->tree_to_vertex[P4EST_CHILDREN * *pt
                                        + p8est_edge_corners[*et][j]];
      }
      P4EST_ASSERT (enode[0] != enode[1]);
      conn->edge_to_tree[ett_offset + egz] = *pt;
      conn->edge_to_edge[ett_offset + egz] =
        *et + (enode[0] < enode[1] ? 0 : P8EST_EDGES);
    }
    ei->edgeid -= real_edges;

    /* determine if this edge is redundant */
    for (egz = 0; egz < ei->trees.elem_count; ++egz) {
      pt = (p4est_topidx_t *) sc_array_index (&ei->trees, egz);
      et = (int8_t *) sc_array_index (&ei->edges, egz);
      einfo.iedge = -1;         /* unused */
      p8est_find_edge_transform_internal (conn, *pt, *et, &einfo,
                                          conn->edge_to_tree +
                                          ett_offset,
                                          conn->edge_to_edge +
                                          ett_offset, ei->trees.elem_count);
      if (eta->elem_count != 0) {
        /* edge is non-redundant */
        break;
      }
    }

    if (eta->elem_count == 0) {
      /* erase all references to this redundant edge */
      for (egz = 0; egz < ei->trees.elem_count; ++egz) {
        pt = (p4est_topidx_t *) sc_array_index (&ei->trees, egz);
        et = (int8_t *) sc_array_index (&ei->edges, egz);
        P4EST_ASSERT (conn->tree_to_edge[P8EST_EDGES * *pt + *et]
                      == ei->edgeid);
        conn->tree_to_edge[P8EST_EDGES * *pt + *et] = -1;
      }
      ei->edgeid = -1;
      ++real_edges;
    }
    else {
      /* accept edge as non-redundant */
      sc_array_reset (eta);
      conn->ett_offset[ett_edge++] = ett_offset;
      ett_offset += (p4est_topidx_t) ei->trees.elem_count;
    }
  }
  sc_array_reset (&edge_to_pz);
  P4EST_ASSERT (ett_edge == conn->num_edges - real_edges);
  P4EST_ASSERT (ett_offset <= ett_count);
  conn->ett_offset[ett_edge] = ett_offset;
  if (real_edges > 0) {
    conn->num_edges -= real_edges;
    conn->ett_offset = P4EST_REALLOC (conn->ett_offset, p4est_topidx_t,
                                      conn->num_edges + 1);
    conn->edge_to_tree =
      P4EST_REALLOC (conn->edge_to_tree, p4est_topidx_t, ett_offset);
    conn->edge_to_edge =
      P4EST_REALLOC (conn->edge_to_edge, int8_t, ett_offset);
  }
  /* clean up storage for all edges */
  for (ez = 0; ez < edge_array.elem_count; ++ez) {
    ei = (p8est_conn_edge_info_t *) sc_array_index (&edge_array, ez);
    P4EST_ASSERT (ei->trees.elem_count >= 1);
    P4EST_ASSERT (ei->trees.elem_count == ei->edges.elem_count);
    P4EST_ASSERT (-1 <= ei->edgeid && ei->edgeid < conn->num_edges);
#ifdef P4EST_ENABLE_DEBUG
    if (ei->trees.elem_count == 1) {
      /* isolated edge does not count */
      P4EST_ASSERT (ei->edgeid == -1);
      pt = (p4est_topidx_t *) sc_array_index (&ei->trees, 0);
      et = (int8_t *) sc_array_index (&ei->edges, 0);
      P4EST_ASSERT (0 <= *pt && *pt < conn->num_trees);
      P4EST_ASSERT (0 <= *et && *et < P8EST_EDGES);
      P4EST_ASSERT (conn->tree_to_edge[P8EST_EDGES * *pt + *et] == -1);
    }
#endif
    sc_array_reset (&ei->trees);
    sc_array_reset (&ei->edges);
  }
  sc_array_reset (&edge_array);
#endif /* P4_TO_P8 */

  /* complete corner identification */
  P4EST_ASSERT (real_corners == 0);
  conn->num_corners = conn->num_vertices;
  conn->ctt_offset = P4EST_ALLOC (p4est_topidx_t, conn->num_corners + 1);
  conn->corner_to_tree = P4EST_ALLOC (p4est_topidx_t, ctt_count);
  conn->corner_to_corner = P4EST_ALLOC (int8_t, ctt_count);
  ctt_corner = 0;
  ctt_offset = 0;
  for (nodeid = 0; nodeid < conn->num_vertices; ++nodeid) {
    nt = node_trees + nodeid;
    nc = node_corners + nodeid;
    zcount = nt->elem_count;
    P4EST_ASSERT (zcount == nc->elem_count);
    if (zcount <= 1) {
      /* isolated corner does not count */
      if (zcount == 1) {
        pt = (p4est_topidx_t *) sc_array_index (nt, 0);
        ct = (int8_t *) sc_array_index (nc, 0);
        P4EST_ASSERT (0 <= *pt && *pt < conn->num_trees);
        P4EST_ASSERT (0 <= *ct && *ct < P4EST_CHILDREN);
        P4EST_ASSERT (conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct] ==
                      nodeid);
        conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct] = -1;
      }
      ++real_corners;
    }
    else {
      /* set up corner connection information */
      for (pz = 0; pz < zcount; ++pz) {
        pt = (p4est_topidx_t *) sc_array_index (nt, pz);
        ct = (int8_t *) sc_array_index (nc, pz);
        P4EST_ASSERT (0 <= *pt && *pt < conn->num_trees);
        P4EST_ASSERT (0 <= *ct && *ct < P4EST_CHILDREN);
        P4EST_ASSERT (conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct]
                      == nodeid);
        if (real_corners > 0) {
          conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct] -= real_corners;
        }
        conn->corner_to_tree[ctt_offset + pz] = *pt;
        conn->corner_to_corner[ctt_offset + pz] = *ct;
      }

      /* determine if this corner is redundant */
      for (pz = 0; pz < zcount; ++pz) {
        pt = (p4est_topidx_t *) sc_array_index (nt, pz);
        ct = (int8_t *) sc_array_index (nc, pz);
        cinfo.icorner = -1;     /* unused */
        (void)
          p4est_find_corner_transform_internal (conn, *pt, *ct, &cinfo,
                                                conn->corner_to_tree +
                                                ctt_offset,
                                                conn->corner_to_corner +
                                                ctt_offset, zcount);
        if (cta->elem_count != 0) {
          /* corner is non-redundant */
          break;
        }
      }

      if (cta->elem_count == 0) {
        /* erase all references to this redundant corner */
        for (pz = 0; pz < zcount; ++pz) {
          pt = (p4est_topidx_t *) sc_array_index (nt, pz);
          ct = (int8_t *) sc_array_index (nc, pz);
          P4EST_ASSERT (conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct]
                        == nodeid - real_corners);
          conn->tree_to_corner[P4EST_CHILDREN * *pt + *ct] = -1;
        }
        ++real_corners;
      }
      else {
        /* accept corner as non-redundant */
        sc_array_reset (cta);
        conn->ctt_offset[ctt_corner++] = ctt_offset;
        ctt_offset += (p4est_topidx_t) zcount;
      }
    }
  }
  P4EST_ASSERT (ctt_corner == conn->num_corners - real_corners);
  P4EST_ASSERT (ctt_offset <= ctt_count);
  conn->ctt_offset[ctt_corner] = ctt_offset;
  if (real_corners > 0) {
    conn->num_corners -= real_corners;
    conn->ctt_offset = P4EST_REALLOC (conn->ctt_offset, p4est_topidx_t,
                                      conn->num_corners + 1);
    conn->corner_to_tree =
      P4EST_REALLOC (conn->corner_to_tree, p4est_topidx_t, ctt_offset);
    conn->corner_to_corner =
      P4EST_REALLOC (conn->corner_to_corner, int8_t, ctt_offset);
  }
  /* clean up storage for all corners */
  for (nodeid = 0; nodeid < conn->num_vertices; ++nodeid) {
    sc_array_reset (node_trees + nodeid);
    sc_array_reset (node_corners + nodeid);
  }
  P4EST_FREE (node_trees);
  P4EST_FREE (node_corners);

  /* and be done */
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
}

void
p4est_connectivity_reduce (p4est_connectivity_t * conn)
{
  conn->num_corners = 0;
  conn->ctt_offset[conn->num_corners] = 0;
  P4EST_FREE (conn->tree_to_corner);
  P4EST_FREE (conn->corner_to_tree);
  P4EST_FREE (conn->corner_to_corner);
  conn->tree_to_corner = NULL;
  conn->corner_to_tree = NULL;
  conn->corner_to_corner = NULL;
#ifdef P4_TO_P8
  conn->num_edges = 0;
  conn->ett_offset[conn->num_edges] = 0;
  P4EST_FREE (conn->tree_to_edge);
  P4EST_FREE (conn->edge_to_tree);
  P4EST_FREE (conn->edge_to_edge);
  conn->tree_to_edge = NULL;
  conn->edge_to_tree = NULL;
  conn->edge_to_edge = NULL;
#endif
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
}

void
p4est_connectivity_permute (p4est_connectivity_t * conn, sc_array_t * inperm,
                            int is_current_to_new)
{
  sc_array_t         *permarray;
  size_t             *perm;
  p4est_topidx_t      ti, ntrees = conn->num_trees;
  p4est_topidx_t      tj, count;
  sc_array_t          array_view;
  int                 j;

  /* we want the permutation to be the current to new map, not
   * the new to current map */
  if (is_current_to_new) {
    permarray = inperm;
    perm = (size_t *) permarray->array;
  }
  else {
    permarray = sc_array_new_size (sizeof (size_t), (size_t) ntrees);
    perm = (size_t *) permarray->array;
    for (ti = 0; ti < ntrees; ti++) {
      size_t              mapti = *((size_t *) sc_array_index (inperm, ti));
      P4EST_ASSERT (mapti < (size_t) ntrees);
      perm[mapti] = (size_t) ti;
    }
  }

  /* first we change the entries in the various tables */

  /* tree_to_tree */
  for (ti = 0; ti < ntrees; ti++) {
    for (j = 0; j < P4EST_FACES; j++) {
      tj = conn->tree_to_tree[P4EST_FACES * ti + j];
      conn->tree_to_tree[P4EST_FACES * ti + j] = (p4est_topidx_t) perm[tj];
    }
  }

#ifdef P4_TO_P8
  /* edge_to_tree */
  if (conn->edge_to_tree != NULL) {
    count = conn->ett_offset[conn->num_edges];
    for (ti = 0; ti < count; ti++) {
      tj = conn->edge_to_tree[ti];
      conn->edge_to_tree[ti] = (p4est_topidx_t) perm[tj];
    }
  }
#endif

  /* corner_to_tree */
  if (conn->corner_to_tree != NULL) {
    count = conn->ctt_offset[conn->num_corners];
    for (ti = 0; ti < count; ti++) {
      tj = conn->corner_to_tree[ti];
      conn->corner_to_tree[ti] = (p4est_topidx_t) perm[tj];
    }
  }

  /* now we reorder the various tables via in-place permutation */

  /* tree_to_vertex */
  sc_array_init_data (&array_view, conn->tree_to_vertex,
                      P4EST_CHILDREN * sizeof (p4est_topidx_t), ntrees);
  sc_array_permute (&array_view, permarray, 1);

  /* tree_to_tree */
  sc_array_init_data (&array_view, conn->tree_to_tree,
                      P4EST_FACES * sizeof (p4est_topidx_t), ntrees);
  sc_array_permute (&array_view, permarray, 1);

  /* tree_to_face */
  sc_array_init_data (&array_view, conn->tree_to_face,
                      P4EST_FACES * sizeof (int8_t), ntrees);
  sc_array_permute (&array_view, permarray, 1);

#ifdef P4_TO_P8
  /* tree_to_edge */
  if (conn->tree_to_edge != NULL) {
    sc_array_init_data (&array_view, conn->tree_to_edge,
                        P8EST_EDGES * sizeof (p4est_topidx_t), ntrees);
    sc_array_permute (&array_view, permarray, 1);
  }
#endif

  /* tree_to_corner */
  if (conn->tree_to_corner != NULL) {
    sc_array_init_data (&array_view, conn->tree_to_corner,
                        P4EST_CHILDREN * sizeof (p4est_topidx_t), ntrees);
    sc_array_permute (&array_view, permarray, 1);
  }

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  if (!is_current_to_new) {
    sc_array_destroy (permarray);
  }
}

#ifdef P4EST_WITH_METIS

static int
reorder_comp (const void *a, const void *b)
{
  const int          *A = (const int *) a;
  const int          *B = (const int *) b;

  if (A[0] < B[0]) {
    return -1;
  }
  else if (B[0] < A[0]) {
    return 1;
  }
  else {
    return (A[1] - B[1]);
  }
}

void
p4est_connectivity_reorder (sc_MPI_Comm comm, int k,
                            p4est_connectivity_t * conn,
                            p4est_connect_type_t ctype)
{
  sc_array_t         *newid = sc_array_new (sizeof (size_t));
  p4est_connectivity_reorder_newid (comm, k, conn, ctype, newid);
  sc_array_destroy (newid);
}

sc_array_t         *
p4est_connectivity_reorder_newid (sc_MPI_Comm comm, int k,
                                  p4est_connectivity_t * conn,
                                  p4est_connect_type_t ctype,
                                  sc_array_t * newid)
{
  const int           n = (int) conn->num_trees;
  int                 metis_n;
  int                *xadj;
  int                *adjncy;
  int                *part;
  int                 totaldeg;
  int                 degree;
  int                 i, j, l;
  int                 rank;
  p4est_corner_info_t ci;
  sc_array_t         *cta = &ci.corner_transforms;
  p4est_corner_transform_t *ct;
#ifdef P4_TO_P8
  p8est_edge_info_t   ei;
  sc_array_t         *eta = &ei.edge_transforms;
  p8est_edge_transform_t *et;
#endif
  int                 volume = -1;
  size_t              zz;
  int                 mpiret;
  size_t             *zp;
  sc_array_t         *sorter;
  int                *ip;
  int                 conntype = p4est_connect_type_int (ctype);
  int                 ncon = 1;

  P4EST_ASSERT (k >= 0);
  P4EST_ASSERT (newid != NULL);
  P4EST_ASSERT (newid->elem_size == sizeof (size_t));

  if (k == 0) {
    mpiret = sc_MPI_Comm_size (comm, &k);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Comm_rank (comm, &rank);
  SC_CHECK_MPI (mpiret);

  /* part will hold the partition number of each tree */
  part = P4EST_ALLOC (int, n);

  if (!rank) {

    xadj = P4EST_ALLOC (int, n + 1);

    switch (conntype) {
    case 1:
      degree = P4EST_FACES;
      break;
    case P4EST_DIM:
      degree = P4EST_INSUL - 1;
      sc_array_init (cta, sizeof (p4est_corner_transform_t));
#ifdef P4_TO_P8
      sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
      break;
#ifdef P4_TO_P8
    case 2:
      degree = P8EST_FACES + P8EST_EDGES;
      sc_array_init (eta, sizeof (p8est_edge_transform_t));
      break;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
    }

    if (degree == P4EST_FACES) {
      /* each tree has the same: metis shouldn't have any trouble with a
       * loop on a face/edge corner that has no neighbor */
      for (i = 0; i < n + 1; i++) {
        xadj[i] = P4EST_FACES * i;
      }
      adjncy = P4EST_ALLOC (int, P4EST_FACES * n);
      for (i = 0; i < n; i++) {
        for (j = 0; j < P4EST_FACES; j++) {
          adjncy[P4EST_FACES * i + j] =
            conn->tree_to_tree[P4EST_FACES * i + j];
        }
      }
    }
    else {
      totaldeg = 0;
      xadj[0] = 0;
      for (i = 0; i < n; i++) {
        totaldeg += P4EST_FACES;
        if (conntype == P4EST_DIM) {
          for (j = 0; j < P4EST_CHILDREN; j++) {
            /* add the number of strict corner neighbors */
            p4est_find_corner_transform (conn, (p4est_topidx_t) i, j, &ci);
            totaldeg += (int) cta->elem_count;
          }
        }
#ifdef P4_TO_P8
        if (conntype >= 2) {
          /* add the number of strict edge neighbors */
          for (j = 0; j < P8EST_EDGES; j++) {
            p8est_find_edge_transform (conn, (p4est_topidx_t) i, j, &ei);
            totaldeg += (int) eta->elem_count;
          }
        }
#endif
        xadj[i + 1] = totaldeg;
      }

      adjncy = P4EST_ALLOC (int, totaldeg);

      l = 0;
      for (i = 0; i < n; i++) {
        for (j = 0; j < P4EST_FACES; j++) {
          adjncy[l++] = (int) conn->tree_to_tree[P4EST_FACES * i + j];
        }
        if (conntype == P4EST_DIM) {
          for (j = 0; j < P4EST_CHILDREN; j++) {
            /* add the number of strict corner neighbors */
            p4est_find_corner_transform (conn, (p4est_topidx_t) i, j, &ci);
            for (zz = 0; zz < cta->elem_count; zz++) {
              ct = p4est_corner_array_index (cta, zz);
              adjncy[l++] = (int) ct->ntree;
            }
          }
        }
#ifdef P4_TO_P8
        if (conntype >= 2) {
          /* add the number of strict edge neighbors */
          for (j = 0; j < P8EST_EDGES; j++) {
            p8est_find_edge_transform (conn, (p4est_topidx_t) i, j, &ei);
            for (zz = 0; zz < eta->elem_count; zz++) {
              et = p8est_edge_array_index (eta, zz);
              adjncy[l++] = (int) et->ntree;
            }
          }
        }
#endif
        P4EST_ASSERT (l == xadj[i + 1]);
      }

      P4EST_ASSERT (l == totaldeg);

      if (conntype == P4EST_DIM) {
        sc_array_reset (cta);
      }
#ifdef P4_TO_P8
      if (conntype >= 2) {
        sc_array_reset (eta);
      }
#endif
    }

    P4EST_GLOBAL_INFO ("Entering metis\n");
    /* now call metis */
    metis_n = n;
    P4EST_EXECUTE_ASSERT_INT
      (METIS_PartGraphRecursive (&metis_n, &ncon, xadj, adjncy, NULL, NULL,
                                 NULL, &k, NULL, NULL, NULL, &volume, part),
       METIS_OK);
    P4EST_ASSERT (metis_n == n);
    P4EST_GLOBAL_INFO ("Done metis\n");

    P4EST_GLOBAL_STATISTICSF ("metis volume %d\n", volume);

    P4EST_FREE (xadj);
    P4EST_FREE (adjncy);
  }

  /* broadcast part to every process: this is expensive, should probably think
   * of a better way to do this */
  sc_MPI_Bcast (part, n, sc_MPI_INT, 0, comm);

  /* now that everyone has part, each process computes the renumbering
   * for itself*/
  sc_array_resize (newid, (size_t) n);
  sorter = sc_array_new_size (2 * sizeof (int), (size_t) n);
  for (i = 0; i < n; i++) {
    ip = (int *) sc_array_index (sorter, i);
    ip[0] = part[i];
    ip[1] = i;
  }
  P4EST_FREE (part);

  /* sort current index by partition given */
  /* this will be the same on every process because the comparison operation
   * does not allow equality between different trees */
  sc_array_sort (sorter, reorder_comp);
  for (i = 0; i < n; i++) {
    ip = (int *) sc_array_index (sorter, i);
    zp = (size_t *) sc_array_index (newid, ip[1]);
    *zp = i;
  }
  sc_array_destroy (sorter);

  p4est_connectivity_permute (conn, newid, 1);

}

#endif /* P4EST_WITH_METIS */

static int
p4est_topidx_compare_2 (const void *A, const void *B)
{
  int                 ret = p4est_topidx_compare (A, B);

  if (!ret) {
    const p4est_topidx_t *a = (const p4est_topidx_t *) A;
    const p4est_topidx_t *b = (const p4est_topidx_t *) B;
    p4est_topidx_t      diff = a[1] - b[1];

    ret = diff ? (diff < 0 ? -1 : 1) : 0;
  }
  return ret;
}

static void
p4est_connectivity_store_corner (p4est_connectivity_t * conn,
                                 p4est_topidx_t t, int c)
{
  p4est_topidx_t      n = ++conn->num_corners;
  p4est_topidx_t     *tc;
  size_t              zz;
  size_t              zk;
  int                 i;
  sc_array_t         *corner_to_tc;

  P4EST_ASSERT (conn->tree_to_corner == NULL ||
                conn->tree_to_corner[P4EST_CHILDREN * t + c] < 0);

  conn->ctt_offset = P4EST_REALLOC (conn->ctt_offset, p4est_topidx_t, n + 1);
  conn->ctt_offset[n] = conn->ctt_offset[n - 1];

  if (conn->tree_to_corner == NULL) {
    conn->tree_to_corner =
      P4EST_ALLOC (p4est_topidx_t, P4EST_CHILDREN * conn->num_trees);
    memset (conn->tree_to_corner, -1,
            P4EST_CHILDREN * conn->num_trees * sizeof (p4est_topidx_t));
  }

  corner_to_tc = sc_array_new (2 * sizeof (p4est_topidx_t));

  conn->tree_to_corner[P4EST_CHILDREN * t + c] = n - 1;
  tc = (p4est_topidx_t *) sc_array_push (corner_to_tc);
  tc[0] = t;
  tc[1] = c;

  for (i = 0; i < P4EST_DIM; i++) {
    int                 f = p4est_corner_faces[c][i];
    p4est_topidx_t      nt = conn->tree_to_tree[P4EST_FACES * t + f];
    int                 nf = conn->tree_to_face[P4EST_FACES * t + f];
    int                 o;
    int                 nc;

    o = nf / P4EST_FACES;
    nf %= P4EST_FACES;

    if (nt == t && nf == f) {
      continue;
    }

    nc = p4est_connectivity_face_neighbor_corner (c, f, nf, o);

    conn->tree_to_corner[P4EST_CHILDREN * nt + nc] = n - 1;
    tc = (p4est_topidx_t *) sc_array_push (corner_to_tc);
    tc[0] = nt;
    tc[1] = nc;
  }
#ifdef P4_TO_P8
  for (i = 0; i < P4EST_DIM; i++) {
    p8est_edge_info_t   ei;
    p8est_edge_transform_t *et;
    int                 e = p8est_corner_edges[c][i];

    sc_array_init (&(ei.edge_transforms), sizeof (p8est_edge_transform_t));
    p8est_find_edge_transform (conn, t, e, &ei);

    for (zz = 0; zz < ei.edge_transforms.elem_count; zz++) {
      p4est_topidx_t      nt;
      int                 ne;
      int                 nc;

      et = p8est_edge_array_index (&(ei.edge_transforms), zz);

      nt = et->ntree;
      ne = et->nedge;
      if (p8est_edge_corners[e][0] == c) {
        nc = p8est_edge_corners[ne][et->nflip];
      }
      else {
        nc = p8est_edge_corners[ne][1 ^ et->nflip];
      }

      conn->tree_to_corner[P4EST_CHILDREN * nt + nc] = n - 1;
      tc = (p4est_topidx_t *) sc_array_push (corner_to_tc);
      tc[0] = nt;
      tc[1] = nc;
    }

    sc_array_reset (&(ei.edge_transforms));
  }
#endif

  sc_array_sort (corner_to_tc, p4est_topidx_compare_2);
  sc_array_uniq (corner_to_tc, p4est_topidx_compare_2);

  zk = corner_to_tc->elem_count;
  conn->ctt_offset[n] += (p4est_topidx_t) zk;
  conn->corner_to_tree = P4EST_REALLOC (conn->corner_to_tree,
                                        p4est_topidx_t, conn->ctt_offset[n]);
  conn->corner_to_corner = P4EST_REALLOC (conn->corner_to_corner,
                                          int8_t, conn->ctt_offset[n]);

  for (zz = 0; zz < zk; zz++) {
    tc = (p4est_topidx_t *) sc_array_index (corner_to_tc, zz);
    conn->corner_to_tree[conn->ctt_offset[n - 1] + zz] = tc[0];
    conn->corner_to_corner[conn->ctt_offset[n - 1] + zz] = (int8_t) tc[1];
  }

  sc_array_destroy (corner_to_tc);
}

#ifdef P4_TO_P8
static void
p8est_connectivity_store_edge (p4est_connectivity_t * conn, p4est_topidx_t t,
                               int e)
{
  p4est_topidx_t      n = ++conn->num_edges;
  p4est_topidx_t     *te;
  size_t              zz;
  size_t              zk;
  int                 i;
  sc_array_t         *edge_to_te;

  P4EST_ASSERT (conn->tree_to_edge == NULL ||
                conn->tree_to_edge[P8EST_EDGES * t + e] < 0);

  conn->ett_offset = P4EST_REALLOC (conn->ett_offset, p4est_topidx_t, n + 1);
  conn->ett_offset[n] = conn->ett_offset[n - 1];

  if (conn->tree_to_edge == NULL) {
    conn->tree_to_edge =
      P4EST_ALLOC (p4est_topidx_t, P8EST_EDGES * conn->num_trees);
    memset (conn->tree_to_edge, -1,
            P8EST_EDGES * conn->num_trees * sizeof (p4est_topidx_t));
  }

  edge_to_te = sc_array_new (2 * sizeof (p4est_topidx_t));

  conn->tree_to_edge[P8EST_EDGES * t + e] = n - 1;
  te = (p4est_topidx_t *) sc_array_push (edge_to_te);
  te[0] = t;
  te[1] = e;

  for (i = 0; i < 2; i++) {
    int                 f = p8est_edge_faces[e][i];
    p4est_topidx_t      nt = conn->tree_to_tree[P4EST_FACES * t + f];
    int                 nf = conn->tree_to_face[P4EST_FACES * t + f];
    int                 o, c[2], nc[2];
    int                 ne;
    int                 ref;
    int                 set;
    int                 j;
    int                 diff;

    o = nf / P4EST_FACES;
    nf %= P4EST_FACES;

    if (t == nt && f == nf) {
      continue;
    }
    ref = p8est_face_permutation_refs[f][nf];
    set = p8est_face_permutation_sets[ref][o];

    for (j = 0; j < 2; j++) {
      c[j] = p8est_edge_corners[e][j];
      nc[j] = p8est_connectivity_face_neighbor_corner_set (c[j], f, nf, set);
    }
    diff = SC_MAX (nc[0], nc[1]) - SC_MIN (nc[0], nc[1]);
    switch (diff) {
    case 1:
      ne = p8est_corner_edges[nc[0]][0];
      break;
    case 2:
      ne = p8est_corner_edges[nc[0]][1];
      break;
    case 4:
      ne = p8est_corner_edges[nc[0]][2];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    conn->tree_to_edge[P8EST_EDGES * nt + ne] = n - 1;
    if (p8est_edge_corners[ne][0] != nc[0]) {
      ne += 12;
    }

    te = (p4est_topidx_t *) sc_array_push (edge_to_te);
    te[0] = nt;
    te[1] = ne;
  }

  sc_array_sort (edge_to_te, p4est_topidx_compare_2);
  sc_array_uniq (edge_to_te, p4est_topidx_compare_2);

  zk = edge_to_te->elem_count;
  conn->ett_offset[n] += (p4est_topidx_t) zk;
  conn->edge_to_tree = P4EST_REALLOC (conn->edge_to_tree,
                                      p4est_topidx_t, conn->ett_offset[n]);
  conn->edge_to_edge = P4EST_REALLOC (conn->edge_to_edge,
                                      int8_t, conn->ett_offset[n]);

  for (zz = 0; zz < zk; zz++) {
    te = (p4est_topidx_t *) sc_array_index (edge_to_te, zz);
    conn->edge_to_tree[conn->ett_offset[n - 1] + zz] = te[0];
    conn->edge_to_edge[conn->ett_offset[n - 1] + zz] = (int8_t) te[1];
  }

  sc_array_destroy (edge_to_te);
}
#endif

static void
p4est_connectivity_join_corners (p4est_connectivity_t * conn,
                                 p4est_topidx_t tree_left,
                                 p4est_topidx_t tree_right,
                                 int corner_left, int corner_right)
{
  p4est_topidx_t      c, c0, c1, swap;
  p4est_topidx_t      startt, endt, n1, it, end0;
  p4est_topidx_t     *swapspace;
  int8_t             *swapspacei;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (tree_left >= 0 && tree_left < conn->num_trees);
  P4EST_ASSERT (tree_right >= 0 && tree_right < conn->num_trees);
  P4EST_ASSERT (corner_left >= 0 && corner_left < P4EST_CHILDREN);
  P4EST_ASSERT (corner_right >= 0 && corner_right < P4EST_CHILDREN);

  /* it could be that the current connectivity did not store corner information,
   * because all of the corners are simple enough that they can be figured out
   * from context.  To simplify things, we're going to only deal with
   * explicitly stored corners. */
  if (conn->tree_to_corner == NULL ||
      conn->tree_to_corner[P4EST_CHILDREN * tree_left + corner_left] < 0) {
    p4est_connectivity_store_corner (conn, tree_left, corner_left);
  }
  if (conn->tree_to_corner == NULL ||
      conn->tree_to_corner[P4EST_CHILDREN * tree_right + corner_right] < 0) {
    p4est_connectivity_store_corner (conn, tree_right, corner_right);
  }

  /* now we know that the two corners are explicitly stored, so it's just a
   * matter of combining their storage and removing references to one of them
   * */
  c0 = conn->tree_to_corner[P4EST_CHILDREN * tree_left + corner_left];
  c1 = conn->tree_to_corner[P4EST_CHILDREN * tree_right + corner_right];

  if (c0 == c1) {
    /* whoops, these two corners were already the same */
    return;
  }
  if (c1 < c0) {
    swap = c0;
    c0 = c1;
    c1 = swap;
  }

  /* remove all reference to c1 */
  startt = conn->ctt_offset[c1];
  endt = conn->ctt_offset[c1 + 1];

  n1 = endt - startt;           /* the number of tree corners that border c1 */
  for (it = startt; it < endt; it++) {  /* get all trees that reference c1 */
    p4est_topidx_t      nt = conn->corner_to_tree[it];  /* nt is a tree the borders c1 */
    int                 ntc = (int) conn->corner_to_corner[it]; /* ntc is nt's numbering for c1 */

    conn->tree_to_corner[P4EST_CHILDREN * nt + ntc] = c0;       /* c1->c0 */
  }

  /* we now have to move the entries in corner_to_tree and corner_to_corner around */
  end0 = conn->ctt_offset[c0 + 1];

  swapspace = P4EST_ALLOC (p4est_topidx_t, n1);
  memcpy (swapspace, conn->corner_to_tree + (size_t) startt,
          n1 * sizeof (p4est_topidx_t));
  memmove (conn->corner_to_tree + (size_t) (end0 + n1),
           conn->corner_to_tree + (size_t) end0,
           (size_t) (startt - end0) * sizeof (p4est_topidx_t));
  memcpy (conn->corner_to_tree + (size_t) end0, swapspace,
          n1 * sizeof (p4est_topidx_t));
  P4EST_FREE (swapspace);

  swapspacei = P4EST_ALLOC (int8_t, n1);
  memcpy (swapspacei, conn->corner_to_corner + (size_t) startt,
          n1 * sizeof (int8_t));
  memmove (conn->corner_to_corner + (size_t) (end0 + n1),
           conn->corner_to_corner + (size_t) end0,
           (size_t) (startt - end0) * sizeof (int8_t));
  memcpy (conn->corner_to_corner + (size_t) end0, swapspacei,
          n1 * sizeof (int8_t));
  P4EST_FREE (swapspacei);

  /* finally, we have to correct ctt_offset */
  for (c = c0 + 1; c <= c1; c++) {
    conn->ctt_offset[c] += n1;
  }

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
}

#ifdef P4_TO_P8
static void
p8est_connectivity_join_edges (p8est_connectivity_t * conn,
                               p4est_topidx_t tree_left,
                               p4est_topidx_t tree_right,
                               int edge_left, int edge_right, int orientation)
{
  int                 i, c_left, c_right;
  p4est_topidx_t      e, e0, e1, swap;
  p4est_topidx_t      startt, endt, n1, it, end0;
  p4est_topidx_t     *swapspace;
  int8_t             *swapspacei;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (tree_left >= 0 && tree_left < conn->num_trees);
  P4EST_ASSERT (tree_right >= 0 && tree_right < conn->num_trees);
  P4EST_ASSERT (edge_left >= 0 && edge_left < P8EST_EDGES);
  P4EST_ASSERT (edge_right >= 0 && edge_right < P8EST_EDGES);
  for (i = 0; i < 2; i++) {
    /* get matching corners */
    c_left = p8est_edge_corners[edge_left][i];
    if (orientation) {
      c_right = p8est_edge_corners[edge_right][1 ^ i];
    }
    else {
      c_right = p8est_edge_corners[edge_right][i];
    }
    /* join the corners */
    p4est_connectivity_join_corners (conn, tree_left, tree_right,
                                     c_left, c_right);
  }

  /* it could be that the current connectivity did not store edge information,
   * because all of the edges are simple enough that they can be figured out
   * from context.  To simplify things, we're going to only deal with
   * explicitly stored edges. */
  if (conn->tree_to_edge == NULL ||
      conn->tree_to_edge[P8EST_EDGES * tree_left + edge_left] < 0) {
    p8est_connectivity_store_edge (conn, tree_left, edge_left);
  }
  if (conn->tree_to_edge == NULL ||
      conn->tree_to_edge[P8EST_EDGES * tree_right + edge_right] < 0) {
    p8est_connectivity_store_edge (conn, tree_right, edge_right);
  }

  /* now we know that the two edges are explicitly stored, so it's just a
   * matter of combining their storage and removing references to one of them
   * */
  e0 = conn->tree_to_edge[P8EST_EDGES * tree_left + edge_left];
  e1 = conn->tree_to_edge[P8EST_EDGES * tree_right + edge_right];

  if (e0 == e1) {
    /* whoops, these two edges were already the same, looks like we did a bunch
     * of work for nothing */
    return;
  }
  if (e1 < e0) {
    swap = e0;
    e0 = e1;
    e1 = swap;
  }

  /* remove all reference to e1 */
  startt = conn->ett_offset[e1];
  endt = conn->ett_offset[e1 + 1];

  n1 = endt - startt;           /* the number of tree edges that border e1 */
  for (it = startt; it < endt; it++) {  /* get all trees that reference e1 */
    p4est_topidx_t      nt = conn->edge_to_tree[it];    /* nt is a tree the borders e1 */
    int                 nte = (int) conn->edge_to_edge[it];     /* nte is nt's numbering for e1,
                                                                   modified by orientation */
    int                 o = nte / P8EST_EDGES;  /* o is that modifying orientation */

    nte %= P8EST_EDGES;         /* okay, now nte is nt's numbering for e1 */
    conn->tree_to_edge[P8EST_EDGES * nt + nte] = e0;    /* e1->e0 */
    /* if edge_left and edge_right have opposite orientations, then the
     * orientation information in edge_to_edge has to be toggled */
    conn->edge_to_edge[it] = P8EST_EDGES * (o ^ orientation) + nte;
  }

  /* we now have to move the entries in edge_to_tree and edge_to_edge around */
  end0 = conn->ett_offset[e0 + 1];

  swapspace = P4EST_ALLOC (p4est_topidx_t, n1);
  memcpy (swapspace, conn->edge_to_tree + (size_t) startt,
          n1 * sizeof (p4est_topidx_t));
  memmove (conn->edge_to_tree + (size_t) (end0 + n1),
           conn->edge_to_tree + (size_t) end0,
           (size_t) (startt - end0) * sizeof (p4est_topidx_t));
  memcpy (conn->edge_to_tree + (size_t) end0, swapspace,
          n1 * sizeof (p4est_topidx_t));
  P4EST_FREE (swapspace);

  swapspacei = P4EST_ALLOC (int8_t, n1);
  memcpy (swapspacei, conn->edge_to_edge + (size_t) startt,
          n1 * sizeof (int8_t));
  memmove (conn->edge_to_edge + (size_t) (end0 + n1),
           conn->edge_to_edge + (size_t) end0,
           (size_t) (startt - end0) * sizeof (int8_t));
  memcpy (conn->edge_to_edge + (size_t) end0, swapspacei,
          n1 * sizeof (int8_t));
  P4EST_FREE (swapspacei);

  /* finally, we have to correct ett_offset */
  for (e = e0 + 1; e <= e1; e++) {
    conn->ett_offset[e] += n1;
  }

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
}
#endif

void
p4est_connectivity_join_faces (p4est_connectivity_t * conn,
                               p4est_topidx_t tree_left,
                               p4est_topidx_t tree_right,
                               int face_left, int face_right, int orientation)
{
#ifdef P4_TO_P8
  int                 ref, set, j;
#endif
  int                 i;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  P4EST_ASSERT (tree_left >= 0 && tree_left < conn->num_trees);
  P4EST_ASSERT (tree_right >= 0 && tree_right < conn->num_trees);
  P4EST_ASSERT (face_left >= 0 && face_left < P4EST_FACES);
  P4EST_ASSERT (face_right >= 0 && face_right < P4EST_FACES);
  P4EST_ASSERT (orientation >= 0 && orientation < P4EST_HALF);
  P4EST_ASSERT (conn->tree_to_tree[P4EST_FACES * tree_left + face_left] ==
                tree_left);
  P4EST_ASSERT (conn->tree_to_tree[P4EST_FACES * tree_right + face_right] ==
                tree_right);
  P4EST_ASSERT (conn->tree_to_face[P4EST_FACES * tree_left + face_left] ==
                (int8_t) face_left);
  P4EST_ASSERT (conn->tree_to_face[P4EST_FACES * tree_right + face_right] ==
                (int8_t) face_right);

#ifdef P4_TO_P8
  /* figure out which edges are next to each other */
  ref = p8est_face_permutation_refs[face_left][face_right];
  set = p8est_face_permutation_sets[ref][orientation];
  for (i = 0; i < 4; i++) {
    int                 c[2], e_left, e_right, e_orient;

    /* get an edge of face_left */
    e_left = p8est_face_edges[face_left][i];

    for (j = 0; j < 2; j++) {
      /* get corners of that edge and their numbers seen from face_right */
      c[j] = p8est_connectivity_face_neighbor_corner_set
        (p8est_edge_corners[e_left][j], face_left, face_right, set);
    }
    /* now from the two corners, we can figure out e_right */
    e_right = p8est_child_corner_edges[c[0]][c[1]];
    P4EST_ASSERT (e_right >= 0);

    /* how are e_left and e_right oriented? 0 for same orientation, 1 for
     * opposite */
    e_orient = (p8est_edge_corners[e_right][0] == c[1]);

    /* now we have two facing edges and their orientation, so we can join them */
    /* these routines will also join the corners */
    p8est_connectivity_join_edges (conn, tree_left, tree_right, e_left,
                                   e_right, e_orient);
  }
#else
  for (i = 0; i < 2; i++) {
    int                 c_left, c_right;

    c_left = p4est_face_corners[face_left][i];
    if (orientation) {          /* if the two faces have opposite orientation */
      c_right = p4est_face_corners[face_right][1 ^ i];
    }
    else {                      /* if the two faces have the same orientation */
      c_right = p4est_face_corners[face_right][i];
    }

    /* join the corners */
    p4est_connectivity_join_corners (conn, tree_left, tree_right, c_left,
                                     c_right);
  }
#endif

  conn->tree_to_tree[P4EST_FACES * tree_left + face_left] = tree_right;
  conn->tree_to_tree[P4EST_FACES * tree_right + face_right] = tree_left;
  conn->tree_to_face[P4EST_FACES * tree_left + face_left] =
    face_right + P4EST_FACES * orientation;
  conn->tree_to_face[P4EST_FACES * tree_right + face_right] =
    face_left + P4EST_FACES * orientation;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
}

#ifdef P4_TO_P8
static int
p8est_edge_compare (const void *a, const void *b)
{
  const p8est_edge_transform_t *A = (const p8est_edge_transform_t *) a;
  const p8est_edge_transform_t *B = (const p8est_edge_transform_t *) b;

  return (A->ntree != B->ntree) ? A->ntree - B->ntree :
    (A->nedge != B->nedge) ? A->nedge - B->nedge :
    (A->naxis[0] != B->naxis[0]) ? A->naxis[0] - B->naxis[0] :
    (A->naxis[1] != B->naxis[1]) ? A->naxis[1] - B->naxis[1] :
    (A->naxis[2] != B->naxis[2]) ? A->naxis[2] - B->naxis[2] :
    (A->nflip != B->nflip) ? A->nflip - B->nflip : A->corners - B->corners;
}
#endif

static int
p4est_corner_compare (const void *a, const void *b)
{
  const p4est_corner_transform_t *A = (const p4est_corner_transform_t *) a;
  const p4est_corner_transform_t *B = (const p4est_corner_transform_t *) b;

  return (A->ntree != B->ntree) ? A->ntree - B->ntree :
    A->ncorner - B->ncorner;
}

int
p4est_connectivity_is_equivalent (p4est_connectivity_t * conn1,
                                  p4est_connectivity_t * conn2)
{
  const size_t        topsize = sizeof (p4est_topidx_t);
  const size_t        int8size = sizeof (int8_t);
  size_t              count;
  p4est_topidx_t      ntrees, t;

  P4EST_ASSERT (p4est_connectivity_is_valid (conn1));
  P4EST_ASSERT (p4est_connectivity_is_valid (conn2));

  /* same pointer or equality are stronger */
  if (conn1 == conn2 || p4est_connectivity_is_equal (conn1, conn2)) {
    return 1;
  }

  ntrees = conn1->num_trees;

  /* clearly must have same number of trees */
  if (conn2->num_trees != ntrees) {
    return 0;
  }

  /* compare tree_to_tree, tree_to_face structure: must be exactly the same */
  count = (size_t) (P4EST_FACES * conn1->num_trees);
  if (memcmp (conn1->tree_to_tree, conn2->tree_to_tree, count * topsize) ||
      memcmp (conn1->tree_to_face, conn2->tree_to_face, count * int8size)) {
    return 0;
  }

  /* test equivalence of edges and corners only through the transforms: the
   * numbering of the edges and corners is not relevant */

#ifdef P4_TO_P8
  {
    p8est_edge_info_t   e1, e2;

    sc_array_init (&e1.edge_transforms, sizeof (p8est_edge_transform_t));
    sc_array_init (&e2.edge_transforms, sizeof (p8est_edge_transform_t));
    for (t = 0; t < ntrees; t++) {
      int                 e;
      size_t              zz;

      for (e = 0; e < P8EST_EDGES; e++) {
        p8est_find_edge_transform (conn1, t, e, &e1);
        p8est_find_edge_transform (conn2, t, e, &e2);
        if (e1.edge_transforms.elem_count != e2.edge_transforms.elem_count) {
          return 0;
        }
        /* sort so memory comparison is correct */
        sc_array_sort (&e1.edge_transforms, p8est_edge_compare);
        sc_array_sort (&e2.edge_transforms, p8est_edge_compare);
        if (e1.edge_transforms.elem_count != e2.edge_transforms.elem_count) {
          return 0;
        }
        for (zz = 0; zz < e1.edge_transforms.elem_count; zz++) {
          p8est_edge_transform_t *t1 = p8est_edge_array_index
            (&e1.edge_transforms, zz);
          p8est_edge_transform_t *t2 = p8est_edge_array_index
            (&e2.edge_transforms, zz);

          if (t1->corners != t2->corners ||
              t1->naxis[0] != t2->naxis[0] ||
              t1->naxis[1] != t2->naxis[1] ||
              t1->naxis[2] != t2->naxis[2] ||
              t1->nedge != t2->nedge ||
              t1->nflip != t2->nflip || t1->ntree != t2->ntree) {
            return 0;
          }
        }
      }
    }
    sc_array_reset (&e1.edge_transforms);
    sc_array_reset (&e2.edge_transforms);
  }
#endif
  {
    p4est_corner_info_t c1, c2;

    sc_array_init (&c1.corner_transforms, sizeof (p4est_corner_transform_t));
    sc_array_init (&c2.corner_transforms, sizeof (p4est_corner_transform_t));
    for (t = 0; t < ntrees; t++) {
      int                 c;
      size_t              zz;

      for (c = 0; c < P4EST_CHILDREN; c++) {
        p4est_find_corner_transform (conn1, t, c, &c1);
        p4est_find_corner_transform (conn2, t, c, &c2);
        if (c1.corner_transforms.elem_count !=
            c2.corner_transforms.elem_count) {
          return 0;
        }
        /* sort so memory comparison is correct */
        sc_array_sort (&c1.corner_transforms, p4est_corner_compare);
        sc_array_sort (&c2.corner_transforms, p4est_corner_compare);

        if (c1.corner_transforms.elem_count !=
            c2.corner_transforms.elem_count) {
          return 0;
        }
        for (zz = 0; zz < c1.corner_transforms.elem_count; zz++) {
          p4est_corner_transform_t *t1 = p4est_corner_array_index
            (&c1.corner_transforms, zz);
          p4est_corner_transform_t *t2 = p4est_corner_array_index
            (&c2.corner_transforms, zz);

          if (t1->ncorner != t2->ncorner || t1->ntree != t2->ntree) {
            return 0;
          }
        }
      }
    }
    sc_array_reset (&c1.corner_transforms);
    sc_array_reset (&c2.corner_transforms);
  }
  return 1;
}

#ifndef P4_TO_P8

int
p4est_connect_type_int (p4est_connect_type_t btype)
{
  switch (btype) {
  case P4EST_CONNECT_FACE:
    return 1;
  case P4EST_CONNECT_CORNER:
    return 2;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

const char         *
p4est_connect_type_string (p4est_connect_type_t btype)
{
  switch (btype) {
  case P4EST_CONNECT_FACE:
    return "FACE";
  case P4EST_CONNECT_CORNER:
    return "CORNER";
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

#endif /* !P4_TO_P8 */
/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
static char        *
p4est_connectivity_getline_upper (FILE * stream)
{
  char               *line = P4EST_ALLOC (char, 1024), *linep = line;
  size_t              lenmax = 1024, len = lenmax;
  int                 c;

  if (line == NULL)
    return NULL;

  for (;;) {
    c = fgetc (stream);
    if (c == EOF && linep == line) {
      P4EST_FREE (linep);
      return NULL;
    }
    c = toupper (c);

    if (--len == 0) {
      char               *linen;

      len = lenmax;
      lenmax *= 2;

      linen = P4EST_REALLOC (linep, char, lenmax);
      if (linen == NULL) {
        P4EST_FREE (linep);
        return NULL;
      }

      line = linen + (line - linep);
      linep = linen;
    }
    if ((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

int
p4est_connectivity_read_inp_stream (FILE * stream,
                                    p4est_topidx_t * num_vertices,
                                    p4est_topidx_t * num_trees,
                                    double *vertices,
                                    p4est_topidx_t * tree_to_vertex)
{
  int                 reading_nodes = 0, reading_elements = 0;
  int                 lines_read = 0, lines_free = 0;
  char               *line;
  p4est_topidx_t      num_nodes = 0;
  p4est_topidx_t      num_elements = 0;
  int                 fill_trees_and_vertices = (vertices != NULL &&
                                                 tree_to_vertex != NULL);

  P4EST_ASSERT ((vertices == NULL && tree_to_vertex == NULL) ||
                (vertices != NULL && tree_to_vertex != NULL));

  for (;;) {
    line = p4est_connectivity_getline_upper (stream);

    if (line == NULL) {
      break;
    }

    ++lines_read;

    /* check for control line */
    if (line[0] == '*') {
      reading_elements = reading_nodes = 0;
      if (strstr (line, "*NODE")) {
        reading_nodes = 1;
        ++lines_free;
        P4EST_FREE (line);
        continue;
      }
      else if (strstr (line, "*ELEMENT")) {
        if (
#ifdef P4_TO_P8
             strstr (line, "TYPE=C3D8")
#else
             strstr (line, "TYPE=C2D4") || strstr (line, "TYPE=CPS4")
             || strstr (line, "TYPE=S4")
#endif
          ) {
          reading_elements = 1;
          ++lines_free;
          P4EST_FREE (line);
          continue;
        }
      }
    }

    if (reading_nodes) {
      if (fill_trees_and_vertices) {
        long long int       node;
        double              x, y, z;
        int                 retval;

        retval = sscanf (line, "%lld, %lf, %lf, %lf", &node, &x, &y, &z);
        if (retval != 4) {
          P4EST_LERROR ("Premature end of file");
          P4EST_FREE (line);
          return 1;
        }

        if (node > *num_vertices) {
          P4EST_LERRORF
            ("Encountered vertex %lld that will not fit in vertices"
             " array of length %lld.  Are the vertices contiguously"
             " numbered?\n", node, (long long int) *num_vertices);
          P4EST_FREE (line);
          return 1;
        }

        vertices[3 * (node - 1) + 0] = x;
        vertices[3 * (node - 1) + 1] = y;
        vertices[3 * (node - 1) + 2] = z;
      }

      ++num_nodes;
    }
    else if (reading_elements) {
      if (fill_trees_and_vertices) {
        long long int       v[P4EST_CHILDREN];
        int                 n;
        int                 retval;

        if (num_elements >= *num_trees) {
          P4EST_LERROR ("Encountered element that will not fit into"
                        " tree_to_vertex array. More elements than expected.\n");
          P4EST_FREE (line);
          return 1;
        }

        /* Note that when we read in the
         * vertices we switch from right-hand
         * vertex ordering to z-order
         */
        retval = sscanf (line, "%*d, %lld, %lld, %lld, %lld"
#ifdef P4_TO_P8
                         ", %lld, %lld, %lld, %lld"
#endif
                         , &v[0], &v[1], &v[3], &v[2]
#ifdef P4_TO_P8
                         , &v[4], &v[5], &v[7], &v[6]
#endif
          );
        if (retval != P4EST_CHILDREN) {
          P4EST_LERROR ("Premature end of file");
          P4EST_FREE (line);
          return 1;
        }

        for (n = 0; n < P4EST_CHILDREN; ++n)
          tree_to_vertex[P4EST_CHILDREN * num_elements + n] = v[n] - 1;
      }

      ++num_elements;
    }

    ++lines_free;
    P4EST_FREE (line);
  }

  *num_vertices = num_nodes;
  *num_trees = num_elements;

  if (num_nodes == 0 || num_elements == 0) {
    P4EST_LERROR ("No elements or nodes found in mesh file.\n");
    return -1;
  }
  else {
    return 0;
  }
}

p4est_connectivity_t *
p4est_connectivity_read_inp (const char *filename)
{
  int                 retval;
  p4est_topidx_t      num_vertices = 0, num_trees = 0, tree;
  int                 face;

  p4est_connectivity_t *conn = NULL;
  FILE               *fid = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Reading connectivity from %s\n", filename);

  fid = fopen (filename, "rb");
  if (fid == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", filename);
    goto dead;
  }

  if (p4est_connectivity_read_inp_stream
      (fid, &num_vertices, &num_trees, NULL, NULL)) {
    P4EST_LERRORF ("Failed to read %s: pass 1\n", filename);
    goto dead;
  }

  rewind (fid);

  conn = p4est_connectivity_new (num_vertices, num_trees,
#ifdef P4_TO_P8
                                 0, 0,
#endif
                                 0, 0);

  if (p4est_connectivity_read_inp_stream (fid, &conn->num_vertices,
                                          &conn->num_trees, conn->vertices,
                                          conn->tree_to_vertex)) {
    P4EST_LERRORF ("Failed to read %s: pass 2\n", filename);
    goto dead;
  }

  /*
   * Fill tree_to_tree and tree_to_face to make sure we have a valid
   * connectivity.
   */
  for (tree = 0; tree < conn->num_trees; ++tree) {
    for (face = 0; face < P4EST_FACES; ++face) {
      conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
      conn->tree_to_face[P4EST_FACES * tree + face] = face;
    }
  }
  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  /* Compute real tree_to_* fields and complete (edge and) corner fields. */
  p4est_connectivity_complete (conn);

  retval = fclose (fid);
  fid = NULL;
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", filename);
    goto dead;
  }

  P4EST_GLOBAL_PRODUCTIONF
    ("New connectivity with %lld trees and %lld vertices\n",
     (long long) conn->num_trees, (long long) conn->num_vertices);

  return conn;

dead:
  /* clean up on error */
  if (fid != NULL) {
    fclose (fid);
  }
  if (conn != NULL) {
    p4est_connectivity_destroy (conn);
  }
  return NULL;
}

/* *INDENT-OFF* */
static p4est_neighbor_transform_t *p4est_neighbor_transform_array_push
  (sc_array_t *array)
{
  return (p4est_neighbor_transform_t *) sc_array_push (array);
}
/* *INDENT-ON* */

static void
p4est_face_transform_to_neighbor_transform (const int ftransform[9],
                                            p4est_neighbor_transform_t * nt)
{
  const int          *my_axis = &ftransform[0];
  const int          *target_axis = &ftransform[3];
  const int          *edge_reverse = &ftransform[6];
#ifndef P4_TO_P8
  int                 ids[] = { 0, 2 };
#else
  int                 ids[] = { 0, 1, 2 };
#endif
  int8_t              sign2;
  p4est_qcoord_t      o_self2, o_neigh2;

  for (int di = 0; di < P4EST_DIM; di++) {
    int                 d = ids[di];

    nt->perm[target_axis[d]] = my_axis[d];
  }
  for (int d = 0; d < P4EST_DIM - 1; d++) {
    nt->sign[target_axis[d]] = edge_reverse[d] ? -1 : 1;
    nt->origin_neighbor[target_axis[d]] = P4EST_ROOT_LEN / 2;
    nt->origin_self[my_axis[d]] = P4EST_ROOT_LEN / 2;
  }
  switch (edge_reverse[2]) {
  case 0:
    sign2 = -1;
    o_self2 = 0;
    o_neigh2 = 0;
    break;
  case 1:
    sign2 = 1;
    o_self2 = 0;
    o_neigh2 = P4EST_ROOT_LEN;
    break;
  case 2:
    sign2 = 1;
    o_self2 = P4EST_ROOT_LEN;
    o_neigh2 = 0;
    break;
  case 3:
    sign2 = -1;
    o_self2 = P4EST_ROOT_LEN;
    o_neigh2 = P4EST_ROOT_LEN;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  nt->sign[target_axis[2]] = sign2;
  nt->origin_self[my_axis[2]] = o_self2;
  nt->origin_neighbor[target_axis[2]] = o_neigh2;
}

#ifdef P4_TO_P8
static void
p8est_edge_transform_to_neighbor_transform (const p8est_edge_transform_t * et,
                                            int8_t iedge,
                                            p4est_neighbor_transform_t * nt)
{
  const int           other_axes[3][2] = { {1, 2}, {0, 2}, {0, 1} };
  int                 iaxis = iedge / 4;
  int                 naxis = et->naxis[0];

  nt->perm[naxis] = iaxis;
  nt->perm[other_axes[naxis][0]] = other_axes[iaxis][0];
  nt->perm[other_axes[naxis][1]] = other_axes[iaxis][1];

  nt->origin_self[iaxis] = P4EST_ROOT_LEN / 2;
  nt->origin_self[other_axes[iaxis][0]] = (iedge & 1) ? P4EST_ROOT_LEN : 0;
  nt->origin_self[other_axes[iaxis][1]] = (iedge & 2) ? P4EST_ROOT_LEN : 0;

  nt->origin_neighbor[naxis] = P4EST_ROOT_LEN / 2;
  nt->origin_neighbor[other_axes[naxis][0]] =
    (et->corners & 1) ? P4EST_ROOT_LEN : 0;
  nt->origin_neighbor[other_axes[naxis][1]] =
    (et->corners & 2) ? P4EST_ROOT_LEN : 0;

  nt->sign[naxis] = et->nflip ? -1 : 1;
  nt->sign[other_axes[naxis][0]] = ((iedge ^ et->corners) & 1) ? 1 : -1;
  nt->sign[other_axes[naxis][1]] = ((iedge ^ et->corners) & 2) ? 1 : -1;
}
#endif

static void
p4est_corner_transform_to_neighbor_transform (p4est_corner_transform_t * ct,
                                              int corner,
                                              p4est_neighbor_transform_t * nt)
{
  for (int d = 0; d < P4EST_DIM; d++) {
    nt->perm[d] = d;
    nt->origin_self[d] = (corner & (1 << d)) ? P4EST_ROOT_LEN : 0;
    nt->origin_neighbor[d] = (ct->ncorner & (1 << d)) ? P4EST_ROOT_LEN : 0;
    nt->sign[d] = ((corner ^ ct->ncorner) & (1 << d)) ? 1 : -1;
  }
}

void
p4est_connectivity_get_neighbor_transforms (p4est_connectivity_t * conn,
                                            p4est_topidx_t tree_id,
                                            p4est_connect_type_t
                                            boundary_type,
                                            int boundary_index,
                                            sc_array_t *
                                            neighbor_transform_array)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 index_lim;
#endif
  int                 dim;

  P4EST_ASSERT (0 <= tree_id && tree_id < conn->num_trees);
  P4EST_ASSERT (P4EST_CONNECT_SELF <= boundary_type
                && boundary_type <= P4EST_CONNECT_FULL);
  P4EST_ASSERT (neighbor_transform_array->elem_size ==
                sizeof (p4est_neighbor_transform_t));
  P4EST_ASSERT (boundary_index >= 0);
  switch (boundary_type) {
  case P4EST_CONNECT_SELF:
#ifdef P4EST_ENABLE_DEBUG
    index_lim = 1;
#endif
    dim = P4EST_DIM;
    break;
  case P4EST_CONNECT_FACE:
#ifdef P4EST_ENABLE_DEBUG
    index_lim = P4EST_FACES;
#endif
    dim = P4EST_DIM - 1;
    break;
  case P4EST_CONNECT_CORNER:
#ifdef P4EST_ENABLE_DEBUG
    index_lim = P4EST_CHILDREN;
#endif
    dim = 0;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
#ifdef P4EST_ENABLE_DEBUG
    index_lim = P8EST_EDGES;
#endif
    dim = 1;
    break;
#endif
  default:
    /* This can only happen for a invalid boundary type. */
    SC_ABORT_NOT_REACHED ();
  }
  P4EST_ASSERT (boundary_index < index_lim);

  /* always add self transformation */
  {
    p4est_neighbor_transform_t *nt = p4est_neighbor_transform_array_push
      (neighbor_transform_array);

    nt->neighbor_type = P4EST_CONNECT_SELF;
    nt->neighbor = tree_id;
    nt->index_self = nt->index_neighbor = 0;
    for (int i = 0; i < P4EST_DIM; i++) {
      nt->origin_self[i] = 0;
      nt->origin_neighbor[i] = 0;
      nt->perm[i] = i;
      nt->sign[i] = 1;
    }
  }
  if (boundary_type == P4EST_CONNECT_SELF) {
    return;
  }

  {
    /* list of trees adjacent to the boundary point */
    int                 nfaces = (dim == P4EST_DIM - 1) ? 1 :
#ifdef P4_TO_P8
      (dim == 1) ? 2 :
#endif
      P4EST_DIM;
    const int          *faces = (dim == P4EST_DIM - 1) ? &boundary_index :
#ifdef P4_TO_P8
      (dim == 1) ? &p8est_edge_faces[boundary_index][0] :
#endif
      &p4est_corner_faces[boundary_index][0];
    const int8_t       *to_face = &conn->tree_to_face[P4EST_FACES * tree_id];

    for (int fi = 0; fi < nfaces; fi++) {
      int                 f = faces[fi];
      int                 ftransform[9];
      int                 ntree =
        p4est_find_face_transform (conn, tree_id, f, ftransform);

      if (ntree >= 0) {
        p4est_neighbor_transform_t *nt = p4est_neighbor_transform_array_push
          (neighbor_transform_array);

        nt->neighbor_type = P4EST_CONNECT_FACE;
        nt->neighbor = ntree;
        nt->index_self = f;
        nt->index_neighbor = to_face[f] % P4EST_FACES;
        p4est_face_transform_to_neighbor_transform (ftransform, nt);
      }
    }
  }
  if (boundary_type == P4EST_CONNECT_FACE) {
    return;
  }

#ifdef P4_TO_P8
  {
    int                 nedges = (dim == 1) ? 1 : 3;
    const int          *edges =
      (dim == 1) ? &boundary_index : &p8est_corner_edges[boundary_index][0];

    for (int ei = 0; ei < nedges; ei++) {
      int                 e = edges[ei];
      p8est_edge_info_t   e_info;
      sc_array_t         *eta = &e_info.edge_transforms;

      sc_array_init (eta, sizeof (p8est_edge_transform_t));
      p8est_find_edge_transform (conn, tree_id, e, &e_info);
      for (size_t iz = 0; iz < eta->elem_count; iz++) {
        p8est_edge_transform_t *et =
          (p8est_edge_transform_t *) sc_array_index (eta, iz);
        p4est_neighbor_transform_t *nt =
          p4est_neighbor_transform_array_push (neighbor_transform_array);

        nt->neighbor_type = P8EST_CONNECT_EDGE;
        nt->index_self = e;
        nt->index_neighbor = et->nedge;
        nt->neighbor = et->ntree;
        p8est_edge_transform_to_neighbor_transform (et, e, nt);
      }
      sc_array_reset (eta);
    }

  }

  if (boundary_type == P8EST_CONNECT_EDGE) {
    return;
  }
#endif

  {
    p4est_corner_info_t c_info;
    sc_array_t         *cta = &c_info.corner_transforms;

    sc_array_init (cta, sizeof (p4est_corner_transform_t));
    p4est_find_corner_transform (conn, tree_id, boundary_index, &c_info);
    for (size_t iz = 0; iz < cta->elem_count; iz++) {
      p4est_neighbor_transform_t *nt =
        p4est_neighbor_transform_array_push (neighbor_transform_array);
      p4est_corner_transform_t *ct =
        (p4est_corner_transform_t *) sc_array_index (cta, iz);

      nt->neighbor = ct->ntree;
      nt->neighbor_type = P4EST_CONNECT_CORNER;
      nt->index_self = boundary_index;
      nt->index_neighbor = ct->ncorner;
      p4est_corner_transform_to_neighbor_transform (ct, boundary_index, nt);
    }

    sc_array_reset (cta);
  }

}

void
p4est_neighbor_transform_coordinates (const p4est_neighbor_transform_t * nt,
                                      const p4est_qcoord_t
                                      self_coords[P4EST_DIM],
                                      p4est_qcoord_t neigh_coords[P4EST_DIM])
{
  p4est_qcoord_t      self_from_origin[P4EST_DIM];

  for (int d = 0; d < P4EST_DIM; d++) {
    self_from_origin[d] = self_coords[d] - nt->origin_self[d];
  }
  for (int d = 0; d < P4EST_DIM; d++) {
    neigh_coords[d] =
      nt->sign[d] * self_from_origin[nt->perm[d]] + nt->origin_neighbor[d];
  }
}

void
p4est_neighbor_transform_coordinates_reverse (const p4est_neighbor_transform_t
                                              * nt,
                                              const p4est_qcoord_t
                                              neigh_coords[P4EST_DIM],
                                              p4est_qcoord_t
                                              self_coords[P4EST_DIM])
{
  p4est_qcoord_t      neigh_from_origin[P4EST_DIM];

  for (int d = 0; d < P4EST_DIM; d++) {
    neigh_from_origin[d] = neigh_coords[d] - nt->origin_neighbor[d];
  }
  for (int d = 0; d < P4EST_DIM; d++) {
    self_coords[nt->perm[d]] =
      nt->sign[d] * neigh_from_origin[d] + nt->origin_self[nt->perm[d]];
  }
}
