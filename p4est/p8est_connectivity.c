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

#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>
#include <p8est.h>

/* *INDENT-OFF* */
const int           p8est_face_corners[6][4] =
{{ 0, 2, 4, 6 },
 { 1, 3, 5, 7 },
 { 0, 1, 4, 5 },
 { 2, 3, 6, 7 },
 { 0, 1, 2, 3 },
 { 4, 5, 6, 7 }};
const int           p8est_face_edges[6][4] =
{{ 4, 6,  8, 10 },
 { 5, 7,  9, 11 },
 { 0, 2,  8,  9 },
 { 1, 3, 10, 11 },
 { 0, 1,  4,  5 },
 { 2, 3,  6,  7 }};
const int           p8est_face_dual[6] = { 1, 0, 3, 2, 5, 4 };
const int           p8est_face_permutations[8][4] =
{{ 0, 1, 2, 3 },                /* no.  0 of 0..23 */
 { 0, 2, 1, 3 },                /* no.  2 of 0..23 */
 { 1, 0, 3, 2 },                /* no.  7 of 0..23 */
 { 1, 3, 0, 2 },                /* no. 10 of 0..23 */
 { 2, 0, 3, 1 },                /* no. 13 of 0..23 */
 { 2, 3, 0, 1 },                /* no. 16 of 0..23 */
 { 3, 1, 2, 0 },                /* no. 21 of 0..23 */
 { 3, 2, 1, 0 }};               /* no. 23 of 0..23 */
const int           p8est_face_permutation_sets[3][4] =
{{ 1, 2, 5, 6 },
 { 0, 3, 4, 7 },
 { 0, 4, 3, 7 }};
const int           p8est_face_permutation_refs[6][6] =
{{ 0, 1, 1, 0, 0, 1 },
 { 2, 0, 0, 1, 1, 0 },
 { 2, 0, 0, 1, 1, 0 },
 { 0, 2, 2, 0, 0, 1 },
 { 0, 2, 2, 0, 0, 1 },
 { 2, 0, 0, 2, 2, 0 }};

const int           p8est_face_edge_permutations[8][4] =
{{ 0, 1, 2, 3 },
 { 0, 1, 3, 2 },
 { 1, 0, 2, 3 },
 { 1, 0, 3, 2 },
 { 2, 3, 0, 1 },
 { 2, 3, 1, 0 },
 { 3, 2, 0, 1 },
 { 3, 2, 1, 0 }};
const int           p8est_face_edge_permutation_sets[3][4] =
{{ 4, 1, 2, 7 },
 { 0, 6, 5, 3 },
 { 0, 5, 6, 3 }};

const int           p8est_edge_faces[12][2] =
{{ 2, 4 },
 { 3, 4 },
 { 2, 5 },
 { 3, 5 },
 { 0, 4 },
 { 1, 4 },
 { 0, 5 },
 { 1, 5 },
 { 0, 2 },
 { 1, 2 },
 { 0, 3 },
 { 1, 3 }};
const int           p8est_edge_corners[12][2] =
{{ 0, 1 },
 { 2, 3 },
 { 4, 5 },
 { 6, 7 },
 { 0, 2 },
 { 1, 3 },
 { 4, 6 },
 { 5, 7 },
 { 0, 4 },
 { 1, 5 },
 { 2, 6 },
 { 3, 7 }};
const int           p8est_edge_edge_corners[12][8] =
{{  0,  1, -1, -1, -1, -1, -1, -1},
 { -1, -1,  0,  1, -1, -1, -1, -1},
 { -1, -1, -1, -1,  0,  1, -1, -1},
 { -1, -1, -1, -1, -1, -1,  0,  1},
 {  0, -1,  1, -1, -1, -1, -1, -1},
 { -1,  0, -1,  1, -1, -1, -1, -1},
 { -1, -1, -1, -1,  0, -1,  1, -1},
 { -1, -1, -1, -1, -1,  0, -1,  1},
 {  0, -1, -1, -1,  1, -1, -1, -1},
 { -1,  0, -1, -1, -1,  1, -1, -1},
 { -1, -1,  0, -1, -1, -1,  1, -1},
 { -1, -1, -1,  0, -1, -1, -1,  1}};
const int           p8est_edge_face_corners[12][6][2] =
{{{ -1, -1 }, { -1, -1 }, {  0,  1 }, { -1, -1 }, {  0,  1 }, { -1, -1 }},
 {{ -1, -1 }, { -1, -1 }, { -1, -1 }, {  0,  1 }, {  2,  3 }, { -1, -1 }},
 {{ -1, -1 }, { -1, -1 }, {  2,  3 }, { -1, -1 }, { -1, -1 }, {  0,  1 }},
 {{ -1, -1 }, { -1, -1 }, { -1, -1 }, {  2,  3 }, { -1, -1 }, {  2,  3 }},
 {{  0,  1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, {  0,  2 }, { -1, -1 }},
 {{ -1, -1 }, {  0,  1 }, { -1, -1 }, { -1, -1 }, {  1,  3 }, { -1, -1 }},
 {{  2,  3 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, {  0,  2 }},
 {{ -1, -1 }, {  2,  3 }, { -1, -1 }, { -1, -1 }, { -1, -1 }, {  1,  3 }},
 {{  0,  2 }, { -1, -1 }, {  0,  2 }, { -1, -1 }, { -1, -1 }, { -1, -1 }},
 {{ -1, -1 }, {  0,  2 }, {  1,  3 }, { -1, -1 }, { -1, -1 }, { -1, -1 }},
 {{  1,  3 }, { -1, -1 }, { -1, -1 }, {  0,  2 }, { -1, -1 }, { -1, -1 }},
 {{ -1, -1 }, {  1,  3 }, { -1, -1 }, {  1,  3 }, { -1, -1 }, { -1, -1 }}};
const int           p8est_edge_face_edges[12][6] =
{{ -1, -1,  0, -1,  0, -1 },
 { -1, -1, -1,  0,  1, -1 },
 { -1, -1,  1, -1, -1,  0 },
 { -1, -1, -1,  1, -1,  1 },
 {  0, -1, -1, -1,  2, -1 },
 { -1,  0, -1, -1,  3, -1 },
 {  1, -1, -1, -1, -1,  2 },
 { -1,  1, -1, -1, -1,  3 },
 {  2, -1,  2, -1, -1, -1 },
 { -1,  2,  3, -1, -1, -1 },
 {  3, -1, -1,  2, -1, -1 },
 { -1,  3, -1,  3, -1, -1 }};

const int           p8est_corner_faces[8][3] =
{{ 0, 2, 4 },
 { 1, 2, 4 },
 { 0, 3, 4 },
 { 1, 3, 4 },
 { 0, 2, 5 },
 { 1, 2, 5 },
 { 0, 3, 5 },
 { 1, 3, 5 }};
const int           p8est_corner_edges[8][3] =
{{ 0, 4,  8 },
 { 0, 5,  9 },
 { 1, 4, 10 },
 { 1, 5, 11 },
 { 2, 6,  8 },
 { 2, 7,  9 },
 { 3, 6, 10 },
 { 3, 7, 11 }};
const int           p8est_corner_face_corners[8][6] =
{{  0, -1,  0, -1,  0, -1 },
 { -1,  0,  1, -1,  1, -1 },
 {  1, -1, -1,  0,  2, -1 },
 { -1,  1, -1,  1,  3, -1 },
 {  2, -1,  2, -1, -1,  0 },
 { -1,  2,  3, -1, -1,  1 },
 {  3, -1, -1,  2, -1,  2 },
 { -1,  3, -1,  3, -1,  3 }};
const int           p8est_corner_edge_corners[8][12] =
{{  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1 },
 {  1, -1, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1 },
 { -1,  0, -1, -1,  1, -1, -1, -1, -1, -1,  0, -1 },
 { -1,  1, -1, -1, -1,  1, -1, -1, -1, -1, -1,  0 },
 { -1, -1,  0, -1, -1, -1,  0, -1,  1, -1, -1, -1 },
 { -1, -1,  1, -1, -1, -1, -1,  0, -1,  1, -1, -1 },
 { -1, -1, -1,  0, -1, -1,  1, -1, -1, -1,  1, -1 },
 { -1, -1, -1,  1, -1, -1, -1,  1, -1, -1, -1,  1 }};
const int           p8est_child_edge_faces[8][12] =
{{ -1,  4,  2, -1, -1,  4,  0, -1, -1,  2,  0, -1 },
 { -1,  4,  2, -1,  4, -1, -1,  1,  2, -1, -1,  1 },
 {  4, -1, -1,  3, -1,  4,  0, -1,  0, -1, -1,  3 },
 {  4, -1, -1,  3,  4, -1, -1,  1, -1,  1,  3, -1 },
 {  2, -1, -1,  5,  0, -1, -1,  5, -1,  2,  0, -1 },
 {  2, -1, -1,  5, -1,  1,  5, -1,  2, -1, -1,  1 },
 { -1,  3,  5, -1,  0, -1, -1,  5,  0, -1, -1,  3 },
 { -1,  3,  5, -1, -1,  1,  5, -1, -1,  1,  3, -1 }};
const int           p8est_child_corner_faces[8][8] =
{{ -1, -1, -1,  4, -1,  2,  0, -1 },
 { -1, -1,  4, -1,  2, -1, -1,  1 },
 { -1,  4, -1, -1,  0, -1, -1,  3 },
 {  4, -1, -1, -1, -1,  1,  3, -1 },
 { -1,  2,  0, -1, -1, -1, -1,  5 },
 {  2, -1, -1,  1, -1, -1,  5, -1 },
 {  0, -1, -1,  3, -1,  5, -1, -1 },
 { -1,  1,  3, -1,  5, -1, -1, -1 }};
const int           p8est_child_corner_edges[8][8] =
{{ -1,  0,  4, -1,  8, -1, -1, -1 },
 {  0, -1, -1,  5, -1,  9, -1, -1 },
 {  4, -1, -1,  1, -1, -1, 10, -1 },
 { -1,  5,  1, -1, -1, -1, -1, 11 },
 {  8, -1, -1, -1, -1,  2,  6, -1 },
 { -1,  9, -1, -1,  2, -1, -1,  7 },
 { -1, -1, 10, -1,  6, -1, -1,  3 },
 { -1, -1, -1, 11, -1,  7,  3, -1 }};
/* *INDENT-ON* */

int
p8est_connectivity_face_neighbor_corner_set (int c, int f, int nf, int set)
{
  int                 fc, nfc;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= set && set < 8);

  fc = p8est_corner_face_corners[c][f];
  P4EST_ASSERT (0 <= fc && fc < P4EST_HALF);
  nfc = p8est_face_permutations[set][fc];
  P4EST_ASSERT (0 <= nfc && nfc < P4EST_HALF);

  return p8est_face_corners[nf][nfc];
}

p4est_connectivity_t *
p8est_connectivity_new_unitcube (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_ett = 0;
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
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    0, 1, 2, 3, 4, 5,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p8est_connectivity_new_rotwrap (void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 1;
  const p4est_topidx_t num_edges = 4;
  const p4est_topidx_t num_corners = 1;
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
  const p4est_topidx_t tree_to_vertex[1 * 8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };
  const p4est_topidx_t tree_to_tree[1 * 6] = {
    0, 0, 0, 0, 0, 0,
  };
  const int8_t        tree_to_face[1 * 6] = {
    1, 0, 2, 3, 11, 10,
  };
  const p4est_topidx_t tree_to_edge[1 * 12] = {
    0, 0, 1, 1, 1, 1, 0, 0, 2, 2, 3, 3,
  };
  const p4est_topidx_t ett_offset[4 + 1] = {
    0, 4, 8, 10, 12,
  };
  const p4est_topidx_t edge_to_tree[12] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int8_t        edge_to_edge[12] = {
    0, 7, 1, 6, 2, 16, 3, 17, 8, 9, 10, 11,
  };
  const p4est_topidx_t tree_to_corner[1 * 8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 8,
  };
  const p4est_topidx_t corner_to_tree[8] = {
    0, 0, 0, 0, 0, 0, 0, 0,
  };
  const int8_t        corner_to_corner[8] = {
    0, 1, 2, 3, 4, 5, 6, 7,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p8est_connectivity_new_drop (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 24;
  const p4est_topidx_t num_trees = 5;
  const p4est_topidx_t num_ett = 2;
  const p4est_topidx_t num_ctt = 1;
  const double        vertices[24 * 3] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    3, 1, 0,
    1, 2, 0,
    3, 2, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1,
    1, 2, 1,
    2, 2, 1,
    1, 0, 2,
    1, 1, 2,
    2, 1, 2,
    2, 2, 2,
    0, 0, 3,
    0, 1, 3,
    3, 1, 3,
    3, 2, 3,
    0, 2, 3,
    1, 2, 2,
  };
  const p4est_topidx_t tree_to_vertex[5 * 8] = {
     0,  1,  2,  3,  7,  8,  9, 10,
     3,  4,  5,  6, 10, 11, 12, 13,
    11,  4, 13,  6, 16, 20, 17, 21,
    15, 16, 23, 17, 19, 20, 22, 21,
     7,  8,  9, 10, 18, 14, 19, 15,
  };
  const p4est_topidx_t tree_to_tree[5 * 6] = {
    0, 0, 0, 0, 0, 4,
    1, 2, 1, 1, 1, 1,
    2, 2, 2, 2, 1, 3,
    3, 2, 3, 3, 3, 3,
    4, 4, 4, 4, 0, 4,
  };
  const int8_t        tree_to_face[5 * 6] = {
    0,  1, 2, 3, 4, 4,
    0, 10, 2, 3, 4, 5,
    0,  1, 2, 3, 7, 1,
    0,  5, 2, 3, 4, 5,
    0,  1, 2, 3, 5, 5,
  };
  const p4est_topidx_t tree_to_edge[5 * 12] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,
    -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,  1, -1, -1, -1,
    -1, -1, -1,  1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ett_offset[2 + 1] = { 0, 2, 4 };
  const p4est_topidx_t edge_to_tree[4] = {
    0, 1, 3, 4
  };
  const int8_t edge_to_edge[4] = {
    11, 8, 8, 15
  };
  const p4est_topidx_t tree_to_corner[5 * 8] = {
    -1, -1, -1, -1, -1, -1, -1,  0,
    -1, -1, -1, -1,  0, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1,  0, -1, -1, -1, -1,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = {
    0, 3
  };
  const p4est_topidx_t corner_to_tree[3] = {
    0, 1, 4
  };
  const int8_t corner_to_corner[3] = {
    7, 4, 3
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      num_ett, num_ctt,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p8est_connectivity_new_twocubes (void)
{
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ett = 0;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[12 * 3] = {
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    0, 0, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[2 * 8] = {
    0, 1, 3, 4, 6, 7, 9, 10,
    1, 2, 4, 5, 7, 8, 10, 11,
  };
  const p4est_topidx_t tree_to_tree[2 * 6] = {
    0, 1, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1,
  };
  const int8_t        tree_to_face[2 * 6] = {
    0, 0, 2, 3, 4, 5,
    1, 1, 2, 3, 4, 5,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p8est_connectivity_new_twowrap (void)
{
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t num_trees = 2;
  const p4est_topidx_t num_ett = 0;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[12 * 3] = {
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    0, 0, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[2 * 8] = {
    3, 9, 0, 6, 4, 10, 1, 7,
    8, 2, 7, 1, 11, 5, 10, 4,
  };
  const p4est_topidx_t tree_to_tree[2 * 6] = {
    0, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 1, 1,
  };
  const int8_t        tree_to_face[2 * 6] = {
    0, 1, 2, 3, 20, 21,
    0, 1, 22, 23, 4, 5,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ett, NULL, NULL,
                                      NULL, &num_ctt, NULL, NULL);
}

p4est_connectivity_t *
p8est_connectivity_new_rotcubes (void)
{
  const p4est_topidx_t num_vertices = 26;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_edges = 3;
  const p4est_topidx_t num_corners = 1;
  const double        vertices[26 * 3] = {
    0, 0, 0,
    1, 0, 2,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    1, -1, 0,
    2, -1, 0,
    1, -1, 1,
    2, -1, 1,
    2, 1, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    0, 0, 1,
    0, 0, 2,
    1, 0, 0,
    1, 1, 2,
    0, 1, 2,
    2.5, 1.5, 2,
    2, 1.5, 2,
    2, 1.5, 2.5,
    2, .5, 2.5,
    2.5, .5, 2,
    2, .5, 2,
  };
  const p4est_topidx_t tree_to_vertex[6 * 8] = {
    0, 17, 3, 4, 15, 11, 13, 14,
    7, 2, 6, 17, 9, 12, 8, 11,
    2, 12, 5, 10, 17, 11, 4, 14,
    19, 13, 18, 14, 16, 15, 1, 11,
    14, 11, 21, 25, 18, 1, 22, 23,
    21, 20, 25, 24, 14, 10, 11, 12,
  };
  const p4est_topidx_t tree_to_tree[6 * 6] = {
    0, 2, 0, 0, 0, 3,
    1, 2, 1, 1, 1, 1,
    2, 5, 1, 2, 2, 0,
    3, 0, 3, 4, 3, 3,
    4, 4, 3, 4, 5, 4,
    4, 5, 5, 5, 5, 2,
  };
  const int8_t        tree_to_face[6 * 6] = {
    0, 5, 2, 3, 4, 13,
    0, 2, 2, 3, 4, 5,
    0, 23, 1, 3, 4, 1,
    0, 17, 2, 8, 4, 5,
    0, 1, 9, 3, 12, 5,
    16, 1, 2, 3, 4, 19,
  };
  const p4est_topidx_t tree_to_edge[6 * 12] = {
    -1, -1, -1, -1, -1, -1, -1, 0, -1, 2, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 2,
    -1, -1, 2, -1, -1, -1, -1, 0, -1, 1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,
    0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, 1, -1, -1, 0, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ett_offset[3 + 1] = { 0, 5, 8, 11 };
  const p4est_topidx_t edge_to_tree[11] = {
    0, 2, 3, 4, 5, 1, 2, 5, 0, 1, 2
  };
  const int8_t        edge_to_edge[11] = {
    7, 7, 23, 12, 18, 7, 9, 15, 9, 11, 2
  };
  const p4est_topidx_t tree_to_corner[6 * 8] = {
    -1, -1, -1, -1, -1, 0, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 0,
    -1, -1, -1, -1, -1, 0, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 0,
    -1, 0, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, 0, -1,
  };
  const p4est_topidx_t ctt_offset[1 + 1] = { 0, 6 };
  const p4est_topidx_t corner_to_tree[6] = { 0, 1, 2, 3, 4, 5 };
  const int8_t        corner_to_corner[6] = { 5, 7, 5, 7, 1, 6 };

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);
}

p4est_connectivity_t *
p8est_connectivity_new_shell (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 18;
  const p4est_topidx_t num_trees =    24;
  const p4est_topidx_t num_edges =    18;
  const p4est_topidx_t num_corners =   0;
  const p4est_topidx_t ctt_offset =    0;
  const double        vertices[18 * 3] = {
    -1, -1,  1,
     0, -1,  1,
     1, -1,  1,
    -1,  0,  1,
     0,  0,  1,
     1,  0,  1,
    -1,  1,  1,
     0,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     0, -1,  2,
     1, -1,  2,
    -1,  0,  2,
     0,  0,  2,
     1,  0,  2,
    -1,  1,  2,
     0,  1,  2,
     1,  1,  2,
  };
  const p4est_topidx_t tree_to_vertex[24 * 8] = {
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
    0, 1, 3, 4,  9, 10, 12, 13,
    1, 2, 4, 5, 10, 11, 13, 14,
    3, 4, 6, 7, 12, 13, 15, 16,
    4, 5, 7, 8, 13, 14, 16, 17,
  };
  const p4est_topidx_t tree_to_tree[24 * 6] = {
    18,  1, 14,  2,  0,  0,
     0, 23, 15,  3,  1,  1,
    16,  3,  0,  4,  2,  2,
     2, 21,  1,  5,  3,  3,
    16,  5,  2,  6,  4,  4,
     4, 21,  3,  7,  5,  5,
    17,  7,  4,  8,  6,  6,
     6, 20,  5,  9,  7,  7,
    17,  9,  6, 10,  8,  8,
     8, 20,  7, 11,  9,  9,
    19, 11,  8, 12, 10, 10,
    10, 22,  9, 13, 11, 11,
    19, 13, 10, 14, 12, 12,
    12, 22, 11, 15, 13, 13,
    18, 15, 12,  0, 14, 14,
    14, 23, 13,  1, 15, 15,
     2, 17,  4, 18, 16, 16,
    16,  8,  6, 19, 17, 17,
     0, 19, 16, 14, 18, 18,
    18, 10, 17, 12, 19, 19,
     9, 21,  7, 22, 20, 20,
    20,  3,  5, 23, 21, 21,
    11, 23, 20, 13, 22, 22,
    22,  1, 21, 15, 23, 23,
  };
  const int8_t        tree_to_face[24 * 6] = {
    6, 0, 3, 2, 4, 5,
    1, 7, 3, 2, 4, 5,
    6, 0, 3, 2, 4, 5,
    1, 7, 3, 2, 4, 5,
    2, 0, 3, 2, 4, 5,
    1, 8, 3, 2, 4, 5,
    2, 0, 3, 2, 4, 5,
    1, 8, 3, 2, 4, 5,
    1, 0, 3, 2, 4, 5,
    1, 0, 3, 2, 4, 5,
    1, 0, 3, 2, 4, 5,
    1, 0, 3, 2, 4, 5,
    9, 0, 3, 2, 4, 5,
    1, 3, 3, 2, 4, 5,
    9, 0, 3, 2, 4, 5,
    1, 3, 3, 2, 4, 5,
    6, 0, 0, 2, 4, 5,
    1, 0, 0, 2, 4, 5,
    6, 0, 3, 6, 4, 5,
    1, 0, 3, 6, 4, 5,
    1, 0, 7, 2, 4, 5,
    1, 7, 7, 2, 4, 5,
    1, 0, 3, 1, 4, 5,
    1, 7, 3, 1, 4, 5,
  };
  const p4est_topidx_t tree_to_edge[24 * 12] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  8,  6,  0,
    -1, -1, -1, -1, -1, -1, -1, -1,  8, -1,  0,  7,
    -1, -1, -1, -1, -1, -1, -1, -1,  6,  0, -1,  9,
    -1, -1, -1, -1, -1, -1, -1, -1,  0,  7,  9, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1,  9, 10,  1,
    -1, -1, -1, -1, -1, -1, -1, -1,  9, -1,  1, 11,
    -1, -1, -1, -1, -1, -1, -1, -1, 10,  1, -1, 12,
    -1, -1, -1, -1, -1, -1, -1, -1,  1, 11, 12, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 12, 13,  2,
    -1, -1, -1, -1, -1, -1, -1, -1, 12, -1,  2, 14,
    -1, -1, -1, -1, -1, -1, -1, -1, 13,  2, -1, 15,
    -1, -1, -1, -1, -1, -1, -1, -1,  2, 14, 15, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 16,  3,
    -1, -1, -1, -1, -1, -1, -1, -1, 15, -1,  3, 17,
    -1, -1, -1, -1, -1, -1, -1, -1, 16,  3, -1,  8,
    -1, -1, -1, -1, -1, -1, -1, -1,  3, 17,  8, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 10,  6,  4,
    -1, -1, -1, -1, -1, -1, -1, -1, 10, -1,  4, 13,
    -1, -1, -1, -1, -1, -1, -1, -1,  6,  4, -1, 16,
    -1, -1, -1, -1, -1, -1, -1, -1,  4, 13, 16, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, 11, 14,  5,
    -1, -1, -1, -1, -1, -1, -1, -1, 11, -1,  5 , 7,
    -1, -1, -1, -1, -1, -1, -1, -1, 14,  5, -1, 17,
    -1, -1, -1, -1, -1, -1, -1, -1,  5,  7, 17, -1,
  };
  const p4est_topidx_t ett_offset[18 + 1] = {
     0,  4 , 8, 12, 16, 20, 24, 28, 32,
    36, 40, 44, 48, 52, 56, 60, 64, 68, 72,
  };
  const p4est_topidx_t edge_to_tree[72] = {
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12, 13, 14, 15,
    16, 17, 18, 19,
    20, 21, 22, 23,
     0,  2, 16, 18,
     1,  3, 21, 23,
     0,  1, 14, 15,
     2,  3,  4,  5,
     4,  6, 16, 17,
     5,  7, 20, 21,
     6,  7,  8,  9,
     8, 10, 17, 19,
     9, 11, 20, 22,
    10, 11, 12, 13,
    12, 14, 18, 19,
    13, 15, 22, 23,
  };
  const int8_t        edge_to_edge[72] = {
    11, 10,  9,  8,
    11, 10,  9,  8,
    11, 10,  9,  8,
    11, 10,  9,  8,
    11, 10,  9,  8,
    11, 10,  9,  8,
    10,  8, 10,  8,
    11,  9, 11,  9,
     9,  8, 11, 10,
    11, 10,  9,  8,
    10,  8,  9,  8,
    11,  9,  9,  8,
    11, 10,  9,  8,
    10,  8, 11,  9,
    11,  9, 10,  8,
    11, 10,  9,  8,
    10,  8, 11, 10,
    11,  9, 11, 10,
  };
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      NULL, &ctt_offset, NULL, NULL);
}

p4est_connectivity_t *
p8est_connectivity_new_sphere (void)
{
/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 16;
  const p4est_topidx_t num_trees =    13;
  const p4est_topidx_t num_edges =    12;
  const p4est_topidx_t num_corners =   0;
  const p4est_topidx_t ctt_offset = 0;
  const double        vertices[16 * 3] = {
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
    -1, -1,  2,
     1, -1,  2,
    -1,  1,  2,
     1,  1,  2,
    -1, -1, -1,
     1, -1, -1,
    -1,  1, -1,
     1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
    -1,  1,  1,
     1,  1,  1,
  };
  const p4est_topidx_t tree_to_vertex[13 * 8] = {
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    0,  1,  2,  3,  4,  5,  6,  7,
    8,  9, 10, 11, 12, 13, 14, 15,
  };
  const p4est_topidx_t tree_to_tree[13 * 6] = {
     5,  3,  4,  1,  6,  0,
     5,  3,  0,  2,  7,  1,
     5,  3,  1,  4,  8,  2,
     2,  0,  1,  4,  9,  3,
     2,  0,  3,  5, 10,  4,
     2,  0,  4,  1, 11,  5,
    11,  9, 10,  7, 12,  0,
    11,  9,  6,  8, 12,  1,
    11,  9,  7, 10, 12,  2,
     8,  6,  7, 10, 12,  3,
     8,  6,  9, 11, 12,  4,
     8,  6, 10,  7, 12,  5,
    11,  9,  6,  8, 10,  7,
  };
  const int8_t        tree_to_face[13 * 6] = {
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  5,  5,
     9,  8,  3,  2,  5,  5,
     6,  0,  3,  6,  5,  5,
     1,  7,  7,  2,  2,  4,
     9,  8,  3,  2,  5,  4,
     6,  0,  3,  6, 15,  4,
     1,  7,  7,  2, 19,  4,
     9,  8,  3,  2, 22,  4,
     6,  0,  3,  6,  6,  4,
    10, 22,  4, 16, 22,  4,
  };
  const p4est_topidx_t tree_to_edge[13 * 12] = {
     0,  2, -1, -1,  8,  9, -1, -1, -1, -1, -1, -1,
     2,  3, -1, -1,  6,  7, -1, -1, -1, -1, -1, -1,
     3,  1, -1, -1, 10, 11, -1, -1, -1, -1, -1, -1,
     7,  5, -1, -1, 11,  9, -1, -1, -1, -1, -1, -1,
     5,  4, -1, -1,  1,  0, -1, -1, -1, -1, -1, -1,
     4,  6, -1, -1, 10,  8, -1, -1, -1, -1, -1, -1,
    -1, -1,  0,  2, -1, -1,  8,  9, -1, -1, -1, -1,
    -1, -1,  2,  3, -1, -1,  6,  7, -1, -1, -1, -1,
    -1, -1,  3,  1, -1, -1, 10, 11, -1, -1, -1, -1,
    -1, -1,  7,  5, -1, -1, 11,  9, -1, -1, -1, -1,
    -1, -1,  5,  4, -1, -1,  1,  0, -1, -1, -1, -1,
    -1, -1,  4,  6, -1, -1, 10,  8, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ett_offset[12 + 1] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,
  };
  const p4est_topidx_t edge_to_tree[48] = {
    0,  4,  6, 10,
    2,  4,  8, 10,
    0,  1,  6,  7,
    1,  2,  7,  8,
    4,  5, 10, 11,
    3,  4,  9, 10,
    1,  5,  7, 11,
    1,  3,  7,  9,
    0,  5,  6, 11,
    0,  3,  6,  9,
    2,  5,  8, 11,
    2,  3,  8,  9,
  };
  const int8_t        edge_to_edge[48] = {
     0, 17,  2, 19,
     1, 16,  3, 18,
     1,  0,  3,  2,
     1,  0,  3,  2,
    13, 12, 15, 14,
    13, 12, 15, 14,
     4, 13,  6, 15,
     5, 12,  7, 14,
     4,  5,  6,  7,
     5, 17,  7, 19,
    16,  4, 18,  6,
    17, 16, 19, 18,
  };
#if 0   /* corner information would be redundant */
  const p4est_topidx_t tree_to_corner[13 * 8] = {
     0,  1,  4,  5, -1, -1, -1, -1,
     4,  5,  6,  7, -1, -1, -1, -1,
     6,  7,  2,  3, -1, -1, -1, -1,
     7,  5,  3,  1, -1, -1, -1, -1,
     3,  1,  2,  0, -1, -1, -1, -1,
     2,  0,  6,  4, -1, -1, -1, -1,
    -1, -1, -1, -1,  0,  1,  4,  5,
    -1, -1, -1, -1,  4,  5,  6,  7,
    -1, -1, -1, -1,  6,  7,  2,  3,
    -1, -1, -1, -1,  7,  5,  3,  1,
    -1, -1, -1, -1,  3,  1,  2,  0,
    -1, -1, -1, -1,  2,  0,  6,  4,
    -1, -1, -1, -1, -1, -1, -1, -1,
  };
  const p4est_topidx_t ctt_offset[8 + 1] = {
    0, 6, 12, 18, 24, 30, 36, 42, 48,
  };
  const p4est_topidx_t corner_to_tree[48] = {
    0,  4,  5,  6, 10, 11,
    0,  3,  4,  6,  9, 10,
    2,  4,  5,  8, 10, 11,
    2,  3,  4,  8,  9, 10,
    0,  1,  5,  6,  7, 11,
    0,  1,  3,  6,  7,  9,
    1,  2,  5,  7,  8, 11,
    1,  2,  3,  7,  8,  9,
  };
  const int8_t        corner_to_corner[48] = {
    0, 3, 1, 4, 7, 5,
    1, 3, 1, 5, 7, 5,
    2, 2, 0, 6, 6, 4,
    3, 2, 0, 7, 6, 4,
    2, 0, 3, 6, 4, 7,
    3, 1, 1, 7, 5, 5,
    2, 0, 2, 6, 4, 6,
    3, 1, 0, 7, 5, 4,
  };
#endif  /* 0 */
/* *INDENT-ON* */

  return p4est_connectivity_new_copy (num_vertices, num_trees,
                                      num_edges, num_corners,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      NULL, &ctt_offset, NULL, NULL);
}

/*
 * Torus section, displaying here treeId's in the first segment.
 *       ______
 *      /      \
 *     /\   3  /\
 *    /  \ __ /  \
 *   |  2 | 4|  0 |
 *   |    |__|    |
 *    \  /    \  /
 *     \/  1   \/
 *      \______/
 *
 * There are 8 edges per segment on the -z face (see drawing below):
 *
 * - 4 edges around -z faces of tree 4,
 *   * tree 0 <-> tree 4
 *   * tree 1 <-> tree 4
 *   * tree 2 <-> tree 4
 *   * tree 3 <-> tree 4
 * - 4 edges:
 *   * tree 0 <-> tree 1
 *   * tree 1 <-> tree 2
 *   * tree 2 <-> tree 3
 *   * tree 3 <-> tree 0
 *
 *  6     7
 *   \_3_/
 *   |   |
 *   2   0
 *   |_1_|
 *  /     \
 * 5       4
 *
 * There are no corners (linking 8 trees).
 */

p8est_connectivity_t *
p8est_connectivity_new_torus (int nSegments)
{

/* *INDENT-OFF* */
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t nTreesPS     =  5; /* number of trees per segment */
  const p4est_topidx_t num_trees    =  nTreesPS*nSegments;
  const p4est_topidx_t nEdgesPS     =  8; /* number of edges per segment */
  const p4est_topidx_t num_edges    =  nEdgesPS*nSegments;
  const p4est_topidx_t num_corners  =  0;
  const p4est_topidx_t num_ctt      =  0; /* corner to tree */
  const p4est_topidx_t num_ett      =  nEdgesPS*nSegments*4; /* edge to tree */
  int i, j, iSegment, nbItems, iTree, iEdge;
  p4est_connectivity_t *conn;

  conn = p4est_connectivity_new (num_vertices, num_trees,
                                 num_edges, num_ett,
                                 num_corners, num_ctt);

  {
    const double vertices[12 * 3] = {
      -1,  -1,  0,
       1,  -1,  0,
      -1,   1,  0,
       1,   1,  0,
      -1,   2,  0,
       1,   2,  0,
      -1,  -1,  1,
       1,  -1,  1,
      -1,   1,  1,
       1,   1,  1,
      -1,   2,  1,
       1,   2,  1,
    };

    for (i=0; i<num_vertices*3; ++i)
      conn->vertices[i] = vertices[i];
  }

  /* tree to vertex */
  {
    const p4est_topidx_t tree_to_vertex[5 * 8] = {
      2, 3, 4, 5, 8, 9, 10, 11,  /* tree 0 */
      2, 3, 4, 5, 8, 9, 10, 11,  /* tree 1 */
      2, 3, 4, 5, 8, 9, 10, 11,  /* tree 2 */
      2, 3, 4, 5, 8, 9, 10, 11,  /* tree 3 */
      0, 1, 2, 3, 6, 7,  8,  9,  /* tree 4  - center*/
    };

    nbItems = nTreesPS*8; /* per segment */

    /* all segments use the same pattern */
    for (iSegment=0; iSegment<nSegments; ++iSegment)
    {
      for (j=0; j<nbItems; ++j)
      {
        conn->tree_to_vertex[j+iSegment*nbItems] = tree_to_vertex[j];
      }
    }
  }
/* *INDENT-ON* */

  /*  tree to tree */
  {
    /*  return global tree id from local tree id and segment id */
#define tGlob(tLoc,iSeg) ( (tLoc) + (iSeg) * nTreesPS )

    /* return global tree id above (+z) */
#define tGlobZp(tLoc,iSeg) ( ( tGlob(tLoc,iSeg) >= num_trees) ? tGlob(tLoc,iSeg) - num_trees : tGlob(tLoc,iSeg) )

    /* return global tree id below (-z) */
#define tGlobZm(tLoc,iSeg) ( ( tGlob(tLoc,iSeg) < 0) ? tGlob(tLoc,iSeg) + num_trees : tGlob(tLoc,iSeg) )

    /* const p4est_topidx_t tree_to_tree[5 * 6] = { */
    /*   3, 1, 4, self, 0, 0,  /\* tree 0 *\/ */
    /*   0, 2, 4, self, 1, 1,  /\* tree 1 *\/ */
    /*   1, 3, 4, self, 2, 2,  /\* tree 2 *\/ */
    /*   2, 0, 4, self, 3, 3,  /\* tree 3 *\/ */
    /*   2, 0, 1, 3, 4, 4,  /\* tree 4 - center *\/ */
    /* }; */

    nbItems = nTreesPS * 6;     /* 5 trees per segment x 6 faces */

    /*  Global tree id */
    for (iSegment = 0; iSegment < nSegments; ++iSegment) {

      iTree = 0;

      conn->tree_to_tree[0 + iTree * 6 + iSegment * nbItems] = tGlob (3, iSegment);     /* tree = 0 modulo 5  */
      conn->tree_to_tree[1 + iTree * 6 + iSegment * nbItems] =
        tGlob (1, iSegment);
      conn->tree_to_tree[2 + iTree * 6 + iSegment * nbItems] =
        tGlob (4, iSegment);
      conn->tree_to_tree[3 + iTree * 6 + iSegment * nbItems] =
        tGlob (iTree, iSegment);
      conn->tree_to_tree[4 + iTree * 6 + iSegment * nbItems] =
        tGlobZm (iTree - nTreesPS, iSegment);
      conn->tree_to_tree[5 + iTree * 6 + iSegment * nbItems] =
        tGlobZp (iTree + nTreesPS, iSegment);

      iTree++;

      conn->tree_to_tree[0 + iTree * 6 + iSegment * nbItems] = tGlob (0, iSegment);     /* tree = 1 modulo 5  */
      conn->tree_to_tree[1 + iTree * 6 + iSegment * nbItems] =
        tGlob (2, iSegment);
      conn->tree_to_tree[2 + iTree * 6 + iSegment * nbItems] =
        tGlob (4, iSegment);
      conn->tree_to_tree[3 + iTree * 6 + iSegment * nbItems] =
        tGlob (iTree, iSegment);
      conn->tree_to_tree[4 + iTree * 6 + iSegment * nbItems] =
        tGlobZm (iTree - nTreesPS, iSegment);
      conn->tree_to_tree[5 + iTree * 6 + iSegment * nbItems] =
        tGlobZp (iTree + nTreesPS, iSegment);

      iTree++;

      conn->tree_to_tree[0 + iTree * 6 + iSegment * nbItems] = tGlob (1, iSegment);     /* tree = 2 modulo 5  */
      conn->tree_to_tree[1 + iTree * 6 + iSegment * nbItems] =
        tGlob (3, iSegment);
      conn->tree_to_tree[2 + iTree * 6 + iSegment * nbItems] =
        tGlob (4, iSegment);
      conn->tree_to_tree[3 + iTree * 6 + iSegment * nbItems] =
        tGlob (iTree, iSegment);
      conn->tree_to_tree[4 + iTree * 6 + iSegment * nbItems] =
        tGlobZm (iTree - nTreesPS, iSegment);
      conn->tree_to_tree[5 + iTree * 6 + iSegment * nbItems] =
        tGlobZp (iTree + nTreesPS, iSegment);

      iTree++;

      conn->tree_to_tree[0 + iTree * 6 + iSegment * nbItems] = tGlob (2, iSegment);     /* tree = 3 modulo 5  */
      conn->tree_to_tree[1 + iTree * 6 + iSegment * nbItems] =
        tGlob (0, iSegment);
      conn->tree_to_tree[2 + iTree * 6 + iSegment * nbItems] =
        tGlob (4, iSegment);
      conn->tree_to_tree[3 + iTree * 6 + iSegment * nbItems] = tGlob (iTree, iSegment); /* self */
      conn->tree_to_tree[4 + iTree * 6 + iSegment * nbItems] =
        tGlobZm (iTree - nTreesPS, iSegment);
      conn->tree_to_tree[5 + iTree * 6 + iSegment * nbItems] =
        tGlobZp (iTree + nTreesPS, iSegment);

      iTree++;

      conn->tree_to_tree[0 + iTree * 6 + iSegment * nbItems] = tGlob (2, iSegment);     /* tree = 4 modulo 5  */
      conn->tree_to_tree[1 + iTree * 6 + iSegment * nbItems] =
        tGlob (0, iSegment);
      conn->tree_to_tree[2 + iTree * 6 + iSegment * nbItems] =
        tGlob (1, iSegment);
      conn->tree_to_tree[3 + iTree * 6 + iSegment * nbItems] =
        tGlob (3, iSegment);
      conn->tree_to_tree[4 + iTree * 6 + iSegment * nbItems] =
        tGlobZm (iTree - nTreesPS, iSegment);
      conn->tree_to_tree[5 + iTree * 6 + iSegment * nbItems] =
        tGlobZp (iTree + nTreesPS, iSegment);

      iTree++;

    }

#undef tGlobZm
#undef tGlobZp
#undef tGlob

  }

  /* tree to face */
  {
    /* const int8_t         tree_to_face[5 * 6] = { */
    /*   1, 0, 7, 3, 4, 5,  /\* tree 0 *\/ */
    /*   1, 0, 8, 3, 4, 5,  /\* tree 1 *\/ */
    /*   1, 0, 0, 3, 4, 5,  /\* tree 2 *\/ */
    /*   1, 0, 3, 3, 4, 5,  /\* tree 3 *\/ */
    /*   2, 8, 8, 2, 4, 5,  /\* tree 4 - center *\/ */
    /* }; */

    /* all segments use the same pattern */
    iTree = 0;                  /* global treeId */
    for (iSegment = 0; iSegment < nSegments; ++iSegment) {

      conn->tree_to_face[0 + iTree * 6] = 1;    /* tree = 0 modulo 5  */
      conn->tree_to_face[1 + iTree * 6] = 0;
      conn->tree_to_face[2 + iTree * 6] = 7;
      conn->tree_to_face[3 + iTree * 6] = 3;
      conn->tree_to_face[4 + iTree * 6] = 5;
      conn->tree_to_face[5 + iTree * 6] = 4;
      iTree++;

      conn->tree_to_face[0 + iTree * 6] = 1;    /* tree = 1 modulo 5  */
      conn->tree_to_face[1 + iTree * 6] = 0;
      conn->tree_to_face[2 + iTree * 6] = 8;
      conn->tree_to_face[3 + iTree * 6] = 3;
      conn->tree_to_face[4 + iTree * 6] = 5;
      conn->tree_to_face[5 + iTree * 6] = 4;
      iTree++;

      conn->tree_to_face[0 + iTree * 6] = 1;    /* tree = 2 modulo 5  */
      conn->tree_to_face[1 + iTree * 6] = 0;
      conn->tree_to_face[2 + iTree * 6] = 0;
      conn->tree_to_face[3 + iTree * 6] = 3;
      conn->tree_to_face[4 + iTree * 6] = 5;
      conn->tree_to_face[5 + iTree * 6] = 4;
      iTree++;

      conn->tree_to_face[0 + iTree * 6] = 1;    /* tree = 3 modulo 5  */
      conn->tree_to_face[1 + iTree * 6] = 0;
      conn->tree_to_face[2 + iTree * 6] = 3;
      conn->tree_to_face[3 + iTree * 6] = 3;
      conn->tree_to_face[4 + iTree * 6] = 5;
      conn->tree_to_face[5 + iTree * 6] = 4;
      iTree++;

      conn->tree_to_face[0 + iTree * 6] = 2;    /* tree = 4 modulo 5  */
      conn->tree_to_face[1 + iTree * 6] = 8;
      conn->tree_to_face[2 + iTree * 6] = 8;
      conn->tree_to_face[3 + iTree * 6] = 2;
      conn->tree_to_face[4 + iTree * 6] = 5;
      conn->tree_to_face[5 + iTree * 6] = 4;
      iTree++;

    }

  }

  /* tree to edge */
  {
    /* helper macro to compute global edge id */
#define eGlob(eLoc,iSeg) ( (eLoc) + (iSeg) * nEdgesPS )

    /* same as above but taking into account torus periodicity */
#define eGlob2(eLoc,iSeg) ( eGlob((eLoc),(iSeg)) < num_edges ? eGlob((eLoc),(iSeg)) : eGlob((eLoc),(iSeg)) - num_edges )

    /* all segments use the same pattern */
    iTree = 0;                  /*  global tree id */
    for (iSegment = 0; iSegment < nSegments; ++iSegment) {
      /* tree = 0 modulo 5  */
      conn->tree_to_edge[0 + iTree * 12] = eGlob (0, iSegment);
      conn->tree_to_edge[1 + iTree * 12] = -1;
      conn->tree_to_edge[2 + iTree * 12] = eGlob2 (0 + nEdgesPS, iSegment);
      conn->tree_to_edge[3 + iTree * 12] = -1;
      conn->tree_to_edge[4 + iTree * 12] = eGlob (7, iSegment);
      conn->tree_to_edge[5 + iTree * 12] = eGlob (4, iSegment);
      conn->tree_to_edge[6 + iTree * 12] = eGlob2 (7 + nEdgesPS, iSegment);
      conn->tree_to_edge[7 + iTree * 12] = eGlob2 (4 + nEdgesPS, iSegment);
      conn->tree_to_edge[8 + iTree * 12] = -1;
      conn->tree_to_edge[9 + iTree * 12] = -1;
      conn->tree_to_edge[10 + iTree * 12] = -1;
      conn->tree_to_edge[11 + iTree * 12] = -1;
      iTree++;

      /* tree = 1 modulo 5  */
      conn->tree_to_edge[0 + iTree * 12] = eGlob (1, iSegment);
      conn->tree_to_edge[1 + iTree * 12] = -1;
      conn->tree_to_edge[2 + iTree * 12] = eGlob2 (1 + nEdgesPS, iSegment);
      conn->tree_to_edge[3 + iTree * 12] = -1;
      conn->tree_to_edge[4 + iTree * 12] = eGlob (4, iSegment);
      conn->tree_to_edge[5 + iTree * 12] = eGlob (5, iSegment);
      conn->tree_to_edge[6 + iTree * 12] = eGlob2 (4 + nEdgesPS, iSegment);
      conn->tree_to_edge[7 + iTree * 12] = eGlob2 (5 + nEdgesPS, iSegment);
      conn->tree_to_edge[8 + iTree * 12] = -1;
      conn->tree_to_edge[9 + iTree * 12] = -1;
      conn->tree_to_edge[10 + iTree * 12] = -1;
      conn->tree_to_edge[11 + iTree * 12] = -1;
      iTree++;

      /* tree = 2 modulo 5  */
      conn->tree_to_edge[0 + iTree * 12] = eGlob (2, iSegment);
      conn->tree_to_edge[1 + iTree * 12] = -1;
      conn->tree_to_edge[2 + iTree * 12] = eGlob2 (2 + nEdgesPS, iSegment);
      conn->tree_to_edge[3 + iTree * 12] = -1;
      conn->tree_to_edge[4 + iTree * 12] = eGlob (5, iSegment);
      conn->tree_to_edge[5 + iTree * 12] = eGlob (6, iSegment);
      conn->tree_to_edge[6 + iTree * 12] = eGlob2 (5 + nEdgesPS, iSegment);
      conn->tree_to_edge[7 + iTree * 12] = eGlob2 (6 + nEdgesPS, iSegment);
      conn->tree_to_edge[8 + iTree * 12] = -1;
      conn->tree_to_edge[9 + iTree * 12] = -1;
      conn->tree_to_edge[10 + iTree * 12] = -1;
      conn->tree_to_edge[11 + iTree * 12] = -1;
      iTree++;

      /* tree = 3 modulo 5  */
      conn->tree_to_edge[0 + iTree * 12] = eGlob (3, iSegment);
      conn->tree_to_edge[1 + iTree * 12] = -1;
      conn->tree_to_edge[2 + iTree * 12] = eGlob2 (3 + nEdgesPS, iSegment);
      conn->tree_to_edge[3 + iTree * 12] = -1;
      conn->tree_to_edge[4 + iTree * 12] = eGlob (6, iSegment);
      conn->tree_to_edge[5 + iTree * 12] = eGlob (7, iSegment);
      conn->tree_to_edge[6 + iTree * 12] = eGlob2 (6 + nEdgesPS, iSegment);
      conn->tree_to_edge[7 + iTree * 12] = eGlob2 (7 + nEdgesPS, iSegment);
      conn->tree_to_edge[8 + iTree * 12] = -1;
      conn->tree_to_edge[9 + iTree * 12] = -1;
      conn->tree_to_edge[10 + iTree * 12] = -1;
      conn->tree_to_edge[11 + iTree * 12] = -1;
      iTree++;

      /* tree = 4 modulo 5  */
      conn->tree_to_edge[0 + iTree * 12] = eGlob (1, iSegment);
      conn->tree_to_edge[1 + iTree * 12] = eGlob (3, iSegment);
      conn->tree_to_edge[2 + iTree * 12] = eGlob2 (1 + nEdgesPS, iSegment);
      conn->tree_to_edge[3 + iTree * 12] = eGlob2 (3 + nEdgesPS, iSegment);
      conn->tree_to_edge[4 + iTree * 12] = eGlob (2, iSegment);
      conn->tree_to_edge[5 + iTree * 12] = eGlob (0, iSegment);
      conn->tree_to_edge[6 + iTree * 12] = eGlob2 (2 + nEdgesPS, iSegment);
      conn->tree_to_edge[7 + iTree * 12] = eGlob2 (0 + nEdgesPS, iSegment);
      conn->tree_to_edge[8 + iTree * 12] = -1;
      conn->tree_to_edge[9 + iTree * 12] = -1;
      conn->tree_to_edge[10 + iTree * 12] = -1;
      conn->tree_to_edge[11 + iTree * 12] = -1;
      iTree++;

    }

#undef eGlob2
#undef eGlob

#define tGlob(tLoc,iSeg) ( (tLoc) + (iSeg) * nTreesPS )

#define tGlob2(tLoc,iSeg) ( ( tGlob(tLoc,iSeg) < 0) ? tGlob(tLoc,iSeg) + num_trees : tGlob(tLoc,iSeg) )

    /* remember, there are 8 edges per segments */
    iEdge = 0;
    for (iSegment = 0; iSegment < nSegments; ++iSegment) {
      conn->edge_to_tree[iEdge + 0] = tGlob2 (0, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (4, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (0 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (4 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (1, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (4, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (1 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (4 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (2, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (4, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (2 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (4 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (3, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (4, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (3 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (4 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (0, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (1, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (0 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (1 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (1, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (2, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (1 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (2 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (2, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (3, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (2 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (3 - nTreesPS, iSegment);
      iEdge += 4;

      conn->edge_to_tree[iEdge + 0] = tGlob2 (3, iSegment);
      conn->edge_to_tree[iEdge + 1] = tGlob2 (0, iSegment);
      conn->edge_to_tree[iEdge + 2] = tGlob2 (3 - nTreesPS, iSegment);
      conn->edge_to_tree[iEdge + 3] = tGlob2 (0 - nTreesPS, iSegment);
      iEdge += 4;
    }

#undef tGlob2
#undef tGlob

    for (i = 0; i <= nEdgesPS * nSegments; ++i) {
      conn->ett_offset[i] = 4 * i;
    }

    /*  edge to edge */
    iEdge = 0;
    for (iSegment = 0; iSegment < nSegments; ++iSegment) {
      conn->edge_to_edge[iEdge + 0] = 0;
      conn->edge_to_edge[iEdge + 1] = 17;
      conn->edge_to_edge[iEdge + 2] = 2;
      conn->edge_to_edge[iEdge + 3] = 19;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 0;
      conn->edge_to_edge[iEdge + 1] = 12;
      conn->edge_to_edge[iEdge + 2] = 2;
      conn->edge_to_edge[iEdge + 3] = 14;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 0;
      conn->edge_to_edge[iEdge + 1] = 4;
      conn->edge_to_edge[iEdge + 2] = 2;
      conn->edge_to_edge[iEdge + 3] = 6;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 0;
      conn->edge_to_edge[iEdge + 1] = 1;
      conn->edge_to_edge[iEdge + 2] = 2;
      conn->edge_to_edge[iEdge + 3] = 3;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 5;
      conn->edge_to_edge[iEdge + 1] = 4;
      conn->edge_to_edge[iEdge + 2] = 7;
      conn->edge_to_edge[iEdge + 3] = 6;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 5;
      conn->edge_to_edge[iEdge + 1] = 4;
      conn->edge_to_edge[iEdge + 2] = 7;
      conn->edge_to_edge[iEdge + 3] = 6;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 5;
      conn->edge_to_edge[iEdge + 1] = 4;
      conn->edge_to_edge[iEdge + 2] = 7;
      conn->edge_to_edge[iEdge + 3] = 6;
      iEdge += 4;

      conn->edge_to_edge[iEdge + 0] = 5;
      conn->edge_to_edge[iEdge + 1] = 4;
      conn->edge_to_edge[iEdge + 2] = 7;
      conn->edge_to_edge[iEdge + 3] = 6;
      iEdge += 4;
    }

  }

  return conn;

}

static int
p8est_find_edge_transform_internal (p4est_connectivity_t * conn,
                                    p4est_topidx_t itree, int iedge,
                                    p8est_edge_info_t * ei,
                                    const p4est_topidx_t * ett,
                                    const int8_t * ete,
                                    p4est_topidx_t edge_trees)
{
  int                 i, j;
  int                 redge, nedge, iflip, nflip;
  int                 pref, pset, fc[2], nc[2];
  int                 face, nface, orient, eorient;
  p4est_topidx_t      etree, ietree, ntree;
  p8est_edge_transform_t *et;
  sc_array_t         *ta = &ei->edge_transforms;
  const int          *fcorners;
  int                 distinct = 1;
  int                 edges[3], edgeorients[3];
  p4est_topidx_t      etrees[3];

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= iedge && iedge < P8EST_EDGES);
  P4EST_ASSERT (ta->elem_size == sizeof (p8est_edge_transform_t));

  etrees[0] = itree;
  edges[0] = iedge;
  edgeorients[0] = 0;

  /* identify touching faces */
  for (i = 0; i < 2; ++i) {
    face = p8est_edge_faces[iedge][i];
    ntree = conn->tree_to_tree[P4EST_FACES * itree + face];
    nface = (int) conn->tree_to_face[P4EST_FACES * itree + face];
    if (ntree != itree || nface != face) {      /* not domain boundary */
      orient = nface / P4EST_FACES;
      nface %= P4EST_FACES;
      fcorners = &(p8est_edge_face_corners[iedge][face][0]);
      P4EST_ASSERT (fcorners[0] >= 0 && fcorners[1] >= 0);
      pref = p8est_face_permutation_refs[face][nface];
      pset = p8est_face_permutation_sets[pref][orient];
      fc[0] = p8est_face_permutations[pset][fcorners[0]];
      fc[1] = p8est_face_permutations[pset][fcorners[1]];

      /* if this is a new edge, add it */
      nc[0] = p8est_face_corners[nface][fc[0]];
      nc[1] = p8est_face_corners[nface][fc[1]];
      nedge = p8est_child_corner_edges[nc[0]][nc[1]];
      P4EST_ASSERT (nedge >= 0);
      eorient = (p8est_edge_corners[nedge][1] == nc[0]);
      for (j = 0; j < distinct; j++) {
        if (ntree == etrees[j] &&
            nedge == edges[j] && eorient == edgeorients[j]) {
          break;
        }
      }
      if (j == distinct) {
        etrees[j] = ntree;
        edges[j] = nedge;
        edgeorients[j] = eorient;
        distinct++;
      }
    }
  }

  /* find orientation of this edge */
  ietree = -1;
  iflip = -1;
  for (etree = 0; etree < edge_trees; ++etree) {
    ntree = ett[etree];
    P4EST_ASSERT (0 <= ntree && ntree < conn->num_trees);
    redge = (int) ete[etree];
    P4EST_ASSERT (redge >= 0 && redge < 2 * P8EST_EDGES);
    nedge = redge % P8EST_EDGES;
    if (nedge == iedge && ntree == itree) {
      iflip = redge / P8EST_EDGES;
      ietree = etree;
      break;
    }
  }
  P4EST_ASSERT (ietree >= 0 && iflip >= 0);

  /* loop through all trees connected through this edge */
  for (etree = 0; etree < edge_trees; ++etree) {
    if (etree == ietree) {
      continue;
    }
    ntree = ett[etree];
    P4EST_ASSERT (0 <= ntree && ntree < conn->num_trees);
    redge = (int) ete[etree];
    P4EST_ASSERT (redge >= 0 && redge < 2 * P8EST_EDGES);
    nedge = redge % P8EST_EDGES;
    nflip = (redge / P8EST_EDGES) ^ iflip;

    for (j = 0; j < distinct; j++) {
      if (ntree == etrees[j] && nedge == edges[j] && nflip == edgeorients[j]) {
        break;
      }
    }
    if (j < distinct) {
      /* already found from self or faces */
      continue;
    }
#if 0
    nows[0] = nows[1] = 0;
    for (i = 0; i < 2; ++i) {
      if (ntree == ntrees[i]) {
        /* check if the edge touches this neighbor contact face */
        P4EST_ASSERT (fcorners[i][0] >= 0);
        nfcorners = p8est_edge_face_corners[nedge][nfaces[i]];
        if (nfcorners[0] >= 0) {
          pref = p8est_face_permutation_refs[faces[i]][nfaces[i]];
          pset = p8est_face_permutation_sets[pref][orients[i]];
          fc[0] = p8est_face_permutations[pset][fcorners[i][0]];
          fc[1] = p8est_face_permutations[pset][fcorners[i][1]];

          if (fc[0] == nfcorners[nflip] && fc[1] == nfcorners[!nflip]) {
            P4EST_ASSERT (!founds[i] && !nows[!i]);
#ifdef P4EST_ENABLE_DEBUG
            founds[i] = 1;
#endif
            nows[i] = 1;
          }
          else if (fc[0] == nfcorners[!nflip] && fc[1] == nfcorners[nflip]) {
            ++flipped;
          }
        }
      }
    }
    if (nows[0] || nows[1]) {
      continue;
    }
#endif

    /* else we have a diagonal edge with ntree */
    et = (p8est_edge_transform_t *) sc_array_push (ta);
    et->ntree = ntree;
    et->nedge = (int8_t) nedge;
    et->naxis[0] = (int8_t) (nedge / 4);
    et->naxis[1] = (int8_t) (nedge < 4 ? 1 : 0);
    et->naxis[2] = (int8_t) (nedge < 8 ? 2 : 1);
    et->nflip = (int8_t) nflip;
    et->corners = (int8_t) (nedge % 4);
  }

  return distinct;
}

#include "p4est_connectivity.c"

int
p8est_connectivity_face_neighbor_face_edge (int fe, int f, int nf, int o)
{
  int                 nfe, pref, pset;

  P4EST_ASSERT (0 <= fe && fe < P4EST_HALF);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= o && o < P4EST_HALF);

  pref = p8est_face_permutation_refs[f][nf];
  pset = p8est_face_edge_permutation_sets[pref][o];
  nfe = p8est_face_edge_permutations[pset][fe];

  P4EST_ASSERT (0 <= nfe && nfe < P4EST_HALF);

  return nfe;
}

int
p8est_connectivity_face_neighbor_edge (int e, int f, int nf, int o)
{
  int                 fe, nfe;

  P4EST_ASSERT (0 <= e && e < P8EST_EDGES);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= o && o < P4EST_HALF);

  fe = p8est_edge_face_edges[e][f];
  P4EST_ASSERT (0 <= fe && fe < P4EST_HALF);

  nfe = p8est_connectivity_face_neighbor_face_edge (fe, f, nf, o);

  return p8est_face_edges[nf][nfe];
}

int
p8est_connectivity_edge_neighbor_edge_corner (int ec, int o)
{
  int                 nec;

  P4EST_ASSERT (0 <= ec && ec < 2);
  P4EST_ASSERT (0 <= o && o < 2);

  nec = ec ^ o;
  P4EST_ASSERT (0 <= nec && nec < 2);

  return nec;
}

int
p8est_connectivity_edge_neighbor_corner (int c, int e, int ne, int o)
{
  int                 ec, nec;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= e && e < P8EST_EDGES);

  ec = p8est_corner_edge_corners[c][e];
  P4EST_ASSERT (0 <= ec && ec < 2);

  nec = p8est_connectivity_edge_neighbor_edge_corner (ec, o);
  P4EST_ASSERT (0 <= nec && nec < 2);

  return p8est_edge_corners[ne][nec];
}

void
p8est_find_edge_transform (p4est_connectivity_t * conn,
                           p4est_topidx_t itree, int iedge,
                           p8est_edge_info_t * ei)
{
  p4est_topidx_t      edge_trees, aedge, ettae;
  sc_array_t         *ta = &ei->edge_transforms;

  P4EST_ASSERT (0 <= itree && itree < conn->num_trees);
  P4EST_ASSERT (0 <= iedge && iedge < P8EST_EDGES);
  P4EST_ASSERT (ta->elem_size == sizeof (p8est_edge_transform_t));

  /* check if this edge exists at all */
  ei->iedge = (int8_t) iedge;
  sc_array_resize (ta, 0);
  if (conn->num_edges == 0) {
    return;
  }
  aedge = conn->tree_to_edge[P8EST_EDGES * itree + iedge];
  if (aedge == -1) {
    return;
  }
  P4EST_ASSERT (0 <= aedge && aedge < conn->num_edges);

  /* retrieve connectivity information for this edge */
  ettae = conn->ett_offset[aedge];
  edge_trees = conn->ett_offset[aedge + 1] - ettae;
  P4EST_ASSERT (0 <= ettae && 1 <= edge_trees);

  /* loop through all edge neighbors and find edge connections */
  P4EST_EXECUTE_ASSERT_INT
    (p8est_find_edge_transform_internal (conn, itree, iedge, ei,
                                         conn->edge_to_tree + ettae,
                                         conn->edge_to_edge + ettae,
                                         edge_trees),
     (int) edge_trees - (int) ta->elem_count);
}

int
p8est_connect_type_int (p8est_connect_type_t btype)
{
  switch (btype) {
  case P8EST_CONNECT_FACE:
    return 1;
  case P8EST_CONNECT_EDGE:
    return 2;
  case P8EST_CONNECT_CORNER:
    return 3;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

const char         *
p8est_connect_type_string (p8est_connect_type_t btype)
{
  switch (btype) {
  case P8EST_CONNECT_FACE:
    return "FACE";
  case P8EST_CONNECT_EDGE:
    return "EDGE";
  case P8EST_CONNECT_CORNER:
    return "CORNER";
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
