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
#include "p4est_bits.c"

int
p8est_quadrant_is_outside_edge (const p4est_quadrant_t * q)
{
  int                 outface[P4EST_DIM];

  outface[0] = (int) (q->x < 0 || q->x >= P4EST_ROOT_LEN);
  outface[1] = (int) (q->y < 0 || q->y >= P4EST_ROOT_LEN);
  outface[2] = (int) (q->z < 0 || q->z >= P4EST_ROOT_LEN);

  return outface[0] + outface[1] + outface[2] == 2;
}

int
p8est_quadrant_is_outside_edge_extra (const p4est_quadrant_t * q, int *edge)
{
  int                 quad_contact[P4EST_FACES];
  int                 face_axis[P4EST_DIM];

  P4EST_ASSERT (q->level <= P4EST_QMAXLEVEL);

  quad_contact[0] = (int) (q->x < 0);
  quad_contact[1] = (int) (q->x >= P4EST_ROOT_LEN);
  quad_contact[2] = (int) (q->y < 0);
  quad_contact[3] = (int) (q->y >= P4EST_ROOT_LEN);
  quad_contact[4] = (int) (q->z < 0);
  quad_contact[5] = (int) (q->z >= P4EST_ROOT_LEN);
  face_axis[0] = quad_contact[0] || quad_contact[1];
  face_axis[1] = quad_contact[2] || quad_contact[3];
  face_axis[2] = quad_contact[4] || quad_contact[5];

  if (face_axis[0] + face_axis[1] + face_axis[2] != 2) {
    return 0;
  }

  if (edge != NULL) {
    if (!face_axis[0]) {
      *edge = 0 + 2 * quad_contact[5] + quad_contact[3];
    }
    else if (!face_axis[1]) {
      *edge = 4 + 2 * quad_contact[5] + quad_contact[1];
    }
    else if (!face_axis[2]) {
      *edge = 8 + 2 * quad_contact[3] + quad_contact[1];
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
    P4EST_ASSERT (p8est_quadrant_touches_edge (q, *edge, 0));
  }

  return 1;
}

int
p4est_quadrant_is_family (const p4est_quadrant_t * q0,
                          const p4est_quadrant_t * q1,
                          const p4est_quadrant_t * q2,
                          const p4est_quadrant_t * q3,
                          const p4est_quadrant_t * q4,
                          const p4est_quadrant_t * q5,
                          const p4est_quadrant_t * q6,
                          const p4est_quadrant_t * q7)
{
  const int8_t        level = q0->level;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (p4est_quadrant_is_extended (q0));
  P4EST_ASSERT (p4est_quadrant_is_extended (q1));
  P4EST_ASSERT (p4est_quadrant_is_extended (q2));
  P4EST_ASSERT (p4est_quadrant_is_extended (q3));
  P4EST_ASSERT (p4est_quadrant_is_extended (q4));
  P4EST_ASSERT (p4est_quadrant_is_extended (q5));
  P4EST_ASSERT (p4est_quadrant_is_extended (q6));
  P4EST_ASSERT (p4est_quadrant_is_extended (q7));

  if (level == 0 || level != q1->level ||
      level != q2->level || level != q3->level ||
      level != q4->level || level != q5->level ||
      level != q6->level || level != q7->level) {
    return 0;
  }

  inc = P4EST_QUADRANT_LEN (level);
  return ((q0->x + inc == q1->x && q0->y == q1->y && q0->z == q1->z) &&
          (q0->x == q2->x && q0->y + inc == q2->y && q0->z == q2->z) &&
          (q1->x == q3->x && q2->y == q3->y && q0->z == q3->z) &&
          (q0->x == q4->x && q0->y == q4->y && q0->z + inc == q4->z) &&
          (q1->x == q5->x && q1->y == q5->y && q4->z == q5->z) &&
          (q2->x == q6->x && q2->y == q6->y && q4->z == q6->z) &&
          (q3->x == q7->x && q3->y == q7->y && q4->z == q7->z));
}

void
p8est_quadrant_edge_neighbor (const p4est_quadrant_t * q,
                              int edge, p4est_quadrant_t * r)
{
  const p4est_qcoord_t qh = P4EST_QUADRANT_LEN (q->level);

  P4EST_ASSERT (0 <= edge && edge < 12);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  switch (edge / 4) {
  case 0:
    r->x = q->x;
    r->y = q->y + (2 * (edge & 0x01) - 1) * qh;
    r->z = q->z + ((edge & 0x02) - 1) * qh;
    break;
  case 1:
    r->x = q->x + (2 * (edge & 0x01) - 1) * qh;
    r->y = q->y;
    r->z = q->z + ((edge & 0x02) - 1) * qh;
    break;
  case 2:
    r->x = q->x + (2 * (edge & 0x01) - 1) * qh;
    r->y = q->y + ((edge & 0x02) - 1) * qh;
    r->z = q->z;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  r->level = q->level;
  P4EST_ASSERT (p4est_quadrant_is_extended (r));
}

void
p8est_quadrant_edge_neighbor_extra (const p4est_quadrant_t * q, p4est_topidx_t
                                    t, int edge, sc_array_t * quads,
                                    sc_array_t * treeids, sc_array_t * nedges,
                                    p4est_connectivity_t * conn)
{
  p4est_quadrant_t    temp;
  p4est_quadrant_t   *qp;
  p4est_topidx_t     *tp;
  int                 face;
  int                *ip;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
  size_t              etree;

  eta = &ei.edge_transforms;

  P4EST_ASSERT (SC_ARRAY_IS_OWNER (quads));
  P4EST_ASSERT (quads->elem_count == 0);
  P4EST_ASSERT (quads->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (treeids));
  P4EST_ASSERT (treeids->elem_count == 0);
  P4EST_ASSERT (treeids->elem_size == sizeof (p4est_topidx_t));
  if (nedges != NULL) {
    P4EST_ASSERT (SC_ARRAY_IS_OWNER (nedges));
    P4EST_ASSERT (nedges->elem_count == 0);
    P4EST_ASSERT (nedges->elem_size == sizeof (int));
  }

  p8est_quadrant_edge_neighbor (q, edge, &temp);
  if (p4est_quadrant_is_inside_root (&temp)) {
    qp = p4est_quadrant_array_push_copy (quads, &temp);
    tp = (p4est_topidx_t *) sc_array_push (treeids);
    *tp = t;
    if (nedges != NULL) {
      ip = (int *) sc_array_push (nedges);
      *ip = (edge ^ 3);
    }
    return;
  }

  if (!p8est_quadrant_is_outside_edge (&temp)) {
    qp = p4est_quadrant_array_push (quads);
    tp = (p4est_topidx_t *) sc_array_push (treeids);

    face = p8est_edge_faces[edge][0];
    p4est_quadrant_face_neighbor (q, face, &temp);
    if (p4est_quadrant_is_inside_root (&temp)) {
      face = p8est_edge_faces[edge][1];
      *tp = p8est_quadrant_face_neighbor_extra (&temp, t, face, qp, NULL,
                                                conn);
      if (*tp == -1) {
        qp = (p4est_quadrant_t *) sc_array_pop (quads);
        tp = (p4est_topidx_t *) sc_array_pop (treeids);
      }
      else if (nedges != NULL) {
        int                 opedge = (edge ^ 1);
        int                 nface =
          conn->tree_to_face[P4EST_FACES * t + face];
        int                 o = nface / P4EST_FACES;
        int                 ref, set;
        int                 c1, c2, nc1, nc2;

        nface = nface % P4EST_FACES;

        P4EST_ASSERT (p8est_edge_faces[opedge][1] == face);
        ref = p8est_face_permutation_refs[face][nface];
        set = p8est_face_permutation_sets[ref][o];

        c1 = p8est_edge_corners[opedge][0];
        nc1 = p8est_corner_face_corners[c1][face];
        P4EST_ASSERT (nc1 >= 0);
        c1 = p8est_face_permutations[set][nc1];
        nc1 = p8est_face_corners[nface][c1];

        c2 = p8est_edge_corners[opedge][1];
        nc2 = p8est_corner_face_corners[c2][face];
        P4EST_ASSERT (nc2 >= 0);
        c2 = p8est_face_permutations[set][nc2];
        nc2 = p8est_face_corners[nface][c2];

        P4EST_ASSERT (nc1 >= 0);
        ip = (int *) sc_array_push (nedges);
        *ip = p8est_child_corner_edges[nc1][nc2];
        if (nc1 > nc2) {
          *ip += P8EST_EDGES;
        }
        P4EST_ASSERT (*ip >= 0);
      }
      return;
    }
    face = p8est_edge_faces[edge][1];
    p4est_quadrant_face_neighbor (q, face, &temp);
    P4EST_ASSERT (p4est_quadrant_is_inside_root (&temp));
    face = p8est_edge_faces[edge][0];
    *tp = p8est_quadrant_face_neighbor_extra (&temp, t, face, qp, NULL, conn);
    if (*tp == -1) {
      qp = (p4est_quadrant_t *) sc_array_pop (quads);
      tp = (p4est_topidx_t *) sc_array_pop (treeids);
    }
    else if (nedges != NULL) {
      int                 opedge = (edge ^ 2);
      int                 nface = conn->tree_to_face[P4EST_FACES * t + face];
      int                 o = nface / P4EST_FACES;
      int                 ref, set;
      int                 c1, c2, nc1, nc2;

      nface = nface % P4EST_FACES;

      P4EST_ASSERT (p8est_edge_faces[opedge][0] == face);
      ref = p8est_face_permutation_refs[face][nface];
      set = p8est_face_permutation_sets[ref][o];

      c1 = p8est_edge_corners[opedge][0];
      nc1 = p8est_corner_face_corners[c1][face];
      P4EST_ASSERT (nc1 >= 0);
      c1 = p8est_face_permutations[set][nc1];
      nc1 = p8est_face_corners[nface][c1];

      c2 = p8est_edge_corners[opedge][1];
      nc2 = p8est_corner_face_corners[c2][face];
      P4EST_ASSERT (nc2 >= 0);
      c2 = p8est_face_permutations[set][nc2];
      nc2 = p8est_face_corners[nface][c2];

      P4EST_ASSERT (nc1 >= 0);
      ip = (int *) sc_array_push (nedges);
      *ip = p8est_child_corner_edges[nc1][nc2];
      if (nc1 > nc2) {
        *ip += P8EST_EDGES;
      }
      P4EST_ASSERT (*ip >= 0);
    }
    return;
  }
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  p8est_find_edge_transform (conn, t, edge, &ei);
  sc_array_resize (quads, eta->elem_count);
  sc_array_resize (treeids, eta->elem_count);
  if (nedges != NULL) {
    sc_array_resize (nedges, eta->elem_count);
  }
  for (etree = 0; etree < eta->elem_count; etree++) {
    qp = p4est_quadrant_array_index (quads, etree);
    tp = (p4est_topidx_t *) sc_array_index (treeids, etree);
    et = p8est_edge_array_index (eta, etree);
    p8est_quadrant_transform_edge (&temp, qp, &ei, et, 1);
    *tp = et->ntree;
    if (nedges != NULL) {
      ip = (int *) sc_array_index (nedges, etree);
      *ip = et->nedge;
      if (et->nflip) {
        *ip += P8EST_EDGES;
      }
    }
  }
  sc_array_reset (eta);
}

int
p8est_quadrant_touches_edge (const p4est_quadrant_t * q, int edge, int inside)
{
  int                 quad_contact[P4EST_FACES];
  int                 axis, side, incount;
  p4est_qcoord_t      lower, upper;

  P4EST_ASSERT (0 <= edge && edge < 12);

  axis = edge / 4;
  if (q->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (p4est_quadrant_is_node (q, inside));
    lower = 0;
    upper = P4EST_ROOT_LEN - (int) inside;
  }
  else {
    if (!inside) {
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      lower = -P4EST_QUADRANT_LEN (q->level);
      upper = P4EST_ROOT_LEN;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_is_valid (q));
      lower = 0;
      upper = P4EST_LAST_OFFSET (q->level);
    }
  }
  quad_contact[0] = (q->x == lower);
  quad_contact[1] = (q->x == upper);
  quad_contact[2] = (q->y == lower);
  quad_contact[3] = (q->y == upper);
  quad_contact[4] = (q->z == lower);
  quad_contact[5] = (q->z == upper);

  incount = 0;
  if (axis != 0) {
    side = edge & 1;
    incount += quad_contact[side];
  }
  if (axis != 1) {
    side = (axis == 0) ? (edge & 1) : ((edge >> 1) & 1);
    incount += quad_contact[2 + side];
  }
  if (axis != 2) {
    side = (edge >> 1) & 1;
    incount += quad_contact[4 + side];
  }
#ifdef P4EST_ENABLE_DEBUG
  upper = P4EST_ROOT_LEN + (p4est_qcoord_t) (q->level == P4EST_MAXLEVEL
                                             && !inside);
  P4EST_ASSERT (axis != 0 || (q->x >= 0 && q->x < upper));
  P4EST_ASSERT (axis != 1 || (q->y >= 0 && q->y < upper));
  P4EST_ASSERT (axis != 2 || (q->z >= 0 && q->z < upper));
#endif

  return incount == 2;
}

void
p8est_coordinates_transform_edge (const p4est_qcoord_t coords_in[],
                                  p4est_qcoord_t coords_out[],
                                  const p8est_edge_info_t * ei,
                                  const p8est_edge_transform_t * et)
{
  int                 iaxis;
  p4est_qcoord_t      my_xyz, *target_xyz[3];

  iaxis = (int) ei->iedge / 4;
  P4EST_ASSERT (0 <= et->naxis[0] && et->naxis[0] < 3);
  P4EST_ASSERT (0 <= et->naxis[1] && et->naxis[1] < 3);
  P4EST_ASSERT (0 <= et->naxis[2] && et->naxis[2] < 3);
  P4EST_ASSERT (et->naxis[0] != et->naxis[1] &&
                et->naxis[0] != et->naxis[2] && et->naxis[1] != et->naxis[2]);
  P4EST_ASSERT (0 <= et->nflip && et->nflip < 2);
  P4EST_ASSERT (coords_in != coords_out);

  for (int d = 0; d < P4EST_DIM; d++) {
    target_xyz[d] = &coords_out[d];
  }

  /* transform coordinate axis parallel to edge */
  my_xyz = coords_in[iaxis];
  if (!et->nflip) {
    *target_xyz[et->naxis[0]] = my_xyz;
  }
  else {
    *target_xyz[et->naxis[0]] = P4EST_ROOT_LEN - my_xyz;
  }

  /* create the other two coordinates */
  switch (et->corners) {
  case 0:
    *target_xyz[et->naxis[1]] = 0;
    *target_xyz[et->naxis[2]] = 0;
    break;
  case 1:
    *target_xyz[et->naxis[1]] = P4EST_ROOT_LEN;
    *target_xyz[et->naxis[2]] = 0;
    break;
  case 2:
    *target_xyz[et->naxis[1]] = 0;
    *target_xyz[et->naxis[2]] = P4EST_ROOT_LEN;
    break;
  case 3:
    *target_xyz[et->naxis[1]] = P4EST_ROOT_LEN;
    *target_xyz[et->naxis[2]] = P4EST_ROOT_LEN;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
p8est_quadrant_transform_edge (const p4est_quadrant_t * q,
                               p4est_quadrant_t * r,
                               const p8est_edge_info_t * ei,
                               const p8est_edge_transform_t * et, int inside)
{
  int                 iaxis;
  p4est_qcoord_t      mh, Rmh;
  p4est_qcoord_t      lshift, rshift;
  p4est_qcoord_t      my_xyz, *target_xyz[3];

  iaxis = (int) ei->iedge / 4;
  P4EST_ASSERT (0 <= et->naxis[0] && et->naxis[0] < 3);
  P4EST_ASSERT (0 <= et->naxis[1] && et->naxis[1] < 3);
  P4EST_ASSERT (0 <= et->naxis[2] && et->naxis[2] < 3);
  P4EST_ASSERT (et->naxis[0] != et->naxis[1] &&
                et->naxis[0] != et->naxis[2] && et->naxis[1] != et->naxis[2]);
  P4EST_ASSERT (0 <= et->nflip && et->nflip < 2);
  P4EST_ASSERT (q != r);

  if (q->level == P4EST_MAXLEVEL) {
    P4EST_ASSERT (!inside);
    P4EST_ASSERT (p8est_quadrant_touches_edge (q, (int) ei->iedge, inside));
    lshift = mh = 0;
    rshift = Rmh = P4EST_ROOT_LEN;
  }
  else {
    P4EST_ASSERT (p8est_quadrant_touches_edge (q, (int) ei->iedge, !inside));
    mh = -P4EST_QUADRANT_LEN (q->level);
    Rmh = P4EST_ROOT_LEN + mh;
    lshift = (inside ? 0 : mh);
    rshift = (inside ? Rmh : P4EST_ROOT_LEN);
  }
  target_xyz[0] = &r->x;
  target_xyz[1] = &r->y;
  target_xyz[2] = &r->z;

  /* transform coordinate axis parallel to edge */
  switch (iaxis) {
  case 0:
    my_xyz = q->x;
    break;
  case 1:
    my_xyz = q->y;
    break;
  case 2:
    my_xyz = q->z;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  if (!et->nflip) {
    *target_xyz[et->naxis[0]] = my_xyz;
  }
  else {
    *target_xyz[et->naxis[0]] = Rmh - my_xyz;
  }

  /* create the other two coordinates */
  switch (et->corners) {
  case 0:
    *target_xyz[et->naxis[1]] = lshift;
    *target_xyz[et->naxis[2]] = lshift;
    break;
  case 1:
    *target_xyz[et->naxis[1]] = rshift;
    *target_xyz[et->naxis[2]] = lshift;
    break;
  case 2:
    *target_xyz[et->naxis[1]] = lshift;
    *target_xyz[et->naxis[2]] = rshift;
    break;
  case 3:
    *target_xyz[et->naxis[1]] = rshift;
    *target_xyz[et->naxis[2]] = rshift;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

#ifdef P4EST_ENABLE_DEBUG
  {
    /* This is the code from the paper. */

    p4est_qcoord_t      qparallel, qt1, qt2;

    qparallel = et->nflip * Rmh + (1 - 2 * et->nflip) * my_xyz;
    P4EST_ASSERT (qparallel == *target_xyz[et->naxis[0]]);

    qt1 = (!(et->nedge & 1)) ? lshift : rshift;
    qt2 = (!(et->nedge & 2)) ? lshift : rshift;
    P4EST_ASSERT (qt1 == *target_xyz[et->naxis[1]]);
    P4EST_ASSERT (qt2 == *target_xyz[et->naxis[2]]);
  }
#endif

  r->level = q->level;
  P4EST_ASSERT (p8est_quadrant_touches_edge (r, (int) et->nedge, inside));
}

void
p8est_quadrant_shift_edge (const p4est_quadrant_t * q,
                           p4est_quadrant_t * r, p4est_quadrant_t * rup,
                           p4est_quadrant_t * rdown, int edge)
{
  int                 outface;
  int                 i, level;
  int                 cid, sid[P4EST_DIM], step[P4EST_DIM];
  p4est_qcoord_t      th;
  p4est_quadrant_t    quad[P4EST_DIM];
  /* *INDENT-OFF* */
  const int           contact[12] = {
    0x14, 0x18, 0x24, 0x28,
    0x11, 0x12, 0x21, 0x22,
    0x05, 0x06, 0x09, 0x0a
  };
  /* *INDENT-ON* */

  P4EST_ASSERT (q != r);
  P4EST_ASSERT (rup == NULL || q != rup);
  P4EST_ASSERT (rdown == NULL || q != rdown);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (edge >= 0 && edge < 12);

  P4EST_QUADRANT_INIT (&quad[0]);
  P4EST_QUADRANT_INIT (&quad[1]);
  P4EST_QUADRANT_INIT (&quad[2]);

  quad[0] = *q;
  quad[1] = *q;
  quad[2] = *q;
  for (;;) {
    th = P4EST_LAST_OFFSET (quad[0].level);
    cid = p4est_quadrant_child_id (&quad[1]);
    switch (edge / 4) {
    case 0:
      sid[0] = 2 * edge;
      sid[1] = 2 * edge + (cid & 0x01);
      sid[2] = 2 * edge + 1;
      step[0] = 0;
      step[1] = 2 * (edge & 0x01) - 1;
      step[2] = (edge & 0x02) - 1;
      break;
    case 1:
      sid[0] = 2 * (edge & 0x02) + (edge & 0x01);
      sid[1] = 2 * (edge & 0x02) + (edge & 0x01) + (cid & 0x02);
      sid[2] = 2 * (edge & 0x02) + (edge & 0x01) + 2;
      step[0] = 2 * (edge & 0x01) - 1;
      step[1] = 0;
      step[2] = (edge & 0x02) - 1;
      break;
    case 2:
      sid[0] = edge - 8;
      sid[1] = edge - 8 + (cid & 0x04);
      sid[2] = edge - 8 + 4;
      step[0] = 2 * (edge & 0x01) - 1;
      step[1] = (edge & 0x02) - 1;
      step[2] = 0;
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    p4est_quadrant_sibling (&quad[1], r, sid[1]);
    if (rup != NULL) {
      p4est_quadrant_sibling (&quad[0], rup, sid[0]);
    }
    if (rdown != NULL) {
      p4est_quadrant_sibling (&quad[2], rdown, sid[2]);
    }
    P4EST_ASSERT (-1 <= step[0] && step[0] <= 1);
    P4EST_ASSERT (-1 <= step[1] && step[1] <= 1);
    P4EST_ASSERT (-1 <= step[2] && step[2] <= 1);

    outface = 0;
    if (step[0] != 0) {
      outface |= ((r->x <= 0) ? 0x01 : 0);
      outface |= ((r->x >= th) ? 0x02 : 0);
    }
    if (step[1] != 0) {
      outface |= ((r->y <= 0) ? 0x04 : 0);
      outface |= ((r->y >= th) ? 0x08 : 0);
    }
    if (step[2] != 0) {
      outface |= ((r->z <= 0) ? 0x10 : 0);
      outface |= ((r->z >= th) ? 0x20 : 0);
    }
    if (outface == contact[edge]) {
      break;
    }
    level = quad[0].level - 1;
    for (i = 0; i < P4EST_DIM; i++) {
      p4est_quadrant_parent (&quad[i], &quad[i]);
      quad[i].x += (p4est_qcoord_t) step[0] * P4EST_QUADRANT_LEN (level);
      quad[i].y += (p4est_qcoord_t) step[1] * P4EST_QUADRANT_LEN (level);
      quad[i].z += (p4est_qcoord_t) step[2] * P4EST_QUADRANT_LEN (level);
    }
    switch (edge / 4) {
    case 0:
      quad[0].x += P4EST_QUADRANT_LEN (level);
      quad[2].x -= P4EST_QUADRANT_LEN (level);
      break;
    case 1:
      quad[0].y += P4EST_QUADRANT_LEN (level);
      quad[2].y -= P4EST_QUADRANT_LEN (level);
      break;
    case 2:
      quad[0].z += P4EST_QUADRANT_LEN (level);
      quad[2].z -= P4EST_QUADRANT_LEN (level);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }

    P4EST_ASSERT (p4est_quadrant_is_extended (&quad[0]));
    P4EST_ASSERT (p4est_quadrant_is_extended (&quad[1]));
    P4EST_ASSERT (p4est_quadrant_is_extended (&quad[2]));
  }

  if (step[0] != 0) {
    if (r->x < 0)
      r->x = 0;
    if (r->x >= P4EST_ROOT_LEN)
      r->x = th;
  }
  if (rup != NULL && rup->x < 0)
    rup->x = 0;
  if (rup != NULL && rup->x >= P4EST_ROOT_LEN)
    rup->x = th;
  if (rdown != NULL && rdown->x < 0)
    rdown->x = 0;
  if (rdown != NULL && rdown->x >= P4EST_ROOT_LEN)
    rdown->x = th;
  if (step[1] != 0) {
    if (r->y < 0)
      r->y = 0;
    if (r->y >= P4EST_ROOT_LEN)
      r->y = th;
  }
  if (rup != NULL && rup->y < 0)
    rup->y = 0;
  if (rup != NULL && rup->y >= P4EST_ROOT_LEN)
    rup->y = th;
  if (rdown != NULL && rdown->y < 0)
    rdown->y = 0;
  if (rdown != NULL && rdown->y >= P4EST_ROOT_LEN)
    rdown->y = th;
  if (step[2] != 0) {
    if (r->z < 0)
      r->z = 0;
    if (r->z >= P4EST_ROOT_LEN)
      r->z = th;
  }
  if (rup != NULL && rup->z < 0)
    rup->z = 0;
  if (rup != NULL && rup->z >= P4EST_ROOT_LEN)
    rup->z = th;
  if (rdown != NULL && rdown->z < 0)
    rdown->z = 0;
  if (rdown != NULL && rdown->z >= P4EST_ROOT_LEN)
    rdown->z = th;
  P4EST_ASSERT (p8est_quadrant_touches_edge (r, edge, 1));
  P4EST_ASSERT (rup == NULL || p8est_quadrant_touches_edge (rup, edge, 1));
  P4EST_ASSERT (rdown == NULL
                || p8est_quadrant_touches_edge (rdown, edge, 1));
}
