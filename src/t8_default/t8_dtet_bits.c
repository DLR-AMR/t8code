/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include "t8_default_common.h"
#include "t8_dtet_bits.h"

typedef int8_t      t8_dtet_cube_id_t;

int                 t8_dtet_cid_type_to_parenttype[8][6] = {
  {0, 1, 2, 3, 4, 5},
  {0, 1, 1, 1, 0, 0},
  {2, 2, 2, 3, 3, 3},
  {1, 1, 2, 2, 2, 1},
  {5, 5, 4, 4, 4, 5},
  {0, 0, 0, 5, 5, 5},
  {4, 3, 3, 3, 4, 4},
  {0, 1, 2, 3, 4, 5}
};

/* In dependence of a type x give the type of
 * the child with Bey number y */
int                 t8_dtet_type_of_child[6][8] = {
  {0, 0, 0, 0, 4, 5, 2, 1},
  {1, 1, 1, 1, 3, 2, 5, 0},
  {2, 2, 2, 2, 0, 1, 4, 3},
  {3, 3, 3, 3, 5, 4, 1, 2},
  {4, 4, 4, 4, 2, 3, 0, 5},
  {5, 5, 5, 5, 1, 0, 3, 4}
};

static              t8_dtet_cube_id_t
t8_dtet_compute_cubeid (const t8_dtet_t * t, int level)
{
  t8_dtet_cube_id_t   id = 0;
  t8_dtet_coord_t     h;

  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  h = T8_DTET_ROOT_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((t->x & h) ? 0x01 : 0);
  id |= ((t->y & h) ? 0x02 : 0);
  id |= ((t->z & h) ? 0x04 : 0);

  return id;
}

int
t8_dtet_is_equal (const t8_dtet_t * t1, const t8_dtet_t * t2)
{
  return (t1->level == t1->level && t1->type == t2->type &&
          t1->x == t1->x && t1->y == t1->y && t1->z == t1->z);
}

void
t8_dtet_parent (const t8_dtet_t * t, t8_dtet_t * parent)
{
  t8_dtet_cube_id_t   cid;
  t8_dtet_coord_t     h;

  T8_ASSERT (t->level > 0);

  parent->eclass = t->eclass;
  parent->level = t->level - 1;

  /* Compute type of parent */
  cid = t8_dtet_compute_cubeid (t, t->level);
  parent->type = t8_dtet_cid_type_to_parenttype[cid][t->type];
  /* Set coordinates of parent */
  h = T8_DTET_ROOT_LEN (t->level);
  parent->x = t->x & ~h;
  parent->y = t->y & ~h;
  parent->z = t->z & ~h;
}

void
t8_dtet_compute_coords (const t8_dtet_t * t,
                        t8_dtet_coord_t coordinates[4][3])
{
  t8_dtet_type_t      type;
  int                 ei, ej;
  int                 i;
  t8_dtet_coord_t     h;

  type = t->type;
  h = T8_DTET_ROOT_LEN (t->level);
  ei = type / 2;
  if (type % 2 == 0) {
    ej = (ei + 2) % 3;
  }
  else {
    ej = (ei + 1) % 3;
  }

  coordinates[0][0] = t->x;
  coordinates[0][1] = t->y;
  coordinates[0][2] = t->z;
  for (i = 0; i < 2; i++) {
    coordinates[1][i] = coordinates[0][i];
    coordinates[2][i] = coordinates[0][i];
    coordinates[3][i] = coordinates[0][i] + h;
  }
  coordinates[1][ei] += h;
  coordinates[2][ei] += h;
  coordinates[2][ej] += h;
}

/* The childid here is the Bey child id,
 * not the Morton child id
 * It is possible that the function is called with
 * elem = child */
void
t8_dtet_child (const t8_dtet_t * elem, int childid, t8_dtet_t * child)
{
  const t8_dtet_t    *t = (const t8_dtet_t *) (elem);
  t8_dtet_t          *c = (t8_dtet_t *) (child);
  t8_dtet_coord_t     t_coordinates[4][3];
  t8_dtet_type_t      type;
  int                 coord2;

  T8_ASSERT (t->level < T8_DTET_MAXLEVEL);
  T8_ASSERT (0 <= childid && childid < 8);

  /* Compute anchor coordinates of child */
  if (childid == 0) {
    c->x = t->x;
    c->y = t->y;
    c->z = t->z;
  }
  else {
    switch (childid) {
    case 1:
    case 4:
    case 5:
      coord2 = 1;
      break;
    case 2:
    case 6:
    case 7:
      coord2 = 2;
      break;
    case 3:
      coord2 = 3;
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    /* i-th anchor coordinate of child is (X_(0,i)+X_(coord2,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    t8_dtet_compute_coords (t, t_coordinates);
    c->x = (t_coordinates[0][0] + t_coordinates[coord2][0]) >> 1;
    c->y = (t_coordinates[0][1] + t_coordinates[coord2][1]) >> 1;
    c->z = (t_coordinates[0][2] + t_coordinates[coord2][2]) >> 1;
  }

  /* Compute type of child */
  type = t->type;
  c->type = t8_dtet_type_of_child[type][childid];

  c->level = t->level + 1;
}

/* The sibid here is the Bey child id of the parent.
 * TODO: Implement this algorithm directly w/o using
 * parent and child */
void
t8_dtet_sibling (const t8_dtet_t * elem, int sibid, t8_dtet_t * sibling)
{
  T8_ASSERT (0 <= sibid && sibid < T8_DTET_CHILDREN);
  T8_ASSERT (((const t8_dtet_t *) (elem))->level > 0);
  t8_dtet_parent (elem, sibling);
  t8_dtet_child (sibling, sibid, sibling);
}

/* Saves the neighbour of T along face "face" in N
 * returns the facenumber of N along which T is its neighbour */
int
t8_dtet_face_neighbour (const t8_dtet_t * t, t8_dtet_t * n, int face)
{
  int                 type_new, type_old;
  int                 sign;
  int                 ret = -1;
  int8_t              level;
  t8_dtet_coord_t     coords[3];

  T8_ASSERT (0 <= face && face < 4);

  n->level = level = t->level;

  coords[0] = t->x;
  coords[1] = t->y;
  coords[2] = t->z;

  type_old = t->type;
  type_new = type_old;
  type_new += 6;                /* We want to compute modulo six and dont want negative numbers */
  if (face == 1 || face == 2) {
    sign = (type_new % 2 == 0 ? 1 : -1);
    sign *= (face % 2 == 0 ? 1 : -1);
    type_new += sign;
    type_new %= 6;
    ret = face;
  }
  else {
    if (face == 0) {
      /* type: 0,1 --> x+1
       *       2,3 --> y+1
       *       4,5 --> z+1 */
      coords[type_old / 2] += T8_DTET_ROOT_LEN (level);
      type_new += (type_new % 2 == 0 ? 4 : 2);
    }
    else {                      /* face == 3 */

      /* type: 1,2 --> z-1
       *       3,4 --> x-1
       *       5,0 --> y-1 */
      coords[((type_new + 3) % 6) / 2] -= T8_DTET_ROOT_LEN (level);
      type_new += (type_new % 2 == 0 ? 2 : 4);
    }
    type_new %= 6;
    ret = 3 - face;
  }

  n->x = coords[0];
  n->y = coords[1];
  n->z = coords[2];
  n->type = type_new;
  return ret;
}

/* we check if t1 and t2 lie in the same subcube and have
 * the same level and parent type */
int
t8_dtet_is_sibling (const t8_dtet_t * t1, const t8_dtet_t * t2)
{
  t8_dtet_coord_t     exclorx, exclory, exclorz;
  t8_dtet_cube_id_t   cid1, cid2;

  if (t1->level == 0) {
    return 0;
  }

  exclorx = t1->x ^ t2->x;
  exclory = t1->y ^ t2->y;
  exclorz = t1->z ^ t2->z;
  cid1 = t8_dtet_compute_cubeid (t1, t1->level);
  cid2 = t8_dtet_compute_cubeid (t2, t2->level);

  return
    (t1->level == t2->level) &&
    ((exclorx & ~T8_DTET_ROOT_LEN (t1->level)) == 0) &&
    ((exclory & ~T8_DTET_ROOT_LEN (t1->level)) == 0) &&
    ((exclorz & ~T8_DTET_ROOT_LEN (t1->level)) == 0) &&
    t8_dtet_cid_type_to_parenttype[cid1][t1->type] ==
    t8_dtet_cid_type_to_parenttype[cid2][t2->type] && t1->type != t2->type
    && 1;
}

int
t8_dtet_is_parent (const t8_dtet_t * t, const t8_dtet_t * c)
{
  t8_dtet_cube_id_t   cid;

  cid = t8_dtet_compute_cubeid (c, c->level);
  return
    (t->level + 1 == c->level) &&
    (t->x == (c->x & ~T8_DTET_ROOT_LEN (c->level))) &&
    (t->y == (c->y & ~T8_DTET_ROOT_LEN (c->level))) &&
    (t->z == (c->z & ~T8_DTET_ROOT_LEN (c->level))) &&
    t->type == t8_dtet_cid_type_to_parenttype[cid][c->type] && 1;
}
