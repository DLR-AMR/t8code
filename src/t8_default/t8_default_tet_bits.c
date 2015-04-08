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
#include "t8_default_tet_bits.h"

typedef int8_t      t8_default_tet_cube_id_t;

int                 t8_tet_cid_type_to_parenttype[8][6] = {
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
int                 t8_tet_type_of_child[6][8] = {
  {0, 0, 0, 0, 4, 5, 2, 1},
  {1, 1, 1, 1, 3, 2, 5, 0},
  {2, 2, 2, 2, 0, 1, 4, 3},
  {3, 3, 3, 3, 5, 4, 1, 2},
  {4, 4, 4, 4, 2, 3, 0, 5},
  {5, 5, 5, 5, 1, 0, 3, 4}
};

static              t8_default_tet_cube_id_t
t8_default_tet_compute_cubeid (const t8_tet_t * t, int level)
{
  t8_default_tet_cube_id_t id = 0;
  t8_tcoord_t         h;

  T8_ASSERT (0 <= level && level <= T8_TET_MAX_LEVEL);
  h = T8_TET_ROOT_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((t->anchor_coordinates[0] & h) ? 0x01 : 0);
  id |= ((t->anchor_coordinates[1] & h) ? 0x02 : 0);
  id |= ((t->anchor_coordinates[2] & h) ? 0x04 : 0);

  return id;
}

void
t8_default_tet_parent (const t8_element_t * elem, t8_element_t * parent)
{
  const t8_tet_t     *t = (const t8_tet_t *) elem;
  t8_tet_t           *p = (t8_tet_t *) parent;
  t8_default_tet_cube_id_t cid;
  t8_tcoord_t         h;
  int                 i;

  T8_ASSERT (t->level > 0);

  p->eclass = t->eclass;
  p->level = t->level - 1;

  /* Compute type of parent */
  cid = t8_default_tet_compute_cubeid (t, t->level);
  t8_default_tet_set_type (p, t8_tet_cid_type_to_parenttype[cid]
                           [t8_default_tet_get_type (t)]);
  /* Set coordinates of parent */
  h = T8_TET_ROOT_LEN (t->level);
  for (i = 0; i < 3; i++) {
    p->anchor_coordinates[i] = t->anchor_coordinates[i] & ~h;
  }
}

void
t8_default_tet_compute_coords (const t8_tet_t * t,
                               t8_tcoord_t coordinates[4][3])
{
  t8_default_tet_type_t type;
  int                 ei, ej;
  int                 i;
  t8_tcoord_t         h;

  type = t8_default_tet_get_type (t);
  h = T8_TET_ROOT_LEN (t->level);
  ei = type / 2;
  if (type % 2 == 0) {
    ej = (ei + 2) % 3;
  }
  else {
    ej = (ei + 1) % 3;
  }
  for (i = 0; i < 2; i++) {
    coordinates[0][i] = t8_default_tet_get_coordinate (t, i);
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
t8_default_tet_child (const t8_element_t * elem, int childid,
                      t8_element_t * child)
{
  const t8_tet_t     *t = (const t8_tet_t *) (elem);
  t8_tet_t           *c = (t8_tet_t *) (child);
  t8_tcoord_t         t_coordinates[4][3], temp_coord;
  t8_default_tet_type_t type;
  int                 coord2, i;

  T8_ASSERT (t->level < T8_TET_MAX_LEVEL);
  T8_ASSERT (0 <= childid && childid < 8);

  /* Compute anchor coordinates of child */
  if (childid == 0) {
    for (i = 0; i < 3; i++) {
      t8_default_tet_set_coordinate (c, i,
                                     t8_default_tet_get_coordinate (t, i));
    }
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
    t8_default_tet_compute_coords (t, t_coordinates);
    for (i = 0; i < 3; i++) {
      temp_coord = (t_coordinates[0][i] + t_coordinates[coord2][i]) >> 1;
      t8_default_tet_set_coordinate (c, i, temp_coord);
    }
  }

  /* Compute type of child */
  type = t8_default_tet_get_type (t);
  t8_default_tet_set_type (c, t8_tet_type_of_child[type][childid]);

  c->level = t->level + 1;
}

/* The sibid here is the Bey child id of the parent.
 * TODO: Implement this algorithm directly w/o using
 * parent and child */
void
t8_default_tet_sibling (const t8_element_t * elem, int sibid,
                        t8_element_t * sibling)
{
  T8_ASSERT (0 <= sibid && sibid < T8_TET_CHILDREN);
  T8_ASSERT (((const t8_tet_t *) (elem))->level > 0);
  t8_default_tet_parent (elem, sibling);
  t8_default_tet_child (sibling, sibid, sibling);
}

/* Saves the neighbour of T along face "face" in N
 * returns the facenumber of N along which T is its neighbour */
int
t8_default_tet_face_neighbour (const t8_tet_t * t, t8_tet_t * n, int face)
{
  int                 type_new, type_old;
  int                 i, sign;
  int                 ret = -1;
  int8_t              level;
  t8_tcoord_t         coords[3];

  T8_ASSERT (0 <= face && face < 4);

  n->level = level = t->level;

  for (i = 0; i < 3; i++) {
    coords[i] = t8_default_tet_get_coordinate (t, i);
  }

  type_old = t8_default_tet_get_type (t);
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
      coords[type_old / 2] += T8_TET_ROOT_LEN (level);
      type_new += (type_new % 2 == 0 ? 4 : 2);
    }
    else {                      /* face == 3 */

      /* type: 1,2 --> z-1
       *       3,4 --> x-1
       *       5,0 --> y-1 */
      coords[((type_new + 3) % 6) / 2] -= T8_TET_ROOT_LEN (level);
      type_new += (type_new % 2 == 0 ? 2 : 4);
    }
    type_new %= 6;
    ret = 3 - face;
  }

  for (i = 0; i < 3; i++) {
    t8_default_tet_set_coordinate (n, i, coords[i]);
  }
  t8_default_tet_set_type (n, type_new);
  return ret;
}
