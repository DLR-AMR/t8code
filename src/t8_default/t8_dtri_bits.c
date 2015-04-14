/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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
#ifndef T8_DTRI_TO_DTET
#include "t8_dtri_connectivity.h"
#include "t8_dtri_bits.h"
#else
#include "t8_dtet_connectivity.h"
#include "t8_dtet_bits.h"
#endif

typedef int8_t      t8_dtri_cube_id_t;

static              t8_dtri_cube_id_t
compute_cubeid (const t8_dtri_t * t, int level)
{
  t8_dtri_cube_id_t   id = 0;
  t8_dtri_coord_t     h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  h = T8_DTRI_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((t->x & h) ? 0x01 : 0);
  id |= ((t->y & h) ? 0x02 : 0);
#ifdef T8_DTRI_TO_DTET
  id |= ((t->z & h) ? 0x04 : 0);
#endif

  return id;
}

void
t8_dtri_parent (const t8_dtri_t * t, t8_dtri_t * parent)
{
  t8_dtri_cube_id_t   cid;
  t8_dtri_coord_t     h;

  T8_ASSERT (t->level > 0);

#ifdef T8_DTRI_TO_DTET
  parent->eclass = t->eclass;
#endif
  parent->level = t->level - 1;

  /* Compute type of parent */
  cid = compute_cubeid (t, t->level);
  parent->type = t8_dtri_cid_type_to_parenttype[cid][t->type];
  /* Set coordinates of parent */
  h = T8_DTRI_LEN (t->level);
  parent->x = t->x & ~h;
  parent->y = t->y & ~h;
#ifdef T8_DTRI_TO_DTET
  parent->z = t->z & ~h;
#endif
}

void
t8_dtri_compute_coords (const t8_dtri_t * t,
                        t8_dtri_coord_t coordinates[T8_DTRI_DIM],
                        const int vertex)
{
  t8_dtri_type_t      type;
  int                 ei;
#ifdef T8_DTRI_TO_DTET
  int                 ej;
#endif
  t8_dtri_coord_t     h;

  T8_ASSERT (0 <= vertex && vertex < T8_DTRI_FACES);

  type = t->type;
  h = T8_DTRI_LEN (t->level);
#ifndef T8_DTRI_TO_DTET
  ei = type;
#else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

  coordinates[0] = t->x;
  coordinates[1] = t->y;
#ifdef T8_DTRI_TO_DTET
  coordinates[2] = t->z;
#endif
  if (vertex == 0) {
    return;
  }
  coordinates[ei] += h;
#ifndef T8_DTRI_TO_DTET
  if (vertex == 2) {
    coordinates[1 - ei] += h;
    return;
  }
#else
  if (vertex == 2) {
    coordinates[ej] += h;
    return;
  }
  if (vertex == 3) {
    coordinates[(ei + 1) % 3] += h;
    coordinates[(ei + 2) % 3] += h;
  }
  /* done 3D */
#endif
}

void
t8_dtri_compute_all_coords (const t8_dtri_t * t,
                            t8_dtri_coord_t
                            coordinates[T8_DTRI_FACES][T8_DTRI_DIM])
{
  t8_dtri_type_t      type;
  int                 ei;
#ifdef T8_DTRI_TO_DTET
  int                 ej;
#endif
  int                 i;
  t8_dtri_coord_t     h;

  type = t->type;
  h = T8_DTRI_LEN (t->level);
#ifndef T8_DTRI_TO_DTET
  ei = type;
#else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

  coordinates[0][0] = t->x;
  coordinates[0][1] = t->y;
#ifdef T8_DTRI_TO_DTET
  coordinates[0][2] = t->z;
#endif
  for (i = 0; i < T8_DTRI_DIM; i++) {
    coordinates[1][i] = coordinates[0][i];
#ifndef T8_DTRI_TO_DTET
    coordinates[2][i] = coordinates[0][i] + h;
#else
    coordinates[2][i] = coordinates[0][i];
    coordinates[3][i] = coordinates[0][i] + h;
#endif
  }
  coordinates[1][ei] += h;
#ifdef T8_DTRI_TO_DTET
  coordinates[2][ei] += h;
  coordinates[2][ej] += h;
#endif
}

/* The childid here is the Morton child id
 * (TODO: define this)
 * It is possible that the function is called with
 * elem = child */
void
t8_dtri_child (const t8_dtri_t * elem, int childid, t8_dtri_t * child)
{
  const t8_dtri_t    *t = (const t8_dtri_t *) elem;
  t8_dtri_t          *c = (t8_dtri_t *) child;
  t8_dtri_coord_t     t_coordinates[T8_DTRI_DIM];
  int                 vertex;
  int                 Bey_cid;

  T8_ASSERT (t->level < T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= childid && childid < T8_DTRI_CHILDREN);

  Bey_cid = t8_dtri_index_to_bey_number[elem->type][childid];

  /* Compute anchor coordinates of child */
  if (Bey_cid == 0) {
    /* TODO: would it be better do drop this if and
     *       capture it with (t->x+t->x)>>1 below? */
    c->x = t->x;
    c->y = t->y;
#ifdef T8_DTRI_TO_DTET
    c->z = t->z;
#endif
  }
  else {
    vertex = t8_dtri_beyid_to_vertex[Bey_cid];
    /* i-th anchor coordinate of child is (X_(0,i)+X_(vertex,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    t8_dtri_compute_coords (t, t_coordinates, vertex);
    c->x = (t->x + t_coordinates[0]) >> 1;
    c->y = (t->y + t_coordinates[1]) >> 1;
#ifdef T8_DTRI_TO_DTET
    c->z = (t->z + t_coordinates[2]) >> 1;
#endif
  }

  /* Compute type of child */
  c->type = t8_dtri_type_of_child[t->type][Bey_cid];

  c->level = t->level + 1;
}

void
t8_dtri_childrenpv (const t8_dtri_t * t, t8_dtri_t * c[T8_DTRI_CHILDREN])
{
  t8_dtri_coord_t     t_coordinates[T8_DTRI_FACES][T8_DTRI_DIM];
  const int8_t        level = t->level + 1;
  int                 i;
  int                 Bey_cid;
  int                 vertex;

  T8_ASSERT (t->level < T8_DTRI_MAXLEVEL);
  t8_dtri_compute_all_coords (t, t_coordinates);
  c[0]->x = t->x;
  c[0]->y = t->y;
#ifdef T8_DTRI_TO_DTET
  c[0]->z = t->z;
#endif
  c[0]->type = t->type;
  c[0]->level = level;
  for (i = 1; i < T8_DTRI_CHILDREN; i++) {
    Bey_cid = t8_dtri_index_to_bey_number[t->type][i];
    vertex = t8_dtri_beyid_to_vertex[Bey_cid];
    /* i-th anchor coordinate of child is (X_(0,i)+X_(vertex,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    c[i]->x = (t->x + t_coordinates[vertex][0]) >> 1;
    c[i]->y = (t->y + t_coordinates[vertex][1]) >> 1;
#ifdef T8_DTRI_TO_DTET
    c[i]->z = (t->z + t_coordinates[vertex][2]) >> 1;
#endif
    c[i]->type = t8_dtri_type_of_child[t->type][Bey_cid];
    c[i]->level = level;
  }
}

/* The sibid here is the Morton child id of the parent.
 * TODO: Implement this algorithm directly w/o using
 * parent and child
 * TODO: CB agrees, make this as non-redundant as possible */
void
t8_dtri_sibling (const t8_dtri_t * elem, int sibid, t8_dtri_t * sibling)
{
  T8_ASSERT (0 <= sibid && sibid < T8_DTRI_CHILDREN);
  T8_ASSERT (((const t8_dtri_t *) elem)->level > 0);
  t8_dtri_parent (elem, sibling);
  t8_dtri_child (sibling, sibid, sibling);
}

/* Saves the neighbour of T along face "face" in N
 * returns the facenumber of N along which T is its neighbour */
int
t8_dtri_face_neighbour (const t8_dtri_t * t, t8_dtri_t * n, int face)
{
  /* TODO: document what happens if outside of root tet */
  int                 type_new, type_old;
#ifdef T8_DTRI_TO_DTET
  int                 sign;
#endif
  int                 ret = -1;
  int8_t              level;
  t8_dtri_coord_t     coords[3];

  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);

  level = t->level;
  type_old = t->type;
  type_new = type_old;
  coords[0] = t->x;
  coords[1] = t->y;
#ifdef T8_DTRI_TO_DTET
  coords[2] = t->z;
#endif

#ifndef T8_DTRI_TO_DTET
  /* 2D */
  ret = 2 - face;
  type_new = 1 - type_old;
  if (face == 0) {
    coords[type_old] += T8_DTRI_LEN (level);
  }
  else if (face == 2) {
    coords[1 - type_old] -= T8_DTRI_LEN (level);
  }
  /* 2D end */
#else
  /* 3D */
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
      coords[type_old / 2] += T8_DTRI_LEN (level);
      type_new += (type_new % 2 == 0 ? 4 : 2);
    }
    else {                      /* face == 3 */

      /* type: 1,2 --> z-1
       *       3,4 --> x-1
       *       5,0 --> y-1 */
      coords[((type_new + 3) % 6) / 2] -= T8_DTRI_LEN (level);
      type_new += (type_new % 2 == 0 ? 2 : 4);
    }
    type_new %= 6;
    ret = 3 - face;
  }
  /* 3D end */
#endif
  n->x = coords[0];
  n->y = coords[1];
#ifdef T8_DTRI_TO_DTET
  n->z = coords[2];
#endif
  n->level = level;
  n->type = type_new;
  return ret;
}

int
t8_dtri_is_inside_root (t8_dtri_t * t)
{
  int                 is_inside;
  is_inside = (t->x >= 0 && t->x < T8_DTRI_ROOT_LEN) && (t->y >= 0) &&
#ifdef T8_DTRI_TO_DTET
    (t->z >= 0) &&
#endif
#ifndef T8_DTRI_TO_DTET
    (t->y - t->x <= 0) && (t->y == t->x ? t->type == 0 : 1) &&
#else
    (t->z - t->x <= 0) &&
    (t->y - t->z <= 0) &&
    (t->z == t->x ? (3 <= t->type && 5 <= t->type) : 1) &&
    (t->y == t->x ? (1 <= t->type && 3 <= t->type) : 1) &&
#endif
    1;
  return is_inside;
}

int
t8_dtri_is_outside (const t8_dtri_t * t, int8_t roottype, int8_t level)
{
  t8_dtri_coord_t     n1, n2, cubex, cubey;
  t8_dtri_coord_t     bitmask;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t     dir3, cubez;
  int                 sign;
#endif

  /* Compute anchor coordinates of the ancestor cube of t */
  bitmask = ~(T8_DTRI_LEN (level) - 1);
  cubex = t->x & bitmask;
  cubey = t->y & bitmask;
#ifdef T8_DTRI_TO_DTET
  cubez = t->z & bitmask;
#endif
#ifndef T8_DTRI_TO_DTET
  /* 2D */
  n1 = (roottype == 0) ? t->x - cubex : t->y - cubey;
  n2 = (roottype == 0) ? t->y - cubey : t->x - cubex;

  return n1 >= T8_DTRI_LEN (level) || n2 < 0 || n2 - n1 > 0
    || (n2 == n1 && t->type == 1 - roottype);
#else
  /* 3D */
  /* *INDENT-OFF* */
  n1 = roottype / 2 == 0 ? t->x - cubex :
       roottype / 2 == 2 ? t->z - cubez : t->y - cubey;
  n2 = (roottype + 3) % 6 == 0 ? t->x - cubex :
       (roottype + 3) % 6 == 2 ? t->z - cubez : t->y - cubey;
  dir3 = roottype % 3 == 2 ? t->x - cubex :
         roottype % 3 == 0 ? t->z - cubez : t->y - cubey;
  sign = (roottype % 2 == 0) ? 1 : -1;
  /* *INDENT-ON* */

  roottype += 6;                /* We need to compute modulo six and want
                                   to avoid negative numbers when substracting from roottype. */
  return n1 >= T8_DTRI_LEN (level) || n2 < 0 || dir3 - n1 > 0 || n2 - dir3 > 0
    || (dir3 == n1 && (t->type == (roottype + sign * 1) % 6
                       || t->type == (roottype + sign * 2) % 6
                       || t->type == (roottype + sign * 3) % 6))
    || (dir3 == n2 && (t->type == (roottype - sign * 1) % 6
                       || t->type == (roottype - sign * 2) % 6
                       || t->type == (roottype - sign * 3) % 6));
#endif
}

int
t8_dtri_is_equal (const t8_dtri_t * t1, const t8_dtri_t * t2)
{
  return (t1->level == t1->level && t1->type == t2->type &&
          t1->x == t1->x && t1->y == t1->y
#ifdef T8_DTRI_TO_DTET
          && t1->z == t1->z
#endif
    );
}

/* we check if t1 and t2 lie in the same subcube and have
 * the same level and parent type */
int
t8_dtri_is_sibling (const t8_dtri_t * t1, const t8_dtri_t * t2)
{
  t8_dtri_coord_t     exclorx, exclory;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t     exclorz;
#endif

  t8_dtri_cube_id_t   cid1, cid2;

  /* TODO: zulassen, dass level 0 element sein eigener sibling ist (?)
   *       JH says no, b/c adopting the convention from p4est a tetrahedron
   *       is never its own sibling */
  if (t1->level == 0) {
    return 0;
  }

  exclorx = t1->x ^ t2->x;
  exclory = t1->y ^ t2->y;
#ifdef T8_DTRI_TO_DTET
  exclorz = t1->z ^ t2->z;
#endif
  cid1 = compute_cubeid (t1, t1->level);
  cid2 = compute_cubeid (t2, t2->level);

  return
    (t1->level == t2->level) &&
    ((exclorx & ~T8_DTRI_LEN (t1->level)) == 0) &&
    ((exclory & ~T8_DTRI_LEN (t1->level)) == 0) &&
#ifdef T8_DTRI_TO_DTET
    ((exclorz & ~T8_DTRI_LEN (t1->level)) == 0) &&
#endif
    t8_dtri_cid_type_to_parenttype[cid1][t1->type] ==
    t8_dtri_cid_type_to_parenttype[cid2][t2->type] && t1->type != t2->type;
}

int
t8_dtri_is_parent (const t8_dtri_t * t, const t8_dtri_t * c)
{
  t8_dtri_cube_id_t   cid;

  cid = compute_cubeid (c, c->level);
  return
    (t->level + 1 == c->level) &&
    (t->x == (c->x & ~T8_DTRI_LEN (c->level))) &&
    (t->y == (c->y & ~T8_DTRI_LEN (c->level))) &&
#ifdef T8_DTRI_TO_DTET
    (t->z == (c->z & ~T8_DTRI_LEN (c->level))) &&
#endif
    t->type == t8_dtri_cid_type_to_parenttype[cid][c->type] && 1;
}

int
t8_dtri_is_ancestor (const t8_dtri_t * t, const t8_dtri_t * c)
{
  t8_dtri_coord_t     n1, n2;
  t8_dtri_coord_t     exclorx;
  t8_dtri_coord_t     exclory;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t     exclorz;
  t8_dtri_coord_t     dir3;
  int                 sign;
#endif
  int8_t              type_t;

  if (t->level >= c->level) {
    return 0;
  }

  exclorx = (t->x ^ c->x) >> (T8_DTRI_MAXLEVEL - t->level);
  exclory = (t->y ^ c->y) >> (T8_DTRI_MAXLEVEL - t->level);
#ifdef T8_DTRI_TO_DTET
  exclorz = (t->z ^ c->z) >> (T8_DTRI_MAXLEVEL - t->level);
#endif

  /* TODO: implement */
  SC_ABORT ("Not implemented");

  if (exclorx == 0 && exclory == 0
#ifdef T8_DTRI_TO_DTET
      && exclorz == 0
#endif
    ) {
    /* t and c have the same cube as ancestor.
     * Now check if t has the correct type to be c's ancestor. */
    type_t = t->type;
#ifndef T8_DTRI_TO_DTET
    /* 2D */
    n1 = (type_t == 0) ? c->x - t->x : c->y - t->y;
    n2 = (type_t == 0) ? c->y - t->y : c->x - t->x;

    return !(n1 >= T8_DTRI_LEN (t->level) || n2 < 0 || n2 - n1 > 0
             || (n2 == n1 && t->type == 1 - type_t));
#else
    /* 3D */
      /* *INDENT-OFF* */
      n1 = type_t / 2 == 0 ? c->x - t->x :
           type_t / 2 == 2 ? c->z - t->z : c->y - t->y;
      n2 = (type_t + 3) % 6 == 0 ? c->x - t->x :
           (type_t + 3) % 6 == 2 ? c->z - t->z : c->y - t->y;
      dir3 = type_t % 3 == 2 ? c->x - t->x :
             type_t % 3 == 0 ? c->z - t->z : c->y - t->y;
      sign = (type_t % 2 == 0) ? 1 : -1;
      /* *INDENT-ON* */

    type_t += 6;                /* We need to compute modulo six and want
                                   to avoid negative numbers when substracting from type_t. */
    return !(n1 >= T8_DTRI_LEN (t->level) || n2 < 0 || dir3 - n1 > 0
             || n2 - dir3 > 0 || (dir3 == n1
                                  && (t->type == (type_t + sign * 1) % 6
                                      || t->type == (type_t + sign * 2) % 6
                                      || t->type == (type_t + sign * 3) % 6))
             || (dir3 == n2
                 && (t->type == (type_t - sign * 1) % 6
                     || t->type == (type_t - sign * 2) % 6
                     || t->type == (type_t - sign * 3) % 6)));
#endif
  }
  else {
    return 0;
  }
}
