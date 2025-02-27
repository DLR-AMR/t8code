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

#ifndef T8_DTRI_TO_DTET
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#else
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
#endif

typedef int8_t t8_dtri_cube_id_t;

/* Compute the cube-id of t's ancestor of level "level" in constant time.
 * If "level" is greater then t->level then the cube-id 0 is returned. */
static t8_dtri_cube_id_t
compute_cubeid (const t8_dtri_t *t, int level)
{
  t8_dtri_cube_id_t id = 0;
  t8_dtri_coord_t h;

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

/* A routine to compute the type of t's ancestor of level "level", if its type at an intermediate level is already 
 * known. If "level" equals t's level then t's type is returned. It is not allowed to call this function with "level" 
 * greater than t->level. This method runs in O(t->level - level).
 */
static t8_dtri_type_t
compute_type_ext (const t8_dtri_t *t, int level, t8_dtri_type_t known_type, int known_level)
{
  int8_t type = known_type;
  t8_dtri_cube_id_t cid;
  int i;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (known_level <= t->level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /* TODO: the type of the root tet is hardcoded to 0
     *       maybe once we want to allow the root tet to have different types */
    return 0;
  }
  for (i = known_level; i > level; i--) {
    cid = compute_cubeid (t, i);
    /* compute type as the type of T^{i+1}, that is T's ancestor of level i+1 */
    type = t8_dtri_cid_type_to_parenttype[cid][type];
  }
  return type;
}

/* A routine to compute the type of t's ancestor of level "level". If "level" equals t's level then t's type is 
 * returned. It is not allowed to call this function with "level" greater than t->level. This method runs in 
 * O(t->level - level).
 */
static t8_dtri_type_t
compute_type (const t8_dtri_t *t, int level)
{
  return compute_type_ext (t, level, t->type, t->level);
}

void
t8_dtri_copy (const t8_dtri_t *t, t8_dtri_t *dest)
{
  if (t == dest) {
    /* Do nothing if they are already the same. */
    return;
  }
  memcpy (dest, t, sizeof (t8_dtri_t));
}

int
t8_dtri_equal (const t8_dtri_t *elem1, const t8_dtri_t *elem2)
{
  return (elem1->level == elem2->level && elem1->type == elem2->type && elem1->x == elem2->x && elem1->y == elem2->y
#ifdef T8_DTRI_TO_DTET
          && elem1->z == elem2->z
#endif
  );
}

int
t8_dtri_compare (const t8_dtri_t *t1, const t8_dtri_t *t2)
{
  int maxlvl;
  t8_linearidx_t id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t1->level, t2->level);
  /* Compute the linear ids of the elements */
  id1 = t8_dtri_linear_id (t1, maxlvl);
  id2 = t8_dtri_linear_id (t2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the triangle with the smaller level is considered smaller */
    T8_ASSERT (t1->level != t2->level || t8_dtri_is_equal (t1, t2));
    return t1->level - t2->level;
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_dtri_parent (const t8_dtri_t *t, t8_dtri_t *parent)
{
  t8_dtri_cube_id_t cid;
  t8_dtri_coord_t h;

  T8_ASSERT (t->level > 0);

  h = T8_DTRI_LEN (t->level);
  /* Compute type of parent */
  cid = compute_cubeid (t, t->level);
  parent->type = t8_dtri_cid_type_to_parenttype[cid][t->type];
  /* Set coordinates of parent */
  parent->x = t->x & ~h;
  parent->y = t->y & ~h;
#ifdef T8_DTRI_TO_DTET
  parent->z = t->z & ~h;
#endif
  parent->level = t->level - 1;
}

void
t8_dtri_ancestor (const t8_dtri_t *t, int level, t8_dtri_t *ancestor)
{
  /* TODO: find out, at which level difference it is faster to use    *
   * the arithmetic computation of ancestor type
   * opposed to iteratively computing the parent type.
   */
  t8_dtri_coord_t delta_x, delta_y, diff_xy;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t delta_z, diff_xz, diff_yz;
  t8_dtet_type_t possible_types[6] = { 1, 1, 1, 1, 1, 1 };
  int i;
#ifdef T8_ENABLE_DEBUG
  int set_type = 0;
#endif
#endif /* T8_DTRI_TO_DTET */

  /* delta_{x,y} = t->{x,y} - ancestor->{x,y}
   * the difference of the coordinates.
   * Needed to compute the type of the ancestor. */
  delta_x = t->x & (T8_DTRI_LEN (level) - 1);
  delta_y = t->y & (T8_DTRI_LEN (level) - 1);
#ifdef T8_DTRI_TO_DTET
  delta_z = t->z & (T8_DTRI_LEN (level) - 1);
#endif

  /* The coordinates of the ancestor. It is necessary
   * to compute the delta first, since ancestor and t
   * could point to the same triangle. */
  ancestor->x = t->x & ~(T8_DTRI_LEN (level) - 1);
  ancestor->y = t->y & ~(T8_DTRI_LEN (level) - 1);
#ifdef T8_DTRI_TO_DTET
  ancestor->z = t->z & ~(T8_DTRI_LEN (level) - 1);
#endif

#ifndef T8_DTRI_TO_DTET
  /* The type of the ancestor depends on delta_x - delta_y */
  diff_xy = delta_x - delta_y;
  if (diff_xy > 0) {
    ancestor->type = 0;
  }
  else if (diff_xy < 0) {
    ancestor->type = 1;
  }
  else {
    T8_ASSERT (diff_xy == 0);
    ancestor->type = t->type;
  }

#else
  /* The sign of each diff reduces the number of possible types
 * for the ancestor. At the end only one possible type is left,
 * this type's entry in the possible_types array will be positive.
 */

  diff_xy = delta_x - delta_y;
  diff_xz = delta_x - delta_z;
  diff_yz = delta_y - delta_z;

  /* delta_x - delta_y */
  if (diff_xy > 0) {
    possible_types[2] = possible_types[3] = possible_types[4] = 0;
  }
  else if (diff_xy < 0) {
    possible_types[0] = possible_types[1] = possible_types[5] = 0;
  }
  else {
    T8_ASSERT (diff_xy == 0);
    if (t->type == 0 || t->type == 1 || t->type == 5) {
      possible_types[2] = possible_types[3] = possible_types[4] = 0;
    }
    else {
      possible_types[0] = possible_types[1] = possible_types[5] = 0;
    }
  }

  /* delta_x - delta_z */
  if (diff_xz > 0) {
    possible_types[3] = possible_types[4] = possible_types[5] = 0;
  }
  else if (diff_xz < 0) {
    possible_types[0] = possible_types[1] = possible_types[2] = 0;
  }
  else {
    T8_ASSERT (diff_xz == 0);
    if (t->type == 0 || t->type == 1 || t->type == 2) {
      possible_types[3] = possible_types[4] = possible_types[5] = 0;
    }
    else {
      possible_types[0] = possible_types[1] = possible_types[2] = 0;
    }
  }

  /* delta_y - delta_z */
  if (diff_yz > 0) {
    possible_types[0] = possible_types[4] = possible_types[5] = 0;
  }
  else if (diff_yz < 0) {
    possible_types[1] = possible_types[2] = possible_types[3] = 0;
  }
  else {
    T8_ASSERT (diff_yz == 0);
    if (t->type == 1 || t->type == 2 || t->type == 3) {
      possible_types[0] = possible_types[4] = possible_types[5] = 0;
    }
    else {
      possible_types[1] = possible_types[2] = possible_types[3] = 0;
    }
  }

  /* Got through possible_types array and find the only entry
   * that is nonzero
   */
  for (i = 0; i < 6; i++) {
    T8_ASSERT (possible_types[i] == 0 || possible_types[i] == 1);
    if (possible_types[i] == 1) {
      ancestor->type = i;
#ifdef T8_ENABLE_DEBUG
      T8_ASSERT (set_type != 1);
      set_type = 1;
#endif
    }
  }
  T8_ASSERT (set_type == 1);
#endif /* T8_DTRI_TO_DTET */
  ancestor->level = level;
}

/* Compute the coordinates of a given vertex of a triangle/tet */
void
t8_dtri_compute_integer_coords (const t8_dtri_t *elem, const int vertex, t8_dtri_coord_t coordinates[T8_DTRI_DIM])
{
  /* Calculate the vertex coordinates of a triangle/tetrahedron in relation to its orientation. Orientations are 
   * described here: https://doi.org/10.1137/15M1040049
   * 1---------------------2
   * |   orientation     /  2
   * |       1         /  / |
   * |               /  /   |
   * |             /  /     |
   * |           /  /       |
   * |         /  /         |
   * |       /  /           |
   * |     /  /             |
   * |   /  /  orientation  |
   * | /  /        0        |
   * 0  /                   |
   *   0--------------------1
   *
   *   y
   *   ^
   *   |
   *   z--> x
   */
  t8_dtri_type_t type;
  int ei;
#ifdef T8_DTRI_TO_DTET
  int ej;
#endif
  t8_dtri_coord_t h;
  T8_ASSERT (0 <= vertex && vertex < T8_DTRI_FACES);

  type = elem->type;
  h = T8_DTRI_LEN (elem->level);
#ifndef T8_DTRI_TO_DTET
  ei = type;
#else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

  coordinates[0] = elem->x;
  coordinates[1] = elem->y;
#ifdef T8_DTRI_TO_DTET
  coordinates[2] = elem->z;
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
t8_dtri_compute_vertex_ref_coords (const t8_dtri_t *elem, const int vertex, double coordinates[T8_DTRI_DIM])
{
  int coords_int[T8_DTRI_DIM];
  T8_ASSERT (0 <= vertex && vertex < T8_DTRI_CORNERS);

  t8_dtri_compute_integer_coords (elem, vertex, coords_int);
  /* Since the integer coordinates are coordinates w.r.t to
   * the embedding into [0,T8_DTRI_ROOT_LEN]^d, we just need
   * to divide them by the root length. */
  coordinates[0] = coords_int[0] / (double) T8_DTRI_ROOT_LEN;
  coordinates[1] = coords_int[1] / (double) T8_DTRI_ROOT_LEN;
#ifdef T8_DTRI_TO_DTET
  coordinates[2] = coords_int[2] / (double) T8_DTRI_ROOT_LEN;
#endif
}

void
t8_dtri_compute_reference_coords (const t8_dtri_t *elem, const double *ref_coords, const size_t num_coords,
#ifndef T8_DTRI_TO_DTET
                                  const size_t skip_coords,
#endif
                                  double *out_coords)
{
  /* Calculate the reference coordinates of a triangle/tetrahedron in
   * relation to its orientation. Orientations are described here:
   * https://doi.org/10.1137/15M1040049
   * 1---------------------2
   * |   orientation     /  2
   * |       1         /  / |
   * |               /  /   |
   * |             /  /     |
   * |           /  /       |
   * |         /  /         |
   * |       /  /           |
   * |     /  /             |
   * |   /  /  orientation  |
   * | /  /        0        |
   * 0  /                   |
   *   0--------------------1
   *
   *   y
   *   ^
   *   |
   *   z--> x
   */
  T8_ASSERT (ref_coords != NULL);

  t8_dtri_type_t type;
  t8_dtri_coord_t h;

  type = elem->type;
  h = T8_DTRI_LEN (elem->level);
#ifndef T8_DTRI_TO_DTET
  const int tri_orientation = type;
#else
  /* These integers define the sequence, in which the ref_coords are added
   * to the out_coords */
  const int tet_orientation0 = type / 2;
  const int tet_orientation1 = (tet_orientation0 + ((type % 2 == 0) ? 1 : 2)) % 3;
  const int tet_orientation2 = (tet_orientation0 + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif
  for (size_t coord = 0; coord < num_coords; ++coord) {
    /* offset defines, how many coordinates to skip in an iteration. */
#ifndef T8_DTRI_TO_DTET
    const size_t offset = (2 + skip_coords) * coord;
    const size_t offset_3d = 3 * coord;
#else
    const size_t offset = 3 * coord;
#endif
    out_coords[offset + 0] = elem->x;
    out_coords[offset + 1] = elem->y;
#ifdef T8_DTRI_TO_DTET
    out_coords[offset + 2] = elem->z;
#endif
#ifndef T8_DTRI_TO_DTET
    out_coords[offset + tri_orientation] += h * ref_coords[offset_3d + 0];
    out_coords[offset + 1 - tri_orientation] += h * ref_coords[offset_3d + 1];
#else
    out_coords[offset + tet_orientation0] += h * ref_coords[offset + 0];
    out_coords[offset + tet_orientation1] += h * ref_coords[offset + 1];
    out_coords[offset + tet_orientation2] += h * ref_coords[offset + 2];

    /* done 3D */
#endif
    /* Since the integer coordinates are coordinates w.r.t to
     * the embedding into [0,T8_DTRI_ROOT_LEN]^d, we just need
     * to divide them by the root length. */
    out_coords[offset + 0] /= (double) T8_DTRI_ROOT_LEN;
    out_coords[offset + 1] /= (double) T8_DTRI_ROOT_LEN;
#ifdef T8_DTRI_TO_DTET
    out_coords[offset + 2] /= (double) T8_DTRI_ROOT_LEN;
#endif
  }
}

/* Compute the coordinates of each vertex of a triangle/tet */
void
t8_dtri_compute_all_coords (const t8_dtri_t *elem, t8_dtri_coord_t coordinates[T8_DTRI_FACES][T8_DTRI_DIM])
{
  /* Calculate the vertex coordinates of a triangle/tetrahedron in
   * relation to its orientation. Orientations are described here:
   * https://doi.org/10.1137/15M1040049
   * 1---------------------2
   * |   orientation     /  2
   * |       1         /  / |
   * |               /  /   |
   * |             /  /     |
   * |           /  /       |
   * |         /  /         |
   * |       /  /           |
   * |     /  /             |
   * |   /  /  orientation  |
   * | /  /        0        |
   * 0  /                   |
   *   0--------------------1
   *
   *   y
   *   ^
   *   |
   *   z--> x
   */
  t8_dtri_type_t type;
  int ei;
#ifdef T8_DTRI_TO_DTET
  int ej;
#endif
  int i;
  t8_dtri_coord_t h;

  type = elem->type;
  h = T8_DTRI_LEN (elem->level);
#ifndef T8_DTRI_TO_DTET
  ei = type;
#else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

  coordinates[0][0] = elem->x;
  coordinates[0][1] = elem->y;
#ifdef T8_DTRI_TO_DTET
  coordinates[0][2] = elem->z;
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
#ifdef T8_ENABLE_DEBUG
  /* We check whether the results are the same as with the
   * t8_dtri_compute_integer_coords function.
   */
  {
    int ivertex;
    t8_dtri_coord_t coords[T8_DTRI_DIM];
    for (ivertex = 0; ivertex < T8_DTRI_FACES; ivertex++) {
      t8_dtri_compute_integer_coords (elem, ivertex, coords);
      T8_ASSERT (coords[0] == coordinates[ivertex][0]);
      T8_ASSERT (coords[1] == coordinates[ivertex][1]);
#ifdef T8_DTRI_TO_DTET
      T8_ASSERT (coords[2] == coordinates[ivertex][2]);
#endif
    }
  }
#endif
}

/* The childid here is the Morton child id
 * (TODO: define this)
 * It is possible that the function is called with
 * elem = child */
void
t8_dtri_child (const t8_dtri_t *t, int childid, t8_dtri_t *child)
{
  t8_dtri_t *c = (t8_dtri_t *) child;
  t8_dtri_coord_t t_coordinates[T8_DTRI_DIM];
  int vertex;
  int Bey_cid;

  T8_ASSERT (t->level < T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= childid && childid < T8_DTRI_CHILDREN);

  Bey_cid = t8_dtri_index_to_bey_number[t->type][childid];

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
    t8_dtri_compute_integer_coords (t, vertex, t_coordinates);
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
t8_dtri_childrenpv (const t8_dtri_t *t, t8_dtri_t *c[T8_DTRI_CHILDREN])
{
  t8_dtri_coord_t t_coordinates[T8_DTRI_FACES][T8_DTRI_DIM];
  const int8_t level = t->level + 1;
  int i;
  int Bey_cid;
  int vertex;
  t8_dtri_type_t t_type = t->type;

  T8_ASSERT (t->level < T8_DTRI_MAXLEVEL);
  t8_dtri_compute_all_coords (t, t_coordinates);
  /* We use t_coordinates[0] for t->x,y,z to ensure that the function is valid
   * if called with t = c[0]. If we would use t->x later it would be the newly
   * computed value c[0]->x. */
  c[0]->x = t_coordinates[0][0]; /* t->x */
  c[0]->y = t_coordinates[0][1]; /* t->y */
#ifdef T8_DTRI_TO_DTET
  c[0]->z = t_coordinates[0][2]; /* t->z */
#endif
  c[0]->type = t_type;
  c[0]->level = level;
  for (i = 1; i < T8_DTRI_CHILDREN; i++) {
    Bey_cid = t8_dtri_index_to_bey_number[t_type][i];
    vertex = t8_dtri_beyid_to_vertex[Bey_cid];
    /* i-th anchor coordinate of child is (X_(0,i)+X_(vertex,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    c[i]->x = (t_coordinates[0][0] + t_coordinates[vertex][0]) >> 1;
    c[i]->y = (t_coordinates[0][1] + t_coordinates[vertex][1]) >> 1;
#ifdef T8_DTRI_TO_DTET
    c[i]->z = (t_coordinates[0][2] + t_coordinates[vertex][2]) >> 1;
#endif
    c[i]->type = t8_dtri_type_of_child[t_type][Bey_cid];
    c[i]->level = level;
#ifdef T8_ENABLE_DEBUG
    {
      /* We check whether the child computed here equals to the child
       * computed in the t8_dtri_child function. */
      t8_dtri_t check_child;
      /* We use c[0] here instead of t, since we explicitly allow t=c[0] as input
       * and thus the values of t may be already overwritten. However the only
       * difference from c[0] to t is in the level. */
      c[0]->level--;
      t8_dtri_child (c[0], i, &check_child);
      T8_ASSERT (check_child.x == c[i]->x && check_child.y == c[i]->y);
      T8_ASSERT (check_child.type == c[i]->type && check_child.level == c[i]->level);
#ifdef T8_DTRI_TO_DTET
      T8_ASSERT (check_child.z == c[i]->z);
#endif
      c[0]->level++;
    }
#endif
  }
}

#ifndef T8_DTRI_TO_DTET
int
t8_dtri_is_familypv (const t8_dtri_t *f[])
{
  const int8_t level = f[0]->level;
  if (level == 0 || level != f[1]->level || level != f[2]->level || level != f[3]->level) {
    return 0;
  }
  /* check whether the types are correct */
  const int type = f[0]->type;
  if (f[1]->type != 0 && f[2]->type != 1 && f[3]->type != type) {
    return 0;
  }
  /* check whether the coordinates are correct
   * triangles 1 and 2 have to have the same coordinates */
  if (f[1]->x != f[2]->x || f[1]->y != f[2]->y) {
    return 0;
  }
  const int dir1 = type;
  const t8_dtri_coord_t inc = T8_DTRI_LEN (level);
  t8_dtri_coord_t coords0[T8_DTRI_CHILDREN];
  t8_dtri_coord_t coords1[T8_DTRI_CHILDREN];
  coords0[0] = f[0]->x;
  coords0[1] = f[0]->y;
  coords1[0] = f[1]->x;
  coords1[1] = f[1]->y;
  return (coords1[dir1] == coords0[dir1] + inc && coords1[1 - dir1] == coords0[1 - dir1] && f[3]->x == f[0]->x + inc
          && f[3]->y == f[0]->y + inc);
}
#endif

/* The sibid here is the Morton child id of the parent.
 * TODO: Implement this algorithm directly w/o using
 * parent and child
 * TODO: CB agrees, make this as non-redundant as possible */
void
t8_dtri_sibling (const t8_dtri_t *elem, int sibid, t8_dtri_t *sibling)
{
  T8_ASSERT (0 <= sibid && sibid < T8_DTRI_CHILDREN);
  T8_ASSERT (((const t8_dtri_t *) elem)->level > 0);
  t8_dtri_parent (elem, sibling);
  t8_dtri_child (sibling, sibid, sibling);
}

/* Saves the neighbour of T along face "face" in N
 * returns the facenumber of N along which T is its neighbour */
int
t8_dtri_face_neighbour (const t8_dtri_t *t, int face, t8_dtri_t *n)
{
  /* TODO: document what happens if outside of root tet */
  int type_new, type_old;
#ifdef T8_DTRI_TO_DTET
  int sign;
#endif
  int ret = -1;
  int8_t level;
  t8_dtri_coord_t coords[3];

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
  type_new += 6; /* We want to compute modulo six and dont want negative numbers */
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
    else { /* face == 3 */

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

void
t8_dtri_nearest_common_ancestor (const t8_dtri_t *t1, const t8_dtri_t *t2, t8_dtri_t *r)
{
  int maxlevel, c_level, r_level;
  t8_dtri_type_t t1_type_at_l, t2_type_at_l;
  uint32_t exclorx, exclory;
#ifdef T8_DTRI_TO_DTET
  uint32_t exclorz;
#endif
  uint32_t maxclor;

  /* Find the level of the nca */
  exclorx = t1->x ^ t2->x;
  exclory = t1->y ^ t2->y;
#ifdef T8_DTRI_TO_DTET
  exclorz = t1->z ^ t2->z;

  maxclor = exclorx | exclory | exclorz;
#else
  maxclor = exclorx | exclory;
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  T8_ASSERT (maxlevel <= T8_DTRI_MAXLEVEL);

  c_level = (int8_t) SC_MIN (T8_DTRI_MAXLEVEL - maxlevel, (int) SC_MIN (t1->level, t2->level));
  r_level = c_level;
  /* c_level is the level of the nca cube surrounding t1 and t2.
   * If t1 and t2 have different types at this level, then we have to shrink this
   * level further. */
  t1_type_at_l = compute_type (t1, c_level);
  t2_type_at_l = compute_type (t2, c_level);
  while (t1_type_at_l != t2_type_at_l) {
    r_level--;
    t1_type_at_l = compute_type_ext (t1, r_level, t1_type_at_l, r_level + 1);
    t2_type_at_l = compute_type_ext (t2, r_level, t2_type_at_l, r_level + 1);
  }
  T8_ASSERT (r_level >= 0);
  /* Construct the ancestor of the first triangle at this level */
  t8_dtri_ancestor (t1, r_level, r);
}

void
t8_dtri_children_at_face (const t8_dtri_t *tri, int face, t8_dtri_t *children[], int num_children, int *child_indices)
{
  int child_ids_local[T8_DTRI_FACE_CHILDREN], i, *child_ids;

  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (num_children == T8_DTRI_FACE_CHILDREN);

  if (child_indices != NULL) {
    child_ids = child_indices;
  }
  else {
    child_ids = child_ids_local;
  }
#ifndef T8_DTRI_TO_DTET
  /* Triangle version */
  /* The first child is '0' for faces 1 and 2 and '1+type' for face 0 */
  child_ids[0] = face == 0 ? 1 + tri->type : 0;
  /* The second child is '1 + type' for face 2 and '3' for faces 0 and 1 */
  child_ids[1] = face == 2 ? 1 + tri->type : 3;
#else
  /* Tetrahedron version */
  for (i = 0; i < T8_DTET_FACE_CHILDREN; i++) {
    child_ids[i] = t8_dtet_face_child_id_by_type[tri->type][face][i];
  }
#endif

  /* Compute the children at the face.
   * We revert the order to compute children[0] last, since the usage
   * allows for tri == children[0].
   */
  for (i = T8_DTRI_FACE_CHILDREN - 1; i >= 0; i--) {
    t8_dtri_child (tri, child_ids[i], children[i]);
  }
}

int
t8_dtri_tree_face (__attribute__((unused)) t8_dtri_t *t, int face)
{
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  /* TODO: Assert if boundary */
#ifndef T8_DTRI_TO_DTET
  /* For triangles of type 0 the face number coincides with the number of the
   * root tree face. Triangles of type 1 cannot lie on the boundary of the
   * tree and thus the return value can be arbitrary. */
  return face;
#else
  /* For tets only tets of type not 3 can have tree boundary faces.
   * All these tets of type not 0 (types 1, 2, 4, and 5) can only have one of
   * their faces as boundary face. */
  switch (t->type) {
  case 0:
    return face;
    break;
  case 1:
    return 0;
    break; /* face 0 of the tet is the boundary face */
  case 2:
    return 1;
    break; /* face 2     "                "          */
  case 3:
    return -1;
    break; /* no face                                */
  case 4:
    return 2;
    break; /* face 1     "                "          */
  case 5:
    return 3;
    break; /* face 3     "                "          */
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
}

int
t8_dtri_root_face_to_face (t8_dtri_t *t, int root_face)
{
  T8_ASSERT (0 <= root_face && root_face < T8_DTRI_FACES);
#ifndef T8_DTRI_TO_DTET
  T8_ASSERT (t->type == 0);
  /* For triangles of type 0 the face number coincides with the number of the
   * root tree face. Triangles of type 1 cannot lie on the boundary of the
   * tree and thus the return value can be arbitrary. */
  return root_face;
#else
  /* For tets only tets of type not 3 can have tree boundary faces.
   * All these tets of type not 0 (types 1, 2, 4, and 5) can only have one of
   * their faces as boundary face. */
  T8_ASSERT (t->type != 3);
  switch (t->type) {
  case 0:
    return root_face;
    break;
  case 1:
    T8_ASSERT (root_face == 0);
    return 0;
    break; /* face 0 of the tet is the boundary face */
  case 2:
    T8_ASSERT (root_face == 1);
    return 2;
    break; /* face 2     "                "          */
  case 4:
    T8_ASSERT (root_face == 2);
    return 1;
    break; /* face 1     "                "          */
  case 5:
    T8_ASSERT (root_face == 3);
    return 3;
    break; /* face 3     "                "          */
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
}

int
t8_dtri_face_child_face (__attribute__((unused))const t8_dtri_t *triangle, int face, __attribute__((unused))int face_child)
{
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);
  T8_ASSERT (0 <= face_child && face_child < T8_DTRI_FACE_CHILDREN);
#ifndef T8_DTRI_TO_DTET
  /* For triangles the face number of the children is the same as the one
   * of the parent */
  return face;
#else
  {
    t8_dtet_type_t face_type;
    int is_middle_face = 0;
    switch (face) {
    case 0:
      return 0;
    case 1:
    case 2:
      /* If the given face_child is at the middle of the face (child 2 for type 0
       * triangles, child 1 for type 1 triangles), then the face number at the
       * tetrahedron is different. */
      face_type = t8_dtet_type_face_to_boundary[triangle->type][face][1];
      if ((face_type == 0 && face_child == 2) || (face_type == 1 && face_child == 1)) {
        is_middle_face = 1;
      }
      return is_middle_face ? face ^ 3 : face; /* face = 1 -> 2 : 1, face = 2 -> 1 : 2 */
    case 3:
      return 3;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }
#endif
}

int
t8_dtri_face_parent_face (const t8_dtri_t *triangle, int face)
{
  int parent_type, child_id, cid;
  T8_ASSERT (0 <= face && face < T8_DTRI_FACES);

  /* By convention, the tree root triangle lies on its (non-existent)
   * parent face. */
  if (triangle->level == 0) {
    return face;
  }
  /* For triangles, a triangle is only at the boundary of its parent if it
   * has the same type as the parent, in which case also the face numbers are the same */
  cid = compute_cubeid (triangle, triangle->level);
  parent_type = t8_dtri_cid_type_to_parenttype[cid][triangle->type];

#ifndef T8_DTRI_TO_DTET
  if (parent_type != triangle->type) {
    return -1;
  }
#endif

  /* The bey child id of the triangle,
   * a triangle of bey id i shares the faces i+1 and i+2 (%3) with its parent */
  child_id = t8_dtri_type_cid_to_beyid[triangle->type][cid];
  t8_dtri_child_id (triangle);
#ifndef T8_DTRI_TO_DTET
  if (face != child_id) {
    return face;
  }
  return -1;
#else
  /* For tets, all tets that have the same type as their parent do
   * share the all faces with number != child_id with their parent */
  if (parent_type == triangle->type && face != child_id) {
    return face;
  }
  /* For the tets whose type does not match their parent's type,
   * we can see from the following table which type/parent_type pairs
   * match with which face/parent_face pairs. */
  /*     p_type type  p_face face |p_type type  p_face face
   *       0      1     0     0   |  3      1     1     2
   *              2     1     2   |         2     0     0
   *              4     2     1   |         4     3     3
   *              5     3     3   |         5     2     1
   *       1      0     0     0   |  4      0     1     2
   *              2     3     3   |         2     2     1
   *              3     2     1   |         3     3     3
   *              5     1     2   |         5     0     0
   *       2      0     2     1   |  5      0     3     3
   *              1     3     3   |         1     2     1
   *              3     0     0   |         3     1     2
   *              4     1     2   |         4     0     0
   */
  /* We see, that parent type i matches all types except i and i+3 (mod 6)
   * We also see, that the faces always match in the following way: 0-0, 1-2, 2-1, 3-3
   * The face values are stored in the table t8_dtet_parent_type_to_face
   */
  if (t8_dtet_parent_type_type_to_face[parent_type][triangle->type] == face) {
    switch (face) {
    case 1:
    case 2:
      return face ^ 3; /* return 1 if face = 2, return 2 if face = 1 */
    default:
      return face; /* 0 if face = 0, 3 if face = 3 */
    }
  }
  return -1;
#endif
}

#ifndef T8_DTRI_TO_DTET
/* This function has only a triangle version. */
void
t8_dtri_transform_face (const t8_dtri_t *trianglein, t8_dtri_t *triangle2, int orientation, int sign,
                        int is_smaller_face)
{
  const t8_dtri_t *triangle1;
  t8_dtri_coord_t h = T8_DTRI_LEN (trianglein->level);
  t8_dtri_coord_t x = trianglein->x;

  T8_ASSERT (0 <= orientation && orientation <= 2);
  triangle2->level = trianglein->level;
  triangle2->type = trianglein->type;
  /*
   * The corners of the triangle are enumerated like this
   *        type 0                    type 1
   *      also root tree
   *         v_2                     v_1  v_2
   *         x                         x--x
   *        /|                         | /
   *       / |                         |/
   *      x--x                         x
   *    v_0  v_1                      v_0
   *
   */

  if (sign) {
    /* The tree faces have the same topological orientation, and
     * thus we have to perform a coordinate switch. */
    /* We use triangl2 as storage, since trianglein and triangle2 are allowed to
     * point to the same triangle */
    triangle1 = triangle2;
    t8_dtri_copy (trianglein, (t8_dtri_t *) triangle1);
    if (trianglein->type == 0) {
      ((t8_dtri_t *) triangle1)->y = trianglein->x - trianglein->y;
    }
    else {
      ((t8_dtri_t *) triangle1)->y = trianglein->x - trianglein->y - h;
    }
  }
  else {
    triangle1 = trianglein;
  }

  if (!is_smaller_face && orientation != 0 && !sign) {
    /* Translate orientation if triangle1 is not on the smaller face.
     *  sign = 0  sign = 1
     *  0 -> 0    0 -> 0
     *  1 -> 2    1 -> 1
     *  2 -> 1    2 -> 2
     */
    orientation = 3 - orientation;
  }
  x = triangle1->x; /* temporary store x coord in case triangle1 = triangle2 */
  switch (orientation) {
  case 0:
    t8_dtri_copy (triangle1, triangle2);
    break;
  case 1:
    triangle2->x = T8_DTRI_ROOT_LEN - h - triangle1->y;
    if (triangle1->type == 0) {
      triangle2->y = x - triangle1->y;
    }
    else {
      triangle2->y = x - triangle1->y - h;
    }
    break;
  case 2:
    if (triangle1->type == 0) {
      triangle2->x = T8_DTRI_ROOT_LEN - h + triangle1->y - x;
    }
    else {
      triangle2->x = T8_DTRI_ROOT_LEN + triangle1->y - x;
    }
    triangle2->y = T8_DTRI_ROOT_LEN - h - x;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
#endif

int
t8_dtri_is_inside_root (t8_dtri_t *t)
{
  int is_inside;

  if (t->level == 0) {
    /* A level 0 simplex is only inside the root simplex if it
     * is the root simplex. */
    return t->type == 0 && t->x == 0 && t->y == 0
#ifdef T8_DTRI_TO_DTET
           && t->z == 0
#endif
      ;
  }
  is_inside = (t->x >= 0 && t->x < T8_DTRI_ROOT_LEN) && (t->y >= 0) &&
#ifdef T8_DTRI_TO_DTET
              (t->z >= 0) &&
#endif
#ifndef T8_DTRI_TO_DTET
              (t->y - t->x <= 0) && (t->y == t->x ? t->type == 0 : 1) &&
#else
              (t->z - t->x <= 0) && (t->y - t->z <= 0) && (t->z == t->x ? (0 <= t->type && t->type < 3) : 1)
              &&                                                     /* types 0, 1, 2 */
              (t->y == t->z ? (0 == t->type || 4 <= t->type) : 1) && /* types 0, 4, 5 */
              /* If the anchor is on the x-y-z diagonal, only type 0 tets are inside root. */
              ((t->x == t->y && t->y == t->z) ? t->type == 0 : 1) &&
#endif
              1;
#ifdef T8_ENABLE_DEBUG
  /* Check if is_inside gives the same result as is_ancestor for the root element. */
  {
    t8_dtri_t root;
    t8_dtri_init_root (&root);
    T8_ASSERT (is_inside == t8_dtri_is_ancestor (&root, t));
  }
#endif
  return is_inside;
}

int
t8_dtri_is_root_boundary (const t8_dtri_t *t, int face)
{
  /* Dependent on the type and the face we have to check
   * different conditions */
#ifndef T8_DTRI_TO_DTET
  switch (t->type) {
  case 0:
    switch (face) {
    case 0:
      return t->x == T8_DTRI_ROOT_LEN - T8_DTRI_LEN (t->level);
    case 1:
      return t->x == t->y;
    case 2:
      return t->y == 0;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  case 1:
    return 0; /* type 1 triangles are never at the boundary */
  default:
    SC_ABORT_NOT_REACHED ();
  }
#else
  switch (t->type) {
  case 0:
    switch (face) {
    case 0:
      return t->x == T8_DTET_ROOT_LEN - T8_DTET_LEN (t->level);
    case 1:
      return t->x == t->z;
    case 2:
      return t->y == t->z;
    case 3:
      return t->y == 0;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  case 1:
    /* type 1 tets are only boundary at face 0 */
    return face == 0 && t->x == T8_DTET_ROOT_LEN - T8_DTET_LEN (t->level);
  case 2:
    /* type 1 tets are only boundary at face 2 */
    return face == 2 && t->x == t->z;
  case 3:
    /* type 3 tets are never at the root boundary */
    return 0;
  case 4:
    /* type 4 tets are only boundary at face 1 */
    return face == 1 && t->y == t->z;
  case 5:
    /* type 5 tets are only boundary at face 3 */
    return face == 3 && t->y == 0;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
  SC_ABORT_NOT_REACHED ();
  return 0; /* Prevent compiler warning */
}

int
t8_dtri_is_equal (const t8_dtri_t *t1, const t8_dtri_t *t2)
{
  return (t1->level == t2->level && t1->type == t2->type && t1->x == t2->x && t1->y == t2->y
#ifdef T8_DTRI_TO_DTET
          && t1->z == t2->z
#endif
  );
}

/* we check if t1 and t2 lie in the same subcube and have
 * the same level and parent type */
int
t8_dtri_is_sibling (const t8_dtri_t *t1, const t8_dtri_t *t2)
{
  t8_dtri_coord_t exclorx, exclory;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t exclorz;
#endif

  t8_dtri_cube_id_t cid1, cid2;

  if (t1->level == 0) {
    return t2->level == 0 && t1->x == t2->x && t1->y == t2->y &&
#ifdef T8_DTRI_TO_DTET
           t1->z == t2->z &&
#endif
           1;
  }

  exclorx = t1->x ^ t2->x;
  exclory = t1->y ^ t2->y;
#ifdef T8_DTRI_TO_DTET
  exclorz = t1->z ^ t2->z;
#endif
  cid1 = compute_cubeid (t1, t1->level);
  cid2 = compute_cubeid (t2, t2->level);

  return (t1->level == t2->level) && ((exclorx & ~T8_DTRI_LEN (t1->level)) == 0)
         && ((exclory & ~T8_DTRI_LEN (t1->level)) == 0) &&
#ifdef T8_DTRI_TO_DTET
         ((exclorz & ~T8_DTRI_LEN (t1->level)) == 0) &&
#endif
         t8_dtri_cid_type_to_parenttype[cid1][t1->type] == t8_dtri_cid_type_to_parenttype[cid2][t2->type];
}

int
t8_dtri_is_parent (const t8_dtri_t *t, const t8_dtri_t *c)
{
  t8_dtri_cube_id_t cid;

  cid = compute_cubeid (c, c->level);
  return (t->level + 1 == c->level) && (t->x == (c->x & ~T8_DTRI_LEN (c->level)))
         && (t->y == (c->y & ~T8_DTRI_LEN (c->level))) &&
#ifdef T8_DTRI_TO_DTET
         (t->z == (c->z & ~T8_DTRI_LEN (c->level))) &&
#endif
         t->type == t8_dtri_cid_type_to_parenttype[cid][c->type] && 1;
}

int
t8_dtri_is_ancestor (const t8_dtri_t *t, const t8_dtri_t *c)
{
  t8_dtri_coord_t n1, n2;
  t8_dtri_coord_t exclorx;
  t8_dtri_coord_t exclory;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t exclorz;
  t8_dtri_coord_t dir3;
  int sign;
#endif
  int8_t type_t;

  if (t->level > c->level) {
    return 0;
  }
  if (t->level == c->level) {
    return t8_dtri_is_equal (t, c);
  }

  exclorx = (t->x ^ c->x) >> (T8_DTRI_MAXLEVEL - t->level);
  exclory = (t->y ^ c->y) >> (T8_DTRI_MAXLEVEL - t->level);
#ifdef T8_DTRI_TO_DTET
  exclorz = (t->z ^ c->z) >> (T8_DTRI_MAXLEVEL - t->level);
#endif

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

    return !(n1 >= T8_DTRI_LEN (t->level) || n2 < 0 || n2 - n1 > 0 || (n2 == n1 && c->type == 1 - type_t));
#else
    /* 3D */
    /* Compute the coordinate directions according to this table:
     *    type(t) 0  1  2  3  4  5
     * n1         x  x  y  y  z  z
     * n2         y  z  z  x  x  y
     * dir3       z  y  x  z  y  x
     */
    /* clang-format off */
    n1 = type_t / 2 == 0 ? c->x - t->x :                      /* type(t) is 0 or 1 */
           type_t / 2 == 2 ? c->z - t->z : c->y - t->y;       /* type(t) is (4 or 5) or (2 or 3) */
    n2 = (type_t + 1) / 2 == 2 ? c->x - t->x :                /* type(t) is 3 or 4 */
           (type_t + 1) / 2 == 1 ? c->z - t->z : c->y - t->y; /* type(t) is (1 or 2) or (0 or 5) */
    dir3 = type_t % 3 == 2 ? c->x - t->x :                    /* type(t) is 2 or 5 */
             type_t % 3 == 0 ? c->z - t->z : c->y - t->y;     /* type(t) is (0 or 3) or (1 or 4) */
    sign = (type_t % 2 == 0) ? 1 : -1;
    /* clang-format on */

    type_t += 6; /* We need to compute modulo six and want
                                   to avoid negative numbers when subtracting from type_t. */

    return !(n1 >= T8_DTRI_LEN (t->level) || n2 < 0 || dir3 - n1 > 0 || n2 - dir3 > 0
             || (dir3 == n2
                 && (c->type == (type_t + sign * 1) % 6 || c->type == (type_t + sign * 2) % 6
                     || c->type == (type_t + sign * 3) % 6))
             || (dir3 == n1
                 && (c->type == (type_t - sign * 1) % 6 || c->type == (type_t - sign * 2) % 6
                     || c->type == (type_t - sign * 3) % 6))
             /* On the x-y-z diagonal only tets of the same type can be
              * ancestor of each other. */
             || (dir3 == n2 && n2 == n1 && type_t - 6 != c->type));
#endif
  }
  else {
    return 0;
  }
}

/* Compute the linear id of the first descendant of a triangle/tet */
static t8_linearidx_t
t8_dtri_linear_id_first_desc (const t8_dtri_t *t, int level)
{
  /* The id of the first descendant is the id of t in a uniform level
   * refinement */
  return t8_dtri_linear_id (t, level);
}

/* Compute the linear id of the last descendant of a triangle/tet */
static t8_linearidx_t
t8_dtri_linear_id_last_desc (const t8_dtri_t *t, int level)
{
  t8_linearidx_t id = 0, t_id;
  int exponent;

  T8_ASSERT (level >= t->level);
  /* The id of the last descendant at level consists of the id of t in
   * the first digits and then the local ids of all last children
   * (3 in 2d, 7 in 3d)
   */
  t_id = t8_dtri_linear_id (t, t->level);
  exponent = level - t->level;
  /* Set the last bits to the local ids of always choosing the last child
   * of t */
  id = (((t8_linearidx_t) 1) << T8_DTRI_DIM * exponent) - 1;
  /* Set the first bits of id to the id of t itself */
  id |= t_id << T8_DTRI_DIM * exponent;
  return id;
}

/* Construct the linear id of a descendant in a corner of t */
static t8_linearidx_t
t8_dtri_linear_id_corner_desc (const t8_dtri_t *t, int corner, int level)
{
  t8_linearidx_t id = 0, t_id, child_id;
  int it;

  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);
  T8_ASSERT (t->level <= level && level <= T8_DTRI_MAXLEVEL);

  switch (corner) {
  case 0:
    return t8_dtri_linear_id_first_desc (t, level);
    break;
  case 1: /* Falls down to case 2 in 3D */
#ifndef T8_DTRI_TO_DTET
    /* For type 0 triangles, the first corner descendant arises from always
     * taking the first child. For type 1 triangles it is always the second child. */
    child_id = t8_dtri_parenttype_beyid_to_Iloc[t->type][1];
    break;
#else
  case 2:
    /* For tets, the first corner descnendant arises from always taking the n-th
     * child, where n is the local child id corresponding to Bey-id corner. */
    child_id = t8_dtet_parenttype_beyid_to_Iloc[t->type][corner];
    break;
#endif
  case T8_DTRI_DIM:
    /* In 2D the 2nd corner and in 3D the 3rd corner correspond to the last
     * descendant of t */
    return t8_dtri_linear_id_last_desc (t, level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* This part is executed for corner 1 (2D) or corner 1 and 2 (3D) */
  t_id = t8_dtri_linear_id (t, t->level);
  for (it = 0; it < level - t->level; it++) {
    /* Store child_id at every position of the new id after t's level and
     * up to level */
    id |= child_id << (T8_DTRI_DIM * it);
  }
  /* At the beginning add the linear id of t */
  id |= t_id << (T8_DTRI_DIM * (level - t->level));
  return id;
}

t8_linearidx_t
t8_dtri_linear_id (const t8_dtri_t *t, int level)
{
  t8_linearidx_t id = 0;
  int8_t type_temp = 0;
  t8_dtri_cube_id_t cid;
  int i;
  int exponent;
  int my_level;

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  my_level = t->level;
  exponent = 0;
  /* If the given level is bigger than t's level
   * we first fill up with the ids of t's descendants at t's
   * origin with the same type as t */
  if (level > my_level) {
    exponent = (level - my_level) * T8_DTRI_DIM;
    type_temp = t->type;
    level = my_level;
  }
  else {
    type_temp = compute_type (t, level);
  }
  for (i = level; i > 0; i--) {
    cid = compute_cubeid (t, i);
    id |= ((t8_linearidx_t) t8_dtri_type_cid_to_Iloc[type_temp][cid]) << exponent;
    exponent += T8_DTRI_DIM; /* multiply with 4 (2d) resp. 8  (3d) */
    type_temp = t8_dtri_cid_type_to_parenttype[cid][type_temp];
  }
  return id;
}

void
t8_dtri_init_linear_id_with_level (t8_dtri_t *t, t8_linearidx_t id, const int start_level, const int end_level,
                                   t8_dtri_type_t parenttype)
{
  int i;
  int offset_coords, offset_index;
  const int children_m1 = T8_DTRI_CHILDREN - 1;
  t8_linearidx_t local_index;
  t8_dtri_cube_id_t cid;
  t8_dtri_type_t type;
  T8_ASSERT (id <= ((t8_linearidx_t) 1) << (T8_DTRI_DIM * end_level));
  /*Ensure, that the function is called with a valid element */
  T8_ASSERT (t->level == start_level);
  T8_ASSERT (t8_dtri_is_valid (t));

  t->level = end_level;

  type = parenttype; /* This is the type of the parent triangle */
  for (i = start_level; i <= end_level; i++) {
    offset_coords = T8_DTRI_MAXLEVEL - i;
    offset_index = end_level - i;
    /* Get the local index of T's ancestor on level i */
    local_index = (id >> (T8_DTRI_DIM * offset_index)) & children_m1;
    /* Get the type and cube-id of T's ancestor on level i */
    cid = t8_dtri_parenttype_Iloc_to_cid[type][local_index];
    type = t8_dtri_parenttype_Iloc_to_type[type][local_index];
    t->x |= (cid & 1) ? 1 << offset_coords : 0;
    t->y |= (cid & 2) ? 1 << offset_coords : 0;
#ifdef T8_DTRI_TO_DTET
    t->z |= (cid & 4) ? 1 << offset_coords : 0;
#endif
  }
  t->type = type;
}

void
t8_dtri_init_linear_id (t8_dtri_t *t, t8_linearidx_t id, int level)
{
  int i;
  int offset_coords, offset_index;
  const int children_m1 = T8_DTRI_CHILDREN - 1;
  t8_linearidx_t local_index;
  t8_dtri_cube_id_t cid;
  t8_dtri_type_t type;
  T8_ASSERT (id <= ((t8_linearidx_t) 1) << (T8_DTRI_DIM * level));

  t->level = level;
  t->x = 0;
  t->y = 0;
#ifdef T8_DTRI_TO_DTET
  t->z = 0;
#endif
  type = 0; /* This is the type of the root triangle */
  for (i = 1; i <= level; i++) {
    offset_coords = T8_DTRI_MAXLEVEL - i;
    offset_index = level - i;
    /* Get the local index of T's ancestor on level i */
    local_index = (id >> (T8_DTRI_DIM * offset_index)) & children_m1;
    /* Get the type and cube-id of T's ancestor on level i */
    cid = t8_dtri_parenttype_Iloc_to_cid[type][local_index];
    type = t8_dtri_parenttype_Iloc_to_type[type][local_index];
    t->x |= (cid & 1) ? 1 << offset_coords : 0;
    t->y |= (cid & 2) ? 1 << offset_coords : 0;
#ifdef T8_DTRI_TO_DTET
    t->z |= (cid & 4) ? 1 << offset_coords : 0;
#endif
  }
  t->type = type;
}

void
t8_dtri_init_root (t8_dtri_t *t)
{
  t->level = 0;
  t->type = 0;
  t->x = 0;
  t->y = 0;
#ifdef T8_DTRI_TO_DTET
  t->z = 0;
#endif
}

/* Stores in s the triangle that is obtained from t by going 'increment' positions
 * along the SFC of a uniform refinement of level 'level'.
 * 'increment' must be greater than -4 (-8) and smaller than +4 (+8).
 * Before calling this function s should store the same entries as t. */
static void
t8_dtri_succ_pred_recursion (const t8_dtri_t *t, t8_dtri_t *s, int level, int increment)
{
  t8_dtri_type_t type_level, type_level_p1;
  t8_dtri_cube_id_t cid;
  int local_index;
  int sign;

  /* We exclude the case level = 0, because the root triangle does
   * not have a successor. */
  T8_ASSERT (1 <= level && level <= t->level);
  T8_ASSERT (-T8_DTRI_CHILDREN < increment && increment < T8_DTRI_CHILDREN);

  if (increment == 0) {
    t8_dtri_copy (t, s);
    return;
  }
  cid = compute_cubeid (t, level);
  type_level = compute_type (t, level);
  local_index = t8_dtri_type_cid_to_Iloc[type_level][cid];
  local_index = (local_index + T8_DTRI_CHILDREN + increment) % T8_DTRI_CHILDREN;
  if (local_index == 0) {
    sign = increment < 0 ? -1 : increment > 0;
    t8_dtri_succ_pred_recursion (t, s, level - 1, sign);
    type_level_p1 = s->type; /* We stored the type of s at level-1 in s->type */
  }
  else {
    type_level_p1 = t8_dtri_cid_type_to_parenttype[cid][type_level];
  }
  type_level = t8_dtri_parenttype_Iloc_to_type[type_level_p1][local_index];
  cid = t8_dtri_parenttype_Iloc_to_cid[type_level_p1][local_index];
  s->type = type_level;
  s->level = level;
  /* Set the x,y(,z) coordinates at level to the cube-id. */
  /* TODO: check if we set the correct bits here! */
  s->x = (cid & 1 ? s->x | 1 << (T8_DTRI_MAXLEVEL - level) : s->x & ~(1 << (T8_DTRI_MAXLEVEL - level)));
  s->y = (cid & 2 ? s->y | 1 << (T8_DTRI_MAXLEVEL - level) : s->y & ~(1 << (T8_DTRI_MAXLEVEL - level)));
#ifdef T8_DTRI_TO_DTET
  s->z = (cid & 4 ? s->z | 1 << (T8_DTRI_MAXLEVEL - level) : s->z & ~(1 << (T8_DTRI_MAXLEVEL - level)));
#endif
}

void
t8_dtri_successor (const t8_dtri_t *t, t8_dtri_t *s, int level)
{
  t8_dtri_copy (t, s);
  t8_dtri_succ_pred_recursion (t, s, level, 1);
}

void
t8_dtri_first_descendant (const t8_dtri_t *t, t8_dtri_t *s, int level)
{
  t8_linearidx_t id;

  T8_ASSERT (level >= t->level);
  /* Compute the linear id of the first descendant */
  id = t8_dtri_linear_id_first_desc (t, level);
  /* The first descendant has exactly this id */
  t8_dtri_init_linear_id (s, id, level);
}

void
t8_dtri_last_descendant (const t8_dtri_t *t, t8_dtri_t *s, int level)
{
  t8_linearidx_t id;

  T8_ASSERT (level >= t->level);
  /* Compute the linear id of t's last descendant */
  id = t8_dtri_linear_id_last_desc (t, level);
  /* Set s to math this linear id */
  t8_dtri_init_linear_id (s, id, level);
}

void
t8_dtri_corner_descendant (const t8_dtri_t *t, t8_dtri_t *s, int corner, int level)
{
  t8_linearidx_t id;
  T8_ASSERT (t->level <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);

  switch (corner) {
  case 0:
    /* The 0-the corner descendant is just the first descendant */
    t8_dtri_first_descendant (t, s, level);
    break;
  case 1:
#ifdef T8_DTRI_TO_DTET /* In 3D corner 1 and 2 are handled the same */
  case 2:
#endif
    /* Compute the linear id of the descendant and construct a triangle with
     * this id */
    id = t8_dtri_linear_id_corner_desc (t, corner, level);
    t8_dtri_init_linear_id (s, id, level);
    break;
  case T8_DTRI_DIM:
    /* The 2nd corner (2D) or 3rd corner (3D) descendant is just the last descendant */
    t8_dtri_last_descendant (t, s, level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
t8_dtri_predecessor (const t8_dtri_t *t, t8_dtri_t *s, int level)
{
  t8_dtri_copy (t, s);
  t8_dtri_succ_pred_recursion (t, s, level, -1);
}

int
t8_dtri_ancestor_id (const t8_dtri_t *t, int level)
{
  t8_dtri_cube_id_t cid;
  t8_dtri_type_t type;

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (level <= t->level);

  cid = compute_cubeid (t, level);
  type = compute_type (t, level);
  return t8_dtri_type_cid_to_Iloc[type][cid];
}

int
t8_dtri_child_id (const t8_dtri_t *t)
{
  return t8_dtri_type_cid_to_Iloc[t->type][compute_cubeid (t, t->level)];
}

int
t8_dtri_get_level (const t8_dtri_t *t)
{
  return t->level;
}

int
t8_dtri_is_valid (const t8_dtri_t *t)
{
  int is_valid;
  t8_dtri_coord_t max_coord;

  /* TODO: depending on the level only certain values for the coordinates are
   *       allowed. Check if the coordinates have these values. */

  /* A triangle/tet is valid if: */
  /* The level is in the valid range */
  is_valid = 0 <= t->level && t->level <= T8_DTRI_MAXLEVEL;
  /* The coordinates are in valid ranges, we allow the x,y,z coordinates
   * to lie in the 3x3 neighborhood of the root cube. */
  max_coord = ((int64_t) 2 * T8_DTRI_ROOT_LEN) - 1;
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= t->x && t->x <= max_coord;
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= t->y && t->y <= max_coord;
#ifdef T8_DTRI_TO_DTET
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= t->z && t->z <= max_coord;
#endif
  /* Its type is in the valid range */
  is_valid = is_valid && 0 <= t->type && t->type < T8_DTRI_NUM_TYPES;

  return is_valid;
}

void
t8_dtri_init (t8_dtri_t *t)
{
  /* Set all values to zero */
  memset (t, 0, sizeof (*t));
}

/* triangles (tets) are packed as x,y (,z) coordinates, type and level */
void
t8_dtri_element_pack (t8_dtri_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                      int *position, sc_MPI_Comm comm)
{
  int mpiret;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(elements[ielem]->x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&elements[ielem]->y, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);

#ifdef T8_DTRI_TO_DTET
    mpiret = sc_MPI_Pack (&elements[ielem]->z, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
#endif

    mpiret = sc_MPI_Pack (&elements[ielem]->type, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&elements[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* triangles (tets) are packed as x,y (,z) coordinates, type and level */
void
t8_dtri_element_pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size)
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  /* coords */
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
#ifdef T8_DTRI_TO_DTET
  int coord_count = 3;
#else
  int coord_count = 2;
#endif
  singlesize += coord_count * datasize;

  /* type, level*/
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += 2 * datasize;

  *pack_size = count * singlesize;
}

/* triangles (tets) are packed as x,y (,z) coordinates, type and level */
void
t8_dtri_element_unpack (void *recvbuf, const int buffer_size, int *position, t8_dtri_t **elements,
                        const unsigned int count, sc_MPI_Comm comm)
{
  int mpiret;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(elements[ielem]->x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(elements[ielem]->y), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
#ifdef T8_DTRI_TO_DTET
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(elements[ielem]->z), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
#endif
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(elements[ielem]->type), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(elements[ielem]->level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}
