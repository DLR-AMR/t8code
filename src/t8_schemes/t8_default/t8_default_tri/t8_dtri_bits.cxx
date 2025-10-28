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
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.hxx>
#else
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
#endif
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>

T8_EXTERN_C_BEGIN ();

typedef int8_t t8_dtri_cube_id_t;

/* Compute the cube-id of t's ancestor of level "level" in constant time.
 * If "level" is greater then t->level then the cube-id 0 is returned. */
static t8_dtri_cube_id_t
compute_cubeid (const t8_dtri_t *element, int level)
{
  t8_dtri_cube_id_t id = 0;
  t8_dtri_coord_t h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  h = T8_DTRI_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((element->x & h) ? 0x01 : 0);
  id |= ((element->y & h) ? 0x02 : 0);
#ifdef T8_DTRI_TO_DTET
  id |= ((element->z & h) ? 0x04 : 0);
#endif

  return id;
}

/* A routine to compute the type of t's ancestor of level "level", if its type at an intermediate level is already
 * known. If "level" equals t's level then t's type is returned. It is not allowed to call this function with "level"
 * greater than t->level. This method runs in O(t->level - level).
 */
static t8_dtri_type_t
compute_type_ext (const t8_dtri_t *element, int level, t8_dtri_type_t known_type, int known_level)
{
  int8_t type = known_type;
  t8_dtri_cube_id_t cid;
  int i;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (known_level <= element->level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /* TODO: the type of the root tet is hardcoded to 0
     *       maybe once we want to allow the root tet to have different types */
    return 0;
  }
  for (i = known_level; i > level; i--) {
    cid = compute_cubeid (element, i);
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
compute_type (const t8_dtri_t *element, int level)
{
  return compute_type_ext (element, level, element->type, element->level);
}

void
t8_dtri_copy (const t8_dtri_t *element, t8_dtri_t *dest)
{
  if (element == dest) {
    /* Do nothing if they are already the same. */
    return;
  }
  memcpy (dest, element, sizeof (t8_dtri_t));
}

int
t8_dtri_equal (const t8_dtri_t *element1, const t8_dtri_t *element2)
{
  return (element1->level == element2->level && element1->type == element2->type && element1->x == element2->x
          && element1->y == element2->y
#ifdef T8_DTRI_TO_DTET
          && element1->z == element2->z
#endif
  );
}

int
t8_dtri_compare (const t8_dtri_t *element1, const t8_dtri_t *element2)
{
  int maxlvl;
  t8_linearidx_t id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (element1->level, element2->level);
  /* Compute the linear ids of the elements */
  id1 = t8_dtri_linear_id (element1, maxlvl);
  id2 = t8_dtri_linear_id (element2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the triangle with the smaller level is considered smaller */
    T8_ASSERT (element1->level != element2->level || t8_dtri_is_equal (element1, element2));
    return element1->level - element2->level;
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_dtri_parent (const t8_dtri_t *element, t8_dtri_t *parent)
{
  t8_dtri_cube_id_t cid;
  t8_dtri_coord_t h;

  T8_ASSERT (element->level > 0);

  h = T8_DTRI_LEN (element->level);
  /* Compute type of parent */
  cid = compute_cubeid (element, element->level);
  parent->type = t8_dtri_cid_type_to_parenttype[cid][element->type];
  /* Set coordinates of parent */
  parent->x = element->x & ~h;
  parent->y = element->y & ~h;
#ifdef T8_DTRI_TO_DTET
  parent->z = element->z & ~h;
#endif
  parent->level = element->level - 1;
}

void
t8_dtri_ancestor (const t8_dtri_t *element, int level, t8_dtri_t *ancestor)
{
  /* TODO: find out, at which level difference it is faster to use
   * the arithmetic computation of ancestor type
   * opposed to iteratively computing the parent type.
   */
  t8_dtri_coord_t delta_x, delta_y, diff_xy;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t delta_z, diff_xz, diff_yz;
  t8_dtet_type_t possible_types[6] = { 1, 1, 1, 1, 1, 1 };
  int i;
#if T8_ENABLE_DEBUG
  int set_type = 0;
#endif
#endif /* T8_DTRI_TO_DTET */

  /* delta_{x,y} = element->{x,y} - ancestor->{x,y}
   * the difference of the coordinates.
   * Needed to compute the type of the ancestor. */
  delta_x = element->x & (T8_DTRI_LEN (level) - 1);
  delta_y = element->y & (T8_DTRI_LEN (level) - 1);
#ifdef T8_DTRI_TO_DTET
  delta_z = element->z & (T8_DTRI_LEN (level) - 1);
#endif

  /* The coordinates of the ancestor. It is necessary
   * to compute the delta first, since ancestor and t
   * could point to the same triangle. */
  ancestor->x = element->x & ~(T8_DTRI_LEN (level) - 1);
  ancestor->y = element->y & ~(T8_DTRI_LEN (level) - 1);
#ifdef T8_DTRI_TO_DTET
  ancestor->z = element->z & ~(T8_DTRI_LEN (level) - 1);
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
    ancestor->type = element->type;
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
    if (element->type == 0 || element->type == 1 || element->type == 5) {
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
    if (element->type == 0 || element->type == 1 || element->type == 2) {
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
    if (element->type == 1 || element->type == 2 || element->type == 3) {
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
#if T8_ENABLE_DEBUG
      T8_ASSERT (set_type != 1);
      set_type = 1;
#endif
    }
  }
  T8_ASSERT (set_type == 1);
#endif /* T8_DTRI_TO_DTET */
  ancestor->level = level;
}

/** Compute the coordinates of the four vertices of a triangle.
 * \param [in] element         Input triangle.
 * \param [out] coordinates An array of 4x3 t8_dtri_coord_t that will be filled with the coordinates of t's vertices.
 */
static void
t8_dtri_compute_all_coords (const t8_dtri_t *element, t8_dtri_coord_t coordinates[T8_DTRI_FACES][T8_DTRI_DIM])
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

  type = element->type;
  h = T8_DTRI_LEN (element->level);
#ifndef T8_DTRI_TO_DTET
  ei = type;
#else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

  coordinates[0][0] = element->x;
  coordinates[0][1] = element->y;
#ifdef T8_DTRI_TO_DTET
  coordinates[0][2] = element->z;
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
#if T8_ENABLE_DEBUG
  /* We check whether the results are the same as with the
   * t8_dtri_compute_integer_coords function.
   */
  {
    int ivertex;
    t8_dtri_coord_t coords[T8_DTRI_DIM];
    for (ivertex = 0; ivertex < T8_DTRI_FACES; ivertex++) {
      t8_default_scheme_tri::element_get_vertex_integer_coords (element, ivertex, coords);
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
t8_dtri_child (const t8_dtri_t *element, int childid, t8_dtri_t *child)
{
  t8_dtri_t *c = (t8_dtri_t *) child;
  t8_dtri_coord_t t_coordinates[T8_DTRI_DIM];
  int vertex;
  int Bey_cid;

  T8_ASSERT (element->level < T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= childid && childid < T8_DTRI_CHILDREN);

  Bey_cid = t8_dtri_index_to_bey_number[element->type][childid];

  /* Compute anchor coordinates of child */
  if (Bey_cid == 0) {
    /* TODO: would it be better do drop this if and
     *       capture it with (t->x+t->x)>>1 below? */
    c->x = element->x;
    c->y = element->y;
#ifdef T8_DTRI_TO_DTET
    c->z = element->z;
#endif
  }
  else {
    vertex = t8_dtri_beyid_to_vertex[Bey_cid];
    /* i-th anchor coordinate of child is (X_(0,i)+X_(vertex,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    t8_default_scheme_tri::element_get_vertex_integer_coords (element, vertex, t_coordinates);
    c->x = (element->x + t_coordinates[0]) >> 1;
    c->y = (element->y + t_coordinates[1]) >> 1;
#ifdef T8_DTRI_TO_DTET
    c->z = (element->z + t_coordinates[2]) >> 1;
#endif
  }

  /* Compute type of child */
  c->type = t8_dtri_type_of_child[element->type][Bey_cid];

  c->level = element->level + 1;
}

void
t8_dtri_childrenpv (const t8_dtri_t *element, t8_dtri_t *children[T8_DTRI_CHILDREN])
{
  t8_dtri_coord_t t_coordinates[T8_DTRI_FACES][T8_DTRI_DIM];
  const int8_t level = element->level + 1;
  int i;
  int Bey_cid;
  int vertex;
  t8_dtri_type_t t_type = element->type;

  T8_ASSERT (element->level < T8_DTRI_MAXLEVEL);
  t8_dtri_compute_all_coords (element, t_coordinates);
  /* We use t_coordinates[0] for t->x,y,z to ensure that the function is valid
   * if called with t = children[0]. If we would use t->x later it would be the newly
   * computed value children[0]->x. */
  children[0]->x = t_coordinates[0][0]; /* t->x */
  children[0]->y = t_coordinates[0][1]; /* t->y */
#ifdef T8_DTRI_TO_DTET
  children[0]->z = t_coordinates[0][2]; /* t->z */
#endif
  children[0]->type = t_type;
  children[0]->level = level;
  for (i = 1; i < T8_DTRI_CHILDREN; i++) {
    Bey_cid = t8_dtri_index_to_bey_number[t_type][i];
    vertex = t8_dtri_beyid_to_vertex[Bey_cid];
    /* i-th anchor coordinate of child is (X_(0,i)+X_(vertex,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    children[i]->x = (t_coordinates[0][0] + t_coordinates[vertex][0]) >> 1;
    children[i]->y = (t_coordinates[0][1] + t_coordinates[vertex][1]) >> 1;
#ifdef T8_DTRI_TO_DTET
    children[i]->z = (t_coordinates[0][2] + t_coordinates[vertex][2]) >> 1;
#endif
    children[i]->type = t8_dtri_type_of_child[t_type][Bey_cid];
    children[i]->level = level;
#if T8_ENABLE_DEBUG
    {
      /* We check whether the child computed here equals to the child
       * computed in the t8_dtri_child function. */
      t8_dtri_t check_child;
      /* We use children[0] here instead of t, since we explicitly allow t=children[0] as input
       * and thus the values of t may be already overwritten. However the only
       * difference from children[0] to t is in the level. */
      children[0]->level--;
      t8_dtri_child (children[0], i, &check_child);
      T8_ASSERT (check_child.x == children[i]->x && check_child.y == children[i]->y);
      T8_ASSERT (check_child.type == children[i]->type && check_child.level == children[i]->level);
#ifdef T8_DTRI_TO_DTET
      T8_ASSERT (check_child.z == children[i]->z);
#endif
      children[0]->level++;
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
t8_dtri_sibling (const t8_dtri_t *element, int sibid, t8_dtri_t *sibling)
{
  T8_ASSERT (0 <= sibid && sibid < T8_DTRI_CHILDREN);
  T8_ASSERT (((const t8_dtri_t *) element)->level > 0);
  t8_dtri_parent (element, sibling);
  t8_dtri_child (sibling, sibid, sibling);
}

/* Saves the neighbour of T along face "face" in N
 * returns the face number of N along which T is its neighbour */
int
t8_dtri_face_neighbour (const t8_dtri_t *element, int face, t8_dtri_t *neigh)
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

  level = element->level;
  type_old = element->type;
  type_new = type_old;
  coords[0] = element->x;
  coords[1] = element->y;
#ifdef T8_DTRI_TO_DTET
  coords[2] = element->z;
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
  neigh->x = coords[0];
  neigh->y = coords[1];
#ifdef T8_DTRI_TO_DTET
  neigh->z = coords[2];
#endif
  neigh->level = level;
  neigh->type = type_new;
  return ret;
}

void
t8_dtri_nearest_common_ancestor (const t8_dtri_t *element1, const t8_dtri_t *element2, t8_dtri_t *nca)
{
  int maxlevel, c_level, r_level;
  t8_dtri_type_t t1_type_at_l, t2_type_at_l;
  uint32_t exclorx, exclory;
#ifdef T8_DTRI_TO_DTET
  uint32_t exclorz;
#endif
  uint32_t maxclor;

  /* Find the level of the nca */
  exclorx = element1->x ^ element2->x;
  exclory = element1->y ^ element2->y;
#ifdef T8_DTRI_TO_DTET
  exclorz = element1->z ^ element2->z;

  maxclor = exclorx | exclory | exclorz;
#else
  maxclor = exclorx | exclory;
#endif
  maxlevel = SC_LOG2_32 (maxclor) + 1;

  T8_ASSERT (maxlevel <= T8_DTRI_MAXLEVEL);

  c_level = (int8_t) SC_MIN (T8_DTRI_MAXLEVEL - maxlevel, (int) SC_MIN (element1->level, element2->level));
  r_level = c_level;
  /* c_level is the level of the nca cube surrounding t1 and t2.
   * If t1 and t2 have different types at this level, then we have to shrink this
   * level further. */
  t1_type_at_l = compute_type (element1, c_level);
  t2_type_at_l = compute_type (element2, c_level);
  while (t1_type_at_l != t2_type_at_l) {
    r_level--;
    t1_type_at_l = compute_type_ext (element1, r_level, t1_type_at_l, r_level + 1);
    t2_type_at_l = compute_type_ext (element2, r_level, t2_type_at_l, r_level + 1);
  }
  T8_ASSERT (r_level >= 0);
  /* Construct the ancestor of the first triangle at this level */
  t8_dtri_ancestor (element1, r_level, nca);
}

void
t8_dtri_children_at_face (const t8_dtri_t *element, int face, t8_dtri_t *children[],
                          __attribute__ ((unused)) int num_children, int *child_indices)
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
  child_ids[0] = face == 0 ? 1 + element->type : 0;
  /* The second child is '1 + type' for face 2 and '3' for faces 0 and 1 */
  child_ids[1] = face == 2 ? 1 + element->type : 3;
#else
  /* Tetrahedron version */
  for (i = 0; i < T8_DTET_FACE_CHILDREN; i++) {
    child_ids[i] = t8_dtet_face_child_id_by_type[element->type][face][i];
  }
#endif

  /* Compute the children at the face.
   * We revert the order to compute children[0] last, since the usage
   * allows for element == children[0].
   */
  for (i = T8_DTRI_FACE_CHILDREN - 1; i >= 0; i--) {
    t8_dtri_child (element, child_ids[i], children[i]);
  }
}

int
t8_dtri_tree_face (__attribute__ ((unused)) t8_dtri_t *element, int face)
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
  switch (element->type) {
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
t8_dtri_root_face_to_face (__attribute__ ((unused)) t8_dtri_t *element, int root_face)
{
  T8_ASSERT (0 <= root_face && root_face < T8_DTRI_FACES);
#ifndef T8_DTRI_TO_DTET
  T8_ASSERT (element->type == 0);
  /* For triangles of type 0 the face number coincides with the number of the
   * root tree face. Triangles of type 1 cannot lie on the boundary of the
   * tree and thus the return value can be arbitrary. */
  return root_face;
#else
  /* For tets only tets of type not 3 can have tree boundary faces.
   * All these tets of type not 0 (types 1, 2, 4, and 5) can only have one of
   * their faces as boundary face. */
  T8_ASSERT (element->type != 3);
  switch (element->type) {
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
t8_dtri_face_child_face (__attribute__ ((unused)) const t8_dtri_t *triangle, int face,
                         __attribute__ ((unused)) int face_child)
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
t8_dtri_is_inside_root (t8_dtri_t *element)
{
  int is_inside;

  if (element->level == 0) {
    /* A level 0 simplex is only inside the root simplex if it
     * is the root simplex. */
    return element->type == 0 && element->x == 0 && element->y == 0
#ifdef T8_DTRI_TO_DTET
           && element->z == 0
#endif
      ;
  }
  is_inside = (element->x >= 0 && element->x < T8_DTRI_ROOT_LEN) && (element->y >= 0) &&
#ifdef T8_DTRI_TO_DTET
              (element->z >= 0) &&
#endif
#ifndef T8_DTRI_TO_DTET
              (element->y - element->x <= 0) && (element->y == element->x ? element->type == 0 : 1) &&
#else
              (element->z - element->x <= 0) && (element->y - element->z <= 0)
              && (element->z == element->x ? (0 <= element->type && element->type < 3) : 1) && /* types 0, 1, 2 */
              (element->y == element->z ? (0 == element->type || 4 <= element->type) : 1) &&   /* types 0, 4, 5 */
              /* If the anchor is on the x-y-z diagonal, only type 0 tets are inside root. */
              ((element->x == element->y && element->y == element->z) ? element->type == 0 : 1) &&
#endif
              1;
#if T8_ENABLE_DEBUG
  /* Check if is_inside gives the same result as is_ancestor for the root element. */
  {
    t8_dtri_t root;
    t8_dtri_init_root (&root);
    T8_ASSERT (is_inside == t8_dtri_is_ancestor (&root, element));
  }
#endif
  return is_inside;
}

int
t8_dtri_is_root_boundary (const t8_dtri_t *element, int face)
{
  /* Dependent on the type and the face we have to check
   * different conditions */
#ifndef T8_DTRI_TO_DTET
  switch (element->type) {
  case 0:
    switch (face) {
    case 0:
      return element->x == T8_DTRI_ROOT_LEN - T8_DTRI_LEN (element->level);
    case 1:
      return element->x == element->y;
    case 2:
      return element->y == 0;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  case 1:
    return 0; /* type 1 triangles are never at the boundary */
  default:
    SC_ABORT_NOT_REACHED ();
  }
#else
  switch (element->type) {
  case 0:
    switch (face) {
    case 0:
      return element->x == T8_DTET_ROOT_LEN - T8_DTET_LEN (element->level);
    case 1:
      return element->x == element->z;
    case 2:
      return element->y == element->z;
    case 3:
      return element->y == 0;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  case 1:
    /* type 1 tets are only boundary at face 0 */
    return face == 0 && element->x == T8_DTET_ROOT_LEN - T8_DTET_LEN (element->level);
  case 2:
    /* type 1 tets are only boundary at face 2 */
    return face == 2 && element->x == element->z;
  case 3:
    /* type 3 tets are never at the root boundary */
    return 0;
  case 4:
    /* type 4 tets are only boundary at face 1 */
    return face == 1 && element->y == element->z;
  case 5:
    /* type 5 tets are only boundary at face 3 */
    return face == 3 && element->y == 0;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
  SC_ABORT_NOT_REACHED ();
  return 0; /* Prevent compiler warning */
}

int
t8_dtri_is_equal (const t8_dtri_t *element1, const t8_dtri_t *element2)
{
  return (element1->level == element2->level && element1->type == element2->type && element1->x == element2->x
          && element1->y == element2->y
#ifdef T8_DTRI_TO_DTET
          && element1->z == element2->z
#endif
  );
}

/* we check if t1 and t2 lie in the same subcube and have
 * the same level and parent type */
int
t8_dtri_is_sibling (const t8_dtri_t *element1, const t8_dtri_t *element2)
{
  t8_dtri_coord_t exclorx, exclory;
#ifdef T8_DTRI_TO_DTET
  t8_dtri_coord_t exclorz;
#endif

  t8_dtri_cube_id_t cid1, cid2;

  if (element1->level == 0) {
    return element2->level == 0 && element1->x == element2->x && element1->y == element2->y &&
#ifdef T8_DTRI_TO_DTET
           element1->z == element2->z &&
#endif
           1;
  }

  exclorx = element1->x ^ element2->x;
  exclory = element1->y ^ element2->y;
#ifdef T8_DTRI_TO_DTET
  exclorz = element1->z ^ element2->z;
#endif
  cid1 = compute_cubeid (element1, element1->level);
  cid2 = compute_cubeid (element2, element2->level);

  return (element1->level == element2->level) && ((exclorx & ~T8_DTRI_LEN (element1->level)) == 0)
         && ((exclory & ~T8_DTRI_LEN (element1->level)) == 0) &&
#ifdef T8_DTRI_TO_DTET
         ((exclorz & ~T8_DTRI_LEN (element1->level)) == 0) &&
#endif
         t8_dtri_cid_type_to_parenttype[cid1][element1->type] == t8_dtri_cid_type_to_parenttype[cid2][element2->type];
}

int
t8_dtri_is_parent (const t8_dtri_t *element, const t8_dtri_t *child)
{
  t8_dtri_cube_id_t cid;

  cid = compute_cubeid (child, child->level);
  return (element->level + 1 == child->level) && (element->x == (child->x & ~T8_DTRI_LEN (child->level)))
         && (element->y == (child->y & ~T8_DTRI_LEN (child->level))) &&
#ifdef T8_DTRI_TO_DTET
         (element->z == (child->z & ~T8_DTRI_LEN (child->level))) &&
#endif
         element->type == t8_dtri_cid_type_to_parenttype[cid][child->type] && 1;
}

int
t8_dtri_is_ancestor (const t8_dtri_t *element, const t8_dtri_t *child)
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

  if (element->level > child->level) {
    return 0;
  }
  if (element->level == child->level) {
    return t8_dtri_is_equal (element, child);
  }

  exclorx = (element->x ^ child->x) >> (T8_DTRI_MAXLEVEL - element->level);
  exclory = (element->y ^ child->y) >> (T8_DTRI_MAXLEVEL - element->level);
#ifdef T8_DTRI_TO_DTET
  exclorz = (element->z ^ child->z) >> (T8_DTRI_MAXLEVEL - element->level);
#endif

  if (exclorx == 0 && exclory == 0
#ifdef T8_DTRI_TO_DTET
      && exclorz == 0
#endif
  ) {
    /* element and child have the same cube as ancestor.
     * Now check if t has the correct type to be c's ancestor. */
    type_t = element->type;
#ifndef T8_DTRI_TO_DTET
    /* 2D */
    n1 = (type_t == 0) ? child->x - element->x : child->y - element->y;
    n2 = (type_t == 0) ? child->y - element->y : child->x - element->x;

    return !(n1 >= T8_DTRI_LEN (element->level) || n2 < 0 || n2 - n1 > 0 || (n2 == n1 && child->type == 1 - type_t));
#else
    /* 3D */
    /* Compute the coordinate directions according to this table:
     *    type(t) 0  1  2  3  4  5
     * n1         x  x  y  y  z  z
     * n2         y  z  z  x  x  y
     * dir3       z  y  x  z  y  x
     */
    /* clang-format off */
    n1 = type_t / 2 == 0 ? child->x - element->x :                      /* type(t) is 0 or 1 */
           type_t / 2 == 2 ? child->z - element->z : child->y - element->y;       /* type(t) is (4 or 5) or (2 or 3) */
    n2 = (type_t + 1) / 2 == 2 ? child->x - element->x :                /* type(t) is 3 or 4 */
           (type_t + 1) / 2 == 1 ? child->z - element->z : child->y - element->y; /* type(t) is (1 or 2) or (0 or 5) */
    dir3 = type_t % 3 == 2 ? child->x - element->x :                    /* type(t) is 2 or 5 */
             type_t % 3 == 0 ? child->z - element->z : child->y - element->y;     /* type(t) is (0 or 3) or (1 or 4) */
    sign = (type_t % 2 == 0) ? 1 : -1;
    /* clang-format on */

    type_t += 6; /* We need to compute modulo six and want
                                   to avoid negative numbers when subtracting from type_t. */

    return !(n1 >= T8_DTRI_LEN (element->level) || n2 < 0 || dir3 - n1 > 0 || n2 - dir3 > 0
             || (dir3 == n2
                 && (child->type == (type_t + sign * 1) % 6 || child->type == (type_t + sign * 2) % 6
                     || child->type == (type_t + sign * 3) % 6))
             || (dir3 == n1
                 && (child->type == (type_t - sign * 1) % 6 || child->type == (type_t - sign * 2) % 6
                     || child->type == (type_t - sign * 3) % 6))
             /* On the x-y-z diagonal only tets of the same type can be
              * ancestor of each other. */
             || (dir3 == n2 && n2 == n1 && type_t - 6 != child->type));
#endif
  }
  else {
    return 0;
  }
}

/* Compute the linear id of the first descendant of a triangle/tet */
static t8_linearidx_t
t8_dtri_linear_id_first_desc (const t8_dtri_t *element, int level)
{
  /* The id of the first descendant is the id of t in a uniform level
   * refinement */
  return t8_dtri_linear_id (element, level);
}

/* Compute the linear id of the last descendant of a triangle/tet */
static t8_linearidx_t
t8_dtri_linear_id_last_desc (const t8_dtri_t *element, int level)
{
  t8_linearidx_t id = 0, t_id;
  int exponent;

  T8_ASSERT (level >= element->level);
  /* The id of the last descendant at level consists of the id of t in
   * the first digits and then the local ids of all last children
   * (3 in 2d, 7 in 3d)
   */
  t_id = t8_dtri_linear_id (element, element->level);
  exponent = level - element->level;
  /* Set the last bits to the local ids of always choosing the last child
   * of t */
  id = (((t8_linearidx_t) 1) << T8_DTRI_DIM * exponent) - 1;
  /* Set the first bits of id to the id of t itself */
  id |= t_id << T8_DTRI_DIM * exponent;
  return id;
}

/* Construct the linear id of a descendant in a corner of t */
static t8_linearidx_t
t8_dtri_linear_id_corner_desc (const t8_dtri_t *element, int corner, int level)
{
  t8_linearidx_t id = 0, t_id, child_id;
  int it;

  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);
  T8_ASSERT (element->level <= level && level <= T8_DTRI_MAXLEVEL);

  switch (corner) {
  case 0:
    return t8_dtri_linear_id_first_desc (element, level);
    break;
  case 1: /* Falls down to case 2 in 3D */
#ifndef T8_DTRI_TO_DTET
    /* For type 0 triangles, the first corner descendant arises from always
     * taking the first child. For type 1 triangles it is always the second child. */
    child_id = t8_dtri_parenttype_beyid_to_Iloc[element->type][1];
    break;
#else
  case 2:
    /* For tets, the first corner descendant arises from always taking the n-th
     * child, where n is the local child id corresponding to Bey-id corner. */
    child_id = t8_dtet_parenttype_beyid_to_Iloc[element->type][corner];
    break;
#endif
  case T8_DTRI_DIM:
    /* In 2D the 2nd corner and in 3D the 3rd corner correspond to the last
     * descendant of element. */
    return t8_dtri_linear_id_last_desc (element, level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* This part is executed for corner 1 (2D) or corner 1 and 2 (3D) */
  t_id = t8_dtri_linear_id (element, element->level);
  for (it = 0; it < level - element->level; it++) {
    /* Store child_id at every position of the new id after t's level and
     * up to level */
    id |= child_id << (T8_DTRI_DIM * it);
  }
  /* At the beginning add the linear id of t */
  id |= t_id << (T8_DTRI_DIM * (level - element->level));
  return id;
}

t8_linearidx_t
t8_dtri_linear_id (const t8_dtri_t *element, int level)
{
  t8_linearidx_t id = 0;
  int8_t type_temp = 0;
  t8_dtri_cube_id_t cid;
  int i;
  int exponent;
  int my_level;

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  my_level = element->level;
  exponent = 0;
  /* If the given level is bigger than t's level
   * we first fill up with the ids of t's descendants at t's
   * origin with the same type as t */
  if (level > my_level) {
    exponent = (level - my_level) * T8_DTRI_DIM;
    type_temp = element->type;
    level = my_level;
  }
  else {
    type_temp = compute_type (element, level);
  }
  for (i = level; i > 0; i--) {
    cid = compute_cubeid (element, i);
    id |= ((t8_linearidx_t) t8_dtri_type_cid_to_Iloc[type_temp][cid]) << exponent;
    exponent += T8_DTRI_DIM; /* multiply with 4 (2d) resp. 8  (3d) */
    type_temp = t8_dtri_cid_type_to_parenttype[cid][type_temp];
  }
  return id;
}

void
t8_dtri_init_linear_id_with_level (t8_dtri_t *element, t8_linearidx_t id, const int start_level, const int end_level,
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
  T8_ASSERT (element->level == start_level);
  T8_ASSERT (t8_dtri_is_valid (element));

  element->level = end_level;

  type = parenttype; /* This is the type of the parent triangle */
  for (i = start_level; i <= end_level; i++) {
    offset_coords = T8_DTRI_MAXLEVEL - i;
    offset_index = end_level - i;
    /* Get the local index of T's ancestor on level i */
    local_index = (id >> (T8_DTRI_DIM * offset_index)) & children_m1;
    /* Get the type and cube-id of T's ancestor on level i */
    cid = t8_dtri_parenttype_Iloc_to_cid[type][local_index];
    type = t8_dtri_parenttype_Iloc_to_type[type][local_index];
    element->x |= (cid & 1) ? 1 << offset_coords : 0;
    element->y |= (cid & 2) ? 1 << offset_coords : 0;
#ifdef T8_DTRI_TO_DTET
    element->z |= (cid & 4) ? 1 << offset_coords : 0;
#endif
  }
  element->type = type;
}

void
t8_dtri_init_linear_id (t8_dtri_t *element, t8_linearidx_t id, int level)
{
  int i;
  int offset_coords, offset_index;
  const int children_m1 = T8_DTRI_CHILDREN - 1;
  t8_linearidx_t local_index;
  t8_dtri_cube_id_t cid;
  t8_dtri_type_t type;
  T8_ASSERT (id <= ((t8_linearidx_t) 1) << (T8_DTRI_DIM * level));

  element->level = level;
  element->x = 0;
  element->y = 0;
#ifdef T8_DTRI_TO_DTET
  element->z = 0;
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
    element->x |= (cid & 1) ? 1 << offset_coords : 0;
    element->y |= (cid & 2) ? 1 << offset_coords : 0;
#ifdef T8_DTRI_TO_DTET
    element->z |= (cid & 4) ? 1 << offset_coords : 0;
#endif
  }
  element->type = type;
}

void
t8_dtri_init_root (t8_dtri_t *element)
{
  element->level = 0;
  element->type = 0;
  element->x = 0;
  element->y = 0;
#ifdef T8_DTRI_TO_DTET
  element->z = 0;
#endif
}

/* Stores in s the triangle that is obtained from t by going 'increment' positions
 * along the SFC of a uniform refinement of level 'level'.
 * 'increment' must be greater than -4 (-8) and smaller than +4 (+8).
 * Before calling this function s should store the same entries as t. */
static void
t8_dtri_succ_pred_recursion (const t8_dtri_t *element, t8_dtri_t *s, int level, int increment)
{
  t8_dtri_type_t type_level, type_level_p1;
  t8_dtri_cube_id_t cid;
  int local_index;
  int sign;

  /* We exclude the case level = 0, because the root triangle does
   * not have a successor. */
  T8_ASSERT (1 <= level && level <= element->level);
  T8_ASSERT (-T8_DTRI_CHILDREN < increment && increment < T8_DTRI_CHILDREN);

  if (increment == 0) {
    t8_dtri_copy (element, s);
    return;
  }
  cid = compute_cubeid (element, level);
  type_level = compute_type (element, level);
  local_index = t8_dtri_type_cid_to_Iloc[type_level][cid];
  local_index = (local_index + T8_DTRI_CHILDREN + increment) % T8_DTRI_CHILDREN;
  if (local_index == 0) {
    sign = increment < 0 ? -1 : increment > 0;
    t8_dtri_succ_pred_recursion (element, s, level - 1, sign);
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
t8_dtri_successor (const t8_dtri_t *element, t8_dtri_t *s, int level)
{
  t8_dtri_copy (element, s);
  t8_dtri_succ_pred_recursion (element, s, level, 1);
}

void
t8_dtri_first_descendant (const t8_dtri_t *element, t8_dtri_t *s, int level)
{
  t8_linearidx_t id;

  T8_ASSERT (level >= element->level);
  /* Compute the linear id of the first descendant */
  id = t8_dtri_linear_id_first_desc (element, level);
  /* The first descendant has exactly this id */
  t8_dtri_init_linear_id (s, id, level);
}

void
t8_dtri_last_descendant (const t8_dtri_t *element, t8_dtri_t *s, int level)
{
  t8_linearidx_t id;

  T8_ASSERT (level >= element->level);
  /* Compute the linear id of t's last descendant */
  id = t8_dtri_linear_id_last_desc (element, level);
  /* Set s to math this linear id */
  t8_dtri_init_linear_id (s, id, level);
}

void
t8_dtri_corner_descendant (const t8_dtri_t *element, t8_dtri_t *s, int corner, int level)
{
  t8_linearidx_t id;
  T8_ASSERT (element->level <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (0 <= corner && corner < T8_DTRI_CORNERS);

  switch (corner) {
  case 0:
    /* The 0-the corner descendant is just the first descendant */
    t8_dtri_first_descendant (element, s, level);
    break;
  case 1:
#ifdef T8_DTRI_TO_DTET /* In 3D corner 1 and 2 are handled the same */
  case 2:
#endif
    /* Compute the linear id of the descendant and construct a triangle with
     * this id */
    id = t8_dtri_linear_id_corner_desc (element, corner, level);
    t8_dtri_init_linear_id (s, id, level);
    break;
  case T8_DTRI_DIM:
    /* The 2nd corner (2D) or 3rd corner (3D) descendant is just the last descendant */
    t8_dtri_last_descendant (element, s, level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

void
t8_dtri_predecessor (const t8_dtri_t *element, t8_dtri_t *s, int level)
{
  t8_dtri_copy (element, s);
  t8_dtri_succ_pred_recursion (element, s, level, -1);
}

int
t8_dtri_ancestor_id (const t8_dtri_t *element, int level)
{
  t8_dtri_cube_id_t cid;
  t8_dtri_type_t type;

  T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
  T8_ASSERT (level <= element->level);

  cid = compute_cubeid (element, level);
  type = compute_type (element, level);
  return t8_dtri_type_cid_to_Iloc[type][cid];
}

int
t8_dtri_child_id (const t8_dtri_t *element)
{
  return t8_dtri_type_cid_to_Iloc[element->type][compute_cubeid (element, element->level)];
}

int
t8_dtri_get_level (const t8_dtri_t *element)
{
  return element->level;
}

int
t8_dtri_is_valid (const t8_dtri_t *element)
{
  int is_valid;
  t8_dtri_coord_t max_coord;

  /* TODO: depending on the level only certain values for the coordinates are
   *       allowed. Check if the coordinates have these values. */

  /* A triangle/tet is valid if: */
  /* The level is in the valid range */
  is_valid = 0 <= element->level && element->level <= T8_DTRI_MAXLEVEL;
  /* The coordinates are in valid ranges, we allow the x,y,z coordinates
   * to lie in the 3x3 neighborhood of the root cube. */
  max_coord = ((int64_t) 2 * T8_DTRI_ROOT_LEN) - 1;
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= element->x && element->x <= max_coord;
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= element->y && element->y <= max_coord;
#ifdef T8_DTRI_TO_DTET
  is_valid = is_valid && -T8_DTRI_ROOT_LEN <= element->z && element->z <= max_coord;
#endif
  /* Its type is in the valid range */
  is_valid = is_valid && 0 <= element->type && element->type < T8_DTRI_NUM_TYPES;

  return is_valid;
}

void
t8_dtri_init (t8_dtri_t *element)
{
  /* Set all values to zero */
  memset (element, 0, sizeof (*element));
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

T8_EXTERN_C_END ();
