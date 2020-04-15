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

#include "t8_dpyramid_bits.h"
#include "t8_dtet_bits.h"
#include <sc_functions.h>

typedef int8_t      t8_dpyramid_cube_id_t;

/*The type of a pyramid depending on the parent pyramid and its local index
 *type = (parent_type, local_index)
 */
const t8_dpyramid_type_t t8_dpyramid_parenttype_Iloc_to_type[2][10] = {
  {6, 3, 6, 0, 6, 0, 3, 6, 7, 6},
  {7, 0, 3, 6, 7, 3, 7, 0, 7, 7}
};

/*The cube Id of a pyramid depending on its parenttype and local index*/
const t8_dpyramid_cube_id_t t8_dpyramid_parenttype_Iloc_to_cid[2][10] = {
  {0, 1, 1, 2, 2, 3, 3, 3, 3, 7},
  {0, 4, 4, 4, 4, 5, 5, 6, 6, 7}
};

/* The local ID of an element in a pyramid. This is important!
 * The local ID is different, if the element is in a tet*/
const int           t8_dpyramid_type_cid_to_Iloc[8][8] = {
    {-1, -1, 3, 5, 1, -1, 7, -1},
    {-1,-1,-1,-1,-1,-1,-1,-1,},
    {-1,-1,-1,-1,-1,-1,-1,-1,},
    {-1, 1, -1, 6, 2, 5, -1, -1},
    {-1,-1,-1,-1,-1,-1,-1,-1,},
    {-1,-1,-1,-1,-1,-1,-1,-1,},
  {0, 2, 4, 7, 3, -1, -1, 9},
  {0, -1, -1, 8, 4, 6, 8, 9}
};

const t8_dpyramid_type_t t8_dpyramid_cid_type_to_parenttype[8][8] = {
  {0, 1, 2, 3, 4, 5, 6, 7},
  {0, 1, 1, 1, 0, 0, 6, -1},
  {2, 2, 2, 3, 3, 3, 6, -1},
  {1, 1, 2, 2, 2, 1, 6, 6},
  {5, 5, 4, 4, 4, 5, 7, 7},
  {0, 0, 0, 5, 5, 5, -1, 7},
  {4, 3, 3, 3, 4, 4 - 1, 7},
  {0, 1, 2, 3, 4, 5, 6, 7}
};

static              t8_dpyramid_cube_id_t
compute_cubeid (const t8_dpyramid_t * p, int level)
{
  t8_dpyramid_cube_id_t id = 0;
  t8_dpyramid_coord_t h;


  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((p->x & h) ? 0x01 : 0);
  id |= ((p->y & h) ? 0x02 : 0);
  id |= ((p->z & h) ? 0x04 : 0);

  return id;
}

int
t8_dpyramid_is_equal (const t8_dpyramid_t * p, const t8_dpyramid_t * q)
{
  if (p->x == q->x && p->y == q->y && p->z == q->z && p->type == q->type
      && p->level == q->level) {
    return 0;
  }
  else {
    return 1;
  }
}

/*Copies a pyramid from source to dest*/
void
t8_dpyramid_copy (const t8_dpyramid_t * source, t8_dpyramid_t * dest)
{
  if (source == dest) {
    return;
  }
  memcpy (dest, source, sizeof (t8_dpyramid_t));
}

int
t8_dpyramid_compare (const t8_dpyramid_t * p1, const t8_dpyramid_t * p2)
{
  int                 maxlvl;
  uint64_t            id1, id2;
  maxlvl = SC_MAX (p1->level, p2->level);

  id1 = t8_dpyramid_linear_id (p1, maxlvl);
  id2 = t8_dpyramid_linear_id (p2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the pyramid with the smaller level
     * is considered smaller */
    return p1->level - p2->level;
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
     id2 */
  return id1 < id2 ? -1 : id1 != id2;
}


int
t8_dpyramid_get_level (const t8_dpyramid_t * p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  return p->level;
}

void
t8_dpyramid_init_linear_id (t8_dpyramid_t * p, int level, uint64_t id)
{
  t8_dpyramid_type_t  type;
  t8_linearidx_t      local_index;
  t8_dpyramid_cube_id_t cid;
  int                 i;
  int                 offset_coords;
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= id && id <= sc_intpow64u (T8_DPYRAMID_CHILDREN, level));

  p->level = level;
  p->x = 0;
  p->y = 0;
  p->z = 0;
  type = 6;                     /*This is the type of the root pyramid */

  for (i = 1; i <= level; i++) {
    offset_coords = T8_DPYRAMID_MAXLEVEL - i;
    local_index = id % T8_DPYRAMID_CHILDREN;

    type = t8_dpyramid_parenttype_Iloc_to_type[type - 6][local_index];
    // Thy types of the tetrahedron children of pyramid are always 0 or 3
    if (0 <= type && type < T8_DTET_NUM_TYPES) {
      /* IDEA: After this step, the pyra and the tet index match.
       * Therefore it is sufficient to compute the tet by the remaining
       * digits of the index. Is there a tet function that computes the
       * tet given a starting level, an end level and an id?
       */
      t8_dtet_init_linear_id_with_level ((t8_dtet_t *) p, (t8_linearidx_t) id,
                                         i, level, type);
      return;
    }
    else {
      cid = t8_dpyramid_parenttype_Iloc_to_cid[T8_DPYRAMID_NUM_TYPES
                                               - type - 1][local_index];
      p->x |= (cid % 2 == 1) ? 1 << offset_coords : 0;
      p->y |= (cid == 2 || cid == 3 || cid == 6
               || cid == 7) ? 1 << offset_coords : 0;
      p->z |= (cid > 3) ? 1 << offset_coords : 0;
    }
    id /= T8_DPYRAMID_CHILDREN;
  }
  p->type = type;
}

/* TODO: What if I am a tet child of a pyramid*/
uint64_t
t8_dpyramid_linear_id (const t8_dpyramid_t * p, int level)
{                               /*
                                   uint64_t                id = 0;
                                   t8_dpyramid_type_t      temp_type = p->type;
                                   t8_dpyramid_cube_id_t   cid;
                                   int i;

                                   T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
                                   for(i = level; i>0; i--)
                                   {
                                   cid = compute_cubeid(p, i);
                                   } */
  return 0;
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                              int level)
{
  if (p->type == 6 || p->type == 7) {
    /*The first descendant of a pyramid has the same anchor coords, but another level */
    t8_dpyramid_copy (p, desc);
    desc->level = level;
  }
  else {

    t8_dtet_first_descendant ((const t8_dtet_t *) p, (t8_dtet_t *) desc,
                              level);
  }
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                             int level)
{
  if (p->type == 6 || p->type == 7) {
    /*The last descendant of a pyramid has a shifted anchor coord and another level */
    t8_dpyramid_copy (p, desc);
    desc->level = level;
    int                 coord_offset =
      T8_DPYRAMID_LEN (p->level) - T8_DPYRAMID_LEN (level);
    desc->x |= coord_offset;
    desc->y |= coord_offset;
    desc->z |= coord_offset;
  }
  else {
    t8_dtet_last_descendant ((const t8_dtet_t *) p, (t8_dtet_t *) desc,
                             level);
  }
}

int
t8_dpyramid_num_vertices (const t8_dpyramid_t * p)
{
  if (p->type < 6) {
    return T8_DTET_CORNERS;
  }
  else {
    return T8_DPYRAMID_CORNERS;
  }
}

int
t8_dpyramid_num_children (const t8_dpyramid_t * p)
{
  if (t8_dpyramid_shape(p) == T8_ECLASS_TET) {
    return T8_DTET_CHILDREN;
  }
  else {
    return T8_DPYRAMID_CHILDREN;
  }
}


int t8_dpyramid_child_id_unknown_parent(const t8_dpyramid_t * p,
                                        t8_dpyramid_t * parent)
{
    t8_dpyramid_cube_id_t   cid = compute_cubeid(p, p->level);
    t8_dpyramid_type_t      type = p->type;
    T8_ASSERT(p->level > 0);
    t8_dpyramid_parent(p, parent);
    if(t8_dpyramid_shape(parent) == T8_ECLASS_PYRAMID){
       return t8_dpyramid_type_cid_to_Iloc[type][cid];
    }
    else{
        return t8_dtet_child_id((const t8_dtet_t *) p);
    }

}

int t8_dpyramid_child_id_known_parent(const t8_dpyramid_t * p,
                                       t8_dpyramid_t * parent)
{
    t8_dpyramid_cube_id_t   cid = compute_cubeid(p, p->level);
    if(t8_dpyramid_shape(parent) == T8_ECLASS_PYRAMID){
       return t8_dpyramid_type_cid_to_Iloc[p->type][cid];
    }
    else{
        return t8_dtet_child_id((const t8_dtet_t *) p);
    }
}

int
t8_dpyramid_child_id (const t8_dpyramid_t * p)
{
  T8_ASSERT (p->level > 0);
  t8_dpyramid_t parent;
  return t8_dpyramid_child_id_unknown_parent(p, &parent);
}

void
t8_dpyramid_child (const t8_dpyramid_t * elem, int child_id,
                   t8_dpyramid_t * child)
{
  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_coord_t h;
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
  if (t8_dpyramid_shape(elem) == T8_ECLASS_TET) {
    t8_dtet_child ((t8_dtet_t *) elem, child_id, (t8_dtet_t *) child);
  }
  else {
    cube_id = t8_dpyramid_parenttype_Iloc_to_cid[elem->type - 6][child_id];
    child->level = elem->level + 1;
    h = T8_DPYRAMID_LEN (child->level);
    child->x = elem->x + ((cube_id & 0x01) ? h : 0);
    child->y = elem->y + ((cube_id & 0x02) ? h : 0);
    child->z = elem->z + ((cube_id & 0x04) ? h : 0);
    child->type =
      t8_dpyramid_parenttype_Iloc_to_type[elem->type - 6][child_id];
  }
}

const t8_dpyramid_type_t t8_dpyramid_type_Iloc_to_parenttype[2][10] = {
  {6, -1, 6, 7, 6, -1, -1, 6, -1, 6},
  {7, -1, -1, 7, 7, -1, 7, -1, 7, 7}
};

const t8_dpyramid_type_t t8_dpyramid_type_cid_to_parenttype[2][8] = {
  {6, 6, 6, 6, 7, -1, -1, 6},
  {7, -1, -1, 6, 7, 7, 7, 7}
};

const int           t8_dpyramid_trailing_zeroes_lookup[37] = { 32, 0, 1,
  26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11,
  0, 13, 4, 7, 17, 0, 25, 22, 31, 15, 29,
  10, 12, 6, 0, 21, 14, 9, 5, 20, 8, 19,
  18
};

/*Get the number of trailing 0 of an 32 bit integer*/
/*Uses the fact, that the values corresponding to the 32 positions are relatively prime to 37
 *thus we can use a lookup-table*/
int
t8_dpyramid_trailing_zeroes (t8_dpyramid_coord_t x)
{
  return t8_dpyramid_trailing_zeroes_lookup[(-x & x) % 37];
}

/* Compute the next possible "significant point" to reach from the ankercoord of p
 * If this point is reachable, return true, else return false */

int
t8_dpyramid_hit_point (const t8_dpyramid_t * p)
{
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (p->type == 0 || p->type == 3);
  t8_dpyramid_coord_t x = p->x, y = p->y, h =
    T8_DPYRAMID_LEN (p->level), n, shift;
  int                 even = (p->z >> (T8_DPYRAMID_MAXLEVEL - p->level)) % 2;
  printf ("hitIN: %i %i %i %i %i\n", p->x, p->y, p->z, p->type, p->level);

  if (even == 0) {
    if ((p->x >> (T8_DPYRAMID_MAXLEVEL - p->level)) % 2 == 0) {
      x += h;
    }
    if ((p->y >> (T8_DPYRAMID_MAXLEVEL - p->level)) % 2 == 0) {
      y += h;
    }
    printf ("intermediate hit-point: %i %i\n", x, y);
    if ((p->z >> (T8_DPYRAMID_MAXLEVEL - p->level + 1)) % 2 == 1) {
      printf ("correction\n");
      if ((y & 11) == 3 && (x & 11) != 3) {
        x = x | (1 << 1);
      }
      if ((y & 11) != 3 && (x & 11) == 3) {
        y = y | (1 << 1);
      }
    }
    if (x <= p->z || y <= p->z) {
      printf ("no hit 0\n");
      return 0;
    }
    printf ("even: %i hitpoint: %i %i\n", even, x, y);
    if (p->x == x && p->y == y) {
      printf ("hit0.1\n");
      return 1;
    }
    else if (p->x + h == x && p->y == y && p->type == 0) {
      printf ("hit0.2\n");
      return 1;
    }
    else if (p->y + h == y && p->x == x && p->type == 3) {
      printf ("hit0.3\n");
      return 1;
    }
    else {
      return 0;
    }
  }
else {  /*
       x = x>>(T8_DPYRAMID_MAXLEVEL - p->level + 1);
       y = y>>(T8_DPYRAMID_MAXLEVEL - p->level + 1);
       if(x%2 != y%2 && x%4 != 0){
       x = x|1;
       y = y|1;
       }
       x = x<<(T8_DPYRAMID_MAXLEVEL - p->level + 1);
       y = y<<(T8_DPYRAMID_MAXLEVEL - p->level + 1); */
    x =(p->x >> (T8_DPYRAMID_MAXLEVEL - p->level + 1)) << (T8_DPYRAMID_MAXLEVEL -
                                                       p->level + 1);
    n = t8_dpyramid_trailing_zeroes (x);
    y = (p->y >> (n)) << (n);
    printf ("shift1: %i %i\n", x, y);

    shift = 1 << (n);
    y = y | shift;
    printf ("zeroes: %i, shift: %i, y: %i\n", n, shift, y);

    printf ("hitPoint: %i %i\n", x, y);
    if (x <= p->z || y <= p->z) {
      return 0;
    }
    if (p->x == x && p->y == y) {
      printf ("hit2.1\n");
      return 1;
    }
    else if ((p->x - h) == x && p->y == y && p->type == 3) {
      printf ("hit2.2\n");
      return 1;
    }
    else if (p->x == x && (p->y - h) == y && p->type == 0) {
      printf ("hit2.3\n");
      return 1;
    }
    else {
      printf ("no hit 2\n");
      return 0;
    }
  }

}

/* Check, if a pyramid in the shape of a tet lies inside a tetrahedron
 * The i first bits give the anchorcoordinate for a possible ancestor of level i
 * for p. */
int
t8_dpyramid_is_inside_tet (const t8_dpyramid_t * p)
{
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (p->type == 0 || p->type == 3);
  int                 i;
  /*the tet is initialized, the ancestor will be computed */
  t8_dtet_t           tet, ancestor;
  tet.x = 0;
  tet.y = 0;
  tet.z = 0;
  for (i = 1; i < p->level; i++) {
    /*Update the coordinate of tet to i first bits */
    tet.x = tet.x | (p->x & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.y = tet.y | (p->y & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.z = tet.z | (p->z & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.level = i;
    /*Compute the ancestor */
    t8_dtet_ancestor ((const t8_dtet_t *) p, i, &ancestor);
    /*Only compare, if the ancestor has type 0 or 3 */
    if (ancestor.type == 0 || ancestor.type == 3) {
      /*Set tet.type to the type of the possible ancestor */
      tet.type = ancestor.type;
      /*Compare */
      if (t8_dtet_is_equal (&tet, &ancestor)) {
        return 0;
      }
    }
  }
  /*No matching tet-ancestor was found, the parent is a pyramid */
  return 1;
}

void
t8_dpyramid_tetparent_type (const t8_dpyramid_t * p, t8_dpyramid_t * parent)
{
  if ((p->z >> (T8_DPYRAMID_MAXLEVEL - p->level)) % 2 == 0) {
    parent->type = 6;
  }
  else {
    parent->type = 7;
  }
}



void
t8_dpyramid_parent (const t8_dpyramid_t * p, t8_dpyramid_t * parent)
{
  T8_ASSERT (p->level >= 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);

  t8_dpyramid_coord_t h = T8_DPYRAMID_LEN (p->level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The parent of a pyramid is a pyramid, maybe of different type */

    int                 cube_id = compute_cubeid(p, p->level);

    parent->type = t8_dpyramid_type_cid_to_parenttype[p->type - 6][cube_id];
    parent->x = p->x & ~h;
    parent->y = p->y & ~h;
    parent->z = p->z & ~h;
    T8_ASSERT (parent->type >= 0);
    parent->level = p->level - 1;
  }
  else if (p->type != 0 && p->type != 3) {
      /* The direct tet-child of a pyra has type 0 or type 3, therefore
       * in this case the parent is a tetrahedron*/
      t8_dtet_parent ((t8_dtet_t *) p, (t8_dtet_t *) parent);
  }
  else if (t8_dpyramid_is_inside_tet (p) == 0) {
      /*Pyramid- / tetparent detection */
      /*If a tetrahedron does not reach a "significant point" its parent is a tet */
      /*Tetcase */;
      t8_dtet_parent ((t8_dtet_t *) p, (t8_dtet_t *) parent);
   }
   else {
      /*p does not lie in larger tet => parent is pyra */

      t8_dpyramid_tetparent_type (p, parent);
      parent->x = p->x & ~h;
      parent->y = p->y & ~h;
      parent->z = p->z & ~h;
      parent->level = p->level - 1;
   }
}

t8_eclass_t
t8_dpyramid_shape (const t8_dpyramid_t * p)
{
  /*The pyramid has the shape of a tetrahedron */
  if (p->type < 6) {
    return T8_ECLASS_TET;
  }
  else {
    return T8_ECLASS_PYRAMID;
  }

}


t8_dpyramid_type_t
compute_type (const t8_dpyramid_t * p, int level)
{
  t8_dpyramid_cube_id_t cid;
  t8_dpyramid_type_t  type = p->type;
  int                 i;
  if (level == p->level) {
    return p->type;
  }
  if (level == 0) {
    /*Type of the root pyra */
    return T8_DPYRAMID_ROOT_TPYE;
  }
  for (i = p->level; i > level; i--) {
    cid = compute_cubeid (p, i);
    type = t8_dpyramid_cid_type_to_parenttype[cid][type];
  }

  return type;
}

void
t8_dpyramid_successor_recursion (const t8_dpyramid_t * elem, t8_dpyramid_t * succ,
                       t8_dpyramid_t * parent, int level)
{
  int                 child_id, num_children;
  t8_dpyramid_copy (elem, succ);
  T8_ASSERT (1 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  succ->level = level;
  T8_ASSERT (succ->type >= 0);
  child_id = t8_dpyramid_child_id_unknown_parent(elem, parent);
  //t8_dpyramid_copy(parent, succ);
  /*Compute number of children*/
  num_children = t8_dpyramid_num_children(parent);
  T8_ASSERT (0 <= child_id
             && child_id < num_children);
  if (child_id == num_children - 1) {
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion(succ, succ, parent, level - 1);
    succ->level = level;
    /* bits auf level auf child 0 setzen */
    succ->x =
     (succ->x >> (T8_DPYRAMID_MAXLEVEL - level + 1))<< (T8_DPYRAMID_MAXLEVEL
                                                        - level + 1);
    succ->y =
      (succ->y >> (T8_DPYRAMID_MAXLEVEL - level + 1)) << (T8_DPYRAMID_MAXLEVEL -
                                                    level + 1);
    succ->z =
      (succ->z >> (T8_DPYRAMID_MAXLEVEL - level + 1)) << (T8_DPYRAMID_MAXLEVEL -
                                                    level + 1);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1*/
    t8_dpyramid_child (parent, child_id + 1, succ);
    succ->level = level;
  }
}

void
t8_dpyramid_successor (const t8_dpyramid_t * elem, t8_dpyramid_t * succ,
                       int level)
{
    t8_dpyramid_t       parent;
    t8_dpyramid_successor_recursion(elem, succ, &parent, level);
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t * p, int vertex, int coords[])
{
  t8_dpyramid_coord_t h;
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_VERTICES);

  if (p->type == 6 || p->type == 7) {
    h = T8_DPYRAMID_LEN (p->level);
    coords[0] = p->x;
    coords[1] = p->y;
    coords[2] = p->z;
    switch (vertex) {
    case 0:
      coords[2] += (p->type == 7) ? h : 0;
      break;
    case 1:
      coords[0] += h;
      coords[2] += (p->type == 7) ? h : 0;
      break;
    case 2:
      coords[1] += h;
      coords[2] += (p->type == 7) ? h : 0;
      break;
    case 3:
      coords[0] += h;
      coords[1] += h;
      coords[2] += (p->type == 7) ? h : 0;
      break;
    case 4:
      coords[0] += (p->type == 6) ? h : 0;
      coords[1] += (p->type == 6) ? h : 0;
      coords[2] += (p->type == 6) ? h : 0;
      break;
    }
  }
  else {
    T8_ASSERT (vertex < T8_DTET_CORNERS);
    t8_dtet_compute_coords ((const t8_dtet_t *) p, vertex, coords);
  }
}
