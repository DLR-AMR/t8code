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
#include "t8_dpyramid_connectivity.h"
#include <sc_functions.h>

typedef int8_t      t8_dpyramid_cube_id_t;


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

int
t8_dpyramid_is_family (const t8_dpyramid_t ** fam)
{
  /*TODO: This can be implemented better! */
  const int level = fam[0]->level;
  int i, ptype;
  t8_dpyramid_coord_t inc = T8_DPYRAMID_LEN(level), x_inc, y_inc;
  if (t8_dpyramid_shape (fam[0]) == T8_ECLASS_TET) {
    return t8_dtet_is_familypv ((const t8_dtet_t **) fam);
  }
  else {
      if(level == 0){
          return 0;
      }
      /*The type of parent is the type of the first child in z-curve-order*/
      ptype = fam[0]->type;
      T8_ASSERT(ptype == 6 || ptype == 7);
      for(i = 1; i<T8_DPYRAMID_CHILDREN; i++){
          /*All elements must have the same level to be a family*/
          if(fam[i]->level != level)
          {
              return 0;
          }
          /*Check if every family-member has the correct type*/
          if(t8_dpyramid_parenttype_Iloc_to_type[ptype][i] != fam[i]->type){
              return 0;
          }
      }

      x_inc = fam[0]->x + inc;
      y_inc = fam[0]->y + inc;
      /*Check the coordinates of the anchor-coordinate*/
      if(ptype == 6){
          return fam[0]->z == fam[1]->z && fam[0]->z == fam[2]->z && fam[0]->z == fam[3]->z &&
               fam[0]->z == fam[4]->z && fam[0]->z == fam[5]->z && fam[0]->z == fam[6]->z &&
               fam[0]->z == fam[7]->z && fam[0]->z == fam[8]->z &&
               fam[0]->z == (fam[9]->z - inc) &&
               fam[0]->x == fam[3]->x && fam[0]->x == fam[4]->x && x_inc == fam[1]->x &&
               x_inc == fam[2]->x && x_inc == fam[5]->x && x_inc == fam[6]->x &&
               x_inc == fam[7]->x &&x_inc == fam[8]->x && x_inc == fam[9]->x &&
               fam[0]->y == fam[1]->y && fam[0]->y == fam[2]->y && y_inc == fam[3]->y &&
               y_inc == fam[4]->y && y_inc == fam[5]->y && y_inc == fam[6]->y &&
               y_inc == fam[7]->y && y_inc == fam[8]->y && y_inc == fam[9]->y;
      }
      else{
           return fam[1]->z == fam[0]->z + inc && fam[1]->z == fam[2]->z &&
                 fam[1]->z == fam[3]->z && fam[1]->z == fam[4]->z && fam[1]->z == fam[5]->z &&
                 fam[1]->z == fam[6]->z && fam[1]->z == fam[7]->z && fam[1]->z == fam[8]->z &&
                 fam[1]->z == fam[9]->z &&
                 fam[0]->x == fam[1]->x && fam[0]->x == fam[2]->x && fam[0]->x == fam[3]->x &&
                 fam[0]->x == fam[4]->x && fam[0]->x == fam[7]->x && fam[0]->x == fam[8]->x &&
                 x_inc == fam[5]->x && x_inc == fam[6]->x && x_inc == fam[9]->x &&
                 fam[0]->y == fam[1]->y && fam[0]->y == fam[2]->y && fam[0]->y == fam[3]->y &&
                 fam[0]->y == fam[4]->y && fam[0]->y == fam[5]->y && fam[0]->y == fam[6]->y &&
                 y_inc == fam[7]->y && y_inc == fam[8]->y && y_inc == fam[9]->y;
      }
    }
}



int
t8_dpyramid_is_root_boundary(const t8_dpyramid_t * p, int face)
{
    T8_ASSERT(0 <= face && face <= T8_DPYRAMID_FACES);
    switch(p->type){
    /*Doublecheck the tet-part*/
    case 0:
        return  (face == 1 && p->x == p->z) ||
                (face == 0 && p->x == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level));
    case 1:
        return  (face == 2 && p->y == p->z) ||
                (face == 0 && p->x == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level));
    case 2:
        return  (face == 2 && p->x == p->z) ||
                (face == 3 && p->y == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level));
    case 3:
        return  (face == 1 && p->y == p->z) ||
                (face == 3 && p->y == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level));
    case 4:
        return 0; /*type 4 tets never touch a root boundary*/
    case 5:
        return 0; /*type 5 tets never touch a root boundary*/
    case 6:
        switch (face) {
        case 0:
            return p->x == p->z;
        case 1:
            return p->x == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level);
        case 2:
            return p->y == T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN(p->level);
        case 3:
            return p->y == p->z;
        case 4:
            return p->z == 0;
        default:
            SC_ABORT_NOT_REACHED();
        }
    case 7:
        return 0;   /*type 7 pyramids are never at the root boundary*/
    default:
        SC_ABORT_NOT_REACHED();
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

int
t8_dpyramid_custom_mod (uint64_t * id, t8_dpyramid_type_t type,
                        t8_linearidx_t pyra, t8_linearidx_t tet)
{
  t8_linearidx_t      test = 0, shift;
  T8_ASSERT (id >= 0);
  int                 remain = -1;
  do {;
    /* Iterate through the local-id. Get the current shift by the type of the
     * current element*/
    shift = t8_dpyramid_parenttype_Iloc_to_type[type][remain + 1] >= 6 ? pyra : tet;
    /*Add the shift to test*/
    test += shift;
    remain++;
  } while (test <= (*id));
  /*test is now larger than id, subtract last shift from test*/
  test -= shift;
  /*Compute the remaining ID*/
  (*id) -= test;
  T8_ASSERT(0 <= remain && remain < T8_DPYRAMID_CHILDREN);
  return remain;
}

void
t8_dpyramid_init_linear_id (t8_dpyramid_t * p, int level, uint64_t id)
{
  t8_dpyramid_type_t  type;
  t8_linearidx_t      local_index, p_sum1 = ((t8_linearidx_t)1)<< (3*level),
    p_sum2 = sc_intpow64u (6, level);
  t8_dpyramid_cube_id_t cid;
  int                 i;
  int                 offset_coords;
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= id && id <= 2 * p_sum1 - p_sum2);
  p->level = level;
  p->x = 0;
  p->y = 0;
  p->z = 0;
  type = 6;                     /*This is the type of the root pyramid */
  for (i = 1; i <= level; i++) {
    offset_coords = T8_DPYRAMID_MAXLEVEL - i;
    p_sum1 >>= 3;
    p_sum2 /= 6;
    // Thy types of the tetrahedron children of pyramid are always 0 or 3
    if (type == 0 || type == 3) {
      t8_dtet_init_linear_id_with_level ((t8_dtet_t *) p, id, i, level, type);
      return;
    }
    /* The local index depends on the alternating number of predecessors
     * caused by switching between pyramids and tetrahedrons, which have
     * a different number of children.*/
    local_index = t8_dpyramid_custom_mod (&id, type, 2 * p_sum1 - p_sum2,
                                          p_sum1);
    cid = t8_dpyramid_parenttype_Iloc_to_cid[type][local_index];
    T8_ASSERT (cid >= 0);
    /* Set the element in its cube*/
    p->x |= (cid % 2 == 1) ? 1 << offset_coords : 0;
    p->y |= (cid == 2 || cid == 3 || cid == 6
             || cid == 7) ? 1 << offset_coords : 0;
    p->z |= (cid > 3) ? 1 << offset_coords : 0;
    /*Compute the type*/
    type = t8_dpyramid_parenttype_Iloc_to_type[type][local_index];
    T8_ASSERT (type >= 0);
  }
  T8_ASSERT(id == 0);
  p->type = type;
}

t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t * p, int level)
{
  T8_ASSERT (level <= T8_DPYRAMID_MAXLEVEL);
  t8_linearidx_t      id = 0, pyra_shift, sum_1 = 1, sum_2 = 1, local_id;
  t8_dpyramid_t       parent, copy;
  int                 i, num_pyra, num_tet;
  t8_dpyramid_copy (p, &copy);
  copy.level = level;
  for (i = level; i > 0; i--) {
    /* Compute the number of pyramids with level maxlvl that are in a pyramid
     * of level i*/
    pyra_shift = (sum_1 << 1) - sum_2;

    /*Compute the parent and the local id of the current element*/
    t8_dpyramid_parent (&copy, &parent);
    local_id = t8_dpyramid_child_id_known_parent (&copy, &parent);

    /* Compute the number of predecessors within the parent that have the
     * shape of a pyramid or a tet*/
    if (t8_dpyramid_shape (&parent) == T8_ECLASS_TET) {
      /* If the parent is a tet, no predecessors are pyramids*/
      num_pyra = 0;
    }
    else {
      /* The number of pyramid-predecessors*/
      num_pyra =
        t8_dpyramid_parenttype_iloc_pyra_w_lower_id[parent.type -
                                                    6][local_id];
    }
    /* The number of tets is the local-id minus the number of pyramid-predecessors*/
    num_tet = local_id - num_pyra;
    /* The Id shifts by the number of predecessor elements*/
    id += num_pyra * pyra_shift + num_tet * sum_1;
    t8_dpyramid_copy (&parent, &copy);
    /* Update the shift*/
    sum_1 = sum_1 << 3;
    sum_2 *= 6;
  }
  T8_ASSERT(p->level >= 0);
  return id;
}

int
t8_dpyramid_face_neighbour(const t8_dpyramid_t *p, int face, t8_dpyramid_t * neigh)
{
    T8_ASSERT(0 <= face && face < T8_DPYRAMID_FACES);
    t8_dpyramid_coord_t len = T8_DPYRAMID_LEN(p->level);
    if(t8_dpyramid_shape(p) == T8_ECLASS_PYRAMID)
    /*pyramid touches tet or pyra*/
    {
        /*Compute the type of the neighbour*/
        if(face == 0 || face == 1){
            neigh->type =  3;
        }
        else if(face == 2 || face == 3){
            neigh->type = 0;
        }
        else
        {
            /*face == 4*/
            neigh->type = (p->type == 6)? 7: 6;
        }
        /*Compute the coords of the neighbour*/
        if(face == 0 || face == 3){
            neigh->x = p->x;
            neigh->y = p->y;
            neigh->z = p->z;
        }
        else if(face == 1){
            neigh->x = p->x + (p->type == 6)?0 : - len;
            neigh->y = p->y + (p->type == 6)?len : 0;
            neigh->z = p->z;
        }
        else if(face == 2){
            neigh->x = p->x+ (p->type == 6)?len : 0;
            neigh->y = p->y+ (p->type == 6)?0 : - len;
            neigh->z = p->z;
        }
        else{
            /*face == 4*/
            neigh->x = p->x;
            neigh->y = p->y;
            neigh->z = p->z + (p->type == 6) ? - len: len;
        }
        neigh->level = p->level;
        t8_debugf("[D] t: %i f: %i, nf: %i\n", p->type, face,
                  t8_dpyramid_type_face_to_nface[p->type - 6][face]);
        return t8_dpyramid_type_face_to_nface[p->type - 6][face];
    }
    else{
        if(p->type != 0 && p->type != 3){
            return t8_dtet_face_neighbour(p, face, neigh);
        }
        if(t8_dpyramid_tet_boundary(p, face)){
            /*tet touches pyra*/
            neigh->x = p->x;
            neigh->y = p->y;
            neigh->z = p->z;
            neigh->level = p->level;
            if(p->type == 0){
                switch(face){
                case 0:
                    neigh->x += len;
                    neigh->type = 7;
                    return 3;
                case 1:
                    neigh->type = 7;
                    return 2;
                case 2:
                    neigh->type = 6;
                    return 2;
                case 3:
                    neigh->y -=len;
                    neigh->type = 6;
                    return 3;
                default:
                    SC_ABORT_NOT_REACHED();
                }
            }
            else{
                /*p->type == 3*/
                switch (face) {
                case 0:
                    neigh->y += len;
                    neigh->type = 7;
                    return 1;
                case 1:
                    neigh->type = 7;
                    return 0;
                case 2:
                    neigh->type = 6;
                    return 0;
                case 3:
                    neigh->x -= len;
                    neigh->type = 6;
                    return 1;
                default:
                    SC_ABORT_NOT_REACHED();
                }
            }
        }
        else{
            /*tet touches tet*/
            return t8_dtet_face_neighbour(p, face, neigh);
        }
    }
}

int
t8_dpyramid_tet_boundary(const t8_dpyramid_t *p, int face)
{
    t8_dpyramid_t anc;
    int level = t8_dpyramid_is_inside_tet(p, p->level, &anc);
    t8_dpyramid_coord_t p_len = T8_DPYRAMID_LEN(p->level),
            a_len = T8_DPYRAMID_LEN(level), len_diff;
    T8_ASSERT(anc.type == 0 || anc.type == 3);
    len_diff = anc.z - p->level;
    if(p->level == 1){
        if(p->type == 0){
            return (p->x == 0 && (face != 1)) ||
                   (p->x != 0 && (face != 3));
        }
        else{
            /*p->type == 3*/
            T8_ASSERT(p->type == 3);
            return (p->y == 0 && (face != 1)) ||
                   (p->y != 0 && (face != 3));
        }
    }
    if(anc.type == 0){
        if(p->type == 0){
            switch(face){
            case 0:
                return  p->x == (anc.x + a_len - p_len);
            case 1:
                return  (p->x == anc.x + len_diff) &&
                        (p->y >= anc.y + len_diff) && (p->y <= anc.y + a_len - p_len);
            case 2:
                return  (p->y == (anc.y + len_diff)) &&
                        (p->x >= (anc.x + len_diff)) && (p->x <= (anc.x + a_len - p_len));
            case 3:
                return  (p->y == anc.y) &&
                        (p->x >= (anc.x + len_diff)) && (p->x <= (anc.x + a_len - p_len));
            default:
                SC_ABORT_NOT_REACHED();
            }
        }
        else{
            return 0;
        }

    }
    else {
        /*anc.type == 3*/
        if(p->type == 3){
            switch(face){
            case 0:
                return  (p->y == anc.y + a_len - p_len) &&
                        (p->x >= anc.x) && (p->x <= anc.x + len_diff);
            case 1:
                return  (p->y == anc.y + len_diff) &&
                        (p->x >= anc.x) && (p->x <= anc.x + len_diff);
            case 2:
                return  (p->x == anc.x + len_diff) &&
                        (p->y >= anc.y) && (p->y <= anc.y + len_diff);
            case 3:
                return  (p->x == anc.x) &&
                        (p->y >= anc.y + len_diff) && (p->y <= anc.y + a_len - p_len);
            default:
                SC_ABORT_NOT_REACHED();
            }
        }
        else{
            return 0;
        }
    }
}

int
t8_dpyramid_is_inside_root(t8_dpyramid_t * p)
{
    if(p->level == 0){
        return p->type == T8_DPYRAMID_ROOT_TPYE && p->x == 0 && p->y == 0 && p->z == 0;
    }
    return (0 <= p->z) && (p->z < T8_DPYRAMID_ROOT_LEN) &&
                    (p->x >= p->z) && (p->x < T8_DPYRAMID_ROOT_LEN) &&
                    (p->y >= p->z) && (p->y < T8_DPYRAMID_ROOT_LEN);
}

int
t8_dpyramid_face_neighbor_inside (const t8_dpyramid_t *p,
                                            t8_dpyramid_t * neigh,
                                            int face, int *neigh_face)
{
    *neigh_face = t8_dpyramid_face_neighbour(p,  face, neigh);
    return t8_dpyramid_is_inside_root(neigh);
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                              int level)
{
  t8_linearidx_t      id;
  T8_ASSERT(level >= p->level);
  T8_ASSERT(0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The first descendant of a pyramid has the same anchor coords, but another level */
    t8_dpyramid_copy(p, desc);
    desc->level = level;
  }
  else {
    id = t8_dpyramid_linear_id (p, level);
    t8_dpyramid_init_linear_id (desc, level, id);
  }
  T8_ASSERT(p->x <= desc->x && p->y <= desc->y && p->z <= desc->z);
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t * p, t8_dpyramid_t * desc,
                             int level)
{
  t8_linearidx_t      id = 0, t_id;
  int                 exponent;
  T8_ASSERT(level >= p->level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_copy (p, desc);
    desc->level = level;
    int                 coord_offset =
      T8_DPYRAMID_LEN (p->level) - T8_DPYRAMID_LEN (level);
    desc->x |= coord_offset;
    desc->y |= coord_offset;
    desc->z |= coord_offset;
  }
  else {
    t_id = t8_dpyramid_linear_id (p, level);
    exponent = level - p->level;
    id = (((t8_linearidx_t) 1) << 3 * exponent) - 1;
    id += t_id;
    t8_dpyramid_init_linear_id (desc, level, id);
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
    t8_debugf("[D] nc-in: %i %i %i %i %i\n", p->x, p->y, p->z, p->type, p->level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_CHILDREN;
  }
  else {
    return T8_DPYRAMID_CHILDREN;
  }
}

int
t8_dpyramid_num_siblings(const t8_dpyramid_t * p)
{
    t8_dpyramid_t parent;
    t8_dpyramid_parent(p, &parent);
    return t8_dpyramid_num_children(&parent);
}

int
t8_dpyramid_child_id_unknown_parent (const t8_dpyramid_t * p,
                                     t8_dpyramid_t * parent)
{
  T8_ASSERT (p->level > 0);
  t8_dpyramid_parent (p, parent);
  return t8_dpyramid_child_id_known_parent (p, parent);

}

int
t8_dpyramid_child_id_known_parent (const t8_dpyramid_t * p,
                                   t8_dpyramid_t * parent)
{
  t8_dpyramid_cube_id_t cid = compute_cubeid (p, p->level);
  if (t8_dpyramid_shape (parent) == T8_ECLASS_PYRAMID) {
    T8_ASSERT (t8_dpyramid_type_cid_to_Iloc[p->type][cid] >= 0);
    return t8_dpyramid_type_cid_to_Iloc[p->type][cid];
  }
  else {
    return t8_dtet_child_id ((const t8_dtet_t *) p);
  }
}

int
t8_dpyramid_child_id (const t8_dpyramid_t * p)
{
  T8_ASSERT (p->level > 0);
  t8_dpyramid_t       parent;
  return t8_dpyramid_child_id_unknown_parent (p, &parent);
}

void
t8_dpyramid_child (const t8_dpyramid_t * elem, int child_id,
                   t8_dpyramid_t * child)
{

  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_coord_t h;
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);

  if (t8_dpyramid_shape (elem) == T8_ECLASS_TET) {
    t8_dtet_child ((t8_dtet_t *) elem, child_id, (t8_dtet_t *) child);
  }
  else {
    cube_id = t8_dpyramid_parenttype_Iloc_to_cid[elem->type][child_id];
    T8_ASSERT (cube_id >= 0);;
    child->level = elem->level + 1;
    h = T8_DPYRAMID_LEN (child->level);
    child->x = elem->x + ((cube_id & 0x01) ? h : 0);
    child->y = elem->y + ((cube_id & 0x02) ? h : 0);
    child->z = elem->z + ((cube_id & 0x04) ? h : 0);
    child->type = t8_dpyramid_parenttype_Iloc_to_type[elem->type][child_id];

  }
  T8_ASSERT (child->type >= 0);
}

void
t8_dpyramid_children (const t8_dpyramid_t * p, int length, t8_dpyramid_t ** c)
{
  int                 i, num_children;
  num_children = t8_dpyramid_num_children (p);
  for (i = num_children - 1; i >= 0; i--) {
    t8_dpyramid_child (p, i, c[i]);
  }

}

/* Check, if a pyramid in the shape of a tet lies inside a tetrahedron
 * The i first bits give the anchorcoordinate for a possible ancestor of level i
 * for p. */
int
t8_dpyramid_is_inside_tet (const t8_dpyramid_t * p, int level, t8_dpyramid_t *anc)
{
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (p->type == 0 || p->type == 3);
  int                 i;
  /*the tet is initialized, the ancestor will be computed */
  t8_dtet_t           tet, ancestor;
  tet.x = 0;
  tet.y = 0;
  tet.z = 0;
  for (i = 1; i < level; i++) {
    /*Update the coordinate of tet to i first bits */
    tet.x = tet.x | (p->x & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.y = tet.y | (p->y & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.z = tet.z | (p->z & (1 << (T8_DPYRAMID_MAXLEVEL - i)));
    tet.level = i;
    /*Compute the tet-ancestor */
    t8_dtet_ancestor ((const t8_dtet_t *) p, i, &ancestor);
    /*Only compare, if the ancestor has type 0 or 3 */
    if (ancestor.type == 0 || ancestor.type == 3) {
      /*Set tet.type to the type of the possible ancestor */
      tet.type = ancestor.type;
      /*Compare */
      if (t8_dtet_is_equal (&tet, &ancestor)) {
          if(anc != NULL){
              t8_dpyramid_copy(&ancestor, anc);
          }
        return i;
      }
    }
  }
  /*No matching tet-ancestor was found, the parent is a pyramid */
  return 0;
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
  T8_ASSERT (p->level > 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);

  t8_dpyramid_coord_t h = T8_DPYRAMID_LEN (p->level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The parent of a pyramid is a pyramid, maybe of different type */

    int                 cube_id = compute_cubeid (p, p->level);

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
  else if (t8_dpyramid_is_inside_tet (p, p->level, NULL) != 0) {
    /*Pyramid- / tetparent detection */
    /*If a tetrahedron does not reach a "significant point" its parent is a tet */
    /*Tetcase */ ;
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
  T8_ASSERT (parent->level >= 0);
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
t8_dpyramid_successor_recursion (const t8_dpyramid_t * elem,
                                 t8_dpyramid_t * succ, t8_dpyramid_t * parent,
                                 int level)
{
  int                 child_id, num_children;
  t8_dpyramid_copy (elem, succ);
  T8_ASSERT (1 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  succ->level = level;
  T8_ASSERT (succ->type >= 0);
  child_id = t8_dpyramid_child_id_unknown_parent (elem, parent);
  /*Compute number of children */
  num_children = t8_dpyramid_num_children (parent);
  T8_ASSERT (0 <= child_id && child_id < num_children);
  if (child_id == num_children - 1) {
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion (succ, succ, parent, level - 1);
    succ->level = level;
    /* bits auf level auf child 0 setzen */
    succ->x =(succ->x >> (T8_DPYRAMID_MAXLEVEL - level + 1))
            << (T8_DPYRAMID_MAXLEVEL - level + 1);
    succ->y =(succ->y >> (T8_DPYRAMID_MAXLEVEL - level + 1))
            << (T8_DPYRAMID_MAXLEVEL - level + 1);
    succ->z =(succ->z >> (T8_DPYRAMID_MAXLEVEL - level + 1))
            << (T8_DPYRAMID_MAXLEVEL -level + 1);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_dpyramid_child (parent, child_id + 1, succ);
    succ->level = level;
  }
}

void
t8_dpyramid_successor (const t8_dpyramid_t * elem, t8_dpyramid_t * succ,
                       int level)
{
  t8_dpyramid_t       parent;
  t8_dpyramid_successor_recursion (elem, succ, &parent, level);
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t * p, int vertex, int coords[])
{
  t8_dpyramid_coord_t h;
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);

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
