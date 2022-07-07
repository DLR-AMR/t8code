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
#include "t8_dpyramid_connectivity.h"
#include <sc_functions.h>
#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>

typedef int8_t t8_dpyramid_cube_id_t;

static t8_dpyramid_cube_id_t
compute_cubeid (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_cube_id_t cube_id = 0;
  t8_dpyramid_coord_t h;

  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }

  cube_id |= ((p->x & h) ? 0x01 : 0);
  cube_id |= ((p->y & h) ? 0x02 : 0);
  cube_id |= ((p->z & h) ? 0x04 : 0);

  return cube_id;
}

/**
 * Starting from a level where the type of \a p is known compute the type
 * of \a p at level \a level
 * 
 * \param [in]  p           Input pyramid
 * \param [in]  level       The level at which we want to know the type of \a p. Must be lower than the level of \a p
 * \param [in]  known_type  The type of \a p at \a known_level
 * \param [in]  known_level The level where we already know the type of \a p
 * \return      t8_dpyramid_type_t The type of \a p at level \a level.
 */
static t8_dpyramid_type_t
compute_type_ext (const t8_dpyramid_t *p, const int level,
                  const t8_dpyramid_type_t known_type, const int known_level)
{
  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_type_t  type = known_type;
  int                 i;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (known_level <= p->level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /*Type of the root pyra */
    return T8_DPYRAMID_ROOT_TPYE;
  }
  for (i = known_level; i > level; i--) {
    cube_id = compute_cubeid (p, i);
    type = t8_dpyramid_cid_type_to_parenttype[cube_id][type];
  }
  return type;
}

t8_dpyramid_type_t
compute_type (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  return compute_type_ext (p, level, p->type, p->level);
}

/* sets the shift last bits of every coordinate to zero*/
static void
pyramid_cut_coords (t8_dpyramid_t *p, const int shift)
{
  T8_ASSERT (0 <= shift && shift <= T8_DPYRAMID_MAXLEVEL);
  p->x = (p->x >> shift) << shift;
  p->y = (p->y >> shift) << shift;
  p->z = (p->z >> shift) << shift;
}

int
t8_dpyramid_ancestor_id (const t8_dpyramid_t *p, const int level)
{
  t8_linearidx_t      id;
  t8_dpyramid_t       helper;
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /*Compute the id of p at the given level */
  id = t8_dpyramid_linear_id (p, level);
  /*Compute the child id of the ancestor of p at the given level */
  t8_dpyramid_init_linear_id (&helper, level, id);
  return t8_dpyramid_child_id (&helper);
}

int
t8_dpyramid_is_family (const t8_dpyramid_t **fam)
{

  const int           level = fam[0]->level;
  int                 i, type_of_first;
  t8_dpyramid_coord_t inc = T8_DPYRAMID_LEN (level), x_inc, y_inc;
  if (t8_dpyramid_shape (fam[0]) == T8_ECLASS_TET) {
    return t8_dtet_is_familypv (fam);
  }
  else {
    if (level == 0) {
      return 0;
    }
    /*The type of parent is the type of the first child in z-curve-order */
    type_of_first = fam[0]->type;
    T8_ASSERT (type_of_first == T8_DPYRAMID_FIRST_TYPE
               || type_of_first == T8_DPYRAMID_SECOND_TYPE);
    for (i = 1; i < T8_DPYRAMID_CHILDREN; i++) {
      /*All elements must have the same level to be a family */
      if (fam[i]->level != level) {
        return 0;
      }
      /*Check if every family-member has the correct type */
      if (t8_dpyramid_parenttype_Iloc_to_type[type_of_first][i] !=
          fam[i]->type) {
        return 0;
      }
    }

    x_inc = fam[0]->x + inc;
    y_inc = fam[0]->y + inc;
    /*Check the coordinates of the anchor-coordinate */
    if (type_of_first == T8_DPYRAMID_FIRST_TYPE) {
      return fam[0]->z == fam[1]->z && fam[0]->z == fam[2]->z
        && fam[0]->z == fam[3]->z && fam[0]->z == fam[4]->z
        && fam[0]->z == fam[5]->z && fam[0]->z == fam[6]->z
        && fam[0]->z == fam[7]->z && fam[0]->z == fam[8]->z
        && fam[0]->z == (fam[9]->z - inc) && fam[0]->x == fam[3]->x
        && fam[0]->x == fam[4]->x && x_inc == fam[1]->x && x_inc == fam[2]->x
        && x_inc == fam[5]->x && x_inc == fam[6]->x && x_inc == fam[7]->x
        && x_inc == fam[8]->x && x_inc == fam[9]->x && fam[0]->y == fam[1]->y
        && fam[0]->y == fam[2]->y && y_inc == fam[3]->y && y_inc == fam[4]->y
        && y_inc == fam[5]->y && y_inc == fam[6]->y && y_inc == fam[7]->y
        && y_inc == fam[8]->y && y_inc == fam[9]->y;
    }
    else {
      return fam[1]->z == fam[0]->z + inc && fam[1]->z == fam[2]->z
        && fam[1]->z == fam[3]->z && fam[1]->z == fam[4]->z
        && fam[1]->z == fam[5]->z && fam[1]->z == fam[6]->z
        && fam[1]->z == fam[7]->z && fam[1]->z == fam[8]->z
        && fam[1]->z == fam[9]->z && fam[0]->x == fam[1]->x
        && fam[0]->x == fam[2]->x && fam[0]->x == fam[3]->x
        && fam[0]->x == fam[4]->x && fam[0]->x == fam[7]->x
        && fam[0]->x == fam[8]->x && x_inc == fam[5]->x && x_inc == fam[6]->x
        && x_inc == fam[9]->x && fam[0]->y == fam[1]->y
        && fam[0]->y == fam[2]->y && fam[0]->y == fam[3]->y
        && fam[0]->y == fam[4]->y && fam[0]->y == fam[5]->y
        && fam[0]->y == fam[6]->y && y_inc == fam[7]->y && y_inc == fam[8]->y
        && y_inc == fam[9]->y;
    }
  }
}

int
t8_dpyramid_is_root_boundary (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= face && face <= T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_coord_t coord_touching_root = T8_DPYRAMID_ROOT_LEN -
    T8_DPYRAMID_LEN (p->level);
  if (!t8_dpyramid_is_inside_root (p)) {
    return 0;
  }
  switch (p->type) {
  case 0:
    return (face == 1 && p->x == p->z) ||
      (face == 0 && p->x == coord_touching_root);
  case 1:
    return (face == 2 && p->y == p->z) ||
      (face == 0 && p->x == coord_touching_root);
  case 2:
    return (face == 2 && p->x == p->z) ||
      (face == 0 && p->y == coord_touching_root);
  case 3:
    return (face == 1 && p->y == p->z) ||
      (face == 0 && p->y == coord_touching_root);
  case 4:
    return 0;                   /*type 4 tets never touch a root boundary */
  case 5:
    return 0;                   /*type 5 tets never touch a root boundary */
  case T8_DPYRAMID_FIRST_TYPE:
    switch (face) {
    case 0:
      return p->x == p->z;
    case 1:
      return p->x == coord_touching_root;
    case 2:
      return p->y == p->z;
    case 3:
      return p->y == coord_touching_root;
    case 4:
      return p->z == 0;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  case T8_DPYRAMID_SECOND_TYPE:
    return 0;                   /*type 7 pyramids are never at the root boundary */
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

/*Copies a pyramid from source to dest*/
void
t8_dpyramid_copy (const t8_dpyramid_t *source, t8_dpyramid_t *dest)
{
  T8_ASSERT (source != NULL && dest != NULL);
  if (source == dest) {
    return;
  }
  memcpy (dest, source, sizeof (t8_dpyramid_t));
}

int
t8_dpyramid_compare (const t8_dpyramid_t *p1, const t8_dpyramid_t *p2)
{
  int                 maxlvl;
  t8_linearidx_t      id1, id2;
  T8_ASSERT (p1->x >= 0 && p1->y >= 0 && p1->z >= 0 &&
             p1->level >= 0 && p1->level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p2->x >= 0 && p2->y >= 0 && p2->z >= 0 &&
             p2->level >= 0 && p2->level <= T8_DPYRAMID_MAXLEVEL);
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
t8_dpyramid_get_level (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  return p->level;
}

/**
 * Computes the local index of a pyramid and updates its current global index. 
 * for further computation of init_linear_id. For each sibling that is a predecessor
 * the number of pyramids or tets on the level used in init_linear_id is substracted from
 * the id. 
 * 
 * \param[in, out] id       The current id that will be updated
 * \param[in] type          The type of the current pyramid
 * \param[in] pyra          Number of pyramids to shift 
 * \param[in] tet           Number of tets to shift
 * \return int              The local-id of the child
 */
static int
t8_dpyramid_update_index (t8_linearidx_t * id, const t8_dpyramid_type_t type,
                          const t8_linearidx_t pyra, const t8_linearidx_t tet)
{
  t8_linearidx_t      test = 0;
  t8_linearidx_t      shift;
  T8_ASSERT (id != NULL);
  T8_ASSERT (*id >= 0);
  int                 remain = -1;
  do {
    /* Iterate through the local-id. Get the current shift by the type of the
     * current element*/
    shift =
      t8_dpyramid_parenttype_Iloc_to_type[type][remain + 1] >=
      T8_DPYRAMID_FIRST_TYPE ? pyra : tet;
    /*Add the shift to test */
    test += shift;
    remain++;
  } while (test <= (*id));
  /*test is now larger than id, subtract last shift from test */
  test -= shift;
  /*Compute the remaining ID */
  (*id) -= test;
  T8_ASSERT (0 <= remain && remain < T8_DPYRAMID_CHILDREN);
  return remain;
}

void
t8_dpyramid_init_linear_id (t8_dpyramid_t *p, const int level,
                            t8_linearidx_t id)
{
  t8_dpyramid_type_t  type;
  t8_linearidx_t      local_index;
  t8_linearidx_t      p_sum1 = ((t8_linearidx_t) 1) << (3 * level);
  t8_linearidx_t      p_sum2 = sc_intpow64u (6, level);
  t8_dpyramid_cube_id_t cube_id;
  int                 i, offset_expo;

  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= id && id <= 2 * p_sum1 - p_sum2);

  p->level = level;
  p->x = 0;
  p->y = 0;
  p->z = 0;
  type = T8_DPYRAMID_ROOT_TPYE; /*This is the type of the root pyramid */
  for (i = 1; i <= level; i++) {
    offset_expo = T8_DPYRAMID_MAXLEVEL - i;
    p_sum1 >>= 3;
    p_sum2 /= 6;
    // Thy types of the tetrahedron children of pyramid are always 0 or 3
    if (type == 0 || type == 3) {
      /* Ensure that a valid tet is used in init_linear_id_with_level */
      p->type = type;
      p->level = i;
#if T8_ENABLE_DEBUG
      ((t8_dtet_t *) p)->eclass_int8 = T8_ECLASS_TET;
#endif
      t8_dtet_init_linear_id_with_level (p, id, i, level, type);
      return;
    }
    /* The local index depends on the alternating number of predecessors
     * caused by switching between pyramids and tetrahedrons, which have
     * a different number of children.*/
    local_index = t8_dpyramid_update_index (&id, type, 2 * p_sum1 - p_sum2,
                                            p_sum1);
    cube_id = t8_dpyramid_parenttype_Iloc_to_cid[type][local_index];
    T8_ASSERT (cube_id >= 0);
    /* Set the element in its cube */
    p->x |= (cube_id % 2 == 1) ? 1 << offset_expo : 0;
    p->y |= (cube_id == 2 || cube_id == 3 || cube_id == 6 || cube_id == 7) ?
      1 << offset_expo : 0;
    p->z |= (cube_id > 3) ? 1 << offset_expo : 0;
    /*Compute the type */
    type = t8_dpyramid_parenttype_Iloc_to_type[type][local_index];
    T8_ASSERT (type >= 0);
  }
  T8_ASSERT (id == 0);
  p->type = type;
}

/* Compute the type of a pyramid at a given level. Starting from its own level,
 * we iterate over the levels and compute the type of this level. If p is a tetrahedron,
 * we compute it in a tetrahedral fashion up unto the last level where p is a tet and
 * continue in a pyramidal fashion*/
int
t8_dpyramid_set_type_at_level (const t8_dpyramid_t *p, const int level)
{
  int                 level_iter;
  int                 start;
  t8_dpyramid_type_t  type = p->type;
  t8_dpyramid_cube_id_t cube_id;
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);

  if (level >= p->level) {
    return p->type;
  }
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*all ancs are pyramids */
    for (level_iter = p->level; level_iter > level; level_iter--) {
      cube_id = compute_cubeid (p, level_iter);
      type =
        t8_dpyramid_type_cid_to_parenttype[type -
                                           T8_DPYRAMID_FIRST_TYPE][cube_id];
    }
    return type;
  }
  else {
    int                 max_tet_lvl;
    t8_dpyramid_t       helper;
    for (level_iter = p->level; level_iter > level && type != 0 && type != 3;
         level_iter--) {
      /*all anc fulfilling the above condition are not pyramids */
      cube_id = compute_cubeid (p, level_iter);
      type = t8_dtet_cid_type_to_parenttype[cube_id][type];
    }
    if (level_iter == level) {
      return type;
    }
    if (level_iter != p->level) {
      start = level_iter;
    }
    else {
      start = p->level;
    }
    t8_dpyramid_copy (p, &helper);
    helper.level = start;
    helper.type = type;
    /*Compute the last tetrahedral level */
    max_tet_lvl = t8_dpyramid_is_inside_tet (&helper, helper.level, NULL);
    if (level >= max_tet_lvl && max_tet_lvl != 0) {
      for (level_iter = start; level_iter > level; level_iter--) {
        /*Up to this level all anc are tets */
        cube_id = compute_cubeid (p, level_iter);
        type = t8_dtet_cid_type_to_parenttype[cube_id][type];
      }
    }
    else {
      if (max_tet_lvl == 0) {
        /*parent is a pyra */
        max_tet_lvl = helper.level;
      }
      for (level_iter = helper.level; level_iter > max_tet_lvl; level_iter--) {
        /*from current level on until max-tet-lvl all ancs are tet */
        cube_id = compute_cubeid (p, level_iter);
        type = t8_dtet_cid_type_to_parenttype[cube_id][type];
      }
      cube_id = compute_cubeid (p, max_tet_lvl);
      type = t8_dtet_type_cid_to_pyramid_parenttype[type][cube_id];
      /*compute type of pyramidal tet parent */
      T8_ASSERT (type == T8_DPYRAMID_FIRST_TYPE
                 || type == T8_DPYRAMID_SECOND_TYPE);
      for (level_iter = max_tet_lvl - 1; level_iter > level; level_iter--) {
        /*remaining ancs are pyramids */
        cube_id = compute_cubeid (p, level_iter);
        type =
          t8_dpyramid_type_cid_to_parenttype[type -
                                             T8_DPYRAMID_FIRST_TYPE][cube_id];
      }
    }
    return type;
  }
}

t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  t8_linearidx_t      id = 0, pyra_shift, sum_1 = 1, sum_2 = 1, local_id;
  t8_dpyramid_t       parent, copy;
  int                 i, num_pyra, num_tet;

  t8_dpyramid_copy (p, &copy);
  copy.type = t8_dpyramid_set_type_at_level (p, level);
  copy.level = level;
  pyramid_cut_coords (&copy, T8_DPYRAMID_MAXLEVEL - level);

  for (i = level; i > 0; i--) {
    /* Compute the number of pyramids with level maxlvl that are in a pyramid
     * of level i*/
    pyra_shift = (sum_1 << 1) - sum_2;

    /*Compute the parent and the local id of the current element */
    t8_dpyramid_parent (&copy, &parent);
    local_id = t8_dpyramid_child_id_known_parent (&copy, &parent);

    /* Compute the number of predecessors within the parent that have the
     * shape of a pyramid or a tet*/
    if (t8_dpyramid_shape (&parent) == T8_ECLASS_TET) {
      /* If the parent is a tet, no predecessors are pyramids */
      num_pyra = 0;
    }
    else {
      /* The number of pyramid-predecessors */
      num_pyra =
        t8_dpyramid_parenttype_iloc_pyra_w_lower_id[parent.type -
                                                    T8_DPYRAMID_FIRST_TYPE]
        [local_id];
    }
    /* The number of tets is the local-id minus the number of pyramid-predecessors */
    num_tet = local_id - num_pyra;
    /* The Id shifts by the number of predecessor elements */
    id += num_pyra * pyra_shift + num_tet * sum_1;
    t8_dpyramid_copy (&parent, &copy);
    /* Update the shift */
    sum_1 = sum_1 << 3;
    sum_2 *= 6;
  }
  T8_ASSERT (p->level >= 0);
  return id;
}

int
t8_dpyramid_face_neighbour (const t8_dpyramid_t *p, const int face,
                            t8_dpyramid_t *neigh)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->level);
  neigh->x = p->x;
  neigh->y = p->y;
  neigh->z = p->z;
  neigh->level = p->level;
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*pyramid touches tet or pyra */
    /*Compute the type of the neighbour */
    if (face == 0 || face == 1) {
      neigh->type = 3;
    }
    else if (face == 2 || face == 3) {
      neigh->type = 0;
    }
    else {
      /*face == 4 */
      neigh->type =
        ((p->type ==
          T8_DPYRAMID_FIRST_TYPE) ? T8_DPYRAMID_SECOND_TYPE :
         T8_DPYRAMID_FIRST_TYPE);
    }
    /*Compute the coords of the neighbour */
    /*Do nothing for face == 0 || face == 2 */
    if (face == 1) {
      neigh->x += ((p->type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
      neigh->y += ((p->type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
    }
    else if (face == 3) {
      neigh->x += ((p->type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
      neigh->y += ((p->type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
    }
    else if (face == 4) {
      neigh->z += ((p->type == T8_DPYRAMID_FIRST_TYPE) ? -length : length);
    }

    return t8_dpyramid_type_face_to_nface[p->type -
                                          T8_DPYRAMID_FIRST_TYPE][face];
  }
  else {
    /*Check if the neighbor is a tet, or a pyra */
    if (p->type != 0 && p->type != 3) {
      /*tets of these types never have a pyra-neighbor */
      return t8_dtet_face_neighbour (p, face, neigh);
    }
    if (t8_dpyramid_tet_boundary (p, face)) {
      /*tet touches pyra, compute the pyra */
      if (p->type == 0) {
        switch (face) {
        case 0:
          neigh->x += length;
          neigh->type = T8_DPYRAMID_SECOND_TYPE;
          return 3;
        case 1:
          neigh->type = T8_DPYRAMID_SECOND_TYPE;
          return 2;
        case 2:
          neigh->type = T8_DPYRAMID_FIRST_TYPE;
          return 2;
        case 3:
          neigh->y -= length;
          neigh->type = T8_DPYRAMID_FIRST_TYPE;
          return 3;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {
        /*p->type == 3 */
        switch (face) {
        case 0:
          neigh->y += length;
          neigh->type = T8_DPYRAMID_SECOND_TYPE;
          return 1;
        case 1:
          neigh->type = T8_DPYRAMID_SECOND_TYPE;
          return 0;
        case 2:
          neigh->type = T8_DPYRAMID_FIRST_TYPE;
          return 0;
        case 3:
          neigh->x -= length;
          neigh->type = T8_DPYRAMID_FIRST_TYPE;
          return 1;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
    }
    else {
      /*tet touches tet */
      return t8_dtet_face_neighbour (p, face, neigh);
    }
  }
}

int
t8_dpyramid_face_parent_face (const t8_dpyramid_t *elem, const int face)
{
  t8_dpyramid_t       parent;
  int                 child_id;
  T8_ASSERT (0 <= elem->level && elem->level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  /*parent is a pyramid */
  if (t8_dpyramid_shape (elem) == T8_ECLASS_PYRAMID) {
    child_id = t8_dpyramid_child_id_unknown_parent (elem, &parent);
    int                 i;
    /*If the pyramid is one of the children in the array, its face-num and the face-num
     * of the parent are the same*/
    for (i = 0; i < 4; i++) {
      if (t8_dpyramid_type_face_to_children_at_face
          [parent.type - T8_DPYRAMID_FIRST_TYPE][face][i]
          == child_id) {
        return face;
      }
    }
    /*No matching face */
    return -1;
  }
  else {
    child_id = t8_dpyramid_child_id_unknown_parent (elem, &parent);
    /*Parent is also a tet, we can call the tet-routine */
    if (t8_dpyramid_shape (&parent) == T8_ECLASS_TET) {
      return t8_dtet_face_parent_face (elem, face);
    }
    /*tet with a pyramid-parent */
    else {
      /*Only tets of type 0 or 3 have a pyra-parent. Parent can have type 6 or 7 */
      if (elem->type == 0 && parent.type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 3 && face == 1) {
          return 0;
        }
        else if (child_id == 5 && face == 0) {
          return 1;
        }
        else
          return -1;
      }
      else if (elem->type == 3 && parent.type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 1 && face == 1) {
          return 2;
        }
        else if (child_id == 6 && face == 0) {
          return 3;
        }
        else
          return -1;
      }
      else if (elem->type == 0 && parent.type == T8_DPYRAMID_SECOND_TYPE) {
        if (child_id == 1 && face == 3) {
          return 1;
        }
        else if (child_id == 7 && face == 2) {
          return 0;
        }
        else
          return -1;
      }
      else if (elem->type == 3 && parent.type == T8_DPYRAMID_SECOND_TYPE) {
        if (child_id == 2 && face == 3) {
          return 3;
        }
        else if (child_id == 5 && face == 2) {
          return 2;
        }
        else
          return -1;
      }
      return -1;
    }
  }
}

int
t8_dpyramid_tet_pyra_face_connection (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (p->type == 0 || p->type == 3);
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /* Depending on its cube-id at its level and its type,
   * 3 faces of a tet connect to a pyramid, one is connecting to a tet*/
  t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->level);
  if ((cube_id == 2 && face != 1) || (cube_id == 6 && face != 2)) {
    return p->type == 0;
  }
  else if ((cube_id == 1 && face != 1) || (cube_id == 5 && face != 2)) {
    return p->type == 3;
  }
  else if (cube_id == 3) {
    return face != 0;
  }
  else if (cube_id == 4) {
    return face != 3;
  }
  else {
    return 0;
  }
}

int
t8_dpyramid_tet_boundary (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_t       anc;
  const int           level = t8_dpyramid_is_inside_tet (p, p->level, &anc);
  int                 bey_id, i, type_temp, valid_touch;
  t8_dpyramid_cube_id_t cube_id;
  if (level == 0) {
    /*Check if the face is connecting to a tet or to a pyra */
    return t8_dpyramid_tet_pyra_face_connection (p, face);
  }
  /*Check if anc-face is connecting to a tet or to a pyra */
  valid_touch = t8_dpyramid_tet_pyra_face_connection (&anc, face);
  if (valid_touch) {
    type_temp = p->type;
    /* If so, check if the tet always lies in the corner of of its parent at this face.
     * Otherwise, the neighbor is a tet*/
    for (i = p->level; i > anc.level; i--) {
      cube_id = compute_cubeid (p, i);
      bey_id = t8_dtet_type_cid_to_beyid[type_temp][cube_id];
      if (t8_dpyramid_face_childid_to_is_inside[face][bey_id] == -1) {
        return 0;
      }
      type_temp = t8_dtet_cid_type_to_parenttype[cube_id][type_temp];
    }
  }
  /*Return if anc-face is connecting to a tet */
  return valid_touch;
}

int
t8_dpyramid_tree_face (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_is_root_boundary (p, face)) {
    if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
      /*If p is a pyramid and touches the boundary, the face-number is the same */
      return face;
    }
    else {
      /*p is a tet and in some occasions p shares a face with its tree */
      if (face == 0 && (p->type == 3 || p->type == 2)) {
        return 3;
      }
      else if (face == 0 && (p->type == 0 || p->type == 1)) {
        return 1;
      }
      else if ((face == 1 && p->type == 3) || (face == 2 && p->type == 1)) {
        return 2;
      }
      else if ((face == 1 && p->type == 0) || (face == 2 && p->type == 2)) {
        return 0;
      }
      else
        return -1;
    }
  }
  return -1;
}

int
t8_dpyramid_is_inside_root (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /*Check, if p is root pyramid */
  if (p->level == 0) {
    return p->type == T8_DPYRAMID_ROOT_TPYE && p->x == 0 && p->y == 0
      && p->z == 0;
  }
  /*Check, if all coordinates are in the limit set by the length of root */
  if ((0 <= p->z) && (p->z < T8_DPYRAMID_ROOT_LEN) &&
      (p->x >= p->z) && (p->x < T8_DPYRAMID_ROOT_LEN) &&
      (p->y >= p->z) && (p->y < T8_DPYRAMID_ROOT_LEN)) {
    if ((p->x == p->z && (p->type == 3 || p->type == 5)) ||
        (p->y == p->z && (p->type == 0 || p->type == 4))) {
      return 0;
    }
    else {
      return 1;
    }
  }
  else {
    return 0;
  }
}

int
t8_dpyramid_face_neighbor_inside (const t8_dpyramid_t *p,
                                  t8_dpyramid_t *neigh,
                                  const int face, int *neigh_face)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /*Compute the face-neighbor, then check if it is inside root */
  *neigh_face = t8_dpyramid_face_neighbour (p, face, neigh);
  return t8_dpyramid_is_inside_root (neigh);
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc,
                              const int level)
{
  t8_linearidx_t      id;
  T8_ASSERT (level >= p->level);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The first descendant of a pyramid has the same anchor coords, but another level */
    t8_dpyramid_copy (p, desc);
    desc->level = level;
  }
  else {
    id = t8_dpyramid_linear_id (p, level);
    t8_dpyramid_init_linear_id (desc, level, id);
  }
  T8_ASSERT (p->x <= desc->x && p->y <= desc->y && p->z <= desc->z);
}

/*compute the descendant at a corner of a tet up to a given level
  You can not use the tetrahedron algorithm here, because it depends on the linear-id
  computation of tet, which if different to the linear id for pyras. */
void
t8_dpyramid_corner_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *d,
                               const int corner, const int level)
{
  int                 child_id, i;
  T8_ASSERT (p->level <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (0 <= corner && corner < T8_DTET_CORNERS);
  if (corner == 0) {
    t8_dpyramid_first_descendant (p, d, level);
  }
  else if (corner == 1 || corner == 2) {
    /*The child at this corner is iterativle the child at child-id up to the
       given level */
    child_id = t8_dtet_parenttype_beyid_to_Iloc[p->type][corner];
    t8_dtet_t           tmp;
    t8_dpyramid_copy (p, &tmp);
    for (i = p->level; i < level; i++) {
      t8_dtet_child (&tmp, child_id, d);
      t8_dpyramid_copy (d, &tmp);
    }
  }
  else {
    /*corner == 3 */
    t8_dpyramid_last_descendant (p, d, level);
  }
}

void
t8_dpyramid_first_descendant_face (const t8_dpyramid_t *p, const int face,
                                   t8_dpyramid_t *first_desc, const int level)
{
  int                 corner;
  t8_dpyramid_coord_t off_set = T8_DPYRAMID_LEN (p->level) -
    T8_DPYRAMID_LEN (level);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->level <= level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    corner = t8_dtet_face_corner[face][0];
    t8_dpyramid_corner_descendant (p, first_desc, corner, level);
  }
  else if (p->level == T8_DPYRAMID_MAXLEVEL) {
    t8_dpyramid_copy (p, first_desc);
  }
  else {
    /*Shift the coordinate of p to compute the descendant */
    if ((p->type == T8_DPYRAMID_FIRST_TYPE
         && (face == 0 || face == 2 || face == 4))
        || (p->type == T8_DPYRAMID_SECOND_TYPE && face != 4)) {
      /*No shifting is needed, fd is the first child with given level */
      t8_dpyramid_child (p, 0, first_desc);
      first_desc->level = level;
    }
    else {
      /*shifting is needed, the fd does not have the same coord as p */
      t8_dpyramid_copy (p, first_desc);
      first_desc->level = level;
      if (p->type == T8_DPYRAMID_FIRST_TYPE && face == 1) {
        first_desc->x |= off_set;
      }
      else if (p->type == T8_DPYRAMID_FIRST_TYPE && face == 3) {
        first_desc->y |= off_set;
      }
      else if (p->type == T8_DPYRAMID_SECOND_TYPE && face == 4) {
        first_desc->z |= off_set;
      }
    }
  }
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc,
                             const int level)
{
  t8_linearidx_t      id = 0, t_id;
  int                 exponent;
  T8_ASSERT (level >= p->level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_copy (p, desc);
    desc->level = level;
    /* Shift the coords to the eights cube. The type of the last descendant is
     * is the type of the input pyramid*/
    t8_dpyramid_coord_t coord_offset =
      T8_DPYRAMID_LEN (p->level) - T8_DPYRAMID_LEN (level);
    desc->x |= coord_offset;
    desc->y |= coord_offset;
    desc->z |= coord_offset;
  }
  else {
    /* Compute current id, shift it to the id of the last desendant and compute it
     * via its id*/
    t_id = t8_dpyramid_linear_id (p, level);
    exponent = level - p->level;
    id = (((t8_linearidx_t) 1) << 3 * exponent) - 1;
    id += t_id;
    t8_dpyramid_init_linear_id (desc, level, id);
  }
}

void
t8_dpyramid_last_descendant_face (const t8_dpyramid_t *p,
                                  const int face, t8_dpyramid_t *last_desc,
                                  const int level)
{
  /*Computation is similar to first-descendant-face */
  int                 corner;
  t8_dpyramid_coord_t off_set =
    T8_DPYRAMID_LEN (p->level) - T8_DPYRAMID_LEN (level);

  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->level <= level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    T8_ASSERT (0 <= face && face < T8_DTET_FACES);
    corner =
      SC_MAX (t8_dtet_face_corner[face][1], t8_dtet_face_corner[face][2]);
    t8_dpyramid_corner_descendant (p, last_desc, corner, level);
  }
  else {
    t8_dpyramid_copy (p, last_desc);
    last_desc->level = level;
    if ((p->type == T8_DPYRAMID_FIRST_TYPE && face != 4) ||
        (p->type == T8_DPYRAMID_SECOND_TYPE
         && (face == 0 || face == 2 || face == 4))) {
      /*No shifting needed */
      t8_dpyramid_last_descendant (p, last_desc, level);
    }
    else if (p->type == T8_DPYRAMID_SECOND_TYPE && face == 1) {
      last_desc->x |= off_set;
      last_desc->z |= off_set;
    }
    else if (p->type == T8_DPYRAMID_SECOND_TYPE && face == 3) {
      last_desc->y |= off_set;
      last_desc->z |= off_set;
    }
    else if (p->type == T8_DPYRAMID_FIRST_TYPE && face == 4) {
      last_desc->x |= off_set;
      last_desc->y |= off_set;
    }
  }
}

int
t8_dpyramid_num_corners (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_CORNERS;
  }
  else {
    return T8_DPYRAMID_CORNERS;
  }
}

int
t8_dpyramid_num_children (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_CHILDREN;
  }
  else {
    return T8_DPYRAMID_CHILDREN;
  }
}

int
t8_dpyramid_num_siblings (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_t       parent;
  if (p->level == 0) {
    /* A level zero pyramid has only itself as sibling. */
    return 1;
  }
  t8_dpyramid_parent (p, &parent);
  return t8_dpyramid_num_children (&parent);
}

int
t8_dpyramid_num_faces (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_FACES;
  }
  else {
    return T8_DPYRAMID_FACES;
  }
}

int
t8_dpyramid_max_num_faces (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_FACES;
  }
  else {
    return T8_DPYRAMID_FACES;
  }
}

void
t8_dpyramid_boundary_face (const t8_dpyramid_t *p, const int face,
                           t8_element_t *boundary)
{
  /* face is face of of p */
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  if (face == 4) {
    /*On the bottom every face is a quad */
    /*Coordinates are scaled, because quad and pyra might have different root-len */
    p4est_quadrant_t   *q = (p4est_quadrant_t *) boundary;
    q->x = ((int64_t) p->x * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->y = ((int64_t) p->y * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->level = p->level;
  }
  else {
    /* Boundary-face is a triangle. */
    /* p-x give t->x for root-face 2,3 and p->y gives t->x for 0,1.
     * t->y is determined by p->z*/
    t8_dtri_t          *t = (t8_dtri_t *) boundary;
    t->level = p->level;
    t->y = p->z * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
    if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
      t->type = 0;
      switch (face) {
      case 0:
        t->x = p->y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 1:
        t->x = p->y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 2:
        t->x = p->x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 3:
        t->x = p->x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }
    else {
      /*Boundary is given by a tet-surface. The cases are ordered by root-face-
       * enumeration*/
      if ((face == 1 && p->type == 0) || (face == 2 && p->type == 2)) {
        t->x = p->y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->type == 0 ? 1 : 0;
      }
      else if (face == 0 && (p->type == 0 || p->type == 1)) {
        t->x = p->y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->type == 0 ? 1 : 0;
      }
      else if ((face == 1 && p->type == 3) || (face == 2 && p->type == 1)) {
        t->x = p->x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->type == 3 ? 1 : 0;
      }
      else {
        t->x = p->x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->type == 3 ? 1 : 0;
      }
    }
    T8_ASSERT (t->type == 0 || t->type == 1);
  }
}

int
t8_dpyramid_extrude_face (const t8_element_t *face, t8_dpyramid_t *p,
                          const int root_face)
{
  T8_ASSERT (0 <= root_face && root_face < T8_DPYRAMID_FACES);

  int                 extruded_face;
  if (root_face == 4) {
    /* Pyramids on the bottom are always type 6 pyramids at the bottom. We need to
     * scale the coordinates, since a quad and a pyra can have different root-len,
     * depending on their maxlvl.*/
    p4est_quadrant_t   *q = (p4est_quadrant_t *) face;
    /*Typecast to int64, we multiply two (possible at max) int32 */
    p->x = ((int64_t) q->x * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->y = ((int64_t) q->y * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->z = 0;
    p->type = T8_DPYRAMID_ROOT_TPYE;
    p->level = q->level;
    return root_face;
  }
  else {
    /*t->y gives the height of the pyramid, t->x gives the p->x or p->y, depending on
     * the root_face. The other coordinate is determined by the root_face.*/
    t8_dtri_t          *t = (t8_dtri_t *) face;
    p->z = ((int64_t) t->y * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    /* level is the same */
    p->level = t->level;
    switch (root_face) {
    case 0:
      p->x = p->z;
      /*Typecast to int64, we multiply two (possible at max) int32 */
      p->y = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 1:
      p->x = T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->level);
      p->y = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 2:
      p->x = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->y = p->z;
      break;
    case 3:
      p->x = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->y = T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->level);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    /*Description of triangles extruding to a pyramid */
    if ((t->y == (t->x & t->y)) && t->type == 0) {
      /*type zero in a pyramid extend to a pyramid */
      p->type = T8_DPYRAMID_FIRST_TYPE;
      extruded_face = root_face;
    }
    else {
      /*type 0 not in a pyramid extend to a tetrahedron */
      p->type = t8_dpyramid_tritype_rootface_to_tettype[t->type][root_face];
      extruded_face =
        t8_dpyramid_tritype_rootface_to_face[t->type][root_face];
    }
  }
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /*return the face-number of the extruded face */
  return extruded_face;
}

int
t8_dpyramid_child_id_unknown_parent (const t8_dpyramid_t *p,
                                     t8_dpyramid_t *parent)
{
  T8_ASSERT (p->level >= 0);
  if (p->level == 0) {
    /* P has no parent. We deliberately set the type of the parent to -1,
     * because a) parent should not be used and setting the type to an invalid
     *            value improves chances of catching a non-allowed use of it in an assertion later.
     *         b) t8_dpyramid_successor caused a compiler warning of parent.type may be used uninitialize
     *            in a subsequent call to t8_dpyramid_shape. This case is actually never executed, but the
     *            compiler doesn't know about it. To prevent this warning, we set the type here.
     */
    parent->type = -1;
    parent->level = -1;
    return -1;
  }
  t8_dpyramid_parent (p, parent);
  return t8_dpyramid_child_id_known_parent (p, parent);

}

int
t8_dpyramid_child_id_known_parent (const t8_dpyramid_t *p,
                                   t8_dpyramid_t *parent)
{
  t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->level);
  if (t8_dpyramid_shape (parent) == T8_ECLASS_PYRAMID) {
    T8_ASSERT (t8_dpyramid_type_cid_to_Iloc[p->type][cube_id] >= 0);
    return t8_dpyramid_type_cid_to_Iloc[p->type][cube_id];
  }
  else {
    return t8_dtet_child_id (p);
  }
}

int
t8_dpyramid_child_id (const t8_dpyramid_t *p)
{
  T8_ASSERT (p->level >= 0);
  t8_dpyramid_t       parent;
  if (p->level == 0)
    return -1;
  return t8_dpyramid_child_id_unknown_parent (p, &parent);
}

void
t8_dpyramid_child (const t8_dpyramid_t *elem, const int child_id,
                   t8_dpyramid_t *child)
{

  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (elem->level + 1);
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
  T8_ASSERT (0 <= elem->level && elem->level <= T8_DPYRAMID_MAXLEVEL);

  if (t8_dpyramid_shape (elem) == T8_ECLASS_TET) {
    t8_dtet_child (elem, child_id, child);
  }
  else {
    /* Compute the cube id and shift the coordinates accordingly */
    cube_id = t8_dpyramid_parenttype_Iloc_to_cid[elem->type][child_id];
    T8_ASSERT (cube_id >= 0);
    child->level = elem->level + 1;
    child->x = elem->x + ((cube_id & 0x01) ? length : 0);
    child->y = elem->y + ((cube_id & 0x02) ? length : 0);
    child->z = elem->z + ((cube_id & 0x04) ? length : 0);
    child->type = t8_dpyramid_parenttype_Iloc_to_type[elem->type][child_id];
  }
  T8_ASSERT (child->type >= 0);
}

void
t8_dpyramid_children (const t8_dpyramid_t *p, t8_dpyramid_t **c)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  int                 i, num_children;
  num_children = t8_dpyramid_num_children (p);
  for (i = num_children - 1; i >= 0; i--) {
    t8_dpyramid_child (p, i, c[i]);
  }
}

void
t8_dpyramid_children_at_face (const t8_dpyramid_t *p, const int face,
                              t8_dpyramid_t *children[],
                              const int num_children, int *child_indices)
{
  T8_ASSERT (num_children == T8_DPYRAMID_FACE_CHILDREN);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    /*Use tet-algo */
    t8_dtet_children_at_face (p, face, children, num_children, child_indices);
  }
  else {
    int                *children_at_face_id,
      children_at_face_id_local[T8_DPYRAMID_FACE_CHILDREN], i;
    if (child_indices != NULL) {
      children_at_face_id = child_indices;
    }
    else {
      children_at_face_id = children_at_face_id_local;
    }
    /*Fill the child ids with the child-ids at the face */
    for (i = 0; i < T8_DPYRAMID_FACE_CHILDREN; i++) {
      children_at_face_id[i] =
        t8_dpyramid_type_face_to_children_at_face[p->type -
                                                  T8_DPYRAMID_FIRST_TYPE]
        [face][i];
    }
    /*Compute the children */
    for (i = T8_DPYRAMID_FACE_CHILDREN - 1; i >= 0; i--) {
      t8_dpyramid_child (p, children_at_face_id[i], children[i]);
    }

  }
}

int
t8_dpyramid_face_child_face (const t8_dpyramid_t *p, const int face,
                             const int face_child)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= face_child && face_child < T8_DPYRAMID_CHILDREN);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return t8_dtet_face_child_face (p, face, face_child);
  }
  else {
    int                 child_face =
      t8_dpyramid_type_face_to_child_face[p->type -
                                          T8_DPYRAMID_FIRST_TYPE][face]
      [face_child];
    T8_ASSERT (child_face >= 0 && child_face <= T8_DPYRAMID_FACES);
    return child_face;
  }
}

t8_element_shape_t
t8_dpyramid_face_shape (const t8_dpyramid_t *pyra, int face)
{
  T8_ASSERT (0 <= face && face <= T8_DPYRAMID_FACES);
  if (t8_dpyramid_shape (pyra) == T8_ECLASS_TET) {
    return T8_ECLASS_TRIANGLE;
  }
  else if (face != 4) {
    return T8_ECLASS_TRIANGLE;
  }
  else {
    return T8_ECLASS_QUAD;
  }
}

int
t8_dpyramid_get_face_corner (const t8_dpyramid_t *pyra, int face, int corner)
{
  T8_ASSERT (0 <= face && face <= T8_DPYRAMID_FACES);
  if (t8_dpyramid_shape (pyra) == T8_ECLASS_TET) {
    return t8_dtet_face_corner[face][corner];
  }
  else {
    int                 corner_number = t8_dpyramid_face_corner[face][corner];
    T8_ASSERT (0 <= corner_number && corner_number < T8_DPYRAMID_FACES);
    return corner_number;
  }
}

int
t8_dpyramid_is_inside_pyra (const t8_dpyramid_t *p,
                            const t8_dpyramid_t *check)
{
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (check->level);
  t8_dpyramid_coord_t diff = p->z - check->z;

  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);

  /*test if p is inside check of type 6 */
  if (((check->x + diff) <= p->x && p->x < (check->x + length)) &&
      ((check->y + diff) <= p->y && p->y < (check->y + length)) &&
      (check->z <= p->z && p->z < (check->z + length))) {
    if ((check->x + diff == p->x && (p->type == 3 || p->type == 1)) ||
        (check->y + diff == p->y && (p->type == 0 || p->type == 2))) {
      /*tet touches face of pyra but is outside of pyra */
      return 0;
    }
    else {
      /*tet is inside pyra of type 6 */
      return T8_DPYRAMID_FIRST_TYPE;
    }
  }
  else if ((check->x <= p->x && p->x <= (check->x + diff)) &&
           (check->y <= p->y && p->y <= (check->y + diff)) &&
           (check->z <= p->z && p->z < (check->z + length))) {
    if ((check->x + diff == p->x && (p->type == 0 || p->type == 2)) ||
        (check->y + diff == p->y && (p->type == 3 || p->type == 1))) {
      /*tet touches face of pyra, but is outside of pyra */
      return 0;
    }
    else {
      /*tet is inside pyra of type 7 */
      return T8_DPYRAMID_SECOND_TYPE;
    }
  }
  else {
    /*tet is inside tet */
    return 0;
  }
}

/* Check, if a pyramid in the shape of a tet lies inside a tetrahedron
 * The i first bits give the anchorcoordinate for a possible ancestor of level i
 * for p.
 * We can store the last tetrahedra ancestor in anc.*/
int
t8_dpyramid_is_inside_tet (const t8_dpyramid_t *p, const int level,
                           t8_dpyramid_t *anc)
{
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (p->type == 0 || p->type == 3);
  int                 i;
  t8_dpyramid_coord_t coord_at_level;
  /*the tet is initialized, the ancestor will be computed */
  t8_dtet_t           tet;
  tet.x = 0;
  tet.y = 0;
  tet.z = 0;
  for (i = 1; i < level; i++) {
    /*Update the coordinate of tet to i first bits */
    coord_at_level = (1 << (T8_DPYRAMID_MAXLEVEL - i));
    tet.x = tet.x | (p->x & coord_at_level);
    tet.y = tet.y | (p->y & coord_at_level);
    tet.z = tet.z | (p->z & coord_at_level);
    tet.level = i;
    if (t8_dpyramid_is_inside_pyra (p, &tet) == 0) {
      /*p is inside a tet */
      if (anc != NULL) {
        t8_dtet_ancestor (p, i, anc);
      }
      return i;
    }
  }
  /*No matching tet-ancestor was found, the parent is a pyramid */
  return 0;
}

void
t8_dpyramid_tetparent_type (const t8_dpyramid_t *p, t8_dpyramid_t *parent)
{
  if ((p->z >> (T8_DPYRAMID_MAXLEVEL - p->level)) % 2 == 0) {
    parent->type = T8_DPYRAMID_FIRST_TYPE;
  }
  else {
    parent->type = T8_DPYRAMID_SECOND_TYPE;
  }
}

void
t8_dpyramid_parent (const t8_dpyramid_t *p, t8_dpyramid_t *parent)
{
  /*t8_debugf("parent: %i %i %i %i %i\n", p->x, p->y, p->z, p->type, p->level); */
  T8_ASSERT (p->level > 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The parent of a pyramid is a pyramid, maybe of different type */

    t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->level);

    parent->type =
      t8_dpyramid_type_cid_to_parenttype[p->type -
                                         T8_DPYRAMID_FIRST_TYPE][cube_id];
    parent->x = p->x & ~length;
    parent->y = p->y & ~length;
    parent->z = p->z & ~length;
    T8_ASSERT (parent->type >= 0);
    parent->level = p->level - 1;
  }
  else if (p->type != 0 && p->type != 3) {
    /* The direct tet-child of a pyra has type 0 or type 3, therefore
     * in this case the parent is a tetrahedron*/
    t8_dtet_parent (p, parent);

  }
  else if (t8_dpyramid_is_inside_tet (p, p->level, NULL) != 0) {
    /*Pyramid- / tetparent detection */
    /*If a tetrahedron does not reach a "significant point" its parent is a tet */
    /*Tetcase */ ;
    t8_dtet_parent (p, parent);
  }
  else {
    /*p does not lie in larger tet => parent is pyra */
    t8_dpyramid_tetparent_type (p, parent);
    parent->x = p->x & ~length;
    parent->y = p->y & ~length;
    parent->z = p->z & ~length;
    parent->level = p->level - 1;
  }
  T8_ASSERT (parent->level >= 0);
}

t8_element_shape_t
t8_dpyramid_shape (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL);
  /*The pyramid has the shape of a tetrahedron */
  if (p->type < T8_DPYRAMID_FIRST_TYPE) {
    return T8_ECLASS_TET;
  }
  else {
    return T8_ECLASS_PYRAMID;
  }

}

static void
t8_dpyramid_successor_recursion (const t8_dpyramid_t *elem,
                                 t8_dpyramid_t *succ, t8_dpyramid_t *parent,
                                 const int level)
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
    int                 shift = T8_DPYRAMID_MAXLEVEL - level + 1;
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion (succ, succ, parent, level - 1);
    succ->level = level;
    /* bits auf level auf child 0 setzen */
    pyramid_cut_coords (succ, shift);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_dpyramid_child (parent, child_id + 1, succ);
    succ->level = level;
  }
}

void
t8_dpyramid_successor (const t8_dpyramid_t *elem, t8_dpyramid_t *succ,
                       const int level)
{
  t8_dpyramid_t       parent;
  t8_dpyramid_successor_recursion (elem, succ, &parent, level);
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t *p, const int vertex,
                            int coords[])
{
  t8_dpyramid_coord_t length;
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    length = T8_DPYRAMID_LEN (p->level);
    coords[0] = p->x;
    coords[1] = p->y;
    coords[2] = p->z;
    switch (vertex) {
    case 0:
      coords[2] += (p->type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 1:
      coords[0] += length;
      coords[2] += (p->type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 2:
      coords[1] += length;
      coords[2] += (p->type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 3:
      coords[0] += length;
      coords[1] += length;
      coords[2] += (p->type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 4:
      coords[0] += (p->type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      coords[1] += (p->type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      coords[2] += (p->type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      break;
    }
  }
  else {
    T8_ASSERT (0 <= vertex && vertex < T8_DTET_CORNERS);
    t8_dtet_compute_coords (p, vertex, coords);
  }
}

void
t8_dpyramid_vertex_reference_coords (const t8_dpyramid_t *elem,
                                     int vertex, double coords[])
{
  int                 coords_int[3];
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);
  t8_dpyramid_compute_coords (elem, vertex, coords_int);
  /*scale the coordinates onto the reference cube */
  coords[0] = coords_int[0] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[1] = coords_int[1] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[2] = coords_int[2] / (double) T8_DPYRAMID_ROOT_LEN;
}

int
t8_dpyramid_is_valid (const t8_dpyramid_t *p)
{
  int                 is_valid;
  const t8_dpyramid_coord_t max_coord =
    ((int64_t) 2 * T8_DPYRAMID_ROOT_LEN) - 1;
  t8_element_shape_t  shape = t8_dpyramid_shape (p);
  /*Check the level */
  is_valid = 0 <= p->level && p->level <= T8_DPYRAMID_MAXLEVEL;
  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->x && p->x <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->y && p->y <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->z && p->z <= max_coord;

  /*The shape can be a pyramid or a tet */
  is_valid = is_valid && (shape == T8_ECLASS_PYRAMID
                          || shape == T8_ECLASS_TET);
  /*Check the type */
  is_valid = is_valid && 0 <= p->type && p->type < T8_DPYRAMID_NUM_TYPES;

  return is_valid;
}

void
t8_dpyramid_debug_print (const t8_dpyramid_t *p)
{
  t8_debugf ("x: %i, y: %i, z: %i, type %i, level: %i\n", p->x, p->y, p->z,
             p->type, p->level);
}
