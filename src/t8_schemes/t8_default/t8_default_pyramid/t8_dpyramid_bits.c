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

  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }
  cube_id |= ((p->pyramid.x & h) ? 0x01 : 0);
  cube_id |= ((p->pyramid.y & h) ? 0x02 : 0);
  cube_id |= ((p->pyramid.z & h) ? 0x04 : 0);

  return cube_id;
}

/**
 * 
 *
 * Starting from a level where the type of \a p is known compute the type
 * of \a p at level \a level
 * 
 * \param [in]  p           Input pyramid
 * \param [in]  level       The level at which we want to know the type of \a p. Must be lower than the level of \a p
 * \param [in]  known_type  The type of \a p at \a known_level
 * \param [in]  known_level The level where we already know the type of \a p
 * \return      t8_dpyramid_type_t The type of \a p at level \a level.
 * 
 * CARFULL: This computation assumes that the shape of the element does not switch between \a known_level
 *          and \a level. 
 */
static t8_dpyramid_type_t
compute_type_same_shape_ext (const t8_dpyramid_t *p, const int level,
                             const t8_dpyramid_type_t known_type,
                             const int known_level)
{
  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_type_t  type = known_type;
  int                 i;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (known_level <= p->pyramid.level);
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

/**
 * Compute the type of a pyramid at \a level.
 * Pay attention, this function assumes that the shape of the element does not switch. 
 * 
 * 
 * \param         p 
 * \param         level 
 * \return        The type of \a p at \a level
 * 
 * CAREFUL: This computation assumes that the shape of the element does not switch between \a known_level
 *          and \a level.
 */
t8_dpyramid_type_t
compute_type_same_shape (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  return compute_type_same_shape_ext (p, level, p->pyramid.type,
                                      p->pyramid.level);
}

/**
 * Set the \a shift last bits of every coordinate to zero. 
 * 
 * \param[in, out]  p     Input pyramid
 * \param[in]       shift Number of bits to set to zero
 */
static void
t8_dpyramid_cut_coordinates (t8_dpyramid_t *p, const int shift)
{
  T8_ASSERT (0 <= shift && shift <= T8_DPYRAMID_MAXLEVEL);
  p->pyramid.x = (p->pyramid.x >> shift) << shift;
  p->pyramid.y = (p->pyramid.y >> shift) << shift;
  p->pyramid.z = (p->pyramid.z >> shift) << shift;
}

int
t8_dpyramid_ancestor_id (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_t       helper;
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_ancestor (p, level, &helper);
  return t8_dpyramid_child_id (&helper);
}

int
t8_dpyramid_is_family (const t8_dpyramid_t **fam)
{

  const int           level = fam[0]->pyramid.level;
  int                 i, type_of_first;
  t8_dpyramid_coord_t inc = T8_DPYRAMID_LEN (level), x_inc, y_inc;
  if (t8_dpyramid_shape (fam[0]) == T8_ECLASS_TET) {
    int                 is_family;
    const t8_dtet_t   **tet_fam =
      T8_ALLOC (const t8_dtet_t *, T8_DTET_CHILDREN);
    for (i = 0; i < T8_DTET_CHILDREN; i++) {
      tet_fam[i] = &fam[i]->pyramid;
    }

    is_family = t8_dtet_is_familypv ((const t8_dtet_t **) tet_fam);
    T8_FREE (tet_fam);
    return is_family;
  }
  else {
    if (level == 0) {
      return 0;
    }
    /*The type of parent is the type of the first child in z-curve-order */
    type_of_first = fam[0]->pyramid.type;
    T8_ASSERT (type_of_first == T8_DPYRAMID_FIRST_TYPE
               || type_of_first == T8_DPYRAMID_SECOND_TYPE);
    for (i = 1; i < T8_DPYRAMID_CHILDREN; i++) {
      /*All elements must have the same level to be a family */
      if (fam[i]->pyramid.level != level) {
        return 0;
      }
      /*Check if every family-member has the correct type */
      if (t8_dpyramid_parenttype_Iloc_to_type[type_of_first][i] !=
          fam[i]->pyramid.type) {
        return 0;
      }
    }

    x_inc = fam[0]->pyramid.x + inc;
    y_inc = fam[0]->pyramid.y + inc;
    /*Check the coordinates of the anchor-coordinate */
    if (type_of_first == T8_DPYRAMID_FIRST_TYPE) {
      return fam[0]->pyramid.z == fam[1]->pyramid.z
        && fam[0]->pyramid.z == fam[2]->pyramid.z
        && fam[0]->pyramid.z == fam[3]->pyramid.z
        && fam[0]->pyramid.z == fam[4]->pyramid.z
        && fam[0]->pyramid.z == fam[5]->pyramid.z
        && fam[0]->pyramid.z == fam[6]->pyramid.z
        && fam[0]->pyramid.z == fam[7]->pyramid.z
        && fam[0]->pyramid.z == fam[8]->pyramid.z
        && fam[0]->pyramid.z == (fam[9]->pyramid.z - inc)
        && fam[0]->pyramid.x == fam[3]->pyramid.x
        && fam[0]->pyramid.x == fam[4]->pyramid.x
        && x_inc == fam[1]->pyramid.x && x_inc == fam[2]->pyramid.x
        && x_inc == fam[5]->pyramid.x && x_inc == fam[6]->pyramid.x
        && x_inc == fam[7]->pyramid.x && x_inc == fam[8]->pyramid.x
        && x_inc == fam[9]->pyramid.x
        && fam[0]->pyramid.y == fam[1]->pyramid.y
        && fam[0]->pyramid.y == fam[2]->pyramid.y
        && y_inc == fam[3]->pyramid.y && y_inc == fam[4]->pyramid.y
        && y_inc == fam[5]->pyramid.y && y_inc == fam[6]->pyramid.y
        && y_inc == fam[7]->pyramid.y && y_inc == fam[8]->pyramid.y
        && y_inc == fam[9]->pyramid.y;
    }
    else {
      return fam[1]->pyramid.z == fam[0]->pyramid.z + inc
        && fam[1]->pyramid.z == fam[2]->pyramid.z
        && fam[1]->pyramid.z == fam[3]->pyramid.z
        && fam[1]->pyramid.z == fam[4]->pyramid.z
        && fam[1]->pyramid.z == fam[5]->pyramid.z
        && fam[1]->pyramid.z == fam[6]->pyramid.z
        && fam[1]->pyramid.z == fam[7]->pyramid.z
        && fam[1]->pyramid.z == fam[8]->pyramid.z
        && fam[1]->pyramid.z == fam[9]->pyramid.z
        && fam[0]->pyramid.x == fam[1]->pyramid.x
        && fam[0]->pyramid.x == fam[2]->pyramid.x
        && fam[0]->pyramid.x == fam[3]->pyramid.x
        && fam[0]->pyramid.x == fam[4]->pyramid.x
        && fam[0]->pyramid.x == fam[7]->pyramid.x
        && fam[0]->pyramid.x == fam[8]->pyramid.x
        && x_inc == fam[5]->pyramid.x && x_inc == fam[6]->pyramid.x
        && x_inc == fam[9]->pyramid.x
        && fam[0]->pyramid.y == fam[1]->pyramid.y
        && fam[0]->pyramid.y == fam[2]->pyramid.y
        && fam[0]->pyramid.y == fam[3]->pyramid.y
        && fam[0]->pyramid.y == fam[4]->pyramid.y
        && fam[0]->pyramid.y == fam[5]->pyramid.y
        && fam[0]->pyramid.y == fam[6]->pyramid.y
        && y_inc == fam[7]->pyramid.y && y_inc == fam[8]->pyramid.y
        && y_inc == fam[9]->pyramid.y;
    }
  }
}

int
t8_dpyramid_is_root_boundary (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= face && face <= T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_coord_t coord_touching_root =
    T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
  if (!t8_dpyramid_is_inside_root (p)) {
    return 0;
  }
  switch (p->pyramid.type) {
  case 0:
    return (face == 1 && p->pyramid.x == p->pyramid.z) ||
      (face == 0 && p->pyramid.x == coord_touching_root);
  case 1:
    return (face == 2 && p->pyramid.y == p->pyramid.z) ||
      (face == 0 && p->pyramid.x == coord_touching_root);
  case 2:
    return (face == 2 && p->pyramid.x == p->pyramid.z) ||
      (face == 0 && p->pyramid.y == coord_touching_root);
  case 3:
    return (face == 1 && p->pyramid.y == p->pyramid.z) ||
      (face == 0 && p->pyramid.y == coord_touching_root);
  case 4:
    return 0;                   /*type 4 tets never touch a root boundary */
  case 5:
    return 0;                   /*type 5 tets never touch a root boundary */
  case T8_DPYRAMID_FIRST_TYPE:
    switch (face) {
    case 0:
      return p->pyramid.x == p->pyramid.z;
    case 1:
      return p->pyramid.x == coord_touching_root;
    case 2:
      return p->pyramid.y == p->pyramid.z;
    case 3:
      return p->pyramid.y == coord_touching_root;
    case 4:
      return p->pyramid.z == 0;
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
  T8_ASSERT (p1->pyramid.x >= 0 && p1->pyramid.y >= 0 && p1->pyramid.z >= 0 &&
             p1->pyramid.level >= 0
             && p1->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p2->pyramid.x >= 0 && p2->pyramid.y >= 0 && p2->pyramid.z >= 0
             && p2->pyramid.level >= 0
             && p2->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  maxlvl = SC_MAX (p1->pyramid.level, p2->pyramid.level);

  id1 = t8_dpyramid_linear_id (p1, maxlvl);
  id2 = t8_dpyramid_linear_id (p2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the pyramid with the smaller level
     * is considered smaller */
    if (p1->pyramid.level == p2->pyramid.level) {
      T8_ASSERT (p1->pyramid.type == p2->pyramid.type);
      return 0;
    }
    else {
      return p1->pyramid.level - p2->pyramid.level;
    }
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
     id2 */
  return id1 < id2 ? -1 : 1;
}

int
t8_dpyramid_get_level (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  return p->pyramid.level;
}

/**
 * Computes the local index of a pyramid and updates its current global index. 
 * for further computation of init_linear_id. For each sibling that is a predecessor
 * the number of pyramids or tets on the level used in init_linear_id is substracted from
 * the id. 
 * 
 * \\param[in, out] id       The current id that will be updated
 * \\param[in] type          The type of the current pyramid
 * \\param[in] pyra          Number of pyramids to shift 
 * \\param[in] tet           Number of tets to shift
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

  p->pyramid.level = level;
  p->pyramid.x = 0;
  p->pyramid.y = 0;
  p->pyramid.z = 0;
  type = T8_DPYRAMID_ROOT_TPYE; /*This is the type of the root pyramid */
  for (i = 1; i <= level; i++) {
    offset_expo = T8_DPYRAMID_MAXLEVEL - i;
    p_sum1 >>= 3;
    p_sum2 /= 6;
    // Thy types of the tetrahedron children of pyramid are always 0 or 3
    if (type == 0 || type == 3) {
      /* Ensure that a valid tet is used in init_linear_id_with_level */
      p->pyramid.type = type;
      p->pyramid.level = i;
#if T8_ENABLE_DEBUG
      ((t8_dtet_t *) p)->eclass_int8 = T8_ECLASS_TET;
#endif
      t8_dtet_init_linear_id_with_level (&(p->pyramid), id, i, level, type);
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
    p->pyramid.x |= (cube_id % 2 == 1) ? 1 << offset_expo : 0;
    p->pyramid.y |= (cube_id == 2 || cube_id == 3 || cube_id == 6
                     || cube_id == 7) ? 1 << offset_expo : 0;
    p->pyramid.z |= (cube_id > 3) ? 1 << offset_expo : 0;
    /*Compute the type */
    type = t8_dpyramid_parenttype_Iloc_to_type[type][local_index];
    T8_ASSERT (type >= 0);
  }
  T8_ASSERT (id == 0);
  p->pyramid.type = type;
}

int
t8_dpyramid_type_at_level (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);

  if (level >= p->pyramid.level) {
    return p->pyramid.type;
  }
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /* The shape does not switch, we can use the compute_type_same_shape function */
    return compute_type_same_shape (p, level);
  }
  else {
    /* The shape may switch */
    T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
    int                 current_level = p->pyramid.level;
    int                 type_at_level = p->pyramid.type;
    /* The shape can not switch for tets that do not have type 0/3 */
    while (type_at_level != 0 && type_at_level != 3 && current_level > level) {
      /* The shape does not switch, we can use compute_type_same_shape_ext */
      current_level--;
      type_at_level =
        compute_type_same_shape_ext (p, current_level, type_at_level,
                                     current_level + 1);
    }
    if (level == current_level) {
      /* We have already reached the desired level and can return. */
      return type_at_level;
    }
    else {
      /* The shape may switch, we can use t8_dpyramid_is_inside_tet, to find out when. */
      T8_ASSERT (type_at_level == 0 || type_at_level == 3);
      T8_ASSERT (current_level > level);
      t8_dpyramid_t       last_tet;
      t8_dpyramid_t       first_pyra;
      t8_dpyramid_t       tmp_pyra;
      t8_dpyramid_copy (p, &tmp_pyra);
      tmp_pyra.pyramid.type = type_at_level;
      tmp_pyra.pyramid.level = current_level;
      /* Compute the last level where the shape is a tetrahedron. If last_tet_level == 0
       * then the parent is already a pyramid. */
      int                 last_tet_level =
        t8_dpyramid_is_inside_tet (&tmp_pyra, current_level, &last_tet);
      if (last_tet_level != 0) {
        /* last_tet was calculated by t8_dpyramid_is_inside_tet */
        current_level = SC_MAX (last_tet_level, level);
        if (current_level > last_tet_level) {
          t8_dtet_ancestor (&(p->pyramid), current_level,
                            &(last_tet.pyramid));
        }
        type_at_level = last_tet.pyramid.type;
      }
      else {
        /* last_tet_level is the current_level */
        t8_dtet_ancestor (&(p->pyramid), current_level, &(last_tet.pyramid));
        type_at_level = last_tet.pyramid.type;
      }
      /* At this point the shape switches. */
      if (current_level > level) {
        current_level--;
        t8_dpyramid_tetparent_type (&last_tet, &first_pyra);
        type_at_level = first_pyra.pyramid.type;
      }
      /* The shape is now a pyramid, and won't switch anymore. */
      while (current_level > level) {
        T8_ASSERT (type_at_level == T8_DPYRAMID_FIRST_TYPE ||
                   type_at_level == T8_DPYRAMID_SECOND_TYPE);
        current_level--;
        type_at_level =
          compute_type_same_shape_ext (p, current_level, type_at_level,
                                       current_level + 1);
      }
      T8_ASSERT (level == current_level);
      return type_at_level;
    }
  }
}

t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_linearidx_t      id = 0, pyra_shift, sum_1 = 1, sum_2 = 1, local_id;
  t8_dpyramid_t       parent, copy;
  int                 i, num_pyra, num_tet;

  t8_dpyramid_copy (p, &copy);
  copy.pyramid.type = t8_dpyramid_type_at_level (p, level);
  copy.pyramid.level = level;
  t8_dpyramid_cut_coordinates (&copy, T8_DPYRAMID_MAXLEVEL - level);

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
        t8_dpyramid_parenttype_iloc_pyra_w_lower_id[parent.pyramid.type -
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
  T8_ASSERT (p->pyramid.level >= 0);
  return id;
}

/**
 * Compute the face-neighbor of p. This function does not allocate memory for the neighbor.
 * 
 * \param[in] p               Input element
 * \param[in] face            A face of \a p
 * \param[in, out] neigh      Allocated memory, will be filled with the data of the neighbor of \a p along the given \a face.
 * \return int                The face of \a neigh that touches \a p
 */
static int
t8_dpyramid_face_neighbour (const t8_dpyramid_t *p, const int face,
                            t8_dpyramid_t *neigh)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->pyramid.level);
  neigh->pyramid.x = p->pyramid.x;
  neigh->pyramid.y = p->pyramid.y;
  neigh->pyramid.z = p->pyramid.z;
  neigh->pyramid.level = p->pyramid.level;
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*pyramid touches tet or pyra */
    /*Compute the type of the neighbour */
    if (face == 0 || face == 1) {
      neigh->pyramid.type = 3;
    }
    else if (face == 2 || face == 3) {
      neigh->pyramid.type = 0;
    }
    else {
      /*face == 4 */
      neigh->pyramid.type =
        ((p->pyramid.type ==
          T8_DPYRAMID_FIRST_TYPE) ? T8_DPYRAMID_SECOND_TYPE :
         T8_DPYRAMID_FIRST_TYPE);
    }
    /*Compute the coords of the neighbour */
    /*Do nothing for face == 0 || face == 2 */
    if (face == 1) {
      neigh->pyramid.x +=
        ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
      neigh->pyramid.y +=
        ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
    }
    else if (face == 3) {
      neigh->pyramid.x +=
        ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
      neigh->pyramid.y +=
        ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
    }
    else if (face == 4) {
      neigh->pyramid.z +=
        ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? -length : length);
    }

    return t8_dpyramid_type_face_to_nface[p->pyramid.type -
                                          T8_DPYRAMID_FIRST_TYPE][face];
  }
  else {
    /*Check if the neighbor is a tet, or a pyra */
    if (p->pyramid.type != 0 && p->pyramid.type != 3) {
      /*tets of these types never have a pyra-neighbor */
      return t8_dtet_face_neighbour (&(p->pyramid), face, &(neigh->pyramid));
    }
    if (t8_dpyramid_tet_boundary (p, face)) {
      /*tet touches pyra, compute the pyra */
      if (p->pyramid.type == 0) {
        switch (face) {
        case 0:
          neigh->pyramid.x += length;
          neigh->pyramid.type = T8_DPYRAMID_SECOND_TYPE;
          return 3;
        case 1:
          neigh->pyramid.type = T8_DPYRAMID_SECOND_TYPE;
          return 2;
        case 2:
          neigh->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
          return 2;
        case 3:
          neigh->pyramid.y -= length;
          neigh->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
          return 3;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {
        /*p->pyramid.type == 3 */
        switch (face) {
        case 0:
          neigh->pyramid.y += length;
          neigh->pyramid.type = T8_DPYRAMID_SECOND_TYPE;
          return 1;
        case 1:
          neigh->pyramid.type = T8_DPYRAMID_SECOND_TYPE;
          return 0;
        case 2:
          neigh->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
          return 0;
        case 3:
          neigh->pyramid.x -= length;
          neigh->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
          return 1;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
    }
    else {
      /*tet touches tet */
      return t8_dtet_face_neighbour (&(p->pyramid), face, &(neigh->pyramid));
    }
  }
}

int
t8_dpyramid_face_parent_face (const t8_dpyramid_t *elem, const int face)
{
  t8_dpyramid_t       parent;
  int                 child_id;
  T8_ASSERT (0 <= elem->pyramid.level
             && elem->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  /*parent is a pyramid */
  if (t8_dpyramid_shape (elem) == T8_ECLASS_PYRAMID) {
    child_id = t8_dpyramid_child_id_unknown_parent (elem, &parent);
    int                 i;
    /*If the pyramid is one of the children in the array, its face-num and the face-num
     * of the parent are the same*/
    for (i = 0; i < 4; i++) {
      if (t8_dpyramid_type_face_to_children_at_face
          [parent.pyramid.type - T8_DPYRAMID_FIRST_TYPE][face][i]
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
      return t8_dtet_face_parent_face (&(elem->pyramid), face);
    }
    /*tet with a pyramid-parent */
    else {
      /*Only tets of type 0 or 3 have a pyra-parent.pyramid. Parent can have type 6 or 7 */
      if (elem->pyramid.type == 0
          && parent.pyramid.type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 3 && face == 1) {
          return 0;
        }
        else if (child_id == 5 && face == 0) {
          return 1;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 3
               && parent.pyramid.type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 1 && face == 1) {
          return 2;
        }
        else if (child_id == 6 && face == 0) {
          return 3;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 0
               && parent.pyramid.type == T8_DPYRAMID_SECOND_TYPE) {
        if (child_id == 1 && face == 3) {
          return 1;
        }
        else if (child_id == 7 && face == 2) {
          return 0;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 3
               && parent.pyramid.type == T8_DPYRAMID_SECOND_TYPE) {
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

/**
 *  Check if anc-face is connecting to a tet or to a pyra 
 * 
 * \param[in] p       Input pyramid
 * \param[in] face    Face of an ancestor of \a p
 * \return            return non-zero if anc-face is connecting to a tet
 */
static int
t8_dpyramid_tet_pyra_face_connection (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (p->pyramid.type == 0 || p->pyramid.type == 3);
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /* Depending on its cube-id at its level and its type,
   * 3 faces of a tet connect to a pyramid, one is connecting to a tet*/
  t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);
  if ((cube_id == 2 && face != 1) || (cube_id == 6 && face != 2)) {
    return p->pyramid.type == 0;
  }
  else if ((cube_id == 1 && face != 1) || (cube_id == 5 && face != 2)) {
    return p->pyramid.type == 3;
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_t       anc;
  const int           level =
    t8_dpyramid_is_inside_tet (p, p->pyramid.level, &anc);
  int                 bey_id, i, type_temp, valid_touch;
  t8_dpyramid_cube_id_t cube_id;
  if (level == 0) {
    /*Check if the face is connecting to a tet or to a pyra */
    return t8_dpyramid_tet_pyra_face_connection (p, face);
  }
  /*Check if anc-face is connecting to a tet or to a pyra */
  valid_touch = t8_dpyramid_tet_pyra_face_connection (&anc, face);
  if (valid_touch) {
    type_temp = p->pyramid.type;
    /* If so, check if the tet always lies in the corner of of its parent at this face.
     * Otherwise, the neighbor is a tet*/
    for (i = p->pyramid.level; i > anc.pyramid.level; i--) {
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_is_root_boundary (p, face)) {
    if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
      /*If p is a pyramid and touches the boundary, the face-number is the same */
      return face;
    }
    else {
      /*p is a tet and in some occasions p shares a face with its tree */
      if (face == 0 && (p->pyramid.type == 3 || p->pyramid.type == 2)) {
        return 3;
      }
      else if (face == 0 && (p->pyramid.type == 0 || p->pyramid.type == 1)) {
        return 1;
      }
      else if ((face == 1 && p->pyramid.type == 3)
               || (face == 2 && p->pyramid.type == 1)) {
        return 2;
      }
      else if ((face == 1 && p->pyramid.type == 0)
               || (face == 2 && p->pyramid.type == 2)) {
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*Check, if p is root pyramid */
  if (p->pyramid.level == 0) {
    return p->pyramid.type == T8_DPYRAMID_ROOT_TPYE && p->pyramid.x == 0
      && p->pyramid.y == 0 && p->pyramid.z == 0;
  }
  /*Check, if all coordinates are in the limit set by the length of root */
  if ((0 <= p->pyramid.z) && (p->pyramid.z < T8_DPYRAMID_ROOT_LEN) &&
      (p->pyramid.x >= p->pyramid.z) && (p->pyramid.x < T8_DPYRAMID_ROOT_LEN)
      && (p->pyramid.y >= p->pyramid.z)
      && (p->pyramid.y < T8_DPYRAMID_ROOT_LEN)) {
    if ((p->pyramid.x == p->pyramid.z
         && (p->pyramid.type == 3 || p->pyramid.type == 5))
        || (p->pyramid.y == p->pyramid.z
            && (p->pyramid.type == 0 || p->pyramid.type == 4))) {
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*Compute the face-neighbor, then check if it is inside root */
  *neigh_face = t8_dpyramid_face_neighbour (p, face, neigh);
  return t8_dpyramid_is_inside_root (neigh);
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc,
                              const int level)
{
  t8_linearidx_t      id;
  T8_ASSERT (level >= p->pyramid.level);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The first descendant of a pyramid has the same anchor coords, but another level */
    t8_dpyramid_copy (p, desc);
    desc->pyramid.level = level;
  }
  else {
    id = t8_dpyramid_linear_id (p, level);
    t8_dpyramid_init_linear_id (desc, level, id);
  }
  T8_ASSERT (p->pyramid.x <= desc->pyramid.x
             && p->pyramid.y <= desc->pyramid.y
             && p->pyramid.z <= desc->pyramid.z);
}

/** Compute the descendant at a corner of a tet up to a given level.
 *  You can not use the tetrahedron algorithm here, because it depends on the linear-id
 *  computation of a tet, which is different to the linear id for pyramids. 
 *  \param[in] p        Input pyramd
 *  \param[in, out] d   Allocated element to fill with the data of element laying in the corner \a corner of \a p
 *  \param[in] corner   A corner of \a p
 *  \param[in] level    The level to compute the corner-descendant at. 
 */
static void
t8_dpyramid_corner_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *d,
                               const int corner, const int level)
{
  int                 child_id, i;
  T8_ASSERT (p->pyramid.level <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (0 <= corner && corner < T8_DTET_CORNERS);
  if (corner == 0) {
    t8_dpyramid_first_descendant (p, d, level);
  }
  else if (corner == 1 || corner == 2) {
    /*The child at this corner is iterativle the child at child-id up to the
       given level */
    child_id = t8_dtet_parenttype_beyid_to_Iloc[p->pyramid.type][corner];
    t8_dpyramid_t       tmp;
    t8_dpyramid_copy (p, &tmp);
    for (i = p->pyramid.level; i < level; i++) {
      t8_dpyramid_child (&tmp, child_id, d);
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
  t8_dpyramid_coord_t off_set = T8_DPYRAMID_LEN (p->pyramid.level) -
    T8_DPYRAMID_LEN (level);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->pyramid.level <= level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    corner = t8_dtet_face_corner[face][0];
    t8_dpyramid_corner_descendant (p, first_desc, corner, level);
  }
  else if (p->pyramid.level == T8_DPYRAMID_MAXLEVEL) {
    t8_dpyramid_copy (p, first_desc);
  }
  else {
    /*Shift the coordinate of p to compute the descendant */
    if ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE
         && (face == 0 || face == 2 || face == 4))
        || (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && face != 4)) {
      /*No shifting is needed, fd is the first child with given level */
      t8_dpyramid_child (p, 0, first_desc);
      first_desc->pyramid.level = level;
    }
    else {
      /*shifting is needed, the fd does not have the same coord as p */
      t8_dpyramid_copy (p, first_desc);
      first_desc->pyramid.level = level;
      if (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && face == 1) {
        first_desc->pyramid.x |= off_set;
      }
      else if (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && face == 3) {
        first_desc->pyramid.y |= off_set;
      }
      else if (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && face == 4) {
        first_desc->pyramid.z |= off_set;
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
  T8_ASSERT (level >= p->pyramid.level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_copy (p, desc);
    desc->pyramid.level = level;
    /* Shift the coords to the eights cube. The type of the last descendant is
     * is the type of the input pyramid*/
    t8_dpyramid_coord_t coord_offset =
      T8_DPYRAMID_LEN (p->pyramid.level) - T8_DPYRAMID_LEN (level);
    desc->pyramid.x |= coord_offset;
    desc->pyramid.y |= coord_offset;
    desc->pyramid.z |= coord_offset;
  }
  else {
    /* Compute current id, shift it to the id of the last desendant and compute it
     * via its id*/
    t_id = t8_dpyramid_linear_id (p, level);
    exponent = level - p->pyramid.level;
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
    T8_DPYRAMID_LEN (p->pyramid.level) - T8_DPYRAMID_LEN (level);

  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->pyramid.level <= level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    T8_ASSERT (0 <= face && face < T8_DTET_FACES);
    corner =
      SC_MAX (t8_dtet_face_corner[face][1], t8_dtet_face_corner[face][2]);
    t8_dpyramid_corner_descendant (p, last_desc, corner, level);
  }
  else {
    t8_dpyramid_copy (p, last_desc);
    last_desc->pyramid.level = level;
    if ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && face != 4) ||
        (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE
         && (face == 0 || face == 2 || face == 4))) {
      /*No shifting needed */
      t8_dpyramid_last_descendant (p, last_desc, level);
    }
    else if (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && face == 1) {
      last_desc->pyramid.x |= off_set;
      last_desc->pyramid.z |= off_set;
    }
    else if (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && face == 3) {
      last_desc->pyramid.y |= off_set;
      last_desc->pyramid.z |= off_set;
    }
    else if (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && face == 4) {
      last_desc->pyramid.x |= off_set;
      last_desc->pyramid.y |= off_set;
    }
  }
}

int
t8_dpyramid_num_corners (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_t       parent;
  if (p->pyramid.level == 0) {
    /* A level zero pyramid has only itself as sibling. */
    return 1;
  }
  t8_dpyramid_parent (p, &parent);
  return t8_dpyramid_num_children (&parent);
}

int
t8_dpyramid_num_faces (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (face == 4) {
    /*On the bottom every face is a quad */
    /*Coordinates are scaled, because quad and pyra might have different root-len */
    p4est_quadrant_t   *q = (p4est_quadrant_t *) boundary;
    q->x = ((int64_t) p->pyramid.x * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->y = ((int64_t) p->pyramid.y * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->level = p->pyramid.level;
  }
  else {
    /* Boundary-face is a triangle. */
    /* p-x give t->x for root-face 2,3 and p->pyramid.y gives t->x for 0,1.
     * t->y is determined by p->pyramid.z*/
    t8_dtri_t          *t = (t8_dtri_t *) boundary;
    t->level = p->pyramid.level;
    t->y = p->pyramid.z * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
    if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
      t->type = 0;
      switch (face) {
      case 0:
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 1:
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 2:
        t->x = p->pyramid.x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      case 3:
        t->x = p->pyramid.x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
    }
    else {
      /*Boundary is given by a tet-surface. The cases are ordered by root-face-
       * enumeration*/
      if ((face == 1 && p->pyramid.type == 0)
          || (face == 2 && p->pyramid.type == 2)) {
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 0 ? 1 : 0;
      }
      else if (face == 0 && (p->pyramid.type == 0 || p->pyramid.type == 1)) {
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 0 ? 1 : 0;
      }
      else if ((face == 1 && p->pyramid.type == 3)
               || (face == 2 && p->pyramid.type == 1)) {
        t->x = p->pyramid.x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 3 ? 1 : 0;
      }
      else {
        t->x = p->pyramid.x * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 3 ? 1 : 0;
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
    p->pyramid.x = ((int64_t) q->x * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->pyramid.y = ((int64_t) q->y * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->pyramid.z = 0;
    p->pyramid.type = T8_DPYRAMID_ROOT_TPYE;
    p->pyramid.level = q->level;
    return root_face;
  }
  else {
    /*t->y gives the height of the pyramid, t->x gives the p->pyramid.x or p->pyramid.y, depending on
     * the root_face. The other coordinate is determined by the root_face.*/
    t8_dtri_t          *t = (t8_dtri_t *) face;
    p->pyramid.z = ((int64_t) t->y * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    /* level is the same */
    p->pyramid.level = t->level;
    switch (root_face) {
    case 0:
      p->pyramid.x = p->pyramid.z;
      /*Typecast to int64, we multiply two (possible at max) int32 */
      p->pyramid.y =
        ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 1:
      p->pyramid.x =
        T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
      p->pyramid.y =
        ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 2:
      p->pyramid.x =
        ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->pyramid.y = p->pyramid.z;
      break;
    case 3:
      p->pyramid.x =
        ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->pyramid.y =
        T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    /*Description of triangles extruding to a pyramid */
    if ((t->y == (t->x & t->y)) && t->type == 0) {
      /*type zero in a pyramid extend to a pyramid */
      p->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
      extruded_face = root_face;
    }
    else {
      /*type 0 not in a pyramid extend to a tetrahedron */
      p->pyramid.type =
        t8_dpyramid_tritype_rootface_to_tettype[t->type][root_face];
      extruded_face =
        t8_dpyramid_tritype_rootface_to_face[t->type][root_face];
    }
  }
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*return the face-number of the extruded face */
  return extruded_face;
}

int
t8_dpyramid_child_id_unknown_parent (const t8_dpyramid_t *p,
                                     t8_dpyramid_t *parent)
{
  T8_ASSERT (p->pyramid.level >= 0);
  if (p->pyramid.level == 0) {
    /* P has no parent.pyramid. We deliberately set the type of the parent to -1,
     * because a) parent should not be used and setting the type to an invalid
     *            value improves chances of catching a non-allowed use of it in an assertion later.
     *         b) t8_dpyramid_successor caused a compiler warning of parent.pyramid.type may be used uninitialize
     *            in a subsequent call to t8_dpyramid_shape. This case is actually never executed, but the
     *            compiler doesn't know about it. To prevent this warning, we set the type here.
     */
    parent->pyramid.type = -1;
    parent->pyramid.level = -1;
    return -1;
  }
  t8_dpyramid_parent (p, parent);
  return t8_dpyramid_child_id_known_parent (p, parent);

}

int
t8_dpyramid_child_id_known_parent (const t8_dpyramid_t *p,
                                   t8_dpyramid_t *parent)
{
  t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);
  if (t8_dpyramid_shape (parent) == T8_ECLASS_PYRAMID) {
    T8_ASSERT (t8_dpyramid_type_cid_to_Iloc[p->pyramid.type][cube_id] >= 0);
    return t8_dpyramid_type_cid_to_Iloc[p->pyramid.type][cube_id];
  }
  else {
    return t8_dtet_child_id (&(p->pyramid));
  }
}

int
t8_dpyramid_child_id (const t8_dpyramid_t *p)
{
  T8_ASSERT (p->pyramid.level >= 0);
  t8_dpyramid_t       parent;
  if (p->pyramid.level == 0)
    return -1;
  return t8_dpyramid_child_id_unknown_parent (p, &parent);
}

void
t8_dpyramid_child (const t8_dpyramid_t *elem, const int child_id,
                   t8_dpyramid_t *child)
{

  t8_dpyramid_cube_id_t cube_id;
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (elem->pyramid.level + 1);
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
  T8_ASSERT (0 <= elem->pyramid.level
             && elem->pyramid.level <= T8_DPYRAMID_MAXLEVEL);

  if (t8_dpyramid_shape (elem) == T8_ECLASS_TET) {
    /*Todo: Set switch_shape_at_level */
    t8_dtet_child (&(elem->pyramid), child_id, &(child->pyramid));
  }
  else {
    /* Compute the cube id and shift the coordinates accordingly */
    cube_id =
      t8_dpyramid_parenttype_Iloc_to_cid[elem->pyramid.type][child_id];
    T8_ASSERT (cube_id >= 0);
    child->pyramid.level = elem->pyramid.level + 1;
    child->pyramid.x = elem->pyramid.x + ((cube_id & 0x01) ? length : 0);
    child->pyramid.y = elem->pyramid.y + ((cube_id & 0x02) ? length : 0);
    child->pyramid.z = elem->pyramid.z + ((cube_id & 0x04) ? length : 0);
    child->pyramid.type =
      t8_dpyramid_parenttype_Iloc_to_type[elem->pyramid.type][child_id];
  }
  T8_ASSERT (child->pyramid.type >= 0);
}

void
t8_dpyramid_children (const t8_dpyramid_t *p, t8_dpyramid_t **c)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
    t8_dtet_t         **tet_children =
      T8_ALLOC (t8_dtet_t *, T8_DTET_FACE_CHILDREN);
    int                 i;
    for (i = 0; i < T8_DTET_FACE_CHILDREN; i++) {
      tet_children[i] = T8_ALLOC (t8_dtet_t, 1);
    }
    t8_dtet_children_at_face (&(p->pyramid), face, tet_children, num_children,
                              child_indices);
    for (i = 0; i < T8_DTET_FACE_CHILDREN; i++) {
      /* TODO: set switch_shape_at_level */
      t8_dtet_copy (tet_children[i], &(children[i]->pyramid));
      T8_FREE (tet_children[i]);
    }
    T8_FREE (tet_children);
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
        t8_dpyramid_type_face_to_children_at_face[p->pyramid.type -
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
    return t8_dtet_face_child_face (&(p->pyramid), face, face_child);
  }
  else {
    int                 child_face =
      t8_dpyramid_type_face_to_child_face[p->pyramid.type -
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

/**
 * Compute if the tetrahedron \a tet lies inside a pyramid  with coordinates given by \a check.
 * Both pyramids of type 6 and 7 are tested, hence the type of \a check does not have to be set.
 * 
 * \param tet     Input pyramid in the shape of a tetrahedron 
 * \param check   Input pyramid, candidate where \a tet could lie in.
 * \return int    the type of the pyramid where tet is inside, or 0 if it does not lie in a pyramid given by the coordinates of \a check.
 */
int
t8_dpyramid_is_inside_pyra (const t8_dpyramid_t *tet,
                            const t8_dpyramid_t *check)
{
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (check->pyramid.level);
  t8_dpyramid_coord_t diff = tet->pyramid.z - check->pyramid.z;
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);

  T8_ASSERT (0 <= tet->pyramid.level
             && tet->pyramid.level <= T8_DPYRAMID_MAXLEVEL);

  /* test if tet is inside the pyramids with coordinates given by check and type 6 */
  if (((check->pyramid.x + diff) <= tet->pyramid.x
       && tet->pyramid.x < (check->pyramid.x + length))
      && ((check->pyramid.y + diff) <= tet->pyramid.y
          && tet->pyramid.y < (check->pyramid.y + length))
      && (check->pyramid.z <= tet->pyramid.z
          && tet->pyramid.z < (check->pyramid.z + length))) {
    if ((check->pyramid.x + diff == tet->pyramid.x
         && (tet->pyramid.type == 3 || tet->pyramid.type == 1))
        || (check->pyramid.y + diff == tet->pyramid.y
            && (tet->pyramid.type == 0 || tet->pyramid.type == 2))) {
      /*tet touches face of pyra but is outside of pyra */
      return 0;
    }
    else {
      /*tet is inside pyra of type 6 */
      return T8_DPYRAMID_FIRST_TYPE;
    }
  }
  /* test if tet is inside the pyramids with coordinates given by check and type 7 */
  else
    if ((check->pyramid.x <= tet->pyramid.x
         && tet->pyramid.x <= (check->pyramid.x + diff))
        && (check->pyramid.y <= tet->pyramid.y
            && tet->pyramid.y <= (check->pyramid.y + diff))
        && (check->pyramid.z <= tet->pyramid.z
            && tet->pyramid.z < (check->pyramid.z + length))) {
    if ((check->pyramid.x + diff == tet->pyramid.x
         && (tet->pyramid.type == 0 || tet->pyramid.type == 2))
        || (check->pyramid.y + diff == tet->pyramid.y
            && (tet->pyramid.type == 3 || tet->pyramid.type == 1))) {
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

/**
 * The i first bits give the anchor coordinate for a possible ancestor of level i
 * for tet.
 * We can store the last tetrahedra ancestor in anc.
 * \param[in] tet     Inpute pyramid in the shape of a tet
 * \param[in] level   the maximal level to check whether \a tet lies in a pyramid
 * \param[in] anc     Can be NULL or an allocated element. If allocated, it will be filled with the data of the last tetrahedral ancestor */
int
t8_dpyramid_is_inside_tet (const t8_dpyramid_t *tet, const int level,
                           t8_dpyramid_t *anc)
{
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);
  T8_ASSERT (tet->pyramid.type == 0 || tet->pyramid.type == 3);
  int                 i;
  t8_dpyramid_coord_t coord_at_level;
  /*the tet is initialized, the ancestor will be computed */
  t8_dpyramid_t       pyra_at_level;    /* Candidate pyramid, where the tet could lie in. */
  pyra_at_level.pyramid.x = 0;
  pyra_at_level.pyramid.y = 0;
  pyra_at_level.pyramid.z = 0;
  for (i = 1; i < level; i++) {
    /*Update the coordinate of tet to i first bits */
    coord_at_level = (1 << (T8_DPYRAMID_MAXLEVEL - i));
    pyra_at_level.pyramid.x =
      pyra_at_level.pyramid.x | (tet->pyramid.x & coord_at_level);
    pyra_at_level.pyramid.y =
      pyra_at_level.pyramid.y | (tet->pyramid.y & coord_at_level);
    pyra_at_level.pyramid.z =
      pyra_at_level.pyramid.z | (tet->pyramid.z & coord_at_level);
    pyra_at_level.pyramid.level = i;
    if (t8_dpyramid_is_inside_pyra (tet, &pyra_at_level) == 0) {
      /*tet is inside a tet */
      if (anc != NULL) {
        t8_dtet_ancestor (&(tet->pyramid), i, &(anc->pyramid));
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
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  if ((p->pyramid.z >> (T8_DPYRAMID_MAXLEVEL - p->pyramid.level)) % 2 == 0) {
    parent->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
  }
  else {
    parent->pyramid.type = T8_DPYRAMID_SECOND_TYPE;
  }
}

void
t8_dpyramid_parent (const t8_dpyramid_t *p, t8_dpyramid_t *parent)
{
  T8_ASSERT (p->pyramid.level > 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->pyramid.level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The parent of a pyramid is a pyramid, maybe of different type */

    t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);

    parent->pyramid.type =
      t8_dpyramid_type_cid_to_parenttype[p->pyramid.type -
                                         T8_DPYRAMID_FIRST_TYPE][cube_id];
    parent->pyramid.x = p->pyramid.x & ~length;
    parent->pyramid.y = p->pyramid.y & ~length;
    parent->pyramid.z = p->pyramid.z & ~length;
    T8_ASSERT (parent->pyramid.type >= 0);
    parent->pyramid.level = p->pyramid.level - 1;
  }
  else if (p->pyramid.type != 0 && p->pyramid.type != 3) {
    /* The direct tet-child of a pyra has type 0 or type 3, therefore
     * in this case the parent is a tetrahedron*/
    /*TODO: Set switch_type_at_level */
    t8_dtet_parent (&(p->pyramid), &(parent->pyramid));
  }
  else if (t8_dpyramid_is_inside_tet (p, p->pyramid.level, NULL) != 0) {
    /*Pyramid- / tetparent detection */
    /*If a tetrahedron does not reach a "significant point" its parent is a tet */
    /*Tetcase */ ;
    t8_dtet_parent (&(p->pyramid), &(parent->pyramid));
  }
  else {
    /*p does not lie in larger tet => parent is pyra */
    t8_dpyramid_tetparent_type (p, parent);
    parent->pyramid.x = p->pyramid.x & ~length;
    parent->pyramid.y = p->pyramid.y & ~length;
    parent->pyramid.z = p->pyramid.z & ~length;
    parent->pyramid.level = p->pyramid.level - 1;
  }
  T8_ASSERT (parent->pyramid.level >= 0);
}

t8_element_shape_t
t8_dpyramid_shape (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level
             && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*The pyramid has the shape of a tetrahedron */
  if (p->pyramid.type < T8_DPYRAMID_FIRST_TYPE) {
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
  succ->pyramid.level = level;
  T8_ASSERT (succ->pyramid.type >= 0);
  child_id = t8_dpyramid_child_id_unknown_parent (elem, parent);
  /*Compute number of children */
  num_children = t8_dpyramid_num_children (parent);
  T8_ASSERT (0 <= child_id && child_id < num_children);
  if (child_id == num_children - 1) {
    int                 shift = T8_DPYRAMID_MAXLEVEL - level + 1;
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion (succ, succ, parent, level - 1);
    succ->pyramid.level = level;
    /* bits auf level auf child 0 setzen */
    t8_dpyramid_cut_coordinates (succ, shift);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_dpyramid_child (parent, child_id + 1, succ);
    succ->pyramid.level = level;
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
    length = T8_DPYRAMID_LEN (p->pyramid.level);
    coords[0] = p->pyramid.x;
    coords[1] = p->pyramid.y;
    coords[2] = p->pyramid.z;
    switch (vertex) {
    case 0:
      coords[2] += (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 1:
      coords[0] += length;
      coords[2] += (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 2:
      coords[1] += length;
      coords[2] += (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 3:
      coords[0] += length;
      coords[1] += length;
      coords[2] += (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE) ? length : 0;
      break;
    case 4:
      coords[0] += (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      coords[1] += (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      coords[2] += (p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0;
      break;
    }
  }
  else {
    T8_ASSERT (0 <= vertex && vertex < T8_DTET_CORNERS);
    t8_dtet_compute_coords (&(p->pyramid), vertex, coords);
  }
}

void
t8_dpyramid_vertex_reference_coords (const t8_dpyramid_t *elem,
                                     const int vertex, double coords[])
{
  int                 coords_int[3];
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);
  t8_dpyramid_compute_coords (elem, vertex, coords_int);
  /*scale the coordinates onto the reference cube */
  coords[0] = coords_int[0] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[1] = coords_int[1] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[2] = coords_int[2] / (double) T8_DPYRAMID_ROOT_LEN;
}

/*Compute the first ancestor in the shape of a pyramid for a tet*/
static void
t8_dpyramid_first_pyra_anc (const t8_dpyramid_t *tet,
                            t8_dpyramid_t *first_pyra_anc)
{
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);
  /*There are no tets on level 0 */
  T8_ASSERT (1 <= tet->pyramid.level
             && tet->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*The ancestor has to have the shape of a pyramid */
  t8_dpyramid_t       last_tet_anc;
  /* t8_dpyramid_is_inside_tet works only for tets of type 0 or type 3 */
  if (tet->pyramid.type != 0 && tet->pyramid.type != 3) {
    /*Find the first tet-anc of type 0 or type 3 */
    t8_dpyramid_type_t  type_at_level = tet->pyramid.type;
    int                 level = tet->pyramid.level;
    while (type_at_level != 0 && type_at_level != 3) {
      level--;
      type_at_level =
        compute_type_same_shape_ext (tet, level, type_at_level, level + 1);
    }
    T8_ASSERT (level > 0);
    T8_ASSERT (type_at_level == 0 || type_at_level == 3);
    t8_dpyramid_t       tmp_tet;

    t8_dtet_ancestor (&(tet->pyramid), level, &(tmp_tet.pyramid));
    T8_ASSERT (tmp_tet.pyramid.type == 0 || tmp_tet.pyramid.type == 3);
    /* With this call tmp_tet has type 0 or type 3 and the first-pyra-anc
     * will be computed using one of the next cases. */
    t8_dpyramid_first_pyra_anc (&tmp_tet, first_pyra_anc);
  }
  else if (t8_dpyramid_is_inside_tet (tet, tet->pyramid.level, &last_tet_anc)
           != 0) {
    /*The parent of last_tet_anc is a pyramid */
    if (last_tet_anc.pyramid.level == 1) {
      first_pyra_anc->pyramid.type = 6;
    }
    else {
      t8_dpyramid_tetparent_type (&last_tet_anc, first_pyra_anc);
    }
    T8_ASSERT (last_tet_anc.pyramid.level >= 1);
    /*Update coordinates */
    t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (last_tet_anc.pyramid.level);
    first_pyra_anc->pyramid.x = last_tet_anc.pyramid.x & ~length;
    first_pyra_anc->pyramid.y = last_tet_anc.pyramid.y & ~length;
    first_pyra_anc->pyramid.z = last_tet_anc.pyramid.z & ~length;
    /*set the level */
    first_pyra_anc->pyramid.level = last_tet_anc.pyramid.level - 1;
  }
  else {
    /* The parent of the tet is already a pyramid */
    t8_dpyramid_parent (tet, first_pyra_anc);
  }
}

/**
 * Smallest level at which an anc of \a tet has the shape of a tetrahedron
 * 
 * \param[in] tet The input element
 * \return The level of the last ancestor with the shape of a tetrahedron  
 */
static int
t8_dpyramid_switches_type_at_level (const t8_dpyramid_t *tet)
{
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);
  t8_dpyramid_type_t  type_at_level = tet->pyramid.type;
  int                 level = tet->pyramid.level;
  int                 last_tet_level;
  t8_dpyramid_t       tmp_tet;
  /* A tetrahedron that has not type 0 or type 3 can not switch the shape, because
   * the tetrahedral children of a pyramid only have type 0 or type 3.*/
  while (type_at_level != 0 && type_at_level != 3) {
    level--;
    type_at_level =
      compute_type_same_shape_ext (tet, level, type_at_level, level + 1);
  }
  T8_ASSERT (type_at_level == 0 || type_at_level == 3);
  t8_dpyramid_copy (tet, &tmp_tet);
  tmp_tet.pyramid.type = type_at_level;
  tmp_tet.pyramid.level = level;
  /* t8_pyramid_is_inside computes the level where the shape switches for 
   * tets of type 0 or 3. */
  last_tet_level = t8_dpyramid_is_inside_tet (&tmp_tet, level, NULL);
  if (last_tet_level != 0) {
    return last_tet_level;
  }
  else {
    /* The parent of tmp_tet is a pyramid, hence we return the current level. */
    return level;
  }
}

/**
 * Compute the ancestor of \a pyra on level \a level.
 * 
 * \param[in]       pyra    Input pyramid
 * \param[in]       level   The level at which we want to compute \a anc
 * \param[in, out]  anc     Allocated input element which will be filled by the data of the anc of \a pyra at level \a level
 */
void
t8_dpyramid_ancestor (const t8_dpyramid_t *pyra, const int level,
                      t8_dpyramid_t *anc)
{
  T8_ASSERT (0 <= level && level <= pyra->pyramid.level);
  /*Set the coordinates and the level of the ancestor */
  t8_dpyramid_copy (pyra, anc);
  if (pyra->pyramid.level == level) {
    return;
  }
  else if (level == pyra->pyramid.level - 1) {
    /* We can reuse the parent computation if we want to go only one level up. */
    t8_dpyramid_parent (pyra, anc);
    return;
  }
  /* The coordinates of the anc are defined by the level. */
  t8_dpyramid_cut_coordinates (anc, T8_DPYRAMID_MAXLEVEL - level);
  anc->pyramid.level = level;
  anc->pyramid.type = t8_dpyramid_type_at_level (pyra, level);
}

void
t8_dpyramid_nearest_common_ancestor (const t8_dpyramid_t *pyra1,
                                     const t8_dpyramid_t *pyra2,
                                     t8_dpyramid_t *nca)
{
  /* If the input elements have different shapes, the nca has to have the
   * shape of a pyramid. The element in the shape of a tet switches the shape. */
  if (t8_dpyramid_shape (pyra1) == T8_ECLASS_PYRAMID &&
      t8_dpyramid_shape (pyra2) == T8_ECLASS_TET) {
    t8_dpyramid_t       first_pyramid_anc;

    t8_dpyramid_first_pyra_anc (pyra2, &first_pyramid_anc);

    /* pyra1 and first_pyramid_anc have the shape of a pyramid now, 
     * we can call the nca again.
     */
    t8_dpyramid_nearest_common_ancestor (pyra1, &first_pyramid_anc, nca);
    return;
  }
  else if (t8_dpyramid_shape (pyra1) == T8_ECLASS_TET &&
           t8_dpyramid_shape (pyra2) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_t       first_pyramid_anc;

    t8_dpyramid_first_pyra_anc (pyra1, &first_pyramid_anc);
    /* pyra2 and first_pyramid_anc have the shape of a pyramid now, 
     * we can call the nca again.
     */
    t8_dpyramid_nearest_common_ancestor (&first_pyramid_anc, pyra2, nca);
    return;
  }
  /* both elements have the shape of a pyramid, hence the nca */
  else if (t8_dpyramid_shape (pyra1) == T8_ECLASS_PYRAMID &&
           t8_dpyramid_shape (pyra2) == T8_ECLASS_PYRAMID) {
    /* The following computations are necessary to find the 
     * first ancestors of pyra1 and pyra2 with the same type. We 
     * have already computed the level at which they have the same
     * coordinate, but the type could be different. */
    int                 level;  /* To iterate over level */
    int                 cube_level;     /* the level of the cube where pyra1 and pyra2 have the same coords */
    int                 real_level;     /* the level of the nca */
    t8_dpyramid_coord_t maxclor;
    t8_dpyramid_type_t  p1_type_at_level;       /* type of pyra1 at level */
    t8_dpyramid_type_t  p2_type_at_level;       /* type of pyra2 at level */
    /* Compute the first level, at which the coordinates differ */
    maxclor = pyra1->pyramid.x ^ pyra2->pyramid.x;
    maxclor |= pyra1->pyramid.y ^ pyra2->pyramid.y;
    maxclor |= pyra1->pyramid.z ^ pyra2->pyramid.z;
    level = SC_LOG2_32 (maxclor) + 1;
    T8_ASSERT (level <= T8_DPYRAMID_MAXLEVEL);
    /* This is the highest possible level. The coordinates are the same,
     * but the types can be different.*/
    cube_level = SC_MIN (T8_DPYRAMID_MAXLEVEL - level,
                         (int) SC_MIN (pyra1->pyramid.level,
                                       pyra2->pyramid.level));
    real_level = cube_level;
    p1_type_at_level = compute_type_same_shape (pyra1, cube_level);
    p2_type_at_level = compute_type_same_shape (pyra2, cube_level);
    /* Iterate over the levels and compute both types at that level.
     * If they are the same, we know the level of the nearest common ancestor. */
    while (p1_type_at_level != p2_type_at_level) {
      real_level--;
      p1_type_at_level =
        compute_type_same_shape_ext (pyra1, real_level, p1_type_at_level,
                                     real_level + 1);
      p2_type_at_level =
        compute_type_same_shape_ext (pyra2, real_level, p2_type_at_level,
                                     real_level + 1);
    }
    T8_ASSERT (real_level >= 0);
    /* Fill the nca */
    t8_dpyramid_copy (pyra1, nca);
    nca->pyramid.level = real_level;
    /* Correct the coordinates of the nca */
    t8_dpyramid_cut_coordinates (nca, T8_DPYRAMID_MAXLEVEL - real_level);
    /* Set the computed type */
    nca->pyramid.type = p1_type_at_level;
    return;
  }
  else {
    /* Both elements are a tet. The ancestor can be at a level befor any of the
     * elementes switches the shape from a tet to a pyra. If one of the tets switches
     * the shape, both tets have to switch the shape. */
    T8_ASSERT (t8_dpyramid_shape (pyra1) == T8_ECLASS_TET);
    T8_ASSERT (t8_dpyramid_shape (pyra2) == T8_ECLASS_TET);
    int                 level;  /* To iterate over level */
    int                 cube_level;     /* the level of the cube where pyra1 and pyra2 have the same coords */
    int                 real_level;     /* the level of the nca */
    t8_dpyramid_coord_t maxclor;
    t8_dpyramid_type_t  p1_type_at_level;       /* type of pyra1 at level */
    t8_dpyramid_type_t  p2_type_at_level;       /* type of pyra2 at level */
    /* Compute the first level, at which the coordinates differ */
    maxclor = pyra1->pyramid.x ^ pyra2->pyramid.x;
    maxclor |= pyra1->pyramid.y ^ pyra2->pyramid.y;
    maxclor |= pyra1->pyramid.z ^ pyra2->pyramid.z;
    level = SC_LOG2_32 (maxclor) + 1;
    T8_ASSERT (level <= T8_DPYRAMID_MAXLEVEL);
    t8_dpyramid_t       pyra1_anc;
    t8_dpyramid_t       pyra2_anc;
    t8_dpyramid_t       last_tet1;
    t8_dpyramid_t       last_tet2;

    /* Cube level is the highest possible level where the nca can be. the coordinates
     * match at the level, but the type can be different.*/
    cube_level = SC_MIN (T8_DPYRAMID_MAXLEVEL - level,
                         (int) SC_MIN (pyra1->pyramid.level,
                                       pyra2->pyramid.level));
    real_level = cube_level;

    /* Get the levels where the elements switch from a tet-shape to a pyra-shape. */
    int                 level_switch_pyra1 =
      t8_dpyramid_switches_type_at_level (pyra1);
    int                 level_switch_pyra2 =
      t8_dpyramid_switches_type_at_level (pyra2);
    t8_dpyramid_ancestor (pyra1, real_level, &pyra1_anc);
    t8_dpyramid_ancestor (pyra2, real_level, &pyra2_anc);

    p1_type_at_level = pyra1_anc.pyramid.type;
    p2_type_at_level = pyra2_anc.pyramid.type;

    /* We iterate over the levels and check if the types of both tets match and 
     * stop at that level.
     * The loop is interupted, if we get to a level where one element switches 
     * the shape.*/
    while (p1_type_at_level != p2_type_at_level &&
           real_level >= level_switch_pyra1 &&
           real_level >= level_switch_pyra2) {
      real_level--;
      p1_type_at_level =
        compute_type_same_shape_ext (pyra1, real_level, p1_type_at_level,
                                     real_level + 1);
      p2_type_at_level =
        compute_type_same_shape_ext (pyra2, real_level, p2_type_at_level,
                                     real_level + 1);
    }
    if (real_level < level_switch_pyra1) {
      /* The first element switches the shape. The type is computed using 
       * assuming a pyramid-parent.pyramid.*/
      t8_dpyramid_t       first_pyra1;
      t8_dtet_ancestor (&(pyra1->pyramid), level_switch_pyra1,
                        &(last_tet1.pyramid));
      t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (level_switch_pyra1);
      first_pyra1.pyramid.x = last_tet1.pyramid.x & ~length;
      first_pyra1.pyramid.y = last_tet1.pyramid.y & ~length;
      first_pyra1.pyramid.z = last_tet1.pyramid.z & ~length;
      t8_dpyramid_tetparent_type (&last_tet1, &first_pyra1);
      first_pyra1.pyramid.level = level_switch_pyra1 - 1;
      t8_dpyramid_nearest_common_ancestor (&first_pyra1, pyra2, nca);
      return;
    }
    else if (real_level < level_switch_pyra2) {

      /* The second element switches the shape. The type is computed using 
       * assuming a pyramid-parent.pyramid.*/
      t8_dpyramid_t       first_pyra2;
      t8_dtet_ancestor (&(pyra2->pyramid), level_switch_pyra2,
                        &(last_tet2.pyramid));
      t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (level_switch_pyra2);
      first_pyra2.pyramid.x = last_tet2.pyramid.x & ~length;
      first_pyra2.pyramid.y = last_tet2.pyramid.y & ~length;
      first_pyra2.pyramid.z = last_tet2.pyramid.z & ~length;
      t8_dpyramid_tetparent_type (&last_tet2, &first_pyra2);
      first_pyra2.pyramid.level = level_switch_pyra2 - 1;
      t8_dpyramid_nearest_common_ancestor (&first_pyra2, pyra1, nca);
      return;
    }
    else {
      /* No anc switches the shape, the nca is a tet */
      T8_ASSERT (p1_type_at_level == p2_type_at_level);
      T8_ASSERT (p1_type_at_level < T8_DPYRAMID_FIRST_TYPE);
      t8_dtet_ancestor (&(pyra1->pyramid), real_level, &(nca->pyramid));
      return;
    }
  }
}

int
t8_dpyramid_is_valid (const t8_dpyramid_t *p)
{
  int                 is_valid;
  const t8_dpyramid_coord_t max_coord =
    ((int64_t) 2 * T8_DPYRAMID_ROOT_LEN) - 1;
  t8_element_shape_t  shape = t8_dpyramid_shape (p);
  /*Check the level */
  is_valid = 0 <= p->pyramid.level
    && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL;
  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.x
    && p->pyramid.x <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.y
    && p->pyramid.y <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.z
    && p->pyramid.z <= max_coord;

  /*The shape can be a pyramid or a tet */
  is_valid = is_valid && (shape == T8_ECLASS_PYRAMID
                          || shape == T8_ECLASS_TET);
  /*Check the type */
  is_valid = is_valid && 0 <= p->pyramid.type
    && p->pyramid.type < T8_DPYRAMID_NUM_TYPES;

  if (p->pyramid.level == 0) {
    is_valid = is_valid && (p->pyramid.type == T8_DPYRAMID_ROOT_TPYE);
  }

  return is_valid;
}

void
t8_dpyramid_debug_print (const t8_dpyramid_t *p)
{
  t8_debugf ("x: %i, y: %i, z: %i, type %i, level: %i\n", p->pyramid.x,
             p->pyramid.y, p->pyramid.z, p->pyramid.type, p->pyramid.level);
}
