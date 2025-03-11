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
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>

typedef int8_t t8_dpyramid_cube_id_t;

static t8_dpyramid_cube_id_t
compute_cubeid (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_cube_id_t cube_id = 0;

  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  const t8_dpyramid_coord_t h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }
  cube_id |= ((p->pyramid.x & h) ? 0x01 : 0);
  cube_id |= ((p->pyramid.y & h) ? 0x02 : 0);
  cube_id |= ((p->pyramid.z & h) ? 0x04 : 0);

  return cube_id;
}

/**
 * Starting from a level where the type of \a p is known compute the type of \a p at level \a level
 * 
 * \param [in]  p           Input pyramid
 * \param [in]  level       The level at which we want to know the type of \a p. Must be lower than the level of \a p
 * \param [in]  known_type  The type of \a p at \a known_level
 * \param [in]  known_level The level where we already know the type of \a p
 * \return      t8_dpyramid_type_t The type of \a p at level \a level.
 * 
 * WARNING: This computation assumes that the shape of the element does not switch between \a known_level
 *          and \a level. 
 */
static t8_dpyramid_type_t
compute_type_same_shape_ext (const t8_dpyramid_t *p, const int level, const t8_dpyramid_type_t known_type,
                             const int known_level)
{
  t8_dpyramid_type_t type = known_type;

  T8_ASSERT (0 <= level && level <= known_level);
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (known_level <= p->pyramid.level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /*Type of the root pyra */
    return T8_DPYRAMID_ROOT_TYPE;
  }
  for (int i = known_level; i > level; i--) {
    const t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, i);
    type = t8_dpyramid_cid_type_to_parenttype[cube_id][type];
  }
  return type;
}

/**
 * Compute the type of a pyramid at \a level.Pay attention, this function assumes that the shape of the element does 
 * not switch. 
 * \param         p 
 * \param         level 
 * \return        The type of \a p at \a level
 * WARNING: This computation assumes that the shape of the element does not switch between \a known_level
 *          and \a level.
 */
t8_dpyramid_type_t
compute_type_same_shape (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  return compute_type_same_shape_ext (p, level, p->pyramid.type, p->pyramid.level);
}

/**
 * Set the \a shift last bits of every coordinate to zero. 
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

/**
 * Compute if the tetrahedron \a tet lies inside a pyramid  with coordinates given by \a check. Both pyramids of type 6
 * and 7 are tested, hence the type of \a check does not have to be set.
 * \param tet     Input pyramid in the shape of a tetrahedron 
 * \param check   Input pyramid, candidate where \a tet could lie in.
 * \return int    the type of the pyramid where tet is inside, or 0 if it does not lie in a pyramid given by the coordinates of \a check.
 */
static int
t8_dpyramid_is_inside_pyra (const t8_dpyramid_t *tet, const t8_dpyramid_t *check)
{
  t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (check->pyramid.level);
  t8_dpyramid_coord_t diff = tet->pyramid.z - check->pyramid.z;
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);

  T8_ASSERT (0 <= tet->pyramid.level && tet->pyramid.level <= T8_DPYRAMID_MAXLEVEL);

  /* test if tet is inside the pyramids with coordinates given by check and type 6 */
  if (((check->pyramid.x + diff) <= tet->pyramid.x && tet->pyramid.x < (check->pyramid.x + length))
      && ((check->pyramid.y + diff) <= tet->pyramid.y && tet->pyramid.y < (check->pyramid.y + length))
      && (check->pyramid.z <= tet->pyramid.z && tet->pyramid.z < (check->pyramid.z + length))) {
    if ((check->pyramid.x + diff == tet->pyramid.x && (tet->pyramid.type == 3 || tet->pyramid.type == 1))
        || (check->pyramid.y + diff == tet->pyramid.y && (tet->pyramid.type == 0 || tet->pyramid.type == 2))) {
      /*tet touches face of pyra but is outside of pyra */
      return 0;
    }
    else {
      /*tet is inside pyra of type 6 */
      return T8_DPYRAMID_FIRST_TYPE;
    }
  }
  /* test if tet is inside the pyramids with coordinates given by check and type 7 */
  else if ((check->pyramid.x <= tet->pyramid.x && tet->pyramid.x <= (check->pyramid.x + diff))
           && (check->pyramid.y <= tet->pyramid.y && tet->pyramid.y <= (check->pyramid.y + diff))
           && (check->pyramid.z <= tet->pyramid.z && tet->pyramid.z < (check->pyramid.z + length))) {
    if ((check->pyramid.x + diff == tet->pyramid.x && (tet->pyramid.type == 0 || tet->pyramid.type == 2))
        || (check->pyramid.y + diff == tet->pyramid.y && (tet->pyramid.type == 3 || tet->pyramid.type == 1))) {
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
 * The i first bits give the anchor coordinate for a possible ancestor of level i for tet.
 * We can store the last tetrahedra ancestor in ancestor.
 * \param[in] tet   Inpute pyramid in the shape of a tet
 * \param[in] level the maximal level to check whether \a tet lies in a pyramid
 * \param[in] ancestor   Can be NULL or an allocated element. If allocated, it will be filled with the data of the last tetrahedral ancestor 
 * \return          0, if the pyramid is insed of a tetrahedron*/
static int
t8_dpyramid_is_inside_tet (const t8_dpyramid_t *tet, const int level, t8_dpyramid_t *ancestor)
{
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);
  T8_ASSERT (tet->pyramid.type == 0 || tet->pyramid.type == 3);
  int i;
  t8_dpyramid_coord_t coord_at_level;
  /*the tet is initialized, the ancestor will be computed */
  t8_dpyramid_t pyra_at_level; /* Candidate pyramid, where the tet could lie in. */
  pyra_at_level.pyramid.x = 0;
  pyra_at_level.pyramid.y = 0;
  pyra_at_level.pyramid.z = 0;
  for (i = 1; i < level; i++) {
    /*Update the coordinate of tet to i first bits */
    coord_at_level = (1 << (T8_DPYRAMID_MAXLEVEL - i));
    pyra_at_level.pyramid.x = pyra_at_level.pyramid.x | (tet->pyramid.x & coord_at_level);
    pyra_at_level.pyramid.y = pyra_at_level.pyramid.y | (tet->pyramid.y & coord_at_level);
    pyra_at_level.pyramid.z = pyra_at_level.pyramid.z | (tet->pyramid.z & coord_at_level);
    pyra_at_level.pyramid.level = i;
    if (t8_dpyramid_is_inside_pyra (tet, &pyra_at_level) == 0) {
      /*tet is inside a tet */
      if (ancestor != NULL) {
        t8_dtet_ancestor (&(tet->pyramid), i, &(ancestor->pyramid));
      }
      return i;
    }
  }
  /*No matching tet-ancestor was found, the parent is a pyramid */
  return 0;
}

/**
 * Smallest level at which an ancestor of \a tet has the shape of a tetrahedron
 * \param[in] tet The input element
 * \return The level of the last ancestor with the shape of a tetrahedron  
 */
static int
t8_dpyramid_compute_switch_shape_at_level (const t8_dpyramid_t *tet)
{
  T8_ASSERT (t8_dpyramid_shape (tet) == T8_ECLASS_TET);
  t8_dpyramid_type_t type_at_level = tet->pyramid.type;
  int level = tet->pyramid.level;
  t8_dpyramid_t tmp_tet;

  /* A tetrahedron that has not type 0 or type 3 can not switch the shape, because
   * the tetrahedral children of a pyramid only have type 0 or type 3.*/
  while (type_at_level != 0 && type_at_level != 3) {
    level--;
    type_at_level = compute_type_same_shape_ext (tet, level, type_at_level, level + 1);
  }
  T8_ASSERT (type_at_level == 0 || type_at_level == 3);
  t8_dpyramid_copy (tet, &tmp_tet);
  tmp_tet.pyramid.type = type_at_level;
  tmp_tet.pyramid.level = level;
  /* t8_pyramid_is_inside computes the level where the shape switches for 
   * tets of type 0 or 3. */
  int last_tet_level = t8_dpyramid_is_inside_tet (&tmp_tet, level, NULL);
  if (last_tet_level != 0) {
    return last_tet_level;
  }
  else {
    /* The parent of tmp_tet is a pyramid, hence we return the current level. */
    return level;
  }
}

/**
 * Sets the field switch_shape_at_level for \a p. \a p has to have the shape of a tetrahedron. switch_shape_at_level
 * is set to the lowest level at which the ancestor of \a p still has the shape of a tetrahedron. switch_shape_at_level
 * is set to -1 for pyramidal shaped elements.
 * \param p       Input element, whose switch_shape_at_level will be set.
 */
static void
t8_dpyramid_set_switch_shape_at_level (t8_dpyramid_t *p)
{
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    p->switch_shape_at_level = t8_dpyramid_compute_switch_shape_at_level (p);
  }
  else {
    p->switch_shape_at_level = -1;
  }
}

int
t8_dpyramid_ancestor_id (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_t helper;
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_ancestor (p, level, &helper);
  return t8_dpyramid_child_id (&helper);
}

int
t8_dpyramid_is_family (t8_dpyramid_t **fam)
{
  const int level = fam[0]->pyramid.level;
  if (t8_dpyramid_shape (fam[0]) == T8_ECLASS_TET) {
    t8_dtet_t *tet_fam[T8_DTET_CHILDREN];
    for (int i = 0; i < T8_DTET_CHILDREN; i++) {
      tet_fam[i] = &fam[i]->pyramid;
    }
    const int is_family = t8_dtet_is_familypv ((const t8_dtet_t **) tet_fam);
    return is_family;
  }
  else {
    if (level == 0) {
      return 0;
    }
    /* The type of parent is the type of the first child in z-curve-order */
    const int type_of_first = fam[0]->pyramid.type;
    T8_ASSERT (type_of_first == T8_DPYRAMID_FIRST_TYPE || type_of_first == T8_DPYRAMID_SECOND_TYPE);
    for (int i = 1; i < T8_DPYRAMID_CHILDREN; i++) {
      /* All elements must have the same level to be a family */
      if (fam[i]->pyramid.level != level) {
        return 0;
      }
      /* Check if every family-member has the correct type */
      if (t8_dpyramid_parenttype_Iloc_to_type[type_of_first][i] != fam[i]->pyramid.type) {
        return 0;
      }
    }

    t8_dpyramid_coord_t inc = T8_DPYRAMID_LEN (level), x_inc, y_inc;
    x_inc = fam[0]->pyramid.x + inc;
    y_inc = fam[0]->pyramid.y + inc;
    /* Check the coordinates of the anchor-coordinate */
    if (type_of_first == T8_DPYRAMID_FIRST_TYPE) {
      return fam[0]->pyramid.z == fam[1]->pyramid.z && fam[0]->pyramid.z == fam[2]->pyramid.z
             && fam[0]->pyramid.z == fam[3]->pyramid.z && fam[0]->pyramid.z == fam[4]->pyramid.z
             && fam[0]->pyramid.z == fam[5]->pyramid.z && fam[0]->pyramid.z == fam[6]->pyramid.z
             && fam[0]->pyramid.z == fam[7]->pyramid.z && fam[0]->pyramid.z == fam[8]->pyramid.z
             && fam[0]->pyramid.z == (fam[9]->pyramid.z - inc) && fam[0]->pyramid.x == fam[3]->pyramid.x
             && fam[0]->pyramid.x == fam[4]->pyramid.x && x_inc == fam[1]->pyramid.x && x_inc == fam[2]->pyramid.x
             && x_inc == fam[5]->pyramid.x && x_inc == fam[6]->pyramid.x && x_inc == fam[7]->pyramid.x
             && x_inc == fam[8]->pyramid.x && x_inc == fam[9]->pyramid.x && fam[0]->pyramid.y == fam[1]->pyramid.y
             && fam[0]->pyramid.y == fam[2]->pyramid.y && y_inc == fam[3]->pyramid.y && y_inc == fam[4]->pyramid.y
             && y_inc == fam[5]->pyramid.y && y_inc == fam[6]->pyramid.y && y_inc == fam[7]->pyramid.y
             && y_inc == fam[8]->pyramid.y && y_inc == fam[9]->pyramid.y;
    }
    else {
      return fam[1]->pyramid.z == fam[0]->pyramid.z + inc && fam[1]->pyramid.z == fam[2]->pyramid.z
             && fam[1]->pyramid.z == fam[3]->pyramid.z && fam[1]->pyramid.z == fam[4]->pyramid.z
             && fam[1]->pyramid.z == fam[5]->pyramid.z && fam[1]->pyramid.z == fam[6]->pyramid.z
             && fam[1]->pyramid.z == fam[7]->pyramid.z && fam[1]->pyramid.z == fam[8]->pyramid.z
             && fam[1]->pyramid.z == fam[9]->pyramid.z && fam[0]->pyramid.x == fam[1]->pyramid.x
             && fam[0]->pyramid.x == fam[2]->pyramid.x && fam[0]->pyramid.x == fam[3]->pyramid.x
             && fam[0]->pyramid.x == fam[4]->pyramid.x && fam[0]->pyramid.x == fam[7]->pyramid.x
             && fam[0]->pyramid.x == fam[8]->pyramid.x && x_inc == fam[5]->pyramid.x && x_inc == fam[6]->pyramid.x
             && x_inc == fam[9]->pyramid.x && fam[0]->pyramid.y == fam[1]->pyramid.y
             && fam[0]->pyramid.y == fam[2]->pyramid.y && fam[0]->pyramid.y == fam[3]->pyramid.y
             && fam[0]->pyramid.y == fam[4]->pyramid.y && fam[0]->pyramid.y == fam[5]->pyramid.y
             && fam[0]->pyramid.y == fam[6]->pyramid.y && y_inc == fam[7]->pyramid.y && y_inc == fam[8]->pyramid.y
             && y_inc == fam[9]->pyramid.y;
    }
  }
}

int
t8_dpyramid_is_root_boundary (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= face && face <= T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  const t8_dpyramid_coord_t coord_touching_root = T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
  if (!t8_dpyramid_is_inside_root (p)) {
    return 0;
  }
  switch (p->pyramid.type) {
  case 0:
    return (face == 1 && p->pyramid.x == p->pyramid.z) || (face == 0 && p->pyramid.x == coord_touching_root);
  case 1:
    return (face == 2 && p->pyramid.y == p->pyramid.z) || (face == 0 && p->pyramid.x == coord_touching_root);
  case 2:
    return (face == 2 && p->pyramid.x == p->pyramid.z) || (face == 0 && p->pyramid.y == coord_touching_root);
  case 3:
    return (face == 1 && p->pyramid.y == p->pyramid.z) || (face == 0 && p->pyramid.y == coord_touching_root);
  case 4:
    return 0; /*type 4 tets never touch a root boundary */
  case 5:
    return 0; /*type 5 tets never touch a root boundary */
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
    return 0; /*type 7 pyramids are never at the root boundary */
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
t8_dpyramid_equal (const t8_dpyramid_t *elem1, const t8_dpyramid_t *elem2)
{
  return t8_dtet_equal (&elem1->pyramid, &elem2->pyramid)
         && elem1->switch_shape_at_level == elem2->switch_shape_at_level;
}

int
t8_dpyramid_compare (const t8_dpyramid_t *p1, const t8_dpyramid_t *p2)
{
  T8_ASSERT (p1->pyramid.x >= 0 && p1->pyramid.y >= 0 && p1->pyramid.z >= 0 && p1->pyramid.level >= 0
             && p1->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p2->pyramid.x >= 0 && p2->pyramid.y >= 0 && p2->pyramid.z >= 0 && p2->pyramid.level >= 0
             && p2->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  const int maxlvl = SC_MAX (p1->pyramid.level, p2->pyramid.level);

  const t8_linearidx_t id1 = t8_dpyramid_linear_id (p1, maxlvl);
  const t8_linearidx_t id2 = t8_dpyramid_linear_id (p2, maxlvl);
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  return p->pyramid.level;
}

/**
 * Computes the local index of a pyramid and updates its current global index. for further computation of 
 * init_linear_id. For each sibling that is a predecessor the number of pyramids or tets on the level used in 
 * init_linear_id is subtracted from the id. 
 * \param[in, out] id       The current id that will be updated
 * \param[in] type          The type of the current pyramid
 * \param[in] pyra          Number of pyramids to shift 
 * \param[in] tet           Number of tets to shift
 * \return int              The local-id of the child
 */
static int
t8_dpyramid_update_index (t8_linearidx_t *id, const t8_dpyramid_type_t type, const t8_linearidx_t pyra,
                          const t8_linearidx_t tet)
{
  t8_linearidx_t test = 0;
  t8_linearidx_t shift;
  T8_ASSERT (id != NULL);
  T8_ASSERT (*id >= 0);
  int remain = -1;
  do {
    /* Iterate through the local-id. Get the current shift by the type of the
     * current element*/
    shift = t8_dpyramid_parenttype_Iloc_to_type[type][remain + 1] >= T8_DPYRAMID_FIRST_TYPE ? pyra : tet;
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
t8_dpyramid_init_linear_id (t8_dpyramid_t *p, const int level, t8_linearidx_t id)
{
  t8_linearidx_t p_sum1 = ((t8_linearidx_t) 1) << (3 * level);
  t8_linearidx_t p_sum2 = sc_intpow64u (6, level);

  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= id && id <= 2 * p_sum1 - p_sum2);

  p->pyramid.level = level;
  p->pyramid.x = 0;
  p->pyramid.y = 0;
  p->pyramid.z = 0;
  t8_dpyramid_type_t type = T8_DPYRAMID_ROOT_TYPE; /*This is the type of the root pyramid */
  for (int i = 1; i <= level; i++) {
    const int offset_expo = T8_DPYRAMID_MAXLEVEL - i;
    p_sum1 >>= 3;
    p_sum2 /= 6;
    // Thy types of the tetrahedron children of pyramid are always 0 or 3
    if (type == 0 || type == 3) {
      /* Ensure that a valid tet is used in init_linear_id_with_level */
      p->pyramid.type = type;
      p->pyramid.level = i;
      p->switch_shape_at_level = i - 1;
      t8_dtet_init_linear_id_with_level (&(p->pyramid), id, i, level, type);
      T8_ASSERT (p->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (p));
      return;
    }
    /* The local index depends on the alternating number of predecessors caused by switching between pyramids and 
     * tetrahedrons, which have a different number of children.*/
    const t8_linearidx_t local_index = t8_dpyramid_update_index (&id, type, 2 * p_sum1 - p_sum2, p_sum1);
    const t8_dpyramid_cube_id_t cube_id = t8_dpyramid_parenttype_Iloc_to_cid[type][local_index];
    T8_ASSERT (cube_id >= 0);
    /* Set the element in its cube */
    p->pyramid.x |= (cube_id % 2 == 1) ? 1 << offset_expo : 0;
    p->pyramid.y |= (cube_id == 2 || cube_id == 3 || cube_id == 6 || cube_id == 7) ? 1 << offset_expo : 0;
    p->pyramid.z |= (cube_id > 3) ? 1 << offset_expo : 0;
    /*Compute the type */
    type = t8_dpyramid_parenttype_Iloc_to_type[type][local_index];
    T8_ASSERT (type >= 0);
  }
  T8_ASSERT (id == 0);
  p->pyramid.type = type;
  /* All ancestors are pyramids, the shape does not change. */
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    p->switch_shape_at_level = level;
    T8_ASSERT (p->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (p));
  }
  else {
    p->switch_shape_at_level = -1;
  }
}

int
t8_dpyramid_type_at_level (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);

  if (level >= p->pyramid.level) {
    return p->pyramid.type;
  }
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID || level >= p->switch_shape_at_level) {
    /* The shape does not switch, we can use the compute_type_same_shape function */
    return compute_type_same_shape (p, level);
  }
  else {
    /* The shape switches */
    T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
    t8_dpyramid_t ancestor;
    t8_dpyramid_ancestor (p, p->switch_shape_at_level, &ancestor);
    t8_dpyramid_parent (&ancestor, &ancestor);
    if (level == ancestor.pyramid.level) {
      /* We have already reached the desired level and can return. */
      return ancestor.pyramid.type;
    }
    else {
      return compute_type_same_shape (&ancestor, level);
    }
  }
}

t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t *p, const int level)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  t8_linearidx_t id = 0, sum_1 = 1, sum_2 = 1;
  t8_dpyramid_t parent, copy;

  t8_dpyramid_copy (p, &copy);
  copy.pyramid.type = t8_dpyramid_type_at_level (p, level);
  copy.pyramid.level = level;
  t8_dpyramid_cut_coordinates (&copy, T8_DPYRAMID_MAXLEVEL - level);

  for (int i = level; i > 0; i--) {
    /* Compute the number of pyramids with level maxlvl that are in a pyramid
     * of level i*/
    const t8_linearidx_t pyra_shift = (sum_1 << 1) - sum_2;

    /*Compute the parent and the local id of the current element */
    t8_dpyramid_parent (&copy, &parent);
    const t8_linearidx_t local_id = t8_dpyramid_child_id (&copy);

    /* Compute the number of predecessors within the parent that have the
     * shape of a pyramid or a tet*/
    int num_pyra;
    if (t8_dpyramid_shape (&parent) == T8_ECLASS_TET) {
      /* If the parent is a tet, no predecessors are pyramids */
      num_pyra = 0;
    }
    else {
      /* The number of pyramid-predecessors */
      num_pyra = t8_dpyramid_parenttype_iloc_pyra_w_lower_id[parent.pyramid.type - T8_DPYRAMID_FIRST_TYPE][local_id];
    }
    /* The number of tets is the local-id minus the number of pyramid-predecessors */
    const int num_tet = local_id - num_pyra;
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
 * \param[in] p               Input element
 * \param[in] face            A face of \a p
 * \param[in, out] neigh      Allocated memory, will be filled with the data of the neighbor of \a p along the given \a face.
 * \return int                The face of \a neigh that touches \a p
 */
static int
t8_dpyramid_face_neighbour (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *neigh)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->pyramid.level);
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
      neigh->pyramid.type
        = ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? T8_DPYRAMID_SECOND_TYPE : T8_DPYRAMID_FIRST_TYPE);
    }
    /*Compute the coords of the neighbour */
    /*Do nothing for face == 0 || face == 2 */
    if (face == 1) {
      neigh->pyramid.x += ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
      neigh->pyramid.y += ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
    }
    else if (face == 3) {
      neigh->pyramid.x += ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? 0 : -length);
      neigh->pyramid.y += ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? length : 0);
    }
    else if (face == 4) {
      neigh->pyramid.z += ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE) ? -length : length);
    }
    return t8_dpyramid_type_face_to_nface[p->pyramid.type - T8_DPYRAMID_FIRST_TYPE][face];
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

/** Compute the pyramid-parent-type of a tetrahedron
 * \param [in] p        Input pyramid
 * \return              The type of the parent.
 */
static int
t8_dpyramid_tetparent_type (const t8_dpyramid_t *p)
{
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  if ((p->pyramid.z >> (T8_DPYRAMID_MAXLEVEL - p->pyramid.level)) % 2 == 0) {
    return T8_DPYRAMID_FIRST_TYPE;
  }
  else {
    return T8_DPYRAMID_SECOND_TYPE;
  }
}

int
t8_dpyramid_face_parent_face (const t8_dpyramid_t *elem, const int face)
{
  t8_dpyramid_t parent;
  T8_ASSERT (0 <= elem->pyramid.level && elem->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  /*parent is a pyramid */
  if (elem->pyramid.level == 0) {
    return face;
  }
  if (t8_dpyramid_shape (elem) == T8_ECLASS_PYRAMID) {
    const int child_id = t8_dpyramid_child_id (elem);
    t8_dpyramid_parent (elem, &parent);

    /*If the pyramid is one of the children in the array, its face-num and the face-num
     * of the parent are the same*/
    for (int i = 0; i < 4; i++) {
      if (t8_dpyramid_type_face_to_children_at_face[parent.pyramid.type - T8_DPYRAMID_FIRST_TYPE][face][i]
          == child_id) {
        return face;
      }
    }
    /*No matching face */
    return -1;
  }
  else {
    const int child_id = t8_dpyramid_child_id (elem);
    /*Parent is also a tet, we can call the tet-routine */
    if (elem->pyramid.level > elem->switch_shape_at_level) {
      return t8_dtet_face_parent_face (&(elem->pyramid), face);
    }
    /*tet with a pyramid-parent */
    else {
      t8_dpyramid_type_t parent_type = t8_dpyramid_tetparent_type (elem);
      /*Only tets of type 0 or 3 have a pyra-parent.pyramid. Parent can have type 6 or 7 */
      if (elem->pyramid.type == 0 && parent_type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 3 && face == 1) {
          return 0;
        }
        else if (child_id == 5 && face == 0) {
          return 1;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 3 && parent_type == T8_DPYRAMID_FIRST_TYPE) {
        if (child_id == 1 && face == 1) {
          return 2;
        }
        else if (child_id == 6 && face == 0) {
          return 3;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 0 && parent_type == T8_DPYRAMID_SECOND_TYPE) {
        if (child_id == 1 && face == 3) {
          return 1;
        }
        else if (child_id == 7 && face == 2) {
          return 0;
        }
        else
          return -1;
      }
      else if (elem->pyramid.type == 3 && parent_type == T8_DPYRAMID_SECOND_TYPE) {
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
 *  Check if ancestor-face is connecting to a tet or to a pyra 
 * \param[in] p       Input pyramid
 * \param[in] face    Face of an ancestor of \a p
 * \return            return non-zero if ancestor-face is connecting to a tet
 */
static int
t8_dpyramid_tet_pyra_face_connection (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (p->pyramid.type == 0 || p->pyramid.type == 3);
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /* Depending on its cube-id at its level and its type,
   * 3 faces of a tet connect to a pyramid, one is connecting to a tet*/
  const t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (p->pyramid.type == 0 || p->pyramid.type == 3);
  t8_dpyramid_t ancestor;
  if (p->pyramid.level == p->switch_shape_at_level) {
    /*Check if the face is connecting to a tet or to a pyra */
    return t8_dpyramid_tet_pyra_face_connection (p, face);
  }
  t8_dpyramid_ancestor (p, p->switch_shape_at_level, &ancestor);
  /*Check if ancestor-face is connecting to a tet or to a pyra */
  int valid_touch = t8_dpyramid_tet_pyra_face_connection (&ancestor, face);
  if (valid_touch) {
    t8_dpyramid_type_t type_temp = p->pyramid.type;
    /* If so, check if the tet always lies in the corner of of its parent at this face.
     * Otherwise, the neighbor is a tet*/
    for (int i = p->pyramid.level; i > ancestor.pyramid.level; i--) {
      const t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, i);
      const int bey_id = t8_dtet_type_cid_to_beyid[type_temp][cube_id];
      if (t8_dpyramid_face_childid_to_is_inside[face][bey_id] == -1) {
        return 0;
      }
      type_temp = t8_dtet_cid_type_to_parenttype[cube_id][type_temp];
    }
  }
  /*Return if ancestor-face is connecting to a tet */
  return valid_touch;
}

int
t8_dpyramid_tree_face (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
      else if ((face == 1 && p->pyramid.type == 3) || (face == 2 && p->pyramid.type == 1)) {
        return 2;
      }
      else if ((face == 1 && p->pyramid.type == 0) || (face == 2 && p->pyramid.type == 2)) {
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*Check, if p is root pyramid */
  if (p->pyramid.level == 0) {
    return p->pyramid.type == T8_DPYRAMID_ROOT_TYPE && p->pyramid.x == 0 && p->pyramid.y == 0 && p->pyramid.z == 0;
  }
  /*Check, if all coordinates are in the limit set by the length of root */
  if ((0 <= p->pyramid.z) && (p->pyramid.z < T8_DPYRAMID_ROOT_LEN) && (p->pyramid.x >= p->pyramid.z)
      && (p->pyramid.x < T8_DPYRAMID_ROOT_LEN) && (p->pyramid.y >= p->pyramid.z)
      && (p->pyramid.y < T8_DPYRAMID_ROOT_LEN)) {
    if ((p->pyramid.x == p->pyramid.z && (p->pyramid.type == 3 || p->pyramid.type == 5))
        || (p->pyramid.y == p->pyramid.z && (p->pyramid.type == 0 || p->pyramid.type == 4))) {
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
t8_dpyramid_face_neighbor_inside (const t8_dpyramid_t *p, t8_dpyramid_t *neigh, const int face, int *neigh_face)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*Compute the face-neighbor, then check if it is inside root */
  *neigh_face = t8_dpyramid_face_neighbour (p, face, neigh);

  int is_inside = t8_dpyramid_is_inside_root (neigh);
  if (is_inside) {
    t8_dpyramid_set_switch_shape_at_level (neigh);
  }
  return is_inside;
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc, const int level)
{
  T8_ASSERT (level >= p->pyramid.level);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The first descendant of a pyramid has the same anchor coords, but another level */
    t8_dpyramid_copy (p, desc);
    desc->pyramid.level = level;
    desc->switch_shape_at_level = -1;
  }
  else {
    t8_linearidx_t id = t8_dpyramid_linear_id (p, level);
    t8_dpyramid_init_linear_id (desc, level, id);
  }
  T8_ASSERT (p->pyramid.x <= desc->pyramid.x && p->pyramid.y <= desc->pyramid.y && p->pyramid.z <= desc->pyramid.z);
}

/** Compute the descendant at a corner of a tet up to a given level. You can not use the tetrahedron algorithm here, 
 * because it depends on the linear-id computation of a tet, which is different to the linear id for pyramids. 
 *  \param[in] p        Input pyramd
 *  \param[in, out] d   Allocated element to fill with the data of element laying in the corner \a corner of \a p
 *  \param[in] corner   A corner of \a p
 *  \param[in] level    The level to compute the corner-descendant at. 
 */
static void
t8_dpyramid_corner_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *d, const int corner, const int level)
{
  T8_ASSERT (p->pyramid.level <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (t8_dpyramid_shape (p) == T8_ECLASS_TET);
  T8_ASSERT (0 <= corner && corner < T8_DTET_CORNERS);
  if (corner == 0) {
    t8_dpyramid_first_descendant (p, d, level);
  }
  else if (corner == 1 || corner == 2) {
    /*The child at this corner is iterativle the child at child-id up to the
       given level */
    const int child_id = t8_dtet_parenttype_beyid_to_Iloc[p->pyramid.type][corner];
    t8_dpyramid_t tmp;
    t8_dpyramid_copy (p, &tmp);
    for (int i = p->pyramid.level; i < level; i++) {
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
t8_dpyramid_first_descendant_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *first_desc, const int level)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->pyramid.level <= level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    const int corner = t8_dtet_face_corner[face][0];
    t8_dpyramid_corner_descendant (p, first_desc, corner, level);
  }
  else if (p->pyramid.level == T8_DPYRAMID_MAXLEVEL) {
    t8_dpyramid_copy (p, first_desc);
  }
  else {
    /*Shift the coordinate of p to compute the descendant */
    if ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && (face == 0 || face == 2 || face == 4))
        || (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && face != 4)) {
      /*No shifting is needed, fd is the first child with given level */
      t8_dpyramid_child (p, 0, first_desc);
      first_desc->pyramid.level = level;
    }
    else {
      /*shifting is needed, the fd does not have the same coord as p */
      t8_dpyramid_copy (p, first_desc);
      const t8_dpyramid_coord_t off_set = T8_DPYRAMID_LEN (p->pyramid.level) - T8_DPYRAMID_LEN (level);
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
t8_dpyramid_last_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc, const int level)
{
  T8_ASSERT (level >= p->pyramid.level);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_copy (p, desc);
    desc->pyramid.level = level;
    /* Shift the coords to the eights cube. The type of the last descendant is
     * is the type of the input pyramid*/
    t8_dpyramid_coord_t coord_offset = T8_DPYRAMID_LEN (p->pyramid.level) - T8_DPYRAMID_LEN (level);
    desc->pyramid.x |= coord_offset;
    desc->pyramid.y |= coord_offset;
    desc->pyramid.z |= coord_offset;
  }
  else {
    /* Compute current id, shift it to the id of the last desendant and compute it
     * via its id*/
    const t8_linearidx_t t_id = t8_dpyramid_linear_id (p, level);
    const int exponent = level - p->pyramid.level;
    t8_linearidx_t id = (((t8_linearidx_t) 1) << 3 * exponent) - 1;
    id += t_id;
    t8_dpyramid_init_linear_id (desc, level, id);
  }
}

void
t8_dpyramid_last_descendant_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *last_desc, const int level)
{
  /*Computation is similar to first-descendant-face */
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p->pyramid.level <= level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    T8_ASSERT (0 <= face && face < T8_DTET_FACES);
    const int corner = SC_MAX (t8_dtet_face_corner[face][1], t8_dtet_face_corner[face][2]);
    t8_dpyramid_corner_descendant (p, last_desc, corner, level);
  }
  else {
    const t8_dpyramid_coord_t off_set = T8_DPYRAMID_LEN (p->pyramid.level) - T8_DPYRAMID_LEN (level);
    t8_dpyramid_copy (p, last_desc);
    last_desc->pyramid.level = level;
    if ((p->pyramid.type == T8_DPYRAMID_FIRST_TYPE && face != 4)
        || (p->pyramid.type == T8_DPYRAMID_SECOND_TYPE && (face == 0 || face == 2 || face == 4))) {
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID || p->switch_shape_at_level == p->pyramid.level) {
    /* The parent is a pyramid. */
    return T8_DPYRAMID_CHILDREN;
  }
  else {
    return T8_DTET_CHILDREN;
  }
}

int
t8_dpyramid_num_faces (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_FACES;
  }
  else {
    return T8_DPYRAMID_FACES;
  }
}

void
t8_dpyramid_boundary_face (const t8_dpyramid_t *p, const int face, t8_element_t *boundary)
{
  /* face is face of of p */
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (face == 4) {
    /*On the bottom every face is a quad */
    /*Coordinates are scaled, because quad and pyra might have different root-len */
    p4est_quadrant_t *q = (p4est_quadrant_t *) boundary;
    q->x = ((int64_t) p->pyramid.x * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->y = ((int64_t) p->pyramid.y * P4EST_ROOT_LEN) / T8_DPYRAMID_ROOT_LEN;
    q->level = p->pyramid.level;
  }
  else {
    /* Boundary-face is a triangle. */
    /* p-x give t->x for root-face 2,3 and p->pyramid.y gives t->x for 0,1.
     * t->y is determined by p->pyramid.z*/
    t8_dtri_t *t = (t8_dtri_t *) boundary;
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
      /*Boundary is given by a tet-surface. The cases are ordered by root-face-enumeration*/
      if ((face == 1 && p->pyramid.type == 0) || (face == 2 && p->pyramid.type == 2)) {
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 0 ? 1 : 0;
      }
      else if (face == 0 && (p->pyramid.type == 0 || p->pyramid.type == 1)) {
        t->x = p->pyramid.y * T8_DTRI_ROOT_BY_DPYRAMID_ROOT;
        t->type = p->pyramid.type == 0 ? 1 : 0;
      }
      else if ((face == 1 && p->pyramid.type == 3) || (face == 2 && p->pyramid.type == 1)) {
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
t8_dpyramid_extrude_face (const t8_element_t *face, t8_dpyramid_t *p, const int root_face)
{
  T8_ASSERT (0 <= root_face && root_face < T8_DPYRAMID_FACES);

  int extruded_face;
  if (root_face == 4) {
    /* Pyramids on the bottom are always type 6 pyramids at the bottom. We need to
     * scale the coordinates, since a quad and a pyra can have different root-len,
     * depending on their maxlvl.*/
    p4est_quadrant_t *q = (p4est_quadrant_t *) face;
    /*Typecast to int64, we multiply two (possible at max) int32 */
    p->pyramid.x = ((int64_t) q->x * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->pyramid.y = ((int64_t) q->y * T8_DPYRAMID_ROOT_LEN) / P4EST_ROOT_LEN;
    p->pyramid.z = 0;
    p->pyramid.type = T8_DPYRAMID_ROOT_TYPE;
    p->pyramid.level = q->level;
    p->switch_shape_at_level = -1;
    return root_face;
  }
  else {
    /*t->y gives the height of the pyramid, t->x gives the p->pyramid.x or p->pyramid.y, depending on
     * the root_face. The other coordinate is determined by the root_face.*/
    t8_dtri_t *t = (t8_dtri_t *) face;
    p->pyramid.z = ((int64_t) t->y * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    /* level is the same */
    p->pyramid.level = t->level;
    switch (root_face) {
    case 0:
      p->pyramid.x = p->pyramid.z;
      /*Typecast to int64, we multiply two (possible at max) int32 */
      p->pyramid.y = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 1:
      p->pyramid.x = T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
      p->pyramid.y = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      break;
    case 2:
      p->pyramid.x = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->pyramid.y = p->pyramid.z;
      break;
    case 3:
      p->pyramid.x = ((int64_t) t->x * T8_DPYRAMID_ROOT_LEN) / T8_DTRI_ROOT_LEN;
      p->pyramid.y = T8_DPYRAMID_ROOT_LEN - T8_DPYRAMID_LEN (p->pyramid.level);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    /*Description of triangles extruding to a pyramid */
    if ((t->y == (t->x & t->y)) && t->type == 0) {
      /*type zero in a pyramid extend to a pyramid */
      p->pyramid.type = T8_DPYRAMID_FIRST_TYPE;
      p->switch_shape_at_level = -1;
      extruded_face = root_face;
    }
    else {
      /*type 0 not in a pyramid extend to a tetrahedron */
      p->pyramid.type = t8_dpyramid_tritype_rootface_to_tettype[t->type][root_face];
      extruded_face = t8_dpyramid_tritype_rootface_to_face[t->type][root_face];
      t8_dpyramid_set_switch_shape_at_level (p);
    }
  }
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*return the face-number of the extruded face */
  return extruded_face;
}

int
t8_dpyramid_child_id (const t8_dpyramid_t *p)
{

  T8_ASSERT (p->pyramid.level >= 0);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    T8_ASSERT (p->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (p));
  }
  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID || p->switch_shape_at_level == p->pyramid.level) {
    if (p->pyramid.level == 0) {
      return 0;
    }
    const t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);
    T8_ASSERT (t8_dpyramid_type_cid_to_Iloc[p->pyramid.type][cube_id] >= 0);
    return t8_dpyramid_type_cid_to_Iloc[p->pyramid.type][cube_id];
  }
  else {
    return t8_dtet_child_id (&(p->pyramid));
  }
}

void
t8_dpyramid_child (const t8_dpyramid_t *elem, const int child_id, t8_dpyramid_t *child)
{
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
  T8_ASSERT (0 <= elem->pyramid.level && elem->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (elem) == T8_ECLASS_TET) {
    t8_dtet_child (&(elem->pyramid), child_id, &(child->pyramid));
    child->switch_shape_at_level = elem->switch_shape_at_level;
  }
  else {
    /* Compute the cube id and shift the coordinates accordingly */
    const t8_dpyramid_cube_id_t cube_id = t8_dpyramid_parenttype_Iloc_to_cid[elem->pyramid.type][child_id];
    const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (elem->pyramid.level + 1);
    T8_ASSERT (cube_id >= 0);
    child->pyramid.level = elem->pyramid.level + 1;
    child->pyramid.x = elem->pyramid.x + ((cube_id & 0x01) ? length : 0);
    child->pyramid.y = elem->pyramid.y + ((cube_id & 0x02) ? length : 0);
    child->pyramid.z = elem->pyramid.z + ((cube_id & 0x04) ? length : 0);
    child->pyramid.type = t8_dpyramid_parenttype_Iloc_to_type[elem->pyramid.type][child_id];
    if (t8_dpyramid_shape (child) == T8_ECLASS_TET) {
      /* Use the level to set switch_shape_at_level, because the function has to be callable
       * which elem = child */
      child->switch_shape_at_level = child->pyramid.level;
    }
    else {
      child->switch_shape_at_level = -1;
    }
  }
  T8_ASSERT (child->pyramid.type >= 0);
}

void
t8_dpyramid_children (const t8_dpyramid_t *p, t8_dpyramid_t **c)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  const int num_children = t8_dpyramid_num_children (p);
  for (int i = num_children - 1; i >= 0; i--) {
    t8_dpyramid_child (p, i, c[i]);
  }
}

void
t8_dpyramid_children_at_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *children[], const int num_children,
                              int *child_indices)
{
  T8_ASSERT (num_children == T8_DPYRAMID_FACE_CHILDREN);
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    /* Use tet-algo */
    t8_dtet_t face_children[T8_DTET_FACE_CHILDREN];
    t8_dtet_t *tet_children[T8_DTET_FACE_CHILDREN]
      = { &face_children[0], &face_children[1], &face_children[2], &face_children[3] };

    t8_dtet_children_at_face (&(p->pyramid), face, tet_children, num_children, child_indices);
    for (int i = 0; i < T8_DTET_FACE_CHILDREN; i++) {
      t8_dtet_copy (tet_children[i], &(children[i]->pyramid));
      children[i]->switch_shape_at_level = p->switch_shape_at_level;
    }
  }
  else {
    int *children_at_face_id, children_at_face_id_local[T8_DPYRAMID_FACE_CHILDREN];
    if (child_indices != NULL) {
      children_at_face_id = child_indices;
    }
    else {
      children_at_face_id = children_at_face_id_local;
    }
    /*Fill the child ids with the child-ids at the face */
    for (int i = 0; i < T8_DPYRAMID_FACE_CHILDREN; i++) {
      children_at_face_id[i]
        = t8_dpyramid_type_face_to_children_at_face[p->pyramid.type - T8_DPYRAMID_FIRST_TYPE][face][i];
    }
    /*Compute the children */
    for (int i = T8_DPYRAMID_FACE_CHILDREN - 1; i >= 0; i--) {
      t8_dpyramid_child (p, children_at_face_id[i], children[i]);
    }
  }
}

int
t8_dpyramid_face_child_face (const t8_dpyramid_t *p, const int face, const int face_child)
{
  T8_ASSERT (0 <= face && face < T8_DPYRAMID_FACES);
  T8_ASSERT (0 <= face_child && face_child < T8_DPYRAMID_CHILDREN);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return t8_dtet_face_child_face (&(p->pyramid), face, face_child);
  }
  else {
    int child_face = t8_dpyramid_type_face_to_child_face[p->pyramid.type - T8_DPYRAMID_FIRST_TYPE][face][face_child];
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
    const int corner_number = t8_dpyramid_face_corner[pyra->pyramid.type - T8_DPYRAMID_FIRST_TYPE][face][corner];
    T8_ASSERT (0 <= corner_number && corner_number < T8_DPYRAMID_FACES);
    return corner_number;
  }
}

void
t8_dpyramid_parent (const t8_dpyramid_t *p, t8_dpyramid_t *parent)
{
  T8_ASSERT (p->pyramid.level > 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);
  const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p->pyramid.level);

  if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
    /*The parent of a pyramid is a pyramid, maybe of different type */
    const t8_dpyramid_cube_id_t cube_id = compute_cubeid (p, p->pyramid.level);

    parent->pyramid.type = t8_dpyramid_type_cid_to_parenttype[p->pyramid.type - T8_DPYRAMID_FIRST_TYPE][cube_id];
    parent->pyramid.x = p->pyramid.x & ~length;
    parent->pyramid.y = p->pyramid.y & ~length;
    parent->pyramid.z = p->pyramid.z & ~length;
    T8_ASSERT (parent->pyramid.type >= 0);
    parent->pyramid.level = p->pyramid.level - 1;
    parent->switch_shape_at_level = -1;
  }
  else if (p->switch_shape_at_level == p->pyramid.level) {
    /*p does not lie in larger tet => parent is pyra */
    parent->pyramid.type = t8_dpyramid_tetparent_type (p);
    parent->pyramid.x = p->pyramid.x & ~length;
    parent->pyramid.y = p->pyramid.y & ~length;
    parent->pyramid.z = p->pyramid.z & ~length;
    parent->pyramid.level = p->pyramid.level - 1;
    parent->switch_shape_at_level = -1;
  }
  else {
    /* The direct tet-child of a pyra has type 0 or type 3, therefore
     * in this case the parent is a tetrahedron*/
    t8_dtet_parent (&(p->pyramid), &(parent->pyramid));
    parent->switch_shape_at_level = p->switch_shape_at_level;
  }
  T8_ASSERT (parent->pyramid.level >= 0);
}

t8_element_shape_t
t8_dpyramid_shape (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL);
  /*The pyramid has the shape of a tetrahedron */
  if (p->pyramid.type < T8_DPYRAMID_FIRST_TYPE) {
    return T8_ECLASS_TET;
  }
  else {
    return T8_ECLASS_PYRAMID;
  }
}

static void
t8_dpyramid_successor_recursion (const t8_dpyramid_t *elem, t8_dpyramid_t *succ, const int level)
{
  T8_ASSERT (1 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_copy (elem, succ);
  succ->pyramid.level = level;
  if (level < succ->switch_shape_at_level) {
    succ->switch_shape_at_level = -1;
  }
  T8_ASSERT (succ->pyramid.type >= 0);
  const int child_id = t8_dpyramid_child_id (elem);
  /*Compute number of children */
  const int num_siblings = t8_dpyramid_num_siblings (elem);
  T8_ASSERT (0 <= child_id && child_id < num_siblings);
  if (child_id == num_siblings - 1) {
    int shift = T8_DPYRAMID_MAXLEVEL - level + 1;
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion (succ, succ, level - 1);
    succ->pyramid.level = level;
    /* bits auf level auf child 0 setzen */
    t8_dpyramid_cut_coordinates (succ, shift);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_dpyramid_parent (succ, succ);
    t8_dpyramid_child (succ, child_id + 1, succ);
    succ->pyramid.level = level;
  }
}

void
t8_dpyramid_successor (const t8_dpyramid_t *elem, t8_dpyramid_t *succ, const int level)
{
  t8_dpyramid_successor_recursion (elem, succ, level);
#if T8_ENABLE_DEBUG
  if (t8_dpyramid_shape (succ) == T8_ECLASS_PYRAMID) {
    T8_ASSERT (succ->switch_shape_at_level < 0);
  }
  else {
    T8_ASSERT (succ->switch_shape_at_level = t8_dpyramid_compute_switch_shape_at_level (succ));
  }
#endif
}

void
t8_dpyramid_compute_integer_coords (const t8_dpyramid_t *elem, const int vertex, int coords[])
{
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);

  if (t8_dpyramid_shape (elem) == T8_ECLASS_PYRAMID) {
    const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (elem->pyramid.level);
    coords[0] = elem->pyramid.x;
    coords[1] = elem->pyramid.y;
    coords[2] = elem->pyramid.z;
    switch (vertex) {
    case 0:
      if (elem->pyramid.type == T8_DPYRAMID_SECOND_TYPE)
        coords[2] += length;
      break;
    case 1:
      coords[0] += length;
      if (elem->pyramid.type == T8_DPYRAMID_SECOND_TYPE)
        coords[2] += length;
      break;
    case 2:
      coords[1] += length;
      if (elem->pyramid.type == T8_DPYRAMID_SECOND_TYPE)
        coords[2] += length;
      break;
    case 3:
      coords[0] += length;
      coords[1] += length;
      if (elem->pyramid.type == T8_DPYRAMID_SECOND_TYPE)
        coords[2] += length;
      break;
    case 4:
      if (elem->pyramid.type == T8_DPYRAMID_FIRST_TYPE) {
        coords[0] += length;
        coords[1] += length;
        coords[2] += length;
      }
      break;
    }
  }
  else {
    T8_ASSERT (t8_dpyramid_shape (elem) == T8_ECLASS_TET);
    T8_ASSERT (0 <= vertex && vertex < T8_DTET_CORNERS);
    t8_dtet_compute_integer_coords (&(elem->pyramid), vertex, coords);
  }
}

void
t8_dpyramid_vertex_reference_coords (const t8_dpyramid_t *elem, const int vertex, double coords[])
{
  int coords_int[3];
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);
  t8_dpyramid_compute_integer_coords (elem, vertex, coords_int);
  /*scale the coordinates onto the reference cube */
  coords[0] = coords_int[0] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[1] = coords_int[1] / (double) T8_DPYRAMID_ROOT_LEN;
  coords[2] = coords_int[2] / (double) T8_DPYRAMID_ROOT_LEN;
}

void
t8_dpyramid_compute_reference_coords (const t8_dpyramid_t *elem, const double *ref_coords, const size_t num_coords,
                                      double *out_coords)
{
  T8_ASSERT (ref_coords != NULL);
  T8_ASSERT (t8_dpyramid_is_valid (elem));
  if (t8_dpyramid_shape (elem) == T8_ECLASS_PYRAMID) {
    const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (elem->pyramid.level);
    size_t coord;
    for (coord = 0; coord < num_coords; ++coord) {
      const size_t offset = coord * 3;
      out_coords[offset + 0] = elem->pyramid.x;
      out_coords[offset + 1] = elem->pyramid.y;
      out_coords[offset + 2] = elem->pyramid.z;

      out_coords[offset + 0] += ref_coords[offset + 0] * length;
      out_coords[offset + 1] += ref_coords[offset + 1] * length;
      out_coords[offset + 2] += ref_coords[offset + 2] * length;
    }
    if (elem->pyramid.type == T8_DPYRAMID_SECOND_TYPE) {
      for (coord = 0; coord < num_coords; ++coord) {
        const size_t offset = coord * 3;
        out_coords[offset + 0] -= ref_coords[offset + 2] * length;
        out_coords[offset + 1] -= ref_coords[offset + 2] * length;
        out_coords[offset + 2] += (1 - 2 * ref_coords[offset + 2]) * length;
      }
    }
    for (coord = 0; coord < num_coords; ++coord) {
      const size_t offset = coord * 3;
      /* Scale the coordinates onto the reference cube */
      out_coords[offset + 0] /= (double) T8_DPYRAMID_ROOT_LEN;
      out_coords[offset + 1] /= (double) T8_DPYRAMID_ROOT_LEN;
      out_coords[offset + 2] /= (double) T8_DPYRAMID_ROOT_LEN;
    }
  }
  else {
    t8_dtet_compute_reference_coords (&(elem->pyramid), ref_coords, num_coords, out_coords);
  }
}

/**
 * Compute the ancestor of \a pyra on level \a level.
 * \param[in]       pyra    Input pyramid
 * \param[in]       level   The level at which we want to compute \a ancestoranc
 * \param[in, out]  ancestor     Allocated input element which will be filled by the data of the ancestor of \a pyra at level \a level
 */
void
t8_dpyramid_ancestor (const t8_dpyramid_t *pyra, const int level, t8_dpyramid_t *ancestor)
{
  T8_ASSERT (0 <= level && level <= pyra->pyramid.level);
  /*Set the coordinates and the level of the ancestor */
  t8_dpyramid_copy (pyra, ancestor);
  if (pyra->pyramid.level == level) {
    return;
  }
  else if (level == pyra->pyramid.level - 1) {
    /* We can reuse the parent computation if we want to go only one level up. */
    t8_dpyramid_parent (pyra, ancestor);
    return;
  }
  /* The coordinates of the ancestor are defined by the level. */
  t8_dpyramid_cut_coordinates (ancestor, T8_DPYRAMID_MAXLEVEL - level);
  ancestor->pyramid.level = level;
  ancestor->pyramid.type = t8_dpyramid_type_at_level (pyra, level);
  if (t8_dpyramid_shape (ancestor) == T8_ECLASS_TET) {
    ancestor->switch_shape_at_level = pyra->switch_shape_at_level;
  }
  else {
    ancestor->switch_shape_at_level = -1;
  }
}

void
t8_dpyramid_nearest_common_ancestor (const t8_dpyramid_t *pyra1, const t8_dpyramid_t *pyra2, t8_dpyramid_t *nca)
{
  /* If the input elements have different shapes, the nca has to have the
   * shape of a pyramid. The element in the shape of a tet switches the shape. */
  if (t8_dpyramid_shape (pyra1) == T8_ECLASS_PYRAMID && t8_dpyramid_shape (pyra2) == T8_ECLASS_TET) {
    t8_dpyramid_t first_pyramid_ancestor;

    t8_dpyramid_ancestor (pyra2, pyra2->switch_shape_at_level - 1, &first_pyramid_ancestor);

    /* pyra1 and first_pyramid_ancestor have the shape of a pyramid now, 
     * we can call the nca again.
     */
    t8_dpyramid_nearest_common_ancestor (pyra1, &first_pyramid_ancestor, nca);
    return;
  }
  else if (t8_dpyramid_shape (pyra1) == T8_ECLASS_TET && t8_dpyramid_shape (pyra2) == T8_ECLASS_PYRAMID) {
    t8_dpyramid_t first_pyramid_ancestor;

    t8_dpyramid_ancestor (pyra1, pyra1->switch_shape_at_level - 1, &first_pyramid_ancestor);

    /* pyra2 and first_pyramid_ancestor have the shape of a pyramid now, 
     * we can call the nca again.
     */
    t8_dpyramid_nearest_common_ancestor (&first_pyramid_ancestor, pyra2, nca);
    return;
  }
  /* both elements have the shape of a pyramid, hence the nca */
  else if (t8_dpyramid_shape (pyra1) == T8_ECLASS_PYRAMID && t8_dpyramid_shape (pyra2) == T8_ECLASS_PYRAMID) {
    /* The following computations are necessary to find the irst ancestors of pyra1 and pyra2 with the same type. We 
     * have already computed the level at which they have the same coordinate, but the type could be different. */
    t8_dpyramid_coord_t maxclor;
    t8_dpyramid_type_t p1_type_at_level; /* type of pyra1 at level */
    t8_dpyramid_type_t p2_type_at_level; /* type of pyra2 at level */
    /* Compute the first level, at which the coordinates differ */
    maxclor = pyra1->pyramid.x ^ pyra2->pyramid.x;
    maxclor |= pyra1->pyramid.y ^ pyra2->pyramid.y;
    maxclor |= pyra1->pyramid.z ^ pyra2->pyramid.z;
    const int level = SC_LOG2_32 (maxclor) + 1;
    T8_ASSERT (level <= T8_DPYRAMID_MAXLEVEL);
    /* This is the highest possible level. The coordinates are the same,
     * but the types can be different.*/
    /* the level of the cube where pyra1 and pyra2 have the same coords */
    const int cube_level
      = SC_MIN (T8_DPYRAMID_MAXLEVEL - level, (int) SC_MIN (pyra1->pyramid.level, pyra2->pyramid.level));
    /* the level of the nca */
    int real_level = cube_level;
    p1_type_at_level = compute_type_same_shape (pyra1, cube_level);
    p2_type_at_level = compute_type_same_shape (pyra2, cube_level);
    /* Iterate over the levels and compute both types at that level.
     * If they are the same, we know the level of the nearest common ancestor. */
    while (p1_type_at_level != p2_type_at_level) {
      real_level--;
      p1_type_at_level = compute_type_same_shape_ext (pyra1, real_level, p1_type_at_level, real_level + 1);
      p2_type_at_level = compute_type_same_shape_ext (pyra2, real_level, p2_type_at_level, real_level + 1);
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
    /* Both elements are a tet. The ancestor can be at a level before any of the
     * elementes switches the shape from a tet to a pyra. If one of the tets switches
     * the shape, both tets have to switch the shape. */
    T8_ASSERT (t8_dpyramid_shape (pyra1) == T8_ECLASS_TET);
    T8_ASSERT (t8_dpyramid_shape (pyra2) == T8_ECLASS_TET);
    t8_dpyramid_coord_t maxclor;
    t8_dpyramid_type_t p1_type_at_level; /* type of pyra1 at level */
    t8_dpyramid_type_t p2_type_at_level; /* type of pyra2 at level */
    /* Compute the first level, at which the coordinates differ */
    maxclor = pyra1->pyramid.x ^ pyra2->pyramid.x;
    maxclor |= pyra1->pyramid.y ^ pyra2->pyramid.y;
    maxclor |= pyra1->pyramid.z ^ pyra2->pyramid.z;
    const int level = SC_LOG2_32 (maxclor) + 1;
    T8_ASSERT (level <= T8_DPYRAMID_MAXLEVEL);
    t8_dpyramid_t pyra1_ancestor;
    t8_dpyramid_t pyra2_ancestor;

    /* Cube level is the highest possible level where the nca can be. the coordinates
     * match at the level, but the type can be different.*/
    /* the level of the cube where pyra1 and pyra2 have the same coords */
    const int cube_level
      = SC_MIN (T8_DPYRAMID_MAXLEVEL - level, (int) SC_MIN (pyra1->pyramid.level, pyra2->pyramid.level));
    int real_level = cube_level; /* the level of the nca */
    t8_dpyramid_ancestor (pyra1, real_level, &pyra1_ancestor);
    t8_dpyramid_ancestor (pyra2, real_level, &pyra2_ancestor);

    p1_type_at_level = pyra1_ancestor.pyramid.type;
    p2_type_at_level = pyra2_ancestor.pyramid.type;
    T8_ASSERT (pyra1->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (pyra1));
    T8_ASSERT (pyra2->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (pyra2));

    /* We iterate over the levels and check if the types of both tets match and stop at that level.
     * The loop is interrupted, if we get to a level where one element switches the shape.*/
    T8_ASSERT (0 < pyra1->switch_shape_at_level && pyra1->switch_shape_at_level <= T8_DPYRAMID_MAXLEVEL);
    T8_ASSERT (0 < pyra2->switch_shape_at_level && pyra2->switch_shape_at_level <= T8_DPYRAMID_MAXLEVEL);
    while (p1_type_at_level != p2_type_at_level && real_level >= pyra1->switch_shape_at_level
           && real_level >= pyra2->switch_shape_at_level) {
      real_level--;
      p1_type_at_level = compute_type_same_shape_ext (pyra1, real_level, p1_type_at_level, real_level + 1);
      p2_type_at_level = compute_type_same_shape_ext (pyra2, real_level, p2_type_at_level, real_level + 1);
    }
    if (real_level < pyra1->switch_shape_at_level || real_level < pyra2->switch_shape_at_level) {
      t8_dpyramid_ancestor (pyra1, pyra1->switch_shape_at_level - 1, &pyra1_ancestor);
      t8_dpyramid_ancestor (pyra2, pyra2->switch_shape_at_level - 1, &pyra2_ancestor);
      t8_dpyramid_nearest_common_ancestor (&pyra1_ancestor, &pyra2_ancestor, nca);
    }
    else {
      T8_ASSERT (p1_type_at_level == p2_type_at_level);
      /* The nearest common ancestor is a tetrahedron. */
      t8_dtet_ancestor (&(pyra1->pyramid), real_level, &(nca->pyramid));
      nca->switch_shape_at_level = pyra1->switch_shape_at_level;
    }
  }
}

int
t8_dpyramid_is_valid (const t8_dpyramid_t *p)
{
  int is_valid;
  const t8_dpyramid_coord_t max_coord = ((int64_t) 2 * T8_DPYRAMID_ROOT_LEN) - 1;
  const t8_element_shape_t shape = t8_dpyramid_shape (p);
  /*Check the level */
  is_valid = 0 <= p->pyramid.level && p->pyramid.level <= T8_DPYRAMID_MAXLEVEL;
  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.x && p->pyramid.x <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.y && p->pyramid.y <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p->pyramid.z && p->pyramid.z <= max_coord;

  /*The shape can be a pyramid or a tet */
  is_valid = is_valid && (shape == T8_ECLASS_PYRAMID || shape == T8_ECLASS_TET);
  /*Check the type */
  is_valid = is_valid && 0 <= p->pyramid.type && p->pyramid.type < T8_DPYRAMID_NUM_TYPES;

  if (p->pyramid.level == 0) {
    is_valid = is_valid && (p->pyramid.type == T8_DPYRAMID_ROOT_TYPE);
  }
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    is_valid = is_valid && (p->switch_shape_at_level > 0);
    is_valid = is_valid && (p->switch_shape_at_level <= T8_DPYRAMID_MAXLEVEL);
    is_valid = is_valid && (p->switch_shape_at_level == t8_dpyramid_compute_switch_shape_at_level (p));
  }
  else {
    is_valid = is_valid && (p->switch_shape_at_level < 0);
  }

  return is_valid;
}
