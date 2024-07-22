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

#include "t8_sele_bits_cxx.hxx"
#include <iostream>

/*********GETTER******/

/** The length of a element at a given level in integer coordinates */
template <t8_eclass_t eclass_T>
t8_element_coord_t
t8_sele_get_len (t8_element_level_t level)
{
  return 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - (level));
}

template <t8_eclass_t eclass_T>
t8_element_coord_t
t8_sele_get_root_len ()
{
  if constexpr (eclass_T == T8_ECLASS_VERTEX) {
    return 0;
  }
  else {
    return 1 << T8_ELEMENT_MAXLEVEL[eclass_T];
  }
}

template <t8_eclass_t eclass_T>
t8_element_shape_t
t8_sele_shape (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (p->type == T8_DPYRAMID_FIRST_PYRA_TYPE || p->type == T8_DPYRAMID_SECOND_PYRA_TYPE) {
      return T8_ECLASS_PYRAMID;
    }
    else {
      return T8_ECLASS_TET;
    }
  }
  return eclass_T;
}

template <t8_eclass_t eclass_T>
int
t8_sele_num_corners (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*TODO*/
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (t8_sele_shape (p) == T8_ECLASS_PYRAMID) {
      return 5;
    }
    return 4;
  }
  return T8_ELEMENT_NUM_CORNERS[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_sele_num_children (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*TODO*/
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (t8_sele_shape (p) == T8_ECLASS_PYRAMID) {
      return 10;
    }
    return 6;
  }
  return T8_ELEMENT_NUM_CHILDREN[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_sele_num_siblings (const t8_standalone_element_t<eclass_T> *p)
{
  if (p->level == 0)
    return 1;
  T8_ASSERT (0 < p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_standalone_element_t<eclass_T> parent;
  t8_sele_parent (p, &parent);
  return t8_sele_num_children (&parent);
}

template <t8_eclass_t eclass_T>
int
t8_sele_num_faces (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    if (t8_sele_shape (p) == T8_ECLASS_PYRAMID) {
      return 5;
    }
    return 4;
  }
  return T8_ELEMENT_NUM_FACES[eclass_T];
}

template <t8_eclass_t eclass_T>
int
t8_sele_max_num_faces (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  return T8_ELEMENT_NUM_FACES[eclass_T];
}

/** Cube helper **/

template <t8_eclass_t eclass_T>
static t8_cube_id_t
t8_sele_compute_cubeid (const t8_standalone_element_t<eclass_T> *p, const int level)
{
  t8_cube_id_t cube_id = 0;

  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  const t8_element_coord_t h = t8_sele_get_len<eclass_T> (level);

  if (level == 0) {
    return 0;
  }
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    cube_id |= ((p->coords[i] & h) ? 1 << i : 0);
  }
  return cube_id;
}

/* For each typebit, consider the coordinate information between level and p->level |10...11|xxxx|0...0| 
 * of both inequality defining dimensions */
template <t8_eclass_t eclass_T>
t8_element_type_t<eclass_T>
t8_sele_compute_type_at_level (const t8_standalone_element_t<eclass_T> *p, int level)
{
  t8_element_type_t<eclass_T> type = 0;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_element_coord_t coord_v0 = p->coords[t8_type_edge_equations<eclass_T>[e][0]];
    t8_element_coord_t coord_v1 = p->coords[t8_type_edge_equations<eclass_T>[e][1]];

    coord_v0 = (coord_v0 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);
    coord_v1 = (coord_v1 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);

    if (coord_v0 == coord_v1) {
      type[e] = p->type[e] | type[e];
    }
    else if (coord_v0 < coord_v1) {
      type |= (1 << e);
    }
    else {
      T8_ASSERT (coord_v0 > coord_v1);
      T8_ASSERT (!(type & (t8_element_type_t<eclass_T>) (1 << e)).all ());
    }
  }
  return type;
}

/**
 * Set the \a shift last bits of every coordinate to zero. 
 * 
 * \param[in, out]  p     Input element
 * \param[in]       shift Number of bits to set to zero
 */
template <t8_eclass_t eclass_T>
static inline void
t8_sele_cut_coordinates (t8_standalone_element_t<eclass_T> *p, const int shift)
{
  T8_ASSERT (0 <= shift && shift <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    p->coords[i] = (p->coords[i] >> shift) << shift;
  }
}

/*Copies a element from source to dest*/
template <t8_eclass_t eclass_T>
void
t8_sele_copy (const t8_standalone_element_t<eclass_T> *source, t8_standalone_element_t<eclass_T> *dest)
{
  T8_ASSERT (source != NULL && dest != NULL);
  if (source == dest) {
    return;
  }
  memcpy (dest, source, sizeof (t8_standalone_element_t<eclass_T>));
}

template <t8_eclass_t eclass_T>
int
t8_sele_compare (const t8_standalone_element_t<eclass_T> *p1, const t8_standalone_element_t<eclass_T> *p2)
{
  const int maxlvl = SC_MAX (p1->level, p2->level);

  const t8_linearidx_t id1 = t8_sele_linear_id (p1, maxlvl);
  const t8_linearidx_t id2 = t8_sele_linear_id (p2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the pyramid with the smaller level
     * is considered smaller */
    if (p1->level == p2->level) {
      return 0;
    }
    else {
      return p1->level - p2->level;
    }
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
     id2 */
  return id1 < id2 ? -1 : 1;
}

template <t8_eclass_t eclass_T>
int
t8_sele_is_valid (const t8_standalone_element_t<eclass_T> *p)
{
  int is_valid;
  const t8_element_coord_t max_coord = ((uint64_t) 2 * (uint64_t) t8_sele_get_root_len<eclass_T> ()) - 1;

  /*Check the level */
  is_valid = 0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T];
  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    is_valid = is_valid && -(int64_t) t8_sele_get_root_len<eclass_T> () <= p->coords[i] && p->coords[i] <= max_coord;
  }

  return is_valid;
}

template <t8_eclass_t eclass_T>
void
t8_sele_debug_print (const t8_standalone_element_t<eclass_T> *p)
{
  t8_debugf ("level: %i\n", p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    t8_debugf ("x_%i: %i \n", i, p->coords[i]);
  }
  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_debugf ("t_%i: %i \n", e, p->type[e]);
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_global_print (const t8_standalone_element_t<eclass_T> *p)
{
  t8_global_productionf ("level: %i", p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    t8_global_productionf ("x_%i: %i \n", i, p->coords[i]);
  }
  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_global_productionf ("t_%i: %i \n", e, p->type[e]);
  }
}

template <t8_eclass_t eclass_T>
static void
t8_sele_root (t8_standalone_element_t<eclass_T> *p)
{
  p->level = 0;
  for (size_t i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    p->coords[i] = 0;
  }
  p->type = 0;
}

/**SFC functionality*/
template <t8_eclass_t eclass_T>
int
t8_sele_child_id (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (p->level >= 0);
  if (p->level == 0) {
    return -1;
  }
  const t8_cube_id_t cube_id = t8_sele_compute_cubeid (p, p->level);
  int8_t child_id;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    child_id = t8_element_type_cubeid_to_Iloc<eclass_T>[p->type.to_ulong ()][cube_id];
  }
  else {
    child_id = cube_id;
  }
  return child_id;
}

template <t8_eclass_t eclass_T>
void
t8_sele_child (const t8_standalone_element_t<eclass_T> *elem, const int child_id,
               t8_standalone_element_t<eclass_T> *child)
{
  T8_ASSERT (0 <= child_id && child_id < T8_ELEMENT_NUM_CHILDREN[eclass_T]);
  T8_ASSERT (0 <= elem->level && elem->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  /* Compute the cube id and shift the coordinates accordingly */
  t8_cube_id_t cube_id;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    cube_id = t8_element_type_Iloc_to_childcubeid<eclass_T>[elem->type.to_ulong ()][child_id];
    child->type = t8_element_type_Iloc_to_childtype<eclass_T>[elem->type.to_ulong ()][child_id];
  }
  else {
    cube_id = child_id;
  }

  const t8_element_coord_t length = t8_sele_get_len<eclass_T> (elem->level + 1);

  // Auslagern also bithelper function?
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    child->coords[i] = elem->coords[i] + ((cube_id & 1 << i) ? length : 0);
  }

  child->level = elem->level + 1;
}

template <t8_eclass_t eclass_T>
int
t8_sele_is_family (t8_standalone_element_t<eclass_T> **fam)
{
  t8_standalone_element_t<eclass_T> parent, compare;
  /* Take the parent of the first element as baseline to compare against */
  t8_sele_parent (fam[0], &parent);
  const int num_children = t8_sele_num_children (&parent);
  for (int childid = 0; childid < num_children; childid++) {
    /* check whether each element has the same parent */
    t8_sele_parent (fam[childid], &compare);
    if (t8_sele_compare (&parent, &compare))
      return 0;
    /* check whether each element is the correct child of the collective parent */
    /* Could be replaced by type comparison as level is already checked in parent comparison */
    t8_sele_child (&parent, childid, &compare);
    if (t8_sele_compare (fam[childid], &compare))
      return 0;
  }
  return 1;
}

template <t8_eclass_t eclass_T>
int
t8_sele_is_inside_root (const t8_standalone_element_t<eclass_T> *p)
{
  t8_standalone_element_t<eclass_T> anc;
  t8_sele_ancestor_equation (p, 0, &anc);

  /**Check that we are in the correct cube*/
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    if (anc.coords[idim])
      return 0;
  }
  //  t8_debugf("anc.type %i\n", anc.type);
  return anc.type == 0;
}

template <t8_eclass_t eclass_T>
void
t8_sele_first_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                          const int level)
{
  T8_ASSERT (level >= p->level);
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*The first descendant of a pyramid has the same anchor coords and type, but another level */
  t8_sele_copy (p, desc);
  desc->level = level;
}

template <t8_eclass_t eclass_T>
void
t8_sele_last_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                         const int level)
{
  T8_ASSERT (level >= p->level);
  t8_sele_copy (p, desc);
  desc->level = level;
  /* Shift the coords to the eights cube. The type of the last descendant
   * is the type of the input element*/
  t8_element_coord_t coord_offset = t8_sele_get_len<eclass_T> (p->level) - t8_sele_get_len<eclass_T> (level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    desc->coords[i] |= coord_offset;
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_ancestor_id (const t8_standalone_element_t<eclass_T> *p, const int level)
{
  t8_standalone_element_t<eclass_T> ancestor;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_sele_ancestor_equation (p, level, &ancestor);
  return t8_sele_child_id (&ancestor);
}

template <t8_eclass_t eclass_T>
void
t8_sele_children (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> **c)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  const int num_children = t8_sele_num_children (p);
  for (int i = num_children - 1; i >= 0; i--) {
    t8_sele_child (p, i, c[i]);
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_parent (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *parent)
{
  T8_ASSERT (p->level > 0);

  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    const t8_cube_id_t cube_id = t8_sele_compute_cubeid (p, p->level);
    parent->type = t8_element_type_cubeid_to_parenttype<eclass_T>[p->type.to_ulong ()][cube_id];
  }

  const t8_element_coord_t length = t8_sele_get_len<eclass_T> (p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    parent->coords[i] = p->coords[i] & ~length;
  }

  parent->level = p->level - 1;
  T8_ASSERT (parent->level >= 0);
}

template <t8_eclass_t eclass_T>
void
t8_sele_successor (const t8_standalone_element_t<eclass_T> *elem, t8_standalone_element_t<eclass_T> *succ, const int level)
{
  T8_ASSERT (1 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_sele_copy (elem, succ);
  succ->level = level;

  const int child_id = t8_sele_child_id (elem);
  /*Compute number of children */
  const int num_siblings = t8_sele_num_siblings (elem);
  T8_ASSERT (0 <= child_id && child_id < num_siblings);
  if (child_id == num_siblings - 1) {
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_sele_successor (succ, succ, level - 1);
    succ->level = level;
    /* set bits from elem-->level to level to 0, emulates first descendant */
    int shift = T8_ELEMENT_MAXLEVEL[eclass_T] - level + 1;
    t8_sele_cut_coordinates (succ, shift);
    /* first descendants */
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_sele_parent (succ, succ);
    t8_sele_child (succ, child_id + 1, succ);
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_ancestor_loop (const t8_standalone_element_t<eclass_T> *el, const int level,
                       t8_standalone_element_t<eclass_T> *anc)
{
  T8_ASSERT (0 <= level && level <= el->level);
  t8_sele_copy (el, anc);
  for (int ilevel = el->level; ilevel < level; ilevel++) {
    t8_sele_parent (anc, anc);
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_ancestor_equation (const t8_standalone_element_t<eclass_T> *el, const int level,
                           t8_standalone_element_t<eclass_T> *anc)
{
  T8_ASSERT (0 <= level && level <= el->level);
  if (el != anc) {
    t8_sele_copy (el, anc);
  }
  if (el->level == level) {
    return;
  }

  /* Set type */
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    anc->type = t8_sele_compute_type_at_level (anc, level);
  }

  /* The coordinates and the type of the anc are defined by the level. */
  t8_sele_cut_coordinates (anc, T8_ELEMENT_MAXLEVEL[eclass_T] - level);

  anc->level = level;
}

template <t8_eclass_t eclass_T>
int
t8_sele_equal (const t8_standalone_element_t<eclass_T> *el1, const t8_standalone_element_t<eclass_T> *el2)
{
  if (el1->level != el2->level)
    return 0;
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    if (el1->coords[idim] != el2->coords[idim])
      return 0;
  }
  return el1->type == el2->type;
}

template <t8_eclass_t eclass_T>
static inline int
t8_sele_cube_ancestor_level (const t8_standalone_element_t<eclass_T> *el1, const t8_standalone_element_t<eclass_T> *el2)
{
  t8_element_coord_t maxexclor = 0;
  int level_inv;
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    maxexclor |= (el1->coords[idim] ^ el2->coords[idim]);
  }

  level_inv = SC_LOG2_32 (maxexclor) + 1;
  T8_ASSERT (level_inv <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  int min_value = SC_MIN (T8_ELEMENT_MAXLEVEL[eclass_T] - level_inv, (int) SC_MIN (el1->level, el2->level));
  return min_value;
}

template <t8_eclass_t eclass_T>
void
t8_sele_nearest_common_ancestor (const t8_standalone_element_t<eclass_T> *el1,
                                 const t8_standalone_element_t<eclass_T> *el2, t8_standalone_element_t<eclass_T> *nca)
{
  t8_standalone_element_t<eclass_T> nca2;
  int cube_anc_level = t8_sele_cube_ancestor_level<eclass_T> (el1, el2);
  int real_level = cube_anc_level;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    t8_element_type_t<eclass_T> anc_type, anc2_type;
    do {
      anc_type = t8_sele_compute_type_at_level (el1, real_level);
      anc2_type = t8_sele_compute_type_at_level (el2, real_level);
      real_level--;
    } while (anc_type != anc2_type);
    real_level++; /* we subtracted once too much*/
  }
  t8_sele_ancestor_equation (el1, real_level, nca);
}

/*******Linear id stuff *************/
template <t8_eclass_t eclass_T>
static t8_linearidx_t
t8_sele_num_descendants_at_leveldiff (const t8_standalone_element_t<eclass_T> *elem, const int leveldiff)
{
  if (leveldiff < 0)
    return 0;
  if constexpr (eclass_T == T8_ECLASS_PYRAMID) {
    t8_linearidx_t two_to_l = 1LL << leveldiff;
    t8_linearidx_t eight_to_l = 1LL << (3 * leveldiff);
    if (t8_sele_shape (elem) == T8_ECLASS_PYRAMID) {
      return ((eight_to_l << 2) - two_to_l) / 3;
    }
    else {
      return ((eight_to_l << 1) + two_to_l) / 3;
    }
  }
  return 1LL << (T8_ELEMENT_DIM[eclass_T] * leveldiff);
}

template <t8_eclass_t eclass_T>
static t8_linearidx_t
t8_sele_num_descendants_of_child_at_leveldiff (const t8_standalone_element_t<eclass_T> *p, int childindex,
                                               int leveldiff)
{
  t8_standalone_element_t<eclass_T> child;
  t8_sele_child (p, childindex, &child);
  t8_linearidx_t num_descendants = t8_sele_num_descendants_at_leveldiff (&child, leveldiff - 1);
  return num_descendants;
}

template <t8_eclass_t eclass_T>
static void
t8_sele_init_linear_id_recursive (t8_standalone_element_t<eclass_T> *p, const int level_diff, t8_linearidx_t id)
{
  T8_ASSERT (0 <= id);
  T8_ASSERT (1 <= level_diff && level_diff <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  if (id == 0) {
    t8_sele_first_descendant (p, p, p->level + level_diff);
    return;
  }

  if (level_diff == 1) {
    T8_ASSERT (id <= T8_ELEMENT_NUM_CHILDREN[eclass_T]);
    t8_sele_child (p, id, p);
    return;
  }

  t8_linearidx_t sum_descendants_of_children_before = 0;
  t8_linearidx_t num_descendants_of_child = 0;
  int childindex;
  /*If needed, can be replaced by binary search in lookuptable */
  for (childindex = 0; childindex < t8_sele_num_children (p); childindex++) {
    num_descendants_of_child = t8_sele_num_descendants_of_child_at_leveldiff (p, childindex, level_diff);
    sum_descendants_of_children_before += num_descendants_of_child;
    if (sum_descendants_of_children_before > id) {
      sum_descendants_of_children_before -= num_descendants_of_child;
      break;
    }
  }
  t8_sele_child (p, childindex, p);
  t8_sele_init_linear_id_recursive (p, level_diff - 1, id - sum_descendants_of_children_before);
}

template <t8_eclass_t eclass_T>
void
t8_sele_init_linear_id (t8_standalone_element_t<eclass_T> *p, const int level, t8_linearidx_t id)
{
  t8_sele_root (p);
  if (level == 0) {
    T8_ASSERT (id == 0);
    return;
  }
  t8_sele_init_linear_id_recursive (p, level, id);
}

template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_sele_linear_id_recursive (t8_standalone_element_t<eclass_T> *p, const t8_linearidx_t id, const int level_diff)
{
  if (p->level == 0)
    return id;

  const int childid = t8_sele_child_id (p);
  t8_sele_parent (p, p);
  t8_linearidx_t parent_id = 0;
  for (int ichild = 0; ichild < childid; ichild++) {
    /* p is now parent, so compute child to get sibling of original p */
    t8_linearidx_t num_child_descendants = t8_sele_num_descendants_of_child_at_leveldiff (p, ichild, level_diff + 1);
    parent_id += num_child_descendants;
  }
  parent_id += id;
  return t8_sele_linear_id_recursive (p, parent_id, level_diff + 1);
}

template <t8_eclass_t eclass_T>
t8_linearidx_t
t8_sele_linear_id (const t8_standalone_element_t<eclass_T> *p, const int level)
{
  t8_standalone_element_t<eclass_T> recursive_start;

  if (level < p->level) {
    t8_sele_ancestor_equation (p, level, &recursive_start);
  }
  else {
    t8_sele_first_descendant (p, &recursive_start, level);
  }

  /* Maybe we can also input p into recursive function and calculate id directly for first desc */
  t8_linearidx_t id = t8_sele_linear_id_recursive (&recursive_start, 0, 0);
  T8_ASSERT (id >= 0);
  return id;
}

/* Vertex information */
template <t8_eclass_t eclass_T>
void
t8_sele_compute_coords (const t8_standalone_element_t<eclass_T> *p, const int vertex, int coords[])
{

  T8_ASSERT (0 <= vertex && vertex < t8_sele_num_corners (p));

  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    int8_t type = p->type.to_ulong ();
    for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
      coords[idim]
        = p->coords[idim]
          + t8_type_vertex_dim_to_binary<eclass_T>[type][vertex][idim] * t8_sele_get_len<eclass_T> (p->level);
    }
  }
  else {
    //Hypercubes
    for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
      coords[idim] = p->coords[idim] + ((vertex & (1 << idim)) >> idim) * t8_sele_get_len<eclass_T> (p->level);
    }
  }
}

/***VISUALISATION****/

template <t8_eclass_t eclass_T>
void
t8_sele_vertex_reference_coords (const t8_standalone_element_t<eclass_T> *elem, const int vertex, double coords[])
{
  if constexpr(eclass_T == T8_ECLASS_VERTEX){
    return;
  }else{
    int coords_int[T8_ELEMENT_DIM[eclass_T]];
    T8_ASSERT (0 <= vertex && vertex < T8_ELEMENT_NUM_CORNERS[eclass_T]);
    t8_sele_compute_coords (elem, vertex, coords_int);
    /*scale the coordinates onto the reference cube */
    for(int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++){
      coords[i] = coords_int[i] / (double) t8_sele_get_root_len<eclass_T> ();
    }
  }
}

/**Face connection**/
template <t8_eclass_t eclass_T>
int
t8_sele_face_neighbor (const t8_standalone_element_t<eclass_T> *elem, t8_standalone_element_t<eclass_T> *neigh,
                       int face, int *neigh_face)
{
  t8_sele_copy (elem, neigh);

  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    int internal_face = t8_sele_lut_face_internal<eclass_T>[elem->type.to_ulong ()][face];
    if (internal_face) {
      /**Determine typebit*/
      int typebit = t8_sele_lut_type_face_to_typebit<eclass_T>[elem->type.to_ulong ()][face];
      /**Change typebit*/
      neigh->type.flip (typebit);
      *neigh_face = t8_sele_lut_type_face_to_neighface<eclass_T>[elem->type.to_ulong ()][face];
      return t8_sele_is_inside_root (neigh);
    }
  }

  /**Face is external, so a cube face*/
  /**Determine facenormal dimension*/
  int facenormal_dim, sign;
  /**Determine sign and adapt coordinate*/
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    sign = t8_sele_lut_type_face_to_sign<eclass_T>[elem->type.to_ulong ()][face];
    facenormal_dim = t8_sele_lut_type_face_to_facenormal_dim<eclass_T>[elem->type.to_ulong ()][face];
    T8_ASSERT (facenormal_dim != -1);
  }
  else {
    sign = face % 2 ? 1 : -1;
    facenormal_dim = face / 2;
  }

  /**Adapt coordinates*/
  t8_element_coord_t length = t8_sele_get_len<eclass_T> (elem->level);
  t8_debugf ("length: %i, sign:%i\n", length, sign);

  t8_debugf ("neigh_coords[%i]: %i\n", facenormal_dim, neigh->coords[facenormal_dim]);
  neigh->coords[facenormal_dim] += length * sign;
  t8_debugf ("neigh_coords[%i]: %i\n", facenormal_dim, neigh->coords[facenormal_dim]);

  /**Adapt typebits*/
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    for (size_t ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
      /**For all neighboring typebits, change typebit*/
      if (t8_type_edge_equations<eclass_T>[ieq][0] == facenormal_dim
          || t8_type_edge_equations<eclass_T>[ieq][1] == facenormal_dim) {
        neigh->type.flip (ieq); /*ASSERT that flip is in correct direction */
      }
    }
  }
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    *neigh_face = face ^ 1;
  }
  else {
    *neigh_face = t8_sele_lut_type_face_to_neighface<eclass_T>[elem->type.to_ulong ()][face];
  }
  /**check inside root*/
  return t8_sele_is_inside_root (neigh);
}

template <t8_eclass_t eclass_T>
int
t8_sele_face_parent_face (const t8_standalone_element_t<eclass_T> *elem, const int face)
{
  if (elem->level == 0)
    return -1;
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    if ((unsigned int) (face % 2) != ((elem->coords[face / 2]) >> (T8_ELEMENT_MAXLEVEL[eclass_T] - elem->level)) % 2) {
      return -1;
    }
    return face;
  }
  else {
    t8_cube_id_t cube_id = t8_sele_compute_cubeid<eclass_T> (elem, elem->level);
    return t8_sele_lut_type_cubeid_face_to_parentface<eclass_T>[elem->type.to_ulong ()][cube_id][face];
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_face_normal_dim (const t8_standalone_element_t<eclass_T> *elem, const int face)
{
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return face / 2;
  }
  else {
    return t8_sele_lut_type_face_to_facenormal_dim<eclass_T>[elem->type.to_ulong ()][face];
  }
}
template <t8_eclass_t eclass_T>
int
t8_sele_face_internal (const t8_standalone_element_t<eclass_T> *elem, const int face)
{
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return 0;
  }
  else {
    return t8_sele_lut_face_internal<eclass_T>[elem->type.to_ulong ()][face];
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_face_is_1_boundary (const t8_standalone_element_t<eclass_T> *elem, const int face)
{
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return face % 2;
  }
  else {
    return t8_sele_lut_type_face_to_is_1_boundary<eclass_T>[elem->type.to_ulong ()][face];
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_is_root_boundary (const t8_standalone_element_t<eclass_T> *elem, int face)
{
  if (!t8_sele_face_internal (elem, face)) {
    int dimid = t8_sele_face_normal_dim (elem, face);
    if (t8_sele_face_is_1_boundary (elem, face)) {
      // a_d must be full of 1s up to level l
      t8_element_coord_t coord_offset = t8_sele_get_root_len<eclass_T> () - t8_sele_get_len<eclass_T> (elem->level);
      if (elem->coords[dimid] != coord_offset) {
        return 0;
      }
      // all edges containing dimid must be fulfilled with x_d-a_d >= x_j-a_j or x_j-a_j <= x_d-a_d
      if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
        for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
          if (t8_type_edge_equations<eclass_T>[ieq][0] == dimid) {
            if (elem->type[ieq]) {
              return 0;
            }
          }
          else if (t8_type_edge_equations<eclass_T>[ieq][1] == dimid) {
            if (!elem->type[ieq]) {
              return 0;
            }
          }
        }
      }
    }
    else {
      //zeroboundary
      // x_d must be full of 0s up to level l
      if (elem->coords[dimid] != 0) {
        return 0;
      }
      // all edges containing dimid must be fulfilled with x_d-a_d <= x_j-a_j or x_j-a_j >= x_d-a_d
      if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
        for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
          if (t8_type_edge_equations<eclass_T>[ieq][0] == dimid) {
            if (!elem->type[ieq]) {
              return 0;
            }
          }
          else if (t8_type_edge_equations<eclass_T>[ieq][1] == dimid) {
            if (elem->type[ieq]) {
              return 0;
            }
          }
        }
      }
    }
  }
  else {
    // internalface
    // get graph edge e (or ieq) = (xi,xj)
    // ai = aj is necessary and sufficient
    if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      int ieq = t8_sele_lut_type_face_to_typebit<eclass_T>[elem->type.to_ulong ()][face];
      if (elem->coords[t8_type_edge_equations<eclass_T>[ieq][0]]
          != elem->coords[t8_type_edge_equations<eclass_T>[ieq][1]]) {
        return 0;
      }
    }
    else {
      SC_ABORT ("Cubes should not have internal faces!\n");
    }
  }
  return 1;
}

template <t8_eclass_t eclass_T>
void
t8_sele_first_descendant_face (const t8_standalone_element_t<eclass_T> *elem, int face,
                               t8_standalone_element_t<eclass_T> *first_desc, int level)
{
  //  t8_debugf("compute first_desc_face for type %i along face %i at level %i for elem:\n", elem->type.to_ulong(), face, level);
  t8_sele_debug_print<eclass_T> (elem);
  first_desc->level = level;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    first_desc->type = elem->type; /**TODO: Check if this is always true! */
  }
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    first_desc->coords[idim] = elem->coords[idim];
  }
  int face_is_1_boundary;
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    face_is_1_boundary = face % 2;
  }
  else {
    face_is_1_boundary = t8_sele_lut_type_face_to_is_1_boundary<eclass_T>[elem->type.to_ulong ()][face];
  }

  if (face_is_1_boundary) {  //the face is a xi=1 boundary
    int facenormal_dim;
    if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      facenormal_dim = face / 2;
    }
    else {
      facenormal_dim = t8_sele_lut_type_face_to_facenormal_dim<eclass_T>[elem->type.to_ulong ()][face];
      T8_ASSERT (facenormal_dim != -1);
    }
    //    t8_debugf("type: %i, face:%i, facenormal_dim: %i\n", elem->type.to_ulong(),face, facenormal_dim);
    t8_element_coord_t coord_offset = t8_sele_get_len<eclass_T> (elem->level) - t8_sele_get_len<eclass_T> (level);

    first_desc->coords[facenormal_dim] += coord_offset;
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_last_descendant_face (const t8_standalone_element_t<eclass_T> *elem, int face,
                              t8_standalone_element_t<eclass_T> *last_desc, int level)
{
  last_desc->level = level;
  t8_element_coord_t coord_offset = t8_sele_get_len<eclass_T> (elem->level) - t8_sele_get_len<eclass_T> (level);
  //  t8_debugf("t8_sele_last_descendant_face with offset %i\n", coord_offset);

  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    last_desc->type = elem->type; /**TODO: Check if this is always true! */
  }

  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    int multiplier = 1;
    if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      if (idim == face / 2) {
        multiplier = face % 2;
      }
    }
    else {
      t8_cube_id_t cube_id = t8_sele_lut_type_face_to_last_facechilds_cubeid<eclass_T>[elem->type.to_ulong ()][face];
      multiplier = (cube_id & (1 << idim)) >> idim;  // = cubeid[idim];
    }
    last_desc->coords[idim] = elem->coords[idim] + multiplier * coord_offset;
  }
  //  t8_debugf("Computed last descendant face:\n");
  //t8_sele_debug_print<eclass_T>(elem);
}

template <t8_eclass_t eclass_T>
void
t8_sele_children_at_face (const t8_standalone_element_t<eclass_T> *elem, int face,
                          t8_standalone_element_t<eclass_T> **children, int num_face_children, int *child_indices)
{
  int allocated_indices = 0;
  if (child_indices == NULL) {
    child_indices = T8_ALLOC_ZERO (int, num_face_children);
    allocated_indices = 1;
  }
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    int face_sign, face_dim;
    face_sign = face % 2;
    face_dim = face / 2;
    for (int ifacechild = 0; ifacechild < num_face_children; ifacechild++) {
      t8_element_coord_t first_part, face_part, last_part;
      /* ifacechild aaaabb, iface = x, then childid = aaaaxbb*/
      first_part = (ifacechild >> face_dim) << (face_dim + 1);
      last_part = ifacechild & ((1 << face_dim) - 1);
      face_part = face_sign << face_dim;
      child_indices[ifacechild] = first_part + face_part + last_part;
    }
  }
  else {
    for (int ifacechild = 0; ifacechild < num_face_children; ifacechild++) {
      child_indices[ifacechild]
        = t8_sele_lut_type_face_facechildid_to_childid<eclass_T>[elem->type.to_ulong ()][face][ifacechild];
      //      t8_debugf("looked up child index %i for type %i, face %i, ifacechild %i\n", child_indices[ifacechild], elem->type.to_ulong(), face, ifacechild);
    }
  }
  for (int ifacechild = 0; ifacechild < num_face_children; ifacechild++) {
    t8_sele_child (elem, child_indices[ifacechild], children[ifacechild]);
  }
  if (allocated_indices) {
    T8_FREE (child_indices);
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_face_child_face (const t8_standalone_element_t<eclass_T> *p, const int face, const int face_child)
{
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return face;
  }
  else {
    int child_id = t8_sele_lut_type_face_facechildid_to_childid<eclass_T>[p->type.to_ulong ()][face][face_child];
    return t8_sele_lut_type_childid_face_to_childface<eclass_T>[p->type.to_ulong ()][child_id][face];
  }
}

template <t8_eclass_t eclass_T>
int
t8_sele_tree_face (const t8_standalone_element_t<eclass_T> *elem, const int face)
{
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return face;
  }
  else {
    return t8_sele_lut_type_face_to_tree_face<eclass_T>[elem->type.to_ulong ()][face];
  }
}

template <t8_eclass_t eclass_T>
void
t8_sele_transform_face (const t8_standalone_element_t<eclass_T> *elem1, t8_standalone_element_t<eclass_T> *elem2,
                        int orientation, int sign, int is_smaller_face)
{
  if constexpr (eclass_T == T8_ECLASS_VERTEX) {
    return;
  }
  int level = elem1->level;
  if constexpr (eclass_T == T8_ECLASS_LINE) {
    t8_sele_copy<eclass_T> (elem1, elem2);
    if (orientation) {
      t8_element_coord_t total_length = 1 << T8_ELEMENT_MAXLEVEL[eclass_T];
      t8_element_coord_t refined_length = 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - level);
      elem2->coords[0] = total_length - refined_length - elem2->coords[0];
    }
    return;
  }
  if constexpr (eclass_T == T8_ECLASS_QUAD) {
    t8_standalone_element_t<eclass_T> tmp;
    if (sign) {
      /* The tree faces have the same topological orientation, and
      * thus we have to perform a coordinate switch. */
      /* We use p as storage, since elem1 and elem2 are allowed to
      * point to the same quad */
      t8_sele_copy<eclass_T> (elem1, &tmp);
      tmp.coords[0] = elem1->coords[1];
      tmp.coords[1] = elem1->coords[0];
    }
    else {
      t8_sele_copy<eclass_T> (elem1, &tmp);
    }

    /*
    * The faces of the root quadrant are enumerated like this:
    *
    *   v_2      v_3
    *     x -->-- x
    *     |       |
    *     ^       ^
    *     |       |
    *     x -->-- x
    *   v_0      v_1
    *
    * Orientation is the corner number of the bigger face that coincides
    * with the corner v_0 of the smaller face.
    */
    /* If this face is not smaller, switch the orientation:
    *  sign = 0   sign = 1
    *  0 -> 0     0 -> 0
    *  1 -> 2     1 -> 1
    *  2 -> 1     2 -> 2
    *  3 -> 3     3 -> 3
    */
    if (!is_smaller_face && (orientation == 1 || orientation == 2) && !sign) {
      orientation = 3 - orientation;
    }
    t8_element_coord_t root_len = t8_sele_get_root_len<eclass_T> ();
    t8_element_coord_t h = t8_sele_get_len<eclass_T> (elem1->level);
    switch (orientation) {
    case 0: /* Nothing to do */
      elem2->coords[0] = tmp.coords[0];
      elem2->coords[1] = tmp.coords[1];
      break;
    case 1:
      elem2->coords[0] = root_len - tmp.coords[1] - h;
      elem2->coords[1] = tmp.coords[0];
      break;
    case 2:
      elem2->coords[0] = tmp.coords[1];
      elem2->coords[1] = root_len - tmp.coords[0] - h;
      break;
    case 3:
      elem2->coords[0] = root_len - tmp.coords[0] - h;
      elem2->coords[1] = root_len - tmp.coords[1] - h;
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    elem2->level = tmp.level;
    return;
  }
  if constexpr (eclass_T == T8_ECLASS_TRIANGLE) {
    T8_ASSERT (0 <= orientation && orientation <= 2);
    t8_element_coord_t h = t8_sele_get_len<eclass_T> (elem1->level);

    // copy elem1 in tmp so that results can directly be set in elem2
    t8_standalone_element_t<eclass_T> tmp;
    t8_sele_copy<eclass_T> (elem1, &tmp);

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
      if (elem1->type == 0) {
        tmp.coords[1] = elem1->coords[0] - elem1->coords[1];
      }
      else {
        tmp.coords[1] = elem1->coords[0] - elem1->coords[1] - h;
      }
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

    elem2->level = tmp.level;
    elem2->type = tmp.type;
    switch (orientation) {
    case 0:
      t8_sele_copy<eclass_T> (&tmp, elem2);
      break;
    case 1:
      elem2->coords[0] = t8_sele_get_root_len<eclass_T> () - h - tmp.coords[1];
      if (tmp.type == 0) {
        elem2->coords[1] = tmp.coords[0] - tmp.coords[1];
      }
      else {
        elem2->coords[1] = tmp.coords[0] - tmp.coords[1] - h;
      }
      break;
    case 2:
      if (tmp.type == 0) {
        elem2->coords[0] = t8_sele_get_root_len<eclass_T> () - h + tmp.coords[1] - tmp.coords[0];
      }
      else {
        elem2->coords[0] = t8_sele_get_root_len<eclass_T> () + tmp.coords[1] - tmp.coords[0];
      }
      elem2->coords[1] = t8_sele_get_root_len<eclass_T> () - h - tmp.coords[0];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }
}

template <t8_eclass_t eclass_T, t8_eclass_t face_eclass_T>
void
t8_sele_boundary_face (const t8_standalone_element_t<eclass_T> *elem, const int root_face,
                       t8_standalone_element_t<face_eclass_T> *boundary)
{
  if constexpr (T8_ELEMENT_DIM[face_eclass_T] >= T8_ELEMENT_DIM[eclass_T]) {
    return;
  }

  boundary->level = elem->level;
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    t8_debugf ("consider dimension idim %i\n", idim);
    int ifacedim;
    if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      int facenormal_dim = root_face / 2;
      if (idim == facenormal_dim) {
        ifacedim = -1;
      }
      else if (idim > facenormal_dim) {
        ifacedim = idim - 1;
      }
      else {
        ifacedim = idim;
      }
    }
    else {
      ifacedim = t8_sele_lut_rootface_dim_to_facedim<eclass_T>[root_face][idim];
    }
    t8_debugf ("found ifacedim %i\n", ifacedim);
    if (ifacedim != -1) {
      if constexpr (face_eclass_T != T8_ECLASS_VERTEX && T8_ELEMENT_DIM[face_eclass_T] < T8_ELEMENT_DIM[eclass_T]) {
        t8_debugf ("set boundary coords %i\n", ifacedim);
        boundary->coords[ifacedim] = elem->coords[idim]
                                     << (T8_ELEMENT_MAXLEVEL[face_eclass_T] - T8_ELEMENT_MAXLEVEL[eclass_T]);
        t8_debugf ("to %i\n", boundary->coords[ifacedim]);
      }
      else {
        SC_ABORT_NOT_REACHED ();
      }
    }
  }
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
      t8_debugf ("consider equation ieq %i\n", ieq);
      int ifaceeq = t8_sele_lut_rootface_eq_to_faceeq<eclass_T>[root_face][ieq];
      t8_debugf ("found ifaceeq %i\n", ifaceeq);

      if (ifaceeq != -1) {
        t8_debugf ("boundarytype size: %li, elemtype size: %li\n", boundary->type.size (), elem->type.size ());
        boundary->type[ifaceeq] = elem->type[ieq];
        //        t8_debugf("set type[%i] = %i\n", ifaceeq, boundary->type[ifaceeq]);
      }
    }
  }
}

template <t8_eclass_t eclass_T, t8_eclass_t face_eclass_T>
int
t8_sele_extrude_face (const t8_standalone_element_t<face_eclass_T> *face, t8_standalone_element_t<eclass_T> *elem,
                      int root_face)
{
  /** Loop over elemdim, get corresponding facedim and set elem coord accordingly 
   * If elemdim is faceboundary, find out if 0 or 1 boundary
   */
  T8_ASSERT (0 <= face->level && face->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  if constexpr (T8_ELEMENT_DIM[face_eclass_T] >= T8_ELEMENT_DIM[eclass_T]) {
    return -1;
  }

  elem->level = face->level;
  t8_element_coord_t h = t8_sele_get_root_len<eclass_T> () - t8_sele_get_len<eclass_T> (elem->level);
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    int ifacedim;
    if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
      int facenormal_dim = root_face / 2;
      if (idim == facenormal_dim) {
        ifacedim = -1;
      }
      else if (idim > facenormal_dim) {
        ifacedim = idim - 1;
      }
      else {
        ifacedim = idim;
      }
    }
    else {
      ifacedim = t8_sele_lut_rootface_dim_to_facedim<eclass_T>[root_face][idim];
    }
    if (ifacedim != -1) {
      if constexpr (face_eclass_T != T8_ECLASS_VERTEX && T8_ELEMENT_DIM[face_eclass_T] < T8_ELEMENT_DIM[eclass_T]) {
        t8_debugf ("shift right by %i\n", (T8_ELEMENT_MAXLEVEL[face_eclass_T] - T8_ELEMENT_MAXLEVEL[eclass_T]));
        elem->coords[idim]
          = face->coords[ifacedim] >> (T8_ELEMENT_MAXLEVEL[face_eclass_T] - T8_ELEMENT_MAXLEVEL[eclass_T]);
      }
      else {
        SC_ABORT_NOT_REACHED ();
      }
    }
    else {
      int is_1_boundary;
      if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
        is_1_boundary = root_face % 2;
      }
      else {
        int root_type = 0;  //TODO: Make other root_types possible
        is_1_boundary = t8_sele_lut_type_face_to_is_1_boundary<eclass_T>[root_type][root_face];
      }
      if (is_1_boundary) {
        t8_debugf ("set to h\n");
        elem->coords[idim] = h;
      }
      else {
        elem->coords[idim] = 0;
      }
    }
    t8_debugf ("set elem->coords[%i] to %i\n", idim, elem->coords[idim]);
  }
  if constexpr (!T8_ELEMENT_NUM_EQUATIONS[eclass_T]) {
    return root_face;
  }
  else {
    t8_element_type_t<eclass_T> root_type = 0;
    elem->type = root_type;
    for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
      int ifaceeq = t8_sele_lut_rootface_eq_to_faceeq<eclass_T>[root_face][ieq];
      if (ifaceeq != -1) {
        elem->type[ieq] = face->type[ifaceeq];
      }
    }
    /** Set those typebits, that are connected to the face_normaldim of root_face*/
    for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; ieq++) {
      int facenormal_dim = t8_sele_lut_type_face_to_facenormal_dim<eclass_T>[root_type.to_ulong ()][root_face];
      if (t8_type_edge_equations<eclass_T>[ieq][0] == facenormal_dim) {
        elem->type[ieq] = !t8_sele_lut_type_face_to_is_1_boundary<eclass_T>[root_type.to_ulong ()][root_face];
      }
      else if (t8_type_edge_equations<eclass_T>[ieq][1] == facenormal_dim) {
        elem->type[ieq] = t8_sele_lut_type_face_to_is_1_boundary<eclass_T>[root_type.to_ulong ()][root_face];
      }
    }
    return t8_sele_lut_type_rootface_to_face<eclass_T>[elem->type.to_ulong ()][root_face];
  }
}