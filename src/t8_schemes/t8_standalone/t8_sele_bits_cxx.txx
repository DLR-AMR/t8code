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

/*********GETTER******/

/** The length of a element at a given level in integer coordinates */
template<t8_eclass_t eclass_T>
t8_element_coord_t  t8_sele_get_len(t8_element_level_t level)
{
  return 1 << (T8_ELEMENT_MAXLEVEL[eclass_T] - (level));
}

template<t8_eclass_t eclass_T>
t8_element_coord_t  t8_sele_get_root_len()
{
  return 1 << T8_ELEMENT_MAXLEVEL[eclass_T];
}

template<t8_eclass_t eclass_T>
t8_element_shape_t
t8_sele_shape (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  if constexpr(eclass_T == T8_ECLASS_PYRAMID){
    if (p->type == T8_DPYRAMID_FIRST_PYRA_TYPE ||
        p->type == T8_DPYRAMID_SECOND_PYRA_TYPE) {
      return T8_ECLASS_PYRAMID;
    }
    else {
      return T8_ECLASS_TET;
    }
  }
  return eclass_T;
}

template<t8_eclass_t eclass_T>
int
t8_sele_num_corners (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*TODO*/
  if constexpr(eclass_T==T8_ECLASS_PYRAMID){
    if(t8_sele_shape(p)==T8_ECLASS_PYRAMID){
      return 5;
    }
    return 4;
  }
  return T8_ELEMENT_NUM_CORNERS[eclass_T];
}

template<t8_eclass_t eclass_T>
int
t8_sele_num_children (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*TODO*/
  if constexpr(eclass_T==T8_ECLASS_PYRAMID){
    if(t8_sele_shape(p)==T8_ECLASS_PYRAMID){
      return 10;
    }
    return 6;
  }
  return T8_ELEMENT_NUM_CHILDREN[eclass_T];
}

template<t8_eclass_t eclass_T>
int
t8_sele_num_siblings (const t8_standalone_element_t<eclass_T> *p)
{
  if (p->level == 0)
    return 1;
  T8_ASSERT (0 < p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_standalone_element_t<eclass_T>       parent;
  t8_sele_parent (p, &parent);
  return t8_sele_num_children (&parent);
}

template<t8_eclass_t eclass_T>
int
t8_sele_num_faces (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  return T8_ELEMENT_NUM_FACES[eclass_T];
}

template<t8_eclass_t eclass_T>
int
t8_sele_max_num_faces (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  return T8_ELEMENT_NUM_FACES[eclass_T];
}

/** Cube helper **/

template<t8_eclass_t eclass_T>
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
template<t8_eclass_t eclass_T>
t8_element_type_t<eclass_T>
t8_sele_compute_type_at_level (const t8_standalone_element_t<eclass_T> *p, int level)
{
  t8_element_type_t<eclass_T>              type = 0;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_element_coord_t coord_v0 = 0;// p->coords[t8_type_edge_equations[e][0]];
    t8_element_coord_t coord_v1 = 0;// p->coords[t8_type_edge_equations[e][1]];

    coord_v0 = (coord_v0 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);
    coord_v1 = (coord_v1 << level) & ((1 << T8_ELEMENT_MAXLEVEL[eclass_T]) - 1);

    if (coord_v0 == coord_v1) {
      type |= (p->type & 1 << e);
    }
    else if (coord_v0 < coord_v1) {
      type |= (1 << e);
    }
    else {
      T8_ASSERT (coord_v0 > coord_v1);
      T8_ASSERT (type && (1 << e) == 0);
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
template<t8_eclass_t eclass_T>
static void
t8_sele_cut_coordinates (t8_standalone_element_t<eclass_T> *p, const int shift)
{
  T8_ASSERT (0 <= shift && shift <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    p->coords[i] = (p->coords[i] >> shift) << shift;
  }
}

/*Copies a element from source to dest*/
template<t8_eclass_t eclass_T>
void
t8_sele_copy (const t8_standalone_element_t<eclass_T> *source, t8_standalone_element_t<eclass_T> *dest)
{
  T8_ASSERT (source != NULL && dest != NULL);
  if (source == dest) {
    return;
  }
  memcpy (dest, source, sizeof (t8_standalone_element_t<eclass_T>));
}

template<t8_eclass_t eclass_T>
int
t8_sele_compare (const t8_standalone_element_t<eclass_T> *p1, const t8_standalone_element_t<eclass_T> *p2)
{
  const int           maxlvl = SC_MAX (p1->level, p2->level);

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

template<t8_eclass_t eclass_T>
int
t8_sele_is_valid (const t8_standalone_element_t<eclass_T> *p)
{
  int                 is_valid;
  const t8_element_coord_t max_coord = ((uint64_t)2 * (uint64_t)t8_sele_get_root_len<eclass_T>) - 1;

  /*Check the level */
  is_valid = 0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T];

  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    is_valid = is_valid && - (int64_t)t8_sele_get_root_len<eclass_T> <= p->coords[i]
      && p->coords[i] <= max_coord;
  }

  return is_valid;
}

template<t8_eclass_t eclass_T>
void
t8_sele_debug_print (const t8_standalone_element_t<eclass_T> *p)
{
  t8_debugf("level: %i\n", p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    t8_debugf("x_%i: %i \n", i, p->coords[i]);
  }
  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_debugf("t_%i: %i \n", e, p->type[e]);
  }
}

template<t8_eclass_t eclass_T>
void
t8_sele_global_print (const t8_standalone_element_t<eclass_T> *p)
{
  t8_global_productionf("level: %i", p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    t8_global_productionf("x_%i: %i \n", i, p->coords[i]);
  }
  for (int e = 0; e < T8_ELEMENT_NUM_EQUATIONS[eclass_T]; e++) {
    t8_global_productionf("t_%i: %i \n", e, p->type[e]);
  }
}

template<t8_eclass_t eclass_T>
static void
t8_sele_root (t8_standalone_element_t<eclass_T> *p)
{
  p->level = 0;
  for (size_t i=0; i<T8_ELEMENT_DIM[eclass_T]; i++) {
    p->coords[i] = 0;
  }
  p->type = 0;
}

/**SFC functionality*/
template<t8_eclass_t eclass_T>
int
t8_sele_child_id (const t8_standalone_element_t<eclass_T> *p)
{
  T8_ASSERT (p->level >= 0);
  if (p->level == 0) {
    return -1;
  }
  const t8_cube_id_t cube_id = t8_sele_compute_cubeid (p, p->level);
  int8_t child_id;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[eclass_T]){
    child_id = t8_element_type_cubeid_to_Iloc<eclass_T>[p->type.to_ulong()][cube_id];
  }else{
    child_id = cube_id;
  }
  return child_id;
}

template<t8_eclass_t eclass_T>
void
t8_sele_child (const t8_standalone_element_t<eclass_T> *elem, const int child_id,
                   t8_standalone_element_t<eclass_T> *child)
{
  T8_ASSERT (0 <= child_id && child_id < T8_ELEMENT_NUM_CHILDREN[eclass_T]);
  T8_ASSERT (0 <= elem->level && elem->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);

  /* Compute the cube id and shift the coordinates accordingly */
  t8_cube_id_t cube_id;
  if constexpr(T8_ELEMENT_NUM_EQUATIONS[eclass_T]){
    cube_id = t8_element_type_Iloc_to_childcubeid<eclass_T>[elem->type.to_ulong()][child_id];
    child->type = t8_element_type_Iloc_to_childtype<eclass_T>[elem->type.to_ulong()][child_id];
  } else {
    cube_id = child_id;
  }

  const t8_element_coord_t length = t8_sele_get_len<eclass_T> (elem->level + 1);

  // Auslagern als bithelper funktion?
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    child->coords[i] = elem->coords[i] + ((cube_id & 1 << i) ? length : 0);
  }

  child->level = elem->level + 1;
}

template<t8_eclass_t eclass_T>
int
t8_sele_is_family (t8_standalone_element_t<eclass_T> **fam)
{
  t8_standalone_element_t<eclass_T>       parent, compare;
  /* Take the parent of the first element as baseline to compare against */
  t8_sele_parent (fam[0], &parent);
  const int           num_children = t8_sele_num_children (&parent);
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

template<t8_eclass_t eclass_T>
int
t8_sele_is_inside_root (const t8_standalone_element_t<eclass_T> *p)
{
  return 1;
}

template<t8_eclass_t eclass_T>
void
t8_sele_first_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                              const int level)
{
  T8_ASSERT (level >= p->level);
  T8_ASSERT (0 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  /*The first descendant of a pyramid has the same anchor coords and type, but another level */
  t8_sele_copy (p, desc);
  desc->level = level;
  //t8_debugf("first_desc\n");
  //t8_sele_debug_print(desc);
}

template<t8_eclass_t eclass_T>
void
t8_sele_last_descendant (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *desc,
                             const int level)
{
  T8_ASSERT (level >= p->level);
  t8_sele_copy (p, desc);
  desc->level = level;
  /* Shift the coords to the eights cube. The type of the last descendant is
   * is the type of the input element*/
  t8_element_coord_t coord_offset =
    t8_sele_get_len<eclass_T> (p->level) - t8_sele_get_len<eclass_T> (level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    desc->coords[i] |= coord_offset;
  }
}

template<t8_eclass_t eclass_T>
int
t8_sele_ancestor_id (const t8_standalone_element_t<eclass_T> *p, const int level)
{
  t8_standalone_element_t<eclass_T>       ancestor;
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_sele_ancestor (p, level, &ancestor);
  return t8_sele_child_id (&ancestor);
}

template<t8_eclass_t eclass_T>
void
t8_sele_children (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> **c)
{
  T8_ASSERT (0 <= p->level && p->level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  const int           num_children = t8_sele_num_children (p);
  for (int i = num_children - 1; i >= 0; i--) {
    t8_sele_child (p, i, c[i]);
  }
}

template<t8_eclass_t eclass_T>
void
t8_sele_parent (const t8_standalone_element_t<eclass_T> *p, t8_standalone_element_t<eclass_T> *parent)
{
  T8_ASSERT (p->level > 0);
  const t8_element_coord_t length = t8_sele_get_len<eclass_T> (p->level);
  const t8_cube_id_t cube_id =
    t8_sele_compute_cubeid (p, p->level);
  for (int i = 0; i < T8_ELEMENT_DIM[eclass_T]; i++) {
    parent->coords[i] = p->coords[i] & ~length;
  }
  parent->level = p->level - 1;
  T8_ASSERT (parent->level >= 0);
}

template<t8_eclass_t eclass_T>
void
t8_sele_successor (const t8_standalone_element_t<eclass_T> *elem,
                       t8_standalone_element_t<eclass_T> *succ, const int level)
{
  T8_ASSERT (1 <= level && level <= T8_ELEMENT_MAXLEVEL[eclass_T]);
  t8_sele_copy (elem, succ);
  succ->level = level;

  const int           child_id = t8_sele_child_id (elem);
  /*Compute number of children */
  const int           num_siblings = t8_sele_num_siblings (elem);
  T8_ASSERT (0 <= child_id && child_id < num_siblings);
  if (child_id == num_siblings - 1) {
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_sele_successor (succ, succ, level - 1);
    succ->level = level;
    /* set bits from elem-->level to level to 0, emulates first descendant */
    int                 shift = T8_ELEMENT_MAXLEVEL[eclass_T] - level + 1;
    t8_sele_cut_coordinates (succ, shift);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_sele_parent (succ, succ);
    t8_sele_child (succ, child_id + 1, succ);
  }
}

template<t8_eclass_t eclass_T>
void
t8_sele_ancestor (const t8_standalone_element_t<eclass_T> *el, const int level,
                      t8_standalone_element_t<eclass_T> *anc)
{
  T8_ASSERT (0 <= level && level <= el->level);
  t8_sele_copy (el, anc);
  if (el->level == level) {
    return;
  }
  else if (level == el->level - 1) {
    /* We can reuse the parent computation if we want to go only one level up. */
    t8_sele_parent (el, anc);
    return;
  }
  /* The coordinates and the type of the anc are defined by the level. */
  t8_sele_cut_coordinates (anc, T8_ELEMENT_MAXLEVEL[eclass_T] - level);

  anc->level = level;
}

template<t8_eclass_t eclass_T>
static int
t8_sele_nca_level (const t8_standalone_element_t<eclass_T> *el1, const t8_standalone_element_t<eclass_T> *el2)
{
  SC_ABORT ("Not implemented.");
  return 0;
}

template<t8_eclass_t eclass_T>
void
t8_sele_nearest_common_ancestor (const t8_standalone_element_t<eclass_T> *el1,
                                     const t8_standalone_element_t<eclass_T> *el2,
                                     t8_standalone_element_t<eclass_T> *nca)
{
  const t8_element_level_t nca_level = t8_sele_nca_level (el1, el2);
  t8_sele_ancestor (el1, nca_level, nca);
}

/*******Linear id stuff *************/
template<t8_eclass_t eclass_T>
static              t8_linearidx_t
t8_sele_num_descendants_at_leveldiff (const t8_standalone_element_t<eclass_T> *p,
                                          const int leveldiff)
{
  return 1 << leveldiff * T8_ELEMENT_DIM[eclass_T];
//  t8_linearidx_t      two_to_l = 1LL << leveldiff;
//  t8_linearidx_t      eight_to_l = 1LL << (3 * leveldiff);
////  //t8_debugf("leveldiff: %d, 2^l: %lu, 8^l: %lu\n", leveldiff, two_to_l, eight_to_l);
//  if (t8_sele_shape (p) == T8_ECLASS_PYRAMID) {
//    return ((eight_to_l << 2) - two_to_l) / 3;
//  }
//  else {
//    return ((eight_to_l << 1) + two_to_l) / 3;
//  }
}

template<t8_eclass_t eclass_T>
static              t8_linearidx_t
t8_sele_num_descendants_of_child_at_leveldiff (const t8_standalone_element_t<eclass_T> *p,
                                                   int childindex,
                                                   int leveldiff)
{
  t8_standalone_element_t<eclass_T>       child;
  t8_sele_child (p, childindex, &child);
  t8_linearidx_t      num_descendants = t8_sele_num_descendants_at_leveldiff (&child, leveldiff - 1);       //????
//  //t8_debugf("%d descendants of child %d with leveldiff %d\n", num_descendants, childindex, leveldiff);
  return num_descendants;
}

template<t8_eclass_t eclass_T>
static void
t8_sele_init_linear_id_recursive (t8_standalone_element_t<eclass_T> *p, const int level_diff,
                                      t8_linearidx_t id)
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

  t8_linearidx_t      sum_descendants_of_children_before = 0;
  t8_linearidx_t      num_descendants_of_child = 0;
  int                 childindex;
  /*If needed, can be replaced by binary search in lookuptable */
  for (childindex = 0; childindex < t8_sele_num_children (p);
       childindex++) {
    num_descendants_of_child =
      t8_sele_num_descendants_of_child_at_leveldiff (p, childindex,
                                                         level_diff);
    sum_descendants_of_children_before += num_descendants_of_child;
    if (sum_descendants_of_children_before > id) {
      sum_descendants_of_children_before -= num_descendants_of_child;
      break;
    }
  }
  t8_sele_child (p, childindex, p);
  t8_sele_init_linear_id_recursive (p, level_diff - 1,
                                        id -
                                        sum_descendants_of_children_before);
}

template<t8_eclass_t eclass_T>
void
t8_sele_init_linear_id (t8_standalone_element_t<eclass_T> *p, const int level,
                            t8_linearidx_t id)
{
  t8_sele_root (p);
  if (level == 0) {
    T8_ASSERT (id == 0);
    return;
  }
  t8_sele_init_linear_id_recursive (p, level, id);
}

template<t8_eclass_t eclass_T>
t8_linearidx_t
t8_sele_linear_id_recursive (t8_standalone_element_t<eclass_T> *p, const t8_linearidx_t id,
                                 const int level_diff)
{
  if (p->level == 0)
    return id;

  const int           childid = t8_sele_child_id (p);
  t8_sele_parent (p, p);
  t8_linearidx_t      parent_id = 0;
  for (int ichild = 0; ichild < childid; ichild++) {
    /* p is now parent, so compute child to get sibling of original p */
    t8_linearidx_t      num_child_descendants =
      t8_sele_num_descendants_of_child_at_leveldiff (p, ichild,
                                                         level_diff + 1);
    parent_id += num_child_descendants;
  }
  parent_id += id;
  return t8_sele_linear_id_recursive (p, parent_id, level_diff + 1);
}

template<t8_eclass_t eclass_T>
t8_linearidx_t
t8_sele_linear_id (const t8_standalone_element_t<eclass_T> *p, const int level)
{
  //t8_debugf("linear_id at level %d\n", level);
  //t8_sele_debug_print(p);
  t8_standalone_element_t<eclass_T>       recursive_start;

  if (level < p->level) {
    //t8_debugf("ancestor\n");
    t8_sele_ancestor (p, level, &recursive_start);
  }
  else {
    //t8_debugf("first_descendant\n");
    t8_sele_first_descendant (p, &recursive_start, level);
  }
  //t8_sele_debug_print(p);

  /* Maybe we can also input p into recursive function and calculate id directly for first desc */
  t8_linearidx_t      id =
    t8_sele_linear_id_recursive (&recursive_start, 0, 0);
  T8_ASSERT (id >= 0);
  return id;
}

/* Vertex information */
template<t8_eclass_t eclass_T>
void
t8_sele_compute_coords (const t8_standalone_element_t<eclass_T> *p, const int vertex,
                            int coords[])
{

  T8_ASSERT (0 <= vertex && vertex < t8_sele_num_corners (p));
  for (int idim = 0; idim < T8_ELEMENT_DIM[eclass_T]; idim++) {
    // Only for hypercubes
    coords[idim] = p->coords[idim] + ((vertex & (1 << idim)) >> idim) * t8_sele_get_len<eclass_T> (p->level);
    //t8_standalone_element_type_t  type = t8_sele_get_type (p);
    //coords[idim] =
    //  p->coords[idim] +
    //  t8_standalone_element_type_vertex_dim_to_binary[type][vertex][idim] *
    //  t8_sele_get_len<eclass_T> (p->level);
  }
}

/***VISUALISATION****/

template<t8_eclass_t eclass_T>
void
t8_sele_vertex_reference_coords (const t8_standalone_element_t<eclass_T> *elem,
                                     const int vertex, double coords[])
{
  int                 coords_int[T8_ELEMENT_DIM[eclass_T]];
  T8_ASSERT (0 <= vertex && vertex < T8_ELEMENT_NUM_CORNERS[eclass_T]);
  t8_sele_compute_coords (elem, vertex, coords_int);
  /*scale the coordinates onto the reference cube */
  coords[0] = coords_int[0] / (double) t8_sele_get_root_len<eclass_T> ();
  coords[1] = coords_int[1] / (double) t8_sele_get_root_len<eclass_T> ();
  coords[2] = coords_int[2] / (double) t8_sele_get_root_len<eclass_T> ();
}
