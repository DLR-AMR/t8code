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

  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  const t8_dpyramid_coord_t h = T8_DPYRAMID_LEN (level);

  if (level == 0) {
    return 0;
  }
  cube_id |= ((p.x & h) ? 0x01 : 0);
  cube_id |= ((p.y & h) ? 0x02 : 0);
  cube_id |= ((p.z & h) ? 0x04 : 0);

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
compute_type_at_level (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_type_t  type;

  return type;
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
  p.x = (p.x >> shift) << shift;
  p.y = (p.y >> shift) << shift;
  p.z = (p.z >> shift) << shift;
}


int
t8_dpyramid_ancestor_id (const t8_dpyramid_t *p, const int level)
{
  t8_dpyramid_t       helper;
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_ancestor (p, level, &helper);
  return t8_dpyramid_child_id (&helper);
}

int
t8_dpyramid_is_family (t8_dpyramid_t **fam)
{
  return false;
}

int
t8_dpyramid_is_root_boundary (const t8_dpyramid_t *p, const int face)
{
  return false;
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
  T8_ASSERT (p1.x >= 0 && p1.y >= 0 && p1.z >= 0 &&
             p1.level >= 0
             && p1.level <= T8_DPYRAMID_MAXLEVEL);
  T8_ASSERT (p2.x >= 0 && p2.y >= 0 && p2.z >= 0
             && p2.level >= 0
             && p2.level <= T8_DPYRAMID_MAXLEVEL);
  const int           maxlvl = SC_MAX (p1.level, p2.level);

  const t8_linearidx_t id1 = t8_dpyramid_linear_id (p1, maxlvl);
  const t8_linearidx_t id2 = t8_dpyramid_linear_id (p2, maxlvl);
  if (id1 == id2) {
    /* The linear ids are the same, the pyramid with the smaller level
     * is considered smaller */
    if (p1.level == p2.level) {
      T8_ASSERT (p1.type == p2.type);
      return 0;
    }
    else {
      return p1.level - p2.level;
    }
  }
  /* return negative if id1 < id2, zero if id1 = id2, positive if id1 >
     id2 */
  return id1 < id2 ? -1 : 1;
}

int
t8_dpyramid_get_level (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  return p.level;
}


void
t8_dpyramid_init_linear_id (t8_dpyramid_t *p, const int level,
                            t8_linearidx_t id)
{
}

t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t *p, const int level)
{
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
}

int
t8_dpyramid_face_parent_face (const t8_dpyramid_t *elem, const int face)
{
}

int
t8_dpyramid_tree_face (const t8_dpyramid_t *p, const int face)
{
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_is_root_boundary (p, face)) {
    if (t8_dpyramid_shape (p) == T8_ECLASS_PYRAMID) {
      /*If p is a pyramid and touches the boundary, the face-number is the same */
      return face;
    }
    else {
      /*p is a tet and in some occasions p shares a face with its tree */
      if (face == 0 && (p.type == 3 || p.type == 2)) {
        return 3;
      }
      else if (face == 0 && (p.type == 0 || p.type == 1)) {
        return 1;
      }
      else if ((face == 1 && p.type == 3)
               || (face == 2 && p.type == 1)) {
        return 2;
      }
      else if ((face == 1 && p.type == 0)
               || (face == 2 && p.type == 2)) {
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
}

void
t8_dpyramid_first_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc,
                              const int level)
{
  T8_ASSERT (level >= p.level);
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  /*The first descendant of a pyramid has the same anchor coords and type, but another level */
  t8_dpyramid_copy (p, desc);
  desc.level = level;
  T8_ASSERT (p.x <= desc.x
             && p.y <= desc.y
             && p.z <= desc.z);
}

void
t8_dpyramid_last_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc,
                             const int level)
{
  T8_ASSERT (level >= p.level);
  t8_dpyramid_copy (p, desc);
  desc.level = level;
  /* Shift the coords to the eights cube. The type of the last descendant is
  * is the type of the input pyramid*/
  t8_dpyramid_coord_t coord_offset =
    T8_DPYRAMID_LEN (p.level) - T8_DPYRAMID_LEN (level);
  desc.x |= coord_offset;
  desc.y |= coord_offset;
  desc.z |= coord_offset;
}

int
t8_dpyramid_num_corners (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (t8_dpyramid_parent(p)) == T8_ECLASS_PYRAMID) {
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
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
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
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  if (t8_dpyramid_shape (p) == T8_ECLASS_TET) {
    return T8_DTET_FACES;
  }
  else {
    return T8_DPYRAMID_FACES;
  }
}

int
t8_dpyramid_child_id (const t8_dpyramid_t *p)
{
  T8_ASSERT (p.level >= 0);
  if (p.level == 0) {
    return -1;
  }
  const t8_dpyramid_cube_id_t cube_id =
    compute_cubeid (p, p.level);
  T8_ASSERT (t8_dpyramid_type_cid_to_Iloc[p.type][cube_id] >= 0);
  return t8_dpyramid_type_cid_to_Iloc[p.type][cube_id];
}

void
t8_dpyramid_child (const t8_dpyramid_t *elem, const int child_id,
                   t8_dpyramid_t *child)
{
  T8_ASSERT (0 <= child_id && child_id < T8_DPYRAMID_CHILDREN);
  T8_ASSERT (0 <= elem.level
             && elem.level <= T8_DPYRAMID_MAXLEVEL);
  /* Compute the cube id and shift the coordinates accordingly */
  const t8_dpyramid_cube_id_t cube_id =
    t8_dpyramid_parenttype_Iloc_to_cid[elem.type][child_id];
  const t8_dpyramid_coord_t length =
    T8_DPYRAMID_LEN (elem.level + 1);
  T8_ASSERT (cube_id >= 0);
  child.level = elem.level + 1;
  child.x = elem.x + ((cube_id & 0x01) ? length : 0);
  child.y = elem.y + ((cube_id & 0x02) ? length : 0);
  child.z = elem.z + ((cube_id & 0x04) ? length : 0);
  child.type =
    t8_dpyramid_parenttype_Iloc_to_type[elem.type][child_id];
  T8_ASSERT (child.type >= 0);
}

void
t8_dpyramid_children (const t8_dpyramid_t *p, t8_dpyramid_t **c)
{
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  const int           num_children = t8_dpyramid_num_children (p);
  for (int i = num_children - 1; i >= 0; i--) {
    t8_dpyramid_child (p, i, c[i]);
  }
}

void
t8_dpyramid_parent (const t8_dpyramid_t *p, t8_dpyramid_t *parent)
{
  T8_ASSERT (p.level > 0);
  T8_ASSERT (T8_DPYRAMID_MAXLEVEL == T8_DTET_MAXLEVEL);
  const t8_dpyramid_coord_t length = T8_DPYRAMID_LEN (p.level);

  const t8_dpyramid_cube_id_t cube_id =
    compute_cubeid (p, p.level);

  parent.type =
    t8_dpyramid_type_cid_to_parenttype[p.type -
                                        T8_DPYRAMID_FIRST_TYPE][cube_id];
  parent.x = p.x & ~length;
  parent.y = p.y & ~length;
  parent.z = p.z & ~length;
  T8_ASSERT (parent.type >= 0);
  parent.level = p.level - 1;
  T8_ASSERT (parent.level >= 0);
}

t8_element_shape_t
t8_dpyramid_shape (const t8_dpyramid_t *p)
{
  T8_ASSERT (0 <= p.level
             && p.level <= T8_DPYRAMID_MAXLEVEL);
  /*The pyramid has the shape of a tetrahedron */
  if (p.type < T8_DPYRAMID_FIRST_TYPE) {
    return T8_ECLASS_TET;
  }
  else {
    return T8_ECLASS_PYRAMID;
  }

}

static void
t8_dpyramid_successor_recursion (const t8_dpyramid_t *elem,
                                 t8_dpyramid_t *succ, const int level)
{
  T8_ASSERT (1 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  t8_dpyramid_copy (elem, succ);
  succ.level = level;
  T8_ASSERT (succ.type >= 0);
  const int           child_id = t8_dpyramid_child_id (elem);
  /*Compute number of children */
  const int           num_siblings = t8_dpyramid_num_siblings (elem);
  T8_ASSERT (0 <= child_id && child_id < num_siblings);
  if (child_id == num_siblings - 1) {
    int                 shift = T8_DPYRAMID_MAXLEVEL - level + 1;
    /* Last-child-case. The successor is the successor of the parent element,
     * but with the given level */
    t8_dpyramid_successor_recursion (succ, succ, level - 1);
    succ.level = level;
    /* bits auf level auf child 0 setzen */
    t8_dpyramid_cut_coordinates (succ, shift);
  }
  else {
    /* Not the last element. Compute child with local ID child_id+1 */
    t8_dpyramid_parent (succ, succ);
    t8_dpyramid_child (succ, child_id + 1, succ);
    succ.level = level;
  }
}

void
t8_dpyramid_successor (const t8_dpyramid_t *elem, t8_dpyramid_t *succ,
                       const int level)
{
  t8_dpyramid_successor_recursion (elem, succ, level);
#ifdef T8_ENABLE_DEBUG
  if (t8_dpyramid_shape (succ) == T8_ECLASS_PYRAMID) {
    T8_ASSERT (succ->switch_shape_at_level < 0);
  }
  else {
    T8_ASSERT (succ->switch_shape_at_level =
               t8_dpyramid_compute_switch_shape_at_level (succ));
  }
#endif
}

void
t8_dpyramid_compute_coords (const t8_dpyramid_t *p, const int vertex,
                            int coords[])
{
  T8_ASSERT (0 <= vertex && vertex < T8_DPYRAMID_CORNERS);
  SC_ABORT("not implemented");
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
  T8_ASSERT (0 <= level && level <= pyra.level);
  t8_dpyramid_copy (pyra, anc);
  if (pyra.level == level) {
    return;
  }
  else if (level == pyra.level - 1) {
    /* We can reuse the parent computation if we want to go only one level up. */
    t8_dpyramid_parent (pyra, anc);
    return;
  }
  /* The coordinates of the anc are defined by the level. */
  t8_dpyramid_cut_coordinates (anc, T8_DPYRAMID_MAXLEVEL - level);
  anc.level = level;
  anc.type = t8_dpyramid_type_at_level (pyra, level);
}

void
t8_dpyramid_nearest_common_ancestor (const t8_dpyramid_t *pyra1,
                                     const t8_dpyramid_t *pyra2,
                                     t8_dpyramid_t *nca)
{
  SC_ABORT("not implemented!\n");
}

int
t8_dpyramid_is_valid (const t8_dpyramid_t *p)
{
  int                 is_valid;
  const t8_dpyramid_coord_t max_coord =
    ((int64_t) 2 * T8_DPYRAMID_ROOT_LEN) - 1;
  const t8_element_shape_t shape = t8_dpyramid_shape (p);
  /*Check the level */
  is_valid = 0 <= p.level
    && p.level <= T8_DPYRAMID_MAXLEVEL;
  /*Check coordinates, we allow a boundary layer around the root-pyramid */
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p.x
    && p.x <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p.y
    && p.y <= max_coord;
  is_valid = is_valid && -T8_DPYRAMID_ROOT_LEN <= p.z
    && p.z <= max_coord;

  /*The shape can be a pyramid or a tet */
  is_valid = is_valid && (shape == T8_ECLASS_PYRAMID
                          || shape == T8_ECLASS_TET);
  /*Check the type */
  is_valid = is_valid && 0 <= p.type
    && p.type < T8_DPYRAMID_NUM_TYPES;

  return is_valid;
}

void
t8_dpyramid_debug_print (const t8_dpyramid_t *p)
{
  t8_debugf
    ("x: %i, y: %i, z: %i, type %i, level: %i\n",
     p.x, p.y, p.z, p.type, p.level);
}
