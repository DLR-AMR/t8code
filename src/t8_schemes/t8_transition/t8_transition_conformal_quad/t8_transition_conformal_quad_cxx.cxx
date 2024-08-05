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

/* Description:
 * This is the low-level structure of 2D quadrilateral elements with transition cells of triangular subelements. */

#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include "t8.h"
#include "t8_transition_conformal_quad_cxx.hxx"

/* *INDENT-OFF* */
/* Connectivity of subelement faces: 
 *     f_0 <-> f_2 
 *     f_1 <-> f_1 (assuming a neighboring transition cell)
 *     f_2 <-> f_0 */
const int subelement_face_dual[3] = { 2, 1, 0 };

/* Connectivity of a subelements location within a transition cell 
 * and the parent quads faces:
 *     location[0] = 0 -> parents face = 1
 *     location[0] = 1 -> parents face = 2
 *     location[0] = 2 -> parents face = 0
 *     location[0] = 3 -> parents face = 3 */
const int subelement_location_to_parent_dual_face[4] = { 1, 2, 0, 3 };

/* Connectivity of a subelements location within a transition cell 
 * and the parent quads faces:
 *     location[0] = 0 (clockwise) -> parents face = 0
 *     location[0] = 1 (clockwise) -> parents face = 3
 *     location[0] = 2 (clockwise) -> parents face = 1
 *     location[0] = 3 (clockwise) -> parents face = 2 */
const int subelement_location_to_parent_face[4] = { 0, 3, 1, 2 };
/* *INDENT-ON* */

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This function is used by other element functions and we thus need to
 * declare it up here */
t8_linearidx_t
t8_element_get_linear_id (const t8_element_t *elem, int level);

int
t8_subelement_scheme_quad_c::t8_element_maxlevel (void) const
{
  return P4EST_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_subelement_scheme_quad_c::t8_element_child_eclass (int childid) const
/* *INDENT-ON* */

{
  T8_ASSERT (0 <= childid && childid < P4EST_CHILDREN);

  return T8_ECLASS_QUAD;
}

int
t8_subelement_scheme_quad_c::t8_element_level (const t8_element_t *elem) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_quad_with_subelements *) pquad_w_sub)->p4q.level;
}

static void
t8_element_copy_surround (const p4est_quadrant_t *q, p4est_quadrant_t *r)
{
  T8_QUAD_SET_TDIM (r, T8_QUAD_GET_TDIM (q));
  if (T8_QUAD_GET_TDIM (q) == 3) {
    T8_QUAD_SET_TNORMAL (r, T8_QUAD_GET_TNORMAL (q));
    T8_QUAD_SET_TCOORD (r, T8_QUAD_GET_TCOORD (q));
  }
}

void
t8_subelement_scheme_quad_c::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  const t8_quad_with_subelements *pquad_w_sub_source = (const t8_quad_with_subelements *) source;
  t8_quad_with_subelements *pquad_w_sub_dest = (t8_quad_with_subelements *) dest;

  const p4est_quadrant_t *q = &pquad_w_sub_source->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_dest->p4q;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (q == r && pquad_w_sub_source->transition_type == pquad_w_sub_dest->transition_type
      && pquad_w_sub_source->subelement_id == pquad_w_sub_dest->subelement_id) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *r = *q;

  t8_element_copy_subelement_values (source, dest);
  t8_element_copy_surround (q, r);
}

int
t8_subelement_scheme_quad_c::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 = (const t8_quad_with_subelements *) elem1;
  const t8_quad_with_subelements *pquad_w_sub_elem2 = (const t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *q = &pquad_w_sub_elem1->p4q;
  const p4est_quadrant_t *r = &pquad_w_sub_elem2->p4q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  int compare = p4est_quadrant_compare (q, r);

  if (compare == 0 && (t8_element_is_subelement (elem1) || t8_element_is_subelement (elem2))) {
    t8_debugf ("Caution, t8_element_compare is used with subelements.\n");
    if (t8_element_is_subelement (elem1) && t8_element_is_subelement (elem2)) {
      /* Caution: The compare function is used for two subelements. */

      if (pquad_w_sub_elem1->transition_type == pquad_w_sub_elem2->transition_type
          && pquad_w_sub_elem1->subelement_id == pquad_w_sub_elem2->subelement_id) {
        /* both subelements are identical */
        return 0;
      }
      /* return != 0 to avoid debug abortion in t8_ghost_add_remote */
      return 1;
    }
    else if (t8_element_is_subelement (elem1)) {
      return -1; /* elem1 is subelement and therefore smaller */
    }
    else if (t8_element_is_subelement (elem2)) {
      return 1; /* elem2 is subelement and therefore smaller */
    }
  }

  /* Note that for subelements, their parent quadrant is compared at this point */
  return compare;
}

void
t8_subelement_scheme_quad_c::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_parent = (t8_quad_with_subelements *) parent;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_parent->p4q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));

  if (t8_element_is_subelement (elem)) {
    pquad_w_sub_parent->p4q = pquad_w_sub_elem->p4q;
  }
  else {
    p4est_quadrant_parent (q, r);
  }

  /* the parent of any element will never be a subelement */
  t8_element_reset_subelement_values (parent);

  t8_element_copy_surround (q, r);
}

void
t8_subelement_scheme_quad_c::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_sibling = (t8_quad_with_subelements *) sibling;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_sibling->p4q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));

  p4est_quadrant_sibling (q, r, sibid);
  t8_element_copy_surround (q, r);
}

int
t8_subelement_scheme_quad_c::t8_element_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_QUAD_SUBELEMENT_FACES : P4EST_FACES);
}

int
t8_subelement_scheme_quad_c::t8_element_max_num_faces (const t8_element_t *elem) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return P4EST_FACES;
}

int
t8_subelement_scheme_quad_c::t8_element_num_children (const t8_element_t *elem) const
{
  /* Note that children of subelements equal the children of the parent quadrant. 
   * Therefore, the number of children of a subelement equals P4EST_CHILDREN */
  T8_ASSERT (t8_element_is_valid (elem));
  return P4EST_CHILDREN;
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_quad_c::t8_element_num_siblings (const t8_element_t *elem) const
/* *INDENT-ON* */

{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  if (pquad_w_sub->transition_type == 0)
    return P4EST_FACES;

  int num_hanging_faces = 0;
  int iface;
  for (
    iface = 0; iface < P4EST_FACES;
    iface++) { /* Count the number of ones of the binary transition type. This number equals the number of hanging faces. */
    num_hanging_faces += (pquad_w_sub->transition_type & (1 << iface)) >> iface;
  }

  return P4EST_CHILDREN + num_hanging_faces;
}

int
t8_subelement_scheme_quad_c::t8_element_num_face_children (const t8_element_t *elem, int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));

  /* if we use this scheme without set_transition, then we are only balanced and two neighbors are possible */
  return 2;
}

int
t8_subelement_scheme_quad_c::t8_element_neighbor_is_sibling (const t8_element_t *elem, const int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));

  if (face == 0 || face == 2) {
    return 1;
  }

  return 0;
}

/* *INDENT-OFF* */
int
t8_subelement_scheme_quad_c::t8_element_get_num_sibling_neighbors_at_face (const t8_element_t *elem,
                                                                           const int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (face == 0 || face == 2);

  return 1;
}

int
t8_subelement_scheme_quad_c::t8_element_get_transition_refine_identifier () const
{
  return T8_TRANSITION_CONFORMAL_QUAD_REFINE_FUNCTION;
}
/* *INDENT-ON* */

int
t8_subelement_scheme_quad_c::t8_element_get_face_corner (const t8_element_t *elem, int face, int corner) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
    /*
     *   2    f_2    3
     *     x -->-- x
     *     |       |
     *     ^       ^
     * f_0 |       | f_1
     *     x -->-- x
     *   0    f_3    1
     */

    T8_ASSERT (0 <= face && face < P4EST_FACES);
    T8_ASSERT (0 <= corner && corner < 4);

    return p4est_face_corners[face][corner];
  }
  else {
    int t8_face_corners_subelement[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
    /*
     *
     *         x - - - - - x 1
     *         | \    f0 / |
     *         |   \ 0 /   |
     *         x - - x  el | f1
     *         |   /   \   |
     *         | /    f2 \ |
     *         x - - x - - x 2
     *               
     * The vertecies of a subelement are enumerated clockwise, starting with the center vertex of the transition cell 
     */

    T8_ASSERT (0 <= face && face < T8_QUAD_SUBELEMENT_FACES);
    T8_ASSERT (0 <= corner && corner < 3);

    return t8_face_corners_subelement[face][corner];
  }
}

int
t8_subelement_scheme_quad_c::t8_element_get_corner_face (const t8_element_t *elem, int corner, int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  T8_ASSERT (0 <= face && face < 2);

  return p4est_corner_faces[corner][face];
}

void
t8_subelement_scheme_quad_c::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  /*
   *
   *         x - - - - - x        x - - x - - x 
   *         |           |        |     |     |
   *         |           |        |  2  |  3  |
   *         |   elem    |   =>   x - - x - - x
   *         |           |        |     |     |
   *         |           |        |  0  |  1  |
   *         x - - - - - x        x - - x - - x
   * 
   */
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_child = (t8_quad_with_subelements *) child;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_child->p4q;

  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);

  /* it should not be possible to construct a child of a subelement */
  // T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (child));
  T8_ASSERT (p4est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P4EST_QMAXLEVEL);

  T8_ASSERT (childid >= 0 && childid < P4EST_CHILDREN);

  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->level = q->level + 1;

  T8_ASSERT (p4est_quadrant_is_parent (q, r));

  t8_element_reset_subelement_values (child);

  t8_element_copy_surround (q, r);
}

/* *INDENT-OFF* */
void
t8_subelement_scheme_quad_c::t8_element_get_sibling_neighbor_in_transition_cell (const t8_element_t *elem,
                                                                                 const int face,
                                                                                 const int num_neighbors,
                                                                                 t8_element_t *neighbor_at_face[],
                                                                                 int *neigh_face[])
{
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_neighbor_is_sibling (elem, face));
  T8_ASSERT (num_neighbors == 1);
  T8_ASSERT (t8_element_is_valid (neighbor_at_face[0]));

  /* If face = 0, then the sibling subelement neighbor is the next subelement in counter clockwise enumeration,
   * if face = 2, then it is the sibling subelement neighbor in  clockwise enumeration. */
  t8_element_copy (elem, neighbor_at_face[0]);

  t8_quad_with_subelements *pquad_w_sub_neighbor_at_face = (t8_quad_with_subelements *) neighbor_at_face[0];

  int num_siblings = t8_element_num_siblings (elem);

  if (face == 0) {
    /* adjust subelement id counter clockwise */
    if (pquad_w_sub_neighbor_at_face->subelement_id == 0) {
      pquad_w_sub_neighbor_at_face->subelement_id += num_siblings - 1;
    }
    else {
      pquad_w_sub_neighbor_at_face->subelement_id -= 1;
    }
  }
  else {
    /* adjust subelement id clockwise */
    if (pquad_w_sub_neighbor_at_face->subelement_id == num_siblings - 1) {
      pquad_w_sub_neighbor_at_face->subelement_id = 0;
    }
    else {
      pquad_w_sub_neighbor_at_face->subelement_id += 1;
    }
  }

  /* return dual face with resprect to neighboring sibling subelement */
  /* Compute the face number as seen from elem.
   *  0 -> 2    2 -> 0
   */
  *neigh_face[0] = subelement_face_dual[face];
}
/* *INDENT-ON* */

void
t8_subelement_scheme_quad_c::t8_element_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  /* if elem is a subelement, then this function will construct the children of its parent p4est quadrant */
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements **pquad_w_sub_children = (t8_quad_with_subelements **) c;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;

  int ichild;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == P4EST_CHILDREN);

#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < P4EST_CHILDREN; i++) {
      T8_ASSERT (t8_element_is_valid (c[i]));
    }
  }
#endif

  /* set coordinates and levels of the children */
  p4est_quadrant_children (q, &pquad_w_sub_children[0]->p4q, &pquad_w_sub_children[1]->p4q,
                           &pquad_w_sub_children[2]->p4q, &pquad_w_sub_children[3]->p4q);

  for (ichild = 0; ichild < P4EST_CHILDREN; ++ichild) {
    t8_element_reset_subelement_values (c[ichild]);
    t8_element_copy_surround (q, &pquad_w_sub_children[ichild]->p4q);
  }
}

int
t8_subelement_scheme_quad_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? pquad_w_sub->subelement_id : p4est_quadrant_child_id (q));
}

int
t8_subelement_scheme_quad_c::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return p4est_quadrant_ancestor_id (q, level);
}

int
t8_subelement_scheme_quad_c::t8_element_is_family (t8_element_t *const *fam) const
{
  /* Note that this test is very rudimentary, especially when there subelements are in fam */
  t8_quad_with_subelements **pquad_w_sub_family = (t8_quad_with_subelements **) fam;

#ifdef T8_ENABLE_DEBUG
  {
    int i;
    int num_siblings = t8_element_num_siblings (fam[0]);
    for (i = 0; i < num_siblings; i++) {
      T8_ASSERT (t8_element_is_valid (fam[i]));
    }
  }
#endif

  /* Subelements can not be refined into other elements of a higher level. 
   * So if the first element of fam is a subelement, we assume that the following num_siblings 
   * many elements are its siblings and therefore form a family. */
  if (pquad_w_sub_family[0]->transition_type != 0) {
    return 1;
  }
  /* If the first element of fam is no subelement we check the following elements of fam */
  else {
    /* If any of the following elements is a subelement, then they can not form a family */
    if (pquad_w_sub_family[1]->transition_type != 0 || pquad_w_sub_family[2]->transition_type != 0
        || pquad_w_sub_family[3]->transition_type != 0) {
      return 0;
    }
    /* If all elements of fam are no subelements, then we can use the p4est check is_family */
    else {
      return p4est_quadrant_is_family (&pquad_w_sub_family[0]->p4q, &pquad_w_sub_family[1]->p4q,
                                       &pquad_w_sub_family[2]->p4q, &pquad_w_sub_family[3]->p4q);
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << P4EST_DIM * level);

  p4est_quadrant_set_morton (q, level, id);
  T8_QUAD_SET_TDIM (q, 2);
}

t8_linearidx_t
t8_subelement_scheme_quad_c::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* Note that the id of a subelement equals the id of its parent quadrant.
   * Therefore, the binary search (for example used in the leaf_face_neighbor function) 
   * will find a random subelement of the transition cell which might not be the desired neighbor of a given element. */

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  return p4est_quadrant_linear_id (q, level);
}

void
t8_subelement_scheme_quad_c::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_desc = (t8_quad_with_subelements *) desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_desc->p4q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  p4est_quadrant_first_descendant (q, r, level);
  T8_QUAD_SET_TDIM (r, 2);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a quad element to a subelement does not change the 
   * p4est quadrant. Therefore, we are constructing the last descendant of the parent 
   * quad element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);
}

void
t8_subelement_scheme_quad_c::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_desc = (t8_quad_with_subelements *) desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_desc->p4q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  p4est_quadrant_last_descendant (q, r, level);
  T8_QUAD_SET_TDIM (r, 2);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a quad element to a subelement does not change the 
   * p4est quadrant. Therefore, we are constructing the last descendant of the parent 
   * quad element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);
}

void
t8_subelement_scheme_quad_c::t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 = (const t8_quad_with_subelements *) elem1;
  t8_quad_with_subelements *pquad_w_sub_elem2 = (t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *q = &pquad_w_sub_elem1->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_elem2->p4q;

  t8_linearidx_t id;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= t8_element_level (elem1) && t8_element_level (elem1) <= P4EST_QMAXLEVEL);

  id = p4est_quadrant_linear_id (q, t8_element_level (elem1));
  T8_ASSERT (id + 1 < ((t8_linearidx_t) 1) << P4EST_DIM * t8_element_level (elem1));
  t8_element_reset_subelement_values (elem2);
  p4est_quadrant_set_morton (r, t8_element_level (elem1), id + 1);
  t8_element_copy_surround (q, r);
}

void
t8_subelement_scheme_quad_c::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                             t8_element_t *nca) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 = (const t8_quad_with_subelements *) elem1;
  const t8_quad_with_subelements *pquad_w_sub_elem2 = (const t8_quad_with_subelements *) elem2;
  t8_quad_with_subelements *pquad_w_sub_nca = (t8_quad_with_subelements *) nca;

  const p4est_quadrant_t *q1 = &pquad_w_sub_elem1->p4q;
  const p4est_quadrant_t *q2 = &pquad_w_sub_elem2->p4q;
  p4est_quadrant_t *r = &pquad_w_sub_nca->p4q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
#if 0
  /* TODO: This assertions throws an error since it expects a 3D hex.
   *       this does not make sense. investigate. */
  T8_ASSERT (t8_element_surround_matches (q1, q2));
#endif

  /* In case of subelements, we use the parent quadrant and construct nca of the parent quadrant */
  t8_element_reset_subelement_values (nca);
  p4est_nearest_common_ancestor (q1, q2, r);
  t8_element_copy_surround (q1, r);
}

t8_element_shape_t
t8_subelement_scheme_quad_c::t8_element_face_shape (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_ECLASS_LINE;
}

void
t8_subelement_scheme_quad_c::t8_element_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                          int num_children, int *child_indices) const
{
#ifdef T8_ENABLE_DEBUG
  {
    int i;
    for (i = 0; i < num_children; i++) {
      T8_ASSERT (t8_element_is_valid (children[i]));
    }
  }
#endif
  /* This function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (num_children == t8_element_num_face_children (elem, face));

  /*
   * Compute the child id of the first and second child at the face.
   *
   *            3
   *
   *      x - - x - - x           This picture shows a refined quadrant
   *      |     |     |           with child_ids and the label for the faces.
   *      | 2   | 3   |           For example for face 2 (bottom face) we see
   * 0    x - - x - - x   1       first_child = 0 and second_child = 1.
   *      |     |     |
   *      | 0   | 1   |
   *      x - - x - - x
   *
   *            2
   */

  T8_ASSERT (num_children == 2);
  int first_child;
  int second_child;
  /* TODO: Think about a short and easy bitwise formula. */
  switch (face) {
  case 0:
    first_child = 0;
    second_child = 2;
    break;
  case 1:
    first_child = 1;
    second_child = 3;
    break;
  case 2:
    first_child = 0;
    second_child = 1;
    break;
  case 3:
    first_child = 2;
    second_child = 3;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* From the child ids we now construct the children at the faces. */
  /* We have to revert the order and compute second child first, since
   * the usage allows for elem == children[0].
   */
  this->t8_element_child (elem, second_child, children[1]);
  this->t8_element_child (elem, first_child, children[0]);
  if (child_indices != NULL) {
    child_indices[0] = first_child;
    child_indices[1] = second_child;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 1);
    return t8_element_face_parent_face (elem, face);
  }
  else {
    /* For quadrants the face enumeration of children is the same as for the parent. */
    return face;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (face >= -1 && face <= P4EST_FACES);

  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  int child_id;

  if (face == -1) {
    return -1;
  }

  /* For subelements we need to adjust the output of this function.
   * A subelements face is a subface of the parent quadrant (the transition cell) if and only if the face number is 1. */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* In this case the face is a subface of the parent. We use the location function in order
       * to determine which of the parents faces intersects the subelements face. */
      int location[3] = {};
      t8_element_get_location_of_subelement (elem, location);

      /* subelements in location are enumerated clockwise (not as quadrant faces) */
      return subelement_location_to_parent_face[location[0]];
    }
    else {
      return -1;
    }
  }

  if (q->level == 0) {
    return face;
  }
  /* Determine whether face is a subface of the parent.
   * This is the case if the child_id matches one of the faces corners */
  child_id = p4est_quadrant_child_id (q);
  if (child_id == p4est_face_corners[face][0] || child_id == p4est_face_corners[face][1]) {
    return face;
  }
  return -1;
}

void
t8_subelement_scheme_quad_c::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation,
                                                        int sign, int is_smaller_face) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 = (const t8_quad_with_subelements *) elem1;
  t8_quad_with_subelements *pquad_w_sub_elem2 = (t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *qin = &pquad_w_sub_elem1->p4q;
  p4est_quadrant_t *p = &pquad_w_sub_elem2->p4q;

  const p4est_quadrant_t *q;
  p4est_qcoord_t h = P4EST_QUADRANT_LEN (qin->level);
  p4est_qcoord_t x = qin->x; /* temp storage for x coordinate in case elem1 = elem 2 */

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= orientation && orientation < P4EST_FACES);

  if (sign) {
    /* The tree faces have the same topological orientation, and
     * thus we have to perform a coordinate switch. */
    /* We use p as storage, since elem1 and elem2 are allowed to
     * point to the same quad */
    q = (const p4est_quadrant_t *) p;
    t8_element_copy_surround (qin, (p4est_quadrant_t *) q);
    ((p4est_quadrant_t *) q)->x = qin->y;
    ((p4est_quadrant_t *) q)->y = x;
    x = q->x; /* temp storage in case elem1 = elem 2 */
  }
  else {
    q = qin;
  }

  p->level = q->level;
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

  switch (orientation) {
  case 0: /* Nothing to do */
    p->x = q->x;
    p->y = q->y;
    break;
  case 1:
    p->x = P4EST_ROOT_LEN - q->y - h;
    p->y = x;
    break;
  case 2:
    p->x = q->y;
    p->y = P4EST_ROOT_LEN - x - h;
    break;
  case 3:
    p->x = P4EST_ROOT_LEN - q->x - h;
    p->y = P4EST_ROOT_LEN - q->y - h;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  T8_QUAD_SET_TDIM (p, 2);

  t8_element_reset_subelement_values (elem2);
}

int
t8_subelement_scheme_quad_c::t8_element_extrude_face (const t8_element_t *face, const t8_eclass_scheme_c *face_scheme,
                                                      t8_element_t *elem, int root_face) const
{
  /* build (extrude) elem from a given face element */
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t *q = &pquad_w_sub->p4q;

  const t8_dline_t *l = (const t8_dline_t *) face;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE (face_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (face_scheme->t8_element_is_valid (elem));
  T8_ASSERT (0 <= root_face && root_face < P4EST_FACES);

  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *        f_2
   *     x -->-- x
   *     |       |
   *     ^       ^
   * f_0 |       | f_1
   *     x -->-- x
   *        f_3
   *
   * The arrows >,^ denote the orientation of the faces.
   * We need to scale the coordinates since a root line may have a different
   * length than a root quad.
   */
  q->level = l->level;
  switch (root_face) {
  case 0:
    q->x = 0;
    q->y = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 1:
    q->x = P4EST_LAST_OFFSET (q->level);
    q->y = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 2:
    q->x = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = 0;
    break;
  case 3:
    q->x = ((int64_t) l->x * P4EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = P4EST_LAST_OFFSET (q->level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  t8_element_reset_subelement_values (elem);
  /* We return the face of q at which we extruded. This is the same number
   * as root_face. */
  return root_face;
}

int
t8_subelement_scheme_quad_c::t8_element_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* If elem is a subelement, then this function should only be called together with 
   * face = 1 since other faces will never intersect a tree face. */
  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 1);

    return t8_element_face_parent_face (elem, face);
  }
  else {
    T8_ASSERT (0 <= face && face < P4EST_FACES);
    /* For quadrants the face and the tree face number are the same. */
    return face;
  }
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_subelement_scheme_quad_c::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                               t8_element_t *first_desc, int level) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_first_desc = (t8_quad_with_subelements *) first_desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *desc = &pquad_w_sub_first_desc->p4q;

  int first_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  /* Get the first corner of q that belongs to face */
  first_face_corner = p4est_face_corners[face][0];
  /* Construce the descendant in that corner */
  p4est_quadrant_corner_descendant (q, desc, first_face_corner, level);
  t8_element_reset_subelement_values (first_desc);
}

/** Construct the last descendant of an element that touches a given face.   */
void
t8_subelement_scheme_quad_c::t8_element_last_descendant_face (const t8_element_t *elem, int face,
                                                              t8_element_t *last_desc, int level) const
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_last_desc = (t8_quad_with_subelements *) last_desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *desc = &pquad_w_sub_last_desc->p4q;

  int last_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (!t8_element_is_subelement (last_desc));
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  /* Get the last corner of q that belongs to face */
  last_face_corner = p4est_face_corners[face][1];
  /* Construce the descendant in that corner */
  p4est_quadrant_corner_descendant (q, desc, last_face_corner, level);
  t8_element_reset_subelement_values (last_desc);
}

void
t8_subelement_scheme_quad_c::t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                       const t8_eclass_scheme_c *boundary_scheme) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  t8_dline_t *l = (t8_dline_t *) boundary;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE (boundary_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));

  if (!t8_element_is_subelement (elem)) {
    T8_ASSERT (0 <= face && face < P4EST_FACES);
    /* The level of the boundary element is the same as the quadrant's level */
    l->level = q->level;
    /*
     * The faces of the quadrant are enumerated like this:
     *        f_2
     *     x ---- x
     *     |      |
     * f_0 |      | f_1
     *     x ---- x
     *        f_3
     *
     * If face = 0 or face = 1 then l->x = q->y
     * if face = 2 or face = 3 then l->x = q->x
     */
    l->x = ((face >> 1 ? q->x : q->y) * ((int64_t) T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN);
  }
  else {
    /* face number 1 is the only face of a subelement that points outward of the transition cell */
    T8_ASSERT (face == 1);
    /* boundary faces of subelements:
     *
     *         x - - - - - x
     *         | \       / |
     *         |   \   /   |
     *         x - - x  e2 | f1
     *         |e1 /   \   |
     *      f1 | /       \ |
     *         x - - x - - x
     *               
     * for a split subelement (e1), the boundary face has a higher level
     * for a non split element (e2), the boundary face has the same level. 
     */

    int location[3]
      = {}; /* location = {location of subelement (face number of transition cell), split, first or second element if split} */
    t8_element_get_location_of_subelement (elem, location);
    int split = location[1];
    int second = location[2];

    if (split) { /* if the subelement lies at a split face */
      l->level = q->level + 1;
      int len = P4EST_QUADRANT_LEN (pquad_w_sub->p4q.level + 1);
      if (second) {             /* second subelement */
        if (location[0] == 0) { /* left face */
          l->x = q->y + len;
        }
        else if (location[0] == 1) { /* upper face */
          l->x = q->x + len;
        }
        else if (location[0] == 2) { /* right face */
          l->x = q->y;
        }
        else { /* lower face */
          l->x = q->x;
        }
      }
      else {                    /* first subelement */
        if (location[0] == 0) { /* left face */
          l->x = q->y;
        }
        else if (location[0] == 1) { /* upper face */
          l->x = q->x;
        }
        else if (location[0] == 2) { /* right face */
          l->x = q->y + len;
        }
        else { /* lower face */
          l->x = q->x + len;
        }
      }
    }
    else { /* if the subelement is not split */
      l->level = q->level;
      if (location[0] == 0) { /* left face */
        l->x = q->y;
      }
      else if (location[0] == 1) { /* upper face */
        l->x = q->x;
      }
      else if (location[0] == 2) { /* right face */
        l->x = q->y;
      }
      else { /* lower face */
        l->x = q->x;
      }
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_boundary (const t8_element_t *elem, int min_dim, int length,
                                                  t8_element_t **boundary) const
{
  SC_ABORT ("Not implemented\n");
#if 0
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif
  int                 iface;

  T8_ASSERT (length ==
             t8_eclass_count_boundary (T8_ECLASS_QUAD, min_dim, per_eclass));

  T8_ASSERT (length == P4EST_FACES);
  for (iface = 0; iface < P4EST_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
#endif
}

int
t8_subelement_scheme_quad_c::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  p4est_qcoord_t coord;

  /* In case of a subelement, we need to change its face number to the face number of the parent quad */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* adjust face of subelement to face of parent */
      face = t8_element_face_parent_face (elem, face);
    }
    else { /* in case of a subelement and face 0 or 2 the face is no subface of the root boundary */
      return 0;
    }
  }

  T8_ASSERT (0 <= face && face < P4EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   */
  coord = face >> 1 ? q->y : q->x;
  /* If face is 0 or 2 check against 0.
   * If face is 1 or 3  check against LAST_OFFSET */
  return coord == (face & 1 ? P4EST_LAST_OFFSET (q->level) : 0);
}

int
t8_subelement_scheme_quad_c::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                              int *neigh_face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < P4EST_FACES);

  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_neigh = (t8_quad_with_subelements *) neigh;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t *n = &pquad_w_sub_neigh->p4q;

  /* In case of a subelement one should construct the face neighbor of the face-corresponding child quadrant
   * of the subelements parent quadrant. Therefore we might want to adjust the level  and adapt the
   * anchor node. */
  if (t8_element_is_subelement (elem)) { /* if elem is a subelement */

    T8_ASSERT (0 <= face && face < T8_QUAD_SUBELEMENT_FACES);

    if (face == 0) { /* in this case the face neighbor of the subelement is a sibling */
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 2) { /* in this case the face neighbor of the subelement is a sibling */
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 1) { /* in this case the face neighbor is no sibling */
      int location[3] = {};
      t8_element_get_location_of_subelement (elem, location);

      /* setting the anchor node of the neighbor element */
      n->x = q->x;
      n->y = q->y;

      /* half the side length of the transition cell of the subelement */
      const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);

      int split = location[1];
      int second = location[2];

      /* we need to take into account whether the subelement is split or not */
      if (split) { /* split */

        /* increase the level by one */
        n->level = q->level + 1;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          if (!second) {
            n->x = q->x - shift;
          }
          else {
            n->x = q->x - shift;
            n->y = q->y + shift;
          }
        }
        else if (location[0] == 2) { /* right face */
          if (!second) {
            n->x = q->x + 2 * shift;
            n->y = q->y + shift;
          }
          else {
            n->x = q->x + 2 * shift;
          }
        }
        else if (location[0] == 3) { /* lower face */
          if (!second) {
            n->x = q->x + shift;
            n->y = q->y - shift;
          }
          else {
            n->y = q->y - shift;
          }
        }
        else { /* upper face */
          if (!second) {
            n->y = q->y + 2 * shift;
          }
          else {
            n->x = q->x + shift;
            n->y = q->y + 2 * shift;
          }
        }
      }

      else { /* not split */
        /* level stays the same */
        n->level = q->level;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          n->x = q->x - 2 * shift;
        }
        else if (location[0] == 2) { /* right face */
          n->x = q->x + 2 * shift;
        }
        else if (location[0] == 3) { /* lower face */
          n->y = q->y - 2 * shift;
        }
        else { /* upper face */
          n->y = q->y + 2 * shift;
        }
      }
    }
  }
  else { /* if elem is no subelement */
    /* Directly construct the face neighbor */
    p4est_quadrant_face_neighbor (q, face, n);
  }

  t8_element_reset_subelement_values (neigh);

  T8_QUAD_SET_TDIM (n, 2);

  /* In the following we set the dual faces of our element at the given face. */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* return dual face with respect to neighboring quad element */
      int location[3] = {};
      t8_element_get_location_of_subelement (elem, location);
      /* if the face is pointing outwards, then we set the face equal to the transition cell face and determine its dual face.
       * Compute the face number as seen from q.
       *  0 -> 1    1 -> 2    2 -> 0    3 -> 3
       */
      *neigh_face = subelement_location_to_parent_dual_face[location[0]];
    }
    else {
      T8_ASSERT (face == 0 || face == 2);
      /* return dual face with resprect to neighboring sibling subelement (note that the constructed neigh is NOT a subelement but the parent quad) */
      /* Compute the face number as seen from q.
       *  0 -> 2    2 -> 0
       */
      *neigh_face = subelement_face_dual[face];
    }
  }
  else {
    /* Compute the face number as seen from q.
     *  0 -> 1    1 -> 0    2 -> 3    3 -> 2
     */
    T8_ASSERT (neigh_face != NULL);
    *neigh_face = p4est_face_dual[face];
  }

  /* return true if neigh is inside the root */
  return p4est_quadrant_is_inside_root (n);
}

void
t8_subelement_scheme_quad_c::t8_element_anchor (const t8_element_t *elem, int coord[3]) const
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
  T8_QUAD_SET_TDIM (q, 2);
}

int
t8_subelement_scheme_quad_c::t8_element_root_len (const t8_element_t *elem) const
{
  return P4EST_ROOT_LEN;
}

int
t8_subelement_scheme_quad_c::t8_element_refines_irregular () const
{
  /* In general, subelements do not refine regularly */
  return 1;
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

void
t8_subelement_scheme_quad_c::t8_element_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                          const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (num_coords > 0);
  T8_ASSERT (ref_coords != NULL);
  T8_ASSERT (out_coords != NULL);
  T8_ASSERT (t8_element_is_valid (elem));
  
  const t8_quad_with_subelements *element = (t8_quad_with_subelements *) elem;
  if (!t8_element_is_subelement (elem)) {
    const t8_element_t * quad = (const t8_element_t *)&element->p4q;
    // We call the default quad reference coord function for the quad coordinates.
    // We could think about doing this in more cases. For that it would
    // be beneficial to store a t8_default_scheme_quad_c pointer as member variable
    // of the t8_subelement_scheme_quad_c class.
    default_quad_scheme.t8_element_reference_coords (quad, ref_coords, num_coords, out_coords);
  }
  else {
    /* This element is a subelement and hence a triangle. */
    T8_ASSERT (t8_element_shape (elem) == T8_ECLASS_TRIANGLE);

    /* We first convert the reference coordinates to barycentric
     * coordinates.
     * Since the reference triangle has vertices
     * (0,0) (1,0) and (1,1)
     * A vector (x,y) in that triangle has barycentric
     * coordinates (0, x-y, y).
     * This is since 0 * (0,0) + (x-y) * (1,0) + y * (1,1) = (0 + x - y + y, y) = (x,y)
     *
     * After that we use t8_geom_triangular_interpolation to compute the actual value.
     */
    constexpr int num_vertices = 3;  // Should use t8_eclass_num_vertices[T8_ECLASS_TRIANGLE]; but is not constexpr
    constexpr int num_vertex_coords = num_vertices * 3;
    double vertex_coords[num_vertex_coords]; // Stores the triangle's vertices
    for (int ivertex = 0;ivertex < num_vertices;++ivertex) {
      t8_element_vertex_reference_coords (elem, ivertex, vertex_coords + 3*ivertex);
    }

    // For each incoming coordinate compute the barycentric coordinates and do the interpolation
    for (size_t icoord = 0;icoord < num_coords;++icoord) {
      // The ref_coords are always 3 dimensional - even though we do not use the 3rd entry.
      // The out_coords are 2 dimensional.
      const int offset_3d = icoord * 3; // offset for 3 dim ref_coords when iterating over points
      const int offset_2d = icoord * 2; // offset for 2 dim out_coords when iterating over points
      const double ref_x = ref_coords[0 + offset_3d];
      const double ref_y = ref_coords[1 + offset_3d];
      const double barycentric_coeff[3] = {0, ref_x - ref_y, ref_y};
      t8_geom_triangular_interpolation (barycentric_coeff, vertex_coords, 3, 2, out_coords + offset_2d);
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_reference_coords (const t8_element_t *t, int vertex,
                                                                 double coords[]) const
{
  int coords_int[2] = {};
  t8_element_vertex_coords (t, vertex, coords_int);

  /* We divide the integer coordinates by the root length of the quad
   * to obtain the reference coordinates. */
  coords[0] = (double) coords_int[0] / (double) P4EST_ROOT_LEN;
  coords[1] = (double) coords_int[1] / (double) P4EST_ROOT_LEN;
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q1 = &pquad_w_sub->p4q;

  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
    int len;

    T8_ASSERT (0 <= vertex && vertex < 4);
    /* Get the length of the quadrant */
    len = P4EST_QUADRANT_LEN (q1->level);

    /* Compute the x and y coordinates of the vertex depending on the
     * vertex number */
    coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
    coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
  }
  else {
    t8_element_vertex_coords_of_subelement (elem, vertex, coords);
  }
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_coords_of_subelement (const t8_element_t *elem, int vertex,
                                                                     int coords[]) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q1 = &pquad_w_sub->p4q;

  int len;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (vertex >= 0 && vertex < T8_QUAD_SUBELEMENT_FACES); /* all subelements are triangles */

  /* get the length of the current quadrant */
  len = P4EST_QUADRANT_LEN (q1->level);

  /* Compute the x and y coordinates of subelement vertices, depending on the transition type, id and vertex number 
   * (faces enumerated clockwise, starting at the center of the transition cell): 
   *
   *               f1                      V1
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - + - - x f2  -->   + - - x 
   *         | 0 / | \ 4 |           V0    V2
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f3
   * 
   * In this example, the below location array would contain the values [2, 1, 1] 
   * (second face, split, first subelement at this face) */

  /* get location information of the given subelement */
  int location[3] = {};
  t8_element_get_location_of_subelement (elem, location);

  /* the face number, the subelement is adjacent to */
  int face_number = location[0];
  /* = 1, if the adjacent face is split and = 0, if not */
  int split = location[1];
  /* = 0, if the subelement is the first (of two) subelements, at the adjacent face and = 1 if it is the second */
  int sub_face_id = location[2];

  /* Check, whether the get_location function provides meaningful location data */
  T8_ASSERT (face_number == 0 || face_number == 1 || face_number == 2 || face_number == 3);
  T8_ASSERT ((split == 0 && sub_face_id == 0) || (split == 1 && (sub_face_id == 0 || sub_face_id == 1)));

  coords[0] = q1->x;
  coords[1] = q1->y;

  /* using the location data to determine vertex coordinates */
  if (vertex == 0) { /* vertex 0 (the first vertex always equals the center of the element) */
    coords[0] += len / 2;
    coords[1] += len / 2;
  }                       /* end of vertex == 0 */
  else if (vertex == 1) { /* vertex 1 */
    if (face_number == 0) {
      if (split && sub_face_id) {
        coords[1] += len / 2;
      }
    }
    else if (face_number == 1) {
      coords[1] += len;
      if (split && sub_face_id) {
        coords[0] += len / 2;
      }
    }
    else if (face_number == 2) {
      coords[0] += len;
      coords[1] += len;
      if (split && sub_face_id) {
        coords[1] -= len / 2;
      }
    }
    else {
      coords[0] += len;
      if (split && sub_face_id) {
        coords[0] -= len / 2;
      }
    }
  }      /* end of vertex == 1 */
  else { /* vertex 2 */
    if (face_number == 0) {
      coords[1] += len;
      if (split && !sub_face_id) {
        coords[1] -= len / 2;
      }
    }
    else if (face_number == 1) {
      coords[0] += len;
      coords[1] += len;
      if (split && !sub_face_id) {
        coords[0] -= len / 2;
      }
    }
    else if (face_number == 2) {
      coords[0] += len;
      if (split && !sub_face_id) {
        coords[1] += len / 2;
      }
    }
    else {
      if (split && !sub_face_id) {
        coords[0] += len / 2;
      }
    }
  } /* end of vertex == 2 */
}

void
t8_subelement_scheme_quad_c::t8_element_to_transition_cell (const t8_element_t *elem, int type, t8_element_t *c[])
{
  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements **pquad_w_sub_subelement = (t8_quad_with_subelements **) c;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;

  /* this function should not be callable by subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (type >= 0 && type <= T8_SUB_QUAD_MAX_TRANSITION_TYPE);

  int num_subelements = t8_element_get_number_of_subelements (type);

#ifdef T8_ENABLE_DEBUG
  {
    int j;
    for (j = 0; j < num_subelements; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif

  /* get the length of a children-quadrant */
  const int8_t level = (int8_t) (q->level);

  T8_ASSERT (p4est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P4EST_QMAXLEVEL);

  /* Setting the parameter values for different subelements. 
   * The different subelement types (up to rotation) are:
   *                               
   *      x - - - - - - x         x - - - - - x        x - - - - - x        x - - - - - x        x - - x - - x        x - - x - - x
   *      |             |         | \   2   / |        | \       / |        | \       / |        | \   |   / |        | \   |   / |
   *      |             |         | 1 \   /   |        |   \   /   |        |   \   /   |        |   \ | /   |        |   \ | /   |
   *      |             |   -->   x - - X   3 |   or   x - - x     |   or   x - - x - - x   or   x - - x - - x   or   x - - x - - x
   *      |             |         | 0 /   \   |        |   / | \   |        |   /   \   |        |   /   \   |        |   / | \   |
   *      | elem        |         | /   4   \ |        | /   |   \ |        | /       \ |        | /       \ |        | /   |   \ |
   *      + - - - - - - x         x - - - - - x        x - - x - - x        x - - - - - x        x - - - - - x        x - - x - - x
   *           
   * Sub_ids are counted clockwise, starting with the (lower) left subelement with id 0.                    
   * Note, that we do not change the p4est quadrant. */

  int sub_id_counter = 0;
  for (sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
    pquad_w_sub_subelement[sub_id_counter]->p4q.x = q->x;
    pquad_w_sub_subelement[sub_id_counter]->p4q.y = q->y;
    pquad_w_sub_subelement[sub_id_counter]->p4q.level = level;

    pquad_w_sub_subelement[sub_id_counter]->transition_type = type;
    pquad_w_sub_subelement[sub_id_counter]->subelement_id = sub_id_counter;

    T8_ASSERT (t8_element_is_valid (c[sub_id_counter]));
    t8_element_copy_surround (q, &pquad_w_sub_subelement[sub_id_counter]->p4q);
  }
}

int
t8_subelement_scheme_quad_c::t8_element_get_number_of_subelements (int transition_type) const
{
  /* we could return 0 for transition type 0 but we will assert this case for safety reasons */
  T8_ASSERT (transition_type != 0);

  /* consider transition_type 13 = 1101 in base two -> there are 4 + (1+1+0+1) = 7 subelements */
  int num_hanging_faces = 0;
  int ichild;
  for (
    ichild = 0; ichild < P4EST_FACES;
    ichild++) { /* Count the number of ones of the binary transition type. This number equals the number of hanging faces. */
    num_hanging_faces += (transition_type & (1 << ichild)) >> ichild;
  }

  /* The number of subelements equals the number of neighbours: */
  return P4EST_FACES + num_hanging_faces;
}

void
t8_subelement_scheme_quad_c::t8_element_get_location_of_subelement (const t8_element_t *elem, int location[]) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  /* this function only works for subelements */
  T8_ASSERT (t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  /* Consider the following transition cell of type 13:
   *            
   *              f0                         1
   *        x - - x - - x              x - - x - - x           
   *        |           |              | \ 2 | 3 / |           faces:                                                      f3   f2   f1   f0
   *        |           |              | 1 \ | / 4 |           binary code:                                                 1    1    0    1   (=13)
   *     f3 x           x f2   -->   1 x - - x - - x 1   -->   rearrange binaries s.t. the faces are enumerated clockwise:  1    1    1    0
   *        |           |              | 0 /   \ 5 |           number subelements at face:                                  2    2    2    1
   *        | elem      |              | /   6   \ |           consider sub_id 3:                                                x -> second subelement on the upper face
   *        + - - - - - x              x - - - - - x
   *              f1                         0
   *           
   * We will use the binary representation to determine the location of the given subelement. 
   * 
   * We need to know: 
   *     i)   the face number of the first vertex (values: {0,1,2,3}).
   *     ii)  whether this face is split in half (values: {0,1}).
   *     iii) if the subelement is the first or second subelement at the face (values: {0,1}).
   * 
   * These information are then saved in the location array which will be used by the element_vertex function, 
   * to automatically determine the vertex coordinates of the given subelement. 
   * 
   * The location array for the above example would be {1,1,1} (upper face, split = true, second subelement at the upper face). */

  /* 1) convert the transition type from a decimal to a binary representation */
  int type = pquad_w_sub->transition_type;
  int binary_array[P4EST_FACES] = {};

  int iface;

  for (
    iface = 0; iface < P4EST_FACES;
    iface++) { /* need an array with 4 elements to store all subelement types of the quad scheme from 1 to 15 ({0,0,0,1} to {1,1,1,1}) */
    binary_array[(P4EST_FACES - 1) - iface] = (type & (1 << iface)) >> iface;
  } /* we now got a binary representation of the transition type, bitwise stored in an array */

  /* 2) rearrange the binary representation to be in clockwise order */
  int binary_array_temp[P4EST_FACES] = {};

  for (iface = 0; iface < P4EST_FACES; iface++) { /* copying the binary array */
    binary_array_temp[iface] = binary_array[iface];
  }

  for (iface = 0; iface < P4EST_FACES; iface++) { /* bringing the entries of binary array into clockwise order */
    binary_array[iface] = binary_array_temp[subelement_location_to_parent_face[iface]];
  }

  /* 3) use the rearranged binary representation, and the sub_id to determine the location of the subelement and store these information in an array */
  /*     3.1) location[0] -> the face_number, the subelement is adjacent to */
  /*     3.2) location[1] -> if the face is split or not */
  /*     3.3) location[2] -> if the subelement is the first or second subelement of the face (always the first, if the face is not split) */
  T8_ASSERT (pquad_w_sub->subelement_id < t8_element_get_number_of_subelements (pquad_w_sub->transition_type));

  int sub_id = pquad_w_sub->subelement_id;
  int sub_face_id;
  int face_number = -1;
  int split;

  int cum_neigh_array[P4EST_FACES] = {};

  /* construct a cumulative array of the number of neighbors from face 0 to face 3 */
  cum_neigh_array[0] = binary_array[0] + 1;
  cum_neigh_array[1] = cum_neigh_array[0] + binary_array[1] + 1;
  cum_neigh_array[2] = cum_neigh_array[1] + binary_array[2] + 1;
  cum_neigh_array[3] = cum_neigh_array[2] + binary_array[3] + 1;

  /* 3.1) we can use the cumulative array to determine the face number of the given subelement */
  if (sub_id < cum_neigh_array[0]) {
    face_number = 0;
  }
  else {
    for (iface = 0; iface < P4EST_FACES - 1; ++iface) {
      if (sub_id >= cum_neigh_array[iface] && sub_id < cum_neigh_array[iface + 1]) {
        face_number = iface + 1;
        break;
      }
    }
  }

  /* make sure that a face_number has been found */
  T8_ASSERT (face_number >= 0);

  /* 3.2) determine, whether the face is split or not */
  if (binary_array[face_number] == 0) {
    split = 0; /* the face is not split */
  }
  else {
    split = 1; /* the face is split */
  }

  /* 3.3) determine, whether the subelement is the first or the second subelement at the face */
  if (sub_id + 1 == cum_neigh_array[face_number] && split == 1) {
    sub_face_id = 1; /* second subelement */
  }
  else {
    sub_face_id = 0; /* first subelement */
  }

  location[0] = face_number;
  location[1] = split;
  location[2] = sub_face_id;
}

void
t8_subelement_scheme_quad_c::t8_element_reset_subelement_values (t8_element *elem) const
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;

  pquad_w_sub->transition_type = 0;
  pquad_w_sub->subelement_id = 0;
}

void
t8_subelement_scheme_quad_c::t8_element_copy_subelement_values (const t8_element *source, t8_element *dest) const
{
  const t8_quad_with_subelements *pquad_w_sub_source = (const t8_quad_with_subelements *) source;
  t8_quad_with_subelements *pquad_w_sub_dest = (t8_quad_with_subelements *) dest;

  pquad_w_sub_dest->transition_type = pquad_w_sub_source->transition_type;
  pquad_w_sub_dest->transition_type = pquad_w_sub_source->transition_type;
  pquad_w_sub_dest->subelement_id = pquad_w_sub_source->subelement_id;
}

int
t8_subelement_scheme_quad_c::t8_element_is_subelement (const t8_element *elem) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  T8_ASSERT (pquad_w_sub->transition_type >= 0);

  /* transition_type == 0 => elem is no subelement.
   * transition_type != 0 => elem is subelement 
   */
  return (pquad_w_sub->transition_type == 0 ? false : true);
}

int
t8_subelement_scheme_quad_c::t8_element_get_transition_type (const t8_element *elem)
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  return pquad_w_sub->transition_type;
}

int
t8_subelement_scheme_quad_c::t8_element_get_subelement_id (const t8_element *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  return pquad_w_sub->subelement_id;
}

t8_element_shape_t
t8_subelement_scheme_quad_c::t8_element_shape (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_ECLASS_TRIANGLE : T8_ECLASS_QUAD);
}

int
t8_subelement_scheme_quad_c::t8_element_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_QUAD_SUBELEMENT_FACES : P4EST_FACES);
}

int
t8_subelement_scheme_quad_c::t8_element_find_neighbor_in_transition_cell (const t8_element_t *elem,
                                                                          const t8_element_t *pseudo_neigh,
                                                                          int elem_face)
{
  /* In this function, we assume pseudo_neigh to be a random subelement of a transition cell that includes
   * the real neighbor of elem at face elem_face. This function will output the subelement_id of the real neighbor of elem. */
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (pseudo_neigh));

  /* we expect neigh to be a element in a transition cell, thus to be a subelement */
  T8_ASSERT (t8_element_is_subelement (pseudo_neigh));

  const t8_quad_with_subelements *pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  const t8_quad_with_subelements *pquad_w_sub_pseudo_neigh = (const t8_quad_with_subelements *) pseudo_neigh;

  /* In the following, all possible neighbor configurations are defined, such that subelement neighbors can be
   * identified in LFN_transitioned. */
  if (pquad_w_sub_elem->transition_type != 0 && (elem_face == 0 || elem_face == 2)) {
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x
     *      | \    elem   / |
     *      |   \       /   |
     *      |N_f0 \   / N_f2|
     *      x - - - + - - - x
     *      |     /   \     |
     *      |   /pseudo \   |
     *      | /   neigh   \ |
     *      x - - - - - - - x
     *
     * Elem and face 0 or face 2 is given and a random sibling subelement neigh is given, too. 
     * We are searching for the subelement id of the real neighbor N_f0 or N_f2, depending on the face number. */
    int shift;
    if (elem_face == 0) {
      shift = -1;
    }
    if (elem_face == 2) {
      shift = 1;
    }
    int num_subelements = t8_element_get_number_of_subelements (pquad_w_sub_elem->transition_type);
    return ((pquad_w_sub_elem->subelement_id + shift) + num_subelements)
           % num_subelements; /* the neighbor is directly before or after elem modulo the number of subelements in the transition cell */
  }
  /* Below are the cases in which the neighbor can not be identified as simple as above. 
   * The idea is to fill a location array with the desired properties of the real neighbor. 
   * Together with the type of the transition cell of pseudo_neigh, we can then identify the sub_id of the right neighbor. */

  if (pquad_w_sub_elem->transition_type != 0 && elem_face == 1) {
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / | \           / |
     *      |   \       /   |   \       /   |
     *      |     \   /     |     \   /     |
     *      x - - - x neigh | elem  x       |
     *      |     /   \     |     / | \     |
     *      |   /pseudo \   |   /   |   \   |
     *      | /   neigh   \ | /     |     \ |
     *      x - - - - - - - x - - - x - - - x
     *
     * A subelement elem is given as well as a random subelement pseudo_neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh. 
     * Note that both transition cells can have different levels. */

    /* get the location of elem */
    int location_elem[3] = {}; /* {face, is_split, number of subelement at face} */
    t8_element_get_location_of_subelement (elem, location_elem);

    /* Initialize the location array of the real neighbor. */
    int location_neigh[3] = { -1, -1, -1 };

    /* the pseudo_neigh tranaition cell has a lower level than the elem transition cell */
    if (pquad_w_sub_pseudo_neigh->p4q.level < pquad_w_sub_elem->p4q.level) {
      if (location_elem[0] == 0) { /* left face of transition cell */
        if (pquad_w_sub_pseudo_neigh->p4q.y == pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 2; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* second subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.y != pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 2; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
      }
      if (location_elem[0] == 1) { /* upper face of transition cell */
        if (pquad_w_sub_pseudo_neigh->p4q.x == pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 3; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* first or second subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.x != pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 3; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
      }
      if (location_elem[0] == 2) { /* right face of transition cell */
        if (pquad_w_sub_pseudo_neigh->p4q.y == pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 0; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.y != pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 0; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* first or second subelement at face */
        }
      }
      if (location_elem[0] == 3) { /* lower face of transition cell */
        if (pquad_w_sub_pseudo_neigh->p4q.x == pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 1; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.x != pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 1; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* second subelement at face */
        }
      }
    }
    /* the pseudo_neigh tranaition cell has not a lower level than the elem transition cell */
    if (pquad_w_sub_pseudo_neigh->p4q.level >= pquad_w_sub_elem->p4q.level) {
      if (location_elem[0] == 0) { /* left face of transition cell */
        location_neigh[0] = 2;     /* face */
        location_neigh[1] = 0;     /* not split */
        location_neigh[2] = 0;     /* first (only) subelement at face */
      }
      if (location_elem[0] == 1) { /* upper face of transition cell */
        location_neigh[0] = 3;     /* face */
        location_neigh[1] = 0;     /* not split */
        location_neigh[2] = 0;     /* first (only) subelement at face */
      }
      if (location_elem[0] == 2) { /* right face of transition cell */
        location_neigh[0] = 0;     /* face */
        location_neigh[1] = 0;     /* not split */
        location_neigh[2] = 0;     /* first (only) subelement at face */
      }
      if (location_elem[0] == 3) { /* lower face of transition cell */
        location_neigh[0] = 1;     /* face */
        location_neigh[1] = 0;     /* not split */
        location_neigh[2] = 0;     /* first (only) subelement at face */
      }
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0 && location_neigh[2] >= 0);

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return t8_element_get_id_from_location (t8_element_get_transition_type (pseudo_neigh), location_neigh);
  }
  if (!t8_element_is_subelement (elem)) {
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / |               |
     *      |   \       /   |               |
     *      |     \   /     |               |
     *      x - - - x neigh |     elem      |
     *      |     /   \     |               |
     *      |   /pseudo \   |               |
     *      | /   neigh   \ |               |
     *      x - - - - - - - x - - - - - - - x
     *
     * Subelement elem is given as well as a random subelement neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh.
     * Note that the transition cell of pseudo_neigh and elem can have different levels. */

    /* Initialize the location array of the real neighbor. */
    int location_neigh[3] = { -1, -1, -1 };

    /* the pseudo_neigh tranaition cell has a lower level than elem */
    if (pquad_w_sub_pseudo_neigh->p4q.level < pquad_w_sub_elem->p4q.level) {
      if (elem_face == 0) { /* left face */
        if (pquad_w_sub_pseudo_neigh->p4q.y == pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 2; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* second subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.y != pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 2; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
      }
      if (elem_face == 1) { /* right face */
        if (pquad_w_sub_pseudo_neigh->p4q.y == pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 0; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.y != pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 0; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* first or second subelement at face */
        }
      }
      if (elem_face == 2) { /* lower face */
        if (pquad_w_sub_pseudo_neigh->p4q.x == pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 1; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.x != pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 1; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* second subelement at face */
        }
      }
      if (elem_face == 3) { /* upper face */
        if (pquad_w_sub_pseudo_neigh->p4q.x == pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 3; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 1; /* first or second subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.x != pquad_w_sub_elem->p4q.x) {
          location_neigh[0] = 3; /* face */
          location_neigh[1] = 1; /* split */
          location_neigh[2] = 0; /* first subelement at face */
        }
      }
    }
    /* the pseudo_neigh tranaition cell has the same level as elem 
     * Note that the level of the trnasition cell can not be higher as the level of elem in this case, 
     * since elem would then be a subelement in a transition cell. */
    if (pquad_w_sub_pseudo_neigh->p4q.level == pquad_w_sub_elem->p4q.level) {
      if (elem_face == 0) {    /* left face */
        location_neigh[0] = 2; /* face */
        location_neigh[1] = 0; /* not split */
        location_neigh[2] = 0; /* first (only) subelement at face */
      }
      if (elem_face == 1) {    /* right face */
        location_neigh[0] = 0; /* face */
        location_neigh[1] = 0; /* not split */
        location_neigh[2] = 0; /* first (only) subelement at face */
      }
      if (elem_face == 2) {    /* lower face */
        location_neigh[0] = 1; /* face */
        location_neigh[1] = 0; /* not split */
        location_neigh[2] = 0; /* first (only) subelement at face */
      }
      if (elem_face == 3) {    /* upper face */
        location_neigh[0] = 3; /* face */
        location_neigh[1] = 0; /* not split */
        location_neigh[2] = 0; /* first (only) subelement at face */
      }
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0 && location_neigh[2] >= 0);

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return t8_element_get_id_from_location (t8_element_get_transition_type (pseudo_neigh), location_neigh);
  }
  return -1; /* return negative if no neighbor element could be found */
}

int
t8_subelement_scheme_quad_c::t8_element_get_id_from_location (int type, int location[])
{
  T8_ASSERT (type >= 0 && type <= T8_SUB_QUAD_MAX_TRANSITION_TYPE);

  int sub_id, subelements_count = 0;
  double type_temp = double (type);  // would work for ints but we use libc pow(double, double)
  int binary_type[4] = {};
  int binary_type_clockwise[4] = {};

  /* get the type as a binary array */
  int iface;
  for (iface = 0; iface < P4EST_FACES; iface++) {
    if (type_temp >= pow (2.0, 4 - (iface + 1))) {
      binary_type[iface] = 1;
      type_temp -= pow (2.0, 4 - (iface + 1));
    }
    else {
      binary_type[iface] = 0;
    }
  }

  for (iface = 0; iface < P4EST_FACES;
       iface++) { /* rearrange the binary type to be in clockwise order of the faces, starting with the left face */
    binary_type_clockwise[iface] = binary_type[subelement_location_to_parent_face[iface]];
  }

  /* count the number of elements up to the given location */
  int element_count;
  for (element_count = 0; element_count <= location[0]; element_count++) {
    if (element_count == location[0]) {
      if (location[1] == 0) {
        subelements_count += 1;
      }
      else {
        if (location[2] == 0) {
          subelements_count += 1;
        }
        else {
          subelements_count += 2;
        }
      }
    }
    else {
      subelements_count += binary_type_clockwise[element_count] + 1;
    }
  }

  /* get the sub_id */
  sub_id = subelements_count - 1;

  return sub_id;
}

int
t8_subelement_scheme_quad_c::t8_element_get_face_number_of_hypotenuse (const t8_element_t *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));

  int location[3] = {};
  t8_element_get_location_of_subelement (elem, location);

  int split = location[1];
  int second = location[2];

  if (!split) { /* if the face is not split, then the hypotenuse is always face number one */
    return 1;
  }
  else {
    if (!second) { /* otherwise, the subelment is mirrored, depending on the value of 'second' */
      return 0;
    }
    else {
      return 2;
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory */
  t8_default_scheme_common_c::t8_element_new (length, elem);

  int elem_count;
  for (elem_count = 0; elem_count < length; elem_count++) {
    t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem[elem_count];
    t8_element_init (1, elem[elem_count]);
    /* set dimension of quad to 2 */
    T8_QUAD_SET_TDIM ((p4est_quadrant_t *) &pquad_w_sub->p4q, 2);
  }
}

void
t8_subelement_scheme_quad_c::t8_element_init (int length, t8_element_t *elem) const
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;

  int elem_count;

  for (elem_count = 0; elem_count < length; elem_count++) {
    /* initialize subelement parameters */
    pquad_w_sub[elem_count].transition_type = 0;
    pquad_w_sub[elem_count].subelement_id = 0;

#ifdef T8_ENABLE_DEBUG
    /* In debugging mode we iterate over all length many elements and 
     * set their quad to the level 0 quad with ID 0. */
    p4est_quadrant_t *quad = &pquad_w_sub[elem_count].p4q;
    /* Set all values to 0 */
    for (int i = 0; i < length; i++) {
      p4est_quadrant_set_morton (quad + i, 0, 0);
      T8_QUAD_SET_TDIM (quad + i, 2);
      T8_ASSERT (p4est_quadrant_is_extended (quad + i));
    }
#endif
  }
}

int
t8_subelement_scheme_quad_c::t8_element_scheme_supports_transitioning (void)
{
  return T8_QUAD_TRANSITION_IS_IMPLEMENTED;
}

int
t8_subelement_scheme_quad_c::t8_element_transition_scheme_is_conformal (void)
{
  return T8_QUAD_TRANSITION_SCHEME_IS_CONFORMAL;
}

int
t8_subelement_scheme_quad_c::t8_element_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{

  if (t8_element_get_subelement_id ((const t8_element_t *) elem1)
      == t8_element_get_subelement_id ((const t8_element_t *) elem2)) {
    return 1;
  }
  else {
    return 0;
  }
}

void
t8_subelement_scheme_quad_c::t8_element_root (t8_element_t *elem) const
{
  SC_ABORT_NOT_REACHED ();
}

void
t8_subelement_scheme_quad_c::t8_element_MPI_Pack (t8_element_t **const elements, const unsigned int count,
                                                  void *send_buffer, const int buffer_size, int *position,
                                                  sc_MPI_Comm comm) const
{
  SC_ABORT_NOT_REACHED ();

  // int mpiret;
  // p4est_quadrant_t **quads = (p4est_quadrant_t **) elements;
  // for (unsigned int ielem = 0; ielem < count; ielem++) {
  //   mpiret = sc_MPI_Pack (&(quads[ielem]->x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
  //   SC_CHECK_MPI (mpiret);
  //   mpiret = sc_MPI_Pack (&quads[ielem]->y, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
  //   SC_CHECK_MPI (mpiret);
  //   mpiret = sc_MPI_Pack (&quads[ielem]->level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
  //   SC_CHECK_MPI (mpiret);
  // }
}

/* each quad is packed as x,y coordinates and the level */
void
t8_subelement_scheme_quad_c::t8_element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  SC_ABORT_NOT_REACHED ();
  // int singlesize = 0;
  // int datasize = 0;
  // int mpiret;

  // /* x,y */
  // mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  // SC_CHECK_MPI (mpiret);
  // singlesize += 2 * datasize;

  // /* level */
  // mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  // SC_CHECK_MPI (mpiret);
  // singlesize += datasize;

  // *pack_size = count * singlesize;
}

/* each quad is packed as x,y coordinates and the level */
void
t8_subelement_scheme_quad_c::t8_element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                                    t8_element_t **elements, const unsigned int count,
                                                    sc_MPI_Comm comm) const
{
  SC_ABORT_NOT_REACHED ();
  // int mpiret;
  // p4est_quadrant_t **quads = (p4est_quadrant_t **) elements;
  // for (unsigned int ielem = 0; ielem < count; ielem++) {
  //   mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->x), 1, sc_MPI_INT, comm);
  //   SC_CHECK_MPI (mpiret);
  //   mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->y), 1, sc_MPI_INT, comm);
  //   SC_CHECK_MPI (mpiret);
  //   mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem]->level), 1, sc_MPI_INT8_T, comm);
  //   SC_CHECK_MPI (mpiret);
  // }
}

#ifdef T8_ENABLE_DEBUG
void
t8_subelement_scheme_quad_c::t8_element_debug_print (const t8_element_t *elem) const
{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  t8_productionf ("\n|------------ t8_element_debug_print: ------------|"
                  "\n|    Transition Type:     %i"
                  "\n|    Subelement ID:       %i"
                  "\n|    Anchor (Morton):     (%i,%i)"
                  "\n|    Anchor (ref coords): (%lf,%lf)"
                  "\n|    Level:               %i"
                  "\n|-------------------------------------------------|\n",
                  pquad_w_sub->transition_type, pquad_w_sub->subelement_id, pquad_w_sub->p4q.x, pquad_w_sub->p4q.y,
                  (double) pquad_w_sub->p4q.x / (double) P4EST_ROOT_LEN,
                  (double) pquad_w_sub->p4q.y / (double) P4EST_ROOT_LEN, pquad_w_sub->p4q.level);

  /* if the element is not valid, abort, but after printing */
  T8_ASSERT (t8_element_is_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_quad_c::t8_element_is_valid (const t8_element_t *elem) const
/* *INDENT-ON* */

{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* the 4pest quadrant AND the subelement values must be valid such that the whole element is valid */
  return (p4est_quadrant_is_extended (q) && t8_element_subelement_values_are_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_quad_c::t8_element_subelement_values_are_valid (const t8_element_t *elem) const
/* *INDENT-ON* */

{
  const t8_quad_with_subelements *pquad_w_sub = (const t8_quad_with_subelements *) elem;

  return ((pquad_w_sub->transition_type >= 0 && pquad_w_sub->transition_type <= T8_SUB_QUAD_MAX_TRANSITION_TYPE)
          || pquad_w_sub->transition_type == 0)
         && ((pquad_w_sub->subelement_id >= 0 && pquad_w_sub->subelement_id <= T8_SUB_QUAD_MAX_SUBELEMENT_ID)
             || pquad_w_sub->subelement_id == 0);
}

void
t8_subelement_scheme_quad_c::t8_element_to_string (const t8_element_t *elem, char *debug_string,
                                                   const int string_size) const
{
  SC_ABORT_NOT_REACHED ();
}
#endif

/* Constructor */
t8_subelement_scheme_quad_c::t8_subelement_scheme_quad_c (void)
{
  eclass = T8_ECLASS_QUAD;
  element_size = sizeof (t8_quad_with_subelements);
  ts_context = sc_mempool_new (element_size);
}

t8_subelement_scheme_quad_c::~t8_subelement_scheme_quad_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
