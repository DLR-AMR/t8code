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

#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_common_cxx.hxx>

#include "t8_subelements_quad_cxx.hxx"

/* *INDENT-OFF* */
const int           subelement_face_dual[3] = { 2, 1, 0};
/* *INDENT-ON* */

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This function is used by other element functions and we thus need to
 * declare it up here */
t8_linearidx_t      t8_element_get_linear_id (const t8_element_t * elem,
                                              int level);

int
t8_subelement_scheme_quad_c::t8_element_maxlevel (void)
{
  return P4EST_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_subelement_scheme_quad_c::t8_element_child_eclass (int childid)
/* *INDENT-ON* */

{
  T8_ASSERT (0 <= childid && childid < P4EST_CHILDREN);

  return T8_ECLASS_QUAD;
}

int
t8_subelement_scheme_quad_c::t8_element_level (const t8_element_t * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_quad_with_subelements *) pquad_w_sub)->p4q.level;
}

static void
t8_element_copy_surround (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  T8_QUAD_SET_TDIM (r, T8_QUAD_GET_TDIM (q));
  if (T8_QUAD_GET_TDIM (q) == 3) {
    T8_QUAD_SET_TNORMAL (r, T8_QUAD_GET_TNORMAL (q));
    T8_QUAD_SET_TCOORD (r, T8_QUAD_GET_TCOORD (q));
  }
}

void
t8_subelement_scheme_quad_c::t8_element_copy (const t8_element_t * source,
                                              t8_element_t * dest)
{
  const t8_quad_with_subelements *pquad_w_sub_source =
    (const t8_quad_with_subelements *) source;
  t8_quad_with_subelements *pquad_w_sub_dest =
    (t8_quad_with_subelements *) dest;

  const p4est_quadrant_t *q = &pquad_w_sub_source->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_dest->p4q;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (q == r &&
      pquad_w_sub_source->dummy_is_subelement ==
      pquad_w_sub_dest->dummy_is_subelement
      && pquad_w_sub_source->subelement_type ==
      pquad_w_sub_dest->subelement_type
      && pquad_w_sub_source->subelement_id ==
      pquad_w_sub_dest->subelement_id) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *r = *q;

  t8_element_copy_subelement_values (source, dest);
  t8_element_copy_surround (q, r);
}

int
t8_subelement_scheme_quad_c::t8_element_compare (const t8_element_t * elem1,
                                                 const t8_element_t * elem2)
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 =
    (const t8_quad_with_subelements *) elem1;
  const t8_quad_with_subelements *pquad_w_sub_elem2 =
    (const t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *q = &pquad_w_sub_elem1->p4q;
  const p4est_quadrant_t *r = &pquad_w_sub_elem2->p4q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  if (p4est_quadrant_compare (q, r) == 0) {
    if (t8_element_test_if_subelement (elem1) == 1
        && t8_element_test_if_subelement (elem2) == 1) {
      SC_ABORT
        ("Both elements are subelements. Specify what the compare function should do in this case");
    }
    else if (t8_element_test_if_subelement (elem1) == 1) {
      return -1;                /* elem1 is subelement and therefore smaller */
    }
    else if (t8_element_test_if_subelement (elem2) == 1) {
      return 1;                 /* elem2 is subelement and therefore smaller */
    }
  }

  /* Note that for subelements, their parent quadrant is compared at this point */
  return p4est_quadrant_compare (q, r);
}

void
t8_subelement_scheme_quad_c::t8_element_parent (const t8_element_t * elem,
                                                t8_element_t * parent)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_parent =
    (t8_quad_with_subelements *) parent;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_parent->p4q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));

  if (pquad_w_sub_elem->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    pquad_w_sub_parent->p4q = pquad_w_sub_elem->p4q;
  }
  else {
    p4est_quadrant_parent (q, r);
  }

  /* resetting the subelement infos and ensuring that a parent element can never be a subelement */
  t8_element_reset_subelement_values (parent);

  t8_element_copy_surround (q, r);
}

void
t8_subelement_scheme_quad_c::t8_element_sibling (const t8_element_t * elem,
                                                 int sibid,
                                                 t8_element_t * sibling)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_sibling =
    (t8_quad_with_subelements *) sibling;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_sibling->p4q;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub_elem->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));

  p4est_quadrant_sibling (q, r, sibid);
  t8_element_copy_surround (q, r);
}

int
t8_subelement_scheme_quad_c::t8_element_num_faces (const t8_element_t * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  T8_ASSERT (t8_element_is_valid (elem));

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    return T8_SUBELEMENT_FACES;
  }
  else {
    return P4EST_FACES;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_max_num_faces (const t8_element_t *
                                                       elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  return P4EST_FACES;
}

int
t8_subelement_scheme_quad_c::t8_element_num_children (const t8_element_t *
                                                      elem)
{
  /* Note that children of subelements equal the children of the parent quadrant. 
   * Therefore, the number of children of a subelement equals P4EST_CHILDREN */
  T8_ASSERT (t8_element_is_valid (elem));
  return P4EST_CHILDREN;
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int   
t8_subelement_scheme_quad_c::t8_element_num_siblings (const t8_element_t *
                                                      elem) const
/* *INDENT-ON* */

{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* TODO: use t8_element_get_number_of_subelements instead (problem with const) */
  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    int                 type = pquad_w_sub->subelement_type;
    int                 num_hanging_faces = 0;
    int                 num_siblings;
    int                 i;

    for (i = 0; i < P4EST_FACES; i++) { /* Count the number of ones of the binary subelement type */
      num_hanging_faces += (type & (1 << i)) >> i;
    }
    /* The number of subelements equals the number of neighbours: */
    num_siblings = P4EST_FACES + num_hanging_faces;
    return num_siblings;
  }
  else {
    return P4EST_CHILDREN;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_num_face_children (const t8_element_t
                                                           * elem, int face)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  return 2;
}

int
t8_subelement_scheme_quad_c::t8_element_get_face_corner (const t8_element_t *
                                                         element, int face,
                                                         int corner)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) element;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);
  /* TODO: check whether this enumeration of the faces is right. It might be f_3 and f_2 switched */
  /*
   *   2    f_2    3
   *     x -->-- x
   *     |       |
   *     ^       ^
   * f_0 |       | f_1
   *     x -->-- x
   *   0    f_3    1
   */

  T8_ASSERT (t8_element_is_valid (element));
  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (0 <= corner && corner < 2);
  return p4est_face_corners[face][corner];
}

int
t8_subelement_scheme_quad_c::t8_element_get_corner_face (const t8_element_t *
                                                         element, int corner,
                                                         int face)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) element;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (element));
  T8_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  T8_ASSERT (0 <= face && face < 2);
  return p4est_corner_faces[corner][face];
}

void
t8_subelement_scheme_quad_c::t8_element_child (const t8_element_t * elem,
                                               int childid,
                                               t8_element_t * child)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_child =
    (t8_quad_with_subelements *) child;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_child->p4q;

  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);

  /* it should not be possible to construct a child of a subelement */
  T8_ASSERT (pquad_w_sub_elem->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

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

void
t8_subelement_scheme_quad_c::t8_element_children (const t8_element_t * elem,
                                                  int length,
                                                  t8_element_t * c[])
{
  /* if elem is a subelement, then this function will construct the children of its parent p4est quadrant */
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements **pquad_w_sub_children =
    (t8_quad_with_subelements **) c;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;

  int                 i;

  T8_ASSERT (t8_element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < P4EST_CHILDREN; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif
  T8_ASSERT (length == P4EST_CHILDREN);

  /* set coordinates and levels of the children */
  p4est_quadrant_children (q, &pquad_w_sub_children[0]->p4q,
                           &pquad_w_sub_children[1]->p4q,
                           &pquad_w_sub_children[2]->p4q,
                           &pquad_w_sub_children[3]->p4q);

  for (i = 0; i < P4EST_CHILDREN; ++i) {
    t8_element_reset_subelement_values (c[i]);
    t8_element_copy_surround (q, &pquad_w_sub_children[i]->p4q);
  }
}

int
t8_subelement_scheme_quad_c::t8_element_child_id (const t8_element_t * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  T8_ASSERT (t8_element_is_valid (elem));

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    /* the child_id of a subelement equals its subelement_id */
    return pquad_w_sub->subelement_id;
  }
  else {
    /* the child_id of a non-subelement element can be determined by its level and anchor node */
    return p4est_quadrant_child_id (q);
  }
}

int
t8_subelement_scheme_quad_c::t8_element_ancestor_id (const t8_element_t *
                                                     elem, int level)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  return p4est_quadrant_ancestor_id (q, level);
}

int
t8_subelement_scheme_quad_c::t8_element_is_family (t8_element_t ** fam)
{
  t8_quad_with_subelements **pquad_w_sub_family =
    (t8_quad_with_subelements **) fam;

#ifdef T8_ENABLE_DEBUG
  int                 i;

  /* TODO: this loop goes from 0 to 3 but with subelements there can be more elements in fam */
  for (i = 0; i < P4EST_CHILDREN; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif

  /* Subelements can not be refined into other elements of a higher level. 
   * So, if the first element of fam is a subelement, we assume that the following num_siblings 
   * many elements are its siblings and therefore form a family. */
  if (pquad_w_sub_family[0]->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    return 1;
  }
  /* If the first element of fam is no subelement we check the following elements of fam */
  else {
    /* If any of the following elements is a subelement, then they can not form a family */
    if (pquad_w_sub_family[1]->dummy_is_subelement ==
        T8_SUB_QUAD_IS_SUBELEMENT
        || pquad_w_sub_family[2]->dummy_is_subelement ==
        T8_SUB_QUAD_IS_SUBELEMENT
        || pquad_w_sub_family[3]->dummy_is_subelement ==
        T8_SUB_QUAD_IS_SUBELEMENT) {
      return 0;
    }
    /* If all elements of fam are no subelements, then we can use the p4est check is_family */
    else {
      return p4est_quadrant_is_family (&pquad_w_sub_family[0]->p4q,
                                       &pquad_w_sub_family[1]->p4q,
                                       &pquad_w_sub_family[2]->p4q,
                                       &pquad_w_sub_family[3]->p4q);
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_set_linear_id (t8_element_t * elem,
                                                       int level,
                                                       t8_linearidx_t id)
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t   *q = &pquad_w_sub->p4q;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << P4EST_DIM * level);

  p4est_quadrant_set_morton (q, level, id);
  T8_QUAD_SET_TDIM (q, 2);
}

t8_linearidx_t
  t8_subelement_scheme_quad_c::t8_element_get_linear_id (const t8_element_t *
                                                         elem, int level)
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t   *q = &pquad_w_sub->p4q;

  /* Note that the id of a subelement equals the id of its parent quadrant.
   * Therefore, the binary search (for example used in the leaf_face_neighbor function) 
   * will find a random subelement of the parent transition cell which might not be the desired neighbor of a given element. */

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  return p4est_quadrant_linear_id (q, level);
}

void
t8_subelement_scheme_quad_c::t8_element_first_descendant (const t8_element_t *
                                                          elem,
                                                          t8_element_t * desc,
                                                          int level)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_desc =
    (t8_quad_with_subelements *) desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_desc->p4q;

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
t8_subelement_scheme_quad_c::t8_element_last_descendant (const t8_element_t *
                                                         elem,
                                                         t8_element_t * desc,
                                                         int level)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_desc =
    (t8_quad_with_subelements *) desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_desc->p4q;

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
t8_subelement_scheme_quad_c::t8_element_successor (const t8_element_t * elem1,
                                                   t8_element_t * elem2,
                                                   int level)
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 =
    (const t8_quad_with_subelements *) elem1;
  t8_quad_with_subelements *pquad_w_sub_elem2 =
    (t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *q = &pquad_w_sub_elem1->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_elem2->p4q;

  t8_linearidx_t      id;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub_elem1->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  id = p4est_quadrant_linear_id (q, level);
  T8_ASSERT (id + 1 < ((t8_linearidx_t) 1) << P4EST_DIM * level);
  t8_element_reset_subelement_values (elem2);
  p4est_quadrant_set_morton (r, level, id + 1);
  t8_element_copy_surround (q, r);
}

void
t8_subelement_scheme_quad_c::t8_element_nca (const t8_element_t * elem1,
                                             const t8_element_t * elem2,
                                             t8_element_t * nca)
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 =
    (const t8_quad_with_subelements *) elem1;
  const t8_quad_with_subelements *pquad_w_sub_elem2 =
    (const t8_quad_with_subelements *) elem2;
  t8_quad_with_subelements *pquad_w_sub_nca =
    (t8_quad_with_subelements *) nca;

  const p4est_quadrant_t *q1 = &pquad_w_sub_elem1->p4q;
  const p4est_quadrant_t *q2 = &pquad_w_sub_elem2->p4q;
  p4est_quadrant_t   *r = &pquad_w_sub_nca->p4q;

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
  t8_subelement_scheme_quad_c::t8_element_face_shape (const t8_element_t *
                                                      elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_ECLASS_LINE;
}

void
t8_subelement_scheme_quad_c::t8_element_children_at_face (const t8_element_t *
                                                          elem, int face,
                                                          t8_element_t *
                                                          children[],
                                                          int num_children,
                                                          int *child_indices)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  int                 first_child, second_child;

#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    for (i = 0; i < num_children; i++) {
      T8_ASSERT (t8_element_is_valid (children[i]));
    }
  }
#endif
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
   *      | 2   | 3   |           For examle for face 2 (bottom face) we see
   * 0    x - - x - - x   1       first_child = 0 and second_child = 1.
   *      |     |     |
   *      | 0   | 1   |
   *      x - - x - - x
   *
   *            2
   */
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
t8_subelement_scheme_quad_c::t8_element_face_child_face (const t8_element_t *
                                                         elem, int face,
                                                         int face_child)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  /* For quadrants the face enumeration of children is the same as for the parent. */
  return face;
}

int
t8_subelement_scheme_quad_c::t8_element_face_parent_face (const t8_element_t *
                                                          elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  int                 child_id;

  /* for subelements qe need to adjust the output of this function.
   * A subelements face is a subface of the parent quadrant if and only if the face number is 1. */
  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    if (face == 1) {
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      if (location[0] == 0) {
        return 0;
      }
      else if (location[0] == 1) {
        return 3;
      }
      else if (location[0] == 2) {
        return 1;
      }
      else if (location[0] == 3) {
        return 2;
      }
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
  if (child_id == p4est_face_corners[face][0]
      || child_id == p4est_face_corners[face][1]) {
    return face;
  }
  return -1;
}

void
t8_subelement_scheme_quad_c::t8_element_transform_face (const t8_element_t *
                                                        elem1,
                                                        t8_element_t * elem2,
                                                        int orientation,
                                                        int sign,
                                                        int is_smaller_face)
{
  const t8_quad_with_subelements *pquad_w_sub_elem1 =
    (const t8_quad_with_subelements *) elem1;
  t8_quad_with_subelements *pquad_w_sub_elem2 =
    (t8_quad_with_subelements *) elem2;

  const p4est_quadrant_t *qin = &pquad_w_sub_elem1->p4q;
  p4est_quadrant_t   *p = &pquad_w_sub_elem2->p4q;

  const p4est_quadrant_t *q;
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (qin->level);
  p4est_qcoord_t      x = qin->x;       /* temp storage for x coordinate in case elem1 = elem 2 */

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub_elem1->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

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
    x = q->x;                   /* temp storage in case elem1 = elem 2 */
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
  case 0:                      /* Nothing to do */
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
t8_subelement_scheme_quad_c::t8_element_extrude_face (const t8_element_t *
                                                      face,
                                                      const t8_eclass_scheme_c
                                                      * face_scheme,
                                                      t8_element_t * elem,
                                                      int root_face)
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t   *q = &pquad_w_sub->p4q;

  const t8_dline_t   *l = (const t8_dline_t *) face;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (face_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (face_scheme->t8_element_is_valid (elem));
  T8_ASSERT (0 <= root_face && root_face < P4EST_FACES);

  /* TODO: check if this enumeration is right (maybe f_2 and f_3 switched) */
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
t8_subelement_scheme_quad_c::t8_element_tree_face (const t8_element_t * elem,
                                                   int face)
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    T8_ASSERT (face != 1);      /* this function does only make sense for subelements at face 1 */
    T8_ASSERT (0 <= face && face < T8_SUBELEMENT_FACES);

    int                 location[3] = { };
    t8_element_get_location_of_subelement (elem, location);

    /* subelements are enumerated clockwise (not as quadrant faces) */
    if (location[1] == 0) {
      return 0;
    }
    else if (location[0] == 1) {
      return 3;
    }
    else if (location[0] == 2) {
      return 1;
    }
    else {
      return 2;
    }
  }
  else {
    T8_ASSERT (0 <= face && face < P4EST_FACES);
    /* For quadrants the face and the tree face number are the same. */
    return face;
  }
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_subelement_scheme_quad_c::t8_element_first_descendant_face (const
                                                               t8_element_t *
                                                               elem, int face,
                                                               t8_element_t *
                                                               first_desc,
                                                               int level)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_first_desc =
    (t8_quad_with_subelements *) first_desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *desc = &pquad_w_sub_first_desc->p4q;

  int                 first_face_corner;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub_elem->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

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
t8_subelement_scheme_quad_c::t8_element_last_descendant_face (const
                                                              t8_element_t *
                                                              elem, int face,
                                                              t8_element_t *
                                                              last_desc,
                                                              int level)
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_last_desc =
    (t8_quad_with_subelements *) last_desc;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *desc = &pquad_w_sub_last_desc->p4q;

  int                 last_face_corner;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub_elem->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);
  T8_ASSERT (pquad_w_sub_last_desc->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (0 <= face && face < P4EST_FACES);
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  /* Get the last corner of q that belongs to face */
  last_face_corner = p4est_face_corners[face][1];
  /* Construce the descendant in that corner */
  p4est_quadrant_corner_descendant (q, desc, last_face_corner, level);
  t8_element_reset_subelement_values (last_desc);
}

void
t8_subelement_scheme_quad_c::t8_element_boundary_face (const t8_element_t *
                                                       elem, int face,
                                                       t8_element_t *
                                                       boundary,
                                                       const
                                                       t8_eclass_scheme_c *
                                                       boundary_scheme)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  t8_dline_t         *l = (t8_dline_t *) boundary;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (boundary_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));
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
  l->x = ((face >> 1 ? q->x : q->y) *
          ((int64_t) T8_DLINE_ROOT_LEN) / P4EST_ROOT_LEN);
}

void
t8_subelement_scheme_quad_c::t8_element_boundary (const t8_element_t * elem,
                                                  int min_dim, int length,
                                                  t8_element_t ** boundary)
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
t8_subelement_scheme_quad_c::t8_element_is_root_boundary (const t8_element_t *
                                                          elem, int face)
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  p4est_qcoord_t      coord;

  /* in case of a subelement we might adjust its face number with its parents face number */
  int                 adjusted_face_in_case_of_subelement = face;

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    if (face == 1) {
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      if (location[0] == 0) {
        adjusted_face_in_case_of_subelement = 0;
      }
      else if (location[0] == 1) {
        adjusted_face_in_case_of_subelement = 3;
      }
      else if (location[0] == 2) {
        adjusted_face_in_case_of_subelement = 1;
      }
      else if (location[0] == 3) {
        adjusted_face_in_case_of_subelement = 2;
      }
    }
    else {                      /* in case of a subelement and face 0 or 2 the face is no subface of the root boundary */
      return false;
    }
  }

  T8_ASSERT (0 <= adjusted_face_in_case_of_subelement
             && adjusted_face_in_case_of_subelement < P4EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   */
  coord = adjusted_face_in_case_of_subelement >> 1 ? q->y : q->x;
  /* If face is 0 or 2 check against 0.
   * If face is 1 or 3  check against LAST_OFFSET */
  return coord ==
    (adjusted_face_in_case_of_subelement & 1 ? P4EST_LAST_OFFSET (q->level) :
     0);
}

int
t8_subelement_scheme_quad_c::t8_element_face_neighbor_inside (const
                                                              t8_element_t *
                                                              elem,
                                                              t8_element_t *
                                                              neigh, int face,
                                                              int *neigh_face)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < P4EST_FACES);

  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements *pquad_w_sub_neigh =
    (t8_quad_with_subelements *) neigh;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;
  p4est_quadrant_t   *n = &pquad_w_sub_neigh->p4q;

  /* In case of a subelement one should construct the face neighbor of the face-corresponding child quadrant
   * of the subelements parent quadrant. Therefore we need to increase the subelements level by one and adapt its
   * anchor node to its specific child_id. */
  if (t8_element_test_if_subelement (elem) == T8_SUB_QUAD_IS_SUBELEMENT) {      /* if elem is a subelement */

    T8_ASSERT (0 <= face && face <= T8_SUBELEMENT_FACES);

    if (face == 0) {
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 2) {
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 1) {
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      /* setting the anchor node of the neighbor element */
      n->x = q->x;
      n->y = q->y;

      /* length of a subelements cathete */
      const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (q->level + 1);

      /* We need to take into account whether the subelement is split or not */
      if (location[1] == 1) {   /* split */

        /* adjust the level of the neighbor of the element */
        n->level = q->level + 1;

        /* adjust the anchor node of the neighbor of the subelement depending on its location at the parent quad */
        if (location[0] == 0) { /* left face */
          if (location[2] == 0) {
            n->x = q->x - shift;
          }
          else {
            n->x = q->x - shift;
            n->y = q->y + shift;
          }
        }
        else if (location[0] == 2) {    /* right face */
          if (location[2] == 0) {
            n->x = q->x + 2 * shift;
            n->y = q->y + shift;
          }
          else {
            n->x = q->x + 2 * shift;
          }
        }
        else if (location[0] == 3) {    /* lower face */
          if (location[2] == 0) {
            n->x = q->x + shift;
            n->y = q->y - shift;
          }
          else {
            n->y = q->y - shift;
          }
        }
        else {                  /* upper face */
          if (location[2] == 0) {
            n->y = q->y + shift;
          }
          else {
            n->x = q->x + shift;
            n->y = q->y + shift;
          }
        }
      }

      else {                    /* not split */
        /* adjust the level of the neighbor of the subelement */
        n->level = q->level;

        /* adjust the anchor node of the neighbor of the subelement depending on its location at the parent quad */
        if (location[0] == 0) { /* left face */
          n->x = q->x - 2 * shift;
        }
        else if (location[0] == 2) {    /* right face */
          n->x = q->x + 2 * shift;
        }
        else if (location[0] == 3) {    /* lower face */
          n->y = q->y - 2 * shift;
        }
        else {                  /* upper face */
          n->y = q->y + 2 * shift;
        }
      }
    }
  }

  else {                        /* if elem is no subelement */
    /* Directly construct the face neighbor */
    p4est_quadrant_face_neighbor (q, face, n);
  }

  t8_element_reset_subelement_values (neigh);

  T8_QUAD_SET_TDIM (n, 2);

  /* In the following we set the dual faces of our element at the given face. 
   * This does only make sense if the neighbor element at face is of the same type as the 
   * current element (either subelement or non subelement) */
  if (t8_element_test_if_subelement (elem) == T8_SUB_QUAD_IS_SUBELEMENT) {
    /* Compute the face number as seen from q.
     *  0 -> 2    2 -> 0    1 -> 1
     */
    *neigh_face = subelement_face_dual[face];
  }
  else {
    /* Compute the face number as seen from q.
     *  0 -> 1    2 -> 3
     *  1 -> 0    3 -> 2
     */
    T8_ASSERT (neigh_face != NULL);
    *neigh_face = p4est_face_dual[face];
  }

  /* return true if neigh is inside the root */
  return p4est_quadrant_is_inside_root (n);
}

void
t8_subelement_scheme_quad_c::t8_element_anchor (const t8_element_t * elem,
                                                int coord[3])
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;
  p4est_quadrant_t   *q = &pquad_w_sub->p4q;

  /* at the moment, this function is only implemented for standard quad elements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));

  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
  T8_QUAD_SET_TDIM (q, 2);
}

int
t8_subelement_scheme_quad_c::t8_element_root_len (const t8_element_t * elem)
{
  return P4EST_ROOT_LEN;
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_coords (const t8_element_t * t,
                                                       int vertex,
                                                       int coords[])
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) t;
  const p4est_quadrant_t *q1 = &pquad_w_sub->p4q;

  T8_ASSERT (t8_element_is_valid (t));

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_NO_SUBELEMENT) {
    int                 len;

    T8_ASSERT (0 <= vertex && vertex < 4);
    /* Get the length of the quadrant */
    len = P4EST_QUADRANT_LEN (q1->level);

    /* Compute the x and y coordinates of the vertex depending on the
     * vertex number */
    coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
    coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
  }
  else {
    t8_element_vertex_coords_of_subelement (t, vertex, coords);
  }
}

void
t8_subelement_scheme_quad_c::t8_element_vertex_coords_of_subelement (const
                                                                     t8_element_t
                                                                     * t,
                                                                     int
                                                                     vertex,
                                                                     int
                                                                     coords[])
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) t;
  const p4est_quadrant_t *q1 = &pquad_w_sub->p4q;

  int                 len;

  T8_ASSERT (t8_element_is_valid (t));
  T8_ASSERT (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT);
  // T8_ASSERT (vertex >= 0 && vertex < 3); /* all subelements are triangles */

  /* get the length of the current quadrant */
  len = P4EST_QUADRANT_LEN (q1->level);

  /* Compute the x and y coordinates of subelement vertices, depending on the subelement type, id and vertex number (faces enumerated clockwise): 
   *
   *               f1                      V0
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - x - - x f2  -->   x - - x 
   *         | 0 / | \ 4 |           V2    V1
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f3
   * 
   * In this example, the below location array would contain the values [2, 1, 1] 
   * (second face, split, first subelement at the second face) */

  /* get location information of the given subelement */
  int                 location[3] = { };
  t8_element_get_location_of_subelement (t, location);

  /* the face number, the subelement is adjacent to */
  int                 face_number = location[0];
  /* = 1, if the adjacent face is split and = 0, if not */
  int                 split = location[1];
  /* = 0, if the subelement is the first (of two) subelements, at the adjacent face and = 1 if it is the second */
  int                 sub_face_id = location[2];

  /* Check, whether the get_location function provides meaningful location data */
  T8_ASSERT (face_number == 0 || face_number == 1 || face_number == 2
             || face_number == 3);
  T8_ASSERT ((split == 0 && sub_face_id == 0)
             || (split == 1 && (sub_face_id == 0 || sub_face_id == 1)));

  /* using the location data to determine vertex coordinates */
  if (vertex == 2) {            /* the third vertex allways equals the center of the element */
    coords[0] = q1->x + len / 2;
    coords[1] = q1->y + len / 2;
  }
  else {                        /* all other verticies can be determined, using the flag parameters sub_face_id, split and vertex, whose values are either 0 or 1 */
    if (face_number == 0) {
      coords[0] = q1->x;
      coords[1] =
        q1->y + len / 2 * sub_face_id + len * vertex -
        len / 2 * split * vertex;
    }
    else if (face_number == 1) {
      coords[0] =
        q1->x + len / 2 * sub_face_id + len * vertex -
        len / 2 * split * vertex;
      coords[1] = q1->y + len;
    }
    else if (face_number == 2) {
      coords[0] = q1->x + len;
      coords[1] =
        q1->y + len - len / 2 * sub_face_id - (len * vertex -
                                               len / 2 * split * vertex);
    }
    else if (face_number == 3) {
      coords[0] =
        q1->x + len - len / 2 * sub_face_id - (len * vertex -
                                               len / 2 * split * vertex);
      coords[1] = q1->y;
    }
  }
}

void
t8_subelement_scheme_quad_c::t8_element_to_subelement (const t8_element_t *
                                                       elem, int type,
                                                       t8_element_t * c[])
{
  const t8_quad_with_subelements *pquad_w_sub_elem =
    (const t8_quad_with_subelements *) elem;
  t8_quad_with_subelements **pquad_w_sub_subelement =
    (t8_quad_with_subelements **) c;

  const p4est_quadrant_t *q = &pquad_w_sub_elem->p4q;

  int                 num_subelements =
    t8_element_get_number_of_subelements (type, elem);

  T8_ASSERT (type >= T8_SUB_QUAD_MIN_SUBELEMENT_TYPE
             && type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE);

  T8_ASSERT (pquad_w_sub_elem->dummy_is_subelement ==
             T8_SUB_QUAD_IS_NO_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));
#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < num_subelements; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif

  /* get the length of a children-quadrant */
  const int8_t        level = (int8_t) (q->level);

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

  int                 sub_id_counter = 0;
  for (sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
    pquad_w_sub_subelement[sub_id_counter]->p4q.x = q->x;
    pquad_w_sub_subelement[sub_id_counter]->p4q.y = q->y;
    pquad_w_sub_subelement[sub_id_counter]->p4q.level = level;
    pquad_w_sub_subelement[sub_id_counter]->dummy_is_subelement =
      T8_SUB_QUAD_IS_SUBELEMENT;
    pquad_w_sub_subelement[sub_id_counter]->subelement_type = type;
    pquad_w_sub_subelement[sub_id_counter]->subelement_id = sub_id_counter;
  }

  int                 i;
  for (i = 0; i < num_subelements; ++i) {
    T8_ASSERT (t8_element_is_valid (c[i]));
    t8_element_copy_surround (q, &pquad_w_sub_subelement[i]->p4q);
  }
}

int
t8_subelement_scheme_quad_c::t8_element_get_number_of_subelements (int
                                                                   subelement_type,
                                                                   const
                                                                   t8_element
                                                                   * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* consider subelement_type 13 = 1101 in base two -> there are 4 + (1+1+0+1) = 7 subelements */
  int                 num_subelements;
  int                 num_hanging_faces = 0;

  int                 i;

  for (i = 0; i < P4EST_FACES; i++) {   /* Count the number of ones of the binary subelement type. This number equals the number of hanging faces. */
    num_hanging_faces += (subelement_type & (1 << i)) >> i;
  }

  /* The number of subelements equals the number of neighbours: */
  num_subelements = P4EST_FACES + num_hanging_faces;
  return num_subelements;
}

void
t8_subelement_scheme_quad_c::t8_element_get_location_of_subelement (const
                                                                    t8_element_t
                                                                    * elem,
                                                                    int
                                                                    location
                                                                    [])
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  /* this function only works for subelements */
  T8_ASSERT (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT);

  T8_ASSERT (t8_element_is_valid (elem));

  /* Consider the following subelement of type 13:
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
   * These informations are then saved in the location array which will be used by the element_vertex function, 
   * to automatically determine the vertex coordinates of the given subelement. 
   * 
   * The location array for the above example would be {1,1,1} (upper face, split = true, second subelement at the upper face). */

  /* 1) convert the subelement type from a decimal to a binary representation */
  int                 type = pquad_w_sub->subelement_type;
  int                 binary_array[P4EST_FACES] = { };

  int                 i;

  for (i = 0; i < P4EST_FACES; i++) {   /* need an array with 4 elements to store all subelement types of the quad scheme from 1 to 15 ({0,0,0,1} to {1,1,1,1}) */
    binary_array[(P4EST_FACES - 1) - i] = (type & (1 << i)) >> i;
  }                             /* we now got a binary represenation of the subelement type, bitwise stored in an array */

  /* 2) rearrange the binary representation to be in clockwise order */
  int                 binary_array_temp[P4EST_FACES] = { };

  int                 j;

  for (j = 0; j < P4EST_FACES; j++) {   /* copying the binary array */
    binary_array_temp[j] = binary_array[j];
  }

  binary_array[0] = binary_array_temp[0];       /* f3 <- f3 */
  binary_array[1] = binary_array_temp[3];       /* f2 <- f0 */
  binary_array[2] = binary_array_temp[1];       /* f1 <- f2 */
  binary_array[3] = binary_array_temp[2];       /* f0 <- f1 */

  /* 3) use the rearranged binary representation, and the sub_id to determine the location of the subelement and store these information in an array */
  /*     3.1) location[0] -> the face_number, the subelement is adjacent to */
  /*     3.2) location[1] -> if the face is split or not */
  /*     3.3) location[2] -> if the subelement is the first or second subelement of the face (always the first, if the face is not split) */
  int                 num_subelements =
    t8_element_get_number_of_subelements (pquad_w_sub->subelement_type, elem);
  T8_ASSERT (pquad_w_sub->subelement_id < num_subelements);

  int                 sub_id = pquad_w_sub->subelement_id;
  int                 sub_face_id;
  int                 face_number;
  int                 split;

  int                 k;

  int                 cum_neigh_array[P4EST_FACES] = { };

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
    for (k = 0; k < P4EST_FACES - 1; ++k) {
      if (sub_id >= cum_neigh_array[k] && sub_id < cum_neigh_array[k + 1]) {
        face_number = k + 1;
        break;
      }
    }
  }

  /* 3.2) determine, whether the face is split or not */
  if (binary_array[face_number] == 0) {
    split = 0;                  /* the face is not split */
  }
  else {
    split = 1;                  /* the face is split */
  }

  /* 3.3) determine, whether the subelement is the first or the second subelement at the face */
  if (sub_id + 1 == cum_neigh_array[face_number] && split == 1) {
    sub_face_id = 1;            /* second subelement */
  }
  else {
    sub_face_id = 0;            /* first subelement */
  }

  location[0] = face_number;
  location[1] = split;
  location[2] = sub_face_id;
}

void
t8_subelement_scheme_quad_c::t8_element_get_element_data (const t8_element_t *
                                                          elem,
                                                          int anchor_node[],
                                                          int level[],
                                                          int
                                                          subelement_data[])
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  anchor_node[0] = pquad_w_sub->p4q.x;
  anchor_node[1] = pquad_w_sub->p4q.y;

  level[0] = pquad_w_sub->p4q.level;

  subelement_data[0] = pquad_w_sub->dummy_is_subelement;
  subelement_data[1] = pquad_w_sub->subelement_type;
  subelement_data[2] = pquad_w_sub->subelement_id;
}

void
t8_subelement_scheme_quad_c::t8_element_reset_subelement_values (t8_element *
                                                                 elem)
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;

  pquad_w_sub->dummy_is_subelement = T8_SUB_QUAD_IS_NO_SUBELEMENT;
  pquad_w_sub->subelement_type = T8_SUB_QUAD_IS_NO_SUBELEMENT;
  pquad_w_sub->subelement_id = T8_SUB_QUAD_IS_NO_SUBELEMENT;
}

void
t8_subelement_scheme_quad_c::t8_element_copy_subelement_values (const
                                                                t8_element *
                                                                source,
                                                                t8_element *
                                                                dest)
{
  const t8_quad_with_subelements *pquad_w_sub_source =
    (const t8_quad_with_subelements *) source;
  t8_quad_with_subelements *pquad_w_sub_dest =
    (t8_quad_with_subelements *) dest;

  pquad_w_sub_dest->dummy_is_subelement =
    pquad_w_sub_source->dummy_is_subelement;
  pquad_w_sub_dest->subelement_type = pquad_w_sub_source->subelement_type;
  pquad_w_sub_dest->subelement_id = pquad_w_sub_source->subelement_id;
}

int
t8_subelement_scheme_quad_c::t8_element_test_if_subelement (const
                                                            t8_element * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) {
    return T8_SUB_QUAD_IS_SUBELEMENT;
  }
  else {
    return T8_SUB_QUAD_IS_NO_SUBELEMENT;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_get_subelement_type (const
                                                             t8_element *
                                                             elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_NO_SUBELEMENT) {
    return 0;
  }
  else {
    return pquad_w_sub->subelement_type;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_get_subelement_id (const
                                                           t8_element * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  return pquad_w_sub->subelement_id;
}

t8_element_shape_t
  t8_subelement_scheme_quad_c::t8_element_shape (const t8_element_t * elem)
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  T8_ASSERT (t8_element_is_valid (elem));

  if (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_NO_SUBELEMENT) {
    return T8_ECLASS_QUAD;
  }
  else {                        /* for quads, all subelements are triangles */
    return T8_ECLASS_TRIANGLE;
  }
}

int
t8_subelement_scheme_quad_c::t8_element_find_neighbor_in_transition_cell
  (const t8_element_t * elem, const t8_element_t * pseudo_neigh,
   int elem_face)
{
  /* In this function, we assume pseudo_neigh to be a random subelement of a transition cell that includes
   * the real neighbor of elem. This function will output the subelement_id of the neighbor of elem. */
  t8_element_is_valid (elem);
  t8_element_is_valid (pseudo_neigh);

  /* we expect neigh to be a element in a transition cell, thus a subelement */
  T8_ASSERT (t8_element_test_if_subelement (pseudo_neigh) ==
             T8_SUB_QUAD_IS_SUBELEMENT);

  const t8_quad_with_subelements *
    pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  const t8_quad_with_subelements *
    pquad_w_sub_pseudo_neigh =
    (const t8_quad_with_subelements *) pseudo_neigh;

  if (pquad_w_sub_elem->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT && (elem_face == 0 || elem_face == 2)) {       /* we search for a neighbor within the transition cell of elem */
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
    int
      shift;
    if (elem_face == 0) {
      shift = -1;
    }
    if (elem_face == 2) {
      shift = 1;
    }
    int
      num_subelements =
      t8_element_get_number_of_subelements (pquad_w_sub_elem->subelement_type,
                                            elem);
    return ((pquad_w_sub_elem->subelement_id + shift) + num_subelements) % num_subelements;     /* the neighbor is directly before or after elem modulo the number of subelements in the transition cell */
  }
  /* Below are the cases in which the neighbor can not be identified as simple as above. 
   * The idea is to fill a location array with the desired properties of the real neighbor. 
   * Togehter with the type of the transition cell of pseudo_neigh, we can identify the sub_id of the right neighbor. */

  if (pquad_w_sub_elem->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT && elem_face == 1) {   /* elem is a subelement and pseudo_neigh is a subelement of a neighboring trnaisiton cell */
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / | \           / |
     *      |   \       /   |   \       /   |
     *      |     \   /     |     \   /     |
     *      x - - - x  N_f1 | elem  x       |
     *      |     /   \     |     / | \     |
     *      |   /pseudo \   |   /   |   \   |
     *      | /   neigh   \ | /     |     \ |
     *      x - - - - - - - x - - - x - - - x
     *
     * A subelement elem and face 1 is given as well as a random subelement pseudo_neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor N_f1. 
     * Note that both transition cells can have different levels. */

    /* get the location of elem */
    int
    location_elem[3] = { };     /* {face, is_split, number of subelement at face} */
    t8_element_get_location_of_subelement (elem, location_elem);

    /* declare the location array of the real neighbor which will be filled depending on the following cases */
    int
    location_neigh[3] = { };

    /* the pseudo_neigh tranaition cell has a lower level than the elem transition cell */
    if (pquad_w_sub_pseudo_neigh->p4q.level < pquad_w_sub_elem->p4q.level) {
      if (location_elem[0] == 0) {      /* left face of transition cell */
        if (pquad_w_sub_pseudo_neigh->p4q.y == pquad_w_sub_elem->p4q.y) {
          location_neigh[0] = 2;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* first or second subelement at face */
        }
        if (pquad_w_sub_pseudo_neigh->p4q.y != pquad_w_sub_elem->p4q.y) {

        }
      }
      if (location_elem[0] == 1) {      /* upper face of transition cell */

      }
      if (location_elem[0] == 2) {      /* right face of transition cell */

      }
      if (location_elem[0] == 3) {      /* lower face of transition cell */

      }
    }
    /* the pseudo_neigh tranaition cell has not a lower level than the elem transition cell */
    if (pquad_w_sub_pseudo_neigh->p4q.level >= pquad_w_sub_elem->p4q.level) {

    }

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return
      t8_element_get_id_from_location (t8_element_get_subelement_type
                                       (pseudo_neigh), location_neigh);
  }
  if (pquad_w_sub_elem->dummy_is_subelement == T8_SUB_QUAD_IS_NO_SUBELEMENT) {  /* elem is no subelement but its neighbor is a subelement from a neighboring transition cell. */
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / |               |
     *      |   \       /   |               |
     *      |     \   /     |               |
     *      x - - - x  N_f1 |     elem      |
     *      |     /   \     |               |
     *      |   /pseudo \   |               |
     *      | /   neigh   \ |               |
     *      x - - - - - - - x - - - - - - - x
     *
     * Subelement Elem and face 1 is given and a random subelement neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor N_f1. */

    /* declare the location array of the real neighbor which will be filled depending on the following cases */
    int
    location_neigh[3] = { };

    return
      t8_element_get_id_from_location (t8_element_get_subelement_type
                                       (pseudo_neigh), location_neigh);
  }
  T8_ASSERT ("A neighbor sub_id should have been found");
}

int
t8_subelement_scheme_quad_c::t8_element_get_id_from_location (int type,
                                                              int location[])
{
  T8_ASSERT (T8_SUB_QUAD_MIN_SUBELEMENT_TYPE <= type
             && type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE);

  int                 sub_id = 0;
  int                 type_temp = type;
  int                 binary_type[4] = { };
  int                 binary_type_clockwise[4] = { };

  /* get the type as a binary array */
  int                 i;
  for (i = 0; i < 4; i++) {
    if (type_temp >= pow (2, 4 - (i + 1))) {
      binary_type[i] = 1;
    }
    else {
      binary_type[i] = 0;
    }
    type_temp -= pow (2, 4 - (i + 1));
  }

  /* rearrange the binary type to be in clockwise order of the faces, starting with the left face */
  binary_type_clockwise[0] = binary_type[0];
  binary_type_clockwise[1] = binary_type[3];
  binary_type_clockwise[2] = binary_type[1];
  binary_type_clockwise[3] = binary_type[2];

  /* count the number of elements up to the given location */
  for (i = 0; i <= location[0]; i++) {
    if (i == location[0]) {
      if (location[1] == 0) {
        sub_id += 1;
      }
      else {
        if (location[2] == 0) {
          sub_id += 1;
        }
        else {
          sub_id += 2;
        }
      }
    }
    else {
      sub_id += binary_type_clockwise[i] + 1;
    }
  }

  /* get the sub_id */
  sub_id -= 1;

  return sub_id;
}

int
t8_subelement_scheme_quad_c::t8_element_adjust_subelement_neighbor_index
  (const t8_element_t * elem, const t8_element_t * neigh, int elem_index,
   int elem_face)
{
  const t8_quad_with_subelements *
    pquad_w_sub_elem = (const t8_quad_with_subelements *) elem;
  const t8_quad_with_subelements *
    pquad_w_sub_neigh = (const t8_quad_with_subelements *) neigh;

  /* The purpose of this function is to solve the following problem:
   * 
   *      x - - - - - - - x
   *      | \    elem   / |
   *      |   \       /   |
   *      |N_f0 \   / N_f2|
   *      x - - - + - - - x
   *      |     /   \     |
   *      |   /       \   |
   *      | /   neigh   \ |
   *      x - - - - - - - x
   * 
   * We are searching for the sibling neighbors N_f0 or N_f2 of a subelement elem but elem_index corresponds to a random sibling subelement neigh
   * instead of elem itself (this is a problem of the index function that is not unique for subelements).
   * Depending on whether we are searching N_f0 or N_f2, we can now use elem_index, and the sub_ids of elem and neigh 
   * in order to adjust elem_index and to get the right neighbor indices by just shifting it by +-1. */

  int
    adjust,
    shift,
    adjusted_index;
  int
    number_of_subelements =
    t8_element_get_number_of_subelements (pquad_w_sub_elem->subelement_type,
                                          elem);

  elem_index -= pquad_w_sub_neigh->subelement_id;       /* now we have the index of the first subelement of this transition cell */

  if (elem_face == 0) {         /* counter clockwise neighbor */
    shift = -1;
  }
  if (elem_face == 2) {         /* clockwise neighbor */
    shift = 1;
  }

  adjust =
    ((pquad_w_sub_elem->subelement_id + shift) +
     number_of_subelements) % number_of_subelements;

  return adjusted_index = elem_index + adjust;
}

void
t8_subelement_scheme_quad_c::t8_element_new (int length, t8_element_t ** elem)
{
  /* allocate memory for a quad with subelements */
  t8_default_scheme_common_c::t8_element_new (length, elem);

  int                 i;
  for (i = 0; i < length; i++) {
    t8_quad_with_subelements *pquad_w_sub =
      (t8_quad_with_subelements *) elem[i];
    t8_element_init (1, elem[i], 0);
    /* set dimension of quad to 2 */
    T8_QUAD_SET_TDIM ((p4est_quadrant_t *) & pquad_w_sub->p4q, 2);
  }
}

void
t8_subelement_scheme_quad_c::t8_element_init (int length, t8_element_t * elem,
                                              int new_called)
{
  t8_quad_with_subelements *pquad_w_sub = (t8_quad_with_subelements *) elem;

  int                 i;

  for (i = 0; i < length; i++) {
    /* initalize subelement parameters */
    pquad_w_sub[i].dummy_is_subelement = T8_SUB_QUAD_IS_NO_SUBELEMENT;
    pquad_w_sub[i].subelement_type = T8_SUB_QUAD_IS_NO_SUBELEMENT;
    pquad_w_sub[i].subelement_id = T8_SUB_QUAD_IS_NO_SUBELEMENT;

#ifdef T8_ENABLE_DEBUG
    /* In debugging mode we iterate over all length many elements and 
     * set their quad to the level 0 quad with ID 0. */
    if (!new_called) {
      p4est_quadrant_t   *quad = &pquad_w_sub[i].p4q;
      p4est_quadrant_set_morton (quad, 0, 0);
      T8_QUAD_SET_TDIM (quad, 2);
      T8_ASSERT (p4est_quadrant_is_extended (quad));
    }
  }
#endif
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_quad_c::t8_element_is_valid (const t8_element_t * elem) const 
/* *INDENT-ON* */
{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;
  const p4est_quadrant_t *q = &pquad_w_sub->p4q;

  /* the 4pest quadrant AND the subelement values must be valid such that the whole element is valid */
  return (p4est_quadrant_is_extended (q)
          && t8_element_subelement_values_are_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_quad_c::t8_element_subelement_values_are_valid (const
                                                                 t8_element_t *
                                                                 elem) const
/* *INDENT-ON* */

{
  const t8_quad_with_subelements *pquad_w_sub =
    (const t8_quad_with_subelements *) elem;

  return (pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_NO_SUBELEMENT
          || pquad_w_sub->dummy_is_subelement == T8_SUB_QUAD_IS_SUBELEMENT) &&
    ((pquad_w_sub->subelement_type >= T8_SUB_QUAD_MIN_SUBELEMENT_TYPE
      && pquad_w_sub->subelement_type <= T8_SUB_QUAD_MAX_SUBELEMENT_TYPE)
     || pquad_w_sub->subelement_type == T8_SUB_QUAD_IS_NO_SUBELEMENT) &&
    ((pquad_w_sub->subelement_id >= T8_SUB_QUAD_MIN_SUBELEMENT_ID
      && pquad_w_sub->subelement_id <= T8_SUB_QUAD_MAX_SUBELEMENT_ID)
     || pquad_w_sub->subelement_id == T8_SUB_QUAD_IS_NO_SUBELEMENT);
}
#endif

/* Constructor */
t8_subelement_scheme_quad_c::t8_subelement_scheme_quad_c (void)
{
  eclass = T8_ECLASS_QUAD;
  element_size = sizeof (t8_pquad_t);
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
