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

#include "t8_different_num_child_line_cxx.hxx"
#include "../../t8_default/t8_default_common_cxx.hxx"
#include "../../t8_default/t8_default_vertex_cxx.hxx"
#include "../../t8_default/t8_dvertex_bits.h"
#include "../../t8_default/t8_dline_bits.h"
#include "../../t8_default/t8_dline.h"

typedef t8_dline_t  t8_default_line_t;

T8_EXTERN_C_BEGIN ();

/* Set the x coordinate of a line element to zero */
static void
t8_set_line_coord_zero (t8_element_t * elem)
{
  ((t8_dline_t *) elem)->x = 0;
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_children_at_face (const t8_element_t * elem, int face,
                               t8_element_t * children[], int num_children,
                               int *child_indices)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Call base class function */
  t8_default_scheme_line_c::t8_element_children_at_face (elem, face, children,
                                                         num_children,
                                                         child_indices);

  /* Explicitely set coordinate to 0 */
  t8_set_line_coord_zero (children[0]);
}

int
t8_different_num_child_scheme_line_c::t8_element_tree_face (const t8_element_t
                                                            * elem, int face)
{
  /* Since every face is a tree boundary and the indices are the same for elements
   * and tree, we just return face */
  return face;
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_transform_face (const t8_element_t * elem1, t8_element_t * elem2,
                             int orientation, int sign, int is_smaller_face)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (orientation == 0 || orientation == 1);
  /* Just copy, since all lines look the same */
  t8_element_copy (elem1, elem2);
}

/** Given a boundary face inside a root tree's face construct
 *  the element inside the root tree that has the given face as a
 *  face. */
/* *INDENT-OFF* */
int
t8_different_num_child_scheme_line_c
::t8_element_extrude_face (const t8_element_t *face,
                           const t8_eclass_scheme_c *face_scheme,
                           t8_element_t *elem, int root_face)
/* *INDENT-ON* */
{
  T8_ASSERT (t8_element_is_valid (face));
  /* Call base class */
  int                 face_number =
    t8_default_scheme_line_c::t8_element_extrude_face (face, face_scheme,
                                                       elem, root_face);
  /* Set coordinate to 0 */
  t8_set_line_coord_zero (elem);
  return face_number;
}

/** Construct the boundary element at a specific face. */
/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_boundary_face (const t8_element_t *elem, int face,
                            t8_element_t *boundary,
                            const t8_eclass_scheme_c *boundary_scheme)
/* *INDENT-ON* */
{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Call base class */
  t8_default_scheme_line_c::t8_element_boundary_face (elem, face, boundary,
                                                      boundary_scheme);
  /* Set coordinate to zero */
  t8_set_line_coord_zero (boundary);
}

/** Construct the first descendant of an element that touches a given face.   */
/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_first_descendant_face (const *elem, int face,
                                    t8_element_t * first_desc, int level)
/* *INDENT-ON* */
{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Call base class */
  t8_default_scheme_line_c::t8_element_first_descendant_face (elem, face,
                                                              first_desc,
                                                              level);
  /* Set coordinate to zero */
  t8_set_line_coord_zero (first_desc);
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_last_descendant_face (const t8_element_t * elem, int face,
                                   t8_element_t * last_desc, int level)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Call base class */
  t8_default_scheme_line_c::t8_element_last_descendant_face (elem, face,
                                                             last_desc,
                                                             level);
  /* Set coordinate to zero */
  t8_set_line_coord_zero (last_desc);
}

/* *INDENT-OFF* */
int
t8_different_num_child_scheme_line_c
::t8_element_is_root_boundary (const t8_element_t * elem, int face)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));

  /* Every element is a root boundary */
  return 1;
}

/* *INDENT-OFF* */
int
t8_different_num_child_scheme_line_c
::t8_element_face_neighbor_inside (const t8_element_t * elem,
                                   t8_element_t * neigh, int face,
                                   int *neigh_face)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DLINE_FACES);

  /* No neighbor is inside */
  return 0;
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_set_linear_id (t8_element_t * elem, int level, t8_linearidx_t id)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);
  T8_ASSERT (id == 0);          /* We only have one element per level */

  t8_dline_init_linear_id ((t8_default_line_t *) elem, level, id);
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_first_descendant (const t8_element_t * elem, t8_element_t * desc,
                               int level)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Copy the element */
  t8_element_copy (elem, desc);
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_last_descendant (const t8_element_t * elem, t8_element_t * desc,
                              int level)
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Copy the element */
  t8_element_copy (elem, desc);
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_vertex_coords (const t8_element_t * t, int vertex, int coords[])
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (t));
  /* The vertex coords are the same as for the level 0 element */
  t8_dline_t          line;
  t8_dline_init_linear_id (&line, 0, 0);
  t8_dline_vertex_coords (&line, vertex, coords);
}

t8_linearidx_t
  t8_different_num_child_scheme_line_c::t8_element_get_linear_id (const
                                                                  t8_element_t
                                                                  * elem,
                                                                  int level)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DLINE_MAXLEVEL);

  /* All elements have Id 0, since there is only one element per level */
  return 0;
}

int
t8_different_num_child_scheme_line_c::t8_element_num_children (const
                                                               t8_element_t *
                                                               elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  return T8_LINE_NUM_CHILD_CHILDREN;
}

int
t8_different_num_child_scheme_line_c::t8_element_child_id (const t8_element_t
                                                           * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  /* Every child is child 0 */
  return 0;
}

/* *INDENT-OFF* */
void
t8_different_num_child_scheme_line_c
::t8_element_children (const t8_element_t * elem, int length,
                       t8_element_t * c[])
/* *INDENT-ON* */

{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == T8_LINE_NUM_CHILD_CHILDREN);

  t8_dline_child ((const t8_dline_t *) elem, 0, (t8_dline_t *) c[0]);
  T8_ASSERT (t8_element_is_valid (c[0]));
}

int
t8_different_num_child_scheme_line_c::t8_element_is_family (t8_element_t **
                                                            fam)
{
#ifdef T8_ENABLE_DEBUG
  int                 i;
  for (i = 0; i < T8_LINE_NUM_CHILD_CHILDREN; i++) {
    T8_ASSERT (t8_element_is_valid (fam[i]));
  }
#endif
  /* Every 'group' of 1 element is a family */
  return 1;
}

#ifdef T8_ENABLE_DEBUG
/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_different_num_child_scheme_line_c::t8_element_is_valid (const t8_element_t * elem) const
/* *INDENT-ON* */
{
  const t8_dline_t   *line = (const t8_dline_t *) elem;
  return t8_dline_is_valid (line) && line->x == 0;
}
#endif

t8_different_num_child_scheme_line_c::~t8_different_num_child_scheme_line_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the line_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
