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

/* In this file we implement a C interface for the member functions of the
 * t8_eclass_scheme_c class.
 * With this interface you can use these member functions from a C file
 * without the need of compiling it with C++.
 */

#include <t8_element.h>
#include <t8_element_cxx.hxx>
#include <t8_element_c_interface.h>

size_t
t8_element_size (const t8_eclass_scheme_c *ts)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_size ();
}

int
t8_element_refines_irregular (const t8_eclass_scheme_c *ts)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_refines_irregular ();
}

int
t8_element_maxlevel (const t8_eclass_scheme_c *ts)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_maxlevel ();
}

int
t8_element_level (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_level (elem);
}

void
t8_element_copy (const t8_eclass_scheme_c *ts, const t8_element_t *source, t8_element_t *dest)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_copy (source, dest);
}

int
t8_element_compare (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_compare (elem1, elem2);
}

int
t8_element_equal (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_equal (elem1, elem2);
}

void
t8_element_parent (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *parent)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_parent (elem, parent);
}

int
t8_element_num_siblings (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_siblings (elem);
}

void
t8_element_sibling (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int sibid, t8_element_t *sibling)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_sibling (elem, sibid, sibling);
}

int
t8_element_num_corners (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_corners (elem);
}

int
t8_element_num_faces (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_faces (elem);
}

int
t8_element_max_num_faces (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_max_num_faces (elem);
}

int
t8_element_num_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_children (elem);
}

int
t8_element_num_face_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_num_face_children (elem, face);
}

int
t8_element_get_face_corner (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, int corner)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_get_face_corner (elem, face, corner);
}

int
t8_element_get_corner_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int corner, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_get_corner_face (elem, corner, face);
}

void
t8_element_child (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int childid, t8_element_t *child)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_child (elem, childid, child);
}

void
t8_element_children (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int length, t8_element_t *c[])
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_children (elem, length, c);
}

int
t8_element_child_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_child_id (elem);
}

int
t8_element_ancestor_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int level)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_ancestor_id (elem, level);
}

int
t8_element_is_family (const t8_eclass_scheme_c *ts, t8_element_t *const *fam)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_is_family (fam);
}

void
t8_element_nca (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_nca (elem1, elem2, nca);
}

t8_element_shape_t
t8_element_face_shape (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_face_shape (elem, face);
}

void
t8_element_children_at_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, t8_element_t *children[],
                             int num_children, int *child_indices)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_children_at_face (elem, face, children, num_children, child_indices);
}

int
t8_element_face_child_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, int face_child)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_face_child_face (elem, face, face_child);
}

int
t8_element_face_parent_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_face_parent_face (elem, face);
}

int
t8_element_tree_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_tree_face (elem, face);
}

void
t8_element_transform_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, t8_element_t *elem2,
                           int orientation, int sign, int is_smaller_face)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_transform_face (elem1, elem2, orientation, sign, is_smaller_face);
}

int
t8_element_extrude_face (const t8_eclass_scheme_c *ts, const t8_element_t *face, const t8_eclass_scheme_c *face_scheme,
                         t8_element_t *elem, int root_face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_extrude_face (face, face_scheme, elem, root_face);
}

void
t8_element_boundary_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face, t8_element_t *boundary,
                          const t8_eclass_scheme_c *boundary_scheme)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_boundary_face (elem, face, boundary, boundary_scheme);
}

void
t8_element_first_descendant_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face,
                                  t8_element_t *first_desc, int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_first_descendant_face (elem, face, first_desc, level);
}

void
t8_element_last_descendant_face (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face,
                                 t8_element_t *last_desc, int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_last_descendant_face (elem, face, last_desc, level);
}

void
t8_element_set_linear_id (const t8_eclass_scheme_c *ts, t8_element_t *elem, const int level, const t8_linearidx_t id,
                          const int multilevel)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_set_linear_id (elem, level, id, multilevel);
}

int
t8_element_is_root_boundary (const t8_eclass_scheme_c *ts, const t8_element_t *elem, int face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_is_root_boundary (elem, face);
}

int
t8_element_face_neighbor_inside (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *neigh, int face,
                                 int *neigh_face)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_face_neighbor_inside (elem, neigh, face, neigh_face);
}

t8_element_shape_t
t8_element_shape (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_shape (elem);
}

t8_linearidx_t
t8_element_get_linear_id (const t8_eclass_scheme_c *ts, const t8_element_t *elem, const int level, const int multilevel)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_get_linear_id (elem, level, multilevel);
}

void
t8_element_first_descendant (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *desc, int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_first_descendant (elem, desc, level);
}

void
t8_element_last_descendant (const t8_eclass_scheme_c *ts, const t8_element_t *elem, t8_element_t *desc, int level)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_last_descendant (elem, desc, level);
}

void
t8_element_successor (const t8_eclass_scheme_c *ts, const t8_element_t *elem1, t8_element_t *elem2, const int level,
                      const int multilevel)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_successor (elem1, elem2, level, multilevel);
}

void
t8_element_vertex_reference_coords (const t8_eclass_scheme_c *ts, const t8_element_t *t, const int vertex,
                                    double coords[])
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_vertex_reference_coords (t, vertex, coords);
}

t8_gloidx_t
t8_element_count_leaves (const t8_eclass_scheme_c *ts, const t8_element_t *t, int level)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_count_leaves (t, level);
}

t8_gloidx_t
t8_element_count_leaves_from_root (const t8_eclass_scheme_c *ts, int level)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_count_leaves_from_root (level);
}

t8_gloidx_t
t8_element_count_elements_from_root (const t8_eclass_scheme_c *ts, int level)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_count_elements_from_root (level);
}

#ifdef T8_ENABLE_DEBUG
int
t8_element_is_valid (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_is_valid (elem);
}

void
t8_element_debug_print (const t8_eclass_scheme_c *ts, const t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_debug_print (elem);
}

void
t8_element_to_string (const t8_eclass_scheme_c *ts, const t8_element_t *elem, char *debug_string, const int string_size)
{
  T8_ASSERT (ts != NULL);
  T8_ASSERT (debug_string != NULL);

  ts->t8_element_to_string (elem, debug_string, string_size);
}
#endif

void
t8_element_new (const t8_eclass_scheme_c *ts, int length, t8_element_t **elems)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_new (length, elems);
}

void
t8_element_destroy (const t8_eclass_scheme_c *ts, int length, t8_element_t **elems)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_destroy (length, elems);
}

void
t8_element_root (const t8_eclass_scheme_c *ts, t8_element_t *elem)
{
  T8_ASSERT (ts != NULL);

  ts->t8_element_root (elem);
}

void
t8_element_MPI_Pack (const t8_eclass_scheme_c *ts, t8_element_t **const elements, const unsigned int count,
                     void *send_buffer, const int buffer_size, int *position, sc_MPI_Comm comm)
{
  T8_ASSERT (ts != NULL);

  return ts->t8_element_MPI_Pack (elements, count, send_buffer, buffer_size, position, comm);
}

void
t8_element_MPI_Pack_size (const t8_eclass_scheme_c *ts, const unsigned int count, sc_MPI_Comm comm, int *pack_size)
{
  T8_ASSERT (ts != NULL);
  return ts->t8_element_MPI_Pack_size (count, comm, pack_size);
}

void
t8_element_MPI_Unpack (const t8_eclass_scheme_c *ts, void *recvbuf, const int buffer_size, int *position,
                       t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm)
{
  T8_ASSERT (ts != NULL);
  return ts->t8_element_MPI_Unpack (recvbuf, buffer_size, position, elements, count, comm);
}
