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
 * t8_scheme class.
 * With this interface you can use these member functions from a C file
 * without the need of compiling it with C++.
 */

#include <t8_element.h>
#include <t8_element.hxx>
#include <t8_element_c_interface.h>
#include <t8_scheme.hxx>
#include <t8_forest/t8_forest_types.h>

size_t
t8_element_get_element_size (const t8_forest_t forest, const t8_eclass_t tree_class)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->get_element_size (tree_class);
}

int
t8_element_refines_irregular (const t8_forest_t forest, const t8_eclass_t tree_class)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->refines_irregular (tree_class);
}

int
t8_element_get_maxlevel (const t8_forest_t forest, const t8_eclass_t tree_class)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->get_maxlevel (tree_class);
}

int
t8_element_get_level (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_level (tree_class, elem);
}

void
t8_element_copy (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *source, t8_element_t *dest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_copy (tree_class, source, dest);
}

int
t8_element_compare (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem1,
                    const t8_element_t *elem2)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_compare (tree_class, elem1, elem2);
}

int
t8_element_is_equal (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem1,
                     const t8_element_t *elem2)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_is_equal (tree_class, elem1, elem2);
}

void
t8_element_get_parent (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                       t8_element_t *parent)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_parent (tree_class, elem, parent);
}

int
t8_element_get_num_siblings (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_num_siblings (tree_class, elem);
}

void
t8_element_get_sibling (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int sibid,
                        t8_element_t *sibling)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_sibling (tree_class, elem, sibid, sibling);
}

int
t8_element_get_num_corners (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_num_corners (tree_class, elem);
}

int
t8_element_get_num_faces (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_num_faces (tree_class, elem);
}

int
t8_element_get_max_num_faces (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_max_num_faces (tree_class, elem);
}

int
t8_element_get_num_children (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_num_children (tree_class, elem);
}

int
t8_element_get_num_face_children (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                  int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_num_face_children (tree_class, elem, face);
}

int
t8_element_get_face_corner (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int face,
                            int corner)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_face_corner (tree_class, elem, face, corner);
}

int
t8_element_get_corner_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                            int corner, int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_corner_face (tree_class, elem, corner, face);
}

void
t8_element_get_child (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int childid,
                      t8_element_t *child)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_child (tree_class, elem, childid, child);
}

void
t8_element_get_children (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int length,
                         t8_element_t *c[])
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_children (tree_class, elem, length, c);
}

int
t8_element_get_child_id (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_child_id (tree_class, elem);
}

int
t8_element_get_ancestor_id (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_ancestor_id (tree_class, elem, level);
}

int
t8_elements_are_family (const t8_forest_t forest, const t8_eclass_t tree_class, t8_element_t *const *fam)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->elements_are_family (tree_class, fam);
}

void
t8_element_get_nca (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem1,
                    const t8_element_t *elem2, t8_element_t *nca)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_nca (tree_class, elem1, elem2, nca);
}

t8_element_shape_t
t8_element_get_face_shape (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_face_shape (tree_class, elem, face);
}

void
t8_element_get_children_at_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                 int face, t8_element_t *children[], int num_children, int *child_indices)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_children_at_face (tree_class, elem, face, children, num_children, child_indices);
}

int
t8_element_face_get_child_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                int face, int face_child)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_face_get_child_face (tree_class, elem, face, face_child);
}

int
t8_element_face_get_parent_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                 int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_face_get_parent_face (tree_class, elem, face);
}

int
t8_element_get_tree_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_tree_face (tree_class, elem, face);
}

void
t8_element_transform_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem1,
                           t8_element_t *elem2, int orientation, int sign, int is_smaller_face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_transform_face (tree_class, elem1, elem2, orientation, sign, is_smaller_face);
}

int
t8_element_extrude_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *face,
                         const t8_eclass_t face_eclass, t8_element_t *elem, int root_face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_extrude_face (tree_class, face, face_scheme, elem, root_face);
}

void
t8_element_construct_boundary_face (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                    int face, t8_element_t *boundary, const t8_eclass_t boundary_face_eclass)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_boundary_face (tree_class, elem, face, boundary, boundary_scheme);
}

void
t8_element_construct_first_descendant_face (const t8_forest_t forest, const t8_eclass_t tree_class,
                                            const t8_element_t *elem, int face, t8_element_t *first_desc, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_first_descendant_face (tree_class, elem, face, first_desc, level);
}

void
t8_element_construct_last_descendant_face (const t8_forest_t forest, const t8_eclass_t tree_class,
                                           const t8_element_t *elem, int face, t8_element_t *last_desc, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_last_descendant_face (tree_class, elem, face, last_desc, level);
}

int
t8_element_is_root_boundary (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_is_root_boundary (tree_class, elem, face);
}

int
t8_element_construct_face_neighbor_inside (const t8_forest_t forest, const t8_eclass_t tree_class,
                                           const t8_element_t *elem, t8_element_t *neigh, int face, int *neigh_face)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_face_neighbor_inside (tree_class, elem, neigh, face, neigh_face);
}

t8_element_shape_t
t8_element_get_shape (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_shape (tree_class, elem);
}

void
t8_element_set_linear_id (const t8_forest_t forest, const t8_eclass_t tree_class, t8_element_t *elem, int level,
                          t8_linearidx_t id)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_set_linear_id (tree_class, elem, level, id);
}

t8_linearidx_t
t8_element_get_linear_id (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_linear_id (tree_class, elem, level);
}

void
t8_element_construct_first_descendant (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                       t8_element_t *desc, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_first_descendant (tree_class, elem, desc, level);
}

void
t8_element_construct_last_descendant (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                                      t8_element_t *desc, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_last_descendant (tree_class, elem, desc, level);
}

void
t8_element_construct_successor (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem1,
                                t8_element_t *elem2)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_construct_successor (tree_class, elem1, elem2);
}

void
t8_element_get_vertex_reference_coords (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *t,
                                        const int vertex, double coords[])
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_vertex_reference_coords (tree_class, t, vertex, coords);
}

void
t8_element_get_reference_coords (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *t,
                                 const double *ref_coords, const size_t num_coords, double coords[])
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_get_reference_coords (tree_class, t, vertex, coords);
}

t8_gloidx_t
t8_element_count_leaves (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *t, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_count_leaves (tree_class, t, level);
}

t8_gloidx_t
t8_element_count_leaves_from_root (const t8_forest_t forest, const t8_eclass_t tree_class, int level)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_count_leaves_from_root (tree_class, level);
}

#ifdef T8_ENABLE_DEBUG
int
t8_element_is_valid (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_is_valid (tree_class, elem);
}

void
t8_element_debug_print (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_debug_print (tree_class, elem);
}

void
t8_element_to_string (const t8_forest_t forest, const t8_eclass_t tree_class, const t8_element_t *elem,
                      char *debug_string, const int string_size)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  T8_ASSERT (debug_string != NULL);

  return forest->scheme->element_to_string (tree_class, elem, debug_string, string_size);
}
#endif

void
t8_element_new (const t8_forest_t forest, const t8_eclass_t tree_class, int length, t8_element_t **elems)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_new (tree_class, length, elems);
}

void
t8_element_init (const t8_forest_t forest, const t8_eclass_t tree_class, int length, t8_element_t *elems)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_init (tree_class, length, elems);
}

void
t8_element_deinit (const t8_forest_t forest, const t8_eclass_t tree_class, int length, t8_element_t *elems)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_deinit (tree_class, length, elems);
}
void
t8_element_destroy (const t8_forest_t forest, const t8_eclass_t tree_class, int length, t8_element_t **elems)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_destroy (tree_class, length, elems);
}

void
t8_element_get_root (const t8_forest_t forest, const t8_eclass_t tree_class, t8_element_t *elem)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->get_root_element (tree_class, elem);
}

void
t8_element_MPI_Pack (const t8_forest_t forest, const t8_eclass_t tree_class, t8_element_t **const elements,
                     const unsigned int count, void *send_buffer, const int buffer_size, int *position,
                     sc_MPI_Comm comm)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_MPI_Pack (tree_class, elements, count, send_buffer, buffer_size, position, comm);
}

void
t8_element_MPI_Pack_size (const t8_forest_t forest, const t8_eclass_t tree_class, const unsigned int count,
                          sc_MPI_Comm comm, int *pack_size)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_MPI_Pack_size (tree_class, count, comm, pack_size);
}

void
t8_element_MPI_Unpack (const t8_forest_t forest, const t8_eclass_t tree_class, void *recvbuf, const int buffer_size,
                       int *position, t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  return forest->scheme->element_MPI_Unpack (tree_class, recvbuf, buffer_size, position, elements, count, comm);
}
