/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_scheme.h>
#include <t8_forest/t8_forest_types.h>

void
t8_scheme_ref (t8_scheme_c *scheme)
{
  T8_ASSERT (scheme != NULL);

  scheme->ref ();
}

void
t8_scheme_unref (t8_scheme_c **pscheme)
{
  T8_ASSERT (pscheme != NULL);

  if ((*pscheme)->unref () < 1) {
    *pscheme = NULL;
  }
}

size_t
t8_element_get_element_size (const t8_scheme_c *scheme, const t8_eclass_t tree_class)
{
  return scheme->get_element_size (tree_class);
}

int
t8_element_refines_irregular (const t8_scheme_c *scheme, const t8_eclass_t tree_class)
{
  return scheme->refines_irregular (tree_class);
}

int
t8_element_get_maxlevel (const t8_scheme_c *scheme, const t8_eclass_t tree_class)
{
  return scheme->get_maxlevel (tree_class);
}

int
t8_element_get_level (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_level (tree_class, elem);
}

void
t8_element_copy (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *source,
                 t8_element_t *dest)
{
  return scheme->element_copy (tree_class, source, dest);
}

int
t8_element_compare (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem1,
                    const t8_element_t *elem2)
{
  return scheme->element_compare (tree_class, elem1, elem2);
}

int
t8_element_is_equal (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem1,
                     const t8_element_t *elem2)
{
  return scheme->element_is_equal (tree_class, elem1, elem2);
}

int
element_is_refinable (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_is_refinable (tree_class, elem);
};

void
t8_element_get_parent (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                       t8_element_t *parent)
{
  return scheme->element_get_parent (tree_class, elem, parent);
}

int
t8_element_get_num_siblings (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_num_siblings (tree_class, elem);
}

void
t8_element_get_sibling (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                        const int sibid, t8_element_t *sibling)
{
  return scheme->element_get_sibling (tree_class, elem, sibid, sibling);
}

int
t8_element_get_num_corners (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_num_corners (tree_class, elem);
}

int
t8_element_get_num_faces (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_num_faces (tree_class, elem);
}

int
t8_element_get_max_num_faces (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_max_num_faces (tree_class, elem);
}

int
t8_element_get_num_children (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_num_children (tree_class, elem);
}

int
t8_get_num_children (const t8_scheme_c *scheme, const t8_eclass_t tree_class)
{
  return scheme->get_max_num_children (tree_class);
}

int
t8_element_get_num_face_children (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                  int face)
{
  return scheme->element_get_num_face_children (tree_class, elem, face);
}

int
t8_element_get_face_corner (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                            const int face, const int corner)
{
  return scheme->element_get_face_corner (tree_class, elem, face, corner);
}

int
t8_element_get_corner_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                            const int corner, const int face)
{
  return scheme->element_get_corner_face (tree_class, elem, corner, face);
}

void
t8_element_get_child (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                      const int childid, t8_element_t *child)
{
  return scheme->element_get_child (tree_class, elem, childid, child);
}

void
t8_element_get_children (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                         const int length, t8_element_t *c[])
{
  return scheme->element_get_children (tree_class, elem, length, c);
}

int
t8_element_get_child_id (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_child_id (tree_class, elem);
}

int
t8_element_get_ancestor_id (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                            const int level)
{
  return scheme->element_get_ancestor_id (tree_class, elem, level);
}

int
t8_elements_are_family (const t8_scheme_c *scheme, const t8_eclass_t tree_class, t8_element_t *const *fam)
{
  return scheme->elements_are_family (tree_class, fam);
}

void
t8_element_get_nca (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem1,
                    const t8_element_t *elem2, t8_element_t *nca)
{
  return scheme->element_get_nca (tree_class, elem1, elem2, nca);
}

t8_element_shape_t
t8_element_get_face_shape (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                           const int face)
{
  return scheme->element_get_face_shape (tree_class, elem, face);
}

void
t8_element_get_children_at_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                 const int face, t8_element_t *children[], const int num_children, int *child_indices)
{
  return scheme->element_get_children_at_face (tree_class, elem, face, children, num_children, child_indices);
}

int
t8_element_face_get_child_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                const int face, const int face_child)
{
  return scheme->element_face_get_child_face (tree_class, elem, face, face_child);
}

int
t8_element_face_get_parent_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                 const int face)
{
  return scheme->element_face_get_parent_face (tree_class, elem, face);
}

int
t8_element_get_tree_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                          const int face)
{
  return scheme->element_get_tree_face (tree_class, elem, face);
}

void
t8_element_transform_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem1,
                           t8_element_t *elem2, const int orientation, const int sign, const int is_smaller_face)
{
  return scheme->element_transform_face (tree_class, elem1, elem2, orientation, sign, is_smaller_face);
}

int
t8_element_extrude_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *face,
                         t8_element_t *elem, const int root_face)
{
  return scheme->element_extrude_face (tree_class, face, elem, root_face);
}

void
t8_element_get_boundary_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                              const int face, t8_element_t *boundary)
{
  return scheme->element_get_boundary_face (tree_class, elem, face, boundary);
}

void
t8_element_get_first_descendant_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                      int face, t8_element_t *first_desc, const int level)
{
  return scheme->element_get_first_descendant_face (tree_class, elem, face, first_desc, level);
}

void
t8_element_get_last_descendant_face (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                     const int face, t8_element_t *last_desc, const int level)
{
  return scheme->element_get_last_descendant_face (tree_class, elem, face, last_desc, level);
}

int
t8_element_is_root_boundary (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                             const int face)
{
  return scheme->element_is_root_boundary (tree_class, elem, face);
}

int
t8_element_get_face_neighbor_inside (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                     t8_element_t *neigh, const int face, int *neigh_face)
{
  return scheme->element_get_face_neighbor_inside (tree_class, elem, neigh, face, neigh_face);
}

t8_element_shape_t
t8_element_get_shape (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_get_shape (tree_class, elem);
}

void
t8_element_set_linear_id (const t8_scheme_c *scheme, const t8_eclass_t tree_class, t8_element_t *elem, int level,
                          t8_linearidx_t id)
{
  return scheme->element_set_linear_id (tree_class, elem, level, id);
}

t8_linearidx_t
t8_element_get_linear_id (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                          const int level)
{
  return scheme->element_get_linear_id (tree_class, elem, level);
}

void
t8_element_get_first_descendant (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                 t8_element_t *desc, int level)
{
  return scheme->element_get_first_descendant (tree_class, elem, desc, level);
}

void
t8_element_get_last_descendant (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                                t8_element_t *desc, const int level)
{
  return scheme->element_get_last_descendant (tree_class, elem, desc, level);
}

void
t8_element_construct_successor (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem1,
                                t8_element_t *elem2)
{
  return scheme->element_construct_successor (tree_class, elem1, elem2);
}

void
t8_element_get_vertex_reference_coords (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *t,
                                        const int vertex, double coords[])
{
  return scheme->element_get_vertex_reference_coords (tree_class, t, vertex, coords);
}

void
t8_element_get_reference_coords (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *t,
                                 const double *ref_coords, const size_t num_coords, double coords[])
{
  return scheme->element_get_reference_coords (tree_class, t, ref_coords, num_coords, coords);
}

t8_gloidx_t
t8_element_count_leaves (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *t,
                         const int level)
{
  return scheme->element_count_leaves (tree_class, t, level);
}

t8_gloidx_t
t8_element_count_leaves_from_root (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const int level)
{
  return scheme->count_leaves_from_root (tree_class, level);
}

#if T8_ENABLE_DEBUG
int
t8_element_is_valid (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_is_valid (tree_class, elem);
}

void
t8_element_debug_print (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem)
{
  return scheme->element_debug_print (tree_class, elem);
}

#endif
void
t8_element_to_string (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t *elem,
                      char *debug_string, const int string_size)
{
  T8_ASSERT (debug_string != NULL);

  return scheme->element_to_string (tree_class, elem, debug_string, string_size);
}

void
t8_element_new (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const int length, t8_element_t **elems)
{
  return scheme->element_new (tree_class, length, elems);
}

void
t8_element_init (const t8_scheme_c *scheme, const t8_eclass_t tree_class, int length, t8_element_t *elems)
{
  return scheme->element_init (tree_class, length, elems);
}

void
t8_element_deinit (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const int length, t8_element_t *elems)
{
  return scheme->element_deinit (tree_class, length, elems);
}
void
t8_element_destroy (const t8_scheme_c *scheme, const t8_eclass_t tree_class, int length, t8_element_t **elems)
{
  return scheme->element_destroy (tree_class, length, elems);
}

void
t8_element_set_to_root (const t8_scheme_c *scheme, const t8_eclass_t tree_class, t8_element_t *elem)
{
  return scheme->set_to_root (tree_class, elem);
}

void
t8_element_MPI_Pack (const t8_scheme_c *scheme, const t8_eclass_t tree_class, t8_element_t **const elements,
                     const unsigned int count, void *send_buffer, const int buffer_size, int *position,
                     sc_MPI_Comm comm)
{
  return scheme->element_MPI_Pack (tree_class, elements, count, send_buffer, buffer_size, position, comm);
}

void
t8_element_MPI_Pack_size (const t8_scheme_c *scheme, const t8_eclass_t tree_class, const unsigned int count,
                          sc_MPI_Comm comm, int *pack_size)
{
  return scheme->element_MPI_Pack_size (tree_class, count, comm, pack_size);
}

void
t8_element_MPI_Unpack (const t8_scheme_c *scheme, const t8_eclass_t tree_class, void *recvbuf, const int buffer_size,
                       int *position, t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm)
{
  return scheme->element_MPI_Unpack (tree_class, recvbuf, buffer_size, position, elements, count, comm);
}
