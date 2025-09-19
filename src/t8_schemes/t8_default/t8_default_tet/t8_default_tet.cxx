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

#include <t8_schemes/t8_default/t8_default_common/t8_default_common.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
#include <t8_schemes/t8_scheme.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

typedef t8_dtet_t t8_default_tet_t;

size_t
t8_default_scheme_tet::get_element_size (void) const
{
  return sizeof (t8_dtet_t);
}

int
t8_default_scheme_tet::get_maxlevel (void) const
{
  return T8_DTET_MAXLEVEL;
}

int
t8_default_scheme_tet::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dtet_get_level ((t8_dtet_t *) elem);
}

void
t8_default_scheme_tet::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  T8_ASSERT (element_is_valid (dest));
  t8_dtet_copy ((const t8_dtet_t *) source, (t8_dtet_t *) dest);
}

int
t8_default_scheme_tet::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dtet_compare ((const t8_dtet_t *) elem1, (const t8_dtet_t *) elem2);
}

int
t8_default_scheme_tet::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dtet_equal ((const t8_dtet_t *) elem1, (const t8_dtet_t *) elem2);
}

void
t8_default_scheme_tet::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t *p = (t8_default_tet_t *) parent;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (parent));
  t8_dtet_parent (t, p);
}

void
t8_default_scheme_tet::element_get_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t *s = (t8_default_tet_t *) sibling;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (sibling));
  t8_dtet_sibling (t, sibid, s);
}

int
t8_default_scheme_tet::element_get_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTET_FACES;
}

int
t8_default_scheme_tet::element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  return T8_DTET_FACES;
}

int
t8_default_scheme_tet::element_get_num_children ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTET_CHILDREN;
}

int
t8_default_scheme_tet::element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem,
                                                      [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DTET_FACE_CHILDREN;
}

int
t8_default_scheme_tet::element_get_face_corner ([[maybe_unused]] const t8_element_t *element, int face,
                                                int corner) const
{
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (0 <= corner && corner < 3);
  return t8_dtet_face_corner[face][corner];
}

void
t8_default_scheme_tet::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t *c = (t8_default_tet_t *) child;
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (child));

  t8_dtet_child (t, childid, c);
}

void
t8_default_scheme_tet::element_get_children (const t8_element_t *elem, [[maybe_unused]] int length,
                                             t8_element_t *c[]) const
{
  T8_ASSERT (length == T8_DTET_CHILDREN);
  T8_ASSERT (element_is_valid (elem));
#if T8_ENABLE_DEBUG
  for (int i = 0; i < T8_DTET_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (c[i]));
  }
#endif
  t8_dtet_childrenpv ((const t8_dtet_t *) elem, (t8_dtet_t **) c);
}

int
t8_default_scheme_tet::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dtet_child_id ((const t8_dtet_t *) elem);
}

int
t8_default_scheme_tet::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  return t8_dtet_ancestor_id ((t8_dtet_t *) elem, level);
}

int
t8_default_scheme_tet::elements_are_family (t8_element_t *const *fam) const
{
#if T8_ENABLE_DEBUG
  for (int i = 0; i < T8_DTET_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (fam[i]));
  }
#endif
  return t8_dtet_is_familypv ((const t8_dtet_t **) fam);
}

void
t8_default_scheme_tet::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  const t8_default_tet_t *t1 = (const t8_default_tet_t *) elem1;
  const t8_default_tet_t *t2 = (const t8_default_tet_t *) elem2;
  t8_default_tet_t *c = (t8_default_tet_t *) nca;

  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  t8_dtet_nearest_common_ancestor (t1, t2, c);
}

t8_element_shape_t
t8_default_scheme_tet::element_get_face_shape ([[maybe_unused]] const t8_element_t *elem,
                                               [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_ECLASS_TRIANGLE;
}

void
t8_default_scheme_tet::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                     int num_children, int *child_indices) const
{
  const t8_dtet_t *t = (const t8_dtet_t *) elem;
  t8_dtet_t **c = (t8_dtet_t **) children;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (num_children == T8_DTET_FACE_CHILDREN);

#if T8_ENABLE_DEBUG
  /* debugging check that all children elements are valid */
  {
    for (int i = 0; i < num_children; i++) {
      T8_ASSERT (element_is_valid (children[i]));
    }
  }
#endif

  t8_dtet_children_at_face (t, face, c, num_children, child_indices);
}

int
t8_default_scheme_tet::element_face_get_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (0 <= face && face < T8_DTET_FACE_CHILDREN);
  return t8_dtet_face_child_face ((const t8_dtet_t *) elem, face, face_child);
}

int
t8_default_scheme_tet::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);

  return t8_dtet_face_parent_face ((const t8_dtet_t *) elem, face);
}

int
t8_default_scheme_tet::element_get_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  return t8_dtet_tree_face ((t8_dtet_t *) elem, face);
}

/* Construct the inner element from a boundary element. */
/* This function is defined here instead of in t8_dri_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtri_bits.c. This would be needed by an implementation, at least
 * for tets. */
int
t8_default_scheme_tet::element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face,
                                             [[maybe_unused]] const t8_scheme *scheme) const
{
  const t8_dtri_t *b = (const t8_dtri_t *) face;
  t8_dtet_t *t = (t8_dtet_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_TRIANGLE, face));
  T8_ASSERT (0 <= root_face && root_face < T8_DTET_FACES);
  t->level = b->level;
  switch (root_face) {
    /* Since the root triangle may have a different scale then the
     * root tetrahedron, we have to rescale the coordinates. */
  case 0:
    t->type = b->type == 0 ? 0 : 1;
    t->x = T8_DTET_ROOT_LEN - T8_DTET_LEN (t->level);
    t->y = ((int64_t) b->y * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    t->z = ((int64_t) b->x * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    break;
  case 1:
    t->type = b->type == 0 ? 0 : 2;
    t->x = t->z = ((int64_t) b->x * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    t->y = ((int64_t) b->y * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    break;
  case 2:
    t->type = b->type == 0 ? 0 : 4;
    t->x = ((int64_t) b->x * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    t->y = t->z = ((int64_t) b->y * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    break;
  case 3:
    t->type = b->type == 0 ? 0 : 5;
    t->x = ((int64_t) b->x * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    t->y = 0;
    t->z = ((int64_t) b->y * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* Compute the face as seen from t */
  return t8_dtet_root_face_to_face (t, root_face);
}

void
t8_default_scheme_tet::element_get_first_descendant_face (const t8_element_t *elem, int face, t8_element_t *first_desc,
                                                          int level) const
{
  int corner;
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  /* Compute the first corner of this face */
  corner = t8_dtet_face_corner[face][0];
  /* Compute the descendant in this corner */
  t8_dtet_corner_descendant ((const t8_dtet_t *) elem, (t8_dtet_t *) first_desc, corner, level);
}

void
t8_default_scheme_tet::element_get_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                                         int level) const
{
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  /* Compute the last corner of this face */
  const int corner = SC_MAX (t8_dtet_face_corner[face][1], t8_dtet_face_corner[face][2]);
  /* Compute the descendant in this corner */
  t8_dtet_corner_descendant ((const t8_dtet_t *) elem, (t8_dtet_t *) last_desc, corner, level);
}

/* Construct the boundary element at a specific face. */
/* This function is defined here instead of in t8_dtet_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtet_bits.c. */
void
t8_default_scheme_tet::element_get_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                  [[maybe_unused]] const t8_scheme *scheme) const
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_dtri_t *b = (t8_dtri_t *) boundary;
  int face_cat;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (scheme->element_is_valid (T8_ECLASS_TRIANGLE, boundary));
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  /* The level of the boundary element is the same as the quadrant's level */
  b->level = t->level;
  /*
   * Depending on t's type and face, b's coordinates and type are defined
   * through t's coordinates. The faces can be divided into 3 categories.
   * category 1: b.x = t.z b.y = t.x - t.y
   * category 2: b.x = t.x b.y = t.z
   * category 3: b.x = t.x b.y = t.y
   *
   * This coordinate then has to be rescaled, since a tet has side length
   * T8_DTET_ROOT_LEN and a triangle T8_DTRI_ROOT_LEN
   *
   */
  face_cat = t8_dtet_type_face_to_boundary[t->type][face][0];
  b->type = t8_dtet_type_face_to_boundary[t->type][face][1];
  T8_ASSERT (face_cat == 1 || face_cat == 2);
  T8_ASSERT (b->type == 0 || b->type == 1);
  switch (face_cat) {
  case 1:
    b->x = t->z * T8_DTRI_ROOT_BY_DTET_ROOT;
    b->y = t->y * T8_DTRI_ROOT_BY_DTET_ROOT;
    break;
  case 2:
    b->x = t->x * T8_DTRI_ROOT_BY_DTET_ROOT;
    b->y = t->z * T8_DTRI_ROOT_BY_DTET_ROOT;
    break;
  case 3:
    b->x = t->x * T8_DTRI_ROOT_BY_DTET_ROOT;
    b->y = t->y * T8_DTRI_ROOT_BY_DTET_ROOT;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

int
t8_default_scheme_tet::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  const t8_dtet_t *t = (const t8_dtet_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  return t8_dtet_is_root_boundary (t, face);
}

int
t8_default_scheme_tet::element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                         int *neigh_face) const
{
  const t8_dtet_t *t = (const t8_dtet_t *) elem;
  t8_dtet_t *n = (t8_dtet_t *) neigh;

  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (neigh_face != NULL);
  *neigh_face = t8_dtet_face_neighbour (t, face, n);
  /* return true if neigh is inside the root */
  return t8_dtet_is_inside_root (n);
}

void
t8_default_scheme_tet::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  T8_ASSERT (id < ((t8_linearidx_t) 1) << 3 * level);
  T8_ASSERT (element_is_valid (elem));

  t8_dtet_init_linear_id ((t8_default_tet_t *) elem, id, level);
}

t8_linearidx_t
t8_default_scheme_tet::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  return t8_dtet_linear_id ((t8_default_tet_t *) elem, level);
}

void
t8_default_scheme_tet::element_construct_successor (const t8_element_t *elem1, t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  T8_ASSERT (0 <= element_get_level (elem1) && element_get_level (elem1) <= T8_DTET_MAXLEVEL);

  t8_dtet_successor ((const t8_default_tet_t *) elem1, (t8_default_tet_t *) elem2, element_get_level (elem1));
}

void
t8_default_scheme_tet::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  t8_dtet_first_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc, level);
}

void
t8_default_scheme_tet::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  t8_dtet_last_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc, level);
}

void
t8_default_scheme_tet::element_get_anchor (const t8_element_t *elem, int anchor[3]) const
{
  t8_dtet_t *tet = (t8_dtet_t *) elem;
  T8_ASSERT (element_is_valid (elem));

  anchor[0] = tet->x;
  anchor[1] = tet->y;
  anchor[2] = tet->z;
}

void
t8_default_scheme_tet::element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtet_compute_integer_coords ((const t8_default_tet_t *) elem, vertex, coords);
}

void
t8_default_scheme_tet::element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                            double coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtet_compute_vertex_ref_coords ((const t8_default_tet_t *) elem, vertex, coords);
}

void
t8_default_scheme_tet::element_get_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                     const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dtet_compute_reference_coords ((const t8_dtet_t *) elem, ref_coords, num_coords, out_coords);
}

/** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
 * Returns false otherwise.
 */
int
t8_default_scheme_tet::refines_irregular () const
{
  /*Tets refine regularly */
  return 0;
}

#if T8_ENABLE_DEBUG
int
t8_default_scheme_tet::element_is_valid (const t8_element_t *element) const

{
  return t8_dtet_is_valid ((const t8_dtet_t *) element);
}

void
t8_default_scheme_tet::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dtet_t *tet = (t8_dtet_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, z: %i, type: %i, level: %i", tet->x, tet->y, tet->z, tet->type,
            tet->level);
}
#endif

void
t8_default_scheme_tet::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a tet */
  t8_default_scheme_common::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#if T8_ENABLE_DEBUG
  {
    for (int i = 0; i < length; i++) {
      set_to_root (elem[i]);
    }
  }
#endif
}

void
t8_default_scheme_tet::element_init ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
{
#if T8_ENABLE_DEBUG
  t8_dtet_t *tets = (t8_dtet_t *) elem;
  for (int i = 0; i < length; i++) {
    t8_dtet_init (tets + i);
  }
#endif
}

void
t8_default_scheme_tet::set_to_root (t8_element_t *elem) const
{
  t8_dtet_t *tet = (t8_dtet_t *) elem;
  tet->level = 0;
  tet->x = 0;
  tet->y = 0;
  tet->z = 0;
  tet->type = 0;
}
/* use macro tri functionality */
void
t8_default_scheme_tet::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                         const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  t8_dtet_element_pack ((t8_dtet_t **) elements, count, send_buffer, buffer_size, position, comm);
}

/* use macro tri functionality */
void
t8_default_scheme_tet::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  t8_dtet_element_pack_size (count, comm, pack_size);
}

/* use macro tri functionality */
void
t8_default_scheme_tet::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position, t8_element_t **elements,
                                           const unsigned int count, sc_MPI_Comm comm) const
{
  t8_dtet_element_unpack (recvbuf, buffer_size, position, (t8_dtet_t **) elements, count, comm);
}

T8_EXTERN_C_END ();
