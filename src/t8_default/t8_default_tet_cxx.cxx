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

#include "t8_default_common_cxx.hxx"
#include "t8_default_tet_cxx.hxx"
#include "t8_dtet_bits.h"
#include "t8_dtri_bits.h"
#include "t8_dtet_connectivity.h"

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

typedef t8_dtet_t   t8_default_tet_t;

int
t8_default_scheme_tet_c::t8_element_maxlevel (void)
{
  return T8_DTET_MAXLEVEL;
}

int
t8_default_scheme_tet_c::t8_element_level (const t8_element_t * elem)
{
  return t8_dtet_get_level ((t8_dtet_t *) elem);
}

void
t8_default_scheme_tet_c::t8_element_copy (const t8_element_t * source,
                                          t8_element_t * dest)
{
  t8_dtet_copy ((const t8_dtet_t *) source, (t8_dtet_t *) dest);
}

int
t8_default_scheme_tet_c::t8_element_compare (const t8_element_t * elem1,
                                             const t8_element_t * elem2)
{
  int                 maxlvl;
  u_int64_t           id1, id2;

  /* Compute the bigger level of the two */
  maxlvl = SC_MAX (t8_default_scheme_tet_c::t8_element_level (elem1),
                   t8_default_scheme_tet_c::t8_element_level (elem2));
  /* Compute the linear ids of the elements */
  id1 = t8_default_scheme_tet_c::t8_element_get_linear_id (elem1, maxlvl);
  id2 = t8_default_scheme_tet_c::t8_element_get_linear_id (elem2, maxlvl);
  /* return negativ if id1 < id2, zero if id1 = id2, positive if id1 > id2 */
  return id1 < id2 ? -1 : id1 != id2;
}

void
t8_default_scheme_tet_c::t8_element_parent (const t8_element_t * elem,
                                            t8_element_t * parent)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *p = (t8_default_tet_t *) parent;

  t8_dtet_parent (t, p);
}

void
t8_default_scheme_tet_c::t8_element_sibling (const t8_element_t * elem,
                                             int sibid,
                                             t8_element_t * sibling)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *s = (t8_default_tet_t *) sibling;

  t8_dtet_sibling (t, sibid, s);
}

int
t8_default_scheme_tet_c::t8_element_num_faces (const t8_element_t * elem)
{
  return T8_DTET_FACES;
}

int
t8_default_scheme_tet_c::t8_element_num_children (const t8_element_t * elem)
{
  return T8_DTET_CHILDREN;
}

int
t8_default_scheme_tet_c::t8_element_num_face_children (const t8_element_t *
                                                       elem, int face)
{
  return T8_DTET_FACE_CHILDREN;
}

void
t8_default_scheme_tet_c::t8_element_child (const t8_element_t * elem,
                                           int childid, t8_element_t * child)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_default_tet_t   *c = (t8_default_tet_t *) child;

  t8_dtet_child (t, childid, c);
}

void
t8_default_scheme_tet_c::t8_element_children (const t8_element_t * elem,
                                              int length, t8_element_t * c[])
{
  T8_ASSERT (length == T8_DTET_CHILDREN);

  t8_dtet_childrenpv ((const t8_dtet_t *) elem, (t8_dtet_t **) c);
}

int
t8_default_scheme_tet_c::t8_element_child_id (const t8_element_t * elem)
{
  return t8_dtet_child_id ((t8_dtet_t *) elem);
}

int
t8_default_scheme_tet_c::t8_element_is_family (t8_element_t ** fam)
{
  return t8_dtet_is_familypv ((const t8_dtet_t **) fam);
}

void
t8_default_scheme_tet_c::t8_element_nca (const t8_element_t * elem1,
                                         const t8_element_t * elem2,
                                         t8_element_t * nca)
{
  const t8_default_tet_t *t1 = (const t8_default_tet_t *) elem1;
  const t8_default_tet_t *t2 = (const t8_default_tet_t *) elem2;
  t8_default_tet_t   *c = (t8_default_tet_t *) nca;

  t8_dtet_nearest_common_ancestor (t1, t2, c);
}

t8_eclass_t
  t8_default_scheme_tet_c::t8_element_face_class (const t8_element_t * elem,
                                                  int face)
{
  return T8_ECLASS_TRIANGLE;
}

void
t8_default_scheme_tet_c::t8_element_children_at_face (const t8_element_t *
                                                      elem, int face,
                                                      t8_element_t *
                                                      children[],
                                                      int num_children)
{
  const t8_dtet_t    *t = (const t8_dtet_t *) elem;
  t8_dtet_t         **c = (t8_dtet_t **) children;
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (num_children == T8_DTET_FACE_CHILDREN);

  t8_dtet_children_at_face (t, face, c, num_children);
}

int
t8_default_scheme_tet_c::t8_element_face_child_face (const t8_element_t *
                                                     elem, int face,
                                                     int face_child)
{
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  T8_ASSERT (0 <= face && face < T8_DTET_FACE_CHILDREN);
  return t8_dtet_face_child_face ((const t8_dtet_t *) elem, face, face_child);
}

int
t8_default_scheme_tet_c::t8_element_tree_face (const t8_element_t * elem,
                                               int face)
{
  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  return t8_dtet_tree_face ((t8_dtet_t *) elem, face);
}

/* Construct the inner element from a boundary element. */
/* This function is defined here instead of in t8_dri_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtri_bits.c. This would be needed by an implementation, at least
 * for tets. */
void
t8_default_scheme_tet_c::t8_element_extrude_face (const t8_element_t * face,
                                                  t8_element_t * elem,
                                                  int root_face)
{
  const t8_dtri_t    *b = (const t8_dtri_t *) face;
  t8_dtet_t          *t = (t8_dtet_t *) elem;

  T8_ASSERT (0 <= root_face && root_face < T8_DTET_FACES);
  t->level = b->level;
  t->eclass = T8_ECLASS_TET;
  switch (root_face) {
    /* Since the root triangle may have a different scale then the
     * root tetrahedron, we have to rescale the coordinates. */
  case 0:
    t->type = b->type == 0 ? 0 : 1;
    t->x = T8_DTET_ROOT_LEN - T8_DTET_LEN (t->level);
    t->y = ((int64_t) b->x * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
    t->z = ((int64_t) b->y * T8_DTET_ROOT_LEN) / T8_DTRI_ROOT_LEN;
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
}

/* Construct the boundary element at a specific face. */
/* This function is defined here instead of in t8_dtet_bits.c since
 * the compile logic does not allow for t8_dtri_t and t8_dtet_t to exist
 * both in t8_dtet_bits.c. */
void
t8_default_scheme_tet_c::t8_element_boundary_face (const t8_element_t * elem,
                                                   int face,
                                                   t8_element_t * boundary)
{
  const t8_default_tet_t *t = (const t8_default_tet_t *) elem;
  t8_dtri_t          *b = (t8_dtri_t *) boundary;
  int                 face_cat;

  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  /* The level of the boundary element is the same as the quadrant's level */
  b->level = t->level;
  /*
   * Depending on t's type and face, b's coordinates and type are defined
   * through t's coordinates. The faces can be divided into 3 category.
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

void
t8_default_scheme_tet_c::t8_element_boundary (const t8_element_t * elem,
                                              int min_dim, int length,
                                              t8_element_t ** boundary)
{
  int                 iface;

  T8_ASSERT (length == T8_DTET_FACES);
  for (iface = 0; iface < T8_DTET_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
}

int
t8_default_scheme_tet_c::t8_element_is_root_boundary (const t8_element_t *
                                                      elem, int face)
{
  const t8_dtet_t    *t = (const t8_dtet_t *) elem;

  return t8_dtet_is_root_boundary (t, face);
}

int
t8_default_scheme_tet_c::t8_element_face_neighbor_inside (const t8_element_t *
                                                          elem,
                                                          t8_element_t *
                                                          neigh, int face)
{
  const t8_dtet_t    *t = (const t8_dtet_t *) elem;
  t8_dtet_t          *n = (t8_dtet_t *) neigh;

  T8_ASSERT (0 <= face && face < T8_DTET_FACES);
  (void) t8_dtet_face_neighbour (t, face, n);
  /* return true if neigh is inside the root */
  return t8_dtet_is_inside_root (n);
}

void
t8_default_scheme_tet_c::t8_element_set_linear_id (t8_element_t * elem,
                                                   int level, uint64_t id)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);
  T8_ASSERT (0 <= id && id < ((uint64_t) 1) << 3 * level);

  t8_dtet_init_linear_id ((t8_default_tet_t *) elem, id, level);
}

uint64_t
  t8_default_scheme_tet_c::t8_element_get_linear_id (const t8_element_t *
                                                     elem, int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  return t8_dtet_linear_id ((t8_default_tet_t *) elem, level);
}

void
t8_default_scheme_tet_c::t8_element_successor (const t8_element_t * elem1,
                                               t8_element_t * elem2,
                                               int level)
{
  T8_ASSERT (0 <= level && level <= T8_DTET_MAXLEVEL);

  t8_dtet_successor ((const t8_default_tet_t *) elem1,
                     (t8_default_tet_t *) elem2, level);
}

void
t8_default_scheme_tet_c::t8_element_first_descendant (const t8_element_t *
                                                      elem,
                                                      t8_element_t * desc)
{
  t8_dtet_first_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc);
}

void
t8_default_scheme_tet_c::t8_element_last_descendant (const t8_element_t *
                                                     elem,
                                                     t8_element_t * desc)
{
  t8_dtet_last_descendant ((t8_dtet_t *) elem, (t8_dtet_t *) desc);
}

void
t8_default_scheme_tet_c::t8_element_anchor (const t8_element_t * elem,
                                            int anchor[3])
{
  t8_dtet_t          *tet = (t8_dtet_t *) elem;

  anchor[0] = tet->x;
  anchor[1] = tet->y;
  anchor[2] = tet->z;
}

int
t8_default_scheme_tet_c::t8_element_root_len (const t8_element_t * elem)
{
  return T8_DTET_ROOT_LEN;
}

void
t8_default_scheme_tet_c::t8_element_vertex_coords (const t8_element_t * t,
                                                   int vertex, int coords[])
{
  t8_dtet_compute_coords ((const t8_default_tet_t *) t, vertex, coords);
}

 /* Constructor */
t8_default_scheme_tet_c::t8_default_scheme_tet_c (void)
{
  eclass = T8_ECLASS_TET;
  element_size = sizeof (t8_dtet_t);
  ts_context = sc_mempool_new (sizeof (element_size));
}

 /* Destructor */
t8_default_scheme_tet_c::~t8_default_scheme_tet_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
