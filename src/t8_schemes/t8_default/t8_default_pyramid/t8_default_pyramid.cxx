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
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_schemes/t8_scheme.hxx>

typedef t8_dpyramid_t t8_default_pyramid_t;

T8_EXTERN_C_BEGIN ();

size_t
t8_default_scheme_pyramid::get_element_size (void) const
{
  return sizeof (t8_dpyramid_t);
}

void
t8_default_scheme_pyramid::element_new (int length, t8_element_t **elem) const
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
t8_default_scheme_pyramid::element_init ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
{
#if T8_ENABLE_DEBUG
  t8_dpyramid_t *pyramid = (t8_dpyramid_t *) elem;
  /* Set all values to 0 */
  for (int i = 0; i < length; i++) {
    t8_dpyramid_init_linear_id (pyramid + i, 0, 0);
    T8_ASSERT (t8_dpyramid_is_valid (pyramid + i));
  }
  pyramid->pyramid.type = 6;
#endif
}

int
t8_default_scheme_pyramid::get_maxlevel (void) const
{
  return T8_DPYRAMID_MAXLEVEL;
}

int
t8_default_scheme_pyramid::element_get_max_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_max_num_faces ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= T8_DPYRAMID_MAXLEVEL);
  return t8_dpyramid_ancestor_id ((const t8_dpyramid_t *) elem, level);
}

int
t8_default_scheme_pyramid::element_get_num_children (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_num_children ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid::element_get_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_num_corners ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid::element_get_num_siblings (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_num_siblings ((const t8_dpyramid_t *) elem);
}

t8_gloidx_t
t8_default_scheme_pyramid::count_leaves_from_root (int level) const
{
  return 2 * sc_intpow64u (8, level) - sc_intpow64u (6, level);
}

int
t8_default_scheme_pyramid::element_get_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_num_faces ((const t8_dpyramid_t *) elem);
}

int
t8_default_scheme_pyramid::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  return t8_dpyramid_compare ((const t8_dpyramid_t *) elem1, (const t8_dpyramid_t *) elem2);
}

int
t8_default_scheme_pyramid::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dpyramid_equal ((const t8_dpyramid_t *) elem1, (const t8_dpyramid_t *) elem2);
}

void
t8_default_scheme_pyramid::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  t8_dpyramid_copy ((const t8_dpyramid_t *) source, (t8_dpyramid_t *) dest);
  T8_ASSERT (element_is_valid (dest));
}

void
t8_default_scheme_pyramid::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= childid && childid < t8_dpyramid_num_children ((const t8_dpyramid_t *) elem));
  t8_dpyramid_child ((t8_dpyramid_t *) elem, childid, (t8_dpyramid_t *) child);
  T8_ASSERT (element_is_valid (child));
}

void
t8_default_scheme_pyramid::element_get_children (const t8_element_t *elem, [[maybe_unused]] int length,
                                                 t8_element_t *c[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_children ((const t8_dpyramid_t *) elem, (t8_dpyramid_t **) c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (element_is_valid (c[i]));
  }
#endif
}

void
t8_default_scheme_pyramid::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                         int num_children, int *child_indices) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_children_at_face ((const t8_dpyramid_t *) elem, face, (t8_dpyramid_t **) children, num_children,
                                child_indices);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < num_children; i++) {
    T8_ASSERT (element_is_valid (children[i]));
  }
#endif
}

int
t8_default_scheme_pyramid::element_face_get_child_face (const t8_element_t *elem, int face, int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_face_child_face ((const t8_dpyramid_t *) elem, face, face_child);
}

t8_element_shape_t
t8_default_scheme_pyramid::element_get_face_shape (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_face_shape ((const t8_dpyramid_t *) elem, face);
}

int
t8_default_scheme_pyramid::element_get_child_id (const t8_element_t *p) const
{
  T8_ASSERT (element_is_valid (p));
  return t8_dpyramid_child_id ((const t8_dpyramid_t *) p);
}

int
t8_default_scheme_pyramid::element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                             int *neigh_face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_face_neighbor_inside ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) neigh, face, neigh_face);
}

int
t8_default_scheme_pyramid::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_face_parent_face ((const t8_dpyramid_t *) elem, face);
}

void
t8_default_scheme_pyramid::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_first_descendant ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) desc, level);
  T8_ASSERT (element_is_valid (desc));
}

void
t8_default_scheme_pyramid::element_get_first_descendant_face (const t8_element_t *elem, int face,
                                                              t8_element_t *first_desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_first_descendant_face ((const t8_dpyramid_t *) elem, face, (t8_dpyramid_t *) first_desc, level);
  T8_ASSERT (element_is_valid (first_desc));
}

int
t8_default_scheme_pyramid::element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  T8_ASSERT (element_is_valid (element));
  return t8_dpyramid_get_face_corner ((const t8_dpyramid_t *) element, face, corner);
}

int
t8_default_scheme_pyramid::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_get_level ((const t8_dpyramid_t *) elem);
}

void
t8_default_scheme_pyramid::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  t8_dpyramid_init_linear_id ((t8_dpyramid_t *) elem, level, id);
  T8_ASSERT (element_is_valid (elem));
}

int
t8_default_scheme_pyramid::elements_are_family (t8_element_t *const *fam) const
{
#if T8_ENABLE_DEBUG
  const int num_siblings = element_get_num_siblings (fam[0]);
  for (int i = 0; i < num_siblings; i++) {
    T8_ASSERT (element_is_valid (fam[i]));
  }
#endif
  return t8_dpyramid_is_family ((t8_dpyramid_t **) fam);
}

int
t8_default_scheme_pyramid::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_is_root_boundary ((const t8_dpyramid_t *) elem, face);
}

void
t8_default_scheme_pyramid::element_get_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                      [[maybe_unused]] const t8_scheme *scheme) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_boundary_face ((const t8_dpyramid_t *) elem, face, boundary);
  T8_ASSERT (
    scheme->element_is_valid (static_cast<t8_eclass_t> (t8_eclass_face_types[T8_ECLASS_PYRAMID][face]), boundary));
}

int
t8_default_scheme_pyramid::element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face,
                                                 [[maybe_unused]] const t8_scheme *scheme) const
{
#if T8_ENABLE_DEBUG
  const t8_eclass_t face_eclass = (t8_eclass_t) t8_eclass_face_types[T8_ECLASS_PYRAMID][root_face];
  T8_ASSERT (scheme->element_is_valid (face_eclass, face));
#endif
  return t8_dpyramid_extrude_face (face, (t8_dpyramid_t *) elem, root_face);
  T8_ASSERT (element_is_valid (elem));
}

t8_element_shape_t
t8_default_scheme_pyramid::element_get_shape (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_shape ((const t8_dpyramid_t *) elem);
}

t8_gloidx_t
t8_default_scheme_pyramid::element_count_leaves (const t8_element_t *element, int level) const
{
  const int element_level = element_get_level (element);
  const int dim = get_dimension ();
  const t8_element_shape_t element_shape = element_get_shape (element);
  if (element_shape == T8_ECLASS_PYRAMID) {
    const int level_diff = level - element_level;
    return element_level > level ? 0 : 2 * sc_intpow64 (8, level_diff) - sc_intpow64 (6, level_diff);
  }
  return count_leaves_from_level (element_level, level, dim);
}

int
t8_default_scheme_pyramid::element_get_tree_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_tree_face ((const t8_dpyramid_t *) elem, face);
}

int
t8_default_scheme_pyramid::element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem,
                                                          [[maybe_unused]] int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DPYRAMID_FACE_CHILDREN;
}

t8_linearidx_t
t8_default_scheme_pyramid::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dpyramid_linear_id ((const t8_dpyramid_t *) elem, level);
}

void
t8_default_scheme_pyramid::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_last_descendant ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) desc, level);
  T8_ASSERT (element_is_valid (desc));
}

void
t8_default_scheme_pyramid::element_get_last_descendant_face (const t8_element_t *elem, int face,
                                                             t8_element_t *last_desc, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_last_descendant_face ((const t8_dpyramid_t *) elem, face, (t8_dpyramid_t *) last_desc, level);
  T8_ASSERT (element_is_valid (last_desc));
}

void
t8_default_scheme_pyramid::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_parent ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) parent);
  T8_ASSERT (element_is_valid (parent));
}

void
t8_default_scheme_pyramid::element_construct_successor (const t8_element_t *elem, t8_element_t *s) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dpyramid_successor ((const t8_dpyramid_t *) elem, (t8_dpyramid_t *) s, element_get_level (elem));
  T8_ASSERT (element_is_valid (s));
}

void
t8_default_scheme_pyramid::element_get_anchor (const t8_element_t *elem, int anchor[3]) const
{
  t8_dpyramid_t *pyra = (t8_dpyramid_t *) elem;

  T8_ASSERT (element_is_valid (elem));
  anchor[0] = pyra->pyramid.x;
  anchor[1] = pyra->pyramid.y;
  anchor[2] = pyra->pyramid.z;
}

void
t8_default_scheme_pyramid::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2,
                                            t8_element_t *nca) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));

  t8_dpyramid_nearest_common_ancestor ((const t8_dpyramid_t *) elem1, (const t8_dpyramid_t *) elem2,
                                       (t8_dpyramid_t *) nca);
  T8_ASSERT (element_is_valid (nca));
}

int
t8_default_scheme_pyramid::refines_irregular () const
{
  /*Pyramids do not refine regularly */
  return true;
}

#if T8_ENABLE_DEBUG
void
t8_default_scheme_pyramid::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dpyramid_t *pyra = (t8_dpyramid_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, z: %i, type: %i, level: %i, switch_shape_at_level: %i",
            pyra->pyramid.x, pyra->pyramid.y, pyra->pyramid.x, pyra->pyramid.type, pyra->pyramid.level,
            pyra->switch_shape_at_level);
}
#endif

void
t8_default_scheme_pyramid::set_to_root (t8_element_t *elem) const
{
  t8_dpyramid_t *pyramid = (t8_dpyramid_t *) elem;
  pyramid->pyramid.level = 0;
  pyramid->pyramid.x = 0;
  pyramid->pyramid.y = 0;
  pyramid->pyramid.z = 0;
  pyramid->pyramid.type = T8_DPYRAMID_ROOT_TYPE;
  pyramid->switch_shape_at_level = -1;
}
/* each pyramid is packed as a tet and the switch_shape_at_level marker */
void
t8_default_scheme_pyramid::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                             const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  t8_default_pyramid_t **pyramids = (t8_default_pyramid_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    t8_dtet_t *p = &pyramids[ielem]->pyramid;
    t8_dtet_element_pack (&p, 1, send_buffer, buffer_size, position, comm);

    SC_CHECK_MPI (sc_MPI_Pack (&pyramids[ielem]->switch_shape_at_level, 1, sc_MPI_INT8_T, send_buffer, buffer_size,
                               position, comm));
  }
}

/* each pyramid is packed as a tet and the switch_shape_at_level marker */
void
t8_default_scheme_pyramid::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;

  t8_dtet_element_pack_size (1, comm, &datasize);
  singlesize += datasize;

  /* switch shape at level */
  SC_CHECK_MPI (sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize));
  singlesize += datasize;

  *pack_size = count * singlesize;
}

/* each pyramid is packed as a tet and the switch_shape_at_level marker */
void
t8_default_scheme_pyramid::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                               t8_element_t **elements, const unsigned int count,
                                               sc_MPI_Comm comm) const
{
  int mpiret;
  t8_default_pyramid_t **pyramids = (t8_default_pyramid_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    t8_dtet *p = &pyramids[ielem]->pyramid;
    t8_dtet_element_unpack (recvbuf, buffer_size, position, &p, 1, comm);

    mpiret
      = sc_MPI_Unpack (recvbuf, buffer_size, position, &pyramids[ielem]->switch_shape_at_level, 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
}

T8_EXTERN_C_END ();
