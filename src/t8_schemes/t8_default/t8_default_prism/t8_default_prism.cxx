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
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>
#include <t8_schemes/t8_scheme.hxx>

typedef t8_dprism_t t8_default_prism_t;

T8_EXTERN_C_BEGIN ();

size_t
t8_default_scheme_prism::get_element_size (void) const
{
  return sizeof (t8_dprism_t);
}

void
t8_default_scheme_prism::element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a tet */
  t8_default_scheme_common::element_new (length, elem);

  /* in debug mode, set sensible default values. */
#ifdef T8_ENABLE_DEBUG
  {
    for (int i = 0; i < length; i++) {
      set_to_root (elem[i]);
    }
  }
#endif
}

void
t8_default_scheme_prism::element_init ([[maybe_unused]] int length, [[maybe_unused]] t8_element_t *elem) const
{
#ifdef T8_ENABLE_DEBUG
  t8_dprism_t *prism = (t8_dprism_t *) elem;
  /* Set all values to 0 */
  for (int i = 0; i < length; i++) {
    t8_dprism_init_linear_id (prism + i, 0, 0);
    T8_ASSERT (t8_dprism_is_valid (prism + i));
  }
#endif
}

int
t8_default_scheme_prism::get_maxlevel (void) const
{
  return T8_DPRISM_MAXLEVEL;
}

int
t8_default_scheme_prism::element_get_level (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_get_level ((const t8_dprism_t *) elem);
}

t8_element_shape_t
t8_default_scheme_prism::element_get_face_shape ([[maybe_unused]] const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);

  return t8_dprism_face_shape (face);
}

void
t8_default_scheme_prism::element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  T8_ASSERT (element_is_valid (source));
  t8_dprism_copy ((const t8_dprism_t *) source, (t8_dprism_t *) dest);
  T8_ASSERT (element_is_valid (dest));
}

int
t8_default_scheme_prism::element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  return t8_dprism_compare ((const t8_dprism_t *) elem1, (const t8_dprism_t *) elem2);
}

int
t8_default_scheme_prism::element_is_equal (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  return t8_dprism_equal ((const t8_dprism_t *) elem1, (const t8_dprism_t *) elem2);
}

void
t8_default_scheme_prism::element_get_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_parent ((const t8_dprism_t *) elem, (t8_dprism_t *) parent);
  T8_ASSERT (element_is_valid (parent));
}

int
t8_default_scheme_prism::element_get_num_children ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DPRISM_CHILDREN;
}

int
t8_default_scheme_prism::element_get_num_face_children ([[maybe_unused]] const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_num_face_children (face);
}

int
t8_default_scheme_prism::element_get_face_corner ([[maybe_unused]] const t8_element_t *element, int face,
                                                  int corner) const
{
  T8_ASSERT (element_is_valid (element));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= corner && corner < T8_DPRISM_CORNERS);

  return t8_dprism_get_face_corner (face, corner);
}

int
t8_default_scheme_prism::element_get_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DPRISM_FACES;
}

int
t8_default_scheme_prism::element_get_child_id (const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_child_id ((const t8_dprism_t *) elem);
}

void
t8_default_scheme_prism::element_get_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_child ((const t8_dprism_t *) elem, childid, (t8_dprism_t *) child);
  T8_ASSERT (element_is_valid (child));
}

int
t8_default_scheme_prism::element_get_max_num_faces ([[maybe_unused]] const t8_element_t *elem) const
{
  T8_ASSERT (element_is_valid (elem));
  return T8_DPRISM_FACES;
}

void
t8_default_scheme_prism::element_get_children (const t8_element_t *elem, int length, t8_element_t *c[]) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (length == T8_DPRISM_CHILDREN);
  t8_dprism_childrenpv ((const t8_dprism_t *) elem, length, (t8_dprism_t **) c);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < length; i++) {
    T8_ASSERT (element_is_valid (c[i]));
  }
#endif
}

int
t8_default_scheme_prism::element_get_ancestor_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_ancestor_id ((t8_dprism_t *) elem, level);
}

void
t8_default_scheme_prism::element_get_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                       int num_children, int *child_indices) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (num_children == t8_dprism_num_face_children (face));
  t8_dprism_children_at_face ((const t8_dprism_t *) elem, face, (t8_dprism_t **) children, num_children, child_indices);
#if T8_ENABLE_DEBUG
  for (int i = 0; i < num_children; i++) {
    T8_ASSERT (element_is_valid (children[i]));
  }
#endif
}

int
t8_default_scheme_prism::element_face_get_child_face ([[maybe_unused]] const t8_element_t *elem, int face,
                                                      [[maybe_unused]] int face_child) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (face_child < t8_dprism_num_face_children (face));
  return t8_dprism_face_child_face (face);
}

int
t8_default_scheme_prism::element_face_get_parent_face (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  return t8_dprism_face_parent_face ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism::element_get_tree_face ([[maybe_unused]] const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  return t8_dprism_tree_face (face);
}

int
t8_default_scheme_prism::element_extrude_face (const t8_element_t *face, t8_element_t *elem, int root_face,
                                               [[maybe_unused]] const t8_scheme *scheme) const
{
#if T8_ENABLE_DEBUG
  const t8_eclass_t face_eclass = (t8_eclass_t) t8_eclass_face_types[T8_ECLASS_PRISM][root_face];
#endif
  T8_ASSERT (scheme->element_is_valid (face_eclass, face));
  T8_ASSERT (0 <= root_face && root_face < T8_DPRISM_FACES);
  t8_dprism_extrude_face (face, elem, root_face);
  T8_ASSERT (element_is_valid (elem));
  /* For the quad-faces of prisms it holds that only the corner-children touch the faces of the parent and
   * their face-numbers coincide. 
   * for the triangular-faces (bottom and top) the faces always have the same number and we can return the
   * root face-number as well. */
  return root_face;
}

int
t8_default_scheme_prism::elements_are_family (t8_element_t *const *fam) const
{
#ifdef T8_ENABLE_DEBUG
  for (int i = 0; i < T8_DPRISM_CHILDREN; i++) {
    T8_ASSERT (element_is_valid (fam[i]));
  }
#endif
  return t8_dprism_is_familypv ((t8_dprism_t **) fam);
}

void
t8_default_scheme_prism::element_get_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  T8_ASSERT (element_is_valid (elem1));
  T8_ASSERT (element_is_valid (elem2));
  const t8_default_prism_t *p1 = (const t8_default_prism_t *) elem1;
  const t8_default_prism_t *p2 = (const t8_default_prism_t *) elem2;
  t8_default_prism_t *c = (t8_default_prism_t *) nca;

  t8_dprism_nearest_common_ancestor (p1, p2, c);
  T8_ASSERT (t8_dprism_is_valid (c));
}

void
t8_default_scheme_prism::element_get_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                    [[maybe_unused]] const t8_scheme *scheme) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  t8_dprism_boundary_face ((const t8_dprism_t *) elem, face, boundary);
  T8_ASSERT (
    scheme->element_is_valid (static_cast<t8_eclass_t> (t8_eclass_face_types[T8_ECLASS_PRISM][face]), boundary));
}

const int t8_dprism_face_corner[5][4] = {
  { 1, 2, 4, 5 },
  { 0, 2, 3, 5 },
  { 0, 1, 3, 4 },
  { 0, 1, 2, -1 }, /*this face is a triangle -> -1 for the 4th corner */
  { 3, 4, 5, -1 }  /*this face is a triangle -> -1 for the 4th corner */
};

void
t8_default_scheme_prism::element_get_first_descendant_face (const t8_element_t *elem, int face,
                                                            t8_element_t *first_desc, int level) const
{
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (element_is_valid (elem));
  const int corner = t8_dprism_face_corner[face][0];
  t8_dprism_corner_descendant ((const t8_dprism_t *) elem, (t8_dprism_t *) first_desc, corner, level);
  T8_ASSERT (element_is_valid (first_desc));
}

void
t8_default_scheme_prism::element_get_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                                           int level) const
{
  int corner;
  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (element_is_valid (elem));
  if (face < 3) {
    corner = t8_dprism_face_corner[face][3];
  }
  else {
    corner = t8_dprism_face_corner[face][2];
  }
  t8_dprism_corner_descendant ((const t8_dprism_t *) elem, (t8_dprism_t *) last_desc, corner, level);
  T8_ASSERT (element_is_valid (last_desc));
}

int
t8_default_scheme_prism::element_is_root_boundary (const t8_element_t *elem, int face) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_is_root_boundary ((const t8_dprism_t *) elem, face);
}

int
t8_default_scheme_prism::element_get_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                           int *neigh_face) const
{
  const t8_dprism_t *p = (const t8_dprism_t *) elem;
  t8_dprism_t *n = (t8_dprism_t *) neigh;

  T8_ASSERT (0 <= face && face < T8_DPRISM_FACES);
  T8_ASSERT (neigh_face != NULL);
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (element_is_valid (neigh));

  *neigh_face = t8_dprism_face_neighbour (p, face, n);
  /* return true if neigh is inside the root */
  return t8_dprism_is_inside_root (n);
}

void
t8_default_scheme_prism::element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (id < ((t8_linearidx_t) 1) << 3 * level);

  t8_dprism_init_linear_id ((t8_default_prism_t *) elem, level, id);

  T8_ASSERT (element_is_valid (elem));
}

void
t8_default_scheme_prism::element_construct_successor (const t8_element_t *elem, t8_element_t *s) const
{
  T8_ASSERT (1 <= element_get_level (elem) && element_get_level (elem) <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (element_is_valid (elem));

  t8_dprism_successor ((const t8_default_prism_t *) elem, (t8_default_prism_t *) s, element_get_level (elem));
  T8_ASSERT (element_is_valid (s));
}

void
t8_default_scheme_prism::element_get_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_first_descendant ((const t8_default_prism_t *) elem, (t8_default_prism_t *) desc, level);
  T8_ASSERT (element_is_valid (desc));
}

void
t8_default_scheme_prism::element_get_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  T8_ASSERT (0 <= level && level <= T8_DPRISM_MAXLEVEL);
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_last_descendant ((const t8_default_prism_t *) elem, (t8_default_prism_t *) desc, level);
  T8_ASSERT (element_is_valid (desc));
}

void
t8_default_scheme_prism::element_get_anchor (const t8_element_t *elem, int anchor[3]) const
{
  t8_dprism_t *prism = (t8_dprism_t *) elem;
  T8_ASSERT (element_is_valid (elem));
  anchor[0] = prism->tri.x / T8_DTRI_ROOT_LEN * T8_DPRISM_ROOT_LEN;
  anchor[1] = prism->tri.y / T8_DTRI_ROOT_LEN * T8_DPRISM_ROOT_LEN;
  anchor[2] = prism->line.x / T8_DLINE_ROOT_LEN * T8_DPRISM_ROOT_LEN;
}

void
t8_default_scheme_prism::element_get_vertex_integer_coords (const t8_element_t *elem, int vertex, int coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_vertex_integer_coords ((const t8_dprism_t *) elem, vertex, coords);
}

void
t8_default_scheme_prism::element_get_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                              double coords[]) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_vertex_ref_coords ((const t8_dprism_t *) elem, vertex, coords);
}

void
t8_default_scheme_prism::element_get_reference_coords (const t8_element_t *elem, const double *ref_coords,
                                                       const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (element_is_valid (elem));
  t8_dprism_compute_reference_coords ((const t8_dprism_t *) elem, ref_coords, num_coords, out_coords);
}

t8_linearidx_t
t8_default_scheme_prism::element_get_linear_id (const t8_element_t *elem, int level) const
{
  T8_ASSERT (element_is_valid (elem));
  return t8_dprism_linear_id ((const t8_dprism_t *) elem, level);
}

int
t8_default_scheme_prism::refines_irregular (void) const
{
  /*prisms refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG

int
t8_default_scheme_prism::element_is_valid (const t8_element_t *elem) const
{
  T8_ASSERT (elem != NULL);
  return t8_dprism_is_valid ((const t8_dprism_t *) elem);
}

void
t8_default_scheme_prism::element_to_string (const t8_element_t *elem, char *debug_string, const int string_size) const
{
  T8_ASSERT (element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_dprism_t *prism = (t8_dprism_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, z: %i, type: %i, level: %i", prism->tri.x, prism->tri.y,
            prism->line.x, prism->tri.type, prism->tri.level);
}

#endif /* T8_ENABLE_DEBUG */

void
t8_default_scheme_prism::set_to_root (t8_element_t *elem) const
{
  t8_dprism_t *prism = (t8_dprism_t *) elem;
  prism->line.level = 0;
  prism->line.x = 0;
  prism->tri.x = 0;
  prism->tri.y = 0;
  prism->tri.type = 0;
  prism->tri.level = 0;
}
/* each prism is packed as x (line.x), y (tri.x), z(tri.y) coordinates, type and the level */
void
t8_default_scheme_prism::element_MPI_Pack (t8_element_t **const elements, const unsigned int count, void *send_buffer,
                                           const int buffer_size, int *position, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_default_prism_t **prisms = (t8_default_prism_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(prisms[ielem]->line.x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&prisms[ielem]->tri.x, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&prisms[ielem]->tri.y, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&prisms[ielem]->tri.type, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);

    T8_ASSERT (prisms[ielem]->line.level == prisms[ielem]->tri.level);
    mpiret = sc_MPI_Pack (&prisms[ielem]->line.level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
}

/* each prism is packed as x (line.x), y (tri.x), z(tri.y) coordinates, type and the level */
void
t8_default_scheme_prism::element_MPI_Pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  /*x,y,z*/
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += 3 * datasize;

  /*type, level*/
  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += 2 * datasize;

  *pack_size = count * singlesize;
}

/* each prism is packed as x (line.x), y (tri.x), z(tri.y) coordinates, type and the level */
void
t8_default_scheme_prism::element_MPI_Unpack (void *recvbuf, const int buffer_size, int *position,
                                             t8_element_t **elements, const unsigned int count, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_default_prism_t **prisms = (t8_default_prism_t **) elements;
  for (unsigned int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(prisms[ielem]->line.x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(prisms[ielem]->tri.x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(prisms[ielem]->tri.y), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(prisms[ielem]->tri.type), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(prisms[ielem]->tri.level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
    prisms[ielem]->line.level = prisms[ielem]->tri.level;
  }
}

T8_EXTERN_C_END ();
