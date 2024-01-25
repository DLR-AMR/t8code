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
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

static int
t8_hilbert_compute_cubeid (const t8_hilbert_t *elem, int level)
{
  if (level == 0) {
    return 0;
  }
  int cube_id = 0;

  T8_ASSERT (0 <= elem->level && elem->level <= T8_HILBERT_MAXLEVEL);
  const t8_hilbert_coord_t h = T8_HILBERT_LEN (level);

  cube_id |= ((elem->x & h) ? 1 : 0);
  cube_id |= ((elem->y & h) ? 2 : 0);
  return cube_id;
}

static void
t8_hilbert_root (t8_hilbert_t *hilbert)
{
  hilbert->x = 0;
  hilbert->y = 0;
  hilbert->level = 0;
  hilbert->type = 0;
}

int
t8_default_scheme_quad_c::t8_element_maxlevel (void) const
{
  return T8_HILBERT_MAXLEVEL;
}

int
t8_default_scheme_quad_c::t8_element_level (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_hilbert_t *) elem)->level;
}

void
t8_default_scheme_quad_c::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  const t8_hilbert_t *source_elem = (const t8_hilbert_t *) source;
  t8_hilbert_t *dest_elem = (t8_hilbert_t *) dest;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (source_elem == dest_elem) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *dest_elem = *source_elem;
}

int
t8_default_scheme_quad_c::t8_element_compare (const t8_element_t *elem1, const t8_element_t *elem2) const
{
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  int maxlevel = SC_MAX (t8_element_level (elem1), t8_element_level (elem2));
  t8_linearidx_t id1 = t8_element_get_linear_id (elem1, maxlevel);
  t8_linearidx_t id2 = t8_element_get_linear_id (elem2, maxlevel);
  return id1 < id2 ? -1 : (id1 > id2 ? 1 : 0);
}

int
t8_default_scheme_quad_c::t8_element_equal (const t8_element_t *e1, const t8_element_t *e2) const
{
  const t8_hilbert_t *elem1 = (const t8_hilbert_t *) e1;
  const t8_hilbert_t *elem2 = (const t8_hilbert_t *) e2;
  return elem1->x == elem2->x && elem1->y == elem2->y && elem1->level == elem2->level && elem1->type == elem2->type;
}

void
t8_default_scheme_quad_c::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_hilbert_t *hilbert_elem = (const t8_hilbert_t *) elem;
  T8_ASSERT (hilbert_elem->level > 0);

  t8_hilbert_t *hilbert_parent = (t8_hilbert_t *) parent;

  const int8_t cube_id = t8_hilbert_compute_cubeid (hilbert_elem, hilbert_elem->level);
  hilbert_parent->type = t8_hilbert_cubeid_type_to_parenttype[cube_id][hilbert_elem->type];

  /* remove least significant coordinate bit at level  */
  const t8_hilbert_coord_t length = T8_HILBERT_LEN (hilbert_elem->level);
  hilbert_parent->x = hilbert_elem->x & ~length;
  hilbert_parent->y = hilbert_elem->y & ~length;

  hilbert_parent->level = hilbert_elem->level - 1;
  T8_ASSERT (hilbert_parent->level >= 0);
}

void
t8_default_scheme_quad_c::t8_element_sibling (const t8_element_t *elem, int sibid, t8_element_t *sibling) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));
  t8_element_parent (elem, sibling);
  t8_element_child (sibling, sibid, sibling);
}

t8_eclass_t
t8_default_scheme_quad_c::t8_element_child_eclass (int childid) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_num_children (const t8_element_t *elem) const
{
  return T8_HILBERT_CHILDREN;
}

void
t8_default_scheme_quad_c::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_hilbert_t *hilbert_elem = (const t8_hilbert_t *) elem;
  t8_hilbert_t *hilbert_child = (t8_hilbert_t *) child;

  T8_ASSERT (0 <= childid && childid < T8_HILBERT_CHILDREN);
  T8_ASSERT (0 <= hilbert_elem->level && hilbert_elem->level <= T8_HILBERT_MAXLEVEL);

  /* Compute the cube id and shift the coordinates accordingly */
  int8_t cube_id;
  cube_id = t8_hilbert_type_Iloc_to_childcubeid[hilbert_elem->type][childid];
  hilbert_child->type = t8_hilbert_type_Iloc_to_childtype[hilbert_elem->type][childid];

  const t8_hilbert_coord_t length = T8_HILBERT_LEN (hilbert_elem->level + 1);

  hilbert_child->x = hilbert_elem->x + ((cube_id & 1) ? length : 0);
  hilbert_child->y = hilbert_elem->y + ((cube_id & 2) ? length : 0);

  hilbert_child->level = hilbert_elem->level + 1;
}

void
t8_default_scheme_quad_c::t8_element_children (const t8_element_t *elem, int length, t8_element_t *children[]) const
{
  T8_ASSERT (length == t8_element_num_children (elem));
  for (int ichild = 0; ichild < length; ichild++) {
    t8_element_child (elem, ichild, children[ichild]);
  }
}

int
t8_default_scheme_quad_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_hilbert_t *hilbert = (const t8_hilbert_t *) elem;
  T8_ASSERT (hilbert->level >= 0);
  if (hilbert->level == 0) {
    return -1;
  }
  const int8_t cube_id = t8_hilbert_compute_cubeid (hilbert, hilbert->level);
  return t8_hilbert_type_cubeid_to_Iloc[hilbert->type][cube_id];
}

int
t8_default_scheme_quad_c::t8_element_is_family (t8_element_t **fam) const
{
  if (t8_element_level (fam[0]) == 0)
    return 0;
  t8_element_t *parent, *parent_compare;
  t8_element_new (1, &parent);
  t8_element_new (1, &parent_compare);
  t8_element_parent (fam[0], parent);
  int num_children = t8_element_num_children (parent);
  bool is_family = true;
  for (int ichild = 1; ichild < num_children; ichild++) {
    t8_element_parent (fam[ichild], parent_compare);
    if (!t8_element_equal (parent, parent_compare)) {
      is_family = false;
      break;
    }
  }
  t8_element_destroy (1, &parent);
  t8_element_destroy (1, &parent_compare);
  return is_family;
}

t8_linearidx_t
t8_hilbert_num_descendants_of_child_at_leveldiff (t8_element_t *elem, const int childid, const int leveldiff)
{
  return 1 << (leveldiff - 1) * T8_HILBERT_DIM;
}

t8_linearidx_t
t8_default_scheme_quad_c::t8_hilbert_linear_id_recursive (t8_element_t *elem, const t8_linearidx_t id,
                                                          const int level_diff) const
{
  if (t8_element_level (elem) == 0)
    return id;

  const int childid = t8_element_child_id (elem);
  t8_element_parent (elem, elem);

  t8_linearidx_t parent_id = 0;
  for (int ichild = 0; ichild < childid; ichild++) {
    t8_linearidx_t num_child_descendants = 1LL << (level_diff * T8_HILBERT_DIM);
    parent_id += num_child_descendants;
  }
  parent_id += id;
  t8_debugf ("lin %li\n", parent_id);
  return t8_hilbert_linear_id_recursive (elem, parent_id, level_diff + 1);
}

void
t8_default_scheme_quad_c::t8_hilbert_init_linear_id_recursive (t8_element_t *elem, const int level_diff,
                                                               t8_linearidx_t id) const
{
  T8_ASSERT (0 <= id);
  T8_ASSERT (1 <= level_diff && level_diff <= T8_HILBERT_MAXLEVEL);

  if (id == 0) {
    t8_element_first_descendant (elem, elem, t8_element_level (elem) + level_diff);
    return;
  }

  if (level_diff == 1) {
    T8_ASSERT (id <= T8_HILBERT_CHILDREN);
    t8_element_child (elem, id, elem);
    return;
  }

  t8_linearidx_t sum_descendants_of_children_before = 0;
  t8_linearidx_t num_descendants_of_child = 0;
  int childindex;
  /*If needed, can be replaced by binary search in lookuptable */
  for (childindex = 0; childindex < t8_element_num_children (elem); childindex++) {
    num_descendants_of_child = (1LL << (level_diff - 1) * 2);
    sum_descendants_of_children_before += num_descendants_of_child;
    if (sum_descendants_of_children_before > id) {
      sum_descendants_of_children_before -= num_descendants_of_child;
      break;
    }
  }
  t8_element_child (elem, childindex, elem);
  t8_hilbert_init_linear_id_recursive (elem, level_diff - 1, id - sum_descendants_of_children_before);
}

void
t8_default_scheme_quad_c::t8_element_set_linear_id (t8_element_t *elem, int level, t8_linearidx_t id) const
{
  t8_hilbert_root ((t8_hilbert_t *) elem);
  if (level == 0) {
    T8_ASSERT (id == 0);
    return;
  }
  t8_hilbert_init_linear_id_recursive (elem, level, id);
}

t8_linearidx_t
t8_default_scheme_quad_c::t8_element_get_linear_id (const t8_element_t *elem, int level) const
{
  t8_element_t *rec_start;
  t8_element_new (1, &rec_start);
  if (level > t8_element_level (elem)) {
    t8_element_first_descendant (elem, rec_start, level);
  }
  else {
    t8_element_copy (elem, rec_start);
    while (t8_element_level (rec_start) > level) {
      t8_element_parent (rec_start, rec_start);
    }
  }

  /* Maybe we can also input p into recursive function and calculate id directly for first desc */
  t8_linearidx_t id = t8_hilbert_linear_id_recursive (rec_start, 0, 0);
  T8_ASSERT (id >= 0);
  t8_element_destroy (1, &rec_start);
  return id;
}

void
t8_default_scheme_quad_c::t8_element_first_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  t8_element_copy (elem, desc);
  while (t8_element_level (desc) < level) {
    t8_element_child (desc, 0, desc);
  }
}

void
t8_default_scheme_quad_c::t8_element_last_descendant (const t8_element_t *elem, t8_element_t *desc, int level) const
{
  t8_element_copy (elem, desc);
  while (t8_element_level (desc) < level) {
    t8_element_child (desc, t8_element_num_children (desc) - 1, desc);
  }
}

void
t8_default_scheme_quad_c::t8_element_successor (const t8_element_t *elem1, t8_element_t *elem2, int level) const
{
  T8_ASSERT (level == t8_element_level (elem1));
  int child_id = t8_element_child_id (elem1);
  if (child_id + 1 == t8_element_num_siblings (elem1)) {
    t8_element_parent (elem1, elem2);
    t8_element_successor (elem2, elem2, level - 1);
    t8_element_child (elem2, 0, elem2);
  }
  else {
    t8_element_sibling (elem1, child_id + 1, elem2);
  }
}
/*****/
int
t8_default_scheme_quad_c::t8_element_ancestor_id (const t8_element_t *elem, int level) const
{
  SC_ABORT ("not implemented");
}

void
t8_default_scheme_quad_c::t8_element_nca (const t8_element_t *elem1, const t8_element_t *elem2, t8_element_t *nca) const
{
  SC_ABORT ("not implemented");
}
/*****/
int
t8_default_scheme_quad_c::t8_element_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_max_num_faces (const t8_element_t *elem) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_num_face_children (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_get_face_corner (const t8_element_t *element, int face, int corner) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_get_corner_face (const t8_element_t *element, int corner, int face) const
{
  SC_ABORT ("not implemented");
}

t8_element_shape_t
t8_default_scheme_quad_c::t8_element_face_shape (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

void
t8_default_scheme_quad_c::t8_element_children_at_face (const t8_element_t *elem, int face, t8_element_t *children[],
                                                       int num_children, int *child_indices) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_face_child_face (const t8_element_t *elem, int face, int face_child) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_face_parent_face (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

void
t8_default_scheme_quad_c::t8_element_transform_face (const t8_element_t *elem1, t8_element_t *elem2, int orientation,
                                                     int sign, int is_smaller_face) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_extrude_face (const t8_element_t *face, const t8_eclass_scheme_c *face_scheme,
                                                   t8_element_t *elem, int root_face) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_tree_face (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_default_scheme_quad_c::t8_element_first_descendant_face (const t8_element_t *elem, int face,
                                                            t8_element_t *first_desc, int level) const
{
  SC_ABORT ("not implemented");
}

/** Construct the last descendant of an element that touches a given face. */
void
t8_default_scheme_quad_c::t8_element_last_descendant_face (const t8_element_t *elem, int face, t8_element_t *last_desc,
                                                           int level) const
{
  SC_ABORT ("not implemented");
}

void
t8_default_scheme_quad_c::t8_element_boundary_face (const t8_element_t *elem, int face, t8_element_t *boundary,
                                                    const t8_eclass_scheme_c *boundary_scheme) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_is_root_boundary (const t8_element_t *elem, int face) const
{
  SC_ABORT ("not implemented");
}

int
t8_default_scheme_quad_c::t8_element_face_neighbor_inside (const t8_element_t *elem, t8_element_t *neigh, int face,
                                                           int *neigh_face) const
{
  SC_ABORT ("not implemented");
}

/***/
void
t8_default_scheme_quad_c::t8_element_anchor (const t8_element_t *elem, int coord[3]) const
{
  t8_hilbert_t *q;

  T8_ASSERT (t8_element_is_valid (elem));
  q = (t8_hilbert_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
}

void
t8_default_scheme_quad_c::t8_element_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                              double coords[]) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= vertex && vertex < t8_element_num_corners (elem));

  double ref_coords[2];
  ref_coords[0] = (vertex & 1 ? 1. : 0.);
  ref_coords[1] = (vertex & 2 ? 1. : 0.);
  t8_element_reference_coords (elem, ref_coords, 1, coords);
}

void
t8_default_scheme_quad_c::t8_element_reference_coords (const t8_element_t *el, const double *ref_coords,
                                                       const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (t8_element_is_valid (el));
  const t8_hilbert_t *elem = (const t8_hilbert_t *) el;
  for (size_t ipoint = 0; ipoint < num_coords; ipoint++) {
    out_coords[ipoint * 3 + 0]
      = (elem->x + ref_coords[ipoint * 3 + 0] * T8_HILBERT_LEN (elem->level)) / T8_HILBERT_ROOT_LEN;
    out_coords[ipoint * 3 + 1]
      = (elem->y + ref_coords[ipoint * 3 + 1] * T8_HILBERT_LEN (elem->level)) / T8_HILBERT_ROOT_LEN;
  }
}

/** memory */
void
t8_default_scheme_quad_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a hilbert_elem */
  t8_default_scheme_common_c::t8_element_new (length, elem);
  t8_element_init (length, elem[0], 0);
}

void
t8_default_scheme_quad_c::t8_element_init (int length, t8_element_t *elem, int new_called) const
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    t8_hilbert_t *current = (t8_hilbert_t *) elem;
    for (int ielem = 0; ielem < length; ielem++, current++) {
      t8_hilbert_root (current);
    }
  }
#endif
}

/** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
 * Returns false otherwise.
 */
int
t8_default_scheme_quad_c::t8_element_refines_irregular () const
{
  /*hilbert_elem refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG
int
t8_default_scheme_quad_c::t8_element_is_valid (const t8_element_t *elem) const
{
  const t8_hilbert_t *hilbert_elem = ((const t8_hilbert_t *) elem);
  if (hilbert_elem->level < 0)
    return 0;
  return 1;
}

void
t8_default_scheme_quad_c::t8_element_to_string (const t8_element_t *elem, char *debug_string,
                                                const int string_size) const
{
  //T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_hilbert_t *hilbert_elem = (t8_hilbert_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, type: %i, level: %i", hilbert_elem->x, hilbert_elem->y,
            hilbert_elem->type, hilbert_elem->level);
}
#endif

/* Constructor */
t8_default_scheme_quad_c::t8_default_scheme_quad_c (void)
{
  eclass = T8_ECLASS_QUAD;
  element_size = sizeof (t8_hilbert_t);
  ts_context = sc_mempool_new (element_size);
}

t8_default_scheme_quad_c::~t8_default_scheme_quad_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}
int
t8_default_scheme_quad_c::t8_element_pack (const t8_element_t *elements, int count, void *send_buffer, int buffer_size,
                                           int *position, sc_MPI_Comm comm) const
{
  int mpiret;
  p4est_quadrant_t *quads = (p4est_quadrant_t *) elements;
  for (int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(quads[ielem].x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads[ielem].y, 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&quads[ielem].level, 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
  return 0;
}

int
t8_default_scheme_quad_c::t8_element_pack_size (int count, sc_MPI_Comm comm, int *pack_size) const
{
  int singlesize = 0;
  int datasize = 0;
  int mpiret;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
  return 0;
}

int
t8_default_scheme_quad_c::t8_element_unpack (void *recvbuf, int buffer_size, int *position, t8_element_t *elements,
                                             int count, sc_MPI_Comm comm) const
{
  int mpiret;
  p4est_quadrant_t *quads = (p4est_quadrant_t *) elements;
  for (int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem].x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem].y), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(quads[ielem].level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
  return 0;
}

T8_EXTERN_C_END ();
