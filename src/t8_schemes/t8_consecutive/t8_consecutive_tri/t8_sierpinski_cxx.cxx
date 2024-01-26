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

#include <t8_schemes/t8_consecutive/t8_consecutive_common/t8_consecutive_common_cxx.hxx>
#include <t8_schemes/t8_consecutive/t8_consecutive_tri/t8_sierpinski_cxx.hxx>
#include <t8_schemes/t8_consecutive/t8_consecutive_tri/t8_sierpinski_connectivity.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

static int
t8_sierpinski_compute_cubeid (const t8_sierpinski_t *elem, int level)
{
  if (level == 0) {
    return 0;
  }
  int cube_id = 0;

  T8_ASSERT (0 <= elem->level && elem->level <= T8_SIERPINSKI_MAXLEVEL);
  const t8_sierpinski_coord_t h = T8_SIERPINSKI_LEN (level);

  cube_id |= ((elem->x & h) ? 1 : 0);
  cube_id |= ((elem->y & h) ? 2 : 0);
  return cube_id;
}
int8_t
t8_sierpinski_implicit_type_at_level (const t8_sierpinski_t *sierpinski, int level)
{
  if (level == 0)
    return T8_SIERPINSKI_ROOT_TYPE;
  T8_ASSERT (level > 0);
  const t8_sierpinski_coord_t length = T8_SIERPINSKI_LEN (level);
  return ((sierpinski->x & length) ^ (sierpinski->y & length)) ? 1 : 0;
}
int
t8_consecutive_scheme_tri_c::t8_element_maxlevel (void) const
{
  return T8_SIERPINSKI_MAXLEVEL;
}

int
t8_consecutive_scheme_tri_c::t8_element_level (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_sierpinski_t *) elem)->level;
}

void
t8_consecutive_scheme_tri_c::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
{
  const t8_sierpinski_t *source_elem = (const t8_sierpinski_t *) source;
  t8_sierpinski_t *dest_elem = (t8_sierpinski_t *) dest;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (source_elem == dest_elem) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *dest_elem = *source_elem;
}

int
t8_consecutive_scheme_tri_c::t8_element_equal (const t8_element_t *e1, const t8_element_t *e2) const
{
  const t8_sierpinski_t *elem1 = (const t8_sierpinski_t *) e1;
  const t8_sierpinski_t *elem2 = (const t8_sierpinski_t *) e2;
  return elem1->x == elem2->x && elem1->y == elem2->y && elem1->level == elem2->level && elem1->type == elem2->type;
}

void
t8_consecutive_scheme_tri_c::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_sierpinski_t *sierpinski_elem = (const t8_sierpinski_t *) elem;
  T8_ASSERT (sierpinski_elem->level > 0);

  t8_sierpinski_t *sierpinski_parent = (t8_sierpinski_t *) parent;

  const int8_t cube_id = t8_sierpinski_compute_cubeid (sierpinski_elem, sierpinski_elem->level);
  const int8_t parent_implicit_type
    = t8_sierpinski_implicit_type_at_level (sierpinski_elem, sierpinski_elem->level - 1);
  sierpinski_parent->type
    = t8_sierpinski_parentimplicittype_cubeid_type_to_parenttype[parent_implicit_type][cube_id][sierpinski_elem->type];

  /* remove least significant coordinate bit at level  */
  const t8_sierpinski_coord_t length = T8_SIERPINSKI_LEN (sierpinski_elem->level);
  sierpinski_parent->x = sierpinski_elem->x & ~length;
  sierpinski_parent->y = sierpinski_elem->y & ~length;

  sierpinski_parent->level = sierpinski_elem->level - 1;
  T8_ASSERT (sierpinski_parent->level >= 0);
}

t8_eclass_t
t8_consecutive_scheme_tri_c::t8_element_child_eclass (int childid) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_tri_c::t8_element_num_children (const t8_element_t *elem) const
{
  return T8_SIERPINSKI_CHILDREN;
}

void
t8_consecutive_scheme_tri_c::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
{
  const t8_sierpinski_t *sierpinski_elem = (const t8_sierpinski_t *) elem;
  t8_sierpinski_t *sierpinski_child = (t8_sierpinski_t *) child;

  T8_ASSERT (0 <= childid && childid < T8_SIERPINSKI_CHILDREN);
  T8_ASSERT (0 <= sierpinski_elem->level && sierpinski_elem->level <= T8_SIERPINSKI_MAXLEVEL);

  /* Compute the cube id and shift the coordinates accordingly */
  int8_t cube_id;
  int8_t parent_implicit_type = t8_sierpinski_implicit_type_at_level (sierpinski_elem, sierpinski_elem->level);
  cube_id = t8_sierpinski_parentimplicittype_parenttype_Iloc_to_childcubeid[parent_implicit_type][sierpinski_elem->type]
                                                                           [childid];
  sierpinski_child->type
    = t8_sierpinski_parentimplicittype_parenttype_Iloc_to_childtype[parent_implicit_type][sierpinski_elem->type]
                                                                   [childid];

  const t8_sierpinski_coord_t length = T8_SIERPINSKI_LEN (sierpinski_elem->level + 1);

  sierpinski_child->x = sierpinski_elem->x + ((cube_id & 1) ? length : 0);
  sierpinski_child->y = sierpinski_elem->y + ((cube_id & 2) ? length : 0);

  sierpinski_child->level = sierpinski_elem->level + 1;
}

int
t8_consecutive_scheme_tri_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_sierpinski_t *sierpinski = (const t8_sierpinski_t *) elem;
  T8_ASSERT (sierpinski->level >= 0);
  if (sierpinski->level == 0) {
    return -1;
  }
  const int8_t cube_id = t8_sierpinski_compute_cubeid (sierpinski, sierpinski->level);
  int8_t parent_implicit_type = t8_sierpinski_implicit_type_at_level (sierpinski, sierpinski->level - 1);
  int childid = t8_sierpinski_parentimplicittype_cubeid_type_to_Iloc[parent_implicit_type][cube_id][sierpinski->type];
  return childid;
}

/***/
void
t8_consecutive_scheme_tri_c::t8_element_anchor (const t8_element_t *elem, int coord[3]) const
{
  t8_sierpinski_t *q;

  T8_ASSERT (t8_element_is_valid (elem));
  q = (t8_sierpinski_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
}

void
t8_consecutive_scheme_tri_c::t8_element_vertex_reference_coords (const t8_element_t *elem, const int vertex,
                                                                 double coords[]) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= vertex && vertex < t8_element_num_corners (elem));

  double ref_coords[2];
  ref_coords[0] = (vertex >= 1) ? 1. : 0.;
  ref_coords[1] = (vertex >= 2) ? 1. : 0.;
  t8_element_reference_coords (elem, ref_coords, 1, coords);
}

void
t8_consecutive_scheme_tri_c::t8_element_reference_coords (const t8_element_t *el, const double *ref_coords,
                                                          const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (t8_element_is_valid (el));
  const t8_sierpinski_t *elem = (const t8_sierpinski_t *) el;
  double a00, a01, a10, a11;
  int implicit_type = t8_sierpinski_implicit_type_at_level (elem, elem->level);
  t8_sierpinski_coord_t x, y;
  x = elem->x;
  y = elem->y;
  if (elem->type) {
    x += T8_SIERPINSKI_LEN (elem->level);
    if (implicit_type) {
      /** 90° */
      a00 = 0;
      a01 = -1;
      a10 = 1;
      a11 = 0;
    }
    else {
      /** 180° */
      a00 = -1;
      a01 = 0;
      a10 = 0;
      a11 = -1;
      y += T8_SIERPINSKI_LEN (elem->level);
    }
  }
  else {
    if (implicit_type) {
      /** 270° */
      a00 = 0;
      a01 = 1;
      a10 = -1;
      a11 = 0;
      y += T8_SIERPINSKI_LEN (elem->level);
    }
    else {
      /** 0° */
      a00 = 1;
      a01 = 0;
      a10 = 0;
      a11 = 1;
    }
  }

  for (size_t ipoint = 0; ipoint < num_coords; ipoint++) {
    out_coords[ipoint * 3 + 0]
      = (x + (a00 * ref_coords[ipoint * 3 + 0] + a01 * ref_coords[ipoint * 3 + 1]) * T8_SIERPINSKI_LEN (elem->level))
        / T8_SIERPINSKI_ROOT_LEN;
    out_coords[ipoint * 3 + 1]
      = (y + (a10 * ref_coords[ipoint * 3 + 0] + a11 * ref_coords[ipoint * 3 + 1]) * T8_SIERPINSKI_LEN (elem->level))
        / T8_SIERPINSKI_ROOT_LEN;
  }
}

/** memory */
void
t8_consecutive_scheme_tri_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a sierpinski_elem */
  t8_consecutive_scheme_common_c::t8_element_new (length, elem);
  t8_element_init (length, elem[0], 0);
}

void
t8_consecutive_scheme_tri_c::t8_element_init (int length, t8_element_t *elem, int new_called) const
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    t8_sierpinski_t *current = (t8_sierpinski_t *) elem;
    for (int ielem = 0; ielem < length; ielem++, current++) {
      t8_element_root ((t8_element_t *) current);
    }
  }
#endif
}

/** Returns true, if there is one element in the tree, that does not refine into 2^dim children.
 * Returns false otherwise.
 */
int
t8_consecutive_scheme_tri_c::t8_element_refines_irregular () const
{
  /*sierpinski_elem refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG
int
t8_consecutive_scheme_tri_c::t8_element_is_valid (const t8_element_t *elem) const
{
  const t8_sierpinski_t *sierpinski_elem = ((const t8_sierpinski_t *) elem);
  if (sierpinski_elem->level < 0)
    return 0;
  return 1;
}

void
t8_consecutive_scheme_tri_c::t8_element_to_string (const t8_element_t *elem, char *debug_string,
                                                   const int string_size) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_sierpinski_t *sierpinski_elem = (t8_sierpinski_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, type: %i, level: %i", sierpinski_elem->x, sierpinski_elem->y,
            sierpinski_elem->type, sierpinski_elem->level);
}
#endif

/* Constructor */
t8_consecutive_scheme_tri_c::t8_consecutive_scheme_tri_c (void)
{
  eclass = T8_ECLASS_TRIANGLE;
  element_size = sizeof (t8_sierpinski_t);
  ts_context = sc_mempool_new (element_size);
}

t8_consecutive_scheme_tri_c::~t8_consecutive_scheme_tri_c ()
{
  /* This destructor is empty since the destructor of the
   * consecutive_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

int
t8_consecutive_scheme_tri_c::t8_element_pack (const t8_element_t *elements, int count, void *send_buffer,
                                              int buffer_size, int *position, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_sierpinski_t *sierpinski = (t8_sierpinski_t *) elements;
  for (int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Pack (&(sierpinski[ielem].x), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&(sierpinski[ielem].y), 1, sc_MPI_INT, send_buffer, buffer_size, position, comm);
    t8_debugf ("packed y %f\n", sierpinski[ielem].y);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Pack (&(sierpinski[ielem].type), 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&(sierpinski[ielem].level), 1, sc_MPI_INT8_T, send_buffer, buffer_size, position, comm);
    SC_CHECK_MPI (mpiret);
  }
  return 0;
}

int
t8_consecutive_scheme_tri_c::t8_element_pack_size (int count, sc_MPI_Comm comm, int *pack_size) const
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

  mpiret = sc_MPI_Pack_size (1, sc_MPI_INT8_T, comm, &datasize);
  SC_CHECK_MPI (mpiret);
  singlesize += datasize;

  *pack_size = count * singlesize;
  return 0;
}

int
t8_consecutive_scheme_tri_c::t8_element_unpack (void *recvbuf, int buffer_size, int *position, t8_element_t *elements,
                                                int count, sc_MPI_Comm comm) const
{
  int mpiret;
  t8_sierpinski_t *sierpinski = (t8_sierpinski_t *) elements;
  for (int ielem = 0; ielem < count; ielem++) {
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(sierpinski[ielem].x), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(sierpinski[ielem].y), 1, sc_MPI_INT, comm);
    t8_debugf ("unpacked y %f\n", sierpinski[ielem].y);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(sierpinski[ielem].type), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Unpack (recvbuf, buffer_size, position, &(sierpinski[ielem].level), 1, sc_MPI_INT8_T, comm);
    SC_CHECK_MPI (mpiret);
  }
  return 0;
}

T8_EXTERN_C_END ();
