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
#include <t8_schemes/t8_consecutive/t8_consecutive_quad/t8_hilbert_cxx.hxx>
#include <t8_schemes/t8_consecutive/t8_consecutive_quad/t8_hilbert_connectivity.h>

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

int
t8_consecutive_scheme_quad_c::t8_element_maxlevel (void) const
{
  return T8_HILBERT_MAXLEVEL;
}

int
t8_consecutive_scheme_quad_c::t8_element_level (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_hilbert_t *) elem)->level;
}

void
t8_consecutive_scheme_quad_c::t8_element_copy (const t8_element_t *source, t8_element_t *dest) const
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
t8_consecutive_scheme_quad_c::t8_element_equal (const t8_element_t *e1, const t8_element_t *e2) const
{
  const t8_hilbert_t *elem1 = (const t8_hilbert_t *) e1;
  const t8_hilbert_t *elem2 = (const t8_hilbert_t *) e2;
  return elem1->x == elem2->x && elem1->y == elem2->y && elem1->level == elem2->level && elem1->type == elem2->type;
}

void
t8_consecutive_scheme_quad_c::t8_element_parent (const t8_element_t *elem, t8_element_t *parent) const
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

t8_eclass_t
t8_consecutive_scheme_quad_c::t8_element_child_eclass (int childid) const
{
  SC_ABORT ("not implemented");
}

int
t8_consecutive_scheme_quad_c::t8_element_num_children (const t8_element_t *elem) const
{
  return T8_HILBERT_CHILDREN;
}

void
t8_consecutive_scheme_quad_c::t8_element_child (const t8_element_t *elem, int childid, t8_element_t *child) const
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

int
t8_consecutive_scheme_quad_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_hilbert_t *hilbert = (const t8_hilbert_t *) elem;
  T8_ASSERT (hilbert->level >= 0);
  if (hilbert->level == 0) {
    return -1;
  }
  const int8_t cube_id = t8_hilbert_compute_cubeid (hilbert, hilbert->level);
  return t8_hilbert_type_cubeid_to_Iloc[hilbert->type][cube_id];
}

/***/
void
t8_consecutive_scheme_quad_c::t8_element_anchor (const t8_element_t *elem, int coord[3]) const
{
  t8_hilbert_t *q;

  T8_ASSERT (t8_element_is_valid (elem));
  q = (t8_hilbert_t *) elem;
  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = 0;
}

void
t8_consecutive_scheme_quad_c::t8_element_vertex_reference_coords (const t8_element_t *elem, const int vertex,
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
t8_consecutive_scheme_quad_c::t8_element_reference_coords (const t8_element_t *el, const double *ref_coords,
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
t8_consecutive_scheme_quad_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory for a hilbert_elem */
  t8_consecutive_scheme_common_c::t8_element_new (length, elem);
  t8_element_init (length, elem[0], 0);
}

void
t8_consecutive_scheme_quad_c::t8_element_init (int length, t8_element_t *elem, int new_called) const
{
#ifdef T8_ENABLE_DEBUG
  if (!new_called) {
    t8_hilbert_t *current = (t8_hilbert_t *) elem;
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
t8_consecutive_scheme_quad_c::t8_element_refines_irregular () const
{
  /*hilbert_elem refine regularly */
  return 0;
}

#ifdef T8_ENABLE_DEBUG
int
t8_consecutive_scheme_quad_c::t8_element_is_valid (const t8_element_t *elem) const
{
  const t8_hilbert_t *hilbert_elem = ((const t8_hilbert_t *) elem);
  if (hilbert_elem->level < 0)
    return 0;
  return 1;
}

void
t8_consecutive_scheme_quad_c::t8_element_to_string (const t8_element_t *elem, char *debug_string,
                                                    const int string_size) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (debug_string != NULL);
  t8_hilbert_t *hilbert_elem = (t8_hilbert_t *) elem;
  snprintf (debug_string, string_size, "x: %i, y: %i, type: %i, level: %i", hilbert_elem->x, hilbert_elem->y,
            hilbert_elem->type, hilbert_elem->level);
}
#endif

/* Constructor */
t8_consecutive_scheme_quad_c::t8_consecutive_scheme_quad_c (void)
{
  eclass = T8_ECLASS_QUAD;
  element_size = sizeof (t8_hilbert_t);
  ts_context = sc_mempool_new (element_size);
}

t8_consecutive_scheme_quad_c::~t8_consecutive_scheme_quad_c ()
{
  /* This destructor is empty since the destructor of the
   * consecutive_common scheme is called automatically and it
   * suffices to destroy the quad_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
