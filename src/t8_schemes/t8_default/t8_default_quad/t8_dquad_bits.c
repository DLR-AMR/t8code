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
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.h>
#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>

static void
t8_element_copy_surround (const p4est_quadrant_t *q, p4est_quadrant_t *r)
{
  T8_QUAD_SET_TDIM (r, T8_QUAD_GET_TDIM (q));
  if (T8_QUAD_GET_TDIM (q) == 3) {
    T8_QUAD_SET_TNORMAL (r, T8_QUAD_GET_TNORMAL (q));
    T8_QUAD_SET_TCOORD (r, T8_QUAD_GET_TCOORD (q));
  }
}

void
t8_dquad_compute_reference_coords (const p4est_quadrant_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords)
{
  const p4est_qcoord_t h = P4EST_QUADRANT_LEN (elem->level);

  for (size_t icoord = 0; icoord < num_coords; ++icoord) {
    const size_t offset_2d = icoord * 2;
    const size_t offset_3d = icoord * 3;
    out_coords[offset_2d + 0] = elem->x + ref_coords[offset_3d + 0] * h;
    out_coords[offset_2d + 1] = elem->y + ref_coords[offset_3d + 1] * h;

    out_coords[offset_2d + 0] /= (double) P4EST_ROOT_LEN;
    out_coords[offset_2d + 1] /= (double) P4EST_ROOT_LEN;
  }
}

void
t8_dquad_child (const p4est_quadrant_t *elem, int childid, p4est_quadrant_t *child)
{
  const p4est_qcoord_t shift = P4EST_QUADRANT_LEN (elem->level + 1);

  T8_ASSERT (p4est_quadrant_is_extended (elem));
  T8_ASSERT (elem->level < P4EST_QMAXLEVEL);
  T8_ASSERT (childid >= 0 && childid < P4EST_CHILDREN);

  child->x = childid & 0x01 ? (elem->x | shift) : elem->x;
  child->y = childid & 0x02 ? (elem->y | shift) : elem->y;
  child->level = elem->level + 1;

  if (elem != child) {
    T8_ASSERT (p4est_quadrant_is_parent (elem, child));
  }
  t8_element_copy_surround (elem, child);
}

void
t8_dquad_parent (const p4est_quadrant_t *elem, p4est_quadrant_t *parent)
{
  T8_ASSERT (p4est_quadrant_is_extended (elem));
  p4est_quadrant_parent (elem, parent);
  t8_element_copy_surround (elem, parent);
}

void
t8_dquad_sibling (const p4est_quadrant_t *elem, int sibid, p4est_quadrant_t *sibling)
{
  T8_ASSERT (p4est_quadrant_is_extended (elem));
  p4est_quadrant_sibling (elem, sibling, sibid);
  t8_element_copy_surround (elem, sibling);
}

void
t8_dquad_copy (const p4est_quadrant_t *source, p4est_quadrant_t *dest)
{
  T8_ASSERT (p4est_quadrant_is_extended (source));
  if (source == dest) {
    /* Do nothing if they are already the same quadrant. */
    return;
  }
  *dest = *source;
  t8_element_copy_surround (source, dest);
}

void
t8_dquad_successor (const p4est_quadrant_t *elem, p4est_quadrant_t *succ, const int level, const int multilevel)
{
  T8_ASSERT (p4est_quadrant_is_extended (elem));
  T8_ASSERT (0 <= elem->level && succ->level <= P4EST_QMAXLEVEL);
  T8_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);
  T8_ASSERT (elem->level <= level);
  /* If level not reached, construct child */
  if (elem->level < level) {
    t8_dquad_child (elem, 0, succ);
  }
  else {
    /* If level reached, construct sibling */
    const int num_siblings = 4;
    int child_id = p4est_quadrant_child_id (elem);
    T8_ASSERT (0 <= child_id && child_id < num_siblings);
    if (child_id < num_siblings - 1) {
      t8_dquad_sibling (elem, child_id + 1, succ);
    }
    else {
      /* elem is last sibling, go up until elem is not last sibling */
      t8_dquad_parent (elem, succ);
      child_id = p4est_quadrant_child_id (succ);
      T8_ASSERT (0 <= child_id && child_id < num_siblings);
      while (child_id == num_siblings - 1) {
        t8_dquad_parent (succ, succ);
        child_id = p4est_quadrant_child_id (succ);
        T8_ASSERT (0 <= child_id && child_id < num_siblings);
      }
      /* Construct next sibling */
      t8_dquad_sibling (succ, child_id + 1, succ);
    }
  }
  if (!multilevel) {
    /* If not multilevel, traverse down the tree */
    while (succ->level < level) {
      t8_dquad_child (succ, 0, succ);
    }
  }
}
