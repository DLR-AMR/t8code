/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_extruded.cxx
 * Implementation of the extruded scheme interface.
 */

#include <t8_schemes/t8_scheme_builder.hxx>
#include <t8_schemes/t8_extruded/t8_extruded.hxx>
#include <t8_schemes/t8_extruded/t8_extruded_scheme.hxx>
#include <t8_eclass/t8_eclass.h>

const t8_scheme *
t8_scheme_new_extruded ()
{
  t8_scheme_builder builder;

  builder.add_eclass_scheme<t8_default_scheme_vertex> ();
  builder.add_eclass_scheme<t8_default_scheme_line> ();
  builder.add_eclass_scheme<t8_default_scheme_quad> ();
  builder.add_eclass_scheme<t8_default_scheme_tri> ();
  builder.add_eclass_scheme<t8_extruded_scheme_hex> ();
  builder.add_eclass_scheme<t8_default_scheme_tet> ();
  builder.add_eclass_scheme<t8_default_scheme_prism> ();
  builder.add_eclass_scheme<t8_default_scheme_pyramid> ();

  return builder.build_scheme ();
}

int
t8_eclass_scheme_is_extruded (const t8_scheme *scheme, const t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_HEX:
    return scheme->check_eclass_scheme_type<t8_extruded_scheme_hex> (T8_ECLASS_HEX);
  default:
    return 0;
  }
}

/** Check whether a single face connection of a hex tree is compatible with the extruded hex scheme.
 * \param [in] face         A face of the hex tree.
 * \param [in] neigh_class  The eclass of the neighbor tree.
 * \param [in] neigh_face   The face of the neighbor tree.
 * \param [in] orientation  The orientation of the face connection.
 * \return                  True if the connection is compatible.
 */
// TODO extruded: static for local helper?
static bool
t8_extruded_hex_face_connection_is_compatible (const int face, const t8_eclass_t neigh_class, const int neigh_face,
                                               const int orientation)
{
  constexpr int num_lateral_faces = t8_extruded_scheme_hex::num_lateral_faces;
  const bool face_is_lateral = face < num_lateral_faces;
  if (neigh_class != T8_ECLASS_HEX) {
    /* Only the bottom and top face may be connected to a quad face of a tree of a different class. */
    return !face_is_lateral && t8_eclass_face_types[neigh_class][neigh_face] == T8_ECLASS_QUAD;
  }
  const bool neigh_face_is_lateral = neigh_face < num_lateral_faces;
  if (face_is_lateral != neigh_face_is_lateral) {
    /* The z-axis of one tree would be in-plane of the other tree. */
    return false;
  }
  if (!face_is_lateral) {
    /* Bottom/top faces are genuine quads, any orientation is fine. */
    return true;
  }
  /* The lateral face coordinates are (in-plane, z). The z-axes are parallel if and only if the face map does not swap
   * the two face coordinates. Depending on whether the faces have the same topological orientation (sign), this is
   * the case for the following orientations. */
  // TODO extruded: correct?
  const int sign
    = t8_eclass_face_orientation[T8_ECLASS_HEX][face] == t8_eclass_face_orientation[T8_ECLASS_HEX][neigh_face];
  if (sign) {
    return orientation == 1 || orientation == 2;
  }
  return orientation == 0 || orientation == 3;
}

bool
t8_cmesh_is_extrusion_compatible (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t ltree = 0; ltree < num_local_trees; ++ltree) {
    // TODO extruded: no hex, no check?
    if (t8_cmesh_get_tree_class (cmesh, ltree) != T8_ECLASS_HEX) {
      continue;
    }
    for (int face = 0; face < t8_eclass_num_faces[T8_ECLASS_HEX]; ++face) {
      int neigh_face;
      int orientation;
      const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, ltree, face, &neigh_face, &orientation);
      if (neigh_tree < 0) {
        /* No neighbor across this face. */
        continue;
      }
      const t8_eclass_t neigh_class = t8_cmesh_get_tree_face_neighbor_eclass (cmesh, ltree, face);
      if (!t8_extruded_hex_face_connection_is_compatible (face, neigh_class, neigh_face, orientation)) {
        t8_debugf ("Face %i of local tree %i is not compatible with the extruded hex scheme.\n", face, ltree);
        return false;
      }
    }
  }
  return true;
}
