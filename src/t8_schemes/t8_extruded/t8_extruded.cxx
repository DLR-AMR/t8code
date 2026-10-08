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
  builder.add_eclass_scheme<t8_extruded_scheme_prism> ();
  builder.add_eclass_scheme<t8_default_scheme_pyramid> ();

  return builder.build_scheme ();
}

int
t8_eclass_scheme_is_extruded (const t8_scheme *scheme, const t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_HEX:
    return scheme->check_eclass_scheme_type<t8_extruded_scheme_hex> (T8_ECLASS_HEX);
  case T8_ECLASS_PRISM:
    return scheme->check_eclass_scheme_type<t8_extruded_scheme_prism> (T8_ECLASS_PRISM);
  default:
    return 0;
  }
}

/** Return the number of lateral faces of an extruded eclass.
 * \param [in] eclass   An eclass.
 * \return              The number of lateral faces if \a eclass is extruded (hex or prism), -1 otherwise.
 */
static int
t8_extruded_num_lateral_faces (const t8_eclass_t eclass)
{
  switch (eclass) {
  case T8_ECLASS_HEX:
    return t8_extruded_scheme_hex::num_lateral_faces;
  case T8_ECLASS_PRISM:
    return t8_extruded_scheme_prism::num_lateral_faces;
  default:
    return -1;
  }
}

/** Check whether a single face connection between two extruded trees is compatible with the extruded schemes.
 * \param [in] eclass       The eclass of the tree, hex or prism.
 * \param [in] face         A face of the tree.
 * \param [in] neigh_class  The eclass of the neighbor tree, hex or prism.
 * \param [in] neigh_face   The face of the neighbor tree.
 * \param [in] orientation  The orientation of the face connection.
 * \return                  True if the connection is compatible.
 */
// TODO extruded: static for local helper?
static bool
t8_extruded_face_connection_is_compatible (const t8_eclass_t eclass, const int face, const t8_eclass_t neigh_class,
                                           const int neigh_face, const int orientation)
{
  const int neigh_num_lateral_faces = t8_extruded_num_lateral_faces (neigh_class);
  T8_ASSERT (neigh_num_lateral_faces > 0);
  const bool face_is_lateral = face < t8_extruded_num_lateral_faces (eclass);
  const bool neigh_face_is_lateral = neigh_face < neigh_num_lateral_faces;
  if (face_is_lateral != neigh_face_is_lateral) {
    /* The z-axis of one tree would be in-plane of the other tree. */
    return false;
  }
  if (!face_is_lateral) {
    /* Bottom/top faces are genuine base elements, any orientation is fine. */
    return true;
  }
  /* The lateral faces are quads with coordinates (in-plane, z) for both hexes and prisms. The z-axes are parallel if
   * and only if the face map does not swap the two face coordinates. Depending on whether the faces have the same
   * topological orientation (sign), this is the case for the following orientations. */
  // TODO extruded: correct?
  const int sign = t8_eclass_face_orientation[eclass][face] == t8_eclass_face_orientation[neigh_class][neigh_face];
  if (sign) {
    return orientation == 1 || orientation == 2;
  }
  return orientation == 0 || orientation == 3;
}

bool
t8_extruded_cmesh_tree_is_compatible (const t8_scheme *scheme, const t8_cmesh_t cmesh, const t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltreeid);
  T8_ASSERT (t8_eclass_scheme_is_extruded (scheme, eclass));
  for (int face = 0; face < t8_eclass_num_faces[eclass]; ++face) {
    int neigh_face;
    int orientation;
    const t8_locidx_t neigh_tree = t8_cmesh_get_face_neighbor (cmesh, ltreeid, face, &neigh_face, &orientation);
    if (neigh_tree < 0) {
      /* No neighbor across this face. */
      continue;
    }
    const t8_eclass_t neigh_class = t8_cmesh_get_tree_face_neighbor_eclass (cmesh, ltreeid, face);
    if (!t8_eclass_scheme_is_extruded (scheme, neigh_class)
        || !t8_extruded_face_connection_is_compatible (eclass, face, neigh_class, neigh_face, orientation)) {
      t8_debugf ("Face %i of local tree %i is not compatible with the extruded schemes.\n", face, ltreeid);
      return false;
    }
  }
  return true;
}
