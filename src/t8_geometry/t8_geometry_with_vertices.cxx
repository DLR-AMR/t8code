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

/* In this file we collect the implementations of the geometry_w_vertices
 * class.
 * */

#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_vec.h>

/* Load the coordinates of the newly active tree to the active_tree_vertices variable. */
void
t8_geometry_with_vertices::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  /* Set active id and eclass */
  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  active_tree = gtreeid;
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  if (0 <= ltreeid && ltreeid < num_local_trees) {
    active_tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  }
  else {
    active_tree_class = t8_cmesh_get_ghost_class (cmesh, ltreeid - num_local_trees);
  }
  /* Load this trees vertices. */
  active_tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  /* Check whether we support this class */
  T8_ASSERT (active_tree_class == T8_ECLASS_VERTEX || active_tree_class == T8_ECLASS_TRIANGLE
             || active_tree_class == T8_ECLASS_TET || active_tree_class == T8_ECLASS_QUAD
             || active_tree_class == T8_ECLASS_HEX || active_tree_class == T8_ECLASS_LINE
             || active_tree_class == T8_ECLASS_PRISM || active_tree_class == T8_ECLASS_PYRAMID);
}

bool
t8_geometry_with_vertices::t8_geom_tree_negative_volume () const
{
  if (t8_eclass_to_dimension[active_tree_class] <= 2) {
    /* Only three dimensional eclass do have a volume */
    return false;
  }
  T8_ASSERT (t8_eclass_to_dimension[active_tree_class]
             == 3);  // Should we include 4 dimensional classes, we need to catch that here and implement it.
  T8_ASSERT (active_tree_class == T8_ECLASS_TET || active_tree_class == T8_ECLASS_HEX
             || active_tree_class == T8_ECLASS_PRISM || active_tree_class == T8_ECLASS_PYRAMID);

  T8_ASSERT (t8_eclass_num_vertices[active_tree_class] >= 4);
  /*
   *      6 ______  7  For Hexes and pyramids, if the vertex 4 is below the 0-1-2-3 plane,
   *       /|     /     the volume is negative. This is the case if and only if
   *    4 /_____5/|     the scalar product of v_4 with the cross product of v_1 and v_2 is
   *      | | _ |_|     smaller 0:
   *      | 2   | / 3   < v_4, v_1 x v_2 > < 0
   *      |/____|/
   *     0      1
   *
   *
   *    For tets/prisms, if the vertex 3 is below/above the 0-1-2 plane, the volume
   *    is negative. This is the case if and only if
   *    the scalar product of v_3 with the cross product of v_1 and v_2 is
   *    greater 0:
   *
   *    < v_3, v_1 x v_2 > > 0
   *
   */

  /* build the vectors v_i as vertices_i - vertices_0 */
  double v_1[3], v_2[3], v_j[3], cross[3], sc_prod;
  int i, j;
  if (active_tree_class == T8_ECLASS_TET || active_tree_class == T8_ECLASS_PRISM) {
    /* In the tet/prism case, the third vector is v_3 */
    j = 3;
  }
  else {
    /* For pyramids and Hexes, the third vector is v_4 */
    j = 4;
  }
  for (i = 0; i < 3; i++) {
    v_1[i] = active_tree_vertices[3 + i] - active_tree_vertices[i];
    v_2[i] = active_tree_vertices[6 + i] - active_tree_vertices[i];
    v_j[i] = active_tree_vertices[3 * j + i] - active_tree_vertices[i];
  }
  /* compute cross = v_1 x v_2 */
  t8_vec_cross (v_1, v_2, cross);
  /* Compute sc_prod = <v_j, cross> */
  sc_prod = t8_vec_dot (v_j, cross);

  for (i = 0; i < 3; ++i) {
    t8_debugf ("Vertex %i: %f %f %f\n", i, active_tree_vertices[3 * i], active_tree_vertices[3 * i + 1],
               active_tree_vertices[3 * i + 2]);
  }

  T8_ASSERT (sc_prod != 0);
  return active_tree_class == T8_ECLASS_TET ? sc_prod > 0 : sc_prod < 0;
}
