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

#include <cmath>
#include <vector>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>

struct t8_healpix_join
{
  t8_gloidx_t tree1;
  int face1;
  t8_gloidx_t tree2;
  int face2;
  int orient;
};

static const struct t8_healpix_join healpix_joins[24] = {
  /* North Ring Connections */
  { 0, 1, 1, 3, 0 },
  { 0, 3, 3, 1, 0 },
  { 1, 1, 2, 3, 0 },
  { 2, 1, 3, 3, 0 },

  /* North-to-Equator Connections */
  { 0, 0, 4, 2, 0 },
  { 0, 2, 5, 3, 1 },
  { 1, 0, 5, 2, 0 },
  { 1, 2, 6, 3, 1 },
  { 2, 0, 6, 2, 0 },
  { 2, 2, 7, 3, 1 },
  { 3, 0, 7, 2, 0 },
  { 3, 2, 4, 3, 1 },

  /* Equatorial Connections */
  { 4, 1, 5, 0, 1 },
  { 5, 1, 6, 0, 1 },
  { 6, 1, 7, 0, 1 },
  { 7, 1, 4, 0, 1 },

  /* Equator-to-South Connections */
  { 4, 0, 8, 3, 1 },
  { 5, 0, 9, 3, 1 },
  { 6, 0, 10, 3, 1 },
  { 7, 0, 11, 3, 1 },

  /* South Ring Connections */
  { 8, 1, 9, 3, 0 },
  { 8, 2, 11, 0, 0 },
  { 9, 1, 10, 3, 0 },
  { 10, 1, 11, 3, 0 }
};
t8_cmesh_t
t8_cmesh_new_healpix (sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  const int ntrees = 12;

  /* Register geometry and retain the pointer */
  t8_cmesh_register_geometry<t8_geometry_healpix> (cmesh);

  for (int itree = 0; itree < ntrees; itree++) {
    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
  }
  /* Explicit O(1) manual face joining */
  for (size_t i = 0; i < 24; ++i) {
    t8_cmesh_set_join (cmesh, healpix_joins[i].tree1, healpix_joins[i].tree2, healpix_joins[i].face1,
                       healpix_joins[i].face2, healpix_joins[i].orient);
  }

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}
