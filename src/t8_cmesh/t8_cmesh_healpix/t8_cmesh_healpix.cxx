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

t8_cmesh_t
t8_cmesh_new_healpix (sc_MPI_Comm comm)
{
  /* Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  const int ntrees = 12;
  const int nverts = 4; /* Number of vertices per cmesh element. */
  t8_eclass_t all_eclasses[ntrees];
  std::vector<double> all_verts;
  all_verts.reserve(ntrees * nverts * 3);

  /* Register geometry and retain the pointer */
  t8_geometry_c *geom = t8_cmesh_register_geometry<t8_geometry_healpix> (cmesh);

  /* Reference coordinates for the 4 corners of a quad element */
  double ref_corners[4][2] = {
    {0.0, 0.0},
    {1.0, 0.0},
    {1.0, 1.0},
    {0.0, 1.0}
  };

  /* Build trees for all 3 layers (upper, middle, lower) */
  for (int itree = 0; itree < ntrees; itree++) {
    const t8_gloidx_t layer = itree / 4;
    const t8_gloidx_t face = itree % 4;

    t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
    all_eclasses[itree] = T8_ECLASS_QUAD;

    /* Associate the custom HEALPix geometry with this tree */
    t8_cmesh_set_tree_geometry (cmesh, itree, geom);

    std::vector<double> verts;
    verts.reserve (nverts * 3);

    for (int i = 0; i < nverts; i++) {
      double coord[3];
      const double xi = std::clamp (ref_corners[i][0], 1e-10, 1.0 - 1e-10);
      const double eta = std::clamp (ref_corners[i][1], 1e-10, 1.0 - 1e-10);

      t8_eval_geom_point (layer, face, xi, eta, coord);
      verts.push_back (coord[0]);
      verts.push_back (coord[1]);
      verts.push_back (coord[2]);
    }

    all_verts.insert (all_verts.end(), verts.begin(), verts.end());
    t8_cmesh_set_tree_vertices (cmesh, itree, verts.data(), nverts);
  }

  /* Compute face connectivity using topological vertices */
  t8_cmesh_set_join_by_vertices (cmesh, 12, all_eclasses, all_verts.data(), nullptr, 0);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}