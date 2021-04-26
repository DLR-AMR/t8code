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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_bspline.hxx>

t8_geometry_bspline::t8_geometry_bspline (int dim, const char *name_in,
                                            t8_geom_load_tree_data_fn
                                            load_tree_data_in,
                                            const Handle_Geom_BSplineSurface &bspline_in)
: bspline(bspline_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;

  load_tree_data = load_tree_data_in;
}

/**
 * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_bspline::t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const
{
  gp_Pnt pnt;
  T8_ASSERT (!bspline.IsNull());
  bspline->D0(ref_coords[0], ref_coords[1], pnt);

  out_coords[0] = pnt.X();
  out_coords[2] = pnt.Z();
  out_coords[1] = pnt.Y();
}

void
t8_geometry_bspline::t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double
                                                *ref_coords,
                                                double *jacobian_out) const
{
  
}

void
t8_geometry_bspline::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid)
{
  if (load_tree_data != NULL) {
    /* Load tree data if a loading function was provided. */
    load_tree_data (cmesh, gtreeid, &tree_data);
  }
  else {
    /* Otherwise it is NULL. */
    tree_data = NULL;
  }
}