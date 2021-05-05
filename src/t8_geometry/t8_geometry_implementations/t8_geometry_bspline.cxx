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

Handle_Geom_BSplineSurface t8_global_bspline;

t8_geometry_bspline::t8_geometry_bspline (int dim, const char *name_in,
                                          const void *user_data_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;
  user_data = user_data_in;
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
  t8_cmesh_attribute_bspline *attribute = (t8_cmesh_attribute_bspline *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                                                                 T8_CMESH_BSPLINE_ATTRIBUTE_KEY,
                                                                                                 gtreeid);
  int face1 = attribute->t8_cmesh_attribute_bspline_get_face();
  const double *parameters = (const double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                                      T8_CMESH_BSPLINE_PARAMETERS_ATTRIBUTE_KEY + face1,
                                                                      gtreeid);
  Handle_Geom_BSplineSurface bspline = t8_global_bspline;
  
  /* Compute opposite face */
  int face2;
  face1 % 2 == 0 ? face2 = face1 + 1 : face2 = face1 - 1;
  
  /* Sequence of normal vectors to do linear interpolations */
  int normal_sequence[3][3] = 
  {
    {1, 2, 0},
    {0, 2, 1},
    {0, 1, 2}
  };
  
  /* Bilinear interpolation for parameters in bspline face1*/  
  double param[2];

  for (int i = 0; i < 2; ++i)
  {
    param[i] = (parameters[0 + i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][0]])
              + parameters[2 + i] * (ref_coords[normal_sequence[(int)(face1 / 2)][0]]))
              * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][1]]);
    param[i] += (parameters[4 + i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][0]])
                +parameters[6 + i] * (ref_coords[normal_sequence[(int)(face1 / 2)][0]]))
                *(ref_coords[normal_sequence[(int)(face1 / 2)][1]]);
  }
  
  /* Check if parameters are valid */
  T8_ASSERT(!bspline.IsNull());
  double u1, u2, v1, v2;
  bspline->Bounds(u1, u2, v1, v2);
  T8_ASSERT(((u1 <= param[0] && param[0] <= u2) || (u2 <= param[0] && param[0] <= u1))
            && ((v1 <= param[1] && param[1] <= v2) || (v2 <= param[1] && param[1] <= v1)));

  /* Compute p1 on face1 with parameters */
  gp_Pnt pnt;
  bspline->D0(param[0], param[1], pnt);
  double p1[] = {pnt.X(), pnt.Y(), pnt.Z()};  
  
  /* Bilinear interpolation for p2 on face2 */
  double p2[3];
  for (int i = 0; i < 3; ++i)
  {
    p2[i] = (active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][face2][0] * 3 + i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][0]])
            +active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][face2][1] * 3 + i] * (ref_coords[normal_sequence[(int)(face1 / 2)][0]]))
            *(1 - ref_coords[normal_sequence[(int)(face1 / 2)][1]]);
    p2[i] += (active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][face2][2] * 3 + i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][0]])
            + active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][face2][3] * 3 + i] * (ref_coords[normal_sequence[(int)(face1 / 2)][0]]))
            * (ref_coords[normal_sequence[(int)(face1 / 2)][1]]);
  }

  /* Linear interpolation between p1 and p2 to get out_coords */
  if (face1 % 2 == 0)
  {
    for (int i = 0; i < 3; ++i)
    {
      out_coords[i] = p1[i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][2]])
                    + p2[i] * (ref_coords[normal_sequence[(int)(face1 / 2)][2]]);
    }
  }
  else
  {
    for (int i = 0; i < 3; ++i)
    {
      out_coords[i] = p2[i] * (1 - ref_coords[normal_sequence[(int)(face1 / 2)][2]])
                    + p1[i] * (ref_coords[normal_sequence[(int)(face1 / 2)][2]]);
    }
  }
}

void
t8_geometry_bspline::t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double
                                                *ref_coords,
                                                double *jacobian_out) const
{
  
}

inline void
t8_geometry_bspline::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                            t8_gloidx_t gtreeid)
{
  /* Set active id and eclass */
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  active_tree = gtreeid;
  active_tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Load this trees vertices. */
  active_tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  T8_ASSERT (t8_eclass_to_dimension[active_tree_class] == dimension);

  /* Check whether we support this class */
  T8_ASSERT (active_tree_class == T8_ECLASS_VERTEX
             || active_tree_class == T8_ECLASS_TRIANGLE
             || active_tree_class == T8_ECLASS_TET
             || active_tree_class == T8_ECLASS_QUAD
             || active_tree_class == T8_ECLASS_HEX
             || active_tree_class == T8_ECLASS_LINE
             || active_tree_class == T8_ECLASS_PRISM);
}

t8_cmesh_attribute_bspline::t8_cmesh_attribute_bspline (const int bspline_index_in,
                                                        const int face_in)
{
  bspline_index = bspline_index_in;
  face =  face_in;
}