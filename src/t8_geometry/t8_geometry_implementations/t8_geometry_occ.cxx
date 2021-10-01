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

#include <t8.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_eclass.h>
#include <t8_geometry/t8_geometry_helpers.h>

#include <gp_Pnt.hxx>

#if T8_WITH_OCC

Handle_Geom_Curve t8_global_occ_curve[10];
Handle_Geom_Surface t8_global_occ_surface[10];

t8_geometry_occ::t8_geometry_occ (int dim, const char *name_in,
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
t8_geometry_occ::t8_geom_evaluate (t8_cmesh_t cmesh,
                                  t8_gloidx_t gtreeid,
                                  const double *ref_coords,
                                  double out_coords[3]) const
{
  const int *edges = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                          T8_CMESH_OCC_CURVE_ATTRIBUTE_KEY,
                                                          gtreeid);
  const int *faces = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                          T8_CMESH_OCC_SURFACE_ATTRIBUTE_KEY,
                                                          gtreeid);
  
  /* Compute coordinates via trilinear interpolation */
  t8_geom_compute_linear_geometry (active_tree_class,
                                   active_tree_vertices, ref_coords,
                                   out_coords);
  
  const int
  t8_edge_vertex_to_tree_vertex[T8_ECLASS_MAX_EDGES][2] =
  {
    {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {4, 6}, {1, 3}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}    /* hex */
  };

  double interpolated_coords[3], param[2], cur_delta[3];
  gp_Pnt pnt;
  
  /* Check each edge for geometry. Currently, only hexes with 12 edges are supported. */
  for (int i_edges = 0; i_edges < 12; ++i_edges)
  {
    if (edges[i_edges] >= 0)
    {
      /* Interpolate coordinates between edge vertices. Due to the indices i_edges of the edges, the edges point in
      * direction of ref_coord i_edges / 4. Therefore, we can use ref_coords[i_edges / 4] for the interpolation.              
      *          6 -------E3------- 7
      *         /|                 /|
      *       E5 |               E7 |
      *       / E10              / E11
      *      /   |              /   |          z y
      *     4 -------E2------- 5    |          |/
      *     |    |             |    |          x-- x
      *     |    2 -------E1---|--- 3
      *    E8   /             E9   /
      *     |  E4              |  E6
      *     | /                | /
      *     |/                 |/
      *     0 -------E0------- 1
      *        
      */
      for (int i_coord_dim = 0; i_coord_dim < 3; ++i_coord_dim)
      {
        interpolated_coords[i_coord_dim] = (1 - ref_coords[i_edges / 4]) * active_tree_vertices[t8_edge_vertex_to_tree_vertex[i_edges][0] * 3 + i_coord_dim]
                                + ref_coords[i_edges / 4] * active_tree_vertices[t8_edge_vertex_to_tree_vertex[i_edges][1] * 3 + i_coord_dim];
      }

      /* Interpolate parameters between edge vertices. Same procedure as above. */
      const double *parameters = (double *) t8_cmesh_get_attribute(cmesh, t8_get_package_id (),
                                                                  T8_CMESH_OCC_CURVE_PARAMETERS_ATTRIBUTE_KEY + i_edges,
                                                                  gtreeid);
      param[0] = (1 - ref_coords[i_edges / 4]) * parameters[0]
                  + ref_coords[i_edges / 4] * parameters[1];
      
      /* Check if calculated parameters are valid */
      T8_ASSERT(!t8_global_occ_curve[edges[i_edges]].IsNull());
      double u1, u2;
      u1 = t8_global_occ_curve[edges[i_edges]]->FirstParameter();
      u2 = t8_global_occ_curve[edges[i_edges]]->LastParameter();
      T8_ASSERT((u1 <= param[0] && param[0] <= u2) || (u2 <= param[0] && param[0] <= u1));

      /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
      t8_global_occ_curve[edges[i_edges]]->D0(param[0], pnt);
      cur_delta[0] = pnt.X() - interpolated_coords[0];
      cur_delta[1] = pnt.Y() - interpolated_coords[1];
      cur_delta[2] = pnt.Z() - interpolated_coords[2];
      
      /* Multiply curve displacement with corresponding ref coords.
      *  The edges are indexed so that all edges which satisfy i_edges % 4 == 0 have to multiplied with the inversed (1 - ref_coord) 
      *  coordinate. All edges which satisfy i_edges % 4 == 1 have to multiplied with one inversed ref_coord and so forth...
      */
      switch (i_edges % 4)
      {
      case 0:
        cur_delta[0] = cur_delta[0] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        cur_delta[1] = cur_delta[1] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        cur_delta[2] = cur_delta[2] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        break;
      case 1:
        cur_delta[0] = cur_delta[0] * ref_coords[(i_edges / 4 + 1) % 3] * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        cur_delta[1] = cur_delta[1] * ref_coords[(i_edges / 4 + 1) % 3] * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        cur_delta[2] = cur_delta[2] * ref_coords[(i_edges / 4 + 1) % 3] * (1 - ref_coords[(i_edges / 4 + 2) % 3]);
        break;
      case 2:
        cur_delta[0] = cur_delta[0] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * ref_coords[(i_edges / 4 + 2) % 3];
        cur_delta[1] = cur_delta[1] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * ref_coords[(i_edges / 4 + 2) % 3];
        cur_delta[2] = cur_delta[2] * (1 - ref_coords[(i_edges / 4 + 1) % 3]) * ref_coords[(i_edges / 4 + 2) % 3];        
        break;
      case 3:
        cur_delta[0] = cur_delta[0] * ref_coords[(i_edges / 4 + 1) % 3] * ref_coords[(i_edges / 4 + 2) % 3];
        cur_delta[1] = cur_delta[1] * ref_coords[(i_edges / 4 + 1) % 3] * ref_coords[(i_edges / 4 + 2) % 3];
        cur_delta[2] = cur_delta[2] * ref_coords[(i_edges / 4 + 1) % 3] * ref_coords[(i_edges / 4 + 2) % 3];        
        break;
      }
      
      /* Add edge displacements to out_coords */
      out_coords[0] += cur_delta[0];
      out_coords[1] += cur_delta[1];
      out_coords[2] += cur_delta[2];
    }
  }
  
  /* Check each face for geometry. Currently, only hexes with 6 faces are supported. */
  for (int i_faces = 0; i_faces < 6; ++i_faces)
  {
    if (faces[i_faces] >= 0)
    {      
      /* Interpolate coordinates between face vertices
      *   
      *               5 ---------------- 7
      *              /|                 /|
      *             / |      f5        / |
      *            /  |     (top)     /  |  <-f3  
      *           /   |              /   | (back)   z y
      *          4 ---------------- 5    |          |/
      *          | f0 |             | f1 |          x-- x
      *          |    2 ------------|--- 3
      *          |   /              |   /
      *   f2->   |  /               |  /
      *  (front) | /        f4      | /
      *          |/      (bottom)   |/
      *          0 ---------------- 1
      *                 
      *  ref_coords[(i_faces / 2 + 1) % 3] and ref_coords[(i_faces / 2 + 2) % 3] always returns a ref_coord perpendicular to face i_faces
      *  
      *  The order in which the interpolations are computed is not consistent with the vertex indices.
      *  Therefore, we need to switch the vertices 1 and 2 when calculating face 2 and 3. This is done with [i_faces / 2 == 1 ? 2 : 1]
      */
      for (int i_coord_dim = 0; i_coord_dim < 3; ++i_coord_dim)
      {
        interpolated_coords[i_coord_dim] = ((1 - ref_coords[(i_faces / 2 + 1) % 3]) * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_faces][0] * 3 + i_coord_dim]
                                + ref_coords[(i_faces / 2 + 1) % 3] * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_faces][i_faces / 2 == 1 ? 2 : 1] * 3 + i_coord_dim])
                                * (1 - ref_coords[(i_faces / 2 + 2) % 3]);
        interpolated_coords[i_coord_dim] += ((1 - ref_coords[(i_faces / 2 + 1) % 3]) * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_faces][i_faces / 2 == 1 ? 1 : 2] * 3 + i_coord_dim]
                                + ref_coords[(i_faces / 2 + 1) % 3] * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i_faces][3] * 3 + i_coord_dim])
                                * ref_coords[(i_faces / 2 + 2) % 3];
      }

      /* Interpolate parameters between face vertices. Same procedure as above. */
      const double *parameters = (double *) t8_cmesh_get_attribute(cmesh, t8_get_package_id (),
                                                      T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + i_faces,
                                                      gtreeid);
      for (int i_param_dim = 0; i_param_dim < 2; ++i_param_dim)
      {
        param[i_param_dim] = (parameters[0 + i_param_dim] * (1 - ref_coords[(i_faces / 2 + 1) % 3])
                  + parameters[(i_faces / 2 == 1 ? 4 : 2) + i_param_dim] * (ref_coords[(i_faces / 2 + 1) % 3]))
                  * (1 - ref_coords[(i_faces / 2 + 2) % 3]);
        param[i_param_dim] += (parameters[(i_faces / 2 == 1 ? 2 : 4) + i_param_dim] * (1 - ref_coords[(i_faces / 2 + 1) % 3])
                    +parameters[6 + i_param_dim] * (ref_coords[(i_faces / 2 + 1) % 3]))
                    *(ref_coords[(i_faces / 2 + 2) % 3]);
      }
      
      /* Check if calculated parameters are valid */
      T8_ASSERT(!t8_global_occ_surface[faces[i_faces]].IsNull());
      double u1, u2, v1, v2;
      t8_global_occ_surface[faces[i_faces]]->Bounds(u1, u2, v1, v2);
      T8_ASSERT(((u1 <= param[0] && param[0] <= u2) || (u2 <= param[0] && param[0] <= u1))
                && ((v1 <= param[1] && param[1] <= v2) || (v2 <= param[1] && param[1] <= v1)));

      /* Compute displacement between vertex interpolation and surface evaluation with interpolated parameters */
      t8_global_occ_surface[faces[i_faces]]->D0(param[0], param[1], pnt);

      /* Compute delta between geometry and interpolated coords, scale them with the appropriate ref_coord 
       * and add them to the out_coords*/
      if (i_faces % 2 == 0)
      {
        out_coords[0] += (pnt.X() - interpolated_coords[0]) * (1 - ref_coords[i_faces / 2]);
        out_coords[1] += (pnt.Y() - interpolated_coords[1]) * (1 - ref_coords[i_faces / 2]);
        out_coords[2] += (pnt.Z() - interpolated_coords[2]) * (1 - ref_coords[i_faces / 2]);
      }
      else
      {
        out_coords[0] += (pnt.X() - interpolated_coords[0]) * ref_coords[i_faces / 2];
        out_coords[1] += (pnt.Y() - interpolated_coords[1]) * ref_coords[i_faces / 2];
        out_coords[2] += (pnt.Z() - interpolated_coords[2]) * ref_coords[i_faces / 2];
      }
    }
  }
}

void
t8_geometry_occ::t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                          t8_gloidx_t gtreeid,
                                          const double
                                          *ref_coords,
                                          double *jacobian_out) const
{
  SC_ABORT ("Not implemented.");
}

inline void
t8_geometry_occ::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid)
{
  t8_geometry_w_vertices::t8_geom_load_tree_data(cmesh, gtreeid);
}

#endif /* T8_WITH_OCC */
