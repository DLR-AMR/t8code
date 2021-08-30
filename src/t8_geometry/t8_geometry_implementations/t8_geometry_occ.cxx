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

  double interpolated_coords[3], param[2], edge_delta[3]= {0, 0, 0}, cur_edge_delta[3], surface_deltas[6][3];
  gp_Pnt pnt;
  
  /* Check each edge for geometry */
  for (int i = 0; i < 12; ++i)
  {
    if (edges[i] >= 0)
    {
      /* Interpolate coordinates between edge vertices. Due to the indices i of the edges, the edges point in
      * direction of ref_coord i % 3. Therefore, we can use ref_coords[i % 3] for the interpolation.              
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
      for (int j = 0; j < 3; ++j)
      {
        interpolated_coords[j] = (1 - ref_coords[i % 3]) * active_tree_vertices[t8_edge_vertex_to_tree_vertex[i][0] * 3 + j]
                                + ref_coords[i % 3] * active_tree_vertices[t8_edge_vertex_to_tree_vertex[i][1] * 3 + j];
      }

      /* Interpolate parameters between edge vertices. Same procedure as above. */
      const double *parameters = (double *) t8_cmesh_get_attribute(cmesh, t8_get_package_id (),
                                                                  T8_CMESH_OCC_CURVE_PARAMETERS_ATTRIBUTE_KEY + i,
                                                                  gtreeid);
      param[0] = (1 - ref_coords[i % 3]) * parameters[0]
                  + ref_coords[i % 3] * parameters[1];
      
      /* Check if calculated parameters are valid */
      T8_ASSERT(!t8_global_occ_curve[edges[i]].IsNull());
      double u1, u2;
      u1 = t8_global_occ_curve[edges[i]]->FirstParameter();
      u2 = t8_global_occ_curve[edges[i]]->LastParameter();
      T8_ASSERT((u1 <= param[0] && param[0] <= u2) || (u2 <= param[0] && param[0] <= u1));

      /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
      t8_global_occ_curve[edges[i]]->D0(param[0], pnt);
      cur_edge_delta[0] = pnt.X() - interpolated_coords[0];
      cur_edge_delta[1] = pnt.Y() - interpolated_coords[1];
      cur_edge_delta[2] = pnt.Z() - interpolated_coords[2];
      
      /* Multiply curve displacement with corresponding ref coords.
      *  The edges are indexed so that all edges which satisfy i % 4 == 0 have to multiplied with the inversed (1 - ref_coord) 
      *  coordinate. All edges which satisfy i % 4 == 1 have to multiplied with one inversed ref_coord and so forth...
      */
      switch (i % 4)
      {
      case 0:
        cur_edge_delta[0] = cur_edge_delta[0] * (1 - ref_coords[(i / 4 + 1) % 3]) * (1 - ref_coords[(i / 4 + 2) % 3]);
        cur_edge_delta[1] = cur_edge_delta[1] * (1 - ref_coords[(i / 4 + 1) % 3]) * (1 - ref_coords[(i / 4 + 2) % 3]);
        cur_edge_delta[2] = cur_edge_delta[2] * (1 - ref_coords[(i / 4 + 1) % 3]) * (1 - ref_coords[(i / 4 + 2) % 3]);
        break;
      case 1:
        cur_edge_delta[0] = cur_edge_delta[0] * ref_coords[(i / 4 + 1) % 3] * (1 - ref_coords[(i / 4 + 2) % 3]);
        cur_edge_delta[1] = cur_edge_delta[1] * ref_coords[(i / 4 + 1) % 3] * (1 - ref_coords[(i / 4 + 2) % 3]);
        cur_edge_delta[2] = cur_edge_delta[2] * ref_coords[(i / 4 + 1) % 3] * (1 - ref_coords[(i / 4 + 2) % 3]);
        break;
      case 2:
        cur_edge_delta[0] = cur_edge_delta[0] * (1 - ref_coords[(i / 4 + 1) % 3]) * ref_coords[(i / 4 + 2) % 3];
        cur_edge_delta[1] = cur_edge_delta[1] * (1 - ref_coords[(i / 4 + 1) % 3]) * ref_coords[(i / 4 + 2) % 3];
        cur_edge_delta[2] = cur_edge_delta[2] * (1 - ref_coords[(i / 4 + 1) % 3]) * ref_coords[(i / 4 + 2) % 3];        
        break;
      case 3:
        cur_edge_delta[0] = cur_edge_delta[0] * ref_coords[(i / 4 + 1) % 3] * ref_coords[(i / 4 + 2) % 3];
        cur_edge_delta[1] = cur_edge_delta[1] * ref_coords[(i / 4 + 1) % 3] * ref_coords[(i / 4 + 2) % 3];
        cur_edge_delta[2] = cur_edge_delta[2] * ref_coords[(i / 4 + 1) % 3] * ref_coords[(i / 4 + 2) % 3];        
        break;
      }
      
      /* Accumulate all edge displacements */
      edge_delta[0] += cur_edge_delta[0];
      edge_delta[1] += cur_edge_delta[1];
      edge_delta[2] += cur_edge_delta[2];
    }
  }
  
  /* Add edge displacements to out_coords */
  out_coords[0] += edge_delta[0];
  out_coords[1] += edge_delta[1];
  out_coords[2] += edge_delta[2];
  
  /* Check each face for geometry */
  for (int i = 0; i < 6; ++i)
  {
    if (faces[i] >= 0)
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
      *  ref_coords[(i / 2 + 1) % 3] and ref_coords[(i / 2 + 2) % 3] always returns a ref_coord perpendicular to face i
      *  
      *  The order in which the interpolations are computed is not consistent with the vertex indices.
      *  Therefore, we need to switch the vertices 1 and 2 when calculating face 2 and 3. This is done with [i / 2 == 1 ? 2 : 1]
      */
      for (int j = 0; j < 3; ++j)
      {
        interpolated_coords[j] = ((1 - ref_coords[(i / 2 + 1) % 3]) * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i][0] * 3 + j]
                                + ref_coords[(i / 2 + 1) % 3] * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i][i / 2 == 1 ? 2 : 1] * 3 + j])
                                * (1 - ref_coords[(i / 2 + 2) % 3]);
        interpolated_coords[j] += ((1 - ref_coords[(i / 2 + 1) % 3]) * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i][i / 2 == 1 ? 1 : 2] * 3 + j]
                                + ref_coords[(i / 2 + 1) % 3] * active_tree_vertices[t8_face_vertex_to_tree_vertex[T8_ECLASS_HEX][i][3] * 3 + j])
                                * ref_coords[(i / 2 + 2) % 3];
      }

      /* Interpolate parameters between face vertices. Same procedure as above. */
      const double *parameters = (double *) t8_cmesh_get_attribute(cmesh, t8_get_package_id (),
                                                      T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + i,
                                                      gtreeid);
      for (int j = 0; j < 2; ++j)
      {
        param[j] = (parameters[0 + j] * (1 - ref_coords[(i / 2 + 1) % 3])
                  + parameters[(i / 2 == 1 ? 4 : 2) + j] * (ref_coords[(i / 2 + 1) % 3]))
                  * (1 - ref_coords[(i / 2 + 2) % 3]);
        param[j] += (parameters[(i / 2 == 1 ? 2 : 4) + j] * (1 - ref_coords[(i / 2 + 1) % 3])
                    +parameters[6 + j] * (ref_coords[(i / 2 + 1) % 3]))
                    *(ref_coords[(i / 2 + 2) % 3]);
      }
      
      /* Check if calculated parameters are valid */
      T8_ASSERT(!t8_global_occ_surface[faces[i]].IsNull());
      double u1, u2, v1, v2;
      t8_global_occ_surface[faces[i]]->Bounds(u1, u2, v1, v2);
      T8_ASSERT(((u1 <= param[0] && param[0] <= u2) || (u2 <= param[0] && param[0] <= u1))
                && ((v1 <= param[1] && param[1] <= v2) || (v2 <= param[1] && param[1] <= v1)));

      /* Compute displacement between vertex interpolation and surface evaluation with interpolated parameters */
      t8_global_occ_surface[faces[i]]->D0(param[0], param[1], pnt);
      surface_deltas[i][0] = pnt.X() - interpolated_coords[0];
      surface_deltas[i][1] = pnt.Y() - interpolated_coords[1];
      surface_deltas[i][2] = pnt.Z() - interpolated_coords[2]; 
    }
  }  

  /* Multiply surface displacements with corresponding ref coords and add them to the out coords */
  for (int i = 0; i < 3; ++i)
  {
    if (faces[i * 2] >= 0)
    {
      for (int j = 0; j < 3; ++j)
      {
        out_coords[j] += surface_deltas[i * 2][j] * (1 - ref_coords[i]);
      }
    }
    if (faces[i * 2 + 1] >= 0)
    {
      for (int j = 0; j < 3; ++j)
      {
        out_coords[j] += surface_deltas[i * 2 + 1][j] * ref_coords[i];
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
  /* Set active id and eclass */
  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  active_tree = gtreeid;
  active_tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Load this trees vertices. */
  active_tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  T8_ASSERT (t8_eclass_to_dimension[active_tree_class] == dimension);

  /* Check whether we support this class */
  T8_ASSERT (active_tree_class == T8_ECLASS_HEX);
}

#endif /* T8_WITH_OCC */
