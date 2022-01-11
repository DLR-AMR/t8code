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
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>
#include <t8_eclass.h>
#include <t8_geometry/t8_geometry_helpers.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopoDS.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#if T8_WITH_OCC

const int
t8_edge_vertex_to_tree_vertex[T8_ECLASS_MAX_EDGES][2] =
{
  {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {4, 6}, {1, 3}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}    /* hex */
};

const int
t8_edge_to_face[T8_ECLASS_MAX_EDGES][2] = 
{
  {2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {0, 5}, {1, 4}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}    /* hex */
};

const int
t8_face_edge_to_tree_edge[T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D]
{
  {8, 10, 4, 5}, {9, 11, 6, 7}, {8, 9, 0, 2}, {10, 11, 1, 3}, {4, 6, 0, 1}, {5, 7, 2, 3}
};

t8_geometry_occ::t8_geometry_occ (int dim, const char *fileprefix, const char *name_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;
  
  BRep_Builder        builder;
  std::string         current_file(fileprefix);
  std::ifstream       is(current_file + ".brep");
  BRepTools::Read(occ_shape, is, builder);
  is.close();
  if(occ_shape.IsNull())
  {
    SC_ABORTF("Could not read brep file or brep file contains no shape \n");
  }
  TopExp::MapShapes(occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes(occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes(occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors(occ_shape, TopAbs_VERTEX, TopAbs_EDGE, occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors(occ_shape, TopAbs_EDGE, TopAbs_FACE, occ_shape_edge2face_map);
}

t8_geometry_occ::t8_geometry_occ (int dim, const TopoDS_Shape occ_shape, const char *name_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;
  if(occ_shape.IsNull())
  {
    SC_ABORTF("Shape is null. \n");
  }
  TopExp::MapShapes(occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes(occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes(occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors(occ_shape, TopAbs_VERTEX, TopAbs_EDGE, occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors(occ_shape, TopAbs_EDGE, TopAbs_FACE, occ_shape_edge2face_map);
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
                                                           T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY,
                                                           gtreeid);
  const int *faces = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                           T8_CMESH_OCC_FACE_ATTRIBUTE_KEY,
                                                           gtreeid);
  
  /* Compute coordinates via trilinear interpolation */
  t8_geom_compute_linear_geometry (active_tree_class,
                                   active_tree_vertices, ref_coords,
                                   out_coords);

  double interpolated_coords[3], param[2], cur_delta[3];
  gp_Pnt pnt;
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;
  
  /* Check each edge for a geometry. Currently, only hexes with 12 edges are supported. */
  for (int i_edges = 0; i_edges < 12; ++i_edges)
  {
    /* We have to check for curves as well as surfaces. Linked curves are stored in the first half of the array, surfaces in the second. 
    *  If a curve is connected to this edge we have to also check, if a surface is connected to at least one of the two adjacent faces.
    *  If there is a face present the edge is only used to get the right parameters while evaluating the surface and we can ignore it. */
    if ((edges[i_edges] > 0 && faces[t8_edge_to_face[i_edges][0]] == 0 && faces[t8_edge_to_face[i_edges][1]] == 0) || 
        edges[i_edges + 12] > 0)
    {
      /* Check if only a surface or a curve is present. Abort if both is true. */ 
      T8_ASSERT(!(edges[i_edges] > 0) != !(edges[i_edges + 12] > 0));

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
                                                                  T8_CMESH_OCC_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edges,
                                                                  gtreeid);
      /* Edges have only one parameter u, surfaces have two, u and v.
       * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
      if (edges[i_edges] > 0)
      {
        /* Linear interpolation between parameters */
        param[0] = (1 - ref_coords[i_edges / 4]) * parameters[0]
                    + ref_coords[i_edges / 4] * parameters[1];

        T8_ASSERT (edges[i_edges] <= occ_shape_edge_map.Size());
        curve = BRep_Tool::Curve(TopoDS::Edge(occ_shape_edge_map.FindKey(edges[i_edges])), first, last);
        
        /* Check if curve are valid */
        T8_ASSERT(!curve.IsNull());

        /* Calculate point on curve with interpolated parameters. */
        curve->D0(param[0], pnt);
      }
      else
      {
        /* Linear interpolation between parameters */
        param[0] = (1 - ref_coords[i_edges / 4]) * parameters[0]
                    + ref_coords[i_edges / 4] * parameters[2];
        param[1] = (1 - ref_coords[i_edges / 4]) * parameters[1]
                    + ref_coords[i_edges / 4] * parameters[3];
        
        T8_ASSERT (edges[i_edges + 12] <= occ_shape_face_map.Size());
        surface = BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(edges[i_edges + 12])));

        /* Check if surface is valid */
        T8_ASSERT(!surface.IsNull());

        /* Compute point on surface with interpolated parameters */
        surface->D0(param[0], param[1], pnt);
      }
      
      /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
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
    if (faces[i_faces] > 0)
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
                                                      T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + i_faces,
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

      T8_ASSERT (faces[i_faces] <= occ_shape_face_map.Size());
      surface = BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(faces[i_faces])));
      
      /* Check if surface is valid */
      T8_ASSERT(!surface.IsNull());

      /* Compute point on surface with interpolated parameters */
      surface->D0(param[0], param[1], pnt);

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

gp_Pnt
t8_geometry_occ::t8_geom_get_occ_point (const int index) const
{
  T8_ASSERT (index <= occ_shape_vertex_map.Size());
  return BRep_Tool::Pnt(TopoDS::Vertex(occ_shape_vertex_map.FindKey(index)));
}

Handle_Geom_Curve
t8_geometry_occ::t8_geom_get_occ_curve (const int index) const
{
  T8_ASSERT (index <= occ_shape_edge_map.Size());
  Standard_Real first, last;
  return BRep_Tool::Curve(TopoDS::Edge(occ_shape_edge_map.FindKey(index)), first, last);
}

Handle_Geom_Surface
t8_geometry_occ::t8_geom_get_occ_surface (const int index) const
{
  T8_ASSERT (index <= occ_shape_face_map.Size());
  return BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(index)));
}

TopTools_IndexedMapOfShape 
t8_geometry_occ::t8_geom_get_occ_shape_vertex_map() const
{
  return occ_shape_vertex_map;
}

TopTools_IndexedMapOfShape 
t8_geometry_occ::t8_geom_get_occ_shape_edge_map() const
{
  return occ_shape_edge_map;
}

TopTools_IndexedMapOfShape
t8_geometry_occ::t8_geom_get_occ_shape_face_map() const
{
  return occ_shape_face_map;
}

int
t8_geometry_occ::t8_geom_get_common_edge (const int vertex1_index, const int vertex2_index) const
{
  auto collection1 = occ_shape_vertex2edge_map.FindFromIndex(vertex1_index);
  auto collection2 = occ_shape_vertex2edge_map.FindFromIndex(vertex2_index);

  for (auto edge1 = collection1.begin(); edge1 != collection1.end(); ++edge1)
  {
    for (auto edge2 = collection2.begin(); edge2 != collection2.end(); ++edge2)
    {
      if (edge1->IsEqual(*edge2))
      {
        return occ_shape_edge2face_map.FindIndex(*edge1);
      }
    }
  }
  return 0;
}

int
t8_geometry_occ::t8_geom_get_common_face (const int edge1_index, const int edge2_index) const
{
  auto collection1 = occ_shape_edge2face_map.FindFromIndex(edge1_index);
  auto collection2 = occ_shape_edge2face_map.FindFromIndex(edge2_index);

  for (auto face1 = collection1.begin(); face1 != collection1.end(); ++face1)
  {
    for (auto face2 = collection2.begin(); face2 != collection2.end(); ++face2)
    {
      if (face1->IsEqual(*face2))
      {
        return occ_shape_face_map.FindIndex(*face1);
      }
    }
  }
  return 0;
}

int
t8_geometry_occ::t8_geom_is_vertex_on_edge (const int vertex_index, const int edge_index) const
{
  auto collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  return collection.Contains(occ_shape_edge_map.FindKey(edge_index));
}

int
t8_geometry_occ::t8_geom_is_edge_on_face (const int edge_index, const int face_index) const
{
  auto collection = occ_shape_edge2face_map.FindFromIndex(edge_index);
  return collection.Contains(occ_shape_face_map.FindKey(face_index));
}

int
t8_geometry_occ::t8_geom_is_vertex_on_face (const int vertex_index, const int face_index) const
{
  auto edge_collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  for (auto edge = edge_collection.begin(); edge != edge_collection.end(); ++edge)
  {
    auto face_collection = occ_shape_edge2face_map.FindFromKey(*edge);
    if (face_collection.Contains(occ_shape_face_map.FindKey(face_index)))
    {
      return 1;
    }
  }
  return 0;
}

void t8_geometry_occ::t8_geom_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                                              const int edge_index, 
                                                              double* edge_param) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_vertex_on_edge(vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  *edge_param = BRep_Tool::Parameter(vertex, edge);
}

void t8_geometry_occ::t8_geom_get_parameters_of_vertex_on_face(const int vertex_index, 
                                                               const int face_index, 
                                                               double* face_params) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_vertex_on_face(vertex_index, face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  uv = BRep_Tool::Parameters(vertex, face);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}

void t8_geometry_occ::t8_geom_edge_parameter_to_face_parameters(const int edge_index, 
                                                                const int face_index, 
                                                                const double edge_param, 
                                                                double* face_params) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_edge_on_face(edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  Handle_Geom2d_Curve curve_on_surface  = BRep_Tool::CurveOnSurface(edge, face, first, last);
  curve_on_surface->D0(edge_param, uv);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_occ.h.
 * Create a new geometry with given dimension. */
t8_geometry_occ_c      *
t8_geometry_occ_new (int dimension, const char *fileprefix, const char *name_in)
{
  t8_geometry_occ *geom = new t8_geometry_occ (dimension, fileprefix, name_in);
  return (t8_geometry_occ_c *) geom;
}

void
t8_geometry_occ_destroy (t8_geometry_occ_c ** geom)
{
#ifdef T8_ENABLE_DEBUG
  t8_geometry_occ_c      *pgeom = *geom;
  T8_ASSERT (dynamic_cast < t8_geometry_occ * >(pgeom) != NULL);
#endif

  delete             *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();

#endif /* T8_WITH_OCC */
