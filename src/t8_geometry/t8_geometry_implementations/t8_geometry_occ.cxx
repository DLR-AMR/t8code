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

#if T8_WITH_OCC

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
#include <Standard_Version.hxx>

#endif /* T8_WITH_OCC */

/* *INDENT-OFF* */
const int           t8_edge_vertex_to_tree_vertex[T8_ECLASS_MAX_EDGES][2] = {
  {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {4, 6}, 
  {1, 3}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}        /* hex */
};

const int           t8_edge_to_face[T8_ECLASS_MAX_EDGES][2] = {
  {2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {0, 5}, 
  {1, 4}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}        /* hex */
};

const int
t8_face_edge_to_tree_edge[T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D] = {
  {8, 10, 4, 5}, {9, 11, 6, 7}, {8, 9, 0, 2}, 
  {10, 11, 1, 3}, {4, 6, 0, 1}, {5, 7, 2, 3}        /* hex */
};
/* *INDENT-ON* */

#if T8_WITH_OCC

t8_geometry_occ::t8_geometry_occ (int dim, const char *fileprefix,
                                  const char *name_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;

  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  BRepTools::Read (occ_shape, is, builder);
  is.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape. "
               "The OCC file may be written with a newer OCC version. "
               "Linked OCC version: %s", OCC_VERSION_COMPLETE);
  }
  TopExp::MapShapes (occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes (occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes (occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_VERTEX, TopAbs_EDGE,
                                       occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_EDGE, TopAbs_FACE,
                                       occ_shape_edge2face_map);
}

t8_geometry_occ::t8_geometry_occ (int dim, const TopoDS_Shape occ_shape,
                                  const char *name_in)
{
  T8_ASSERT (0 <= dim && dim <= 3);

  name = name_in;
  dimension = dim;
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  TopExp::MapShapes (occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes (occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes (occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_VERTEX, TopAbs_EDGE,
                                       occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_EDGE, TopAbs_FACE,
                                       occ_shape_edge2face_map);
}

void
t8_geometry_occ::t8_geom_evaluate (t8_cmesh_t cmesh,
                                   t8_gloidx_t gtreeid,
                                   const double *ref_coords,
                                   double out_coords[3]) const
{
  /* Compute coordinates via trilinear interpolation */
  t8_geom_compute_linear_geometry (active_tree_class,
                                   active_tree_vertices, ref_coords,
                                   out_coords);

  const int           num_edges = t8_eclass_num_edges[active_tree_class];
  const int           num_faces = t8_eclass_num_faces[active_tree_class];
  double              interpolated_coords[3],
    interpolated_curve_param, interpolated_surface_params[2], cur_delta[3];
  gp_Pnt              pnt;
  double              interpolation_coeffs[3],
    temp_face_vertices[T8_ECLASS_MAX_CORNERS_2D * 3],
    temp_edge_vertices[2 * 3];
  Handle_Geom_Curve   curve;
  Handle_Geom_Surface surface;
  Standard_Real       first, last;

  /* Other tree classes are not implemented yet. */
  T8_ASSERT (active_tree_class == T8_ECLASS_HEX);

  /* Check each edge for a geometry. */
  for (int i_edges = 0; i_edges < num_edges; ++i_edges) {
    /* We have to check for curves as well as surfaces. Linked curves are stored 
     * in the first half of the array, surfaces in the second. 
     * If a curve is connected to this edge we have to also check, 
     * if a surface is connected to at least one of the two adjacent faces.
     * If there is a face present the edge is only used to get the right parameters 
     * while evaluating the surface and we can ignore it. */
    if (edges[i_edges] > 0 || edges[i_edges + num_edges] > 0) {
      /* Check if only a surface or a curve is present. Abort if both is true. */
      T8_ASSERT (!(edges[i_edges] > 0) != !(edges[i_edges + num_edges] > 0));

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
      const int           edge_direction = i_edges / 4;
      /* Save the edge vertices temporarily. */
      t8_geom_get_edge_vertices (active_tree_class,
                                 active_tree_vertices,
                                 i_edges, 3, temp_edge_vertices);
      /* Interpolate between them. */
      t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                    temp_edge_vertices,
                                    3, 1, interpolated_coords);

      /* Interpolate parameters between edge vertices. Same procedure as above. */
      const double       *parameters =
        (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                           T8_CMESH_OCC_EDGE_PARAMETERS_ATTRIBUTE_KEY
                                           + i_edges,
                                           gtreeid);
      T8_ASSERT (parameters != NULL);
      /* Edges have only one parameter u, surfaces have two, u and v.
       * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
      if (edges[i_edges] > 0) {
        /* Linear interpolation between parameters */
        t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                      parameters, 1, 1,
                                      &interpolated_curve_param);

        /* *INDENT-OFF* */
        T8_ASSERT (edges[i_edges] <= occ_shape_edge_map.Size ());
        curve =
          BRep_Tool::Curve (TopoDS:: Edge (occ_shape_edge_map.FindKey (edges[i_edges])), 
                                                                       first, last);
        /* *INDENT-ON* */

        /* Check if curve are valid */
        T8_ASSERT (!curve.IsNull ());

        /* Calculate point on curve with interpolated parameters. */
        curve->D0 (interpolated_curve_param, pnt);
      }
      else {
        /* Linear interpolation between parameters */
        t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                      parameters, 2, 1,
                                      interpolated_surface_params);

        T8_ASSERT (edges[i_edges + num_edges] <= occ_shape_face_map.Size ());
        surface =
          BRep_Tool::Surface (TopoDS::Face
                              (occ_shape_face_map.FindKey
                               (edges[i_edges + num_edges])));

        /* Check if surface is valid */
        T8_ASSERT (!surface.IsNull ());

        /* Compute point on surface with interpolated parameters */
        surface->D0 (interpolated_surface_params[0],
                     interpolated_surface_params[1], pnt);
      }

      /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
      cur_delta[0] = pnt.X () - interpolated_coords[0];
      cur_delta[1] = pnt.Y () - interpolated_coords[1];
      cur_delta[2] = pnt.Z () - interpolated_coords[2];

      /* Multiply curve displacement with corresponding ref coords.
       * The edges are indexed so that all edges which satisfy i_edges % 4 == 0 
       * have to multiplied with the inversed (1 - ref_coord) 
       * coordinate. All edges which satisfy i_edges % 4 == 1 have to multiplied with one 
       * inversed ref_coord and so forth... */

      /* *INDENT-OFF* */
      double scaling_factor = 0;
      switch (i_edges % 4) {
      case 0:
        scaling_factor = (1 - ref_coords[(edge_direction + 1) % 3]) 
                       * (1 - ref_coords[(edge_direction + 2) % 3]);
        break;
      case 1:
        scaling_factor = ref_coords[(edge_direction + 1) % 3] 
                       * (1 - ref_coords[(edge_direction + 2) % 3]);
        break;
      case 2:
        scaling_factor = (1 - ref_coords[(edge_direction + 1) % 3]) 
                       * ref_coords[(edge_direction + 2) % 3];
        break;
      case 3:
        scaling_factor = ref_coords[(edge_direction + 1) % 3] 
                       * ref_coords[(edge_direction + 2) % 3];
        break;
      }
      /* *INDENT-ON* */

      /* Add edge displacements to out_coords */
      out_coords[0] += cur_delta[0] * scaling_factor;
      out_coords[1] += cur_delta[1] * scaling_factor;
      out_coords[2] += cur_delta[2] * scaling_factor;
    }
  }

  /* Iterate over each face to calculate the displacements generated by each face */
  for (int i_faces = 0; i_faces < num_faces; ++i_faces) {
    /* Check if face has a linked surface */
    if (faces[i_faces] > 0) {
      /* Allocate some variables and save the normal direction of the face and the face vertices 
       * in a separate array for later usage. */
      const int           face_normal_direction = i_faces / 2;
      t8_geom_get_face_vertices (T8_ECLASS_HEX,
                                 active_tree_vertices,
                                 i_faces, 3, temp_face_vertices);

      double              face_displacement_from_edges[3] = { 0 },
        surface_parameter_displacement_from_edges[2] = { 0 },
        surface_parameters_from_curve[2] = { 0 };

      /* Retrieve surface parameters of nodes */
      const double       *surface_parameters =
        (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                           T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY
                                           + i_faces,
                                           gtreeid);
      T8_ASSERT (surface_parameters != NULL);
      /* Iterate over each edge of face */
      for (int i_face_edge = 0; i_face_edge < 4; ++i_face_edge) {
        /* Check if curve is present */
        if (edges[t8_face_edge_to_tree_edge[i_faces][i_face_edge]] > 0) {
          /* Calculating some indices */
          const int           edge_direction =
            t8_face_edge_to_tree_edge[i_faces][i_face_edge] / 4;
          int                 orthogonal_direction_of_edge_on_face = 0;
          switch (edge_direction + face_normal_direction) {
          case 1:
            orthogonal_direction_of_edge_on_face = 2;
            break;
          case 3:
            orthogonal_direction_of_edge_on_face = 0;
            break;
          case 2:
            orthogonal_direction_of_edge_on_face = 1;
            break;
          }

          /* Retrieve parameters of nodes und curve */
          const double       *curve_parameters =
            (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                               T8_CMESH_OCC_EDGE_PARAMETERS_ATTRIBUTE_KEY
                                               +
                                               t8_face_edge_to_tree_edge
                                               [i_faces][i_face_edge],
                                               gtreeid);
          T8_ASSERT (curve_parameters != NULL);
          /* Interpolate linearly between the parameters of the two nodes on the curve */
          t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                        curve_parameters, 1, 1,
                                        &interpolated_curve_param);
          /* Do the same interpolation but with the surface parameters of the same two nodes as above */
          double              interpolated_surface_parameters_on_edge[2];
          double              edge_parameters_on_face[4];
          t8_geom_get_face_vertices ((t8_eclass_t)
                                     t8_eclass_face_types[active_tree_class]
                                     [i_faces], surface_parameters,
                                     i_face_edge, 2, edge_parameters_on_face);
          t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                        edge_parameters_on_face, 2, 1,
                                        interpolated_surface_parameters_on_edge);
          /* Do the same interpolation but this time between the coordinates of the same two nodes as above */
          double              interpolated_edge_coordinates[3];
          double              edge_vertices_on_face[6];
          t8_geom_get_face_vertices ((t8_eclass_t)
                                     t8_eclass_face_types[active_tree_class]
                                     [i_faces], temp_face_vertices,
                                     i_face_edge, 3, edge_vertices_on_face);
          t8_geom_linear_interpolation (&ref_coords[edge_direction],
                                        edge_vertices_on_face, 3, 1,
                                        interpolated_edge_coordinates);

          /* Retrieve the curve of the edge */
          T8_ASSERT (edges[t8_face_edge_to_tree_edge[i_faces][i_face_edge]] <=
                     occ_shape_edge_map.Size ());
          curve =
            BRep_Tool::Curve (TopoDS::Edge (occ_shape_edge_map.FindKey (edges
                                                                        [t8_face_edge_to_tree_edge
                                                                         [i_faces]
                                                                         [i_face_edge]])), first, last);
          /* Check if curve is valid */
          T8_ASSERT (!curve.IsNull ());
          /* Calculate point on curve with interpolated parameters */
          curve->D0 (interpolated_curve_param, pnt);

          /* Calculate the displacement generated by the presence of the curve */
          if (i_face_edge % 2 == 0) {
            face_displacement_from_edges[0]
              += (1 - ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.X () - interpolated_edge_coordinates[0]);
            face_displacement_from_edges[1]
              += (1 - ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.Y () - interpolated_edge_coordinates[1]);
            face_displacement_from_edges[2]
              += (1 - ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.Z () - interpolated_edge_coordinates[2]);
          }
          else {
            face_displacement_from_edges[0]
              += (ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.X () - interpolated_edge_coordinates[0]);
            face_displacement_from_edges[1]
              += (ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.Y () - interpolated_edge_coordinates[1]);
            face_displacement_from_edges[2]
              += (ref_coords[orthogonal_direction_of_edge_on_face])
              * (pnt.Z () - interpolated_edge_coordinates[2]);
          }

          /* Convert the interpolated parameter of the curve into the corresponding parameters on the surface */
          t8_geometry_occ::t8_geom_edge_parameter_to_face_parameters (edges
                                                                      [t8_face_edge_to_tree_edge
                                                                       [i_faces]
                                                                       [i_face_edge]], faces[i_faces], interpolated_curve_param, surface_parameters_from_curve);

          /* Calculate the displacement between the interpolated parameters on the surface 
           * and the parameters on the surface converted from the parameter of the curve
           * and scale them with the corresponding ref coord */
          if (i_face_edge % 2 == 0) {
            for (int dim = 0; dim < 2; ++dim) {
              surface_parameter_displacement_from_edges[dim]
                += (1 - ref_coords[orthogonal_direction_of_edge_on_face])
                * (surface_parameters_from_curve[dim]
                   - interpolated_surface_parameters_on_edge[dim]);
            }
          }
          else {
            for (int dim = 0; dim < 2; ++dim) {
              surface_parameter_displacement_from_edges[dim]
                += ref_coords[orthogonal_direction_of_edge_on_face]
                * (surface_parameters_from_curve[dim]
                   - interpolated_surface_parameters_on_edge[dim]);
            }
          }
        }
      }

      /* Do a bilinear interpolation between the nodes of the face and the parameters of the face */
      interpolation_coeffs[0] = ref_coords[(face_normal_direction + 1) % 3];
      interpolation_coeffs[1] = ref_coords[(face_normal_direction + 2) % 3];
      /* The normal vectors of the faces 0, 1, 4, 5 of a hex point in the same direction
       * as the corresponding axis (normal(f0) = (1, 0, 0)). The faces 2 and 3 are oriented
       * in the other direction (normal(f2) = (0, -1, 0)). Therefore we have to switch two 
       * vertices to change the orientation. Here we switch vertex 1 and 2. */
      double              temp_surface_parameters[8];
      memcpy (temp_surface_parameters, surface_parameters,
              sizeof (double) * 8);

      if (face_normal_direction == 1) {
        double              temp;
        /* Switch vertex coordinates */
        for (int dim = 0; dim < 3; ++dim) {
          temp = temp_face_vertices[1 * 3 + dim];
          temp_face_vertices[1 * 3 + dim] = temp_face_vertices[2 * 3 + dim];
          temp_face_vertices[2 * 3 + dim] = temp;
        }
        /* Switch vertex parameters */
        for (int dim = 0; dim < 2; ++dim) {
          temp = temp_surface_parameters[1 * 2 + dim];
          temp_surface_parameters[1 * 2 + dim] =
            temp_surface_parameters[2 * 2 + dim];
          temp_surface_parameters[2 * 2 + dim] = temp;
        }
      }
      /* The actual bilinear interpolations */
      t8_geom_linear_interpolation (interpolation_coeffs,
                                    temp_face_vertices,
                                    3, 2, interpolated_coords);
      t8_geom_linear_interpolation (interpolation_coeffs,
                                    temp_surface_parameters,
                                    2, 2, interpolated_surface_params);

      for (int dim = 0; dim < 3; ++dim) {
        /* Correct the interpolated coordinates with the displacement generated by the linked edges */
        interpolated_coords[dim] += face_displacement_from_edges[dim];
      }

      for (int dim = 0; dim < 2; ++dim) {
        /* Correct the interpolated parameters with the parameter displacement generated by the edges */
        interpolated_surface_params[dim] +=
          surface_parameter_displacement_from_edges[dim];
      }

      /* *INDENT-OFF* */
      /* Retrieve the surface of the edge */
      T8_ASSERT (faces[i_faces] <= occ_shape_face_map.Size ());
      surface =
        BRep_Tool::Surface (TopoDS::Face (occ_shape_face_map.FindKey (faces[i_faces])));
      /* *INDENT-ON* */

      /* Check if surface is valid */
      T8_ASSERT (!surface.IsNull ());

      /* Compute point on surface with interpolated parameters */
      surface->D0 (interpolated_surface_params[0],
                   interpolated_surface_params[1], pnt);

      /* Compute the displacement between surface and interpolated coords, scale them with the appropriate ref_coord 
       * and add them to the out_coords. */
      if (i_faces % 2 == 0) {
        out_coords[0]
          += (pnt.X () - interpolated_coords[0])
          * (1 - ref_coords[face_normal_direction]);
        out_coords[1]
          += (pnt.Y () - interpolated_coords[1])
          * (1 - ref_coords[face_normal_direction]);
        out_coords[2]
          += (pnt.Z () - interpolated_coords[2])
          * (1 - ref_coords[face_normal_direction]);
      }
      else {
        out_coords[0]
          += (pnt.X () - interpolated_coords[0])
          * ref_coords[face_normal_direction];
        out_coords[1]
          += (pnt.Y () - interpolated_coords[1])
          * ref_coords[face_normal_direction];
        out_coords[2]
          += (pnt.Z () - interpolated_coords[2])
          * ref_coords[face_normal_direction];
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
  double              h = 1e-9;
  double              in1[3], in2[3];
  double              out1[3], out2[3];
  for (int dim = 0; dim < 3; ++dim) {
    memcpy (in1, ref_coords, sizeof (double) * 3);
    memcpy (in2, ref_coords, sizeof (double) * 3);

    if (ref_coords[dim] < h) {
      in2[dim] += ref_coords[dim] + h;
    }
    else if (ref_coords[dim] > 1 - h) {
      in1[dim] -= h;
    }
    else {
      in1[dim] -= 0.5 * h;
      in2[dim] += 0.5 * h;
    }
    t8_geometry_occ::t8_geom_evaluate (cmesh, gtreeid, in1, out1);
    t8_geometry_occ::t8_geom_evaluate (cmesh, gtreeid, in2, out2);
    for (int dim2 = 0; dim2 < 3; ++dim2) {
      jacobian_out[dim * 3 + dim2] = (out2[dim2] - out1[dim2]) / h;
    }
  }
}

/* Our indent skript has huge problems with c++ */
/* *INDENT-OFF* */
inline void
t8_geometry_occ::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                         t8_gloidx_t gtreeid)
{
  t8_geometry_w_vertices::t8_geom_load_tree_data(cmesh, gtreeid);
  edges = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY,
                                                gtreeid);
  faces = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                T8_CMESH_OCC_FACE_ATTRIBUTE_KEY,
                                                gtreeid);
  T8_ASSERT(edges != NULL);
  T8_ASSERT(faces != NULL);
}

const gp_Pnt
t8_geometry_occ::t8_geom_get_occ_point (const int index) const
{
  T8_ASSERT (index <= occ_shape_vertex_map.Size());
  return BRep_Tool::Pnt(TopoDS::Vertex(occ_shape_vertex_map.FindKey(index)));
}

const Handle_Geom_Curve
t8_geometry_occ::t8_geom_get_occ_curve (const int index) const
{
  T8_ASSERT (index <= occ_shape_edge_map.Size());
  Standard_Real first, last;
  return BRep_Tool::Curve(TopoDS::Edge(occ_shape_edge_map.FindKey(index)), 
                          first, last);
}

const Handle_Geom_Surface
t8_geometry_occ::t8_geom_get_occ_surface (const int index) const
{
  T8_ASSERT (index <= occ_shape_face_map.Size());
  return BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(index)));
}

const TopTools_IndexedMapOfShape 
t8_geometry_occ::t8_geom_get_occ_shape_vertex_map() const
{
  return occ_shape_vertex_map;
}

const TopTools_IndexedMapOfShape 
t8_geometry_occ::t8_geom_get_occ_shape_edge_map() const
{
  return occ_shape_edge_map;
}

const TopTools_IndexedMapOfShape
t8_geometry_occ::t8_geom_get_occ_shape_face_map() const
{
  return occ_shape_face_map;
}

int
t8_geometry_occ::t8_geom_get_common_edge (const int vertex1_index, 
                                          const int vertex2_index) const
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
t8_geometry_occ::t8_geom_get_common_face (const int edge1_index, 
                                          const int edge2_index) const
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
t8_geometry_occ::t8_geom_is_vertex_on_edge (const int vertex_index, 
                                            const int edge_index) const
{
  auto collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  return collection.Contains(occ_shape_edge_map.FindKey(edge_index));
}

int
t8_geometry_occ::t8_geom_is_edge_on_face (const int edge_index, 
                                          const int face_index) const
{
  auto collection = occ_shape_edge2face_map.FindFromIndex(edge_index);
  return collection.Contains(occ_shape_face_map.FindKey(face_index));
}

int
t8_geometry_occ::t8_geom_is_vertex_on_face (const int vertex_index, 
                                            const int face_index) const
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

void 
t8_geometry_occ::t8_geom_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                                         const int edge_index, 
                                                         double* edge_param) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_vertex_on_edge(vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  *edge_param = BRep_Tool::Parameter(vertex, edge);
}

void 
t8_geometry_occ::t8_geom_get_parameters_of_vertex_on_face(const int vertex_index, 
                                                          const int face_index, 
                                                          double* face_params) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_vertex_on_face(vertex_index, 
                                                       face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  uv = BRep_Tool::Parameters(vertex, face);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}

void 
t8_geometry_occ::t8_geom_edge_parameter_to_face_parameters(const int edge_index, 
                                                           const int face_index, 
                                                           const double edge_param, 
                                                           double* face_params) const
{
  T8_ASSERT(t8_geometry_occ::t8_geom_is_edge_on_face(edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  Handle_Geom2d_Curve curve_on_surface  = BRep_Tool::CurveOnSurface(edge, face, 
                                                                    first, last);
  curve_on_surface->D0(edge_param, uv);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}
/* *INDENT-ON* */

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_occ.h.
 * Create a new geometry with given dimension. */
t8_geometry_occ_c  *
t8_geometry_occ_new (int dimension, const char *fileprefix,
                     const char *name_in)
{
  t8_geometry_occ    *geom =
    new t8_geometry_occ (dimension, fileprefix, name_in);
  return (t8_geometry_occ_c *) geom;
}

void
t8_geometry_occ_destroy (t8_geometry_occ_c ** geom)
{
#ifdef T8_ENABLE_DEBUG
  t8_geometry_occ_c  *pgeom = *geom;
  T8_ASSERT (dynamic_cast < t8_geometry_occ * >(pgeom) != NULL);
#endif

  delete             *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();

#endif /* T8_WITH_OCC */
