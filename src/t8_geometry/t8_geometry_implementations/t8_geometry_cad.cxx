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
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.h>
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
#include <Standard_Version.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>

/* The lookup table contains the coordinate of each edge of a tetrahedron,
 * which is used for the interpolation.
 * For example: On edges 0, 1, and 2 the x coordinates is used to interpolate */
const int t8_interpolation_coefficient_tet_edge[6] = { 0, 0, 0, 2, 2, 1 };
/* The lookup table contains the coordinates of each face of a tetrahedron.
 * For example: face 0 is described by coordinates z and y. */
const int t8_face_ref_coords_tet[4][2] = { { 2, 1 }, { 0, 1 }, { 0, 1 }, { 0, 2 } };

t8_geometry_cad::t8_geometry_cad (std::string fileprefix, std::string name_in): t8_geometry_with_vertices (name_in)
{
  BRep_Builder builder;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  if (is.is_open () == false) {
    SC_ABORTF ("Cannot find the file %s.brep.\n", fileprefix.c_str ());
  }
  BRepTools::Read (cad_shape, is, builder);
  is.close ();
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape. "
               "The cad file may be written with a newer cad version. "
               "Linked cad version: %s",
               OCC_VERSION_COMPLETE);
  }
  TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
  TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
  TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_VERTEX, TopAbs_EDGE, cad_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_EDGE, TopAbs_FACE, cad_shape_edge2face_map);
}

t8_geometry_cad::t8_geometry_cad (const TopoDS_Shape cad_shape, std::string name_in)
  : t8_geometry_with_vertices (name_in)
{
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
  TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
  TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_VERTEX, TopAbs_EDGE, cad_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_EDGE, TopAbs_FACE, cad_shape_edge2face_map);
}

t8_geometry_cad::t8_geometry_cad (): t8_geometry_with_vertices ("t8_geom_cad")
{
  cad_shape.Nullify ();
}

void
t8_geometry_cad::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                   const size_t num_coords, double *out_coords) const
{
  switch (active_tree_class) {
  case T8_ECLASS_TRIANGLE:
    t8_geometry_cad::t8_geom_evaluate_cad_tri (cmesh, gtreeid, ref_coords, num_coords, out_coords);
    break;
  case T8_ECLASS_QUAD:
    t8_geometry_cad::t8_geom_evaluate_cad_quad (cmesh, gtreeid, ref_coords, num_coords, out_coords);
    break;
  case T8_ECLASS_HEX:
    t8_geometry_cad::t8_geom_evaluate_cad_hex (cmesh, gtreeid, ref_coords, num_coords, out_coords);
    break;
  case T8_ECLASS_PRISM:
    t8_geometry_cad::t8_geom_evaluate_cad_prism (cmesh, gtreeid, ref_coords, num_coords, out_coords);
    break;
  case T8_ECLASS_TET:
    t8_geometry_cad::t8_geom_evaluate_cad_tet (cmesh, gtreeid, ref_coords, num_coords, out_coords);
    break;
  default:
    SC_ABORTF ("Error: Curved cad geometry for element type %s not yet implemented. \n",
               t8_eclass_to_string[active_tree_class]);
  }
}

void
t8_geometry_cad::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                            const size_t num_coords, double *jacobian_out) const
{
  double h = 1e-9;
  double in1[3], in2[3];
  double out1[3], out2[3];
  for (size_t icoord = 0; icoord < num_coords; icoord++) {

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
      t8_geometry_cad::t8_geom_evaluate (cmesh, gtreeid, in1, 1, out1);
      t8_geometry_cad::t8_geom_evaluate (cmesh, gtreeid, in2, 1, out2);
      for (int dim2 = 0; dim2 < 3; ++dim2) {
        jacobian_out[9 * icoord + dim * 3 + dim2] = (out2[dim2] - out1[dim2]) / h;
      }
    }
  }
}

inline void
t8_geometry_cad::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_geometry_with_vertices::t8_geom_load_tree_data (cmesh, gtreeid);
  edges = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, ltreeid);
  faces = (const int *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, ltreeid);
  T8_ASSERT (edges != NULL);
  T8_ASSERT (faces != NULL);
}

void
t8_geometry_cad::t8_geom_evaluate_cad_tri (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                           const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (active_tree_class == T8_ECLASS_TRIANGLE);

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  const int num_edges = t8_eclass_num_edges[active_tree_class];
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;
  gp_Pnt pnt;
  double displacement;
  double scaling_factor;
  double scaled_displacement;
  /* Allocate storage for later usage. Storage depents on size of the batch. */
  double *ref_intersection = T8_ALLOC (double, 2 * num_coords);
  double *glob_intersection = T8_ALLOC (double, 3 * num_coords);
  double interpolated_curve_parameter;
  double *converted_edge_surface_parameters = T8_ALLOC (double, 2 * num_coords);
  double *interpolated_edge_surface_parameters = T8_ALLOC (double, 2 * num_coords);
  double *interpolated_surface_parameters = T8_ALLOC (double, 2 * num_coords);

  /*
   * Reference Space    |      Global Space
   *      [0,1]^2       |         [x,y]^2
   *                    |
   *              2     |            2
   *            / |     |           / \
   *           /  |     |          /    \
   *          /   |     |         /       \              o ref_coords in reference space
   *         /    |     |        /          \            @ glob_ref_coords in global space
   *        E1     E0   |       E1           E0
   *       /      |     |      /                \
   *      /       o     |     /                  @
   *     /        |     |    /                   |   y
   *    /         |     |   /                   /    |
   *   0----E2----1     |  0---------E2--------1     x--x
   * 
   */

  /* Linear mapping from ref_coords to out_coords for each reference point */
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

  /* Check if face has a linked geometry */
  if (*faces > 0) {
#if T8_ENABLE_DEBUG
    for (int i_edge = 0; i_edge < num_edges; i_edge++) {
      /* If face carries a surface, edges can't carry surfaces too */
      T8_ASSERT (edges[i_edge + num_edges] == 0);
    }
#endif /* T8_ENABLE_DEBUG */
    /* Retrieve surface parameters */
    const double *face_parameters = (double *) t8_cmesh_get_attribute (
      cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY, ltreeid);
    T8_ASSERT (face_parameters != NULL);

    /* Retrieve surface_parameter for each reference point in global space by triangular interpolation from ref_coords to global space */
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_2d = i_coord * 2;
      /* We use the out coords as buffer for the surface parameters to
       * avoid allocating extra dynamic memory */
      t8_geom_triangular_interpolation (ref_coords + offset_2d, face_parameters, 2, 2,
                                        interpolated_surface_parameters + offset_2d);
    }
    /* Check every edge and search for edge displacement */
    for (int i_edge = 0; i_edge < num_edges; i_edge++) {
      if (edges[i_edge] > 0) {
        /* Calculate the intersections of straight lines from the opposite vertex of the current edge,
         * through each reference point, onto the current edge */
        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          t8_geom_get_ref_intersection (i_edge, ref_coords + offset_2d, ref_intersection + offset_2d);
        }
        /* Converting ref_intersections to global_intersections by interpolation */
        t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_intersection, num_coords,
                                         glob_intersection);

        /* Get parameters of the current edge if the edge is curved */
        const double *edge_parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (edge_parameters != NULL);

        const int num_face_nodes = t8_eclass_num_vertices[active_tree_class];
        /* Calculate the curve parameter at the intersection point for each reference coordinate */
        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          /* If the current edge is edge 0, we use the y coordinate of the ref intersection.
           * For all other edges we can use the x coordinate. */
          t8_geom_linear_interpolation (&ref_intersection[(i_edge == 0) + offset_2d], edge_parameters, 1, 1,
                                        &interpolated_curve_parameter);
          /* Convert the interpolated edge parameter of each reference point to surface parameters */
          t8_geometry_cad::t8_geom_edge_parameter_to_face_parameters (edges[i_edge], *faces, num_face_nodes,
                                                                      interpolated_curve_parameter, face_parameters,
                                                                      converted_edge_surface_parameters + offset_2d);
        }

        double edge_surface_parameters[4];
        /* Get the surface parameters at the vertices of the current edge */
        t8_geom_get_face_vertices (active_tree_class, face_parameters, i_edge, 2, edge_surface_parameters);
        /* Interpolate between the surface parameters of the current edge with the ref_intersection of each reference point */
        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          t8_geom_linear_interpolation (&ref_intersection[(i_edge == 0) + offset_2d], edge_surface_parameters, 2, 1,
                                        interpolated_edge_surface_parameters + offset_2d);
        }

        /* Determine the scaling factor for each reference point by calculating the distances from the opposite vertex
         * to the glob_intersections and to the reference points */
        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          const int offset_3d = i_coord * 3;
          scaling_factor = t8_geom_get_triangle_scaling_factor (i_edge, active_tree_vertices,
                                                                glob_intersection + offset_3d, out_coords + offset_3d);

          /* Calculate the parameter displacement for each reference point and add it to the surface parameters */
          for (int dim = 0; dim < 2; ++dim) {
            displacement = converted_edge_surface_parameters[dim + offset_2d]
                           - interpolated_edge_surface_parameters[dim + offset_2d];
            scaled_displacement = displacement * scaling_factor;
            interpolated_surface_parameters[dim + offset_2d] += scaled_displacement;
          }
          /* Retrieve surface */
          T8_ASSERT (*faces <= cad_shape_face_map.Size ());
          surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (*faces)));
          /* Check if surface is valid */
          T8_ASSERT (!surface.IsNull ());

          /* Evaluate surface and save result */
          surface->D0 (interpolated_surface_parameters[offset_2d], interpolated_surface_parameters[offset_2d + 1], pnt);

          for (int dim = 0; dim < 3; ++dim) {
            out_coords[dim + offset_3d] = pnt.Coord (dim + 1);
          }
        }
      }
      /* Retrieve surface */
      T8_ASSERT (*faces <= cad_shape_face_map.Size ());
      surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (*faces)));
      /* Check if surface is valid */
      T8_ASSERT (!surface.IsNull ());

      /* Evaluate surface and save result */
      surface->D0 (interpolated_surface_parameters[0], interpolated_surface_parameters[1], pnt);

      for (int dim = 0; dim < 3; ++dim) {
        out_coords[dim] = pnt.Coord (dim + 1);
      }
    }
  }
  else {
    /* No surface is linked to the face. We therefore search for edge displacements
     * for the triangle cell and add them.
     * Iterate over each edge and check if it is linked. */
    for (int i_edge = 0; i_edge < num_edges; i_edge++) {
      if (edges[i_edge] > 0 || edges[i_edge + num_edges] > 0) {
        /* Get parameters of the current edge if the edge is curved */
        const double *parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (parameters != NULL);

        /* Calculate the intersections of straight lines from the opposite vertex of the current edge,
         * through each reference point, onto the current edge */
        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          t8_geom_get_ref_intersection (i_edge, ref_coords + offset_2d, ref_intersection + offset_2d);
        }
        /* Converting ref_intersections to global_intersections */
        t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_intersection, num_coords,
                                         glob_intersection);

        for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
          const int offset_2d = i_coord * 2;
          const int offset_3d = i_coord * 3;
          if (edges[i_edge] > 0) {
            /* Interpolate between the curve parameters of the current edge with the ref_intersection of each reference point */
            t8_geom_linear_interpolation (&ref_intersection[(i_edge == 0) + offset_2d], parameters, 1, 1,
                                          &interpolated_curve_parameter);
            /* Retrieve curve */
            T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());
            curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);
            /* Check if curve is valid */
            T8_ASSERT (!curve.IsNull ());

            /* Calculate point on curve with interpolated parameters. */
            curve->D0 (interpolated_curve_parameter, pnt);
          }
          else {
            /* Interpolate between the surface parameters of the current edge with the ref_intersection of each reference point */
            t8_geom_linear_interpolation (&ref_intersection[(i_edge == 0) + offset_2d], parameters, 2, 1,
                                          interpolated_surface_parameters + offset_2d);
            T8_ASSERT (edges[i_edge + num_edges] <= cad_shape_face_map.Size ());
            surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (edges[i_edge + num_edges])));
            /* Check if surface is valid */
            T8_ASSERT (!surface.IsNull ());

            /* Compute point on surface with interpolated parameters */
            surface->D0 (interpolated_surface_parameters[offset_2d], interpolated_surface_parameters[offset_2d + 1],
                         pnt);
          }
          /* Determine the scaling factor by calculating the distances from the opposite vertex
          * to the glob_intersection and to the reference point */
          scaling_factor = t8_geom_get_triangle_scaling_factor (i_edge, active_tree_vertices,
                                                                glob_intersection + offset_3d, out_coords + offset_3d);

          /* Calculate displacement between points on curve and point on linear curve.
           * Then scale it and add the scaled displacement to the result. */
          for (int dim = 0; dim < 3; ++dim) {
            displacement = pnt.Coord (dim + 1) - glob_intersection[dim + offset_3d];
            scaled_displacement = displacement * scaling_factor;
            out_coords[dim + offset_3d] += scaled_displacement;
          }
        }
      }
    }
  }
  T8_FREE (ref_intersection);
  T8_FREE (glob_intersection);
  T8_FREE (converted_edge_surface_parameters);
  T8_FREE (interpolated_edge_surface_parameters);
  T8_FREE (interpolated_surface_parameters);
}

void
t8_geometry_cad::t8_geom_evaluate_cad_quad (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                            const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (active_tree_class == T8_ECLASS_QUAD);

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  const int num_edges = t8_eclass_num_edges[active_tree_class];
  gp_Pnt pnt;
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;

  /* Check if face has a linked geometry */
  if (*faces > 0) {
#if T8_ENABLE_DEBUG
    /* Check, that edges do not carry a surface as well */
    for (int i_edge = 0; i_edge < num_edges; ++i_edge) {
      T8_ASSERT (edges[i_edge + num_edges] <= 0);
    }
#endif /* T8_ENABLE_DEBUG */

    /* Retrieve surface parameters */
    const double *face_parameters = (double *) t8_cmesh_get_attribute (
      cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY, ltreeid);
    T8_ASSERT (face_parameters != NULL);

    /* Interpolate between surface parameters */
    for (size_t coord = 0; coord < num_coords; ++coord) {
      const int offset_3d = coord * 3;
      const int offset_2d = coord * 2;
      /* We use the out coords as buffer for the surface parameters to
       * avoid allocating dynamic memory */
      t8_geom_linear_interpolation (ref_coords + offset_2d, face_parameters, 2, 2, out_coords + offset_3d);
    }

    /* Iterate over each edge to search for parameter displacements */
    for (int i_edge = 0; i_edge < num_edges; ++i_edge) {
      if (edges[i_edge] > 0) {
        /* The edges of a quad point in direction of ref_coord (1 - i_edge >> 1).
         *
         *     2 -------E3------- 3
         *     |                  |
         *     |                  |
         *    E0                 E1
         *     |                  |
         *     |                  |    y
         *     |                  |    |
         *     0 -------E2------- 1    x-- x
         *        
         */
        const int edge_orthogonal_direction = (i_edge >> 1);
        const int edge_direction = 1 - edge_orthogonal_direction;
        const int num_face_nodes = t8_eclass_num_vertices[active_tree_class];
        double temp_edge_parameters[2];
        double temp_face_parameters[2];
        /* Retrieve edge parameters and interpolate */
        const double *edge_parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (edge_parameters != NULL);
        T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());

        curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);

        /* Check if curve is valid */
        T8_ASSERT (!curve.IsNull ());

        for (size_t coord = 0; coord < num_coords; ++coord) {
          const int offset_3d = coord * 3;
          const int offset_2d = coord * 2;
          /* Interpolate between curve parameters and surface parameters of the same nodes */
          t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_2d], edge_parameters, 1, 1,
                                        temp_edge_parameters);

          /* Convert curve parameter to surface parameters */
          t8_geometry_cad::t8_geom_edge_parameter_to_face_parameters (
            edges[i_edge], *faces, num_face_nodes, temp_edge_parameters[0], face_parameters, temp_edge_parameters);

          /* Interpolate between the surface parameters of the current edge */
          double edge_surface_parameters[4];
          t8_geom_get_face_vertices (active_tree_class, face_parameters, i_edge, 2, edge_surface_parameters);
          t8_geom_linear_interpolation (&ref_coords[edge_direction], edge_surface_parameters, 2, 1,
                                        temp_face_parameters);

          /* Calculate parameter displacement and add it to the surface parameters */
          for (int dim = 0; dim < 2; ++dim) {
            const double displacement = temp_edge_parameters[dim] - temp_face_parameters[dim];
            if (i_edge & 1) {
              out_coords[dim + offset_3d] += displacement * (1 - ref_coords[edge_orthogonal_direction]);
            }
            else {
              out_coords[dim + offset_3d] += displacement * ref_coords[edge_orthogonal_direction];
            }
          }
        }
      }
    }

    /* Retrieve surface */
    T8_ASSERT (*faces <= cad_shape_face_map.Size ());
    surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (*faces)));

    /* Check if surface is valid */
    T8_ASSERT (!surface.IsNull ());

    /* Evaluate surface and save result */
    for (size_t coord = 0; coord < num_coords; ++coord) {
      const int offset_3d = coord * 3;
      surface->D0 (out_coords[offset_3d], out_coords[offset_3d + 1], pnt);
      for (int dim = 0; dim < 3; ++dim) {
        out_coords[dim + offset_3d] = pnt.Coord (dim + 1);
      }
    }
  }
  else {
    /* No surface is linked to the face. We therefore use a bilinear
     * interpolation for the quad cell and add the edge displacements.
     * Iterate over each edge and check if it is linked. */
    for (size_t coord = 0; coord < num_coords; ++coord) {
      const int offset_3d = coord * 3;
      const int offset_2d = coord * 2;
      t8_geom_linear_interpolation (ref_coords + offset_2d, active_tree_vertices, 3, 2, out_coords + offset_3d);
    }
    for (int i_edge = 0; i_edge < num_edges; ++i_edge) {
      if (edges[i_edge] > 0 || edges[i_edge + num_edges] > 0) {
        /* An edge can only be linked to a curve or a surface, not both */
        T8_ASSERT (!(edges[i_edge] > 0) || !(edges[i_edge + num_edges] > 0));

        /* The edges of a quad point in direction of ref_coord (1 - i_edge >> 1).
         *
         *     2 -------E3------- 3
         *     |                  |
         *     |                  |
         *    E0                 E1
         *     |                  |
         *     |                  |    y
         *     |                  |    |
         *     0 -------E2------- 1    x-- x
         *        
         */
        const int edge_orthogonal_direction = (i_edge >> 1);
        const int edge_direction = 1 - edge_orthogonal_direction;
        double temp_edge_vertices[2 * 3];
        double temp_parameters[2];
        double temp_coords[3];

        /* Save the edge vertices temporarily. */
        t8_geom_get_face_vertices (active_tree_class, active_tree_vertices, i_edge, 3, temp_edge_vertices);

        /* Get pointer to edge parameters */
        const double *parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (parameters != NULL);
        /* Curves have only one parameter u, surfaces have two, u and v.
         * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
        if (edges[i_edge] > 0) {
          /* Get curve */
          T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());
          /* Infinite indent loop */
          /* *INDENT-OFF* */
          curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);

          /* Check if curve are valid */
          T8_ASSERT (!curve.IsNull ());
          for (size_t coord = 0; coord < num_coords; ++coord) {
            const int offset_3d = coord * 3;
            const int offset_2d = coord * 2;
            /* Linear interpolation between edge vertices */
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_2d], temp_edge_vertices, 3, 1,
                                          temp_coords);

            /* Linear interpolation between parameters */
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_2d], parameters, 1, 1, temp_parameters);
            /* Calculate point on curve with interpolated parameters. */
            curve->D0 (temp_parameters[0], pnt);

            /* Calculate scaling factor for edge */
            for (int dim = 0; dim < 3; ++dim) {
              const double displacement = pnt.Coord (dim + 1) - temp_coords[dim];

              if (i_edge % 2 == 0) {
                out_coords[dim + offset_3d] += displacement * (1 - ref_coords[edge_orthogonal_direction + offset_2d]);
              }
              else {
                out_coords[dim + offset_3d] += displacement * ref_coords[edge_orthogonal_direction + offset_2d];
              }
            }
          }
        }
        else {
          /* Get surface */
          T8_ASSERT (edges[i_edge + num_edges] <= cad_shape_face_map.Size ());
          surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (edges[i_edge + num_edges])));

          /* Check if surface is valid */
          T8_ASSERT (!surface.IsNull ());

          for (size_t coord = 0; coord < num_coords; ++coord) {
            const int offset_3d = coord * 3;
            const int offset_2d = coord * 2;
            /* Linear interpolation between edge vertices */
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_2d], temp_edge_vertices, 3, 1,
                                          temp_coords);

            /* Linear interpolation between parameters */
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_2d], parameters, 2, 1, temp_parameters);
            /* Calculate point on sirface with interpolated parameters. */
            surface->D0 (temp_parameters[0 + offset_2d], temp_parameters[1 + offset_2d], pnt);

            /* Calculate scaling factor for edge */
            for (int dim = 0; dim < 3; ++dim) {
              const double displacement = pnt.Coord (dim + 1) - temp_coords[dim];

              if (i_edge % 2 == 0) {
                out_coords[dim + offset_3d] += displacement * (1 - ref_coords[edge_orthogonal_direction + offset_2d]);
              }
              else {
                out_coords[dim + offset_3d] += displacement * ref_coords[edge_orthogonal_direction + offset_2d];
              }
            }
          }
        }
      }
    }
  }
}

void
t8_geometry_cad::t8_geom_evaluate_cad_tet (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                           const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (active_tree_class == T8_ECLASS_TET);

  /* Compute linear coordinates via triangular interpolation (barycentric coordinates) */
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  const int num_edges = t8_eclass_num_edges[active_tree_class];
  const int num_faces = t8_eclass_num_faces[active_tree_class];
  gp_Pnt pnt;
  double temp_edge_vertices[2 * 3], temp_face_vertices[T8_ECLASS_MAX_CORNERS_2D * 3], interpolated_curve_param,
    interpolated_surface_params[2], cur_delta[3];
  double interpolated_surface_parameters[2], interpolated_coords[3];
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;

  for (size_t coord = 0; coord < num_coords; ++coord) {
    const int offset_3d = coord * 3;

    /* Check each edge for a geometry. */
    for (int i_edge = 0; i_edge < num_edges; ++i_edge) {

      /* We have to check for curves as well as surfaces. Linked curves are stored 
      * in the first half of the array, surfaces in the second. 
      * If a curve is connected to this edge we have to also check, 
      * if a surface is connected to at least one of the two adjacent faces. */
      if (edges[i_edge] > 0 || edges[i_edge + num_edges] > 0) {

        /* Check if only a surface or a curve is present. Abort if both is true. */
        T8_ASSERT (!(edges[i_edge] > 0 && edges[i_edge + num_edges] > 0));

        /* 
        *             _0
        *          _- / \
        *       E0   /   \
        *    _-     E1    \
        *  1 --__  /       E2
        *   \    -/-__      \
        *   E3   /    E4__   \
        *     \ /         --__\
        *      2------ E5 -----3
        */

        /* Save the vertices of the current edge */
        t8_geom_get_edge_vertices (active_tree_class, active_tree_vertices, i_edge, 3, temp_edge_vertices);

        /* Get the interpolation coefficients for the current edge */
        double interpolation_coeff = ref_coords[t8_interpolation_coefficient_tet_edge[i_edge] + offset_3d];

        /* Interpolate between the edge vertices with the interpolation_coeff of the current edge */
        t8_geom_linear_interpolation (&interpolation_coeff, temp_edge_vertices, 3, 1, interpolated_coords);

        /* Retrieve edge parameters of the current edge */
        const double *parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (parameters != NULL);

        /* Interpolate between the parameters of the current edge. Same procedure as above.
        * Curves have only one parameter u, surfaces have two, u and v.
        * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
        if (edges[i_edge] > 0) { /* Check for linked curves */

          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&interpolation_coeff, parameters, 1, 1, &interpolated_curve_param);
          T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());

          /* Retrieve the curve and check if curve is valid */
          curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);
          T8_ASSERT (!curve.IsNull ());

          /* Calculate point on curve with the interpolated parameter */
          curve->D0 (interpolated_curve_param, pnt);
        }
        else { /* Check for linked surfaces */

          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&interpolation_coeff, parameters, 2, 1, interpolated_surface_params);
          T8_ASSERT (edges[i_edge + num_edges] <= cad_shape_face_map.Size ());

          /* Retrieve the surface and check if surface is valid */
          surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (edges[i_edge + num_edges])));
          T8_ASSERT (!surface.IsNull ());

          /* Compute point on surface with interpolated parameters */
          surface->D0 (interpolated_surface_params[0], interpolated_surface_params[1], pnt);
        }

        /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
        cur_delta[0] = pnt.X () - interpolated_coords[0];
        cur_delta[1] = pnt.Y () - interpolated_coords[1];
        cur_delta[2] = pnt.Z () - interpolated_coords[2];

        /* Save the scaling factor for both neighbouring faces of the current edge.
        * The scaling factors scale the edge displacement orthogonal along the adjacent faces. */
        double scaling_factor_neigh_face_1 = t8_geom_get_scaling_factor_of_edge_on_face_tet (
          i_edge, t8_edge_to_face[active_tree_class][i_edge][0], ref_coords + offset_3d);
        double scaling_factor_neigh_face_2 = t8_geom_get_scaling_factor_of_edge_on_face_tet (
          i_edge, t8_edge_to_face[active_tree_class][i_edge][1], ref_coords + offset_3d);
        double scaling_factor = scaling_factor_neigh_face_1 * scaling_factor_neigh_face_2;

        /* out_coord correction with scaling */
        out_coords[0 + offset_3d] += cur_delta[0] * (scaling_factor * scaling_factor);
        out_coords[1 + offset_3d] += cur_delta[1] * (scaling_factor * scaling_factor);
        out_coords[2 + offset_3d] += cur_delta[2] * (scaling_factor * scaling_factor);
      }
    }

    /* Iterate over each face to calculate the displacements generated */
    for (int i_faces = 0; i_faces < num_faces; ++i_faces) {

      /* Check if face has a linked surface */
      if (faces[i_faces] > 0) {
        /* 
        * The faces of a tetrahedron are numerated in such a way,
        * that face X is the opposing face to the corner node X. 
        * 
        *             _0
        *          _- / \
        *       _-   /   \
        *    _-  F3 /     \
        *  1 --__  / F2    \
        *   \    -/-__ F1   \
        *    \   /    --__   \
        *     \ / F0      --__\
        *      2---------------3
        */

        /* Allocate some variables and save the face vertices for later usage */
        t8_geom_get_face_vertices (active_tree_class, active_tree_vertices, i_faces, 3, temp_face_vertices);
        double face_displacement_from_edges[3] = { 0 };

        /* Save the face intersection of a ray passing trough the reference coordinate
        * and the opposite vertex of the face */
        double face_intersection[3] = { 0 };
        t8_geom_get_tet_face_intersection (i_faces, ref_coords + offset_3d, face_intersection);

        /* Turn 3D face_intersection into 2D coordinates on current face 
        * for parameter interpolation */
        double face_intersection_2d[2];
        face_intersection_2d[0] = face_intersection[t8_face_ref_coords_tet[i_faces][0]];
        face_intersection_2d[1] = face_intersection[t8_face_ref_coords_tet[i_faces][1]];

        /* Retrieve surface_parameters of the linked face */
        const double *surface_parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + i_faces, ltreeid);
        T8_ASSERT (surface_parameters != NULL);

        /* Iterate over each edge of face to search for additional edge displacement */
        for (int i_face_edge = 0; i_face_edge < 3; ++i_face_edge) {
          /* Save the tree edge */
          int i_tree_edge = t8_face_edge_to_tree_edge[active_tree_class][i_faces][i_face_edge];

          /* Check if curve is present */
          if (edges[i_tree_edge] > 0) {

            /* Retrieve parameters of nodes on curve */
            const double *curve_parameters = (double *) t8_cmesh_get_attribute (
              cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_tree_edge, ltreeid);
            T8_ASSERT (curve_parameters != NULL);

            /* Get the interpolation coefficients for the current edge */
            const double *interpolation_coeff = &face_intersection[t8_interpolation_coefficient_tet_edge[i_tree_edge]];

            /* Interpolate linearly between the parameters of the two nodes on the curve */
            t8_geom_linear_interpolation (interpolation_coeff, curve_parameters, 1, 1, &interpolated_curve_param);

            /* Do the same interpolation but this time between the coordinates of the same two nodes as above */
            double interpolated_edge_coordinates[3];
            double edge_vertices_on_face[6];
            t8_geom_get_edge_vertices (active_tree_class, active_tree_vertices, i_tree_edge, 3, edge_vertices_on_face);
            t8_geom_linear_interpolation (interpolation_coeff, edge_vertices_on_face, 3, 1,
                                          interpolated_edge_coordinates);

            /* Retrieve the curve of the edge and check if it is valid */
            T8_ASSERT (edges[i_tree_edge] <= cad_shape_edge_map.Size ());
            curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_tree_edge])), first, last);
            T8_ASSERT (!curve.IsNull ());

            /* Calculate point on curve with interpolated parameter */
            curve->D0 (interpolated_curve_param, pnt);

            /* Calculate the same scaling factors for the neighbouring faces
            * as in the evaluation of the edges of tetrahedral tree */
            double scaling_factor
              = t8_geom_get_scaling_factor_of_edge_on_face_tet (i_tree_edge, i_faces, face_intersection);

            /* Calculate the displacement generated by the presence of the curve */
            for (int dim = 0; dim < 3; ++dim) {
              face_displacement_from_edges[dim]
                += (scaling_factor * scaling_factor) * (pnt.Coord (dim + 1) - interpolated_edge_coordinates[dim]);
            }
          }
        }

        /* Interpolate to get the interpolated_coords on the current face in global space
        * and the interpolated_surface_parameters on the current face */
        t8_geom_triangular_interpolation (face_intersection, active_tree_vertices, 3, 3, interpolated_coords);
        t8_geom_triangular_interpolation (face_intersection_2d, surface_parameters, 2, 2,
                                          interpolated_surface_parameters);

        for (int dim = 0; dim < 3; ++dim) {
          /* Correct the interpolated coordinates with the displacement generated by the linked edges */
          interpolated_coords[dim] += face_displacement_from_edges[dim];
        }

        /* Retrieve the surface and check if it is valid */
        T8_ASSERT (faces[i_faces] <= cad_shape_face_map.Size ());
        surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (faces[i_faces])));
        T8_ASSERT (!surface.IsNull ());

        /* Compute point on surface with interpolated surface parameters */
        surface->D0 (interpolated_surface_parameters[0], interpolated_surface_parameters[1], pnt);

        /* Compute the scaling factor. The scaling happens along the straight from
        * the opposite vertex of the face to the face_intersection. */
        double dist_ref_coords;
        double dist_face_intersection;

        /* Save the opposite vertex of the face in reference space.
        * Reminder: Opposite vertex of a face has the same index as the face. */
        const double *ref_opposite_vertex = t8_element_corner_ref_coords[active_tree_class][i_faces];

        dist_ref_coords = sqrt (
          (ref_opposite_vertex[0] - ref_coords[0 + offset_3d]) * (ref_opposite_vertex[0] - ref_coords[0 + offset_3d])
          + (ref_opposite_vertex[1] - ref_coords[1 + offset_3d]) * (ref_opposite_vertex[1] - ref_coords[1 + offset_3d])
          + (ref_opposite_vertex[2] - ref_coords[2 + offset_3d])
              * (ref_opposite_vertex[2] - ref_coords[2 + offset_3d]));
        dist_face_intersection
          = sqrt ((ref_opposite_vertex[0] - face_intersection[0]) * (ref_opposite_vertex[0] - face_intersection[0])
                  + (ref_opposite_vertex[1] - face_intersection[1]) * (ref_opposite_vertex[1] - face_intersection[1])
                  + (ref_opposite_vertex[2] - face_intersection[2]) * (ref_opposite_vertex[2] - face_intersection[2]));
        double scaling_factor = dist_ref_coords / dist_face_intersection;

        /* Compute the displacement between surface and interpolated coords and scale it with the scaling factor */
        out_coords[0 + offset_3d] += (pnt.X () - interpolated_coords[0]) * (scaling_factor * scaling_factor);
        out_coords[1 + offset_3d] += (pnt.Y () - interpolated_coords[1]) * (scaling_factor * scaling_factor);
        out_coords[2 + offset_3d] += (pnt.Z () - interpolated_coords[2]) * (scaling_factor * scaling_factor);
      }
    }
  }
}

void
t8_geometry_cad::t8_geom_evaluate_cad_hex (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                           const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (active_tree_class == T8_ECLASS_HEX);

  /* Compute coordinates via trilinear interpolation */
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  const int num_edges = t8_eclass_num_edges[active_tree_class];
  const int num_faces = t8_eclass_num_faces[active_tree_class];
  double *interpolated_coords = T8_ALLOC (double, 3 * num_coords);
  double interpolated_curve_param, interpolated_surface_params[2], cur_delta[3];
  gp_Pnt pnt;
  double interpolation_coeffs[2], temp_face_vertices[T8_ECLASS_MAX_CORNERS_2D * 3], temp_edge_vertices[2 * 3];
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;

  for (size_t coord = 0; coord < num_coords; ++coord) {
    const int offset_3d = coord * 3;

    /* Check each edge for a geometry. */
    for (int i_edge = 0; i_edge < num_edges; ++i_edge) {
      /* We have to check for curves as well as surfaces. Linked curves are stored 
      * in the first half of the array, surfaces in the second. 
      * If a curve is connected to this edge we have to also check, 
      * if a surface is connected to at least one of the two adjacent faces. */
      if (edges[i_edge] > 0 || edges[i_edge + num_edges] > 0) {
        /* Check if only a surface or a curve is present. Abort if both is true. */
        T8_ASSERT (!(edges[i_edge] > 0) != !(edges[i_edge + num_edges] > 0));

        /* Interpolate coordinates between edge vertices. Due to the indices i_edge of the edges, the edges point in
        * direction of ref_coord i_edge >> 2. Therefore, we can use ref_coords[i_edge >> 2] for the interpolation.
        *          6 -------E3------- 7
        *         /|                 /|
        *       E6 |               E7 |
        *       / E10              / E11
        *      /   |              /   |          z y
        *     4 -------E2------- 5    |          |/
        *     |    |             |    |          x-- x
        *     |    2 -------E1---|--- 3
        *    E8   /             E9   /
        *     |  E4              |  E5
        *     | /                | /
        *     |/                 |/
        *     0 -------E0------- 1
        *        
        */
        const int edge_direction = i_edge / 4;
        /* Save the edge vertices temporarily. */
        t8_geom_get_edge_vertices (active_tree_class, active_tree_vertices, i_edge, 3, temp_edge_vertices);
        /* Interpolate between them. */
        t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], temp_edge_vertices, 3, 1,
                                      interpolated_coords + offset_3d);
        /* Interpolate parameters between edge vertices. Same procedure as above. */
        const double *parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (parameters != NULL);

        /* Curves have only one parameter u, surfaces have two, u and v.
        * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
        if (edges[i_edge] > 0) {
          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], parameters, 1, 1,
                                        &interpolated_curve_param);

          T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());
          curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);

          /* Check if curve are valid */
          T8_ASSERT (!curve.IsNull ());

          /* Calculate point on curve with interpolated parameters. */
          curve->D0 (interpolated_curve_param, pnt);
        }
        else {
          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], parameters, 2, 1,
                                        interpolated_surface_params);

          T8_ASSERT (edges[i_edge + num_edges] <= cad_shape_face_map.Size ());
          surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (edges[i_edge + num_edges])));

          /* Check if surface is valid */
          T8_ASSERT (!surface.IsNull ());

          /* Compute point on surface with interpolated parameters */
          surface->D0 (interpolated_surface_params[0], interpolated_surface_params[1], pnt);
        }

        /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
        cur_delta[0] = pnt.X () - interpolated_coords[offset_3d];
        cur_delta[1] = pnt.Y () - interpolated_coords[offset_3d + 1];
        cur_delta[2] = pnt.Z () - interpolated_coords[offset_3d + 2];

        /* Multiply curve displacement with corresponding ref coords.
        * The edges are indexed so that all edges which satisfy i_edge % 4 == 0 
        * have to multiplied with the inversed (1 - ref_coord) 
        * coordinate. All edges which satisfy i_edge % 4 == 1 have to multiplied with one 
        * inversed ref_coord and so forth...
        * An exception are edge 5 and 6, which have to be switched because they do not follow that rule.
        * Edges which are located at ref_coord[i] = 0 have to be multiplied with (1 - ref_coord[i]) and if the
        * edge is located at ref_coord[i] = 1 it has to be multiplied with ref_coord[i]. */

        double scaling_factor = 0;
        const int temp_i_edge = i_edge == 5 ? 6 : i_edge == 6 ? 5 : i_edge;
        switch (temp_i_edge % 4) {
        case 0:
          scaling_factor = (1 - ref_coords[((edge_direction + 1) % 3) + offset_3d])
                           * (1 - ref_coords[((edge_direction + 2) % 3) + offset_3d]);
          break;
        case 1:
          scaling_factor = ref_coords[((edge_direction + 1) % 3) + offset_3d]
                           * (1 - ref_coords[((edge_direction + 2) % 3) + offset_3d]);
          break;
        case 2:
          scaling_factor = (1 - ref_coords[((edge_direction + 1) % 3) + offset_3d])
                           * ref_coords[((edge_direction + 2) % 3) + offset_3d];
          break;
        case 3:
          scaling_factor
            = ref_coords[((edge_direction + 1) % 3) + offset_3d] * ref_coords[((edge_direction + 2) % 3) + offset_3d];
          break;
        }

        /* Add edge displacements to out_coords */
        out_coords[offset_3d + 0] += cur_delta[0] * scaling_factor;
        out_coords[offset_3d + 1] += cur_delta[1] * scaling_factor;
        out_coords[offset_3d + 2] += cur_delta[2] * scaling_factor;
      }
    }

    /* Iterate over each face to calculate the displacements generated by each face */
    for (int i_faces = 0; i_faces < num_faces; ++i_faces) {
      /* Check if face has a linked surface */
      if (faces[i_faces] > 0) {
        /* Allocate some variables and save the normal direction of the face and the face vertices 
        * in a separate array for later usage. */
        const int face_normal_direction = i_faces / 2;
        t8_geom_get_face_vertices (T8_ECLASS_HEX, active_tree_vertices, i_faces, 3, temp_face_vertices);

        /* Retrieve surface parameters of nodes */
        const double *surface_parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + i_faces, ltreeid);
        T8_ASSERT (surface_parameters != NULL);

        double face_displacement_from_edges[3] = { 0 };
        double surface_parameter_displacement_from_edges[2] = { 0 };
        double surface_parameters_from_curve[2];

        /* Iterate over each edge of face */
        for (int i_face_edge = 0; i_face_edge < 4; ++i_face_edge) {
          /* Check if curve is present */
          if (edges[t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge]] > 0) {
            /* Calculating some indices */
            const int edge_direction = t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge] / 4;
            int orthogonal_direction_of_edge_on_face = 0;
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
            /* Convert global tree id to local tree id, for receiving cmesh attributes. */
            t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
            /* Retrieve parameters of nodes und curve */
            const double *curve_parameters
              = (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (),
                                                   T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY
                                                     + t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge],
                                                   ltreeid);
            T8_ASSERT (curve_parameters != NULL);

            /* Interpolate linearly between the parameters of the two nodes on the curve */
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], curve_parameters, 1, 1,
                                          &interpolated_curve_param);
            /* Do the same interpolation but with the surface parameters of the same two nodes as above */
            double interpolated_surface_parameters_on_edge[2];
            double edge_parameters_on_face[4];
            t8_geom_get_face_vertices ((t8_eclass_t) t8_eclass_face_types[active_tree_class][i_faces],
                                       surface_parameters, i_face_edge, 2, edge_parameters_on_face);
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], edge_parameters_on_face, 2, 1,
                                          interpolated_surface_parameters_on_edge);
            /* Do the same interpolation but this time between the coordinates of the same two nodes as above */
            double interpolated_edge_coordinates[3];
            double edge_vertices_on_face[6];
            t8_geom_get_face_vertices ((t8_eclass_t) t8_eclass_face_types[active_tree_class][i_faces],
                                       temp_face_vertices, i_face_edge, 3, edge_vertices_on_face);
            t8_geom_linear_interpolation (&ref_coords[edge_direction + offset_3d], edge_vertices_on_face, 3, 1,
                                          interpolated_edge_coordinates);

            /* Retrieve the curve of the edge */
            T8_ASSERT (edges[t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge]]
                       <= cad_shape_edge_map.Size ());
            curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (
                                        edges[t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge]])),
                                      first, last);
            /* Check if curve is valid */
            T8_ASSERT (!curve.IsNull ());
            /* Calculate point on curve with interpolated parameters */
            curve->D0 (interpolated_curve_param, pnt);

            /* Calculate the displacement generated by the presence of the curve */
            if (i_face_edge % 2 == 0) {
              for (int dim = 0; dim <= 2; ++dim) {
                face_displacement_from_edges[dim] += (1 - ref_coords[orthogonal_direction_of_edge_on_face + offset_3d])
                                                     * (pnt.Coord (dim + 1) - interpolated_edge_coordinates[dim]);
              }
            }
            else {
              for (int dim = 0; dim <= 2; ++dim) {
                face_displacement_from_edges[dim] += (ref_coords[orthogonal_direction_of_edge_on_face + offset_3d])
                                                     * (pnt.Coord (dim + 1) - interpolated_edge_coordinates[dim]);
              }
            }
            /* Convert the interpolated parameter of the curve into the corresponding parameters on the surface */
            const int num_face_nodes = t8_eclass_num_vertices[T8_ECLASS_QUAD];
            t8_geometry_cad::t8_geom_edge_parameter_to_face_parameters (
              edges[t8_face_edge_to_tree_edge[T8_ECLASS_HEX][i_faces][i_face_edge]], faces[i_faces], num_face_nodes,
              interpolated_curve_param, surface_parameters, surface_parameters_from_curve);

            /* Calculate the displacement between the interpolated parameters on the surface 
            * and the parameters on the surface converted from the parameter of the curve
            * and scale them with the corresponding ref coord */
            if (i_face_edge % 2 == 0) {
              for (int dim = 0; dim < 2; ++dim) {
                surface_parameter_displacement_from_edges[dim]
                  += (1 - ref_coords[orthogonal_direction_of_edge_on_face + offset_3d])
                     * (surface_parameters_from_curve[dim] - interpolated_surface_parameters_on_edge[dim]);
              }
            }
            else {
              for (int dim = 0; dim < 2; ++dim) {
                surface_parameter_displacement_from_edges[dim]
                  += ref_coords[orthogonal_direction_of_edge_on_face + offset_3d]
                     * (surface_parameters_from_curve[dim] - interpolated_surface_parameters_on_edge[dim]);
              }
            }
          }
        }

        /* Do a bilinear interpolation between the nodes of the face and the parameters of the face */
        interpolation_coeffs[0] = ref_coords[((face_normal_direction + 1) % 3) + offset_3d];
        interpolation_coeffs[1] = ref_coords[((face_normal_direction + 2) % 3) + offset_3d];
        /* The normal vectors of the faces 0, 1, 4, 5 of a hex point in the same direction
        * as the corresponding axis (normal(f0) = (1, 0, 0)). The faces 2 and 3 are oriented
        * in the other direction (normal(f2) = (0, -1, 0)). Therefore we have to switch two 
        * vertices to change the orientation. Here we switch vertex 1 and 2. */
        double temp_surface_parameters[8];
        memcpy (temp_surface_parameters, surface_parameters, sizeof (double) * 8);

        if (face_normal_direction == 1) {
          double temp;
          /* Switch vertex coordinates */
          for (int dim = 0; dim < 3; ++dim) {
            temp = temp_face_vertices[1 * 3 + dim];
            temp_face_vertices[1 * 3 + dim] = temp_face_vertices[2 * 3 + dim];
            temp_face_vertices[2 * 3 + dim] = temp;
          }
          /* Switch vertex parameters */
          for (int dim = 0; dim < 2; ++dim) {
            temp = temp_surface_parameters[1 * 2 + dim];
            temp_surface_parameters[1 * 2 + dim] = temp_surface_parameters[2 * 2 + dim];
            temp_surface_parameters[2 * 2 + dim] = temp;
          }
        }
        /* The actual bilinear interpolations */
        t8_geom_linear_interpolation (interpolation_coeffs, temp_face_vertices, 3, 2, interpolated_coords + offset_3d);
        t8_geom_linear_interpolation (interpolation_coeffs, temp_surface_parameters, 2, 2, interpolated_surface_params);

        for (int dim = 0; dim < 3; ++dim) {
          /* Correct the interpolated coordinates with the displacement generated by the linked edges */
          interpolated_coords[dim + offset_3d] += face_displacement_from_edges[dim];
        }

        for (int dim = 0; dim < 2; ++dim) {
          /* Correct the interpolated parameters with the parameter displacement generated by the edges */
          interpolated_surface_params[dim] += surface_parameter_displacement_from_edges[dim];
        }

        /* Retrieve the surface of the edge */
        T8_ASSERT (faces[i_faces] <= cad_shape_face_map.Size ());
        surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (faces[i_faces])));

        /* Check if surface is valid */
        T8_ASSERT (!surface.IsNull ());

        /* Compute point on surface with interpolated parameters */
        surface->D0 (interpolated_surface_params[0], interpolated_surface_params[1], pnt);

        /* Compute the displacement between surface and interpolated coords, scale them with the appropriate ref_coord 
        * and add them to the out_coords. */
        if (i_faces % 2 == 0) {
          out_coords[offset_3d]
            += (pnt.X () - interpolated_coords[offset_3d]) * (1 - ref_coords[face_normal_direction + offset_3d]);
          out_coords[offset_3d + 1]
            += (pnt.Y () - interpolated_coords[offset_3d + 1]) * (1 - ref_coords[face_normal_direction + offset_3d]);
          out_coords[offset_3d + 2]
            += (pnt.Z () - interpolated_coords[offset_3d + 2]) * (1 - ref_coords[face_normal_direction + offset_3d]);
        }
        else {
          out_coords[offset_3d]
            += (pnt.X () - interpolated_coords[offset_3d]) * ref_coords[face_normal_direction + offset_3d];
          out_coords[offset_3d + 1]
            += (pnt.Y () - interpolated_coords[offset_3d + 1]) * ref_coords[face_normal_direction + offset_3d];
          out_coords[offset_3d + 2]
            += (pnt.Z () - interpolated_coords[offset_3d + 2]) * ref_coords[face_normal_direction + offset_3d];
        }
      }
    }
  }
  T8_FREE (interpolated_coords);
}

void
t8_geometry_cad::t8_geom_evaluate_cad_prism (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                             const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (active_tree_class == T8_ECLASS_PRISM);
  /* The array contains the coordinate [x,y,z] to interpolate for each edge of a prism. 
   * For example: On edge 0 the interpolation coordinate is y. */
  const int t8_interpolation_coefficient_prism_edge[9] = { 1, 0, 0, 1, 0, 0, 2, 2, 2 };
  /* The array contains the coordinates [x,y,z] to interpolate for each face of a prism. 
   * For example: On face 0 the interpolation coordinates are y and z. */
  const int t8_interpolation_coefficients_prism_face[5][2] = { { 1, 2 }, { 0, 2 }, { 0, 2 }, { 0, 1 }, { 0, 1 } };

  /* Compute coordinates in global space from ref_coords in order to shift them afterwards */
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  double interpolated_curve_param, interpolated_surface_params[2], cur_delta[3];
  gp_Pnt pnt;
  double interpolated_coords[3], interpolation_coeffs[3], temp_face_vertices[T8_ECLASS_MAX_CORNERS_2D * 3],
    temp_edge_vertices[2 * 3];
  Handle_Geom_Curve curve;
  Handle_Geom_Surface surface;
  Standard_Real first, last;

  /* Check each edge for a geometry. */
  for (int i_edge = 0; i_edge < T8_DPRISM_EDGES; ++i_edge) {
    /* We have to check for curves as well as surfaces. Linked curves are stored 
     * in the first half of the array, surfaces in the second. 
     * If a curve is connected to this edge we have to also check, 
     * if a surface is connected to at least one of the two adjacent faces. */
    if (edges[i_edge] > 0 || edges[i_edge + T8_DPRISM_EDGES] > 0) {
      /* Check if only a surface or a curve is present. Abort if both is true. */
      T8_ASSERT (!(edges[i_edge] > 0) != !(edges[i_edge + T8_DPRISM_EDGES] > 0));
      /*
       *     z     y
       *     |  _-                _-4
       *     |_-                _- / \
       *     0----x           _-  /   \
       *                    _-   /     \
       *                  E6   E5      E3 
       *                _-     /         \ 
       *              _-      /           \ 
       *            _-      _3-----E4----_-5
       *           1      _-           _-
       *          / \   E8           _-
       *         /   \ -           E7
       *        /   _-\          _-
       *      E2  _-  E0       _-
       *      / _-      \    _-
       *     /_-         \ _-
       *    0-----E1------2
       */

      /* Save the edge vertices temporarily. */
      t8_geom_get_edge_vertices (active_tree_class, active_tree_vertices, i_edge, 3, temp_edge_vertices);

      /* Loop for batch processing of reference points */
      for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
        const int offset_3d = i_coord * 3;
        /* Interpolate between the edge vertices */
        t8_geom_linear_interpolation (&ref_coords[t8_interpolation_coefficient_prism_edge[i_edge] + offset_3d],
                                      temp_edge_vertices, 3, 1, interpolated_coords);

        /* Get the parameters of the curve/surface linked to the current edge */
        const double *parameters = (double *) t8_cmesh_get_attribute (
          cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_edge, ltreeid);
        T8_ASSERT (parameters != NULL);

        /* Curves have only one parameter u, surfaces have two, u and v.
        * Therefore, we have to distinguish if the edge has a curve or surface linked to it. */
        if (edges[i_edge] > 0) {
          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&ref_coords[t8_interpolation_coefficient_prism_edge[i_edge] + offset_3d],
                                        parameters, 1, 1, &interpolated_curve_param);

          T8_ASSERT (edges[i_edge] <= cad_shape_edge_map.Size ());
          curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_edge])), first, last);

          /* Check if curve are valid */
          T8_ASSERT (!curve.IsNull ());

          /* Compute point on curve with interpolated parameters. */
          curve->D0 (interpolated_curve_param, pnt);
        }
        else {
          /* Linear interpolation between parameters */
          t8_geom_linear_interpolation (&ref_coords[t8_interpolation_coefficient_prism_edge[i_edge] + offset_3d],
                                        parameters, 2, 1, interpolated_surface_params);

          T8_ASSERT (edges[i_edge + T8_DPRISM_EDGES] <= cad_shape_face_map.Size ());
          surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (edges[i_edge + T8_DPRISM_EDGES])));

          /* Check if surface is valid */
          T8_ASSERT (!surface.IsNull ());

          /* Compute point on surface with interpolated parameters */
          surface->D0 (interpolated_surface_params[0], interpolated_surface_params[1], pnt);
        }

        /* Compute displacement between vertex interpolation and curve evaluation with interpolated parameters */
        cur_delta[0] = pnt.X () - interpolated_coords[0];
        cur_delta[1] = pnt.Y () - interpolated_coords[1];
        cur_delta[2] = pnt.Z () - interpolated_coords[2];

        /* Save the scaling factor for both neighbouring faces of the current edge.
         * The scaling factors scale the edge displacement orthogonal along the adjacent faces. */
        double scaling_factor_neigh_face_1 = t8_geom_get_scaling_factor_of_edge_on_face_prism (
          i_edge, t8_edge_to_face[active_tree_class][i_edge][0], ref_coords + offset_3d);
        double scaling_factor_neigh_face_2 = t8_geom_get_scaling_factor_of_edge_on_face_prism (
          i_edge, t8_edge_to_face[active_tree_class][i_edge][1], ref_coords + offset_3d);
        double scaling_factor = scaling_factor_neigh_face_1 * scaling_factor_neigh_face_2;

        /* Add edge displacements to out_coords */
        out_coords[offset_3d + 0] += cur_delta[0] * scaling_factor;
        out_coords[offset_3d + 1] += cur_delta[1] * scaling_factor;
        out_coords[offset_3d + 2] += cur_delta[2] * scaling_factor;
      }
    }
  }
  /* Iterate over each face to calculate the displacements generated by each face */
  for (int i_faces = 0; i_faces < T8_DPRISM_FACES; ++i_faces) {
    /* Check if face has a linked surface */
    if (faces[i_faces] > 0) {
      /* Save the face vertices for later usage */
      t8_geom_get_face_vertices (active_tree_class, active_tree_vertices, i_faces, 3, temp_face_vertices);

      /* Retrieve surface parameters of nodes */
      const double *surface_parameters = (double *) t8_cmesh_get_attribute (
        cmesh, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + i_faces, ltreeid);
      T8_ASSERT (surface_parameters != NULL);

      /* Loop for batch processing of reference points */
      for (size_t coord = 0; coord < num_coords; ++coord) {
        const int offset_3d = coord * 3;

        double face_displacement_from_edges[3] = { 0 };

        /* Iterate over each edge of face */
        for (int i_face_edge = 0;
             i_face_edge < t8_eclass_num_vertices[t8_eclass_face_types[active_tree_class][i_faces]]; ++i_face_edge) {

          const int i_tree_edge = t8_face_edge_to_tree_edge[active_tree_class][i_faces][i_face_edge];
          const int interpolation_coeff = t8_interpolation_coefficient_prism_edge[i_tree_edge];

          /* Check if curve is present */
          if (edges[i_tree_edge] > 0) {
            /* Convert global tree id to local tree id, for receiving cmesh attributes. */
            t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
            /* Retrieve parameters of nodes of the curve */
            const double *curve_parameters = (double *) t8_cmesh_get_attribute (
              cmesh, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + i_tree_edge, ltreeid);
            T8_ASSERT (curve_parameters != NULL);

            /* Interpolate linearly between the parameters of the two nodes on the curve */
            t8_geom_linear_interpolation (&ref_coords[interpolation_coeff + offset_3d], curve_parameters, 1, 1,
                                          &interpolated_curve_param);

            /* Do the same interpolation but this time between the coordinates of the same two nodes as above */
            double interpolated_edge_coordinates[3];
            double edge_vertices_on_face[6];
            t8_geom_get_edge_vertices (active_tree_class, active_tree_vertices, i_tree_edge, 3, edge_vertices_on_face);
            t8_geom_linear_interpolation (&ref_coords[interpolation_coeff + offset_3d], edge_vertices_on_face, 3, 1,
                                          interpolated_edge_coordinates);

            /* Retrieve the curve of the edge */
            T8_ASSERT (edges[i_tree_edge] <= cad_shape_edge_map.Size ());
            curve = BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (edges[i_tree_edge])), first, last);
            /* Check if curve is valid */
            T8_ASSERT (!curve.IsNull ());
            /* Calculate point on curve with interpolated parameters */
            curve->D0 (interpolated_curve_param, pnt);

            /* Compute the scaling_factor of the edge displacement on the current face */
            double scaling_factor = t8_geom_get_scaling_factor_of_edge_on_face_prism (
              t8_face_edge_to_tree_edge[active_tree_class][i_faces][i_face_edge], i_faces, ref_coords + offset_3d);

            /* Save the calculated and scaled displacement */
            for (int dim = 0; dim <= 2; ++dim) {
              face_displacement_from_edges[dim]
                += (pnt.Coord (dim + 1) - interpolated_edge_coordinates[dim]) * scaling_factor;
            }
          }
        }

        interpolation_coeffs[0] = ref_coords[t8_interpolation_coefficients_prism_face[i_faces][0] + offset_3d];
        interpolation_coeffs[1] = ref_coords[t8_interpolation_coefficients_prism_face[i_faces][1] + offset_3d];

        /* Do a bilinear interpolation between the nodes of the face and the parameters of the face */
        if (i_faces <= 2) {
          t8_geom_linear_interpolation (interpolation_coeffs, temp_face_vertices, 3, 2, interpolated_coords);
          t8_geom_linear_interpolation (interpolation_coeffs, surface_parameters, 2, 2, interpolated_surface_params);
        }
        else {
          t8_geom_triangular_interpolation (interpolation_coeffs, temp_face_vertices, 3, 2, interpolated_coords);
          t8_geom_triangular_interpolation (interpolation_coeffs, surface_parameters, 2, 2,
                                            interpolated_surface_params);
        }

        for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
          /* Correct the interpolated coordinates with the displacement generated by the linked edges */
          interpolated_coords[dim] += face_displacement_from_edges[dim];
        }

        /* Retrieve the surface of the edge */
        T8_ASSERT (faces[i_faces] <= cad_shape_face_map.Size ());
        surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (faces[i_faces])));

        /* Check if surface is valid */
        T8_ASSERT (!surface.IsNull ());

        /* Compute point on surface with interpolated parameters */
        surface->D0 (interpolated_surface_params[0], interpolated_surface_params[1], pnt);

        /* Compute the displacement between surface and interpolated coords, scale them with the appropriate scaling_factor 
         * and add them to the out_coords. */
        double scaling_factor = t8_geom_get_scaling_factor_face_through_volume_prism (i_faces, ref_coords + offset_3d);

        out_coords[offset_3d] += (pnt.X () - interpolated_coords[0]) * scaling_factor;
        out_coords[offset_3d + 1] += (pnt.Y () - interpolated_coords[1]) * scaling_factor;
        out_coords[offset_3d + 2] += (pnt.Z () - interpolated_coords[2]) * scaling_factor;
      }
    }
  }
}

int
t8_geometry_cad::t8_geom_is_line (const int curve_index) const
{
  const Handle_Geom_Curve curve = t8_geom_get_cad_curve (curve_index);
  const GeomAdaptor_Curve curve_adaptor (curve);
  return curve_adaptor.GetType () == GeomAbs_Line;
}

int
t8_geometry_cad::t8_geom_is_plane (const int surface_index) const
{
  const Handle_Geom_Surface surface = t8_geom_get_cad_surface (surface_index);
  const GeomAdaptor_Surface surface_adaptor (surface);
  return surface_adaptor.GetType () == GeomAbs_Plane;
}

const gp_Pnt
t8_geometry_cad::t8_geom_get_cad_point (const int index) const
{
  T8_ASSERT (index <= cad_shape_vertex_map.Size ());
  return BRep_Tool::Pnt (TopoDS::Vertex (cad_shape_vertex_map.FindKey (index)));
}

const Handle_Geom_Curve
t8_geometry_cad::t8_geom_get_cad_curve (const int index) const
{
  T8_ASSERT (index <= cad_shape_edge_map.Size ());
  Standard_Real first, last;
  return BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (index)), first, last);
}

const Handle_Geom_Surface
t8_geometry_cad::t8_geom_get_cad_surface (const int index) const
{
  T8_ASSERT (index <= cad_shape_face_map.Size ());
  return BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (index)));
}

const TopTools_IndexedMapOfShape
t8_geometry_cad::t8_geom_get_cad_shape_vertex_map () const
{
  return cad_shape_vertex_map;
}

const TopTools_IndexedMapOfShape
t8_geometry_cad::t8_geom_get_cad_shape_edge_map () const
{
  return cad_shape_edge_map;
}

const TopTools_IndexedMapOfShape
t8_geometry_cad::t8_geom_get_cad_shape_face_map () const
{
  return cad_shape_face_map;
}

int
t8_geometry_cad::t8_geom_get_common_edge (const int vertex1_index, const int vertex2_index) const
{
  const TopTools_ListOfShape collection1 = cad_shape_vertex2edge_map.FindFromIndex (vertex1_index);
  const TopTools_ListOfShape collection2 = cad_shape_vertex2edge_map.FindFromIndex (vertex2_index);

  for (auto edge1 = collection1.begin (); edge1 != collection1.end (); ++edge1) {
    for (auto edge2 = collection2.begin (); edge2 != collection2.end (); ++edge2) {
      if (edge1->IsEqual (*edge2)) {
        return cad_shape_edge2face_map.FindIndex (*edge1);
      }
    }
  }
  return 0;
}

int
t8_geometry_cad::t8_geom_get_common_face (const int edge1_index, const int edge2_index) const
{
  const TopTools_ListOfShape collection1 = cad_shape_edge2face_map.FindFromIndex (edge1_index);
  const TopTools_ListOfShape collection2 = cad_shape_edge2face_map.FindFromIndex (edge2_index);

  for (auto face1 = collection1.begin (); face1 != collection1.end (); ++face1) {
    for (auto face2 = collection2.begin (); face2 != collection2.end (); ++face2) {
      if (face1->IsEqual (*face2)) {
        return cad_shape_face_map.FindIndex (*face1);
      }
    }
  }
  return 0;
}

int
t8_geometry_cad::t8_geom_is_vertex_on_edge (const int vertex_index, const int edge_index) const
{
  const TopTools_ListOfShape collection = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
  return collection.Contains (cad_shape_edge_map.FindKey (edge_index));
}

int
t8_geometry_cad::t8_geom_is_edge_on_face (const int edge_index, const int face_index) const
{
  const TopTools_ListOfShape collection = cad_shape_edge2face_map.FindFromIndex (edge_index);
  return collection.Contains (cad_shape_face_map.FindKey (face_index));
}

int
t8_geometry_cad::t8_geom_is_vertex_on_face (const int vertex_index, const int face_index) const
{
  const TopTools_ListOfShape edge_collection = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
  for (auto edge = edge_collection.begin (); edge != edge_collection.end (); ++edge) {
    const TopTools_ListOfShape face_collection = cad_shape_edge2face_map.FindFromKey (*edge);
    if (face_collection.Contains (cad_shape_face_map.FindKey (face_index))) {
      return 1;
    }
  }
  return 0;
}

void
t8_geometry_cad::t8_geom_get_parameter_of_vertex_on_edge (const int vertex_index, const int edge_index,
                                                          double *edge_param) const
{
  T8_ASSERT (t8_geometry_cad::t8_geom_is_vertex_on_edge (vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));
  *edge_param = BRep_Tool::Parameter (vertex, edge);
}

void
t8_geometry_cad::t8_geom_get_parameters_of_vertex_on_face (const int vertex_index, const int face_index,
                                                           double *face_params) const
{
  T8_ASSERT (t8_geometry_cad::t8_geom_is_vertex_on_face (vertex_index, face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));
  uv = BRep_Tool::Parameters (vertex, face);
  face_params[0] = uv.X ();
  face_params[1] = uv.Y ();
}

void
t8_geometry_cad::t8_geom_edge_parameter_to_face_parameters (const int edge_index, const int face_index,
                                                            const int num_face_nodes, const double edge_param,
                                                            const double *surface_params, double *face_params) const
{
  T8_ASSERT (t8_geometry_cad::t8_geom_is_edge_on_face (edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));
  Handle_Geom2d_Curve curve_on_surface = BRep_Tool::CurveOnSurface (edge, face, first, last);
  Handle_Geom_Surface surface = BRep_Tool::Surface (face);
  curve_on_surface->D0 (edge_param, uv);
  face_params[0] = uv.X ();
  face_params[1] = uv.Y ();

  /* Check for right conversion of edge to surface parameter and correct if needed */
  /* Checking u parameter */
  if (surface_params != NULL) {
    double parametric_bounds[4];
    surface->Bounds (parametric_bounds[0], parametric_bounds[1], parametric_bounds[2], parametric_bounds[3]);
    if (surface->IsUClosed ()) {
      for (int i_face_node = 0; i_face_node < num_face_nodes; ++i_face_node) {
        if (surface_params[i_face_node * 2] == parametric_bounds[0]) {
          if (face_params[0] == parametric_bounds[1]) {
            face_params[0] = parametric_bounds[0];
          }
        }
        else if (surface_params[i_face_node * 2] == parametric_bounds[1]) {
          if (face_params[0] == parametric_bounds[0]) {
            face_params[0] = parametric_bounds[1];
          }
        }
      }
    }
    /* Checking v parameter */
    if (surface->IsVClosed ()) {
      for (int i_face_node = 0; i_face_node < num_face_nodes; ++i_face_node) {
        if (surface_params[i_face_node * 2 + 1] == parametric_bounds[0]) {
          if (face_params[1] == parametric_bounds[1]) {
            face_params[1] = parametric_bounds[0];
          }
        }
        else if (surface_params[i_face_node * 2 + 1] == parametric_bounds[1]) {
          if (face_params[1] == parametric_bounds[0]) {
            face_params[1] = parametric_bounds[1];
          }
        }
      }
    }
  }
}

void
t8_geometry_cad::t8_geom_get_face_parametric_bounds (const int surface_index, double *bounds) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (surface_index);
  cad_surface->Bounds (bounds[0], bounds[1], bounds[2], bounds[3]);
}

void
t8_geometry_cad::t8_geom_get_edge_parametric_bounds (const int edge_index, double *bounds) const
{
  const Handle_Geom_Curve cad_edge = t8_geom_get_cad_curve (edge_index);
  bounds[0] = cad_edge->FirstParameter ();
  bounds[1] = cad_edge->LastParameter ();
}

int
t8_geometry_cad::t8_geom_is_edge_closed (int edge_index) const
{
  const Handle_Geom_Curve cad_edge = t8_geom_get_cad_curve (edge_index);
  return cad_edge->IsClosed ();
}

int
t8_geometry_cad::t8_geom_is_surface_closed (int geometry_index, int parameter) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (geometry_index);
  switch (parameter) {
  case 0:
    return cad_surface->IsUClosed ();
    break;
  case 1:
    return cad_surface->IsVClosed ();
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
}

/* This part should be callable from C */
T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_cad.h.
 * Create a new geometry. */
t8_geometry_cad_c *
t8_geometry_cad_new (const char *fileprefix, const char *name_in)
{
  t8_geometry_cad *geom = new t8_geometry_cad (fileprefix, name_in);
  return (t8_geometry_cad_c *) geom;
}

void
t8_geometry_cad_destroy (t8_geometry_cad_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_CAD);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();
