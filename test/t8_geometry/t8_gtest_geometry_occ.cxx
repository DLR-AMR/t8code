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

#include <gtest/gtest.h>
#include <t8_cmesh/t8_cmesh.c>
#include <t8_forest/t8_forest_general.h>
#include <t8_vtk.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>

#if T8_WITH_OCC
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#endif

/* In this file we collect tests for t8code's OpenCASCADE geometry module.
 * These tests are
 *  - linked_edges: Checks if the geometry mapping works. 
 *                  We define an OpenCASCADE curve and link it to an edge of a hexahedron. 
 *                  After that we probe whether the correct coordinates inside the hexahedron are returned.
 *                  We repeat this check for all 12 edges.
 *  - linked_faces: Checks if the geometry mapping works. 
 *                  We define an OpenCASCADE surface and link it to a face of a hexahedron. 
 *                  After that we probe whether the correct coordinates inside the hexahedron are returned.
 *                  We repeat this check for all 6 faces.
 *  - jacobian:     Checks the resulting jacobian of an identity.
 */

#if T8_WITH_OCC
/** Euler rotation around intrinsic zxz. 
 * \param [in] pos_vec                Position vector of three dimensional points to rotate.
 * \param [in] rot_vec                Three dimensional rotation vector around z, x and z in rad.
 * \param [out] res_vec               Vector for resulting points.
 * \param [in] rot_origin             Origin for Euler rotation.
 * \param [in] n_points               Number of points in pos_vec.
 */
static void
t8_euler_rotation (double *pos_vec,
                   double *rot_vec,
                   double *res_vec, double *rot_origin, int n_points)
{
  for (int point_iterator = 0; point_iterator < n_points; ++point_iterator) {
    res_vec[point_iterator * 3 + 0] =
      (cos (rot_vec[0]) * cos (rot_vec[2]) -
       sin (rot_vec[0]) * cos (rot_vec[1]) * sin (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 0] - rot_origin[0]);
    res_vec[point_iterator * 3 + 0] +=
      (-cos (rot_vec[0]) * sin (rot_vec[2]) -
       sin (rot_vec[0]) * cos (rot_vec[1]) * cos (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 1] - rot_origin[1]);
    res_vec[point_iterator * 3 + 0] +=
      (sin (rot_vec[0]) * sin (rot_vec[1])) *
      (pos_vec[point_iterator * 3 + 2] - rot_origin[2]);
    res_vec[point_iterator * 3 + 0] += rot_origin[0];

    res_vec[point_iterator * 3 + 1] =
      (sin (rot_vec[0]) * cos (rot_vec[2]) +
       cos (rot_vec[0]) * cos (rot_vec[1]) * sin (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 0] - rot_origin[0]);
    res_vec[point_iterator * 3 + 1] +=
      (-sin (rot_vec[0]) * sin (rot_vec[2]) +
       cos (rot_vec[0]) * cos (rot_vec[1]) * cos (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 1] - rot_origin[1]);
    res_vec[point_iterator * 3 + 1] +=
      (-cos (rot_vec[0]) * sin (rot_vec[1])) *
      (pos_vec[point_iterator * 3 + 2] - rot_origin[2]);
    res_vec[point_iterator * 3 + 1] += rot_origin[1];

    res_vec[point_iterator * 3 + 2] =
      (sin (rot_vec[1]) * sin (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 0] - rot_origin[0]);
    res_vec[point_iterator * 3 + 2] +=
      (sin (rot_vec[1]) * cos (rot_vec[2])) *
      (pos_vec[point_iterator * 3 + 1] - rot_origin[1]);
    res_vec[point_iterator * 3 + 2] +=
      (cos (rot_vec[1])) * (pos_vec[point_iterator * 3 + 2] - rot_origin[2]);
    res_vec[point_iterator * 3 + 2] += rot_origin[2];
  }
}

/** Constructs an occ surface for testing purposes. Surface is build between vertex 0, 1, 4 and 5 of a unit hexahedron.
 * Saves the surface in the geometry shape.
 * \return                            The geometry.
 */
t8_geometry_occ    *
t8_create_occ_surface_geometry ()
{
  Handle_Geom_Surface surface;
  TopoDS_Shape        shape;
  TColgp_Array2OfPnt  point_array (1, 3, 1, 3);

  point_array (1, 1) = gp_Pnt (0, 0, 0);
  point_array (2, 1) = gp_Pnt (-0.2, 0.2, 0.5);
  point_array (3, 1) = gp_Pnt (0, 0, 1);

  point_array (1, 2) = gp_Pnt (0.5, 0.2, -0.2);
  point_array (2, 2) = gp_Pnt (0.5, 0, 0.5);
  point_array (3, 2) = gp_Pnt (0.5, -0.2, 1.2);

  point_array (1, 3) = gp_Pnt (1, 0, 0);
  point_array (2, 3) = gp_Pnt (1.2, -0.2, 0.5);
  point_array (3, 3) = gp_Pnt (1, 0, 1);

  surface = GeomAPI_PointsToBSplineSurface (point_array).Surface ();
  shape = BRepBuilderAPI_MakeFace (surface, 1e-6).Face ();
  t8_geometry_occ    *geometry = new t8_geometry_occ (3, shape, "occ dim=3");
  return geometry;
}

/** Constructs an occ curve for testing purpsoes. Curve is build between vertex 0, 1, 4 and 5 of a unit hexahedron.
 * Saves the curve in the geometry shape.
 * \return                            The occ geometry.
 */
t8_geometry_occ    *
t8_create_occ_curve_geometry ()
{
  Handle_Geom_Curve   curve;
  TopoDS_Shape        shape;
  TColgp_Array1OfPnt  point_array (1, 5);

  point_array (1) = gp_Pnt (0.00, 0.00, 0.00);
  point_array (2) = gp_Pnt (0.24, 0.20, 0.20);
  point_array (3) = gp_Pnt (0.50, 0.00, 0.40);
  point_array (4) = gp_Pnt (0.75, -0.20, 0.20);
  point_array (5) = gp_Pnt (1.00, 0.00, 0.00);

  curve = GeomAPI_PointsToBSpline (point_array).Curve ();
  shape = BRepBuilderAPI_MakeEdge (curve).Edge ();
  t8_geometry_occ    *geometry = new t8_geometry_occ (3, shape, "occ dim=3");
  return geometry;
}
#endif /* T8_WITH_OCC */

/** Constructs a cmesh with an occ geometry linked hypercube.
 * \param [in] rot_vec                The rotation vector to rotate the cube before linking a geometry to it.
 * \param [in] face                   The index of the face to link a surface to. -1 for no face.
 * \param [in] edge                   The index of the edge to link a curve to. -1 for no edge.
 * \param [in] parameters             Parameters of the curve/surface.
 * \return                            A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t
t8_create_occ_hypercube (double *rot_vec,
                         int face, int edge, double *parameters)
{
#if T8_WITH_OCC
  if (edge >= 0 && face >= 0) {
    SC_ABORTF ("Please specify only an edge or a face.");
  }

  t8_cmesh_t          cmesh;
  t8_geometry_occ    *geometry = NULL;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);

  double              rotated_vertices[24], vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };
  double              rotation_origin[3] = { 0.5, 0.5, 0.5 };

  t8_euler_rotation (vertices, rot_vec, rotated_vertices, rotation_origin, 8);
  t8_cmesh_set_tree_vertices (cmesh, 0, rotated_vertices, 24);

  int                 faces[6] = { 0 };
  int                 edges[24] = { 0 };
  T8_ASSERT (face < 0 || edge < 0);
  if (face >= 0) {
    faces[face] = 1;
    geometry = t8_create_occ_surface_geometry ();
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (),
                            T8_CMESH_OCC_FACE_PARAMETERS_ATTRIBUTE_KEY + face,
                            parameters, 8 * sizeof (double), 0);
  }
  else if (edge >= 0) {
    edges[edge] = 1;
    geometry = t8_create_occ_curve_geometry ();
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (),
                            T8_CMESH_OCC_EDGE_PARAMETERS_ATTRIBUTE_KEY + edge,
                            parameters, 2 * sizeof (double), 0);
  }
  else {
    /* Even if we do not want to link any geometry to the edges or faces, 
     * we have to create a geometry. Hence an occ geometry can only be created
     * with an actual shape, we just create a geometry with a curve and do not
     * link the curve to any edge. */
    geometry = t8_create_occ_curve_geometry ();
  }
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (),
                          T8_CMESH_OCC_FACE_ATTRIBUTE_KEY, faces,
                          6 * sizeof (int), 0);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (),
                          T8_CMESH_OCC_EDGE_ATTRIBUTE_KEY, edges,
                          24 * sizeof (int), 0);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  return cmesh;

#else /* !T8_WITH_OCC */
  SC_ABORTF ("OCC not linked");
#endif /* T8_WITH_OCC */
}

/** Tests the occ geometry functions for hexahedra.
 * \param [in] rot_vec                The rotation vector to rotate the hypercube.
 * \param [in] face                   The face to test. -1 for no face.
 * \param [in] edge                   The edge to test. -1 for no edge.
 * \param [in] parameters             The parameters of the curve/surface.
 * \param [in] test_ref_coords        List of 8 coordinates to test.
 * \param [in] test_return_coords     List of 8 expected output coordinates.
 * \param [in] comm                   The mpi communicator to use.
 * \return                            Returns 1 if passed, 0 if failed.
 */
void
t8_test_geometry_occ (double *rot_vec,
                      int face,
                      int edge,
                      double *parameters,
                      double *test_ref_coords, double *test_return_coords)
{
#if T8_WITH_OCC
  double              out_coords[3];
  double              rotated_test_ref_coords[24];
  double              rotation_origin[3] = { 0.5, 0.5, 0.5 };
  double              inversed_rot_vec[3];
  double              tol =
    T8_PRECISION_EPS > 1e-10 ? T8_PRECISION_EPS : 1e-10;
  t8_cmesh_t          cmesh =
    t8_create_occ_hypercube (rot_vec, face, edge, parameters);

  for (int i_coord = 0; i_coord < 3; ++i_coord) {
    inversed_rot_vec[2 - i_coord] = -rot_vec[i_coord];
  }
  t8_euler_rotation (test_ref_coords, inversed_rot_vec,
                     rotated_test_ref_coords, rotation_origin, 8);
  for (int i_coord = 0; i_coord < 8; ++i_coord) {
    t8_geometry_evaluate (cmesh, 0, rotated_test_ref_coords + i_coord * 3,
                          out_coords);
    EXPECT_NEAR (out_coords[0], test_return_coords[0 + i_coord * 3], tol);
    EXPECT_NEAR (out_coords[1], test_return_coords[1 + i_coord * 3], tol);
    EXPECT_NEAR (out_coords[2], test_return_coords[2 + i_coord * 3], tol);
  }
  t8_cmesh_destroy (&cmesh);

#else /* !T8_WITH_OCC */
  SC_ABORTF ("OCC not linked");
#endif /* T8_WITH_OCC */
}

#if T8_WITH_OCC
TEST (t8_gtest_geometry_occ, linked_faces)
{
  double              test_ref_coords[24] = {
    0.1, 0.1, 0.1,
    0.8, 0.1, 0.1,
    0.15, 0.9, 0.1,
    0.9, 0.9, 0.3,
    0.3, 0.1, 0.7,
    0.9, 0.25, 0.95,
    0.1, 0.9, 0.9,
    0.95, 0.85, 0.8
  };
  double              surface_test_return_coords[24] = {
    0.0396282769, 0.1897542602, 0.0396282769,
    0.8553975402, 0.1510451803, -0.0012778561,
    0.1434278361, 0.9117760771, 0.0909403721,
    0.9149739120, 0.8893780561, 0.2953610950,
    0.2190065733, 0.1000000000, 0.7809934267,
    0.9318450385, 0.1898146343, 0.9989190836,
    0.0932920308, 0.9000000000, 0.9067079692,
    0.9673042609, 0.8312979801, 0.8063743210
  };
  double              surface_rot_vecs[18] = {
    M_PI / 2, 0, 0,             // Face 0
    M_PI * 3 / 2, 0, 0,         // Face 1
    0, 0, 0,                    // Face 2
    M_PI, 0, 0,                 // Face 3
    0, M_PI * 3 / 2, 0,         // Face 4
    0, M_PI / 2, 0              // Face 5
  };
  double              surface_parameters[48] = {
    0, 1, 0, 0, 1, 1, 1, 0,     // Face 0
    0, 0, 0, 1, 1, 0, 1, 1,     // Face 1
    0, 0, 0, 1, 1, 0, 1, 1,     // Face 2
    0, 1, 0, 0, 1, 1, 1, 0,     // Face 3
    1, 0, 1, 1, 0, 0, 0, 1,     // Face 4
    0, 0, 0, 1, 1, 0, 1, 1      // Face 5
  };
  for (int i_faces = 0; i_faces < 6; ++i_faces) {
    t8_test_geometry_occ (surface_rot_vecs + i_faces * 3,
                          i_faces,
                          -1,
                          surface_parameters + i_faces * 8,
                          test_ref_coords, surface_test_return_coords);
  }
}

TEST (t8_gtest_geometry_occ, linked_edges)
{
  double              test_ref_coords[24] = {
    0.1, 0.1, 0.1,
    0.8, 0.1, 0.1,
    0.15, 0.9, 0.1,
    0.9, 0.9, 0.3,
    0.3, 0.1, 0.7,
    0.9, 0.25, 0.95,
    0.1, 0.9, 0.9,
    0.95, 0.85, 0.8
  };
  double              curve_test_return_coords[24] = {
    0.0955204602, 0.2235162028, 0.1217553783,
    0.7995278713, -0.0659838746, 0.2083328730,
    0.1494299582, 0.9170222805, 0.1069555502,
    0.8999105642, 0.8892289094, 0.3015732294,
    0.2987855815, 0.1481519479, 0.7726155646,
    0.8999520880, 0.2442297729, 0.9508428015,
    0.0999446970, 0.9015248914, 0.9002685849,
    0.9499697383, 0.8472575225, 0.7998496263
  };

  double              curve_rot_vecs[36] = {
    0, 0, 0,                    // Edge 0
    0, -M_PI / 2, 0,            // Edge 1
    0, M_PI / 2, 0,             // Edge 2
    0, M_PI, 0,                 // Edge 3
    M_PI / 2, 0, 0,             // Edge 4
    0, M_PI, -M_PI / 2,         // Edge 5
    -M_PI / 2, 0, 0,            // Edge 6
    0, M_PI, M_PI / 2,          // Edge 7
    -M_PI / 2, M_PI / 2, M_PI / 2,      // Edge 8
    M_PI / 2, M_PI / 2, -M_PI / 2,      // Edge 9
    M_PI / 2, -M_PI / 2, 0,     // Edge 10
    -M_PI / 2, -M_PI / 2, 0,    // Edge 11
  };
  double              curve_parameters[24] = {
    0, 1,                       // Edge 0
    0, 1,                       // Edge 1
    0, 1,                       // Edge 2
    0, 1,                       // Edge 3
    1, 0,                       // Edge 4
    0, 1,                       // Edge 5
    0, 1,                       // Edge 6
    1, 0,                       // Edge 7
    1, 0,                       // Edge 8
    0, 1,                       // Edge 9
    1, 0,                       // Edge 10
    0, 1                        // Edge 11
  };
  for (int i_edges = 0; i_edges < 12; ++i_edges) {
    t8_test_geometry_occ (curve_rot_vecs + i_edges * 3,
                          -1,
                          i_edges,
                          curve_parameters + i_edges * 2,
                          test_ref_coords, curve_test_return_coords);
  }
}
#endif /* T8_WITH_OCC */

#if T8_WITH_OCC
TEST (t8_gtest_geometry_occ, jacobian)
{
  t8_cmesh_t          cmesh;
  double              jacobian[9], rot_vec[3] = { 0, 0, 0 }, ref_coords[3] =
    { 0.5, 0.5, 0.5 };
  double              jacobian_expect[9] = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  };
  cmesh = t8_create_occ_hypercube (rot_vec, -1, -1, NULL);
  t8_geometry_jacobian (cmesh, 0, ref_coords, jacobian);
  for (int i = 0; i < 9; ++i) {
    EXPECT_FLOAT_EQ (jacobian[i], jacobian_expect[i]);
  }
  t8_cmesh_destroy (&cmesh);
}
#endif /* T8_WITH_OCC */
