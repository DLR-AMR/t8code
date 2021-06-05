/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_cmesh_vtk.h>

#if T8_WITH_OCC
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <gp_Pnt.hxx>
#include <NCollection_Array2.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <Geom_BSplineSurface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <StlAPI.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <BRepSweep_Prism.hxx>

#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Circ.hxx>
#include <gp_Vec.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Geom_Surface.hxx>
#endif

typedef enum
{
  T8_GEOM_ZERO = 0,
  T8_GEOM_SINCOS = T8_GEOM_ZERO,
  T8_GEOM_CYLINDER,
  T8_GEOM_MOEBIUS,
  T8_GEOM_CIRCLE,
  T8_GEOM_3D,
  T8_GEOM_MOVING,
  T8_GEOM_OCC_SURFACE_CUBES,
  T8_GEOM_OCC_SURFACE_CYLINDER,
  T8_GEOM_COUNT
} t8_analytic_geom_type;

static void
t8_analytic_sincos (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                    const double *ref_coords, double out_coords[3],
                    const void *tree_data, const void *user_data)
{
  double              x = ref_coords[0];
  if (gtreeid == 1) {
    /* Translate ref coordinates by +1 in x direction. */
    x += 1;
  }
  out_coords[0] = x;
  out_coords[1] = ref_coords[1];
  out_coords[2] =
    0.2 * sin (ref_coords[0] * 2 * M_PI) * cos (ref_coords[1] * 2 * M_PI);
}

static void
t8_analytic_cylinder (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                      const double *ref_coords, double out_coords[3],
                      const void *tree_data, const void *user_data)
{
  out_coords[0] = cos (ref_coords[0] * 2 * M_PI);
  out_coords[1] = ref_coords[1];
  out_coords[2] = sin (ref_coords[0] * 2 * M_PI);
}

static void
t8_analytic_moebius (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                     const double *ref_coords, double out_coords[3],
                     const void *tree_vertices, const void *user_data)
{
  const double       *tree_v = (const double *) tree_vertices;
  double              t;
  double              phi;

  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_eclass_t         tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store
   * in out_coords. */
  t8_geom_compute_linear_geometry (tree_class, tree_v, ref_coords,
                                   out_coords);

  /* At first, we map x from [0,1] to [-.5,.5]
   * and y to [0, 2*PI] */
  t = out_coords[0] - .5;
  phi = out_coords[1] * 2 * M_PI;

  /* We now apply the parametrization for the moebius strip. */
  out_coords[0] = (1 - t * sin (phi / 2)) * cos (phi);
  out_coords[1] = (1 - t * sin (phi / 2)) * sin (phi);
  out_coords[2] = t * cos (phi / 2);
}

static void
t8_analytic_circle (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                    const double *ref_coords, double out_coords[3],
                    const void *tree_vertices, const void *user_data)
{
  const double       *tree_v = (const double *) tree_vertices;
  double              x;
  double              y;

  t8_locidx_t         ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  t8_eclass_t         tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store
   * in out_coords. */
  t8_geom_compute_linear_geometry (tree_class, tree_v, ref_coords,
                                   out_coords);

  /* We now remap the coords to match the square [-1,1]^2 */
  x = out_coords[0] * 2 - 1;
  y = out_coords[1] * 2 - 1;

  /* An now we apply the formula that projects the square to the circle. */
  out_coords[0] = x * sqrt (1 - y * y / 2);
  out_coords[1] = y * sqrt (1 - x * x / 2);
  out_coords[2] = 0;
}

/* This geometry rotates with time around the origin.
 * The rotation direction is reversed after 2 seconds.
 * Additionally, the z coordinate is modifyied according to the
 * sincos function and multiplied with the current time.
 * To use this, a double variable time has to be set as the geometries
 * user data.
 */
static void
t8_analytic_moving (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                    const double *ref_coords, double out_coords[3],
                    const void *tree_data, const void *user_data)
{
  double              x = ref_coords[0] - .5;
  double              y = ref_coords[1] - .5;
  T8_ASSERT (tree_data == NULL);
  T8_ASSERT (user_data != NULL);
  const double        time = *(double *) user_data;
  double              radius_sqr = x * x + y * y;
  double              phi = radius_sqr * (time > 2 ? 4 - time : time);

  /* Change gridlines by applying a 4th order polynomial mapping
   * [0,1]^2 -> [0,1]^2.
   * And then map this to [-0.5,-0.5]^2 */
  //x = x*x*x*x - 3.5*x*x*x + 3.5 * x;
  int                 sign = x < 0 ? 1 : -1;
  double              rho = 0.5 - time / 10;
  x = sign * (1 - exp (-fabs (-x) / rho)) / (2 * (1 - exp (-0.5 / rho)));
  sign = y < 0 ? 1 : -1;
  y = sign * (1 - exp (-fabs (-y) / rho)) / (2 * (1 - exp (-0.5 / rho)));
  //y = y*y*y - 3.5*y*y + 2.5 * y - 0.5;

  /* Rotate the x-y axis and add sincon in z axis. */
  out_coords[0] = x * (cos (phi)) - y * sin (phi);
  out_coords[1] = y * (cos (phi)) + x * sin (phi);
  out_coords[2] = 0;
  //sin(2*M_PI * time) * 0.2 * sin (out_coords[0] * 2 * M_PI) * cos (out_coords[1] * 2 * M_PI);
}

static void
t8_analytic_3D_cube (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                     const double *ref_coords, double out_coords[3],
                     const void *tree_data, const void *user_data)
{
  out_coords[0] = ref_coords[0];
  out_coords[1] = ref_coords[1];
  out_coords[2] = ref_coords[2] * (0.8 +
                                   0.2 * sin (ref_coords[0] * 2 * M_PI) *
                                   cos (ref_coords[1] * 2 * M_PI));
}

/* This adapt callback function will refine all elements at the
 * domain boundary up to a given maximum refinement level. */
static int
t8_geom_adapt_boundary (t8_forest_t forest, t8_forest_t forest_from,
                        t8_locidx_t ltree_id, t8_locidx_t lelement_id,
                        t8_eclass_scheme_c * ts, int num_elements,
                        t8_element_t * elements[])
{
  t8_cmesh_t          cmesh = t8_forest_get_cmesh (forest_from);
  /* Get the number of faces of the element. */
  int                 num_faces = ts->t8_element_num_faces (elements[0]);
  int                 iface;
  /* Get the maximum level from the forest's user data 
   * (must be set before using the callback). */
  int                 maxlevel = *(int *) t8_forest_get_user_data (forest);

  /* We do not refine more then the given maximum level. */
  if (ts->t8_element_level (elements[0]) >= maxlevel) {
    return 0;
  }

  /* Check for each face of the element whether it lies on the 
   * domain boundary. If so, the element is refined. */
  for (iface = 0; iface < num_faces; ++iface) {
    if (ts->t8_element_is_root_boundary (elements[0], iface)) {
      /* This element's face is at its tree boundary. Check whether
         the tree's face is at the domain boundary. */
      int                 tree_face =
        ts->t8_element_tree_face (elements[0], iface);
      t8_locidx_t         lctreeid =
        t8_forest_ltreeid_to_cmesh_ltreeid (forest_from, ltree_id);
      if (t8_cmesh_tree_face_is_boundary (cmesh, lctreeid, tree_face)) {
        /* The tree's face is at the domain boundary, we refine the element. */
        return 1;
      }
    }
  }
  /* All other elements remain unchanged. */
  return 0;
}

static void
t8_analytic_geom (int level, t8_analytic_geom_type geom_type)
{
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  char                vtuname[BUFSIZ];
  t8_geometry_c      *geometry;
  int                 uniform_level;
  double              time = 0; /* used for moving geometry */
  int                 sreturn;

  t8_cmesh_init (&cmesh);
  /* Depending on the geometry type, add the tree, set the geometry
   * and set the output file name. */
  switch (geom_type) {
  case T8_GEOM_SINCOS:
    t8_global_productionf
      ("Creating uniform level %i forest with a sinus/cosinus geometry.\n",
       level);
    /* Sin/cos geometry. Has two quad trees. */
    geometry =
      new t8_geometry_analytic (2, "analytic sinus cosinus dim=2",
                                t8_analytic_sincos, NULL, NULL, NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
    snprintf (vtuname, BUFSIZ, "forest_analytic_sincos_lvl_%i", level);
    break;
  case T8_GEOM_CYLINDER:
    t8_global_productionf
      ("Creating uniform level %i forest with a cylinder geometry.\n", level);
    /* Cylinder geometry. Has one quad tree that is periodic in x direction. */
    geometry =
      new t8_geometry_analytic (2, "analytic geometry cylinder",
                                t8_analytic_cylinder, NULL, NULL, NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
    snprintf (vtuname, BUFSIZ, "forest_analytic_cylinder_lvl_%i", level);
    break;
  case T8_GEOM_MOEBIUS:
    t8_global_productionf
      ("Creating uniform level %i forest with a moebius geometry.\n", level);
    {
      /* Moebius geometry on hybrid unit square. */
      t8_cmesh_t          hybrid_square =
        t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
      t8_cmesh_set_derive (cmesh, hybrid_square);
      geometry =
        new t8_geometry_analytic (2, "analytic moebius", t8_analytic_moebius,
                                  NULL, t8_geom_load_tree_data_vertices,
                                  NULL);
      snprintf (vtuname, BUFSIZ, "forest_analytic_moebius_lvl_%i", level);
    }
    break;
  case T8_GEOM_CIRCLE:
    t8_global_productionf ("Creating forest with a circle geometry.\n");
    t8_global_productionf
      ("This forest will get refined at the boundary to level %i.\n", level);
    {
      /* Circle geometry on triangulated unit square. */
      t8_cmesh_t          tri_square =
        t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0,
                                0);
      t8_cmesh_set_derive (cmesh, tri_square);
      geometry =
        new t8_geometry_analytic (2, "analytic circle", t8_analytic_circle,
                                  NULL, t8_geom_load_tree_data_vertices,
                                  NULL);
      snprintf (vtuname, BUFSIZ, "forest_analytic_circle_lvl_%i", level);
    }
    break;
  case T8_GEOM_3D:
    t8_global_productionf
      ("Creating uniform level %i forest with a 3D function graph geometry.\n",
       level);
    /* Cube geometry with sincos on top. Has one hexahedron tree. */
    geometry =
      new t8_geometry_analytic (3, "cube geom", t8_analytic_3D_cube, NULL,
                                NULL, NULL);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    snprintf (vtuname, BUFSIZ, "forest_analytic_3D_lvl_%i", level);
    break;
  case T8_GEOM_MOVING:
    t8_global_productionf
      ("Creating uniform level %i forest with a moving geometry.\n", level);
    /* Quad geometry that rotates with time. */
    geometry =
      new t8_geometry_analytic (2, "analytic moving dim=2",
                                t8_analytic_moving, NULL, NULL, &time);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    snprintf (vtuname, BUFSIZ, "forest_moving_lvl_%i", level);
    break;
  case T8_GEOM_OCC_SURFACE_CUBES:
    {
      #if T8_WITH_OCC
      t8_global_productionf
      ("Creating uniform level %i forests with a occ surface geometry.\n",
       level);

      /* Create a OCC surface */
      Handle_Geom_Surface occ_surface;
      TColgp_Array2OfPnt point_array(1, 5, 1, 3);
      
      point_array(1, 1) = gp_Pnt(-0.2, 0.1, 1.2);
      point_array(2, 1) = gp_Pnt(0.5, 0, 1.4);
      point_array(3, 1) = gp_Pnt(1.0, -0.2, 1.1);
      point_array(4, 1) = gp_Pnt(1.5, 0, 1.0);
      point_array(5, 1) = gp_Pnt(2.1, -0.2, 0.9);

      point_array(1, 2) = gp_Pnt(0.0, 0.5, 1.0);
      point_array(2, 2) = gp_Pnt(0.5, 0.5, 1.2);
      point_array(3, 2) = gp_Pnt(1.0, 0.5, 1.0);
      point_array(4, 2) = gp_Pnt(1.5, 0.5, 0.8);
      point_array(5, 2) = gp_Pnt(2.2, 0.5, 0.6);

      point_array(1, 3) = gp_Pnt(0.0, 1, 0.8);
      point_array(2, 3) = gp_Pnt(0.5, 1, 0.8);
      point_array(3, 3) = gp_Pnt(1.0, 0.85, 0.9);
      point_array(4, 3) = gp_Pnt(1.5, 1, 1.1);
      point_array(5, 3) = gp_Pnt(2.0, 1.1, 1.2);

      occ_surface = GeomAPI_PointsToBSplineSurface(point_array).Surface();
      t8_global_occ_surface[0] = occ_surface;

      geometry = new t8_geometry_occ (3, "occ surface dim=3", NULL);
      
      /* Create tree 0*/
      t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
      double vertices0[24] = {
        -0.5, 0, 0,
        1, 0, 0,
        -0.25, 1.25, 0,
        1, 1, 0,
        0, 0, 1,
        1, 0, 1,
        0, 1, 1,
        1, 1, 1
      };
      t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 24);

      /* Give tree 1 information about its surface and the parameters of the vertices*/
      double parameters0[8] = {0, 0,
                              0.5, 0,
                              0, 1,
                              0.5, 1};
      int faces[6] = {-1, -1, -1, -1, -1, 0};
      t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id(), T8_CMESH_OCC_SURFACE_ATTRIBUTE_KEY, 
                              faces, 6 * sizeof(int), 0);
      t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 5, 
                              parameters0, 8 * sizeof(double), 0);

      /* Create tree 1 */
      t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_HEX);
      double vertices1[24] = {
        1, 0, 0,
        2, 0, 0,
        1, 1, 0,
        2, 1, 0,
        1, 0, 1,
        2, 0, 1,
        1, 1, 1,
        2, 1, 1
      };
      t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 24);

      /* Give tree 1 information about its surface and the parameters of the vertices */
      double parameters1[8] = {0.5, 0,
                              1, 0,
                              0.5, 1,
                              1, 1};
      int face1 = 5;
      t8_cmesh_attribute_occ_surface occ_surface_attribute1(0, face1);
      t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id(), T8_CMESH_OCC_SURFACE_ATTRIBUTE_KEY, 
                              faces, 6 * sizeof(int), 0);
      t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 5, 
                              parameters1, 8 * sizeof(double), 0);

      /* Join tree 0 and tree 1 together */
      t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
      
      snprintf (vtuname, BUFSIZ, "forest_occ_surface_cubes_lvl_%i", level);
      break;
      #else /* !T8_WITH_OCC */
      SC_ABORTF("OCC not linked");
      #endif /* T8_WITH_OCC */
    }
  case T8_GEOM_OCC_SURFACE_CYLINDER:
    {
      #if T8_WITH_OCC
      t8_global_productionf
      ("Creating uniform level %i forests with a occ cylinder geometry.\n",
       level);

      /* Create occ cylinder surfaces */
      double radius_inner = 0.5;
      double radius_outer = 1.0;
      gp_Pnt origin(0, 0, 0);
      gp_Dir z_dir(0, 0, 1);
      gp_Ax2 axis(origin, z_dir);
      gp_Vec height(0, 0, 1);
      gp_Circ circle_outer(axis, radius_outer);
      gp_Circ circle_inner(axis, radius_inner);
      BRepBuilderAPI_MakeEdge make_outer_edge(circle_outer);
      TopoDS_Edge edge_outer = make_outer_edge.Edge();
      TopoDS_Face face_outer = TopoDS::Face(BRepPrimAPI_MakePrism(edge_outer, height));
      Handle_Geom_Surface cylinder_outer = BRep_Tool::Surface(face_outer);
      BRepBuilderAPI_MakeEdge make_inner_edge(circle_inner);
      TopoDS_Edge edge_inner = make_inner_edge.Edge();
      TopoDS_Face face_inner = TopoDS::Face(BRepPrimAPI_MakePrism(edge_inner, height));
      Handle_Geom_Surface cylinder_inner = BRep_Tool::Surface(face_inner);
      t8_global_occ_surface[0] = cylinder_outer;
      t8_global_occ_surface[1] = cylinder_inner;

      geometry = new t8_geometry_occ (3, "occ surface dim=3", NULL);      
      
      int num = 4;
      double *vertices, *parameters;
      vertices = T8_ALLOC(double, num * 24);
      parameters = T8_ALLOC(double, num * 8);
      int faces[6] = {0, 1, -1, -1, -1, -1};
      for (int i = 0; i < num; ++i)
      {
        t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_HEX);
        vertices[i * 24 + 0] = cos((i + 1) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 1] = sin((i + 1) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 2] = 0;
        vertices[i * 24 + 3] = cos((i + 1) * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 4] = sin((i + 1) * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 5] = 0;
        vertices[i * 24 + 6] = cos(i * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 7] = sin((i) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 8] = 0;
        vertices[i * 24 + 9] = cos(i * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 10] = sin(i * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 11] = 0;
        vertices[i * 24 + 12] = cos((i + 1) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 13] = sin((i + 1) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 14] = 1;
        vertices[i * 24 + 15] = cos((i + 1) * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 16] = sin((i + 1) * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 17] = 1;
        vertices[i * 24 + 18] = cos(i * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 19] = sin((i) * 2 * M_PI / num) * radius_outer;
        vertices[i * 24 + 20] = 1;
        vertices[i * 24 + 21] = cos(i * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 22] = sin(i * 2 * M_PI / num) * radius_inner;
        vertices[i * 24 + 23] = 1;
        t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 24, 24);
        parameters[i * 8 + 0] = (i + 1) * 2 * M_PI / num;
        parameters[i * 8 + 1] = 0;
        parameters[i * 8 + 2] = i * 2 * M_PI / num;
        parameters[i * 8 + 3] = 0;
        parameters[i * 8 + 4] = (i + 1) * 2 * M_PI / num;
        parameters[i * 8 + 5] = -1;
        parameters[i * 8 + 6] = i * 2 * M_PI / num;
        parameters[i * 8 + 7] = -1;
        t8_cmesh_set_attribute (cmesh, i, t8_get_package_id(), T8_CMESH_OCC_SURFACE_ATTRIBUTE_KEY, 
                                faces, 6 * sizeof(int), 1);
        t8_cmesh_set_attribute (cmesh, i, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 0, 
                                parameters + i * 8, 8 * sizeof(double), 1);
        t8_cmesh_set_attribute (cmesh, i, t8_get_package_id(), T8_CMESH_OCC_SURFACE_PARAMETERS_ATTRIBUTE_KEY + 1, 
                                parameters + i * 8, 8 * sizeof(double), 1);
      }
      
      T8_FREE(vertices);
      T8_FREE(parameters);
      snprintf (vtuname, BUFSIZ, "forest_geometry_cylinder_lvl_%i", level);
      break;
      #else /* !T8_WITH_OCC */
      SC_ABORTF("OCC not linked");
      #endif /* T8_WITH_OCC */
    }
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* Register the geometry */
  t8_cmesh_register_geometry (cmesh, geometry);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* The initial uniform refinement level is the input level except
   * when geom_type is T8_GEOM_CIRCLE. In that case we start with level
   * 2 and refine recursively only along the boundary. */
  uniform_level = geom_type == T8_GEOM_CIRCLE ? SC_MIN (2, level) : level;
  /* Create a uniform forest */
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level,
                           0, sc_MPI_COMM_WORLD);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_forest_t         forest_adapt;
    /* Create a forest that is only refined at the tree boundaries. 
     * We pass the input level as user pointer and use it in the adapt 
     * callback to stop refinement after this level. */
    forest_adapt =
      t8_forest_new_adapt (forest, t8_geom_adapt_boundary, 1, 1, &level);
    forest = forest_adapt;
  }

  /* Write to vtk */
  t8_forest_write_vtk (forest, vtuname);
  /* Output */
  t8_global_productionf ("Wrote forest to vtu files %s.*\n", vtuname);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_global_productionf
      ("\tNote that this mesh is heavily distorted and we do not\n");
    t8_global_productionf
      ("\trecommend using such a mesh in a production code.\n");
    t8_global_productionf
      ("\tThis example is for demonstrative purposes only.\n");
  }
  if (geom_type == T8_GEOM_MOVING) {
    /* Moving geometry, we start a time simulation and write out the mesh
     * after each time step. */
    int                 timestep = 0;
    const int           num_timesteps = 100;
    const int           end_time = 4;
    char                vtuname_with_timestep[BUFSIZ];

    for (timestep = 0; timestep < num_timesteps; ++timestep) {
      /* Modify the time. Note that the user_data pointer of our
       * geometry points to this entry, which changes the shape of the tree. */
      time += ((double) end_time) / num_timesteps;
      /* At the time step to the output filename */
      sreturn =
        snprintf (vtuname_with_timestep, BUFSIZ, "%s_%04i", vtuname,
                  timestep);
      if (sreturn >= BUFSIZ) {
        /* The vtu name message was truncated */
        /* Note: gcc >= 7.1 prints a warning if we
         * do not check the return value of snprintf. */
        t8_debugf ("Warning: Truncated vtu name to '%s'\n",
                   vtuname_with_timestep);
      }

      t8_forest_write_vtk (forest, vtuname_with_timestep);
    }
  }

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 level;
  int                 parsed, helpme;
  int                 geom_type;
  int                 sreturn;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the analytic geometry capabitlities of t8code.\n"
                      "You can choose from different geometries on which to build a uniform forest.\n"
                      "Usage: %s\n", usage);

  if (sreturn >= BUFSIZ) {
    /* The help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated help message to '%s'\n", help);
  }

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The refinement level of the mesh.");
  sc_options_add_int (opt, 'g', "geometry", &geom_type, 0,
                      "Specify the geometry to use.\n"
                      "\t\t0 - The graph of sin(x) * cos (y) with two 2D quad trees.\n"
                      "\t\t1 - A cylinder with one 2D quad tree.\n"
                      "\t\t2 - A moebius strip on a hybrid mesh with 4 triangles and 2 quads.\n"
                      "\t\t3 - A square of two triangles that is mapped into a circle.\n"
                      "\t\t    The mesh will not be uniform. Instead it is refined at the domain boundary.\n"
                      "\t\t4 - A cube that is distorted in z-direction with one 3D cube tree.\n"
                      "\t\t5 - A moving mesh consisting of a single 2D quad tree.\n");

  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && T8_GEOM_ZERO <= geom_type
           && geom_type < T8_GEOM_COUNT) {
    t8_analytic_geom (level, (t8_analytic_geom_type) geom_type);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
