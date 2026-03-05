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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#if T8_ENABLE_OCC
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#endif
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>

#if T8_ENABLE_OCC
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <gp_Pnt.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Circ.hxx>
#include <gp_Vec.hxx>
#include <BRep_Tool.hxx>
#endif

typedef enum {
  T8_GEOM_ZERO = 0,
  T8_GEOM_SINCOS = T8_GEOM_ZERO,
  T8_GEOM_CYLINDER,
  T8_GEOM_MOEBIUS,
  T8_GEOM_TWO_GEOMETRIES,
  T8_GEOM_CIRCLE,
  T8_GEOM_3D,
  T8_GEOM_MOVING,
  T8_GEOM_ANALYTIC_QUAD_TO_SPHERE,
  T8_GEOM_CAD_TRIANGLE,
  T8_GEOM_CAD_CURVE_CUBE,
  T8_GEOM_CAD_SURFACE_CUBES,
  T8_GEOM_CAD_SURFACE_CYLINDER,
  T8_GEOM_COUNT
} t8_example_geom_type;

/** This geometry maps a point (x,y) in R^2 
 * to the point (x,y, 0.2 * sin(2PI X) * cos(2PI Y)).
 * It should only be used for 2 dimensional forests.
 * 
 * This geometry does not provide a jacobian.
 */
struct t8_geometry_sincos: public t8_geometry
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry_sincos (): t8_geometry ("t8_sincos_geometry")
  {
  }

  /**
   * Maps points (x,y) in R^2 
   * to the point (x,y, 0.2 * sin(2PI X) * cos(2PI Y)).
   * It is specifically designed to work on two tree cmeshes and 
   * models the rectangle [0,2] x [0,1].
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_2d = 2 * i_coord;
      const int offset_3d = 3 * i_coord;
      double x = ref_coords[offset_2d];
      if (gtreeid == 1) {
        /* Translate ref coordinates by +1 in x direction for the second tree. */
        x += 1;
      }
      out_coords[offset_3d] = x;
      out_coords[offset_3d + 1] = ref_coords[offset_2d + 1];
      out_coords[offset_3d + 2]
        = 0.2 * sin (ref_coords[offset_2d] * 2 * M_PI) * cos (ref_coords[offset_2d + 1] * 2 * M_PI);
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /** Check if  the currently active tree has a negative volume. In this case return zero. */
  bool
  t8_geom_tree_negative_volume () const
  {
    return 0;
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf (
        "t8_geometry_sincos is not compatible with tree type %s\n It is only compatible with quad elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }
};

/** This geometry maps the unit square \f$ [0,1]^2 \f$ to the moebius strip.
 * The unit square can be modelled with any cmesh (consisting of any number of trees).
 * 
 * It inherits from the w_vertices geometry since we use the tree's vertex coordinates.
 * This geometry does not provide a jacobian.
 */
struct t8_geometry_moebius: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry_moebius (): t8_geometry_with_vertices ("t8_moebius_geometry")
  {
  }

  /**
   * Maps points in \f$ [0,1]^2 \f$ to the moebius band.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    double t;
    double phi;

    /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store in out_coords. */
    t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_3d = i_coord * 3;
      /* At first, we map x from [0,1] to [-.5,.5]
      * and y to [0, 2*PI] */
      t = out_coords[offset_3d] - .5;
      phi = out_coords[offset_3d + 1] * 2 * M_PI;
      /* We now apply the parametrization for the moebius strip. */
      out_coords[offset_3d] = (1 - t * sin (phi / 2)) * cos (phi);
      out_coords[offset_3d + 1] = (1 - t * sin (phi / 2)) * sin (phi);
      out_coords[offset_3d + 2] = t * cos (phi / 2);
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf (
        "t8_geometry_moebius is not compatible with tree type %s\n It is only compatible with quad elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }

  /* Load tree data is inherited from vertices geometry. */
};

/** This geometry maps the unit square to a cylinder.
 * It should only be used for cmeshes with a single quad tree.
 * 
 * This geometry does not provide a jacobian.
 */
struct t8_geometry_cylinder: public t8_geometry
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry_cylinder (): t8_geometry ("t8_cylinder_geometry")
  {
  }

  /**
   * Map a reference point in the unit square to a cylinder.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_3d = i_coord * 3;
      const int offset_2d = i_coord * 2;
      out_coords[offset_3d] = cos (ref_coords[offset_2d] * 2 * M_PI);
      out_coords[offset_3d + 1] = ref_coords[offset_2d + 1];
      out_coords[offset_3d + 2] = sin (ref_coords[offset_2d] * 2 * M_PI);
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /** Check if  the currently active tree has a negative volume. In this case return zero. */
  bool
  t8_geom_tree_negative_volume () const
  {
    return 0;
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf (
        "t8_geometry_cylinder is not compatible with tree type %s\n It is only compatible with quad elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }
};

/**
 * This geometry map a unit square \f$ [0,1]^2 \f$ cmesh to a circle with midpoint 0
 * and radius 1.
 * This geometry massively distorts elements near the boundary and should not be
 * used for actual numerical experiments.
 * 
 * This geometry does not provide a jacobian.
 */
struct t8_geometry_circle: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry_circle (): t8_geometry_with_vertices ("t8_circle_geometry")
  {
  }

  /**
   * Map a reference point in the unit square to a circle.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    double x;
    double y;

    /* Compute the linear coordinates (in [0,1]^2) of the reference vertex and store in out_coords. */

    /* No idea why, but indent insert a lot of newlines here */
    t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);

    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_3d = i_coord * 3;
      /* We now remap the coords to match the square [-1,1]^2 */
      x = out_coords[offset_3d] * 2 - 1;
      y = out_coords[offset_3d + 1] * 2 - 1;

      /* An now we apply the formula that projects the square to the circle. */
      out_coords[offset_3d] = x * sqrt (1 - y * y / 2);
      out_coords[offset_3d + 1] = y * sqrt (1 - x * x / 2);
      out_coords[offset_3d + 2] = 0;
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf (
        "t8_geometry_circle is not compatible with tree type %s\n It is only compatible with quad elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }

  /* Load tree data is inherited from vertices geometry. */
};

/* This geometry rotates \f$ [0,1]^2 \f$ with time around the origin.
 * The rotation direction is reversed after 2 seconds.
 * Additionally, the z coordinate is modified according to the
 * sincos function and multiplied with the current time.
 * To use this, a pointer to a double variable time is added to the geometry.
 * This variable can be modified from outside.
 * 
 * The geometry can only be used with single tree cmeshes (unit square).
 */

struct t8_geometry_moving: public t8_geometry
{
 public:
  /* Basic constructor that sets the name and the time pointer. */
  t8_geometry_moving (const double *time): t8_geometry ("t8_moving_geometry"), ptime (time)
  {
  }

  /**
   * Map a reference point in the unit square to a square distorted with time.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    double x, y, radius_sqr, phi, rho;
    int sign;
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_3d = i_coord * 3;
      const int offset_2d = i_coord * 2;
      x = ref_coords[offset_2d] - .5;
      y = ref_coords[offset_2d + 1] - .5;
      const double time = *ptime;
      radius_sqr = x * x + y * y;
      phi = radius_sqr * (time > 2 ? 4 - time : time);

      /* Change gridlines by applying a 4th order polynomial mapping
      * [0,1]^2 -> [0,1]^2.
      * And then map this to [-0.5,-0.5]^2 */
      sign = x < 0 ? 1 : -1;
      rho = 0.5 - time / 10;
      x = sign * (1 - exp (-fabs (-x) / rho)) / (2 * (1 - exp (-0.5 / rho)));
      sign = y < 0 ? 1 : -1;
      y = sign * (1 - exp (-fabs (-y) / rho)) / (2 * (1 - exp (-0.5 / rho)));

      /* Rotate the x-y axis and add sincos in z axis. */
      out_coords[offset_3d] = x * (cos (phi)) - y * sin (phi);
      out_coords[offset_3d + 1] = y * (cos (phi)) + x * sin (phi);
      out_coords[offset_3d + 2] = 0;
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /** Check if  the currently active tree has a negative volume. In this case return zero. */
  bool
  t8_geom_tree_negative_volume () const
  {
    return 0;
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf (
        "t8_geometry_moving is not compatible with tree type %s\n It is only compatible with quad elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }

 protected:
  const double *ptime; /* Time pointer to outside time variable */
};

/** Map the unit cube \f$ [0,1]^3 \f$ onto a cube that is distorted
 * in z direction.
 * Can be used with 1 tree unit cube cmesh only.
 */
struct t8_geometry_cube_zdistorted: public t8_geometry
{
 public:
  /* Basic constructor that sets the name. */
  t8_geometry_cube_zdistorted (): t8_geometry ("t8_cube_zdistorted_geometry")
  {
  }
  /**
   * Map a reference point in the unit cube to a cube distorted in the z axis.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                    const size_t num_coords, double *out_coords) const
  {
    for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
      const int offset_3d = i_coord * 3;
      out_coords[offset_3d] = ref_coords[offset_3d];
      out_coords[offset_3d + 1] = ref_coords[offset_3d + 1];
      out_coords[offset_3d + 2]
        = ref_coords[offset_3d + 2]
          * (0.8 + 0.2 * sin (ref_coords[offset_3d] * 2 * M_PI) * cos (ref_coords[offset_3d + 1] * 2 * M_PI));
    }
  }

  /* Jacobian, not implemented. */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /** Check if  the currently active tree has a negative volume. In this case return zero. */
  bool
  t8_geom_tree_negative_volume () const
  {
    return 0;
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only hex elements are supported by this geometry.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_HEX) {
      t8_productionf (
        "t8_geometry_cube_zdistorted is not compatible with tree type %s\n It is only compatible with hex elements.\n",
        t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_UNDEFINED;
  }
};

/* This adapt callback function will refine all elements at the
 * domain boundary up to a given maximum refinement level. */
static int
t8_geom_adapt_boundary (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, const t8_eclass_t tree_class,
                        [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                        [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                        t8_element_t *elements[])
{
  t8_cmesh_t cmesh = t8_forest_get_cmesh (forest_from);
  /* Get the number of faces of the element. */
  const int num_faces = scheme->element_get_num_faces (tree_class, elements[0]);
  int iface;
  /* Get the maximum level from the forest's user data 
   * (must be set before using the callback). */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);

  /* We do not refine more then the given maximum level. */
  if (scheme->element_get_level (tree_class, elements[0]) >= maxlevel) {
    return 0;
  }

  /* Check for each face of the element whether it lies on the 
   * domain boundary. If so, the element is refined. */
  for (iface = 0; iface < num_faces; ++iface) {
    if (scheme->element_is_root_boundary (tree_class, elements[0], iface)) {
      /* This element's face is at its tree boundary. Check whether
         the tree's face is at the domain boundary. */
      int tree_face = scheme->element_get_tree_face (tree_class, elements[0], iface);
      t8_locidx_t lctreeid = t8_forest_ltreeid_to_cmesh_ltreeid (forest_from, ltree_id);
      if (t8_cmesh_tree_face_is_boundary (cmesh, lctreeid, tree_face)) {
        /* The tree's face is at the domain boundary, we refine the element. */
        return 1;
      }
    }
  }
  /* All other elements remain unchanged. */
  return 0;
}

void
quad_to_sphere_callback ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                         const double *ref_coords, const size_t num_coords, double *out_coords,
                         [[maybe_unused]] const void *tree_data, [[maybe_unused]] const void *user_data)
{
  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset = 3 * i_coord;

    const double radius = 1.0;
    const double latitude = 2 * M_PI * ref_coords[offset + 0];
    const double longitude = ref_coords[offset + 1] * M_PI;

    out_coords[offset + 0] = radius * sin (longitude) * cos (latitude);
    out_coords[offset + 1] = radius * sin (longitude) * sin (latitude);
    out_coords[offset + 2] = radius * cos (longitude);
  }
}

static void
t8_analytic_geom (int level, t8_example_geom_type geom_type)
{
  t8_forest_t forest;
  t8_cmesh_t cmesh;
  char vtuname[BUFSIZ];
  int uniform_level;
  double time = 0; /* used for moving geometry */
  int sreturn;

  t8_cmesh_init (&cmesh);
  /* Depending on the geometry type, add the tree, set the geometry
   * and set the output file name. */
  switch (geom_type) {
  case T8_GEOM_SINCOS:
    t8_global_productionf ("Creating uniform level %i forest with a sine/cosine geometry.\n", level);
    /* Sin/cos geometry. Has two quad trees. */
    t8_cmesh_register_geometry<t8_geometry_sincos> (cmesh);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
    snprintf (vtuname, BUFSIZ, "forest_sincos_lvl_%i", level);
    break;
  case T8_GEOM_CYLINDER:
    t8_global_productionf ("Creating uniform level %i forest with a cylinder geometry.\n", level);
    /* Cylinder geometry. Has one quad tree that is periodic in x direction. */
    t8_cmesh_register_geometry<t8_geometry_cylinder> (cmesh);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
    snprintf (vtuname, BUFSIZ, "forest_cylinder_lvl_%i", level);
    break;
  case T8_GEOM_MOEBIUS:
    t8_global_productionf ("Creating uniform level %i forest with a moebius geometry.\n", level);
    {
      /* Moebius geometry on hybrid unit square. */
      t8_cmesh_t hybrid_square = t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
      t8_cmesh_set_derive (cmesh, hybrid_square);
      t8_cmesh_register_geometry<t8_geometry_moebius> (cmesh);
      snprintf (vtuname, BUFSIZ, "forest_moebius_lvl_%i", level);
    }
    break;
  case T8_GEOM_TWO_GEOMETRIES:
    t8_global_productionf ("Creating uniform level %i forest with a cylinder and a sine cosine geometry.\n", level);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    /* Tree 0 is connected to itself to form a cylinder */
    t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
    /* Cylinder geometry on tree 0. Sincos geometry on tree 1. */
    {
      t8_geometry *geometry_cylinder = t8_cmesh_register_geometry<t8_geometry_cylinder> (cmesh);
      t8_geometry *geometry_sincos = t8_cmesh_register_geometry<t8_geometry_sincos> (cmesh);
      t8_cmesh_set_tree_geometry (cmesh, 0, geometry_cylinder);
      t8_cmesh_set_tree_geometry (cmesh, 1, geometry_sincos);
    }
    snprintf (vtuname, BUFSIZ, "forest_cylinder_and_sincos_lvl_%i", level);
    break;
  case T8_GEOM_CIRCLE:
    t8_global_productionf ("Creating forest with a circle geometry.\n");
    t8_global_productionf ("This forest will get refined at the boundary to level %i.\n", level);
    {
      /* Circle geometry on triangulated unit square. */
      t8_cmesh_t tri_square = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
      t8_cmesh_set_derive (cmesh, tri_square);
      t8_cmesh_register_geometry<t8_geometry_circle> (cmesh);
      snprintf (vtuname, BUFSIZ, "forest_circle_lvl_%i", level);
    }
    break;
  case T8_GEOM_3D:
    t8_global_productionf ("Creating uniform level %i forest with a 3D function graph geometry.\n", level);
    /* Cube geometry with sincos on top. Has one hexahedron tree. */
    t8_cmesh_register_geometry<t8_geometry_cube_zdistorted> (cmesh);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    snprintf (vtuname, BUFSIZ, "forest_cube_3D_lvl_%i", level);
    break;
  case T8_GEOM_MOVING:
    t8_global_productionf ("Creating uniform level %i forest with a moving geometry.\n", level);
    /* Quad geometry that rotates with time. */
    t8_cmesh_register_geometry<t8_geometry_moving> (cmesh, &time);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    snprintf (vtuname, BUFSIZ, "forest_moving_lvl_%i", level);
    break;
  case T8_GEOM_ANALYTIC_QUAD_TO_SPHERE:
    t8_global_productionf ("Wrapping a quad around a sphere.\n");
    t8_cmesh_register_geometry<t8_geometry_analytic> (cmesh, "geom_quad_to_sphere", quad_to_sphere_callback, nullptr,
                                                      nullptr, nullptr, nullptr, nullptr);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_join (cmesh, 0, 0, 1, 0, 0);

    snprintf (vtuname, BUFSIZ, "forest_quad_to_sphere");
    break;
  case T8_GEOM_CAD_TRIANGLE: {
#if T8_ENABLE_OCC
    t8_global_productionf ("Creating uniform level %i forests with an cad triangle geometry.\n", level);

    /* Constructing a triangle with one curved edge (f1) */
    Handle_Geom_BSplineCurve cad_curve;
    TColgp_Array1OfPnt point_array (1, 3);
    TopoDS_Shape shape;

    /* Define knots along the bsplines. */
    point_array (1) = gp_Pnt (0.0, 0.0, 0.0);
    point_array (2) = gp_Pnt (0.4, 1.3, 0.0);
    point_array (3) = gp_Pnt (1.0, 2.0, 0.0);

    /* Generate bsplines from arrays. */
    cad_curve = GeomAPI_PointsToBSpline (point_array).Curve ();

    /* Fill shape with bsplines so that we can create a geometry with this shape. */
    shape = BRepBuilderAPI_MakeEdge (cad_curve).Edge ();

    /* Create a cad geometry. */
    t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, shape);

    /* The arrays indicate which face/edge carries a geometry. 
       * 0 means no geometry and any other number indicates the position of the geometry 
       * in the global geometry array. Here edge 1 carries the created cad_curve. */
    int faces[1] = { 0 };
    int edges[6] = { 0, 1, 0, 0, 0, 0 };
    /* Create tree 0 */
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
    double vertices[9] = { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 2.0, 0.0 };
    t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);

    /* The valid parameter range for bsplines is [0, 1]. Therefore, we define the parameter range accordingly. */
    double parameters_edge[2] = { 0, 1 };

    /* Give the tree information about its curves and the parameters of the vertices. 
       * Each parameter set is given to the tree via its attribute key + the edge or face index it corresponds with. */
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces, 1 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges, 6 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + 1,
                            parameters_edge, 2 * sizeof (double), 0);

    snprintf (vtuname, BUFSIZ, "forest_cad_triangle_lvl_%i", level);
    break;
#else  /* !T8_ENABLE_OCC */
    SC_ABORTF ("OCC not linked");
#endif /* T8_ENABLE_OCC */
  }
  case T8_GEOM_CAD_CURVE_CUBE: {
#if T8_ENABLE_OCC
    t8_global_productionf ("Creating uniform level %i forests with cad curve geometries.\n", level);

    /* Create two cad bsplines which oscillate along the x-axis. 
       * For this we need to define two arrays from which we create the bsplines. */
    Handle_Geom_Curve cad_curve0;
    Handle_Geom_Curve cad_curve1;
    TColgp_Array1OfPnt point_array0 (1, 5);
    TColgp_Array1OfPnt point_array1 (1, 5);
    TopoDS_Shape shape;

    /* Define knots along the bsplines. */
    point_array0 (1) = gp_Pnt (0, 0, 0);
    point_array0 (2) = gp_Pnt (0.25, 0.1, 0.1);
    point_array0 (3) = gp_Pnt (0.5, 0, 0);
    point_array0 (4) = gp_Pnt (0.75, -0.1, -0.1);
    point_array0 (5) = gp_Pnt (1, 0, 0);

    point_array1 (1) = gp_Pnt (0, 1, 1);
    point_array1 (2) = gp_Pnt (0.25, 0.9, 1.1);
    point_array1 (3) = gp_Pnt (0.5, 1, 1);
    point_array1 (4) = gp_Pnt (0.9, 1.1, 0.9);
    point_array1 (5) = gp_Pnt (1, 1, 1);

    /* Generate bsplines from arrays. */
    cad_curve0 = GeomAPI_PointsToBSpline (point_array0).Curve ();
    cad_curve1 = GeomAPI_PointsToBSpline (point_array1).Curve ();

    /* Fill shape with bsplines so that we can create a geometry with this shape. */
    shape = BRepBuilderAPI_MakeEdge (cad_curve0).Edge ();
    shape = BRepAlgoAPI_Fuse (shape, BRepBuilderAPI_MakeEdge (cad_curve1).Edge ());

    /* Create a cad geometry. */
    t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, shape);

    /* The arrays indicate which face/edge carries a geometry. 
     * 0 means no geometry and any other number indicates the position of the geometry 
     * in the global geometry array. Here edge 0 carries cad_curve0 and edge 3 carries cad_curve1.
     * We add them in the next step. */
    int faces[6] = { 0, 0, 0, 0, 0, 0 };
    int edges[24] = { 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    /* Create tree 0 */
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    double vertices[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
    t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);

    /* The valid parameter range for bsplines is [0, 1]. We defined the bsplines in such a way, 
       * that parameter 0 and 1 resemble the two vertices of the connected edge. */
    double parameters[2] = { 0, 1 };

    /* Give the tree information about its curves and the parameters of the vertices. 
       * Each parameter set is given to the tree via its attribute key + the edge or face index it corresponds with. */
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces, 6 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges, 24 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + 0, parameters,
                            2 * sizeof (double), 0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_PARAMETERS_ATTRIBUTE_KEY + 3, parameters,
                            2 * sizeof (double), 0);

    snprintf (vtuname, BUFSIZ, "forest_cad_curve_cube_lvl_%i", level);
    break;
#else  /* !T8_ENABLE_OCC */
    SC_ABORTF ("OCC not linked");
#endif /* T8_ENABLE_OCC */
  }
  case T8_GEOM_CAD_SURFACE_CUBES: {
#if T8_ENABLE_OCC
    t8_global_productionf ("Creating uniform level %i forests with a cad surface geometry.\n", level);

    /* Create a cad bspline surface with 2D array of knots */
    Handle_Geom_Surface cad_surface;
    TColgp_Array2OfPnt point_array (1, 5, 1, 3);
    TopoDS_Shape shape;

    /* Filling the 2D surface array with knots. The resulting surface resembles 
     * a surface at the top (face 5) of the trees.
     * Some of the knots have the same position as the vertices of the trees. 
     * These knots are marked with the tree id and vertex index. 
     * We also marked the direction of the u- and v-parameter.
     *
     *  x--> u-parameter
     *  |
     *  v v-parameter
     *
     *     point_array  1       2       3       4       5
     *
     *         1      t0_v6--------t0_v7&t1_v6--------t1_v7
     *                  |               |               |
     *                  |               |               |
     *         2        | tree 0 face 5 | tree 1 face 5 |
     *                  |               |               |
     *                  |               |               |
     *         3      t0_v4--------t0_v5&t1_v4--------t1_v5
     *
     * z-dir
     *    X--> x-dir
     *    |
     *    v
     *    y-dir
     */
    point_array (1, 1) = gp_Pnt (-0.2, -0.2, 1.2);  // t0_v6
    point_array (2, 1) = gp_Pnt (0.5, 0.0, 1.0);
    point_array (3, 1) = gp_Pnt (1.0, -0.2, 0.8);  // t0_v7 & t1_v6
    point_array (4, 1) = gp_Pnt (1.5, 0.0, 1.0);
    point_array (5, 1) = gp_Pnt (2.2, -0.2, 1.2);  // t1_v7

    point_array (1, 2) = gp_Pnt (0.0, 0.5, 1.0);
    point_array (2, 2) = gp_Pnt (0.5, 0.5, 1.0);
    point_array (3, 2) = gp_Pnt (1.0, 0.5, 0.8);
    point_array (4, 2) = gp_Pnt (1.5, 0.5, 1.0);
    point_array (5, 2) = gp_Pnt (2.0, 0.5, 1.0);

    point_array (1, 3) = gp_Pnt (-0.2, 1.2, 1.2);  // t0_v4
    point_array (2, 3) = gp_Pnt (0.5, 1.0, 1.0);
    point_array (3, 3) = gp_Pnt (1.0, 1.2, 0.8);  // t0_v5 & t1_v4
    point_array (4, 3) = gp_Pnt (1.5, 1.0, 1.0);
    point_array (5, 3) = gp_Pnt (2.2, 1.2, 1.2);  // t1_v5

    /* Generate bspline surface from array and fill shape with it
     * so that we can create a geometry with this shape. */
    cad_surface = GeomAPI_PointsToBSplineSurface (point_array).Surface ();
    shape = BRepBuilderAPI_MakeFace (cad_surface, 1e-6).Face ();

    /* The arrays indicate which face/edge carries a geometry. 
     * 0 means no geometry and any other number indicates the position of the geometry 
     * in the global geometry array. Here face 5 carries the surface, we add it in the next step. 
     * There are no geometries linked to the edges, hence all entries are 0. */
    int faces[6] = { 0, 0, 0, 0, 0, 1 };
    int edges[24] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    /* Create cad geometry. */
    t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, shape);

    /* Create tree 0 */
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
    double vertices0[24] = {
      0.0,  0.0,  0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, -0.2, 1.2, 1.2,  // Point (1, 3) from array
      1.0,  1.2,  0.8,                                                               // Point (3, 3) from array
      -0.2, -0.2, 1.2,                                                               // Point (1, 1) from array
      1.0,  -0.2, 0.8                                                                // Point (3, 1) from array
    };
    t8_cmesh_set_tree_vertices (cmesh, 0, vertices0, 8);

    /* The valid parameter range for bspline surfaces is [0,1]^2. We defined the bspline surface in such a way, 
       * that parameters 0, 0.5 and 1 resemble the vertices of the connected surface. */
    double parameters0[8] = { 0, 0, 0.5, 0, 0, 1, 0.5, 1 };

    /* Give tree 0 information about its surface and the parameters of the vertices. 
     * Each parameter set is given to the tree via its attribute key + the edge or face index it corresponds with. 
     */
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces, 6 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges, 24 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 5, parameters0,
                            8 * sizeof (double), 0);

    /* Create tree 1 */
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_HEX);
    double vertices1[24] = {
      1.0, 0.0,  0.0, 2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 1.0, 0.5, 0.8, /* Point (3, 3) from array */
      2.2, 1.2,  1.2,                                                             /* Point (5, 3) from array */
      1.0, -0.2, 0.8,                                                             /* Point (3, 1) from array */
      2.2, -0.2, 1.2                                                              /* Point (5, 1) from array */
    };
    t8_cmesh_set_tree_vertices (cmesh, 1, vertices1, 8);

    /* The valid parameter range for bspline surfaces is [0,1]^2. We defined the bspline surface in such a way,
     * that parameters 0, 0.5 and 1 resemble the vertices of the connected surface. */
    double parameters1[8] = { 0.5, 0, 1, 0, 0.5, 1, 1, 1 };

    /* Give tree 1 information about its surface and the parameters of the vertices. 
     *  Each parameter set is given to the tree via its attribute key + the edge or face index it corresponds with. 
     *  We can use the same edges and faces array, because we link the surface to the same face on tree 1. */
    t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces, 6 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges, 24 * sizeof (int),
                            0);
    t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 5, parameters1,
                            8 * sizeof (double), 0);

    /* Join tree 0 and tree 1 together */
    t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);

    snprintf (vtuname, BUFSIZ, "forest_cad_surface_cubes_lvl_%i", level);
    break;
#else  /* !T8_ENABLE_OCC */
    SC_ABORTF ("OCC not linked");
#endif /* T8_ENABLE_OCC */
  }
  case T8_GEOM_CAD_SURFACE_CYLINDER: {
#if T8_ENABLE_OCC
    t8_global_productionf ("Creating uniform level %i forests with an cad cylinder geometry.\n", level);

    /* Create cad cylinder surfaces. We use an outer radius of 0.5 to get a diameter of 1. */
    double radius_inner = 0.25;
    double radius_outer = 0.5;

    /* Define origin, z-axis and height vector for creating and extruding circles. */
    gp_Pnt origin (0, 0, 0);
    gp_Dir z_dir (0, 0, 1);
    gp_Ax2 axis (origin, z_dir);
    gp_Vec height (0, 0, 1);

    /* Create inner and outer cylinder mantles. */
    gp_Circ circle_outer (axis, radius_outer);
    gp_Circ circle_inner (axis, radius_inner);
    BRepBuilderAPI_MakeEdge make_outer_edge (circle_outer);
    TopoDS_Edge edge_outer = make_outer_edge.Edge ();
    TopoDS_Face face_outer = TopoDS::Face (BRepPrimAPI_MakePrism (edge_outer, height));
    Handle_Geom_Surface cylinder_outer = BRep_Tool::Surface (face_outer);
    BRepBuilderAPI_MakeEdge make_inner_edge (circle_inner);
    TopoDS_Edge edge_inner = make_inner_edge.Edge ();
    TopoDS_Face face_inner = TopoDS::Face (BRepPrimAPI_MakePrism (edge_inner, height));
    Handle_Geom_Surface cylinder_inner = BRep_Tool::Surface (face_inner);
    TopoDS_Shape shape = BRepBuilderAPI_MakeFace (cylinder_outer, 1e-6).Face ();

    /* Fill shape with mantles so that we can create a geometry with this shape. */
    shape = BRepAlgoAPI_Fuse (shape, BRepBuilderAPI_MakeFace (cylinder_inner, 1e-6).Face ());

    /* The arrays indicate which face/edge carries a geometry. 
     * 0 means no geometry and any other number indicates the position of the geometry 
     * in the global geometry array. Here face 0 carries the outer cylinder and face 1 carries the inner cylinder.
     * We add them in the next step. The edges do not have any geometries, hence all entries are 0. */
    int faces[6] = { 1, 2, 0, 0, 0, 0 };
    int edges[24] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    /* Create a cad geometry. */
    t8_cmesh_register_geometry<t8_geometry_cad> (cmesh, shape);

    /* Create corresponding trees and parameters. 
     * Here we create num trees by a coordinate transformation from cylinder to cartesian coordinates. */
    int num = 4;
    double *vertices, *parameters;
    vertices = T8_ALLOC (double, num * 24);
    parameters = T8_ALLOC (double, num * 8);
    for (int i = 0; i < num; ++i) {
      t8_cmesh_set_tree_class (cmesh, i, T8_ECLASS_HEX);
      /* Coordinate transformation. */
      vertices[i * 24 + 0] = cos ((i + 1) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 1] = sin ((i + 1) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 2] = 0;
      vertices[i * 24 + 3] = cos ((i + 1) * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 4] = sin ((i + 1) * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 5] = 0;
      vertices[i * 24 + 6] = cos (i * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 7] = sin ((i) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 8] = 0;
      vertices[i * 24 + 9] = cos (i * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 10] = sin (i * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 11] = 0;
      vertices[i * 24 + 12] = cos ((i + 1) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 13] = sin ((i + 1) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 14] = 1;
      vertices[i * 24 + 15] = cos ((i + 1) * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 16] = sin ((i + 1) * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 17] = 1;
      vertices[i * 24 + 18] = cos (i * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 19] = sin ((i) * 2 * M_PI / num) * radius_outer;
      vertices[i * 24 + 20] = 1;
      vertices[i * 24 + 21] = cos (i * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 22] = sin (i * 2 * M_PI / num) * radius_inner;
      vertices[i * 24 + 23] = 1;
      t8_cmesh_set_tree_vertices (cmesh, i, vertices + i * 24, 8);

      /* Create corresponding parameters for the cylinders. 
       * The parameter range of the cylinders is u ∈ [0, 2 * M_PI] and v ∈ ]inf, -inf[ */
      parameters[i * 8 + 0] = (i + 1) * 2 * M_PI / num;
      parameters[i * 8 + 1] = 0;
      parameters[i * 8 + 2] = i * 2 * M_PI / num;
      parameters[i * 8 + 3] = 0;
      parameters[i * 8 + 4] = (i + 1) * 2 * M_PI / num;
      parameters[i * 8 + 5] = -1;
      parameters[i * 8 + 6] = i * 2 * M_PI / num;
      parameters[i * 8 + 7] = -1;

      /* Give the trees information about their surfaces and the parameters of the vertices. 
       * Each parameter set is given to the tree via its attribute key + face index it corresponds with. 
       * We can use the same edges and faces array, because we link the surface to the same faces on every tree.*/
      t8_cmesh_set_attribute (cmesh, i, t8_get_package_id (), T8_CMESH_CAD_FACE_ATTRIBUTE_KEY, faces, 6 * sizeof (int),
                              1);
      t8_cmesh_set_attribute (cmesh, i, t8_get_package_id (), T8_CMESH_CAD_EDGE_ATTRIBUTE_KEY, edges, 24 * sizeof (int),
                              1);
      t8_cmesh_set_attribute (cmesh, i, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 0,
                              parameters + i * 8, 8 * sizeof (double), 0);
      t8_cmesh_set_attribute (cmesh, i, t8_get_package_id (), T8_CMESH_CAD_FACE_PARAMETERS_ATTRIBUTE_KEY + 1,
                              parameters + i * 8, 8 * sizeof (double), 0);
    }

    T8_FREE (vertices);
    T8_FREE (parameters);
    snprintf (vtuname, BUFSIZ, "forest_geometry_cylinder_lvl_%i", level);
    break;
#else  /* !TT8_ENABLE_OCC */
    SC_ABORTF ("OCC not linked");
#endif /* T8_ENABLE_OCC */
  }
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* The initial uniform refinement level is the input level except
   * when geom_type is T8_GEOM_CIRCLE. In that case we start with level
   * 2 and refine recursively only along the boundary. */
  uniform_level = geom_type == T8_GEOM_CIRCLE ? SC_MIN (2, level) : level;
  /* Create a uniform forest */
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), uniform_level, 0, sc_MPI_COMM_WORLD);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_forest_t forest_adapt;
    /* Create a forest that is only refined at the tree boundaries. 
     * We pass the input level as user pointer and use it in the adapt 
     * callback to stop refinement after this level. */
    forest_adapt = t8_forest_new_adapt (forest, t8_geom_adapt_boundary, 1, 1, &level);
    forest = forest_adapt;
  }

  /* Write to vtk. We use the extended vtk function to export a curved vtk mesh.
   * This is only viable if you link to vtk. */
  t8_forest_write_vtk_ext (forest, vtuname, 1, 1, 1, 1, 0, 1, 0, 0, NULL);
  /* Output */
  t8_global_productionf ("Wrote forest to vtu files %s.*\n", vtuname);
  if (geom_type == T8_GEOM_CIRCLE) {
    t8_global_productionf ("\tNote that this mesh is heavily distorted and we do not\n");
    t8_global_productionf ("\trecommend using such a mesh in a production code.\n");
    t8_global_productionf ("\tThis example is for demonstrative purposes only.\n");
  }
  if (geom_type == T8_GEOM_MOVING) {
    /* Moving geometry, we start a time simulation and write out the mesh
     * after each time step. */
    int timestep = 0;
    const int num_timesteps = 100;
    const double end_time = 4;
    char vtuname_with_timestep[BUFSIZ];

    for (timestep = 0; timestep < num_timesteps; ++timestep) {
      /* Modify the time. Note that a pointer of our
       * geometry points to this entry, which changes the shape of the tree. */
      time += end_time / num_timesteps;
      /* At the time step to the output filename */
      sreturn = snprintf (vtuname_with_timestep, BUFSIZ, "%s_%04i", vtuname, timestep);
      if (sreturn >= BUFSIZ) {
        /* The vtu name message was truncated */
        /* Note: gcc >= 7.1 prints a warning if we
         * do not check the return value of snprintf. */
        t8_debugf ("Warning: Truncated vtu name to '%s'\n", vtuname_with_timestep);
      }

      t8_forest_write_vtk (forest, vtuname_with_timestep);
      t8_debugf ("Wrote vtu file %s\n", vtuname_with_timestep);
    }
  }

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_options_t *opt;
  char usage[BUFSIZ];
  char help[BUFSIZ];
  int level;
  int parsed, helpme;
  int geom_type;
  int sreturn;

  /* brief help message */
  snprintf (usage, BUFSIZ,
            "\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  sreturn = snprintf (help, BUFSIZ,
                      "Demonstrates the some of the geometry capabitlities of t8code.\n"
                      "You can choose from different geometries on which to build a uniform forest.\n"
                      "Usage: %s\n",
                      usage);

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
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 2, "The uniform refinement level of the mesh. Default: 2");
  sc_options_add_int (opt, 'g', "geometry", &geom_type, -1,
                      "Specify the geometry to use.\n"
                      "\t\t0 - The graph of sin(x) * cos(y) with two 2D quad trees.\n"
                      "\t\t1 - A cylinder with one 2D quad tree.\n"
                      "\t\t2 - A moebius strip on a hybrid mesh with 4 triangles and 2 quads.\n"
                      "\t\t3 - A mesh of two trees with different geometries each.\n\t\t    Using the cylinder for the "
                      "first tree, the sin/cos for the second.\n"
                      "\t\t4 - A square of two triangles that is mapped into a circle.\n"
                      "\t\t    The mesh will not be uniform. Instead it is refined at the domain boundary.\n"
                      "\t\t5 - A cube that is distorted in z-direction with one 3D cube tree.\n"
                      "\t\t6 - A moving mesh consisting of a single 2D quad tree.\n"
                      "\t\t7 - A quad morphed into a sphere.\n"
                      "\t\t8 - A curved triangle with an cad curve.\n"
                      "\t\t9 - A cube with two cad curves as edges.\n"
                      "\t\t10 - Two cubes with one cad surface as face.\n"
                      "\t\t11 - A hollow cylinder with a cad surface on the in- and outside.\n");

  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= level && T8_GEOM_ZERO <= geom_type && geom_type < T8_GEOM_COUNT) {
    t8_analytic_geom (level, (t8_example_geom_type) geom_type);
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
