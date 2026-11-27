/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
#include <t8_geometry/t8_geometry_base.hxx>
/**
* The t8_geometry_healpix struct implements a spherical geometry mapping
* based on the HEALPix (Hierarchical Equal Area isoLatitude Pixelization)
* base-mesh configuration. It inherits from t8_geometry and
* is designed to map planar quadrilateral faces onto a spherical surface.
*/
struct t8_geometry_healpix: public t8_geometry
{
 private:
  double sphere_outer_radius;
  double sphere_inner_radius = 0;

 public:
  /**
   * @brief Constructs a HEALPix geometry with the specified sphere inner and outer radius.
   *
   * @param outer_radius The outer radius of the sphere.
   * @param inner_radius The inner radius of the sphere.
   */
  t8_geometry_healpix (double outer_radius, double inner_radius): t8_geometry ("t8_geometry_healpix")
  {
    if (inner_radius > outer_radius) {
      t8_errorf ("\n[t8_healpix] Warning inner radius is greater than outer radius \n[t8_healpix] defaulting to a "
                 "radius of 1\n");
    }

    sphere_outer_radius = outer_radius;
    sphere_inner_radius = inner_radius;
  }
  /**
   * @brief Default constructor.
   *
   * Constructs a HEALPix geometry with a given sphere radius.
   */
  t8_geometry_healpix (double radius): t8_geometry ("t8_geometry_healpix")
  {
    sphere_outer_radius = radius;
  }
  /**
   * @brief Default constructor.
   *
   * Constructs a HEALPix geometry with a sphere radius of 1.
   */
  t8_geometry_healpix (): t8_geometry ("t8_geometry_healpix")
  {
    sphere_outer_radius = 1;
  }

  /**
   * @brief Destructor.
   */
  virtual ~t8_geometry_healpix ()
  {
  }

  /**
   * Map the faces of a cube to a spherical surface.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   *
   * This routine expects an input mesh of twelve quadrilaterals in a configuration of the healpix base-mesh
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const override;

  /**
   * Jacobian, not implemented.
   * \param[in] cmesh       The cmesh in which the point lies.
   * \param[in] gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param[in] ref_coords  Array of tree dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
   * \param[in] num_coords  The number of points to map.
   * \param[in] jacobian    The Jacobian matrix to be filled.
   */
  void
  t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                             [[maybe_unused]] const double *ref_coords, [[maybe_unused]] const size_t num_coords,
                             [[maybe_unused]] double *jacobian) const override
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only quad elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const override
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf ("t8_geometry_healpix is not compatible with tree type %s\n"
                      "It is only compatible with quad elements.\n",
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
  t8_geom_get_type () const override;

  /**
   * Override to suppress warning.
   * \return                True if the currently loaded tree has a negative volume.
   */
  virtual bool
  t8_geom_tree_negative_volume () const override
  {
    return 0;
  };
};
