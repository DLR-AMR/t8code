/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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
#include <t8_geometry/t8_geometry_with_vertices.hxx>


struct t8_geometry_healpix: public t8_geometry_with_vertices
{
 public:
  /* Basic constructor that sets the dimension and the name. */
  t8_geometry_healpix (): t8_geometry_with_vertices ("t8_geometry_healpix")
  {
    std::cout<<"done1"<<std::endl;
  }

  /* The destructor. */
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
   * This routine expects an input mesh of six hexaeders arranged into a cube.
   *
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

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
                             [[maybe_unused]] double *jacobian) const
  {
    SC_ABORT_NOT_REACHED ();
  }

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only hex elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    if (active_tree_class != T8_ECLASS_QUAD) {
      t8_productionf ("t8_geometry_healpix is not compatible with tree type %s\n"
                      "It is only compatible with hex elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return false;
    }
    return true;
  }

  /* Load tree data is inherited from t8_geometry_with_vertices. */
};
