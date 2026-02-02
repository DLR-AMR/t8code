/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_geometry_cad.hxx
 * This geometry implements OpenCASCADE geometries. It enables the option to link different
 * 1 and 2 dimensional cad geometries to the edges and faces of refinement trees.
 * The geometry of the refinement tree is extended into the volume accordingly.
 */

#ifndef T8_GEOMETRY_CAD_HXX
#define T8_GEOMETRY_CAD_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.h>
#include <t8_cad/t8_cad.hxx>
#include <memory>
#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>

/**
 * This geometry uses OpenCASCADE CAD geometries to curve
 * the trees to the actual shape of the underlying domain.
 */
struct t8_geometry_cad: public t8_geometry_with_vertices
{
 public:
  /**
   * Constructor of the cad geometry. The geometry
   * is currently viable with quad/hex and triangle trees. Tets will be supported soon.
   * The geometry uses as many vertices as the tree type has, as well as
   * additional geometry information, which is extracted from a .brep file.
   * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
   * Since the internals of this geometry are finely tuned to the .brep file
   * it is recommended to only use it with the \ref t8_cmesh_from_msh_file function.
   * \param [in] fileprefix Prefix of a .brep file from which to extract an cad geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_cad (std::string fileprefix, std::string name = "t8_geom_cad");

  /**
   * Constructor of the cad geometry. The geometry
   * is currently viable with quad/hex and triangle trees. Tets will be supported soon.
   * The geometry uses as many vertices as the tree type has, as well as
   * additional geometry information, which is given via the \a cad_shape.
   * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
   * This constructor can be used in short scripts or in combination with a
   * mesh generator, to omit the file IO of the
   * \ref t8_geometry_cad (std::string fileprefix,  std::string name) constructor.
   * \param [in] cad_shape  cad shape geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_cad (const TopoDS_Shape cad_shape, std::string name = "t8_geom_cad");

  /**
   * Constructor of the cad geometry for testing purposes. Sets an invalid cad_shape.
   */
  t8_geometry_cad ();

  /** The destructor. */
  virtual ~t8_geometry_cad ()
  {
    /* Nothing to do. */
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_CAD;
  };

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  virtual void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \a num_coords x dimension x 3. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian  (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   */
  virtual void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *jacobian) const;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the vertex coordinates of the
   * tree into an internal buffer (as is done in the linear geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * This geometry supports all element types, hence it returns true.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const
  {
    return true;
  }

  /**
   * Getter function for the CAD manager.
   * 
   * \return The CAD manager of the geometry.
  */
  std::shared_ptr<t8_cad>
  get_cad_manager () const
  {
    return cad_manager;
  }

 private:
  /**
   * Maps points in the reference space \f$ [0,1]^2 \f$ to \f$ \mathbb{R}^3 \f$. Only for triangle trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of 2 entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void
  t8_geom_evaluate_cad_tri (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                            double *out_coords) const;

  /**
   * Maps points in the reference space \f$ [0,1]^2 \f$ to \f$ \mathbb{R}^3 \f$. Only for quad trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_quad (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *out_coords) const;

  /**
   * Map a point in the reference space \f$ [0,1]^3 \f$ to \f$ \mathbb R^3 \f$. Only for tet trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  The number of points to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void
  t8_geom_evaluate_cad_tet (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                            double *out_coords) const;

  /**
   * Map a point in the reference space \f$ \f$ [0,1]^3 \f$ \f$ to \f$ \mathbb{R}^3 \f$. Only for hex trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_hex (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                            double *out_coords) const;

  /**
   * Maps points in the reference space \f$ \f$ [0,1]^3 \f$ \f$ to \f$ \mathbb{R}^3 \f$. Only for prism trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_prism (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                              double *out_coords) const;

  /**
   * Evaluate a point on a CAD curve.
   * \param [in] curve_index The index of the curve in the CAD shape.
   * \param [in] param       The parameter on the curve to evaluate.
   * \return The point on the curve.
   */
  gp_Pnt
  process_curve (const int curve_index, const double param) const;

  /**
   * Evaluate a point on a CAD surface.
   * \param [in] surface_index The index of the surface in the CAD shape.
   * \param [in] params        The parameter on the surface to evaluate.
   * \return The point on the surface.
   */
  gp_Pnt
  process_surface (const int surface_index, const double *params) const;

  const int *edges; /**< The linked edges of the currently active tree. */
  const int *faces; /**< The linked faces of the currently active tree. */

  std::shared_ptr<t8_cad> cad_manager; /**< The CAD manager of the geometry. */
};

#endif /* !T8_GEOMETRY_CAD_HXX */
