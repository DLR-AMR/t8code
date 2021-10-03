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

/** \file t8_geometry_occ.hxx
 * This gometry implements OpenCASCADE geometries. It enables the option to link different 
 * 1 and 2 dimensional occ geometries to the edges and faces of refinement trees. 
 * The geometry of the refinement tree gets bent accordingly.
 */

#ifndef T8_GEOMETRY_OCC_HXX
#define T8_GEOMETRY_OCC_HXX

#define T8_WITH_OCC 1
#if T8_WITH_OCC

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_cmesh/t8_cmesh_types.h>

#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Standard_Handle.hxx>



extern Handle_Geom_Curve t8_global_occ_curve[];
extern Handle_Geom_Surface t8_global_occ_surface[];

/**
 * Definition of an occ geometry function.
 * This function maps reference coordinates to physical
 * coordinates.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 * \param [in]  tree_data   The data of the current tree as loaded by a \ref t8_geom_load_tree_data_fn.
 * \param [in]  user_data   The user data pointer stored in the geometry.
 */
typedef void        (*t8_geom_occ_fn) (t8_cmesh_t cmesh,
                                      t8_gloidx_t gtreeid,
                                      const double *ref_coords,
                                      double out_coords[3],
                                      const void *tree_data,
                                      const void *user_data);

struct t8_geometry_occ:public t8_geometry_w_vertices
{
public:

  /**
   * Constructor with analytical and jacobian functions.
   * \param [in] dimension  The dimension of this geometry.
   * \param [in] name       The name to give this geometry.
   * \param [in] load_tree_data The function that is used to load a tree's data.
   */
  t8_geometry_occ (int dimension, const char *name);

  /** The destructor. 
   * Clears the allocated memory.
   */
                      virtual ~ t8_geometry_occ ()
  {
    /* Nothing to do */
  }

  /**
   * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  virtual void        t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const;

  /**
   * Not yet implemented.
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   * \note The jacobian will be
   *            (1)              (1 0)             (1 0 0)
   * dim 1: J = (0)   dim 2: J = (0 1)  dim 3: J = (0 1 0)
   *            (0)              (0 0)             (0 0 1)
   */
  virtual void        t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const;


  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the vertex coordinates of the 
   * tree into an internal buffer (as is done in the linear geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void        t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid);

  /** Push a new occ curve into the occ_curves array. The position gets saved in index.
   * \param [in]  curve      The occ curve.
   * \param [out] index      The index of the curve in the occ_curves array.
   */
  virtual void        t8_geom_push_occ_curve (Handle_Geom_Curve curve, 
                                              int index) const;

  /** Push a new occ surface into the occ_surfaces array. The position gets saved in index.
   * \param [in]  surface    The occ surface.
   * \param [out] index      The index of the surface in the occ_surfaces array.
   */
  virtual void        t8_geom_push_occ_surface (Handle_Geom_Surface surface, 
                                              int index) const;

private:

  const void         *tree_data;        /** Tree data pointer that can be set in \a load_tree_data and 
                                           is passed onto \a analytical_function and \a jacobian. */

  sc_array_t         *occ_curves;       /** Occ curve geometry pointer. Curves can be pushed with t8_push_occ_curve(). */

  sc_array_t         *occ_surfaces;     /** Occ surface geometry pointer. Curves can be pushed with t8_push_occ_surface(). */
};

#endif /* T8_WITH_OCC */

#endif /* !T8_GEOMETRY_OCC_HXX! */
