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
 * This geometry implements OpenCASCADE geometries. It enables the option to link different 
 * 1 and 2 dimensional occ geometries to the edges and faces of refinement trees. 
 * The geometry of the refinement tree is extended into the volume accordingly.
 */

#ifndef T8_GEOMETRY_OCC_HXX
#define T8_GEOMETRY_OCC_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cad/t8_cad_geom.hxx>

#if T8_WITH_OCC

#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>

#endif /* T8_WITH_OCC */

/** The vertices of each edge of a hexahedron. Used in the occ geometry. */
extern const int    t8_edge_vertex_to_tree_vertex[T8_ECLASS_MAX_EDGES][2];

/** The faces connected to each edge. */
extern const int    t8_edge_to_face[T8_ECLASS_MAX_EDGES][2];

/** The edges of a face to the edges of a tree */
extern const int
     t8_face_edge_to_tree_edge[T8_ECLASS_MAX_FACES][T8_ECLASS_MAX_EDGES_2D];

#if T8_WITH_OCC

/**
 * Definition of an occ geometry function.
 * This function maps reference coordinates to physical
 * coordinates regarding the occ geometries linked to the cells edges and faces.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 * \param [in]  tree_data   The data of the current tree as loaded by a \ref t8_geom_load_tree_data_fn.
 */
typedef void        (*t8_geom_occ_fn) (t8_cmesh_t cmesh,
                                       t8_gloidx_t gtreeid,
                                       const double *ref_coords,
                                       double out_coords[3],
                                       const void *tree_data);

/* *INDENT-OFF* */
class t8_geometry_occ:public t8_geometry_w_vertices, public t8_cad_geom
{
public:

  /**
   * Constructor of the occ geometry class.
   * \param [in] dimension  The dimension of this geometry.
   * \param [in] fileprefix Prefix of a .brep file from which to extract an occ geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_occ (int dimension, const char *fileprefix, const char *name);

  /**
   * Constructor of the occ geometry class.
   * \param [in] dimension  The dimension of this geometry.
   * \param [in] shape      Occ shape geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_occ (int dimension, const TopoDS_Shape shape, 
                   const char *name);

  /** The destructor. 
   * Clears the allocated memory.
   */
  virtual ~t8_geometry_occ ()
  {
    /* Nothing to do. */
  }

  /**
   * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  virtual void
  t8_geom_evaluate (t8_cmesh_t cmesh,
                    t8_gloidx_t gtreeid,
                    const double *ref_coords,
                    double out_coords[3]) const;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   */
  virtual void
  t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
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
  virtual void 
  t8_geom_load_tree_data (t8_cmesh_t cmesh,
                          t8_gloidx_t gtreeid);

private:
  const int                 *edges;  /**< The linked edges of the currently active tree. */
  const int                 *faces;  /**< The linked faces of the currently active tree. */
};
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */

#endif /* !T8_GEOMETRY_OCC_HXX! */
