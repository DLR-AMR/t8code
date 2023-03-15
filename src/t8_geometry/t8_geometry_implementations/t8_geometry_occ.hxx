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
struct t8_geometry_occ:public t8_geometry_w_vertices
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
   * \param [in] occ_shape  Occ shape geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_occ (int dimension, const TopoDS_Shape occ_shape, 
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
  
  /** Get an occ point from the occ_shape.
   * \param [in] index      The index of the point in the occ_shape.
   * \return                The occ point.
   */
  const gp_Pnt
  t8_geom_get_occ_point (const int index) const;
  
  /** Get an occ curve from the occ_shape.
   * \param [in] index      The index of the curve in the occ_shape.
   * \return                The occ curve.
   */
  const Handle_Geom_Curve
  t8_geom_get_occ_curve (const int index) const;

  /** Get an occ surface from the occ_shape.
   * \param [in] index      The index of the surface in the occ_shape.
   * \return                The occ surface.
   */
  const Handle_Geom_Surface
  t8_geom_get_occ_surface (const int index) const;

  /** Get the occ_shape_vertex2edge_map.
   * \return                The occ_shape_vertex_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_occ_shape_vertex_map() const;

  /** Get the occ_shape_edge2face_map.
   * \return                The occ_shape_edge_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_occ_shape_edge_map() const;

  /** Get the occ_shape_face_map.
   * \return                The occ_shape_face_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_occ_shape_face_map() const;

  /** Check if two occ points share a common occ edge.
   * \param [in]  vertex1_index  The index of the first occ point.
   * \param [in]  vertex2_index  The index of the second occ point.
   * \return                    Index of the shared edge. 0 if there is no shared edge.
   */
  int
  t8_geom_get_common_edge (const int vertex1_index, 
                           const int vertex2_index) const;

  /** Check if two occ edges share a common occ face.
   * \param [in]  edge1_index    The index of the first occ edge.
   * \param [in]  edge2_index    The index of the second occ edge.
   * \return                    Index of the shared face. 0 if there is no shared face.
   */
  int
  t8_geom_get_common_face (const int edge1_index, 
                           const int edge2_index) const;

  /** Check if a occ vertex lies on an occ edge.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  edge_index     The index of the occ edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_geom_is_vertex_on_edge (const int vertex_index, 
                             const int edge_index) const;

  /** Check if a occ vertex lies on an occ edge.
   * \param [in]  edge_index     The index of the occ vertex.
   * \param [in]  face_index     The index of the occ edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_geom_is_edge_on_face (const int edge_index, 
                           const int face_index) const;

  /** Check if a occ vertex lies on an occ face.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  face_index     The index of the occ face.
   * \return                    1 if vertex lies on face, otherwise 0.
   */
  int
  t8_geom_is_vertex_on_face (const int vertex_index, 
                             const int face_index) const;

  /** Retrieves the parameter of an occ vertex on an occ edge.
   *  The vertex has to lie on the edge.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  edge_index     The index of the occ edge.
   * \param [out] edge_param     The parameter of the vertex on the edge.
   */
  void
  t8_geom_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                          const int edge_index, 
                                          double* edge_param) const;

  /** Retrieves the parameters of an occ vertex on a occ face.
   *  The vertex has to lie on the face.
   * \param [in]  vertex_index   The index of the occ vertex.
   * \param [in]  face_index     The index of the occ face.
   * \param [out] face_params    The parameters of the vertex on the face.
   */                                   
  void 
  t8_geom_get_parameters_of_vertex_on_face(const int vertex_index, 
                                           const int face_index, 
                                           double* face_params) const;

  /** Converts the parameter of an occ edge to the corresponding parameters on an occ face.
   *  The edge has to lie on the face.
   * \param [in]  edge_index     The index of the occ edge.
   * \param [in]  face_index     The index of the occ face.
   * \param [in]  edge_param     The parameter on the edge.
   * \param [out] face_params    The corresponding parameters on the face.
   */   
  void
  t8_geom_edge_parameter_to_face_parameters(const int edge_index, 
                                            const int face_index, 
                                            const double edge_param, 
                                            double* face_params) const;

private:
  /**
   * Map a point in the reference space $$[0,1]^2$$ to $$\mathbb R^3$$. Only for quad trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of 2 entries, specifying a point in [0,1]^2.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void
  t8_geom_evaluate_occ_quad (t8_cmesh_t cmesh,
                             t8_gloidx_t gtreeid,
                             const double *ref_coords,
                             double out_coords[3]) const;

  /**
   * Map a point in the reference space $$[0,1]^3$$ to $$\mathbb R^3$$. Only for hex trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of 3 entries, specifying a point in [0,1]^3.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void
  t8_geom_evaluate_occ_hex (t8_cmesh_t cmesh,
                            t8_gloidx_t gtreeid,
                            const double *ref_coords,
                            double out_coords[3]) const;

  const int                                  *edges;                      /**< The linked edges of the currently active tree. */
  const int                                  *faces;                      /**< The linked faces of the currently active tree. */
  TopoDS_Shape                                occ_shape;                  /**< Occ geometry */
  TopTools_IndexedMapOfShape                  occ_shape_vertex_map;       /**< Map of all TopoDS_Vertex in shape. */
  TopTools_IndexedMapOfShape                  occ_shape_edge_map;         /**< Map of all TopoDS_Edge in shape. */
  TopTools_IndexedMapOfShape                  occ_shape_face_map;         /**< Map of all TopoDS_Face in shape. */
  TopTools_IndexedDataMapOfShapeListOfShape   occ_shape_vertex2edge_map;  /**< Maps all TopoDS_Vertex of shape to all its connected TopoDS_Edge */
  TopTools_IndexedDataMapOfShapeListOfShape   occ_shape_edge2face_map;    /**< Maps all TopoDS_Edge of shape to all its connected TopoDS_Face */
};
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */

#endif /* !T8_GEOMETRY_OCC_HXX! */
