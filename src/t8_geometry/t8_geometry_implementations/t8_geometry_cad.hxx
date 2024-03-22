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

/** \file t8_geometry_cad.hxx
 * This geometry implements OpenCASCADE geometries. It enables the option to link different 
 * 1 and 2 dimensional cad geometries to the edges and faces of refinement trees. 
 * The geometry of the refinement tree is extended into the volume accordingly.
 */

#ifndef T8_GEOMETRY_CAD_HXX
#define T8_GEOMETRY_CAD_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.h>

#if T8_WITH_OCC

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
   * Constructor of the cad geometry with a given dimension. The geometry
   * is currently viable with quad/hex and triangle trees. Tets will be supported soon.
   * The geometry uses as many vertices as the tree type has, as well as
   * additional geometry information, which is extracted from a .brep file.
   * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
   * Since the internals of this geometry are finely tuned to the .brep file
   * it is recommended to only use it with the \ref t8_cmesh_readmshfile function.
   * \param [in] dim        The dimension of this geometry.
   * \param [in] fileprefix Prefix of a .brep file from which to extract an cad geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_cad (int dim, std::string fileprefix, std::string name = "t8_geom_cad");

  /**
   * Constructor of the cad geometry with a given dimension. The geometry
   * is currently viable with quad/hex and triangle trees. Tets will be supported soon.
   * The geometry uses as many vertices as the tree type has, as well as
   * additional geometry information, which is given via the \a cad_shape.
   * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
   * This constructor can be used in short scripts or in combination with a
   * mesh generator, to omit the file IO of the 
   * \ref t8_geometry_cad (int dim, std::string fileprefix,  std::string name) constructor.
   * \param [in] dim        The dimension of this geometry.
   * \param [in] cad_shape  cad shape geometry.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_cad (int dim, const TopoDS_Shape cad_shape, std::string name = "t8_geom_cad");

  /**
   * Constructor of the cad geometry for testing purposes. Sets an invalid cad_shape.
   * \param [in] dim        The dimension of this geometry.
   */
  t8_geometry_cad (int dim);

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
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  virtual void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
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

  /** Check if a cad_curve is a line.
   * \param [in] curve_index      The index of the cad_curve.
   * \return                      1 if curve is a line, 0 if curve is not a line.
   */
  int
  t8_geom_is_line (const int curve_index) const;

  /** Check if a cad_surface is a plane.
   * \param [in] surface_index      The index of the cad_surface.
   * \return                        1 if surface is a plane linear, 0 if surface is not a plane.
   */
  int
  t8_geom_is_plane (const int surface_index) const;

  /** Get an cad point from the cad_shape.
   * \param [in] index      The index of the point in the cad_shape.
   * \return                The cad point.
   */
  const gp_Pnt
  t8_geom_get_cad_point (const int index) const;

  /** Get an cad curve from the cad_shape.
   * \param [in] index      The index of the curve in the cad_shape.
   * \return                The cad curve.
   */
  const Handle_Geom_Curve
  t8_geom_get_cad_curve (const int index) const;

  /** Get an cad surface from the cad_shape.
   * \param [in] index      The index of the surface in the cad_shape.
   * \return                The cad surface.
   */
  const Handle_Geom_Surface
  t8_geom_get_cad_surface (const int index) const;

  /** Get the cad_shape_vertex2edge_map.
   * \return                The cad_shape_vertex_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_cad_shape_vertex_map () const;

  /** Get the cad_shape_edge2face_map.
   * \return                The cad_shape_edge_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_cad_shape_edge_map () const;

  /** Get the cad_shape_face_map.
   * \return                The cad_shape_face_map.
   */
  const TopTools_IndexedMapOfShape
  t8_geom_get_cad_shape_face_map () const;

  /** Check if two cad points share a common cad edge.
   * \param [in]  vertex1_index  The index of the first cad point.
   * \param [in]  vertex2_index  The index of the second cad point.
   * \return                    Index of the shared edge. 0 if there is no shared edge.
   */
  int
  t8_geom_get_common_edge (const int vertex1_index, const int vertex2_index) const;

  /** Check if two cad edges share a common cad face.
   * \param [in]  edge1_index    The index of the first cad edge.
   * \param [in]  edge2_index    The index of the second cad edge.
   * \return                    Index of the shared face. 0 if there is no shared face.
   */
  int
  t8_geom_get_common_face (const int edge1_index, const int edge2_index) const;

  /** Check if a cad vertex lies on an cad edge.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  edge_index     The index of the cad edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_geom_is_vertex_on_edge (const int vertex_index, const int edge_index) const;

  /** Check if a cad vertex lies on an cad edge.
   * \param [in]  edge_index     The index of the cad vertex.
   * \param [in]  face_index     The index of the cad edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  t8_geom_is_edge_on_face (const int edge_index, const int face_index) const;

  /** Check if a cad vertex lies on an cad face.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  face_index     The index of the cad face.
   * \return                    1 if vertex lies on face, otherwise 0.
   */
  int
  t8_geom_is_vertex_on_face (const int vertex_index, const int face_index) const;

  /** Retrieves the parameter of an cad vertex on an cad edge.
   *  The vertex has to lie on the edge.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  edge_index     The index of the cad edge.
   * \param [out] edge_param     The parameter of the vertex on the edge.
   */
  void
  t8_geom_get_parameter_of_vertex_on_edge (const int vertex_index, const int edge_index, double *edge_param) const;

  /** Retrieves the parameters of an cad vertex on a cad face.
   *  The vertex has to lie on the face.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  face_index     The index of the cad face.
   * \param [out] face_params    The parameters of the vertex on the face.
   */
  void
  t8_geom_get_parameters_of_vertex_on_face (const int vertex_index, const int face_index, double *face_params) const;

  /** Converts the parameters of an cad edge to the corresponding parameters on an cad face.
   * The edge has to lie on the face.
   * For the conversion of edge parameters of mesh elements to topological face parameters of a closed surface, it is additionally
   * checked, whether the conversion was correct, to prevent disorted elements. 
   * \param [in]  edge_index     The index of the cad edge, which parameters should be converted to face parameters.
   * \param [in]  face_index     The index of the cad face, on to which the edge parameters should be converted.
   * \param [in]  num_face_nodes The number of the face nodes of the evaluated element. Only needed for closed surface check, otherwise NULL.
   * \param [in]  edge_param     The parameter on the edge.
   * \param [in]  surface_param  The parameters of the surface nodes.
   *                             When provided, there are additional checks for closed geometries.
   *                             If there are no surface parameter
   *                             to pass in to the function, you can pass NULL.
   * \param [out] face_params    The corresponding parameters on the face.
   */
  void
  t8_geom_edge_parameter_to_face_parameters (const int edge_index, const int face_index, const int num_face_nodes,
                                             const double edge_param, const double *surface_params,
                                             double *face_params) const;

  /** Finds the parametric bounds of an cad face.
   * \param [in]  face_index   The index of the cad face.
   * \param [out] bounds          The parametric bounds of the cad face.
   */
  void
  t8_geom_get_face_parametric_bounds (const int surface_index, double *bounds) const;

  /** Finds the parametric bounds of an cad edge.
   * \param [in]  edge_index   The index of the cad edge.
   * \param [out] bounds       The parametric bounds of the cad edge.
   */
  void
  t8_geom_get_edge_parametric_bounds (const int edge_index, double *bounds) const;

  /** Checks if an edge is closed in its U parameter.
   * \param [in]  edge_index   The index of the closed edge.
   * \return                   1 if edge is closed in U. 0 if edge is not closed in U.
   */
  int
  t8_geom_is_edge_closed (int edge_index) const;

  /** Checks if a surface is closed in its U parameter or V parameter.
   * \param [in]  geometry_index   The index of the closed geometry.
   * \param [in]  parameter        The parameter, which should be check for closeness.
   *                               0 stands for the U parameter and 1 for the V parameter. 
   * \return                       1 if geometry is closed in U. 0 if geometry is not closed in U.
   */
  int
  t8_geom_is_surface_closed (int geometry_index, int parameter) const;

 private:
  /**
   * Maps points in the reference space \f$ [0,1]^2 \f$ to \f$ \mathbb{R}^3 \f$. Only for triangle trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of 2 entries, specifying a point in \f$ [0,1]^2 \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  void
  t8_geom_evaluate_cad_tri (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                            double *out_coords) const;

  /**
   * Maps points in the reference space \f$ [0,1]^2 \f$ to \f$ \mathbb{R}^3 \f$. Only for quad trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_quad (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *out_coords) const;

  /**
   * Map a point in the reference space $$[0,1]^3$$ to $$\mathbb R^3$$. Only for tet trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
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
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_hex (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                            double *out_coords) const;

  /**
   * Maps points in the reference space \f$ \f$ [0,1]^3 \f$ \f$ to \f$ \mathbb{R}^3 \f$. Only for prism trees.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  void
  t8_geom_evaluate_cad_prism (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                              double *out_coords) const;

  const int *edges;                                /**< The linked edges of the currently active tree. */
  const int *faces;                                /**< The linked faces of the currently active tree. */
  TopoDS_Shape cad_shape;                          /**< cad geometry */
  TopTools_IndexedMapOfShape cad_shape_vertex_map; /**< Map of all TopoDS_Vertex in shape. */
  TopTools_IndexedMapOfShape cad_shape_edge_map;   /**< Map of all TopoDS_Edge in shape. */
  TopTools_IndexedMapOfShape cad_shape_face_map;   /**< Map of all TopoDS_Face in shape. */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_vertex2edge_map; /**< Maps all TopoDS_Vertex of shape to all its connected TopoDS_Edge */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_edge2face_map; /**< Maps all TopoDS_Edge of shape to all its connected TopoDS_Face */
};

#endif /* T8_WITH_OCC */

#endif /* !T8_GEOMETRY_CAD_HXX */
