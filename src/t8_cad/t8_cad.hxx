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

/** \file t8_cad.hxx
 * This file implements the t8_cad class. It manages OpenCASCADE shapes and implements
 * helper functions for working with the shapes.
 */

#ifndef T8_CAD_HXX
#define T8_CAD_HXX

#include <Standard_Handle.hxx>
#include <gp_Pnt.hxx>
#include <TopExp.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Curve.hxx>

/**
 * This class manages OpenCASCADE shapes and implements helper functions for working with the shapes.
*/
class t8_cad {
 public:
  /**
    * Constructor of the cad shape.
    * The shape is initialized based on a .brep file with the given prefix.
    * The internal structure extracts and stores geometric information such as
    * vertices, edges, and faces from this file. The number and type of vertices
    * should match the tree type (quad/hex/tri), and the cad data must be valid.
    * This constructor is intended for general use, including file-based mesh creation.
    *
    * \param [in] fileprefix  Prefix of a .brep file from which to extract cad geometry.
    */
  t8_cad (std::string fileprefix);

  /**
    * Constructor of the cad shape.
    * The shape is initialized directly from an existing TopoDS_Shape.
    * This constructor is especially useful for use in scripts, testing,
    * or integration with mesh generators that already provide geometry in memory.
    * It avoids file I/O and allows full control over the CAD input.
    *
    * \param [in] cad_shape  cad shape geometry object.
    */
  t8_cad (const TopoDS_Shape cad_shape);

  /**
   * Constructor of the cad shape for testing purposes. Sets an invalid cad_shape.
   */
  t8_cad ();

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
   * \param [in]  surface_params The parameters of the surface nodes.
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
   * \param [in]  surface_index   The index of the cad face.
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
  TopoDS_Shape cad_shape;                          /**< cad geometry */
  TopTools_IndexedMapOfShape cad_shape_vertex_map; /**< Map of all TopoDS_Vertex in shape. */
  TopTools_IndexedMapOfShape cad_shape_edge_map;   /**< Map of all TopoDS_Edge in shape. */
  TopTools_IndexedMapOfShape cad_shape_face_map;   /**< Map of all TopoDS_Face in shape. */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_vertex2edge_map; /**< Maps all TopoDS_Vertex of shape to all its connected TopoDS_Edge */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_edge2face_map; /**< Maps all TopoDS_Edge of shape to all its connected TopoDS_Face */
};

#endif /* !T8_CAD_HXX */
