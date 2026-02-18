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

/** \file t8_cad_handle.hxx
 * This file implements the t8_cad_handle class. It manages OpenCASCADE shapes and implements
 * helper functions for working with the shapes.
 */

#ifndef T8_CAD_HANDLE_HXX
#define T8_CAD_HANDLE_HXX

#include <gp_Pnt.hxx>
#include <TopExp.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <optional>
#include <span>

#include <t8_cmesh/t8_cmesh.h>
/**
 * This class manages OpenCASCADE shapes and implements helper functions for working with the shapes.
*/
class t8_cad_handle {
 public:
  /**
    * Constructor of the cad shape handle.
    * The shape is initialized based on a .brep file with the given prefix.
    * The internal structure extracts and stores geometric information such as
    * vertices, edges, and faces from this file. The number and type of vertices
    * should match the tree type (quad/hex/tri), and the cad data must be valid.
    * This constructor is intended for general use, including file-based mesh creation.
    *
    * \param [in] fileprefix  Prefix of a .brep file from which to extract cad geometry.
    */
  t8_cad_handle (std::string fileprefix);
  /**
    * Constructor of the cad shape.
    * The shape is initialized directly from an existing TopoDS_Shape.
    * This constructor is especially useful for use in scripts, testing,
    * or integration with mesh generators that already provide geometry in memory.
    * It avoids file I/O and allows full control over the CAD input.
    *
    * \param [in] cad_shape  cad shape geometry object.
    */
  t8_cad_handle (const TopoDS_Shape cad_shape);

  /**
   * Constructor of the cad shape for testing purposes. Sets an invalid cad_shape.
   */
  t8_cad_handle ();

  /**
 * Destructor of the cad shape handle.
 */
  ~t8_cad_handle ();

  /**
   * Loads a cad shape from a .brep file with the given prefix and maps it.
   * \param [in] fileprefix  Prefix of a .brep file from which to extract cad geometry.
   */
  void
  load_cad_from_file (const std::string fileprefix);

  /**
   * Loads a cad shape from an existing TopoDS_Shape and maps it.
   * \param [in] new_cad_shape  The input cad shape.
   */
  void
  load_cad_from_shape (const TopoDS_Shape &cad_shape);

  /** Check if a cad_curve is a line.
   * \param [in] curve_index      The index of the cad_curve.
   * \return                      1 if curve is a line, 0 if curve is not a line.
   */
  int
  is_line (const int curve_index) const;

  /** Check if a cad_surface is a plane.
   * \param [in] surface_index      The index of the cad_surface.
   * \return                        1 if surface is a plane linear, 0 if surface is not a plane.
   */
  int
  is_plane (const int surface_index) const;

  /** Get a cad vertex from the cad_shape.
   * \param [in] index      The index of the vertex in the cad_shape.
   * \return                The cad vertex.
   */
  const TopoDS_Vertex
  get_cad_vertex (const int index) const;

  /** Get a cad edge from the cad_shape.
   * \param [in] index      The index of the edge in the cad_shape.
   * \return                The cad edge.
   */
  const TopoDS_Edge
  get_cad_edge (const int index) const;

  /** Get a cad face from the cad_shape.
   * \param [in] index      The index of the face in the cad_shape.
   * \return                The cad face.
   */
  const TopoDS_Face
  get_cad_face (const int index) const;

  /** Get a cad point from the cad_shape.
   * \param [in] index      The index of the point in the cad_shape.
   * \return                The cad point.
   */
  const gp_Pnt
  get_cad_point (const int index) const;

  /** Get a cad curve from the cad_shape.
   * \param [in] index      The index of the curve in the cad_shape.
   * \return                The cad curve.
   */
  const Handle_Geom_Curve
  get_cad_curve (const int index) const;

  /** Get a cad surface from the cad_shape.
   * \param [in] index      The index of the surface in the cad_shape.
   * \return                The cad surface.
   */
  const Handle_Geom_Surface
  get_cad_surface (const int index) const;

  /** Get the cad_shape_vertex2edge_map.
   * \return                The cad_shape_vertex_map.
   */
  const TopTools_IndexedMapOfShape
  get_cad_shape_vertex_map () const;

  /** Get the cad_shape_edge2face_map.
   * \return                The cad_shape_edge_map.
   */
  const TopTools_IndexedMapOfShape
  get_cad_shape_edge_map () const;

  /** Get the cad_shape_face_map.
   * \return                The cad_shape_face_map.
   */
  const TopTools_IndexedMapOfShape
  get_cad_shape_face_map () const;

  /** Check if two cad points share a common cad edge.
   * \param [in]  vertex1_index  The index of the first cad point.
   * \param [in]  vertex2_index  The index of the second cad point.
   * \return                    Index of the shared edge. 0 if there is no shared edge.
   */
  int
  get_common_edge_of_vertices (const int vertex1_index, const int vertex2_index) const;

  /** Check if two cad edges share a common cad face.
   * \param [in]  edge1_index    The index of the first cad edge.
   * \param [in]  edge2_index    The index of the second cad edge.
   * \return                    Index of the shared face. 0 if there is no shared face.
   */
  int

  get_common_face_of_edges (const int edge1_index, const int edge2_index) const;

  /** Check if a cad vertex and cad edge share a common face.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  edge_index     The index of the cad edge.
   * \return                     Index of the shared face. 0 if there is no shared face.
   */
  int
  get_common_face_of_vertex_and_edge (const int vertex_index, const int edge_index) const;

  /** Check if two cad vertices share a common cad face.
   * \param [in]  vertex1_index The index of the first cad edge.
   * \param [in]  vertex2_index The index of the second cad edge.
   * \return                    Index of the shared face. 0 if there is no shared face.
   */
  int
  get_common_face_of_vertices (const int vertex1_index, const int vertex2_index) const;

  /** Check if a cad vertex lies on an cad edge.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  edge_index     The index of the cad edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  is_vertex_on_edge (const int vertex_index, const int edge_index) const;

  /** Check if a cad vertex lies on an cad edge.
   * \param [in]  edge_index     The index of the cad vertex.
   * \param [in]  face_index     The index of the cad edge.
   * \return                    1 if vertex lies on edge, otherwise 0.
   */
  int
  is_edge_on_face (const int edge_index, const int face_index) const;

  /** Check if a cad vertex lies on an cad face.
   * \param [in]  vertex_index   The index of the cad vertex.
   * \param [in]  face_index     The index of the cad face.
   * \return                    1 if vertex lies on face, otherwise 0.
   */
  int
  is_vertex_on_face (const int vertex_index, const int face_index) const;

  /** Returns true if \a vertex_index is on a seam of \a edge_index.
   * A seam is a vertex which connects a curve to itself.
   *
   * \param [in] vertex_index   The index of the cad vertex.
   * \param [in] edge_index     The index of the cad edge.
   * \return true if the vertex is a seam. false otherwise.
   */
  bool
  vertex_is_seam (const int vertex_index, const int edge_index) const;

  /** Returns true if \a vertex_index is on a seam of \a face_index.
   * A seam is an edge which connects a surface to itself.
   *
   * \param [in] vertex_index   The index of the cad vertex.
   * \param [in] face_index     The index of the cad face.
   * \return true if the edge is a seam. false otherwise.
   */
  bool
  vertex_is_on_seam_edge (const int vertex_index, const int face_index) const;

  /** Returns true if \a edge_index is a seam of \a face_index.
   * A seam is an edge which connects a surface to itself.
   *
   * \param [in] edge_index   The index of the cad edge.
   * \param [in] face_index   The index of the cad face.
   * \return true if the edge is a seam. false otherwise.
   */
  bool
  edge_is_seam (const int edge_index, const int face_index) const;

  /** Retrieves the parameter of an cad vertex on an cad edge.
   * The vertex has to lie on the edge.
   * \warning If the edge is closed in any direction and the vertex is the closing bound,
   * it is random which side of the closed edge the parameter is from. The parameter
   * is correct, but for curved elements it has to be checked if the parameter has to be
   * converted onto the other bound.
   * \param [in]  vertex_index            The index of the cad vertex.
   * \param [in]  edge_index              The index of the cad edge.
   * \param [out] edge_param              The parameter of the vertex on the edge.
   * \param [in]  reference_edge_param    Reference parameters on the edge.
   */
  void
  get_parameter_of_vertex_on_edge (const int vertex_index, const int edge_index, double *edge_param,
                                   std::optional<double> reference_edge_param = std::nullopt) const;

  /** Retrieves the parameters of an cad vertex on a cad face.
   * The vertex has to lie on the face.
   * If the vertex is in a seam, the vertex closer to the \a reference_face_params is chosen.
   * \warning If the face is closed in any direction and the vertex is on the seam and no reference
   * parameters are provided it is random which side of the closed face the parameters are from.
   * The parameters are correct, but for curved elements it has to be checked if the parameters have to be
   * converted onto the other bound.
   * \param [in]  vertex_index              The index of the cad vertex.
   * \param [in]  face_index                The index of the cad face.
   * \param [out] face_params               The parameters of the vertex on the face.
   * \param [in]  reference_face_params     Reference parameters on the surface.
   */
  void
  get_parameters_of_vertex_on_face (const int vertex_index, const int face_index, double face_params[2],
                                    std::optional<std::span<const double, 2> > reference_face_params
                                    = std::nullopt) const;

  /** Converts the parameters of a cad edge to the corresponding parameters on a cad face.
   * The edge has to lie on the face.
   * If the edge is a seam, the edge closer to the \a reference_face_params is chosen.
   * \warning If the face is closed in any direction and the edge is the seam and no reference
   * parameters are provided it is random which side of the closed face the parameters are from.
   * The parameters are correct, but for curved elements it has to be checked if the parameters have to be
   * converted onto the other side of the seam.
   * \param [in]  edge_index                The index of the cad edge, which parameters should be converted to face parameters.
   * \param [in]  face_index                The index of the cad face, on to which the edge parameters should be converted.
   * \param [in]  edge_param                The parameter on the edge.
   * \param [out] face_params_out           The corresponding parameters on the face.
   * \param [in]  reference_face_params     Reference parameters on the surface.
   */
  void
  edge_parameter_to_face_parameters (const int edge_index, const int face_index, const double edge_param,
                                     double face_params_out[2],
                                     std::optional<std::span<const double, 2> > reference_face_params
                                     = std::nullopt) const;

  /** Finds the parametric bounds of an cad face.
   * \param [in]  surface_index   The index of the cad face.
   * \param [out] bounds          The parametric bounds of the cad face.
   */
  void
  get_face_parametric_bounds (const int surface_index, double *bounds) const;

  /** Finds the parametric bounds of an cad edge.
   * \param [in]  edge_index   The index of the cad edge.
   * \param [out] bounds       The parametric bounds of the cad edge.
   */
  void
  get_edge_parametric_bounds (const int edge_index, double *bounds) const;

  /** Checks if an edge is closed.
   * \param [in]  edge_index   The index of the closed edge.
   * \return                   true if edge is closed
   */
  int
  edge_is_closed (int edge_index) const;

  /**
   * Increase the reference count of the cad handle.
   */
  inline void
  ref ()
  {
    t8_refcount_ref (&rc);
  }

  /**
   * Decrease the reference count of the cad handle.
   * If the reference count reaches zero, the cad handle is deleted.
   */
  inline void
  unref ()
  {
    if (t8_refcount_unref (&rc)) {
      t8_debugf ("Deleting the cad_handle.\n");
      delete this;
    }
  }

  /** Checks if a surface is closed in any direction.
   * \param [in]  geometry_index   The index of the closed geometry.
   * \return                       true if geometry is closed in any direction.
   */

  bool
  surface_is_closed (int geometry_index) const;

  /** Checks if a surface is closed in its U direction or V direction.
   * \param [in]  geometry_index   The index of the closed geometry.
   * \param [in]  direction        The direction, which should be check for closeness.
   *                               0 stands for the U direction and 1 for the V direction.
   * \return                       true if geometry is closed in the given direction.
   */
  bool
  surface_is_closed (int geometry_index, int direction) const;

 private:
  /**
   * Map the cad shape to extract vertices, edges, faces, and their relationships.
   * \param [in] cad_shape_in  The input cad shape to be mapped.
   */
  void
  map_cad_shape (const TopoDS_Shape &cad_shape_in);

  TopoDS_Shape cad_shape;                          /**< cad geometry */
  TopTools_IndexedMapOfShape cad_shape_vertex_map; /**< Map of all TopoDS_Vertex. */
  TopTools_IndexedMapOfShape cad_shape_edge_map;   /**< Map of all TopoDS_Edge. */
  TopTools_IndexedMapOfShape cad_shape_face_map;   /**< Map of all TopoDS_Face. */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_vertex2edge_map; /**< Maps all TopoDS_Vertex to all its connected TopoDS_Edge */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_edge2face_map; /**< Maps all TopoDS_Edge to all its connected TopoDS_Face */
  /** The reference count of the cad handle. TODO: Replace by shared_ptr when cmesh becomes a class. */
  TopTools_IndexedDataMapOfShapeListOfShape
    cad_shape_vertex2face_map; /**< Maps all TopoDS_Vertex to all its connected TopoDS_Face */
  t8_refcount_t rc;
};

#endif /* !T8_CAD_HANDLE_HXX */
